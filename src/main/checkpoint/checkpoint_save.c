/* checkpoint_save.c */
/* Wolfgang Tichy, 8/2019 */


#include "nmesh.h"
#include "checkpoint.h"


/******************************************************************/
/* some functions to save nmesh data for checkpoints  */
/******************************************************************/

/******************************************************************/
/* functions to save pars */
/******************************************************************/
/* save pars */
int checkpoint_save_pars(tMesh *mesh, char *fname)
{
  FILE *fp;
  char *list, *saveptr, *name;
  int IObufsz = Geti(Par("fwrite_bufsize"));
  char *IObuf; /* larger buffer for write */

  /* only Rank0 writes the file */
  if(!Rank0) return 0;

  /* open destination file */
  fp = fopen_buf(fname, "wb", &IObuf,IObufsz);
  if(!fp) errorexits("failed opening %s", fname);

  fprintf(fp, "# parameters listed in checkpoint_save_pars\n");

  /* put vals of checkpoint_save_pars in list,
     duplicate Gets(Par("checkpoint_save_pars")) because strtok_r
     will modify list */
  list = strdup( Gets(Par("checkpoint_save_pars")) );

  /* loop over contents of list, and print pars */
  for(name=strtok_r(list, " ", &saveptr); name!=NULL;
      name=strtok_r(NULL, " ", &saveptr))
  {
    fprintf(fp, "%s = %s\n", name, Gets(Par(name)));
  }

  free(list);
  fclose_buf(fp, &IObuf);
  return 0;
}


/******************************************************************/
/* functions to save patches */
/******************************************************************/
/* save patch info */
int checkpoint_save_patches(tMesh *mesh, char *fname)
{
  FILE *fp;
  int p;
  int IObufsz = Geti(Par("fwrite_bufsize"));
  char *IObuf; /* larger buffer for write */

  /* only Rank0 writes the file */
  if(!Rank0) return 0;

  /* open destination file */
  fp = fopen_buf(fname, "wb", &IObuf,IObufsz);
  if(!fp) errorexits("failed opening %s", fname);

  /* header with some mesh info */
  fprintf(fp, "mesh->\n");
  fprintf(fp, " time = %.19g\n", mesh->time);
  fprintf(fp, " iteration = %d\n", mesh->iteration);
  fprintf(fp, " dt = %.19g\n", mesh->dt);
  fprintf(fp, "\n");

  /* write data of each patch */
  forpatches(mesh, p)
  {
    tPat *pat = mesh->pat[p];

    checkpoint_write_pat(fp, pat);
    fprintf(fp, "\n");
  }

  fclose_buf(fp, &IObuf);
  return 0;
}

/* write non-pointer part of tPat */
void checkpoint_write_pat(FILE *fp, tPat *pat)
{
  int f;

  fprintf(fp, "patch%d:\n", pat->p);
  fprintf(fp, "pat->\n");

  for(f=0; f<6; f++)
    fprintf(fp, " bbox[%d] = %.19g\n", f, pat->bbox[f]);

  //fprintf(fp,   " nmax = %d\n", pat->nmax);

  //for(d=0; d<3; d++)
  //  fprintf(fp, " periodic[%d] = %d\n", d, pat->periodic[d]);

  //for(d=0; d<3; d++)
  //  fprintf(fp, " rnode->n[%d] = %d\n", d, pat->rnode->n[d]);

  //for(d=0; d<3; d++)
  //  fprintf(fp, " rnode->pt_typ[%d] = %d\n", d, pat->rnode->pt_typ[d]);

  //printCI(pat);
  checkpoint_write_CI(fp, pat->CI);
}


/* write non-pointer part of tCoordInfo */
void checkpoint_write_CI(FILE *fp, tCoordInfo *CI)
{
  int d, f, i, useF=0;

  fprintf(fp, " CI->\n");
  /*
  for(f=0; f<6; f++)
    fprintf(fp,   "  iSurf[%d] = %d\n", f, CI->iSurf[f]);

  for(f=0; f<6; f++)
    for(d=0; d<3; d++)
      fprintf(fp, "  idSurfdX[%d][%d] = %d\n", f,d, CI->idSurfdX[f][d]);
  */
  for(f=0; f<6; f++)
    fprintf(fp,   "  s[%d] = %.19g\n", f, CI->s[f]);

  for(d=0; d<3; d++)
    fprintf(fp,   "  xc[%d] = %.19g\n", d, CI->xc[d]);

  for(i=0; i<4; i++)
    fprintf(fp,   "  co[%d] = %.19g\n", i, CI->co[i]);

  fprintf(fp,     "  dom = %d\n", CI->dom);
  fprintf(fp,     "  type = %d\n", CI->type);

  for(f=0; f<6; f++) if(CI->FSurf[f]) { useF = 1; break; }
  fprintf(fp,     "  use_FSurf = %d\n", useF);

  /* this signifies end of patch info, so this needs to be last */
  fprintf(fp,     "  label = %d\n", CI->label);
}


/******************************************************************/
/* functions to save elms */
/******************************************************************/
/* save elm info */
int checkpoint_save_elms(tMesh *mesh, char *fname)
{
  FILE *fp = NULL;
  int IObufsz = Geti(Par("fwrite_bufsize"));
  char *IObuf; /* larger buffer for write */

  /* only Rank0 writes the file */
  if(Rank0)
  {
    /* open destination file */
    fp = fopen_buf(fname, "wb", &IObuf,IObufsz);
    if(!fp) errorexits("failed opening %s", fname);

    fprintf(fp, "number of elms, followed by all elms, their n[3] and "
                "optionally their py_typ[3]\n\n");
    fprintf(fp, "nelms = %lu\n\n", mesh->eidlim[nMPI_size()-1]);
  }

  /* write all nodes */
  checkpoint_write_elms(mesh, fp);

  /* Rank0 needs to close file */
  if(Rank0)
    fclose_buf(fp, &IObuf);

  return 0;
}

/* write out info about all elms */
int checkpoint_write_elms(tMesh *mesh, FILE *fp)
{
  int size = nMPI_size();
  int rk;
  for(rk=0; rk<size; rk++)
  {
    int n_def[3] = {0};      /* defaults for n */
    int pt_typ_def[3] = {0}; /* and pt_typ     */
    ulong nelm0s;
    tElm0 *elm0 = amr_alloc_get_elm0array_of_rank(mesh, rk, &nelm0s);
    ulong i;

    for(i=0; i<nelm0s; i++)
    {
      int n[3], pt_typ[3];
      int d, l, write_n, write_pt_typ;
      tEloc eloc[1];

      /* get n, pt_typ from elm and decide if we write them as well */
      write_n = write_pt_typ = 0;
      for(d=0; d<3; d++)
      {
        n[d]      = elm0[i].n[d];
        pt_typ[d] = elm0[i].pt_typ[d];

        if(n[d]      != n_def[d])      { write_n = 1; }
        if(pt_typ[d] != pt_typ_def[d]) { write_n = 1; write_pt_typ = 1; }
      }

      /* update defaults */
      for(d=0; d<3; d++)
      {
        n_def[d]      = n[d];
        pt_typ_def[d] = pt_typ[d];
      }

      /* get loc strings in eloc */
      eloc_from_eploc(eloc, elm0[i].eploc);

      /* only rank0 writes */
      if(Rank0)
      {
        /* write elm name */
        fprintf(fp, "%d_", eloc->p);
        for(l = 0; l < eloc->l; l++) fputc(eloc->loc[l], fp);
        fprintf(fp, "\n");

        /* add info about n, pt_typ */
        for(d=0; d<3; d++)
        {
          if(write_n)                 fprintf(fp, "%d", n[d]);
          if(write_pt_typ)            fprintf(fp, " %d", pt_typ[d]);
          if(write_n || write_pt_typ) fprintf(fp, "\n");
        }
        fprintf(fp, "\n");
      }
    } /* end rk-loop */
    free(elm0);
  }
  return 0;
}


/******************************************************************/
/* functions to save variables */
/******************************************************************/

/* save all EvoVars */
int checkpoint_save_EvoVars(tMesh *mesh, char *fname)
{
  tVarList *vl = vlalloc(mesh);
  int vi;

  /* loop over all vars and put all important EvoVars into vl */
  for(vi=0; vi<mesh->nvdb; vi++)
  {
    int vt = MeshVarType(mesh, vi);
    if( (vt==EVOVAR) || (vt==DATAVAR) )
      if(!var_added_by_evolve_init_evosys(mesh, vi))
        vlpushone(vl, vi);
  }

  /* We make variables.bin portable by always writing the double arrays
     with the variable data in little endian format. */
  /* write var list vl in little endian format */
  checkpoint_save_VL(mesh, fname, vl, 0);
  vlfree(vl);
  return 0;
}

/* save all amr_elm_nbinfo vars */
int checkpoint_save_nbinfoVars(tMesh *mesh, char *fname)
{
  tVarList *vl = vlalloc(mesh);

  vlpush(vl, Ind("amr_elm_nbinfo0"));

  /* nbinfo.bin is not portable, because we write an array of structs.
     Therefore we simply save it in native byte order. NOTE: If a checkpoint
     is moved between machines nbinfo.bin should be deleted! */
  /* write var list vl in native format */
  checkpoint_save_VL(mesh, fname, vl, 1);
  vlfree(vl);
  return 0;
}

/* write a varlist vl */
int checkpoint_save_VL(tMesh *mesh, char *fname, tVarList *vl,
                       int write_native)
{
  FILE *fp;
  int rk;
  int IObufsz = Geti(Par("fwrite_bufsize"));
  char *IObuf; /* larger buffer for write */

  /* MPI motivated loop to assign work */
  for(rk=0; rk<nMPI_size(); rk++)
  {
    /* do work when it is my turn */
    if(rk == nMPI_rank())
    {
      /* open destination file */
      if(Rank0) fp = fopen_buf(fname, "wb", &IObuf,IObufsz);
      else      fp = fopen_buf(fname, "ab", &IObuf,IObufsz);
      if(!fp) errorexits("failed opening %s", fname);

      /* write var list vl in native or little endian format */
      checkpoint_write_vl(fp, vl, write_native);

      fclose_buf(fp, &IObuf);
      fs_sync(mesh); /* make sure every MPI proc flushes buffers to disk */
    }
    /* wait until everyone is here */
    nMPI_barrier();
  } /* end rk-loop */

  return 0;
}

/* output varlist on each node */
void checkpoint_write_vl(FILE *fp, tVarList *vl, int write_native)
{
  tMesh *mesh = vl->mesh;
  char name[256];
  int vli;

  /* write all var names */
  if(Rank0)
  {
    fprintf(fp, "number of variables and variable list in this file:\n");
    fprintf(fp, "%d\n", vl->n);
    forvl(vl, vli)
    {
      int vi = Vind(vl, vli);
      fprintf(fp, "%s\n", VarName(vi));
    }
    fprintf(fp, "\n\n");
    fprintf(fp, "variable-list-data:\n\n");
  }

  /* loop over all leaf nodes */
  formylnodes_noomp(mesh)
  {
    tNode *node = MyLnode;

    /* node name and number of points np on this node */
    nodename(node, name,255);
    fprintf(fp, "{\n");
    fprintf(fp, "%s\n", name);
    fprintf(fp, "%d\n", node->np);

    /* write data for all in vl on this node */
    forvl(vl, vli)
    {
      /* do something only if this proc has dat */
      if(node->dat)
      {
        int vi = Vind(vl, vli);
        tArray *va = node->dat->v[vi];
        if(va)
        {
          double *v = Arrd_(va);

          if(va->N == node->np) /* print only index in varlist */
            fprintf(fp, "%d\n", vli);
          else /* print index in varlist and array dims */
            fprintf(fp, "%d %d %d %d\n", vli, va->n[0],va->n[1],va->n[2]);

          /* write var array in raw binary */
          if(write_native) fwrite(v, sizeof(double), va->N, fp);
          else             fwrite_little(v, sizeof(double), va->N, fp);
          fprintf(fp, "\n");
          //if(Node_eid(node)==28) printf("v[]=%g\n", v[0]);
        }
      }
    }
    fprintf(fp, "}\n");
    fprintf(fp, "\n");
  } /* end node-loop */
}
