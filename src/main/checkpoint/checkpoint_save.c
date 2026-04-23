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
int checkpoint_save_patches(tMesh *mesh, char *dir, char *fname)
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

    checkpoint_write_pat(fp, pat, dir);
    fprintf(fp, "\n");
  }

  fclose_buf(fp, &IObuf);
  return 0;
}

/* write non-pointer part of tPat */
void checkpoint_write_pat(FILE *fp, tPat *pat, const char *dir)
{
  int f;

  fprintf(fp, "patch%d:\n", pat->p);
  fprintf(fp, "pat->\n");

  for(f=0; f<6; f++)
    fprintf(fp, " bbox[%d] = %.19g\n", f, pat->bbox[f]);

  fprintf(fp,   " npg = %d\n", pat->npg);
  fprintf(fp,   " pg0 = %d\n", pat->pg0);

  //for(d=0; d<3; d++)
  //  fprintf(fp, " periodic[%d] = %d\n", d, pat->periodic[d]);

  //printCI(pat);
  checkpoint_write_CI(fp, pat, dir);
}


/* write non-pointer part of tCoordInfo */
void checkpoint_write_CI(FILE *fp, tPat *pat, const char *dir)
{
  tCoordInfo *CI = pat->CI;
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
  checkpoint_write_CI_Fcoef(pat, dir);

  /* this signifies end of patch info, so this needs to be last */
  fprintf(fp,     "  label = %d\n", CI->label);
}

/* write arrays with coeffs in CI->Fcoef */
void checkpoint_write_CI_Fcoef(tPat *pat, const char *dir)
{
  tCoordInfo *CI = pat->CI;
  char *Fcname = cmalloc(strlen(dir)+128);
  int f;

  for(f=0; f<6; f++)
    if(CI->Fcoef[f])
    {
      /* set output file name in Fcname */
      checkpoint_set_CI_Fcoef_filename(pat, f, dir, Fcname);

      /* write Fcoef into file Fcname */
      //PRF;printf(": writing %s\n", Fcname);
      array_write(NULL, CI->Fcoef[f], Fcname);
    }
  free(Fcname);
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
                "optionally their pt_typ[3]\n\n");
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
  int rank = nMPI_rank();
  /* n_def[d]=0 forces initial write of n[d] as n[d] can never be 0 */
  int n_def[]      = {0,0,0};             /* defaults for n */
  int pt_typ_def[] = {P_LGL,P_LGL,P_LGL}; /* and pt_typ     */
  int rk;
  for(rk=0; rk<size; rk++)
    if(rank==rk || rank==0)
    {
      ulong nelm0s = amr_nelms_on_rank(mesh, rk);
      tElm0 *elm0 = checked_calloc(nelm0s, sizeof(elm0[0]));

      /* send elm0 from rank rk to 0 */
      amr_send_elm0array_of_src_to_dest(mesh, rk, 0, nelm0s, elm0);

      if(rank==0)
      {
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
      }
      free(elm0);
    } /* end rk-loop-if */
  return 0;
}


/******************************************************************/
/* functions to save dat->info */
/******************************************************************/
/* save dat->info */
int checkpoint_save_datinfo(tMesh *mesh, char *fname)
{
  FILE *fp = NULL;
  int IObufsz = Geti(Par("fwrite_bufsize"));
  char *IObuf; /* larger buffer for write */
  int rk;

  /* ranks write one after the other */
  for(rk=0; rk<nMPI_size(); rk++)
  {
    /* do work when it is my turn */
    if(rk == nMPI_rank())
    {
      /* only Rank0 writes the file header */
      if(rk==0)
      {
        /* open destination file */
        fp = fopen_buf(fname, "wb", &IObuf,IObufsz);
        if(!fp) errorexits("failed opening %s", fname);
        fprintf(fp, "number of elms, followed by all elms, and some of their"
                    "dat->info\n\n");
        fprintf(fp, "nelms = %lu\n\n", mesh->eidlim[nMPI_size()-1]);
      }
      else
      {
        fp = fopen_buf(fname, "ab", &IObuf,IObufsz);
        if(!fp) errorexits("failed opening %s", fname);
      }

      /* write all elm and their dat->info */
      errorexit("this is not finished...");
      formyelms(mesh)
      {
        tElm *elm = MyElm;
        tNodeInfo *info = elm->dat->info;
        tEloc eloc[1];
        int l;

        /* get loc strings in eloc */
        eloc_from_eploc(eloc, elm->eploc);

        /* write elm name */
        fprintf(fp, "%d_", eloc->p);
        for(l = 0; l < eloc->l; l++) fputc(eloc->loc[l], fp);
        fprintf(fp, "\n");

        fprintf(fp, "%d ", info->use_fv);
        // ...
        fprintf(fp, "\n");
      }
      fclose_buf(fp, &IObuf);
    }
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

      fclose_buf_file_sync(mesh, fp, &IObuf, Par("file_sync"));
      fs_sync(mesh); /* make sure every MPI proc flushes buffers to disk */
    }
    /* wait until everyone is here */
    MCK( nMPI_barrier() );
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


/******************************************************************/
/* functions to save CRCs */
/******************************************************************/

/* Rank0 saves all CRCs from nmesh_CRCs */
int checkpoint_save_CRCs(tMesh *mesh, char *fname)
{
  ulong CRC[5]; //nmesh CRCs for pars,pats,elms,nbinfo,vars
  FILE *fp;

  /* get current CRCs */
  nmesh_CRCs(mesh, 5, CRC);

  if(Rank0)
  {
    fp = fopen(fname, "w");
    if(fp)
    {
      //fprintf("CRCs for several things in nmesh's memory\n");
      fprintf(fp, "pars:\t%lu\n", CRC[0]);
      fprintf(fp, "pats:\t%lu\n", CRC[1]);
      fprintf(fp, "elms:\t%lu\n", CRC[2]);
      fprintf(fp, "nbinfo:\t%lu\n", CRC[3]);
      fprintf(fp, "vars:\t%lu\n", CRC[4]);
      fclose(fp);
    }
  }
  return 0;
}
