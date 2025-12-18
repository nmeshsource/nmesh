/* output_mesh.c */
/* Wolfgang Tichy, May 2023 */

#include "nmesh.h"
#include "output.h"


/*************************************************************************/
/* funcs to write elms, i.e. leaf nodes */
/*************************************************************************/

/* write essential info */
void write_elm0(FILE *fp, tElm *e, int mode, const char *s)
{
  char nns[100];
  fprintf(fp, "%s", elmname(e, nns,99));
  fprintf(fp, ": eid%lu [%g,%g]x[%g,%g]x[%g,%g]",
         Elm_eid(e),
         e->bbox[0],e->bbox[1], e->bbox[2],e->bbox[3], e->bbox[4],e->bbox[5]);
  fprintf(fp, " n=%dx%dx%d=%d", e->n[0],e->n[1],e->n[2], e->np);
  switch(mode)
  {
  case 1:
    fprintf(fp, " datrank=%d", e->datrank);
    break;
  //default:
  }
  fprintf(fp, "%s", s);
}


/* write fnb on one face */
void write_elm_fnb(FILE *fp, tElm *e, int face, int sorted)
{
  int j;
  int nlocs = 100;
  int nnb = e->nfnb[face];
  char **fnbstr;

  if(nnb)
    fnbstr = rows_calloc_matrix(nnb,nlocs, sizeof(fnbstr[0][0]));
  else
    fnbstr = NULL;

  for(j=0; j<nnb; j++)
  {
    if(e->fnb[face][j]) elmname(e->fnb[face][j], fnbstr[j],nlocs);
    else                sprintf(fnbstr[j], "nil");
  }
  //for(j=0; j<nnb; j++) printf("%s\n",fnbstr[j]);

  if( sorted && (nnb>1) )
    qsort(fnbstr, nnb, sizeof(fnbstr[0]), qsort_compar_strlist);

  fprintf(fp, " {");
  for(j=0; j<nnb; j++)
  {
    fprintf(fp, "%s", fnbstr[j]);
    if(j<nnb-1) fprintf(fp, " ");
  }
  fprintf(fp, "}");
  //rows_free(fnbstr, nnb);
  for(j=0; j<nnb; j++) free(fnbstr[j]);
  free(fnbstr);
}


/* write essential info + some more */
void write_elm(FILE *fp, tElm *e, int mode)
{
  tDat *dat = e->dat;
  int i;
  //union { const tElm *elm; tElm0 *elm0; } e2e0;
  //e2e0.elm = e;
  //printelm0(e2e0.elm0, "");
  write_elm0(fp, e, mode, "");

  switch(mode)
  {
  case 1:
    fprintf(fp, " dat:%s\n", dat ? "yes" : "no");
    break;
  default:
    fprintf(fp, "\n");
  }

  fprintf(fp, " fnb =");
  for(i=0; i<6; i++)
    write_elm_fnb(fp, e, i, 1);
  fprintf(fp, "\n");
}

/* open file and write all my node elms into it */
void write_mylnodes(tMesh *mesh, const char *info, int mode)
{
  int size = nMPI_size();
  int rank = nMPI_rank();
  int rk;
  int outd = Par("outdir");
  char *outdir = Gets(outd);
  char fname[1000];
  FILE *fp;

  /* MPI motivated loop to assign work */
  for(rk=0; rk<size; rk++)
  {
    /* do work when it is my turn */
    if(rk == rank)
    {
      sprintf(fname, "%s/%s", outdir, "myelms.txt");
      fp = fopen(fname, "a");

      if(rk==0) fprintf(fp, "%s\n", info);

      formylnodes_noomp(mesh)
      {
        tNode *elm = MyLnode;
        write_elm(fp, elm, mode);
      }
      if(rk==size-1) fprintf(fp, "\n");
      fclose(fp);
    }
    /* wait until everyone is here */
    MCK( nMPI_barrier() );
  } /* end rk-loop */
}

/* open file and write neighbor node elms into it */
void write_nblnodes(tMesh *mesh, const char *info, int mode)
{
  int rk;
  int outd = Par("outdir");
  char *outdir = Gets(outd);
  char fname[1000];
  FILE *fp;

  /* MPI motivated loop to assign work */
  for(rk=0; rk<nMPI_size(); rk++)
  {
    /* do work when it is my turn */
    if(rk == nMPI_rank())
    {
      ulong ei;

      sprintf(fname, "%s/%s.%d", outdir, "nbelms", rk);
      fp = fopen(fname, "a");

      fprintf(fp, "%s\n", info);

      for(ei=0; ei < mesh->nnbelm; ei++)
      {
        tNode *elm = mesh->nbelm[ei];
        write_elm(fp, elm, mode);
      }
      fprintf(fp, "\n");
      fclose(fp);
    }
    /* wait until everyone is here */
    MCK( nMPI_barrier() );
  } /* end rk-loop */
}

/*************************************************************************/
/* funcs to write all dat->info */
/*************************************************************************/

/* write a single dat->info from inside a elm->dat */
void write_elm_dat_info(FILE *fp, tElm *elm, int mode)
{
  tDat *dat = elm->dat;
  if(dat)
  {
    tNodeInfo *info = dat->info;
    int i;
    fprintf(fp, " evo_troubled=%d trbl_score=%d use_fv=%d nlim=%d",
            info->evo_troubled, info->trbl_score, info->use_fv, info->nlim);
    if(mode==1)
    {
      fprintf(fp, " trbl_ref=");
      fwrite_little(info->trbl_ref, sizeof(info->trbl_ref[0]),1 , fp);
      fprintf(fp, "\n");
      fprintf(fp, " load_timer_running=%d load_TimeIn_s=%g load_start=",
              info->load_timer_running, info->load_TimeIn_s);
      fwrite_little(info->load_start, sizeof(info->load_start[0]),1 , fp);
      fprintf(fp, "\n");
      fprintf(fp, " desrank=%d ", info->desrank);
    }
    else
    {
      fprintf(fp, "\n");
    }
    fprintf(fp, " nnbinfo=");
    for(i=0; i<6; i++) fprintf(fp, "%d ", info->nnbinfo[i]);
    fprintf(fp, " unlimited=%d\n", info->unlimited);
  }
}

/* open file and write the dat->info from all elms into it */
void write_elm_dat_infos(tMesh *mesh, const char *header, int mode)
{
  int size = nMPI_size();
  int rank = nMPI_rank();
  int rk;
  int outd = Par("outdir");
  char *outdir = Gets(outd);
  char fname[1000];
  FILE *fp;

  /* MPI motivated loop to assign work */
  for(rk=0; rk<size; rk++)
  {
    /* do work when it is my turn */
    if(rk == rank)
    {
      sprintf(fname, "%s/%s", outdir, "elm_dat_infos.txt");
      fp = fopen(fname, "a");

      if(rk==0 && header) fprintf(fp, "%s\n", header);

      formyelms_noomp(mesh)
      {
        tNode *elm = MyElm;
        write_elm0(fp, elm, mode, ":\n");
        write_elm_dat_info(fp, elm, mode);
      }
      if(rk==size-1) fprintf(fp, "\n");
      fclose(fp);
    }
    /* wait until everyone is here */
    MCK( nMPI_barrier() );
  } /* end rk-loop */
}

/*************************************************************************/
/* funcs to write all amr_elm_nbinfo? vars */
/*************************************************************************/

/* write the amr_elm_nbinfo? for one elm */
void write_amr_elm_nbinfo(FILE *fp, tElm *elm, int elm_nbinfo0)
{
  tDat *dat = elm->dat;
  if(dat)
  {
    int f;
    for(f=0; f<6; f++)
    {
      int i_nbinfo = elm_nbinfo0 + f;
      tArray *nbinfo = VarA(elm, i_nbinfo);
      int Neplocs, ni;

      /* if there is no nbinfo do nothing */
      if(!nbinfo) continue;
      Neplocs = array_Neplocs(nbinfo);
      if(!Neplocs) continue;

      fprintf(fp, " f%d: %d: ", f, Neplocs);

      /* loop over eplocs in nbinfo and write them */
      for(ni=0; ni<Neplocs; ni++)
      {
        tEploc *eploc = &(nbinfo->eploc[ni]);
        ulong eid = eploc->eid;
        char nbname[100];
        eplocname(eploc, nbname,100);
        fprintf(fp, " %s:eid%lu", nbname, eid);
      }
      fprintf(fp, "\n");
    }
  }
}

/* open file and write amr_elm_nbinfo? from all elms into it */
void write_amr_elm_nbinfos(tMesh *mesh, const char *header, int mode)
{
  int elm_nbinfo0 = Ind("amr_elm_nbinfo0");
  int size = nMPI_size();
  int rank = nMPI_rank();
  int rk;
  int outd = Par("outdir");
  char *outdir = Gets(outd);
  char fname[1000];
  FILE *fp;

  /* MPI motivated loop to assign work */
  for(rk=0; rk<size; rk++)
  {
    /* do work when it is my turn */
    if(rk == rank)
    {
      sprintf(fname, "%s/%s", outdir, "amr_elm_nbinfos.txt");
      fp = fopen(fname, "a");

      if(rk==0 && header) fprintf(fp, "%s\n", header);

      formyelms_noomp(mesh)
      {
        tNode *elm = MyElm;
        write_elm0(fp, elm, mode, ":\n");
        write_amr_elm_nbinfo(fp, elm, elm_nbinfo0);
      }
      if(rk==size-1) fprintf(fp, "\n");
      fclose(fp);
    }
    /* wait until everyone is here */
    MCK( nMPI_barrier() );
  } /* end rk-loop */
}
