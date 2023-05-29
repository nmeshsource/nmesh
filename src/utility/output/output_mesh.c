/* output_mesh.c */
/* Wolfgang Tichy, May 2023 */

#include "nmesh.h"
#include "output.h"



/* write essential info */
void write_lnode0(FILE *fp, tNode *e, const char *s)
{
  char nns[100];
  fprintf(fp, "%s", nodename(e, nns,99));
  fprintf(fp, ": eid%lu [%g,%g]x[%g,%g]x[%g,%g]",
         Node_eid(e),
         e->bbox[0],e->bbox[1], e->bbox[2],e->bbox[3], e->bbox[4],e->bbox[5]);
  fprintf(fp, " n=%dx%dx%d=%d", e->n[0],e->n[1],e->n[2], e->np);
  fprintf(fp, " datrank=%d", e->datrank);
  fprintf(fp, "%s", s);
}

/* write essential info + some more */
void write_lnode(FILE *fp, tNode *e, int mode)
{
  char nns[100];
  tDat *dat = e->dat;
  int i;
  //union { const tNode *elm; tNode0 *elm0; } e2e0;
  //e2e0.elm = e;
  //printelm0(e2e0.elm0, "");
  write_lnode0(fp, e, "");

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
  {
    int j;
    fprintf(fp, " {");
    {
      for(j=0; j<e->nfnb[i]; j++)
      {
        //if(e->fnb[i][j]) printeploc(e->fnb[i][j]->eploc);
        if(e->fnb[i][j]) fprintf(fp, "%s", nodename(e->fnb[i][j], nns,99));
        else             fprintf(fp, "nil");
        if(j<e->nfnb[i]-1) fprintf(fp, " ");
      }
    }
    fprintf(fp, "}");
  }
  fprintf(fp, "\n");
}

/* open file and write all my node elms into it */
void write_mylnodes(tMesh *mesh, const char *info, int mode)
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
      sprintf(fname, "%s/%s.%d", outdir, "myelms", rk);
      fp = fopen(fname, "a");

      fprintf(fp, "%s\n", info);

      formylnodes_noomp(mesh)
      {
        tNode *elm = MyLnode;
        write_lnode(fp, elm, mode);
      }
      fprintf(fp, "\n");
      fclose(fp);
    }
    /* wait until everyone is here */
    nMPI_barrier();
  } /* end rk-loop */
}

/* open file and write neighbor node elms into it */
/*
void write_nblnodes(tMesh *mesh)
{
  int ei;
  for(ei=0; ei < mesh->nnbelm; ei++)
  {
    tNode *elm = mesh->nbelm[ei];
    write_lnode(fp, elm);
  }
}
*/
