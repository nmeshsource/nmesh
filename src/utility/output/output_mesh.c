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


/* write fnb on one face */
void write_lnode_fnb(FILE *fp, tNode *e, int face, int sorted)
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
    if(e->fnb[face][j]) nodename(e->fnb[face][j], fnbstr[j],nlocs);
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
void write_lnode(FILE *fp, tNode *e, int mode)
{
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
    write_lnode_fnb(fp, e, i, 1);
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
        write_lnode(fp, elm, mode);
      }
      if(rk==size-1) fprintf(fp, "\n");
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
