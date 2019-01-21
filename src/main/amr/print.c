/* print.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "amr.h"


/* print information about various objects */

void printmesh(tMesh *g)
{
  int i;

  printf("mesh=%p: nboxes=%d, nvariables=%d, dt=%g\n",
	 g, g->nboxes, g->nvariables, g->dt);
  for (i = 0; i < g->nboxes; i++)
    printbox(g->box[i]);
}

void printpatch(tPat *box)
{
  printf("p=%d:\n", box->b);
}

void printnode(tNode *n) 
{
  printf("l=%d", n->l);
  printf("\n");
}

/*
void printbface(tBface *bface)
{
  int np;
  printf("b%d f%d fi%d  ob%d ofi%d  oXi=%d oYi=%d oZi=%d  ",
          bface->b,bface->f,bface->fi, bface->ob,bface->ofi,
          bface->oXi,bface->oYi,bface->oZi);
  printf("bits=%u%u%u%u%u%u,%u,%u%u",
         bface->touch,
         bface->sameX, bface->sameY, bface->sameZ,  bface->same_fpts,
         bface->fpts_off_face,
         bface->setnormalderiv, bface->innerbound,  bface->outerbound);
  if(bface->fpts==NULL) np = 0;
  else                  np = bface->fpts->npoints[bface->b];
  printf("  (%d points)\n", np);
  //prPointList(bface->fpts);
  //prPointList_ijk_inbox(bface->fpts, bface->b);
}

void printbfaces(tPat *pat)
{
  int i;
  printf("pat->b=%d\n", pat->b);
  for(i=0; i<pat->nbfaces; i++)
  {
    tBface *bface = pat->bface[i];
    printf("bface[%d]: ", i);
    printbface(bface);
  }
}
*/

/* print a variable in a node */
void printvar_innode(tNode *node, char *name)
{
  double *v = 0; //node->v[Ind(name)];
  int i,j,k;
  int *n = node->n;

  printf("%s, %p, Ind=%d:\n", name, v, Ind(name));
  if(v)
  {
    //forallijk(i,j,k)
    {
      printf("%g ", v[Index(i,j,k)]);
    }
  }    
  printf("%s ends here.\n", name);
}


/* print a variable list */
void printVarList(tVarList *vl)
{
  int i;
  tMesh *mesh = vl->mesh;

  for(i = 0; i < vl->n; i++)
    printf("%s\n", VarName(vl->index[i]));      
}

