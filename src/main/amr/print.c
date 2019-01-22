/* print.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "amr.h"


/* print information about various objects */

void printmesh(tMesh *m)
{
  int i;

  printf("mesh=%p: npats=%d, nvdb=%d, dt=%g\n",
	 m, m->npats, m->nvdb, m->dt);
  //forallpats
  // printpatch(pat);
}

void printpatch(tPat *pat)
{
  printf("p=%d:\n", pat->p);
}

void printnode(tNode *n) 
{
  printf("l=%d", n->l);
  printf("\n");
}


/* print a variable in a node */
void printvar_innode(tNode *node, char *name)
{
  tMesh *mesh = node->pat->mesh;
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
  tMesh *mesh = vl->mesh;
  int i;

  for(i = 0; i < vl->n; i++)
    printf("%s\n", VarName(vl->index[i]));      
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
