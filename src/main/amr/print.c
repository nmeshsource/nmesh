/* print.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "amr.h"


/* print information about various objects */

void printmesh(tMesh *m)
{
  int p;

  printf("mesh=%p: npats=%d npdb=%d nvdb=%d dt=%g\n",
	 m, m->npats, m->npdb, m->nvdb, m->dt);
  forpatches(m, p)
    printpatch(m->pat[p]);
}

void printpatch(tPat *pat)
{
  tNlist *el;
  tNode *node;

  printf("p%d: [%g,%g]x[%g,%g]x[%g,%g] nD=%d\n",
         pat->p, pat->bbox[0], pat->bbox[1], pat->bbox[2], pat->bbox[3],
         pat->bbox[4],pat->bbox[5], pat->nD);
  printf("root node:\n");
  printnode(pat->rnode);
  printf("leaf nodes:\n");
  forlnodes(pat, node)
  {
    printnode(node);
  } endforlnodes;


  printf("leaf nodes again:\n");
  fornodelist(pat->lns, el)
  {
    printf("%p: prev=%p next=%p\n", el, el->prev, el->next);
    printnode(el->node);
  }


}

void printnode(tNode *n) 
{
  int i;

  printf("ijk%d: p%d [%g,%g]x[%g,%g]x[%g,%g] np=%dx%dx%d=%d",
         n->ijk, n->pat->p, n->bbox[0], n->bbox[1], n->bbox[2], n->bbox[3],
         n->bbox[4],n->bbox[5], n->n[0], n->n[1], n->n[2], n->np);
  printf(" leaf=%d l=%d", n->leaf, n->l);
  printf(" patface=");
  for(i=0; i<6; i++) printf("%d", n->patface[i]);
  printf("\n");
  printf(" node = %p   parent = %p\n", n, n->parent);
  printf(" nb =");
  for(i=0; i<6; i++) printf(" %p", n->nb[i]);
  printf("\n");
  printf(" child =");
  for(i=0; i<8; i++) printf(" %p", n->child[i]);
  printf("\n");
  printf(" D =");
  for(i=0; i<3; i++) printf(" %p", n->D[i]);
  printf("\n");
  printf(" datrank=%d  dat=%p\n", n->datrank, n->dat);

  printf("\n");
}

void printnodelist(tNlist *nl)
{
  tNlist *el;
  fornodelist(nl, el)
  {
    printf("%p: prev=%p next=%p\n", el, el->prev, el->next);
    printnode(el->node);
  }
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
