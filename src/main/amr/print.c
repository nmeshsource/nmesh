/* print.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "amr.h"


/* print information about various objects */

void printmesh(tMesh *m)
{
  int p;
  tNode *node;

  printf("mesh=%p: npats=%d npdb=%d nvdb=%d nmyln=%d dt=%g\n",
	 m, m->npats, m->npdb, m->nvdb, m->nmyln, m->dt);
  forpatches(m, p)
    printpatch(m->pat[p]);
  printf("leaf nodes:\n");
  forlnodes(m, node)
  {
    printnode(node);
  } endforlnodes;
}

void printpatch(tPat *pat)
{
  printf("p%d: [%g,%g]x[%g,%g]x[%g,%g] nmax=%d\n",
         pat->p, pat->bbox[0], pat->bbox[1], pat->bbox[2], pat->bbox[3],
         pat->bbox[4],pat->bbox[5], pat->nmax);
  printf("root node:\n");
  printnode(pat->rnode);
/*
  printf("leaf nodes again:\n");
  fornodelist(pat->lns, el)
  {
    printf("%p: prev=%p next=%p\n", el, el->prev, el->next);
    printnode(el->node);
  }
*/
}

void printnode(tNode *n) 
{
  int i, j;

  printf("ijk%d: nid%ld p%d [%g,%g]x[%g,%g]x[%g,%g] np=%dx%dx%d=%d\n",
         n->ijk, n->nid, n->pat->p, n->bbox[0], n->bbox[1], n->bbox[2],
         n->bbox[3], n->bbox[4],n->bbox[5], n->n[0], n->n[1], n->n[2], n->np);
  printf(" l%d leaf=%d: ", n->l, n->leaf);
  printf(" datrank=%d  dat: %s ", n->datrank, n->dat ? "yes" : "no");
  printf(" patface=");
  for(i=0; i<6; i++) printf("%d", n->patface[i]);
  printf("\n");
  printf(" nb =");
  for(i=0; i<6; i++) printf(" %ld", get_node_nid(n->nb[i]));
  printf("   parent->nid=%ld\n", get_node_nid(n->parent));
  for(i=0; i<6; i++)
  {
    printf(" |fnb[%d]:", i);
    for(j=0; j<n->nfnb[i]; j++) printf(" %ld", get_node_nid(n->fnb[i][j]));
    printf("|");
  }
  printf("\n");
  printf(" child =");
  for(i=0; i<8; i++) printf(" %ld", get_node_nid(n->child[i]));
  //printf("\n");
  //printf(" Dt =");
  //for(i=0; i<3; i++) printf(" %p", n->Dt[i]);
  //printf("\n");

  printf("\n");
}

void printnode_and_neighbors(tNode *n)
{
  tNode *n0, *n0p, *n1;
  int dir, ni;

  printnode(n);
  for(dir=0; dir<3; dir++)
  {
    printf("dir%d l%d neighbors cover:\n", dir, n->l);
    ni=0;
    for(n0=n;   n0; n0=n0->nb[dir*2], ni--) n0p = n0;
    ni++;
    for(n1=n0p; n1; n1=n1->nb[dir*2+1], ni++)
    {
      printf(" [%+5g,%+5g]x[%+5g,%+5g]x[%+5g,%+5g]  ni=%d, nid%ld\n",
             n1->bbox[0], n1->bbox[1], n1->bbox[2], n1->bbox[3],
             n1->bbox[4], n1->bbox[5], ni, n1->nid);
    }
  }
}

void printnodelist_and_neighbors(tNlist *nl)
{
  tNlist *el;
  fornodelist(first_nodelist(nl), el)
  {
    printf("nid%ld: ", get_node_nid(el->node));
    if(el->prev) printf("prev=%ld ", get_node_nid(el->prev->node));
    if(el->next) printf("next=%ld\n", get_node_nid(el->next->node));
    else         printf("\n");
    printnode_and_neighbors(el->node);
  }
}

void printnodelist(tNlist *nl)
{
  tNlist *el;
  fornodelist(first_nodelist(nl), el)
  {
    printf("nid%ld: ", get_node_nid(el->node));
    if(el->prev) printf("prev=%ld ", get_node_nid(el->prev->node));
    if(el->next) printf("next=%ld\n", get_node_nid(el->next->node));
    else         printf("\n");
    //printf("%p: prev=%p next=%p\n", el, el->prev, el->next);
    //printnode(el->node);
  }
  if(!nl) printf("<empty nodelist>\n");
} 


/* print a variable in a node */
void printvar_innode(tNode *node, char *name)
{
  tMesh *mesh = node->pat->mesh;
  int vi = Ind(name);
  tArray *va = node->dat->v[vi];

  printf("%s, Ind=%d:\n", name, vi);
  printarray(va);
}

/* print an array */
void printarray(tArray *A)
{
  int i,j,k;
  printf("->n[] = {%d,%d,%d}\n", A->n[0],A->n[1],A->n[2]);
  for(k=0; k<A->n[2]; k++)
  {
    for(j=0; j<A->n[1]; j++)
    {
      for(i=0; i<A->n[0]; i++)
        printf(" %g", A->a[Ind_n(i,j,k, A->n)]);
      printf("\n");
    }
    printf("\n");
  }
}

/* print an array */
void printarray_matrix0(tArray *A)
{
  int i,J;
  int nJ = A->n[1] * A->n[2];
  int dJ = A->n[0];
  printf("->n[] = {%d,%d,%d} => %dx%d\n", A->n[0],A->n[1],A->n[2], A->n[0],nJ);
  for(i=0; i<A->n[0]; i++)
  {
    for(J=0; J<nJ; J++)
      printf(" %g", A->a[i + dJ*J]);
    printf("\n");
  }
}

/* print an array */
void printarray_matrix1(tArray *A)
{
  int j,J;
  int nJ = A->n[0] * A->n[2];
  int dJ = A->n[1];
  printf("->n[] = {%d,%d,%d} => %dx%d\n", A->n[0],A->n[1],A->n[2], A->n[1],nJ);
  for(j=0; j<A->n[1]; j++)
  {
    for(J=0; J<nJ; J++)
      printf(" %g", A->a[j + dJ*J]);
    printf("\n");
  }
}

/* print an array */
void printarray_matrix2(tArray *A)
{
  int k,J;
  int nJ = A->n[1] * A->n[0];
  int dJ = A->n[2];
  printf("->n[] = {%d,%d,%d} => %dx%d\n", A->n[0],A->n[1],A->n[2], A->n[2],nJ);
  for(k=0; k<A->n[2]; k++)
  {
    for(J=0; J<nJ; J++)
      printf(" %g", A->a[k + dJ*J]);
    printf("\n");
  }
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
