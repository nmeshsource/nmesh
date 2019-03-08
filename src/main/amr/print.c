/* print.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "amr.h"

/* mesh made in main */
extern tMesh *main_mesh;


/* print information about various objects */

void printmesh(tMesh *m)
{
  int p;
  tNode *node;

  if(m==main_mesh) printf("mesh=main_mesh: ");
  else             printf("mesh=%p: ", m);
  printf("npats=%d npdb=%d nvdb=%d nln=%ld myln->nm=%d dt=%g\n",
	 m->npats, m->npdb, m->nvdb, m->nln, m->myln->nm, m->dt);
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

void printCI(tPat *pat)
{
  tCoordInfo *CI = pat->CI;
  int i, j, p;

  printf("p%d: CI->dom=%d CI->type=%d\n",
         pat->p, CI->dom, CI->type);

  /* s, xc */
  printf("CI->s[] = [ ");
  for(i=0;i<6;i++) printf("%g ", CI->s[i]);
  printf("]\n");
  printf("CI->xc[] = [ ");
  for(i=0;i<4;i++) printf("%g ", CI->xc[i]);
  printf("]\n");

  /* iFS, iSurf */
  for(i=p=0; i<6; i++) if(CI->iFS[i]) p++;
  if(p)
  {
    printf("CI->iFS[] = [ ");
    for(i=0;i<6;i++) printf("%d ", CI->iFS[i]);
    printf("]\n");
    printf("CI->iSurf[] = [ ");
    for(i=0;i<6;i++) printf("%d ", CI->iSurf[i]);
    printf("]\n");
    printf("idSurfdX[][] = [");
    for(i=0;i<6;i++)
    {
      for(j=1;j<4;j++) printf(" %d", CI->idSurfdX[i][j]);
      printf(";");
    }
    printf("]\n");
  }

  /* FSurf */
  printf("CI->FSurf[] = [ ");
  for(i=0;i<6;i++) printf("%p ", CI->FSurf[i]);
  printf("]\n");
  for(i=p=0; i<6; i++) if(CI->dFSurfdC[i]) p++;
  if(p)
  {
    printf("dFSurfdC[] = [");
    for(i=0;i<6;i++) printf(" %p", CI->dFSurfdC[i]);
    printf("]\n");
  }
}

void printnode(tNode *n) 
{
  int i, j;
  char s[100];

  printf("nid%ld: ijk%d %s p%d [%g,%g]x[%g,%g]x[%g,%g] np=%dx%dx%d=%d\n",
          n->nid, n->ijk, node_location_str(n, s,100),
          n->pat->p, n->bbox[0], n->bbox[1], n->bbox[2],
         n->bbox[3], n->bbox[4],n->bbox[5], n->n[0], n->n[1], n->n[2], n->np);
  printf(" l%d leaf=%d: ", n->l, n->leaf);
  printf(" datrank=%d  dat: %s ", n->datrank, n->dat ? "yes" : "no");
  printf(" patface=");
  for(i=0; i<6; i++) printf("%d", n->patface[i]);
  printf("\n");
  printf(" nb =");
  for(i=0; i<6; i++) printf(" %ld", get_node_nid(n->nb[i]));
  printf("   parent->nid=%ld\n", get_node_nid(n->parent));
  printf(" fnb =");
  for(i=0; i<6; i++)
  {
    //printf(" %d:{", i);
    printf(" {");
    for(j=0; j<n->nfnb[i]; j++) printf(" %ld", get_node_nid(n->fnb[i][j]));
    printf(" }");
  }
  printf("\n");
  if(!n->leaf)
  {
    printf(" child =");
    for(i=0; i<8; i++) printf(" %ld", get_node_nid(n->child[i]));
    printf("\n");
  }
  //printf("\n");
  //printf(" Dt =");
  //for(i=0; i<3; i++) printf(" %p", n->Dt[i]);
  //printf("\n");
}

void printnode_and_neighbors(tNode *n)
{
  tNode *n0, *n0p, *n1;
  int dir, ni;
  char s[100];

  printnode(n);
  for(dir=0; dir<3; dir++)
  {
    printf("dir%d l%d neighbors cover:\n", dir, n->l);
    ni=0;
    for(n0=n;   n0; n0=n0->nb[dir*2], ni--) n0p = n0;
    ni++;
    for(n1=n0p; n1; n1=n1->nb[dir*2+1], ni++)
    {
      printf(" [%+5g,%+5g]x[%+5g,%+5g]x[%+5g,%+5g]  ni=%d, nid%ld %s\n",
             n1->bbox[0], n1->bbox[1], n1->bbox[2], n1->bbox[3],
             n1->bbox[4], n1->bbox[5], ni, n1->nid,
             node_location_str(n1, s,99));
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
void printvar_innode(tNode *node, int vi)
{
  int f;
  tMesh *mesh = node->pat->mesh;
  char *name = VarName(vi);
  int zones =  MeshVarSurfacezones(mesh, vi);
  int type =  MeshVarType(mesh, vi);
  tDat *dat = node->dat;
  tArray *va = dat ? dat->v[vi] : NULL;

  printf("%s Ind=%d type=%d zones=%d Array", name, vi, type, zones);
  printarray(va);

  if(va)
    for(f=0; f<6; f++)
    {
      tSurface *sf = dat->s[f][vi];
      tArray *msa  = sf ? sf->mysurf : NULL;
      tArray *asa  = sf ? sf->ajsurf : NULL;
      tArray **nsa = sf ? sf->nbsurf : NULL;
      if(msa)
      {
        printf("f%d allocd=%d mysurf", f, sf->allocd_mysurf);
        printarray(msa);
      }
      if(asa)
      {
        printf("f%d allocd=%d ajsurf", f, sf->allocd_ajsurf);
        printarray(asa);
      }
      if(nsa)
      {
        int ni;
        for(ni=0; ni<sf->nnbsurf; ni++)
        {
          printf("f%d allocd=%d nbsurf[%d]", f, sf->allocd_nbsurf[ni], ni);
          printarray(nsa[ni]);
        }
      }
    }
}

void printvar_ajsurfdiff(tNode *node, int vi)
{
  int f;
  tMesh *mesh = node->pat->mesh;
  char *name = VarName(vi);
  int zones =  MeshVarSurfacezones(mesh, vi);
  int type =  MeshVarType(mesh, vi);
  tDat *dat = node->dat;
  tArray *va = dat ? dat->v[vi] : NULL;

  printf("%s Ind=%d type=%d zones=%d\n", name, vi, type, zones);

  if(va)
    for(f=0; f<6; f++)
    {
      tSurface *sf = dat->s[f][vi];
      tArray *msa  = sf ? sf->mysurf : NULL;
      tArray *asa  = sf ? sf->ajsurf : NULL;
      tArray **nsa = sf ? sf->nbsurf : NULL;
      if(msa)
      {
        printf("f%d allocd=%d mysurf\n", f, sf->allocd_mysurf);
      }
      if(nsa)
      {
        int ni;
        for(ni=0; ni<sf->nnbsurf; ni++)
        {
          printf("f%d allocd=%d nbsurf[%d]\n", f, sf->allocd_nbsurf[ni], ni);
        }
      }
      if(asa)
      {
        tArray *diff = alloc_array(asa->n);
        int k;
        printf("f%d allocd=%d ajsurf\n", f, sf->allocd_ajsurf);
        forarray(diff, k) diff->d[k] = asa->d[k] - msa->d[k];
        printarray(diff);
        free_array(diff);
      }
    }
}

/* print an array */
void printarray_sel(tArray *A, int dbl)
{
  int i,j,k;
  if(!A)
  {
    printf(" <NULL array>\n\n");
    return;
  }
  printf("->n[] = {%d,%d,%d}\n", A->n[0],A->n[1],A->n[2]);
  for(k=0; k<A->n[2]; k++)
  {
    for(j=0; j<A->n[1]; j++)
    {
      for(i=0; i<A->n[0]; i++)
      {
        if(dbl) printf(" %g", A->d[Ind_n(i,j,k, A->n)]);
        else    printf(" %d", A->i[Ind_n(i,j,k, A->n)]);
      }
      printf("\n");
    }
    printf("\n");
  }
}
void printarray(tArray *A)
{
  printarray_sel(A, 1);
}
void printarray_int(tArray *A)
{
  printarray_sel(A, 0);
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
      printf(" %g", A->d[i + dJ*J]);
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
      printf(" %g", A->d[j + dJ*J]);
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
      printf(" %g", A->d[k + dJ*J]);
    printf("\n");
  }
}

/* print one bface */
void printthisbface(tBface *bface, char *s)
{
  if(!bface) return;
  printf(" %s: p%d f%d op%d", s, bface->pat->p,bface->f, bface->op);
  if(bface->brct_isset)
    printf(" [%g,%g]x[%g,%g]",
           bface->brct[0], bface->brct[1], bface->brct[2], bface->brct[3]);
  else
    printf(" [ , ]x[ , ]");
  printf(" ioX={%d,%d,%d} bits=%1d%1d%1d\n",
         bface->ioX[0],bface->ioX[1],bface->ioX[2],
         bface->face2, bface->innerbound,  bface->outerbound);
}

/* print one bface and its pair */
void printbface(tBface *bface)
{
  printthisbface(bface, "A");
  printthisbface(bface->obface, "B");
}

/* print all patch bfaces */
void printbfaces(tPat *pat)
{
  tBface *bf;

  printf("pat->p=%d\n", pat->p);
  forbfaces(pat, bf)
    printbface(bf);
}

/* print all bfaces in mesh */
void printallbfaces(tMesh *mesh)
{
  int p;
  forpatches(mesh, p)
    printbfaces(mesh->pat[p]);
}

/* print 3-vec */
void pr3v(char *s, double x[3])
{
  printf("%s=%g %g %g  ", s, x[0],x[1],x[2]);
}

void prbbox(double *bb, int dim)
{
  int d;
  printf("[%g,%g]", bb[0], bb[1]);
  for(d=1; d<dim; d++)
    printf("x[%g,%g]", bb[2*d], bb[2*d+1]);
  printf(" ");
}
