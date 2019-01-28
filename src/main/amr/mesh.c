/* mesh.c */
/* Wolfgang Tichy, 1/2019 */


#include "nmesh.h"
#include "amr.h"

#define PR 0




/* make an empty mesh, into which we an then initialize or into which
   we can e.g. copy the contents of an existing mesh */
tMesh *make_empty_mesh(int pr)
{
  tMesh *mesh;

  /* print info */
  if(pr) prdivider(0);

  mesh = alloc_mesh(0);

  if(pr) { PRFs(":\n"); printmesh(mesh); }

  /* return pointer to newly created mesh */
  return mesh;
}

/* add apatch to the mesh */
int add_patch(tMesh *mesh, double bbox[6], int nroot[3], int nD)
{
  tPat *pat;
  int p = mesh->npats;
  int i;

  /* make room for new patch in mesh and then add an empty patch */
  realloc_patlist_in_mesh(mesh, p + 1);
  pat = alloc_patch(mesh, p, nD);
  mesh->pat[p] = pat;

  /* set bbox */
  for(i=0; i<6; i++) pat->bbox[i] = bbox[i];

  /* set diff matrices */
  pat->D[5][1]->a[4] = 4;

  /* setup root node */
  pat->rnode = make_root_node(pat, nroot, 0);
  pat->lns = alloc_nodelist(pat->rnode);

  return 0;
}



/* a function just for testing */
int setup_test_mesh(tMesh *mesh)
{
  double bbox[6] = { -4,4, -2,2, -1,1 };
  int n1max = 55;
  int n[3] = { 9,9,9 };
  tNlist *nlist, *el;
  tNode *nd;
  int i;

  PRFs(":\n");

//tNode *tnode = alloc_node();
//mesh->pat[0]->rnode = 0;
  //realloc_patlist_in_mesh(mesh, 1);
  add_patch(mesh, bbox, n, n1max);

  nd = mesh->pat[0]->rnode;
  nlist = make8_child_nodes(nd, n);

  mesh->pat[0]->lns = first_replace1_in_nodelist(mesh->pat[0]->lns, nlist);
  //printnodelist(nlist);
  printmesh(mesh);

  el = mesh->pat[0]->lns;
  for(i=1; i<=1; i++) el = el->next;
  nd = el->node;
  nlist = make8_child_nodes(nd, n);
  mesh->pat[0]->lns = first_replace1_in_nodelist(el, nlist);

  el = mesh->pat[0]->lns;
  for(i=1; i<=8+2; i++) el = el->next;
  nd = el->node;
  nlist = make8_child_nodes(nd, n);
  mesh->pat[0]->lns = first_replace1_in_nodelist(el, nlist);

  //printnodelist(nlist);
  printmesh(mesh);

  printnodelist_and_neighbors(mesh->pat[0]->lns);


double  A[6] = { 1,2,
                 3,4,
                 5,6 };
/*
double B0[24] = { 1,2,3,4,5,6,7,8,9,10,11,12,
                  13,14,15,16,17,18,19,20,21,22,23,24 };
*/
double B0[24], B1[24], B2[24];
double AB[36];
int j,k;
int nB0[] = {2,3,4};
int nB1[] = {3,2,4};
int nB2[] = {4,3,2};

for(i=0; i<12; i++) { B0[2*i] = 2*i+1; B0[2*i+1] = 2*i+2; }

for(k=0; k<4; k++)
for(j=0; j<3; j++)
for(i=0; i<2; i++)
B2[Ind_n(k,j,i, nB2)] = B1[Ind_n(j,i,k, nB1)] = B0[Ind_n(i,j,k, nB0)];

n[0]=2; n[1]=2; n[2]=1;
tArray *Aa = alloc_array(n);
Aa->a = A;
//Aa->n[0]=3; Aa->n[1]=2; Aa->n[1]=1;

tArray *B0a = alloc_array(nB0);
B0a->a = B0;
tArray *B1a = alloc_array(nB1);
B1a->a = B1;
tArray *B2a = alloc_array(nB2);
B2a->a = B2;


int nC0[] = {2,3,4};
tArray *C0a = alloc_array(nC0);
double *C0 = C0a->a;
//printarray_matrix0(Aa);
printarray(Aa);
//printarray_matrix0(B0a);
printarray(B0a);
mm_array0(Aa,B0a, C0a);
//printarray_matrix0(C0a);
printarray(C0a);


Yo(1);
int nC1[] = {3,2,4};
tArray *C1a = alloc_array(nC1);
double *C1 = C1a->a;
//printarray_matrix1(B1a);
printarray(B1a);
mm_array1(Aa,B1a, C1a);
//printarray_matrix0(Ca1);
set_const_array(C0a, 0.);
for(k=0; k<4; k++)
for(j=0; j<2; j++)
for(i=0; i<3; i++)
C0[Ind_n(j,i,k, nC0)] = C1[Ind_n(i,j,k, nC1)];
printarray_matrix0(C0a);
printarray(C1a);

Yo(2);
int nC2[] = {4,3,2};
tArray *C2a = alloc_array(nC2);
double *C2 = C2a->a;
//printarray_matrix2(B2a);
printarray(B2a);
mm_array2(Aa,B2a, C2a);
//printarray_matrix0(Ca1);
set_const_array(C0a, 0.);
for(k=0; k<2; k++)
for(j=0; j<3; j++)
for(i=0; i<4; i++)
C0[Ind_n(k,j,i, nC0)] = C2[Ind_n(i,j,k, nC2)];
printarray_matrix0(C0a);
printarray(C2a);


abort();
  return 0;
}
