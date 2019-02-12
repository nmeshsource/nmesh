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
int add_patch(tMesh *mesh, double bbox[6], int nroot[3], int nmax)
{
  tNlist *nlist;
  tPat *pat;
  int p = mesh->npats;
  int i, ni, dir;
  double *winterp = dmalloc(nmax);

  /* make room for new patch in mesh and then add an empty patch */
  realloc_patlist_in_mesh(mesh, p + 1);
  pat = alloc_patch(mesh, p, nmax);
  mesh->pat[p] = pat;

  /* set bbox */
  for(i=0; i<6; i++) pat->bbox[i] = bbox[i];

  /* set diff, and other matrices */
  for(ni=1; ni<=nmax; ni++)
  {
    for(dir=0; dir<3; dir++)
    {
      double *Xb = pat->Xb[ni][dir]->a;
      double *Winteg = pat->Winteg[ni][dir]->a;
      double *DT = pat->Dt[ni][dir]->a;
      double *AT = pat->At[ni][dir]->a;
      double *ST = pat->St[ni][dir]->a;

      LGL_x_winteg(ni, Xb, Winteg);
      Lagrange_winterp(ni, Xb, winterp);
      LGL_DT(ni, Xb, winterp, DT);
      LGL_AT_ST_matrices(ni, Xb, Winteg, AT, ST);
    }
  }

  /* setup root node */
  pat->rnode = make_root_node(pat, nroot, 0);
  /* add root node to global mesh->lns list */
  nlist = alloc_nodelist(pat->rnode);
  append_nodelist_to_mesh_lns_myln(mesh, nlist);

  free(winterp);
  return 0;
}

void test_array_thingies(tMesh *mesh)
{
  int n[3];
  int i,j,k;
  double  A[6] = { 1,2,
                   3,4,
                   5,6 };
  /*
  double B0[24] = { 1,2,3,4,5,6,7,8,9,10,11,12,
                    13,14,15,16,17,18,19,20,21,22,23,24 };
  */
  double B0[24], B1[24], B2[24];
  //double AB[36];
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
  //Aa->a = A;
  point_array_a_to_data(Aa, A);
  //Aa->n[0]=3; Aa->n[1]=2; Aa->n[1]=1;

  tArray *B0a = alloc_array(nB0);
  //B0a->a = B0;
  point_array_a_to_data(B0a, B0);
  tArray *B1a = alloc_array(nB1);
  //B1a->a = B1;
  point_array_a_to_data(B1a, B1);
  tArray *B2a = alloc_array(nB2);
  //B2a->a = B2;
  point_array_a_to_data(B2a, B2);

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

  //free_array(Aa);
  //free_array(B0a);
  //free_array(B1a);
  //free_array(B2a);
  free(Aa);
  free(B0a);
  free(B1a);
  free(B2a);
  free_array(C0a);
  free_array(C1a);
  free_array(C2a);
}


/* a function just for testing */
int setup_test_mesh(tMesh *mesh)
{
  double bbox[6] = { -4,4, -2,2, -1,1 };
  int n1max = 55;
  int n[3] = { 5,5,5 };
  tNlist *el, *el2;
  int i;

  PRFs(":\n");


//tNode *tnode = alloc_node();
//mesh->pat[0]->rnode = 0;
  //realloc_patlist_in_mesh(mesh, 1);
  add_patch(mesh, bbox, n, n1max);

  enablevar(mesh, Ind("X"));
//  tNlist *nlist;
//  tNode *nd;
//  nd = mesh->pat[0]->rnode;
//  nlist = make8_child_nodes(nd, n);
//  replace1_in_mesh_lns_myln(mesh->lns, nlist);
  make8children_in_mesh_lns_myln(mesh->lns, n);

  //printnodelist(nlist);
  printmesh(mesh);

  el = mesh->lns;
  for(i=1; i<=1; i++) el = el->next;
//  nd = el->node;
//  nlist = make8_child_nodes(nd, n);
//  replace1_in_mesh_lns_myln(el, nlist);
  make8children_in_mesh_lns_myln(el, n);

  el = mesh->lns;
  for(i=1; i<=8+2; i++) el = el->next;
//  nd = el->node;
//  nlist = make8_child_nodes(nd, n);
//  replace1_in_mesh_lns_myln(el, nlist);
  el = make8children_in_mesh_lns_myln(el, n);

  //printnodelist(nlist);
  printmesh(mesh);
  printnodelist_and_neighbors(mesh->lns);

  destroy8siblings_in_mesh_lns_myln(el);
  printmesh(mesh);
  printnodelist_and_neighbors(mesh->lns);

  //test_array_thingies(mesh);
  //abort();

  //printarray(mesh->lns->next->node->St[1]);
  printarray_matrix0(mesh->lns->next->node->St[1]);
  printarray_matrix0(mesh->lns->next->node->Dt[1]);

  printarray(mesh->lns->next->node->Xb[1]);
  printarray(mesh->lns->next->node->Winteg[1]);

//  el = mesh->lns;
//  printnodelist(el);

  el = alloc_nodelist(mesh->pat[0]->rnode->child[1]);
  printnode_and_neighbors(el->node);

Yo(1);
  printnodelist(el);


Yo(2);
  el2 = all_descendants_along_face(el, 0, &i);
Yo(3);
  printf("i=%d\n",i);
  printnodelist(el2);

//  tNlist *find_patch_neighbors(tNode *node, int face)
//  el = 
free_nodelist(el);
free_nodelist(el2);


  return 0;
}
