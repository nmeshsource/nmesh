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
  double bbox[6] = { -.5,1, -2.1,2, 4,5 };
  int n1max = 55;
  int n[3] = { 9,9,9 };
  tNlist *nlist;

  PRFs(":\n");

//tNode *tnode = alloc_node();
//mesh->pat[0]->rnode = 0;
  //realloc_patlist_in_mesh(mesh, 1);
  add_patch(mesh, bbox, n, n1max);

  nlist = make8_child_nodes(mesh->pat[0]->rnode, n);
  mesh->pat[0]->lns = replace1_in_nodelist(mesh->pat[0]->lns, nlist);
  //printnodelist(nlist);
  printmesh(mesh);

  nlist = make8_child_nodes(mesh->pat[0]->lns->next->node, n);
  //mesh->pat[0]->lns = 
  replace1_in_nodelist(mesh->pat[0]->lns->next, nlist);
  //printnodelist(nlist);
  printmesh(mesh);


  return 0;
}
