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
int add_patch(tMesh *mesh, int nroot[3], int nmax[3])
{
  tPat *pat;
  int p = mesh->npats;
  int nD = max3(nmax[0],nmax[1],nmax[2]) * 2;

  /* make room for new patch in mesh and then add an empty patch */
  realloc_patlist_in_mesh(mesh, p + 1);
  pat = alloc_patch(mesh, p, nD);
  mesh->pat[p]  = pat;

  /* set diff matrices */

  /* setup root node */
  pat->rnode = make_root_node(pat, nroot, 0);
  pat->lns = alloc_nodelist(pat->rnode);

  return 0;
}



/* a function just for testing */
int setup_test_mesh(tMesh *mesh)
{
  int nmax[3] = { 5,5,5 };
  int n[3]    = { 3,3,3 };
  PRFs(":\n");

  //realloc_patlist_in_mesh(mesh, 1);
  add_patch(mesh, n, nmax);
  make8_child_nodes(mesh->pat[0]->rnode, nmax);
  printmesh(mesh);


  return 0;
}
