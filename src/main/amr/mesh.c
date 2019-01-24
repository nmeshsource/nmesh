/* mesh.c */
/* Wolfgang Tichy, April 2005 */


#include "nmesh.h"
#include "amr.h"

#define PR 0




/* initialize mesh called in main() */

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
int add_patch(tMesh *mesh, int nmax[3])
{
  int p = mesh->npats;
  int nD = max3(nmax[0],nmax[1],nmax[2]) * 2;

  realloc_patlist_in_mesh(mesh, p + 1);
  mesh->pat[p]  = alloc_patch(mesh, p, nD);
  return 0;
}


/* a function just for testing */
int setup_test_mesh(tMesh *mesh)
{
  int nmax[3] = { 5,5,5 };
  PRFs(":\n");

  //realloc_patlist_in_mesh(mesh, 1);
  add_patch(mesh, nmax);
  printmesh(mesh);


  return 0;
}
