/* mesh.c */
/* Wolfgang Tichy, April 2005 */


#include "nmesh.h"
#include "MemoryMan.h"

#define PR 0




/* initialize mesh 
   called in main()
   here we sort out the various parameter options before calling
   for example make_mesh_box
*/
// ...

/* make an empty mesh, into which we an then initialize or into which
   we can e.g. copy the contents of an existing mesh */
tMesh *make_empty_mesh(int pr)
{
  tMesh *mesh;

  /* print info */
  if(pr) 
  {
    prdivider(0);
    PRF;
  }

  mesh = alloc_mesh();
  //if(pr) printf("g->nboxes=%d  g->box=%p  nvariables=%d\n",
  //               g->nboxes, g->box, nvariables);

  if(pr) printmesh(mesh);

  /* return pointer to newly created mesh */
  return mesh;
}
