/* evolve.c */
/* Wolfgang Tichy, 2/2019 */

#include "nmesh.h"
#include "evolve.h"

#define PR 1




/* register a list of variable lists in evosys */
void evolve_register_u(tMesh *mesh, pVLList *u)
{
  mesh->evosys->u = u;
}

/* register a RHS with evosys */
void evolve_register_u(tMesh *mesh,
                       void (*rhs)(pVLList *rhs, pVLList *u, double time))
{
  mesh->evosys->setrhs = rhs;
}




int evolve_mesh(tMesh *mesh)                  
{}

/*  */
int evolve(tNode *node)
{
  tMesh *mesh = node->pat->mesh;
  tEvoSys *evosys = mesh->evosys;

  /* do nothing if we have no vars to evolve */
  if(!evosys->u) return 0;

  if(!evosys->setrhs) errorexit("no RHS!");

  /* if there are no aux vars add them */
  if(!evosys->u_p)
  {
  
  }
  
  return 0;
}
