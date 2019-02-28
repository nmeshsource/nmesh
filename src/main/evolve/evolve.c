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
{
  tEvoSys *evosys = mesh->evosys;
  int i, myid;

  /* do nothing if we have no vars to evolve */
  if(!evosys->u) return 0;

  if(!evosys->setrhs) errorexit("no RHS!");

  /* if there are no aux vars add them */
  if(!evosys->u_p)
  {
    printf("Adding variables for RK4 evolution:\n");
    forList(evosys->u, i)
    {
      tVarList *u   = ListEntry(evosys->u, i);
      ListEntry(evosys->w  , i) = AddDuplicate(u, "_w", -1,-1);
      ListEntry(evosys->rhs, i) = AddDuplicate(u, "_r", 1,0);
      ListEntry(evosys->u_p, i) = AddDuplicate(u, "_p", 1,0);
      //ListEntry(evosys->s0 , i) = AddDuplicate(u, "_s0", 1,0);
      //ListEntry(evosys->s1 , i) = AddDuplicate(u, "_s1", 1,0);
    }
  }

  /* evolve each node */
  formylnodes(mesh, myid)
  {
    tNode *node = GetMyNode(mesh, myid);
    evolve(node);
  }

}

/*  */
int evolve(tNode *node)
{
  tMesh *mesh = node->pat->mesh;
  tEvoSys *evosys = mesh->evosys;

  
  return 0;
}


junk:
      tVarList *w   = ListEntry(evosys->w, i);
      tVarList *rhs = ListEntry(evosys->rhs, i);
      tVarList *u_p = ListEntry(evosys->u_p, i);
      tVarList *s0  = ListEntry(evosys->s[0], i);
      tVarList *s1  = ListEntry(evosys->s[1], i);
