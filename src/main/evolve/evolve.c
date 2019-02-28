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
void evolve_register_rhs(tMesh *mesh,
                         void (*setrhs)(tNode *node, pVLList *rhs, pVLList *u))
{
  mesh->evosys->setrhs = setrhs;
}

/*
void std_setrhs(tNode *node, pVLList *rhs, pVLList *u)
{
  tMesh *mesh = node->pat->mesh;
  tEvoSys *evosys = mesh->evosys;

  //forList(evosys->u, i)
  geom_rhs(node, ListEntry(evosys->rhs,0), ListEntry(evosys->u,0));
  matter_rhs(node, ListEntry(evosys->rhs,1), ListEntry(evosys->u,1));
}
*/


/* evolve the entire leaf node mesh one time step forward */
int evolve_myln(tMesh *mesh)
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

    /* FIXME: for now all nodes use the same time step */
    node->dt = mesh->dt;
    node->time = mesh->time;

    evolve(node);
  }
  return 0;
}

/* evolve one node */
void evolve(tNode *node)
{
  evolve_RK4(node);
}



/* free extra VarLists */
int evolve_finalize(tMesh *mesh)
{
  tEvoSys *evosys = mesh->evosys;
  int i;

  /* do nothing if we have no vars to evolve */
  if(!evosys->u) return 0;

  /* free memory */
  printf("Freeing extra variable lists for evolution:\n");
  free_pVLList(evosys->u); /* free list only, not content */
  freeall_pVLList(evosys->w, vlfree,0); /* free list and its content */
  freeall_pVLList(evosys->rhs, vlfree,0);
  freeall_pVLList(evosys->u_p, vlfree,0);
  for(i=0; i<NUTEMP; i++)
    freeall_pVLList(evosys->s[i], vlfree,0);

  return 0;
}
