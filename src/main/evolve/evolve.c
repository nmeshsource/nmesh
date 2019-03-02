/* evolve.c */
/* Wolfgang Tichy, 2/2019 */

#include "nmesh.h"
#include "evolve.h"

#define PR 0

/* The functions below are a little complex because they deal with lists of
   variable lists and a list of RHSs (one for each VarList). The was done to
   be able to couple two systems such as e.g. Z4 and matter. Each have their
   own vars, but to compute the Z4 RHS one needs the stress-energy tensor T
   (a src) that depends on the matter, while the matter RHS needs the ADM
   metric (another src) that depends on Z4. */


/* register a list of variable lists and its RHS and source functions
   in evosys */
void evolve_register_subsys_u_rhs_src(tMesh *mesh, tVarList *u,
                                      FuncPointer rhs, FuncPointer src)
{
  tEvoSys *evosys = mesh->evosys;

  /* allocate lists in evosys */
  if(!evosys->u)      evosys->u      = alloc_pVLList();
  if(!evosys->setrhs) evosys->setrhs = alloc_FuncPointerList();
  if(!evosys->setsrc) evosys->setsrc = alloc_FuncPointerList();

  /* add u, rhs, src to lists in evosys */
  push_pVLList(evosys->u, u);
  push_FuncPointerList(evosys->setrhs, rhs);
  push_FuncPointerList(evosys->setsrc, src);
}

/* free extra VarLists and other Lists */
int evolve_finalize(tMesh *mesh)
{
  tEvoSys *evosys = mesh->evosys;
  int i;

  /* do nothing if we have no vars to evolve */
  if(!evosys->u) return 0;

  PRFs(":\n");

  /* free memory in varlists */
  printf("Freeing extra variable lists for evolution:\n");
  freeall_pVLList(evosys->u, vlfree,0); /* free list and its content */
  freeall_pVLList(evosys->w, vlfree,0); /* free list and its content */
  freeall_pVLList(evosys->rhs, vlfree,0);
  freeall_pVLList(evosys->u_p, vlfree,0);
  for(i=0; i<NEVOTEMP; i++)
    freeall_pVLList(evosys->s[i], vlfree,0);

  /* free Lists */
  printf("Freeing rhs lists for evolution:\n");
  //free_pVLList(evosys->u); /* free list only, not content */
  free_FuncPointerList(evosys->setrhs);
  free_FuncPointerList(evosys->setsrc);

  /* now set all of evosys to zero */
//evolve_print_evosys(mesh);
  memset(evosys, 0, sizeof(evosys[0]));
evolve_print_evosys(mesh);

  return 0;
}

/* print evosys */
void evolve_print_evosys(tMesh *mesh)
{
  tEvoSys *evosys = mesh->evosys;
  int i;

  PRFs(":\n");
  if(evosys->u)
  {
    pr_pVLList(evosys->u);
    forList(evosys->u, i)
    {
      printf("%d: ", i);
      prvarlist(ListEntry(evosys->u,i));
    }
  }
  if(evosys->setrhs)
  {
    forList(evosys->setrhs, i)
      if(ListEntry(evosys->setrhs,i)) printf("%d: setrhs: yes\n", i);
    forList(evosys->setsrc, i)
      if(ListEntry(evosys->setsrc,i)) printf("%d: setsrc: yes\n", i);
  }
}


/* set RHS of all evo subsystems. This first also calls the setsrc
   functions in case some sources in the RHSs have to be set. */
void evolve_setrhs(tNode *node, pVLList *rhs, pVLList *u)
{
  tMesh *mesh = node->pat->mesh;
  tEvoSys *evosys = mesh->evosys;
  int i;

  /* set all sources */
  forList(u, i)
    if(ListEntry(evosys->setsrc,i))
      ListEntry(evosys->setsrc,i)(node, ListEntry(u,i));

  /* set all RHSs */
  forList(u, i)
    if(ListEntry(evosys->setrhs,i))
      ListEntry(evosys->setrhs,i)(node, ListEntry(rhs,i),
                                        ListEntry(u,i));
}



/* evolve the entire leaf node mesh one time step forward */
int evolve_myln(tMesh *mesh)
{
  tEvoSys *evosys = mesh->evosys;
  int evolve_method = Par("evolve_method");
  void (*evolve)(tNode *node); /* func pointer for evo method */
  int i, myid;

  /* do nothing if we have no vars to evolve */
  if(!evosys->u) return 0;

  if(PR) PRFs(":\n");

  if(!evosys->setrhs) errorexit("no RHS!");

  /* if there are no aux vars add them */
  if(!evosys->rhs)
  {
    /* add lists */
    evosys->w   = alloc_pVLList();
    evosys->rhs = alloc_pVLList();
    evosys->u_p = alloc_pVLList();

    printf("Adding variables for RK4 evolution:\n");
    forList(evosys->u, i)
    {
      tVarList *u   = ListEntry(evosys->u, i);

      push_pVLList(evosys->w,   AddDuplicateEnable(u, "_w", 1,-1));
      push_pVLList(evosys->rhs, AddDuplicateEnable(u, "_r", 1,0));
      push_pVLList(evosys->u_p, AddDuplicateEnable(u, "_p", 1,0));
      //push_pVLList(evosys->s[0], AddDuplicateEnable(u, "_s0", 1,0));
    }
    //printf("evosys->w = %p\n", evosys->w);
  }

  /* select evo method */
  if(Getv(evolve_method, "RK"))
    evolve = evolve_RK4;
  else if(Getv(evolve_method, "Euler"))
    evolve = evolve_Euler;

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
