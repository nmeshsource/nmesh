/* evosys.c */
/* Wolfgang Tichy, 7/2019 */

#include "nmesh.h"
#include "evolve.h"

#define PR 0

/* The functions below are a little complex because they deal with lists of
   variable lists and a list of RHSs (one for each VarList). This was done to
   be able to couple two systems such as e.g. Z4 and matter. Each have their
   own vars, but to compute the Z4 RHS one needs the stress-energy tensor T
   (a src) that depends on the matter, while the matter RHS needs the ADM
   metric (another src) that depends on Z4. */


/* register a variable list u and its RHS in evosys. The func. pointers will
   be called in the order they appear here. rhs is the most important one
   and contains the RHS for the evo eqn of u. The others are there to apply
   limiters or to set source terms before rhs is called */
void evolve_register_subsys(tMesh *mesh, tVarList *u,
                FuncPointer prelim, FuncPointer limdata, FuncPointer limiter,
                FuncPointer presurf, FuncPointer setsrc, FuncPointer volrhs,
                FuncPointer surfrhs)
{
  tEvoSys *evosys = mesh->evosys;
  int b;

  /* allocate lists in evosys */
  if(!evosys->u)       evosys->u       = alloc_pVLList();
  for(b=0; b<NEVOFUNCBINS; b++)
    if(!evosys->f[b]) evosys->f[b] = alloc_FuncPointerList();

  /* add u, rhs, src, ... to lists in evosys */
  push_pVLList(evosys->u, u);
  push_FuncPointerList(evosys->f[PRELIM], prelim);
  push_FuncPointerList(evosys->f[LIMDATA], limdata);
  push_FuncPointerList(evosys->f[LIMITER], limiter);
  push_FuncPointerList(evosys->f[PRESURF], presurf);
  push_FuncPointerList(evosys->f[SETSRC], setsrc);
  push_FuncPointerList(evosys->f[VOLRHS], volrhs);
  push_FuncPointerList(evosys->f[SURFRHS], surfrhs);
}

/* register a list of variable lists and its RHS, source functions and
   limiters in evosys */
void evolve_register_subsys_u_rhs_lim(tMesh *mesh, tVarList *u,
                                      FuncPointer volrhs, FuncPointer surfrhs,
                                      FuncPointer limdata,
                                      FuncPointer limiter)
{
  evolve_register_subsys(mesh, u, NULL,limdata,limiter,
                         NULL,NULL,volrhs,surfrhs);
}

/* free extra VarLists and other Lists */
int evolve_free_evosys(tMesh *mesh)
{
  tEvoSys *evosys = mesh->evosys;
  int i, b;

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
  for(b=0; b<NEVOFUNCBINS; b++)
    free_FuncPointerList(evosys->f[b]);


  /* now set all of evosys to zero */
  //evolve_print_evosys(mesh);
  memset(evosys, 0, sizeof(evosys[0]));
  //evolve_print_evosys(mesh);

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
    //pr_pVLList(evosys->u);
    forList(evosys->u, i)
    {
      printf("%d: ", i);
      prvarlist(ListEntry(evosys->u,i));
      if(ListEntry(evosys->f[PRELIM],i))  printf("%d: PRELIM:  yes\n", i);
      if(ListEntry(evosys->f[LIMDATA],i)) printf("%d: LIMDATA: yes\n", i);
      if(ListEntry(evosys->f[LIMITER],i)) printf("%d: LIMITER: yes\n", i);
      if(ListEntry(evosys->f[PRESURF],i)) printf("%d: PRESURF: yes\n", i);
      if(ListEntry(evosys->f[SETSRC],i))  printf("%d: SETSRC:  yes\n", i);
      if(ListEntry(evosys->f[VOLRHS],i))  printf("%d: VOLRHS:  yes\n", i);
      if(ListEntry(evosys->f[SURFRHS],i)) printf("%d: SURFRHS: yes\n", i);
    }
  }
  printf("evolve function pointers are called in this order:\n");
  if(evosys->f[LIMDATA])
  {
    forList(evosys->f[PRELIM], i)
      if(ListEntry(evosys->f[PRELIM],i))  printf("%d: PRELIM\n", i);
    forList(evosys->f[LIMDATA], i)
      if(ListEntry(evosys->f[LIMDATA],i)) printf("%d: LIMDATA\n", i);
    forList(evosys->f[LIMITER], i)
      if(ListEntry(evosys->f[LIMITER],i)) printf("%d: LIMITER\n", i);
  }
  if(evosys->f[VOLRHS])
  {
    forList(evosys->f[PRESURF], i)
      if(ListEntry(evosys->f[PRESURF],i)) printf("%d: PRESURF\n", i);
    forList(evosys->f[SETSRC], i)
      if(ListEntry(evosys->f[SETSRC],i))  printf("%d: SETSRC\n", i);
    forList(evosys->f[VOLRHS], i)
      if(ListEntry(evosys->f[VOLRHS],i))  printf("%d: VOLRHS\n", i);
    forList(evosys->f[SURFRHS], i)
      if(ListEntry(evosys->f[SURFRHS],i)) printf("%d: SURFRHS\n", i);
  }
}


/* make some vars and put them in evosys */
int evolve_init_evosys(tMesh *mesh)
{
  tEvoSys *evosys = mesh->evosys;
  int i;

  /* do nothing if we have no vars to evolve */
  if(!evosys->u) return 0;

  if(PR) PRFs(":\n");

  if(!evosys->f[VOLRHS]) errorexit("no RHS!");

  /* if there are no aux vars add them */
  if(!evosys->rhs)
  {
    /* free surfaces since we are adding EvoVars with more surfaces */
    evolve_free_communication_structs(mesh);

    /* add lists */
    evosys->w   = alloc_pVLList();
    evosys->rhs = alloc_pVLList();
    evosys->u_p = alloc_pVLList();

    printf("Adding variables for RK evolution:\n");
    forList(evosys->u, i)
    {
      tVarList *u   = ListEntry(evosys->u, i);

      push_pVLList(evosys->w,   AddDuplicateEnable(u, "_w", AUXVAR,-1));
      push_pVLList(evosys->rhs, AddDuplicateEnable(u, "_r", AUXVAR,0));
      push_pVLList(evosys->u_p, AddDuplicateEnable(u, "_p", AUXVAR,0));
      //push_pVLList(evosys->s[0], AddDuplicateEnable(u, "_s0", AUXVAR,0));
    }
    //printf("evosys->w = %p\n", evosys->w);

    /* now that we have more vars re-init surfaces */
    evolve_init_communication_structs(mesh);
  }
  return 0;
}

/* check if var vi is in list w */
int var_in_pVLList(pVLList *w, int vi)
{
  if(w)
  {
    int i;
    forList(w, i)
    {
      tVarList *vl = ListEntry(w, i);
      if(!vl) continue;
      if(vlindex(vl, vi) >= 0) return 1;
    }
  }
  return 0;
}

/* return 1 if var vi is one of the vars that evolve_init_evosys has added */
int var_added_by_evolve_init_evosys(tMesh *mesh, int vi)
{
  tEvoSys *evosys = mesh->evosys;
  int i;

  if(var_in_pVLList(evosys->w, vi)) return 1;

  if(var_in_pVLList(evosys->rhs, vi)) return 1;

  if(var_in_pVLList(evosys->u_p, vi)) return 1;

  for(i=0; i<NEVOTEMP; i++)
    if(var_in_pVLList(evosys->s[i], vi)) return 1;

  return 0;
}

/* init structs that are used for node to neighbor node communication */
int evolve_init_communication_structs(tMesh *mesh)
{
  /* init node to nb surface exchanges */
  MPIexchange_init_all_myln(mesh);

  /* int node to nb indc exchanges */
  init_all_myln_myindc_in_evosys(mesh);
  return 0;
}

/* free structs that are used for node to neighbor node communication */
int evolve_free_communication_structs(tMesh *mesh)
{
  /* free indc */
  free_all_myln_myindc_in_evosys(mesh);

  /* free surfaces */
  MPIexchange_free_all_myln(mesh);
  return 0;
}

/* init all indc on all nodes in the mesh for all varlists in evosys */
void init_all_myln_myindc_in_evosys(tMesh *mesh)
{
  tEvoSys *evosys = mesh->evosys;
  int i;

  PRFs(":\n");

  /* check if evo vars in u need myindc */
  if(evosys->u)
  {
    forList(evosys->u, i)
    {
      tVarList *vl = ListEntry(evosys->u,i);

      if(ListEntry(evosys->f[LIMDATA],i))
      {
        /* NOTE: ListEntry(evosys->f[LIMDATA],i)(NULL, vl)
                 must return number of data vals we need */
        int nvals = ListEntry(evosys->f[LIMDATA],i)(NULL, vl);
        if(nvals>0)
          init_all_myln_myindc_for_vl(mesh, vl, nvals);
      }
    } /* end forList */
  }

  /* check if evo vars in w need myindc */
  if(evosys->w)
  {
    forList(evosys->w, i)
    {
      tVarList *vl = ListEntry(evosys->w,i);

      if(ListEntry(evosys->f[LIMDATA],i))
      {
        /* NOTE: ListEntry(evosys->f[LIMDATA],i)(NULL, vl)
                 must return number of data vals we need */
        int nvals = ListEntry(evosys->f[LIMDATA],i)(NULL, vl);
        if(nvals>0)
          init_all_myln_myindc_for_vl(mesh, vl, nvals);
      }
    } /* end forList */
  }
}

/* free all indc on all nodes in the mesh for varlists in evosys */
void free_all_myln_myindc_in_evosys(tMesh *mesh)
{
  tEvoSys *evosys = mesh->evosys;
  int i;

  PRFs(":\n");

  /* check if evo vars in u need myindc */
  if(evosys->u)
  {
    forList(evosys->u, i)
    {
      tVarList *vl = ListEntry(evosys->u,i);
      free_all_myln_indc_for_vl(mesh, vl);
    } /* end forList */
  }

  /* check if evo vars in w need myindc */
  if(evosys->w)
  {
    forList(evosys->w, i)
    {
      tVarList *vl = ListEntry(evosys->w,i);
      free_all_myln_indc_for_vl(mesh, vl);
    } /* end forList */
  }
}


/*************************************************************************/
/* NOTE: functions below do not work yet !!! */
/*************************************************************************/

/*************************************************************************/
/* functions to evolve on just one node
   will work only once request_all_vl_surfaces and such in
   main/amr/surface.c start working
*/
/*************************************************************************/

/* request all surfaces on node for all vars in u */
void evolve_request_surfaces(tNode *node, pVLList *u)
{
  tMesh *mesh = node->pat->mesh;
  tVarList *allu = vlalloc(mesh);
  int i;

  if(PR) PRFs(":\n");

  /* 1. make varlist with all in u */
  forList(u, i)
    vlpushvl(allu, ListEntry(u,i));

  /* 2. now start surface requests */
  set_all_vl_mysurf(node, allu);
  request_all_vl_surfaces(node, allu);

  /* we don't need allu anymore */
  vlfree(allu);
}

/* free all surfaces */
void evolve_free_surfaces(tNode *node, pVLList *u)
{
  tMesh *mesh = node->pat->mesh;
  tVarList *allu = vlalloc(mesh);
  int i;

  if(PR) PRFs(":\n");

  /* 1. make varlist with all in u */
  forList(u, i)
    vlpushvl(allu, ListEntry(u,i));

  /* 2. now free */
  free_all_vl_surfaces(node, allu);

  /* we don't need allu anymore */
  vlfree(allu);
}
