/* evolve.c */
/* Wolfgang Tichy, 2/2019 */

#include "nmesh.h"
#include "evolve.h"

#define PR 0

/* The functions below are a little complex because they deal with lists of
   variable lists and a list of RHSs (one for each VarList). This was done to
   be able to couple two systems such as e.g. Z4 and matter. Each have their
   own vars, but to compute the Z4 RHS one needs the stress-energy tensor T
   (a src) that depends on the matter, while the matter RHS needs the ADM
   metric (another src) that depends on Z4. */


/* evolve the entire leaf node mesh one time step forward */
int evolve_myln(tMesh *mesh)
{
  tEvoSys *evosys = mesh->evosys;
  int evolve_method = Par("evolve_method");
  void (*Evolve)(tNode *node) = NULL;      /* func pointer for evo method */
  void (*Evolve_mesh)(tMesh *mesh) = NULL; /* func pointer for evo method */
  int allnodes = 1;
  tVarList *allu = NULL;
  tVarList *u_ch = NULL;
  int i;

  /* do nothing if we have no vars to evolve */
  if(!evosys->u) return 0;

  if(PR) PRFs(":\n");

  /* make aux vars if needed */
  evolve_init_evosys(mesh);

  /* select evo method */
  if(Getv(evolve_method, "RK4"))
  {
    Evolve = evolve_RK4;
    Evolve_mesh = evolve_RK4_mesh;
  }
  else if(Getv(evolve_method, "Euler"))
  {
    Evolve = evolve_Euler;
    Evolve_mesh = evolve_Euler_mesh;
  }
  else if(Getv(evolve_method, "sspRK3"))
  {
    Evolve = NULL;
    Evolve_mesh = evolve_sspRK3_mesh;
  }
  else
    errorexits("unknown value:   evolve_method = %s", Gets(evolve_method));

  /* make varlist allu with all in evosys->u */
  allu = vlalloc(mesh);
  forList(evosys->u, i) vlpushvl(allu, ListEntry(evosys->u,i));

  /* set u_ch = allu */
  if(Getb(Par("evolve_compute_change")))
  {
    u_ch = AddDuplicateEnable(allu, "_change", AUXVAR, 0);
    vlcopy(u_ch, allu); /* u_ch = allu */
  }

  /* how we evolve the mesh */
  if(allnodes)
  {
    Evolve_mesh(mesh);
  }
  else /* evolve each node on its own */
  {
    /* evolve each node */
    formylnodes(mesh)
    {
      tNode *node = MyLnode;

      /* FIXME: for now all nodes use the same time step */
      node->dt = mesh->dt;
      node->time = mesh->time;

      Evolve(node);
    }
  }

  /* set u_ch = allu - u_ch after evo step */
  if(Getb(Par("evolve_compute_change")))
    vladd(u_ch, 1.,allu, -1.,u_ch); /* u_ch = allu - u_ch */

  /* we don't need allu, u_ch anymore */
  vlfree(u_ch);
  vlfree(allu);

  return 0;
}


/* set RHS of all evo subsystems. This first also calls the setsrc
   functions in case some sources in the RHSs have to be set. */
/* Version for entire mesh: */
void evolve_setrhs_mesh(tMesh *mesh, pVLList *rhs, pVLList *u)
{
  tEvoSys *evosys = mesh->evosys;
  int i;

  if(PR) PRFs(":\n");

  /* set things before surface exchange, e.g. cons2prim */
  forList(u, i)
  {
    tVarList *vlr = ListEntry(rhs,i);
    tVarList *vlu = ListEntry(u,i);

    if(ListEntry(evosys->presurf,i))
    {
      formylnodes(mesh)
      {
        tNode *node = MyLnode;
        ListEntry(evosys->presurf,i)(node, vlr, vlu);
      }
    }
  }

  /* set time on all nodes */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;

    node->time = mesh->time;
    node->dt   = mesh->dt;
  }

  /* do surface exchange on entire mesh */
  MPIexchange_set_all_myln_localdata(mesh);
  MPIexchange_request_all_myln_data(mesh);

  /* set all sources */
  forList(u, i)
    if(ListEntry(evosys->setsrc,i))
      ListEntry(evosys->setsrc,i)(mesh, ListEntry(rhs,i), ListEntry(u,i));

  /* set all volume RHSs */
  forList(u, i)
    if(ListEntry(evosys->volrhs,i))
      ListEntry(evosys->volrhs,i)(mesh, ListEntry(rhs,i), ListEntry(u,i));

  /* get surfaces so that we can compute fluxes */
  MPIexchange_get_all_myln_data(mesh);

  /* add all surface RHSs */
  forList(u, i)
    if(ListEntry(evosys->surfrhs,i))
      ListEntry(evosys->surfrhs,i)(mesh, ListEntry(rhs,i), ListEntry(u,i));
}

/* apply limiters to evo subsystems. */
/* Version for entire mesh: */
void evolve_limiter_mesh(tMesh *mesh, pVLList *u)
{
  tEvoSys *evosys = mesh->evosys;
  int i;

  if(PR) PRFs(":\n");

  /* loop over list of varlists and call funcs that we need before limiters */
  forList(u, i)
  {
    tVarList *vl = ListEntry(u,i);

    if(ListEntry(evosys->prelim,i))
    {
      /* call prelim functions */
      formylnodes(mesh)
      {
        tNode *node = MyLnode;
        ListEntry(evosys->prelim,i)(node, vl);
      }
    }
  }

  /* loop over list of varlists */
  forList(u, i)
  {
    tVarList *vl = ListEntry(u,i);

    /* set limiter data in indicators (indc) */
    if(ListEntry(evosys->limdata,i))
    {
      /* set data limiter needs in myindc arrays of each node */
      formylnodes(mesh)
      {
        tNode *node = MyLnode;
        ListEntry(evosys->limdata,i)(node, vl);
      }

      /* initiate indc exchange */
      request_all_myln_indc_exchange_for_vl(mesh, vl);
      /* After this we could do work. MPI is now busy sending buffers */

      /* now get the indicators and wait for MPI buffers if necessary */
      get_all_myln_indc_for_vl(mesh, vl);
    }

    /* apply limiter */
    if(ListEntry(evosys->limiter,i))
    {
      /* use limiter on each node */
      formylnodes(mesh)
      {
        tNode *node = MyLnode;
        int ret = ListEntry(evosys->limiter,i)(node, vl);

        /* increase nlim if limiter was active, otherwise reset nlim */
        if(ret) node->dat->nlim += 1;
        else    node->dat->nlim = 0;
      }
    }
  } /* end loop over list of varlists */
}


/* apply filters to all evo subsystems */
int evolve_filter_evosys_mesh(tMesh *mesh)
{
  tEvoSys *evosys = mesh->evosys;
  pVLList *u = evosys->u;
  int i;
  int filter_on = Getb(Par("evolve_filter"));

  if(filter_on)
  {
    double af = Getd(Par("evolve_filter_alp"));
    double sf = Getd(Par("evolve_filter_s"));

    if(PR) PRFs(":\n");

    /* loop over list of varlists and filter each varlist */
    forList(u, i)
    {
      tVarList *vl = ListEntry(u,i);
      expfilter_vl(vl, af, sf);
    }
  }

  return 0;
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

/* set RHS of all evo subsystems. This first also calls the setsrc
   functions in case some sources in the RHSs have to be set. */
/* Version for just one node: */
void evolve_setrhs(tNode *node, pVLList *rhs, pVLList *u, int request_surfs)
{
  tMesh *mesh = node->pat->mesh;
  tEvoSys *evosys = mesh->evosys;
  int i;

  if(PR) PRFs(":\n");

  /* request all surfaces on node for all vars in u */
  if(request_surfs)
    evolve_request_surfaces(node, u);

  /* set all sources */
  forList(u, i)
    if(ListEntry(evosys->setsrc,i))
      ListEntry(evosys->setsrc,i)(node, ListEntry(rhs,i), ListEntry(u,i));

  /* set all vol. RHSs */
  forList(u, i)
    if(ListEntry(evosys->volrhs,i))
      ListEntry(evosys->volrhs,i)(node, ListEntry(rhs,i), ListEntry(u,i));

  //Test:  get_all_surfaces(node);

  /* add all surf. RHSs */
  forList(u, i)
    if(ListEntry(evosys->surfrhs,i))
      ListEntry(evosys->surfrhs,i)(node, ListEntry(rhs,i), ListEntry(u,i));

  /* do not free all surface info, because we call evolve_setrhs repeatedly */
  //if(0) evolve_free_surfaces(node, u);
}
