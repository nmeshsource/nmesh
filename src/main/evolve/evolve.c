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

  /* since this is the start of an evo step reset evo_troubled flag here: */
  trouble_reset_evo_troubled_mesh(mesh);

  /* how we evolve the mesh */
  if(allnodes)
  {
    int trouble_score;

    /* make one full evo step */
    Evolve_mesh(mesh);

    /* get trouble score, and redo step with fv or switch back to dg */
    trouble_score = evolve_set_trouble_score_mesh(mesh);
    if(trouble_score>0)
    {
      //prdivider(0);
      PRF;printf(": trouble_score=%d (bad) => switch troubled nodes & "
                 "redo evo step\n", trouble_score);
      //node = node_from_nodename(mesh, "0_465");
      //if(node && node->dat) printelm_nodeinfo(node);
      //if(node) { Yo(89); GRHD_cons2prim_needs_cons_fix(node, ListEntry(evosys->u_p,0)); }
      /* go back to u_p and switch to fv */
      evolve_prepare_do_over_mesh(mesh);
      //if(node && node->dat) printelm_nodeinfo(node);
      //if(node) { Yo(94); GRHD_cons2prim_needs_cons_fix(node, ListEntry(evosys->u,0)); }
      //tPoint pt[] =  {{.node=node, .ijk=151}};
      //printvarlist_atpoint(pt, ListEntry(evosys->u,0));
      /* Now all new fv nodes have newly interpolated evo vars, so we need
         to limit them again. ALSO, we need to run PRELIM again to update
         fields (like the ADM metric) that are not in evosys->u to the new
         evosys->u = evosys->u_p set in evolve_prepare_do_over_mesh. */
      //evolve_limiter_mesh(mesh, evosys->u, 1); //but only if trbl_score>0
      evolve_limiter_mesh(mesh, evosys->u, 0); //in all nodes
      //if(node && node->dat) printelm_nodeinfo(node);
      //if(node) GRHD_cons2prim_needs_cons_fix(node, ListEntry(evosys->u,0));

      /* redo evo step */
      //PRF;printf(": redo evo step\n");
      Evolve_mesh(mesh);
      //prdivider(0);
    }
    else if(trouble_score<=-NOTROUBLES)
    {
      //prdivider(0);
      PRF;printf(": trouble_score=%d (great) => switch nontroubled nodes\n",
                 trouble_score);
      /* switch to dg */
      evolve_switch_nontroubled_nodes_mesh(mesh);
      /* now some aux vars (and others) are not set */
      //prdivider(0);
    }

    /* we limit the final u only here */
    evolve_limiter_mesh(mesh, evosys->u, 0);

    /* update some vars by calling funcs in PRESURF, SETSRC
       often PRESURF does cons2prim, SETSRC sets stress-energy */
    if(trouble_score<=-NOTROUBLES)
      evolve_setsrc_again_nontroubled_nodes_mesh(mesh, evosys->rhs,
                                                 evosys->u);
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

  if(PR) PRFs(":\n");

  /* loop over all nodes before MPI requests */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    int i, troubled;

    /* set time on all nodes */
    node->time = mesh->time;
    node->dt   = mesh->dt;

    /* set things before surface exchange, e.g. cons2prim */
    troubled = 0;
    forList(u, i)
    {
      tVarList *vlr = ListEntry(rhs,i);
      tVarList *vlu = ListEntry(u,i);
      if(ListEntry(evosys->f[PRESURF],i))
        troubled |= ListEntry(evosys->f[PRESURF],i)(node, vlr, vlu);
    }
    node->dat->info->evo_troubled |= troubled;
  }

  /* do surface exchange on entire mesh */
  MPIexchange_set_all_myln_localdata(mesh);
  MPIexchange_request_all_myln_data(mesh);

  /* Now work on things that do not depend on surface data:
     I.e. we overlap the communication started by
     MPIexchange_request_all_myln_data with other calculations. */

  /* loop over all nodes after MPI exchange has been requested */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    int i, troubled;

    troubled = 0;
    forList(u, i)
    {
      /* set all sources */
      if(ListEntry(evosys->f[SETSRC],i))
        troubled |= ListEntry(evosys->f[SETSRC],i)(node, ListEntry(rhs,i),
                                                   ListEntry(u,i));
      /* set all volume RHSs */
      if(ListEntry(evosys->f[VOLRHS],i))
        troubled |= ListEntry(evosys->f[VOLRHS],i)(node, ListEntry(rhs,i),
                                                   ListEntry(u,i));
    }
    node->dat->info->evo_troubled |= troubled;
  }

  /* After we have done all we can without the surface data, we now wait
     until we get all the surface data: */

  /* get surfaces so that we can compute fluxes */
  MPIexchange_get_all_myln_data(mesh);

  /* loop over all nodes, after MPI data has been received */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    int i, troubled;

    /* add all surface RHSs */
    troubled = 0;
    forList(u, i)
      if(ListEntry(evosys->f[SURFRHS],i))
        troubled |= ListEntry(evosys->f[SURFRHS],i)(node, ListEntry(rhs,i),
                                                    ListEntry(u,i));
    node->dat->info->evo_troubled |= troubled;
  }
}

/* Apply limiters to evo subsystems.
   We do it only if trbl_score >= badlim or trbl_score<=goodlim.
   If opt=0 badlim=INT_MIN so that we ALWAYS do it.
   If opt=1 badlim=1       so that we do it only for switched nodes. */
/* Version for entire mesh: */
void evolve_limiter_mesh(tMesh *mesh, pVLList *u, int opt)
{
  tEvoSys *evosys = mesh->evosys;
  int j;
  tVarList *vl;
  int goodlim = -NOTROUBLES;
  int badlim;

  if(opt==1) badlim = 1;
  else       badlim = INT_MIN;

  if(PR) PRFs(":\n");

  /* loop over all nodes before MPI requests */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    int trbl_score = node->dat->info->trbl_score;
    int i, troubled;

    if(trbl_score >= badlim || trbl_score<=goodlim)
    {
      troubled = 0;
      forList(u, i)
      {
        /* call funcs that we need before limiters */
        if(ListEntry(evosys->f[PRELIM],i))
          troubled |= ListEntry(evosys->f[PRELIM],i)(node, ListEntry(u,i));

        /* set data limiter needs in myindc arrays of each node */
        if(ListEntry(evosys->f[LIMDATA],i))
          troubled |= ListEntry(evosys->f[LIMDATA],i)(node, ListEntry(u,i));
      }
      node->dat->info->evo_troubled |= troubled;
    } /* end if */
  }

  /* create varlist that needs MPI exchange */
  vl = vlalloc(mesh);
  forList(u, j)
  {
    if(ListEntry(evosys->f[LIMDATA],j))
      vlpushvl(vl, ListEntry(u,j)); //add ListEntry(u,j) to vl
  }
  /* exchange data in indicators (indc) if needed */
  if(VLn(vl)>0)
  {
    /* initiate indc exchange */
    request_all_myln_indc_exchange_for_vl(mesh, vl);
    /* After this we could do work. MPI is now busy sending buffers */

    /* now get the indicators and wait for MPI buffers if necessary */
    get_all_myln_indc_for_vl(mesh, vl);
  }
  vlfree(vl);

  /* loop over all nodes after MPI exchange */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    int trbl_score = node->dat->info->trbl_score;
    int i;

    if(trbl_score >= badlim || trbl_score<=goodlim)
    {
      forList(u, i)
      {
        /* apply limiter */
        if(ListEntry(evosys->f[LIMITER],i))
        {
          int ret = ListEntry(evosys->f[LIMITER],i)(node, ListEntry(u,i));

          /* increase nlim if limiter was active, otherwise reset nlim */
          if(ret) node->dat->info->nlim += 1;
          else    node->dat->info->nlim = 0;

          /* also count a ret!=0 as trouble */
          node->dat->info->evo_troubled |= ret;
        }
      }
    } /* end if */
  }
}

/* update some vars by calling funcs in PRESURF, SETSRC
   often PRESURF does cons2prim, SETSRC sets stress-energy,
   PRELIM sets ADM metric */
void evolve_setsrc_again_nontroubled_nodes_mesh(tMesh *mesh,
                                                pVLList *rhs, pVLList *u)
{
  tEvoSys *evosys = mesh->evosys;

  if(PR) PRFs(":\n");

  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    int i;

    /* set some things again, e.g. cons2prim */
    /* this is only needed on nodes that had a very neg. trouble score and
       were thus converted to dg */
    if(node->dat->info->trbl_score <= -NOTROUBLES)
      forList(u, i)
      {
        tVarList *vlr = ListEntry(rhs,i);
        tVarList *vlu = ListEntry(u,i);

        if(ListEntry(evosys->f[PRESURF],i))
          ListEntry(evosys->f[PRESURF],i)(node, vlr, vlu);

        if(ListEntry(evosys->f[SETSRC],i))
          ListEntry(evosys->f[SETSRC],i)(node, vlr, vlu);

        /* we do not need to run PRELIM, since evolve_limiter_mesh will
           run right before evolve_setsrc_mesh is called */
        /* if(ListEntry(evosys->f[PRELIM],i))
             ListEntry(evosys->f[PRELIM],i)(node, vlu);*/
      }
  }
}


/* apply filters to all evo subsystems */
int evolve_filter_evosys_mesh(tMesh *mesh)
{
  tEvoSys *evosys = mesh->evosys;
  pVLList *u = evosys->u;
  int i;
  int evolve_filter         = Par("evolve_filter");
  int evolve_filter_varlist = Par("evolve_filter_varlist");
  int filter_all_evovars = Getb(evolve_filter);
  int filter_varlist     = GetLen(evolve_filter_varlist);

  if(PR) PRFs(":\n");

  if(filter_all_evovars)
  {
    double af = Getd(Par("evolve_filter_alp"));
    double sf = Getd(Par("evolve_filter_s"));

    if(PR) { PRF;printf(": filtering all evolution vars in evosys->u\n"); }

    /* loop over list of varlists and filter each varlist */
    forList(u, i)
    {
      tVarList *vl = ListEntry(u,i);
      expfilter_vl(vl, af, sf);
    }
  }

  if(filter_varlist)
  {
    double af = Getd(Par("evolve_filter_alp"));
    double sf = Getd(Par("evolve_filter_s"));
    char *list = strdup(Gets(evolve_filter_varlist));
    char *name, *saveptr;
    tVarList *vl = vlalloc(mesh);

    if(PR) { PRF;printf(": filtering varlist\n%s\n", list); }

    /* loop over list and make a varlist from it */
    for(name=strtok_r(list, " ", &saveptr); name;
        name=strtok_r(0,    " ", &saveptr))
    {
      int vi = MeshVarIndLax(mesh, name);
      if(vi>=0) vlpush(vl, vi);
    }

    if(vl->n) expfilter_vl(vl, af, sf);

    vlfree(vl);
    free(list);
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
    if(ListEntry(evosys->f[SETSRC],i))
      ListEntry(evosys->f[SETSRC],i)(node, ListEntry(rhs,i), ListEntry(u,i));

  /* set all vol. RHSs */
  forList(u, i)
    if(ListEntry(evosys->f[VOLRHS],i))
      ListEntry(evosys->f[VOLRHS],i)(node, ListEntry(rhs,i), ListEntry(u,i));

  //Test:  get_all_surfaces(node);

  /* add all surf. RHSs */
  forList(u, i)
    if(ListEntry(evosys->f[SURFRHS],i))
      ListEntry(evosys->f[SURFRHS],i)(node, ListEntry(rhs,i), ListEntry(u,i));

  /* do not free all surface info, because we call evolve_setrhs repeatedly */
  //if(0) evolve_free_surfaces(node, u);
}
