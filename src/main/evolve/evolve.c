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
    vlcopy(NULL, u_ch, allu); /* u_ch = allu */
  }

  /* how we evolve the mesh */
  if(allnodes)
  {
    int trouble_score;
    //tNode *node = node_from_nodename(mesh, "0_366");
    //tPoint pt[] =  {{.node=node, .ijk=17}};

    /* make one full evo step */
    Evolve_mesh(mesh);

    /* get global trouble score */
    trouble_score = evolve_set_trouble_score_mesh(mesh);

    /* If trouble_score>0 at least some nodes are troubled.
       In this case we switch them to fv and redo the evo step. */
    if(trouble_score>0)
    {
      PRF;printf(": trouble_score=%d (bad) => switch troubled nodes & "
                 "redo evo step\n", trouble_score);
      //if(node) { Yo(80); printelm_nodeinfo(node); }
      //if(node) { GRHD_cons2prim_needs_cons_fix(node, ListEntry(evosys->u_p,0)); }
      //printvarlist_atpoint(pt, ListEntry(evosys->u_p,0), "");
      //printvarlist_atpoint(pt, ListEntry(evosys->u_p,1), "");
      /* go back to u_p and switch to fv */
      evolve_prepare_do_over_mesh(mesh);
      //if(node) { Yo(84); printelm_nodeinfo(node); }
      //if(node) { GRHD_cons2prim_needs_cons_fix(node, ListEntry(evosys->u,0)); }
      //printvarlist_atpoint(pt, ListEntry(evosys->u,0), "");
      //printvarlist_atpoint(pt, ListEntry(evosys->u,1), "");
      /* Now all new fv nodes have newly interpolated evo vars, so we need
         to limit them again. ALSO, we need to run PRELIM again to update
         fields (like the ADM metric) that are not in evosys->u to the new
         evosys->u = evosys->u_p set in evolve_prepare_do_over_mesh. */
      //evolve_limiter_mesh(mesh, evosys->u, 1); //but only if trbl_score>0
      evolve_limiter_mesh(mesh, evosys->u, 0); //in all nodes
      refine_set_datinfo_unlimited_mesh(mesh, 0); //mark elms as now limited
      //if(node) { Yo(90); printelm_nodeinfo(node); }
      //if(node) { GRHD_cons2prim_needs_cons_fix(node, ListEntry(evosys->u,0)); }
      //printvarlist_atpoint(pt, ListEntry(evosys->u,0), "");
      //printvarlist_atpoint(pt, ListEntry(evosys->u,1), "");

      /* redo evo step */
      //PRF;printf(": redo evo step\n");
      Evolve_mesh(mesh);
      /* It may be good to update trouble_score here. But we assume that
         the nodes with trbl_score<0 in the 1st step still have <0 when
         we redo the step. Also, if we call evolve_set_trouble_score_mesh
         again, trbl_score of some nodes may further decrease ... */
    }

    /* Since this is the end of the evo step (that we now consider
       successful), reset evo_troubled flag here. If any evo trouble
       happens aftrwards (e.g. when we call the limiter later) this should
       be taken into account when evolve_set_trouble_score_mesh is called
       the next time. */
    trouble_reset_evo_troubled_mesh(mesh);

    /* If trouble_score>0 or trouble_score<=-NOTROUBLES, some nodes may have
       trbl_score<=-NOTROUBLES. We now switch these nodes back to dg. */
    if(trouble_score>0 || trouble_score<=-NOTROUBLES)
    {
      if(trouble_score<=-NOTROUBLES) //all nodes have trbl_score<=0
      {
        PRF;printf(": trouble_score=%d (great) => switch nontroubled nodes\n",
                   trouble_score);
      }
      /* switch all nodes with negative enough trbl_score to dg */
      evolve_switch_nontroubled_nodes_mesh(mesh);
      /* now some aux vars (and others) are not set */
    }

    /* we limit the final u only here */
    evolve_limiter_mesh(mesh, evosys->u, 0);
    //if(node) { Yo(100); GRHD_cons2prim_needs_cons_fix(node, ListEntry(evosys->u_p,0)); }

    /* update some vars by calling funcs in PRESURF*, SETSRC*
       often PRESURF does cons2prim, SETSRC sets stress-energy */
    if(trouble_score>0 || trouble_score<=-NOTROUBLES)
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
    vladd(NULL, u_ch, 1.,allu, -1.,u_ch); /* u_ch = allu - u_ch */

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
  pVLList *x = evosys->x; /* extra vars needed for LDG */
  int have_XRHS; /* is set to 1 if we have a XVOLRHS or a XSURFRHS for x */
  int ie;

  if(PR) PRFs(":\n");

  /* loop over all nodes before MPI requests */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    int i, troubled;

    /* time PRESURF */
    loadtimer_start(node);

    /* set time on all nodes */
    node->time = mesh->time;
    node->dt   = mesh->dt;

    /* set things before surface exchange, e.g. cons2prim */
    troubled = 0;
    forList(u, i)
    {
      tEvoVars evv[1] = {{
        .vlu = ListEntry(u,i),
        .vlr = ListEntry(rhs,i),
        .vlx = ListEntry(x,i)    }};
      if(ListEntry(evosys->f[PRESURF0],i))
        troubled |= ListEntry(evosys->f[PRESURF0],i)(node, evv);
      if(ListEntry(evosys->f[PRESURF],i))
        troubled |= ListEntry(evosys->f[PRESURF],i)(node, evv);
      if(ListEntry(evosys->f[PRESURF2],i))
        troubled |= ListEntry(evosys->f[PRESURF2],i)(node, evv);
    }
    node->dat->info->evo_troubled |= troubled;

    /* add time spend on PRESURF */
    loadtimer_stop(node);
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

    /* time SETSRC & VOLRHS */
    loadtimer_start(node);

    troubled = 0;
    forList(u, i)
    {
      tEvoVars evv[1] = {{
        .vlu = ListEntry(u,i),
        .vlr = ListEntry(rhs,i),
        .vlx = ListEntry(x,i)    }};

      /* set all early sources */
      if(ListEntry(evosys->f[SETSRC0],i))
        troubled |= ListEntry(evosys->f[SETSRC0],i)(node, evv);
      /* set all sources */
      if(ListEntry(evosys->f[SETSRC],i))
        troubled |= ListEntry(evosys->f[SETSRC],i)(node, evv);

      /* check if we have to set extra vars for LDG */
      if(ListEntry(evosys->f[XVOLRHS],i))
        troubled |= ListEntry(evosys->f[XVOLRHS],i)(node, evv);
      if(ListEntry(evosys->f[XSURFRHS],i))
      {
        get_all_surfaces(node); //get surfaces of u, don't need surf of x yet
        //FIXME: make get_all_surfaces_vl that does it just for a varlist
        troubled |= ListEntry(evosys->f[XSURFRHS],i)(node, evv);
      }

      /* set all volume RHSs */
      if(ListEntry(evosys->f[VOLRHS],i))
        troubled |= ListEntry(evosys->f[VOLRHS],i)(node, evv);
    }
    node->dat->info->evo_troubled |= troubled;

    /* add time spend on SETSRC & VOLRHS */
    loadtimer_stop(node);
  }

  /* check if there is a single XVOLRHS */
  have_XRHS = 0;
  forList(u, ie)
    if( (ListEntry(evosys->f[XVOLRHS],ie)) ||
        (ListEntry(evosys->f[XSURFRHS],ie)) ) { have_XRHS = 1;  break; }

  /* After we have done all we can without the surface data, we now wait
     until we get all the surface data: */
  /* get surfaces of u so that we can compute fluxes that depend on u */
  MPIexchange_get_all_myln_data(mesh);

  /* We may also need to get the surface data for the extra vars in x.
     Note, the surface data for u has already arrived. */
  if(have_XRHS)
  {
    /* For now we just do the entire surface exchange again. (FIXME) */
    MPIexchange_set_all_myln_localdata(mesh);
    MPIexchange_request_all_myln_data(mesh);
    MPIexchange_get_all_myln_data(mesh);
  }

  /* Now we have all surface data */

  /* loop over all nodes, after MPI data has been received */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    int i, troubled;

    /* time SURFRHS */
    loadtimer_start(node);

    /* add all surface RHSs */
    troubled = 0;
    forList(u, i)
    {
      tEvoVars evv[1] = {{
        .vlu = ListEntry(u,i),
        .vlr = ListEntry(rhs,i),
        .vlx = ListEntry(x,i)    }};
      if(ListEntry(evosys->f[SURFRHS],i))
        troubled |= ListEntry(evosys->f[SURFRHS],i)(node, evv);
    }
    node->dat->info->evo_troubled |= troubled;

    /* add time spend on SURFRHS */
    loadtimer_stop(node);
  }
}

/* parse options for evolve_limiter_mesh */
int evolve_call_limiter(tElm *elm, int opt)
{
  int trbl_score;

  switch(opt)
  {
  case 0: /* always call limiter */
    return 1;
  case 1: /* call limiter only for switched elms */
    trbl_score = elm->dat->info->trbl_score;
    if(trbl_score>=1 || trbl_score<=-NOTROUBLES) return 1;
    else                                         return 0;
  case 2: /* call limiter only for elms that are still unlimited */
    if(elm->dat->info->unlimited) return 1;
    else                          return 0;
  default:
    return 1;
  }
}

/* Apply limiters to evo subsystems.
   If opt=1
     we do it only if trbl_score >= 1 or trbl_score<=-NOTROUBLES
     i.e. only for switched nodes.
   If opt=0
     we ALWAYS do it. */
/* Version for entire mesh: */
void evolve_limiter_mesh(tMesh *mesh, pVLList *u, int opt)
{
  tEvoSys *evosys = mesh->evosys;
  int j;
  tVarList *vl;

  if(!u) return;

  if(PR) PRFs(":\n");

  /* loop over all nodes before MPI requests */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    int i, troubled;

    /* time PRELIM & LIMDATA */
    loadtimer_start(node);

    if(evolve_call_limiter(node, opt))
    {
      troubled = 0;
      forList(u, i)
      {
        tEvoVars evv[1] = {{.vlu=ListEntry(u,i), .vlr=NULL, .vlx=NULL}};
        /* call funcs that we need before limiters */
        if(ListEntry(evosys->f[PRELIM0],i))
          troubled |= ListEntry(evosys->f[PRELIM0],i)(node, evv);
        if(ListEntry(evosys->f[PRELIM],i))
          troubled |= ListEntry(evosys->f[PRELIM],i)(node, evv);

        /* set data limiter needs in myindc arrays of each node */
        if(ListEntry(evosys->f[LIMDATA],i))
          troubled |= ListEntry(evosys->f[LIMDATA],i)(node, evv);
      }
      node->dat->info->evo_troubled |= troubled;
    } /* end if */

    /* add time spend on PRELIM & LIMDATA */
    loadtimer_stop(node);
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
    int i;

    /* time LIMITER */
    loadtimer_start(node);

    if(evolve_call_limiter(node, opt))
    {
      forList(u, i)
      {
        tEvoVars evv[1] = {{.vlu=ListEntry(u,i), .vlr=NULL, .vlx=NULL}};
        /* apply limiter */
        if(ListEntry(evosys->f[LIMITER],i))
        {
          int ret = ListEntry(evosys->f[LIMITER],i)(node, evv);

          /* increase nlim if limiter was active, otherwise reset nlim */
          if(ret) node->dat->info->nlim += 1;
          else    node->dat->info->nlim = 0;

          /* also count a ret!=0 as trouble */
          node->dat->info->evo_troubled |= ret;
        }
      }
    } /* end if */

    /* add time spend on LIMITER */
    loadtimer_stop(node);
  }
}


/* update some vars by calling funcs in PRESURF*, SETSRC*
   often PRESURF does cons2prim, SETSRC sets stress-energy,
   PRELIM sets ADM metric */
void evolve_setsrc_again_nontroubled_nodes_mesh(tMesh *mesh,
                                                pVLList *rhs, pVLList *u)
{
  tEvoSys *evosys = mesh->evosys;
  pVLList *x = evosys->x; /* extra vars needed for LDG */

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
        tEvoVars evv[1] = {{
          .vlu = ListEntry(u,i),
          .vlr = ListEntry(rhs,i),
          .vlx = ListEntry(x,i)    }};

        if(ListEntry(evosys->f[PRESURF0],i))
          ListEntry(evosys->f[PRESURF0],i)(node, evv);

        if(ListEntry(evosys->f[PRESURF],i))
          ListEntry(evosys->f[PRESURF],i)(node, evv);

        if(ListEntry(evosys->f[PRESURF2],i))
          ListEntry(evosys->f[PRESURF2],i)(node, evv);

        if(ListEntry(evosys->f[SETSRC0],i))
          ListEntry(evosys->f[SETSRC0],i)(node, evv);

        if(ListEntry(evosys->f[SETSRC],i))
          ListEntry(evosys->f[SETSRC],i)(node, evv);

        /* we do not need to run PRELIM, since evolve_limiter_mesh will
           run right before evolve_setsrc_mesh is called */
        /* if(ListEntry(evosys->f[PRELIM0],i))
             ListEntry(evosys->f[PRELIM0],i)(node, evv);
           if(ListEntry(evosys->f[PRELIM],i))
             ListEntry(evosys->f[PRELIM],i)(node, evv); */
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
    int   dnf = Geti(Par("evolve_filter_dn"));

    if(PR) { PRF;printf(": filtering all evolution vars in evosys->u\n"); }

    /* loop over list of varlists and filter each varlist */
    forList(u, i)
    {
      tVarList *vl = ListEntry(u,i);
      expfilter_mesh_vl1(vl, af, sf, dnf);
    }
  }

  if(filter_varlist)
  {
    double af = Getd(Par("evolve_filter_alp"));
    double sf = Getd(Par("evolve_filter_s"));
    int   dnf = Geti(Par("evolve_filter_dn"));
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

    if(vl->n) expfilter_mesh_vl1(vl, af, sf, dnf);

    vlfree(vl);
    free(list);
  }

  return 0;
}


/* funcs to time the EVOLVE time bin, and also to output times */

/* start special EVOLVE_time timer */
int EVOLVE_timer_start(tMesh *mesh)
{
  if(Getb(Par("evolve_output_timers"))) timer_start("EVOLVE-time", 1);
  return 0;
}
/* stop special EVOLVE_time timer */
int EVOLVE_timer_stop(tMesh *mesh)
{
  if(Getb(Par("evolve_output_timers"))) timer_stop("EVOLVE-time", 1);
  return 0;
}

/* output EVOLVE timer and the loadtimers */
int evolve_output_timers(tMesh *mesh)
{
  if(Getb(Par("evolve_output_timers")))
  {
    static int firstcall = 1;
    char *outdir = Gets(Par("outdir"));
    char f[100], s[1000];
    FILE *fp;
    double EVOLVE_time = timer_get_dtime("EVOLVE-time", 1);
    double rank_loadtime;

    /* open file */
    snprintf(f, 100, "%%s/evolve_timers.%%0%dd", (int) log10(nMPI_size())+1);
    snprintf(s, 1000, f, outdir, nMPI_rank());
    fp = fopen(s, "a");
    if(!fp) errorexits("could not open %s", s);
    if(firstcall)
    {
      fprintf(fp, "#    PhysTime     EVOLVE-time   rank_loadtime"
                  "   elm0_loadtime   elm1_loadtime   ...\n");
      firstcall = 0;
    }

    /* calc total load_Time on this rank */
    rank_loadtime = 0.;
    formyelms(mesh)
    {
      tElm *elm = MyElm;
      tDat *dat = elm->dat;
      if(dat)
      {
        tNodeInfo *info = dat->info;
        rank_loadtime += info->load_TimeIn_s;
      }
    }
    /* ouput timer data */
    fprintf(fp, "%13g  %13gs  %13gs", mesh->time, EVOLVE_time, rank_loadtime);
    formyelms(mesh)
    {
      tElm *elm = MyElm;
      tDat *dat = elm->dat;
      if(dat)
        fprintf(fp, "  %13gs", dat->info->load_TimeIn_s);
    }
    fprintf(fp, "\n");
    fclose(fp);
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
  {
    tEvoVars evv[1] = {{.vlu=ListEntry(u,i), .vlr=ListEntry(rhs,i)}};
    if(ListEntry(evosys->f[SETSRC],i))
      ListEntry(evosys->f[SETSRC],i)(node, evv);
  }

  /* set all vol. RHSs */
  forList(u, i)
  {
    tEvoVars evv[1] = {{.vlu=ListEntry(u,i), .vlr=ListEntry(rhs,i)}};
    if(ListEntry(evosys->f[VOLRHS],i))
      ListEntry(evosys->f[VOLRHS],i)(node, evv);
  }

  //Test:  get_all_surfaces(node);

  /* add all surf. RHSs */
  forList(u, i)
  {
    tEvoVars evv[1] = {{.vlu=ListEntry(u,i), .vlr=ListEntry(rhs,i)}};
    if(ListEntry(evosys->f[SURFRHS],i))
      ListEntry(evosys->f[SURFRHS],i)(node, evv);
  }

  /* do not free all surface info, because we call evolve_setrhs repeatedly */
  //if(0) evolve_free_surfaces(node, u);
}
