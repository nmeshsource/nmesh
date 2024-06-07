/* evolve.c */
/* Wolfgang Tichy, 5/2022 */

#include "nmesh.h"
#include "evolve.h"

#define PR 0


/* use amr vars */
extern tAMR amr[1];
extern tDGglobals DGglobals[1];


/* Determine and set trouble score in a node,
   i.e. set node->dat->info->trbl_score .
   Note evosys->f[TROUBLE] returns a score of 1, 0, or -1:
    1 trouble  =>  switch to more robust method (e.g. fv)
    0 node ok with current method  =>  do nothing
   -1 all ok   =>  switch to more accurate but delicate method (e.g. dg)
   examples: if shock in dg return 1
             if shock in fv return 0
             if no shock in dg return 0
             if no shock in fv return -1
   Here we also accumulate the 1 and -1 scores and return them */
int evolve_set_trouble_score(tNode *node)
{
  tMesh *mesh = node->pat->mesh;
  tEvoSys *evosys = mesh->evosys;
  pVLList *u_p = evosys->u_p;
  pVLList *u   = evosys->u;
  pVLList *r   = evosys->rhs;
  int max_trouble = 1073741824;
  int troubled = 0; /* default is to assume that we change nothing */
  int tr_max = -max_trouble; /* init to low value */
  int i;

  if(node->dat == NULL) errorexit("node->dat is NULL");

  //if(PR) PRFs(":\n");

  /* check all evo systems for trouble and put max into tr_max */
  forList(u, i)
  {
    tVarList *vlu_p = ListEntry(u_p,i);
    tVarList *vlu   = ListEntry(u,i);
    tVarList *vlr   = ListEntry(r,i);

    /* run TROUBLE func */
    if(tr_max<=0) /* need to check only if there no trouble yet */
    {
      if(ListEntry(evosys->f[TROUBLE],i))
      {
        int tr = 0;
        tr = ListEntry(evosys->f[TROUBLE],i)(node, vlr, vlu, vlu_p);
        if(tr>tr_max) tr_max = tr;
      }
    }
  }

  /* set troubled flag to 1, 0, or -1 */
  if(tr_max>0) /* new trouble found, switch to fv */
    troubled = 1;
  else if(tr_max==0 || tr_max==-max_trouble) /* ok, keep as is */
    troubled = 0;
  else /* all is very good */
    troubled = -1;
  //pr_nodename(node);
  //printf(" tr_max=%d troubled=%d", tr_max, troubled);

  /* set node->dat->info->trbl_score */
  if(troubled>0) /* i.e. there is trouble now */
  {
    /* if there was no trouble, we set the trouble score to 1 */
    if(node->dat->info->trbl_score<=0) node->dat->info->trbl_score  = 1;
    /* if there was trouble before, we continuously increase the score */
    else                            node->dat->info->trbl_score += 1;
  }
  else if(troubled==0) /* is ok, keep node as is */
  {
    node->dat->info->trbl_score = 0;
  }
  else /* all is good, can switch back to dg */
  {
    /* if there was trouble, we set the trouble score to -1 */
    if(node->dat->info->trbl_score>=0) node->dat->info->trbl_score  = -1;
    /* if there was no trouble, we continuously lower the score */
    else                            node->dat->info->trbl_score -= 1;
  }

  /* make sure trouble score does not become too large or small */
  if(node->dat->info->trbl_score < -max_trouble)
    node->dat->info->trbl_score = -max_trouble;
  if(node->dat->info->trbl_score > max_trouble)
    node->dat->info->trbl_score = max_trouble;

  //printf(" ->trbl_score=%d\n", node->dat->info->trbl_score);

  return node->dat->info->trbl_score;
}

/* determine and set trouble score in each node,
   i.e. set node->dat->info->trbl_score */
int evolve_set_trouble_score_mesh(tMesh *mesh)
{
  tEvoSys *evosys = mesh->evosys;
  pVLList *u_p = evosys->u_p;
  int Max_trb, max_trb=INT_MIN;
  int Min_trb, min_trb=INT_MAX;
  if(PR) PRFs(":\n");

  /* collect min/max values of all nodes and their of neighbors */
  evolve_collect_u_p_data_mesh(mesh, u_p);

  /* loop over all nodes, check for trouble, and node-info trouble score */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    evolve_set_trouble_score(node);
  }

  /* now find rank-local max of node->dat->info->trbl_score */
  formylnodes_noomp(mesh)
  {
    tNode *node = MyLnode;
    int trb = node->dat->info->trbl_score;
    if(trb>max_trb) max_trb = trb; /* max trouble */
    if(trb<min_trb) min_trb = trb; /* min trouble */
  }

  /* Max over all ranks */
  Max_trb = max_trb;
  MCK( nMPI_Allreduce(&max_trb, &Max_trb, 1, nMPI_INT, nMPI_MAX) );

  /* Min over all ranks */
  Min_trb = min_trb;
  MCK( nMPI_Allreduce(&min_trb, &Min_trb, 1, nMPI_INT, nMPI_MIN) );

  /* if there is bad trouble somewhere */
  if(Max_trb>0)
    return Max_trb; /* returns max. trouble value of all nodes */
  else
    return Min_trb; /* no trouble, return smallest to signal need for dg */
}


/* switch from dg to fv based on node->dat->info->trbl_score flag */
void evolve_switch_troubled_nodes_mesh(tMesh *mesh)
{
  tRef ref[1]; /* for ref info */
  int firstit;
  int ptUNI[] = { P_UNIFORM, P_UNIFORM, P_UNIFORM };
  int ptLGL[] = { P_LGL, P_LGL, P_LGL };
  if(PR) PRFs(":\n");

  /* free surfaces & indc since they will change now anyway */
  evolve_free_communication_structs(mesh);

  /* Set ref to uniform. All nodes with trbl_score>0 should have trbl_ref
     set to uniform already. */
  /* loop over all nodes, and flag all troubled nodes for uniform grid */
  firstit = 1;
  ref->method = REF_METH_DONOTHING;
  formylnodes_noomp(mesh)
  {
    tNode *node = MyLnode;

    /* make sure time on node is correct */
    node->time = mesh->time;

    /* mark troubled nodes as to be refined */
    if(node->dat->info->trbl_score > 0)
    {
      tRef *info_ref = node->dat->info->trbl_ref;
      //if(elmname_is(node, "0_700")) {PRFs(": trbl_ref: ");printref(info_ref);}

      if(!firstit)
        if(info_ref->method != ref->method)
          errorexit("nodes with trbl_score > 0 must all have same trbl_ref");
      firstit = 0;

      /* set ref and node->rflag, since ref is global don't use OpenMP! */
      ref[0] = info_ref[0];
      node->rflag = ref->method;
    }
    else
    {
      node->rflag = 0;
    }
  }

  /* get ref->method from my proc to all others */
  refine_synchronize_ref_method(ref);

  /* do p-refinement to desired n and point type */
  prefine_nodes_if_rflag(mesh, ref);
  /* now some aux vars (and others) are not set */
  /* this will be fixed by evolve_setsrc_again_nontroubled_nodes_mesh */

  /* clear rflag on all leaf nodes */
  refine_set_rflag_forall_nodes(mesh, 0);

  /* switch on fv on all uniform nodes */
  refine_set_use_fv_if_pt_typ(mesh, ptUNI, 1);

  /* switch off fv on all LGL nodes */
  refine_set_use_fv_if_pt_typ(mesh, ptLGL, 0);

  /* now that nodes are changed re-init surfaces & indc */
  evolve_init_communication_structs(mesh);

  /* balance load, now that some nodes use a different method */
  //FIXME: do something better than simple_load_balance
}

/* switch from fv to dg based on node->dat->info->trbl_score flag */
void evolve_switch_nontroubled_nodes_mesh(tMesh *mesh)
{
  tRef ref[1]; /* for ref info */
  int firstit, order_sav, force_sav;
  int ptUNI[] = { P_UNIFORM, P_UNIFORM, P_UNIFORM };
  int ptLGL[] = { P_LGL, P_LGL, P_LGL };
  tVarList *vl_extrap;

  if(PR) PRFs(":\n");

  /* free surfaces & indc since they will change now anyway */
  evolve_free_communication_structs(mesh);

  /* Set ref to LGL. All nodes with trouble<0 should have trbl_ref
     set to LGL already. */
  /* loop over all nodes, and flag all non-troubled nodes for LGL grid */
  firstit = 1;
  ref->method = REF_METH_DONOTHING;
  formylnodes_noomp(mesh)
  {
    tNode *node = MyLnode;

    /* mark non-troubled nodes as to be refined */
    if(node->dat->info->trbl_score <= -NOTROUBLES)
    {
      tRef *info_ref = node->dat->info->trbl_ref;
      //if(elmname_is(node, "0_700")) {PRFs(": trbl_ref: ");printref(info_ref);}

      if(!firstit)
        if(info_ref->method != ref->method)
          errorexit("nodes with trbl_score<=-NOTROUBLES must all have same trbl_ref");
      firstit = 0;

      /* set ref and node->rflag, since ref is global don't use OpenMP! */
      ref[0] = info_ref[0];
      node->rflag = ref->method;
    }
    else
    {
      node->rflag = 0;
    }
  }

  /* get ref->method from my proc to all others */
  refine_synchronize_ref_method(ref);

  /* set varlist where we use extrap to face before interp to dg */
  if(DGglobals->fv2d_interp_use_extrap1)
    vl_extrap = DGglobals->fv2d_interp_use_extrap1_vl;
  else
    vl_extrap = NULL;

  /* do p-refinement to desired n and point type: */
  /* We have just determined that these fv nodes are so untroubled that we
     want to switch to them to dg. I.e. there will be no shocks. Thus we
     will use high order Lagrange interpolation (and not WENO6).
  */
  /* 1. select interpolation scheme and order */
  force_sav = amr->force_interp_scheme;         //backup flag
  order_sav = Geti(amr->Lagrange_interp_order); //backup order par
  amr->force_interp_scheme = INTERP_LAGRANGE;   //set internal amr flag
  Seti(amr->Lagrange_interp_order, 12);         //set 12th order Lag interp
  /* 2. now call p-refinement from amr */
  rec1d_fv_uface_to_uin_1_if_rflag(mesh, vl_extrap, 0); //extrap to face
  prefine_nodes_if_rflag(mesh, ref);
  /* now some aux vars (and others) are not set */
  /* this will be fixed by evolve_setsrc_again_nontroubled_nodes_mesh */
  /* 3. restore par amr_Lagrange_interp_order and amr->force_interp_scheme */
  amr->force_interp_scheme = force_sav;
  Seti(amr->Lagrange_interp_order, order_sav);

  /* clear rflag on all leaf nodes */
  refine_set_rflag_forall_nodes(mesh, 0);

  /* switch on fv on all uniform nodes */
  refine_set_use_fv_if_pt_typ(mesh, ptUNI, 1);

  /* switch off fv on all LGL nodes */
  refine_set_use_fv_if_pt_typ(mesh, ptLGL, 0);

  /* now that nodes are changed re-init surfaces & indc */
  evolve_init_communication_structs(mesh);

  /* balance load, now that some nodes use a different method */
  //FIXME: do something better than simple_load_balance
}

/* set u = u_p, and then switch to fv */
void evolve_prepare_do_over_mesh(tMesh *mesh)
{
  tEvoSys *evosys = mesh->evosys;
  double  t = mesh->time;
  double dt = mesh->dt;
  pVLList *u_p = evosys->u_p;
  pVLList *u   = evosys->u;

  if(PR) PRFs(":\n");

  /* go back to t-dt and set u = u_p */
  copy_pVLList(u, u_p, vlcopy,0);
  mesh->time = t-dt;

  /* loop over all nodes, and set time */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    node->time = mesh->time;
  }

  /* switch nodes based on trouble flag */
  evolve_switch_troubled_nodes_mesh(mesh);
}

/* Set myindc for u_p to get min/max of u_p needed for a RDMP trouble indicator.
   Here we use limdata_MRS to get min,max,average. */
void evolve_collect_u_p_data_mesh(tMesh *mesh, pVLList *u_p)
{
  tEvoSys *evosys = mesh->evosys;
  int i;

  if(PR) PRFs(":\n");

  /* loop over list of varlists */
  forList(u_p, i)
  {
    tVarList *vl = ListEntry(u_p,i);

    /* set limiter data in indicators (indc) if trouble indicator needs it */
    if(ListEntry(evosys->f[TROUBLE],i))
    {
      /* set data RDMP indicator needs in myindc arrays of each node */
      formylnodes(mesh)
      {
        tNode *node = MyLnode;
        limdata_MRS(node, vl); /* sets min,max,average */
      }

      /* initiate indc exchange */
      request_all_myln_indc_exchange_for_vl(mesh, vl);
      /* After this we could do work. MPI is now busy sending buffers */

      /* now get the indicators and wait for MPI buffers if necessary */
      get_all_myln_indc_for_vl(mesh, vl);
    }
  } /* end forList */
}

/* Relaxed discrete maximum principle (RDMP) indicator as in 2109.11645
   This indicator is originally from 1406.7416 (Eqs 23-25), where a
   discrete maximum principle (DMP) that is relaxed by delta is introduced.
   We implement it as in 2109.11645 (Eqs 65-66). */
int evolve_RDMP_trouble(tNode *node, tVarList *vlu, tVarList *vlu_p,
                        double deltafac, double delta0, double epsilon)
{
  //tMesh *mesh = node->pat->mesh;
  tDat *dat;
  int nvars = VLn(vlu);
  int vli, f, ni;
  double min_u_p[nvars], max_u_p[nvars];
  double delta;
  double lower_lim, upper_lim;
  int troubled;

  dat = node->dat;
  if(!dat) return 0;

  /* find abs max and min of u_p and compare with u */
  troubled = 0;
  forvl(vlu_p, vli)
  {
    int iu = Vind(vlu, vli);
    double *u = Vard(node, iu);
    int ijk;
    int iu_p = Vind(vlu_p, vli);
    double min = dat->ic[iu_p]->myindc->d[1];    /* get my min(u_p) */
    double max = dat->ic[iu_p]->myindc->d[2];    /* get my max(u_p) */
    //double wbar = dat->ic[iu_p]->myindc->d[0];   /* get my average */

    /* find min and max of u_p in neighbors + this node */
    max_u_p[vli] = max;
    min_u_p[vli] = min;
    for(f=0; f<6; f++)
      for(ni=0; ni<node->nfnb[f]; ni++)
      {
        double min_nb = dat->ic[iu_p]->nbindc[f][ni]->d[1];
        double max_nb = dat->ic[iu_p]->nbindc[f][ni]->d[2];
        if(min_nb < min_u_p[vli]) min_u_p[vli] = min_nb;
        if(max_nb > max_u_p[vli]) max_u_p[vli] = max_nb;
      }

    /* set delta */
    delta = deltafac * max2(delta0, epsilon*(max_u_p[vli] - min_u_p[vli]));

    /* set lower and upper allowed limits for u */
    lower_lim = min_u_p[vli] - delta;
    upper_lim = max_u_p[vli] + delta;

    /* check if u is within range */
    forpoints(node, ijk)
      if(u[ijk] < lower_lim || u[ijk] > upper_lim)
      {
        troubled = 1; /* u is troubled at one point */
        break;
      }
    if(troubled) /* if this u is troubled we don't need to check others */
      break;
  }

  /* return whether there is trouble */
  return troubled;
}

/* Compute Persson trouble indicator for array u, with point type pt_typ,
   using only the on first ncoeffs coeffs.
   Ref: Per-Olof Persson and Jaime Peraire. Sub-cell shock capturing for
   discontinuous Galerkin methods.
   In 44th AIAA Aerospace Sciences Meeting and Exhibit.
   American Institute of Aeronautics and Astronautics, Inc., 2006. */
int evolve_Persson_array_trouble(tArray *u, double u_scale, int pt_typ[3],
                                 int ncoeffs[3], double alpha)
{
  tArray *At[3];
  tArray *ca;
  int i,j,k, n_max;
  double c2_sum, c2_hi, se, se_lim;
  int troubled;

  /* The Persson's indicator se calculated below is always negative.
     So if alpha<0 se_lim=-alpha*log10(n_max)>0. Thus se>=se_lim
     will never be true, and we always return troubled=0 for alpha<0. */
  if(alpha < 0.) return 0;

  /* get ana. matrices */
  At3_pt_typ_n(pt_typ, Arrn(u), At);

  /* get coeffs of var iu */
  ca = alloc_array(Arrn(u));
  basis_array_analysis3_At(At, u, ca);

  /* set the 0th coeff (that corresponds to the node average) to u_scale,
     so that we get the same result, no matter how high the node average */
  if(u_scale >= 0.) Arrd(ca)[0] = u_scale * 2.*sqrt(2.);

  /* compute sum of squares of the ncoeffs coeffs */
  c2_sum = 0.;
  forijk(i,j,k, ncoeffs)
  {
    int ijk = Ind_n(i,j,k, Arrn(u));
    double co = Arrd(ca)[ijk];
    c2_sum += co*co;
  }

  /* compute sum of squares of highest coeffs */
  c2_hi = 0.;
  forijk(i,j,k, ncoeffs)
    if( ((i==ncoeffs[0]-1) && (ncoeffs[0]>1)) ||
        ((j==ncoeffs[1]-1) && (ncoeffs[1]>1)) ||
        ((k==ncoeffs[2]-1) && (ncoeffs[2]>1)) )
    {
      int ijk = Ind_n(i,j,k, Arrn(u));
      double co = Arrd(ca)[ijk];
      c2_hi += co*co;
    }

  /* Persson's indicator */
  se = log10( c2_hi / (c2_sum + DBL_MIN) + LOGARGFLOOR );

  /* find max of number of points in all 3 dirs for ncoeffs coeffs */
  n_max = max3(ncoeffs[0], ncoeffs[1], ncoeffs[2]);

  /* 2109.11645 says that in 1D there is trouble if
     se >= -alpha * log10(n)  inside a node
     for a dg node alpha=4, for a fv node alpha=5 */
  se_lim = -alpha * log10(n_max);

  //PRFs(": ");pr_nodename(node);printf(" iu=%d: ", iu);
  //printf("c2_nm1_sum=%g ", c2_nm1_sum);
  //printf("c2_hi=%g c2_sum=%g  ", c2_hi, c2_sum);
  //printf(" %g %g\n", se, se_lim);

  /* set troubled flag */
  if(se >= se_lim)
  {
    //char ns[100];

    //PRFs(": ");printf("%s iu=%d: ", nodename(node,ns,99), iu);
    //printf(" %g %g\n", se, se_lim);
    troubled = 1;
  }
  else
  {
    troubled = 0;
  }

  if(PR)
  {
    printf(" c2_hi=%g c2_sum=%g  ", c2_hi, c2_sum);
    printf(" %g %g => %d", se, se_lim, troubled);
    printf("\n");
  }

  free_array(ca);
  return troubled;
}

/* Compute Persson trouble for var iu in node based on first ncoeffs
   coeffs */
int evolve_Persson_trouble_ncoeffs(tNode *node, int iu, double u_scale,
                                   int ncoeffs[3],
                                   double alpha, double alpha_fv)
{
  //tMesh *mesh = node->pat->mesh;
  int fv = node->dat->info->use_fv;
  double alpha_a;

  /* 2109.11645 says that in 1D there is trouble if
     se >=     -alpha * log10(n)  inside a dg node
     se >= -(alpha+1) * log10(n)  inside a fv node, in both cases alpha=4 */
  if(fv) alpha_a = alpha_fv;
  else   alpha_a = alpha;

  if(PR) //(nodename_is(node, "0_701"))
  {
    pr_nodename(node);
    printf(" %d", fv);
  }

  return evolve_Persson_array_trouble(VarA(node,iu), u_scale, node->pt_typ,
                                      ncoeffs, alpha_a);
}

/* Compute Persson trouble for var iu in node on dg-grid based on first
   ncoeffs coeffs */
int evolve_Persson_trouble_ncoeffs_dg(tNode *node, int iu, double u_scale,
                                      int ncoeffs[3],
                                      double alpha, double alpha_fv)
{
  int fv = node->dat->info->use_fv;
  int troubled;

  if(PR) //(nodename_is(node, "0_701"))
  {
    pr_nodename(node);
    printf(" %d", fv);
  }

  if(fv)
  {
    /* use trbl_ref to find n_dg and pt_typ_dg*/
    tRef *ref = node->dat->info->trbl_ref;
    int n_dg[3], pt_typ_dg[3];
    tArray *u_interp;
    int npts, fv_extrap;

    /* since evolve_Persson_array_trouble returns 0 for alpha_fv<0 anyway,
       we do this here already */
    if(alpha_fv < 0.) return 0;

    /* decide if we extrap to faces 1st before we call interp_to_pt_typ */
    fv_extrap = 0;
    if(DGglobals->fv2d_interp_use_extrap1)
      if( vlindex(DGglobals->fv2d_interp_use_extrap1_vl, iu) >= 0 )
        fv_extrap = 1; // but always 0 could also be an option...

    /* get num. and type of points on dg grid that we would switch to */
    hp_refine_set_n_pt_typ(node, ref, n_dg, pt_typ_dg);

    /* Pick a good value for npts:
       Since We want to test for smoothness in here, we use LAGRANGE as
       it is very sensitive to shocks. Yet in order to not incur too large
       of a numerical error (when interpolating on a uniform grid), we
       need to limit the number of interpolation points npts. */
    npts = 12; // should this be a GRHD par that's passed in???

    /* interpolate u to dg grid */
    u_interp = alloc_array(n_dg);
    if(fv_extrap) rec1d_uface_to_uin_1_var(node, iu, 0); //extrap back to face
    interp_to_pt_typ(node, iu, pt_typ_dg, npts,INTERP_LAGRANGE,1., u_interp);
    if(fv_extrap) rec1d_uface_to_uin_1_var(node, iu, 1); //undo extrap
    troubled = evolve_Persson_array_trouble(u_interp, u_scale, pt_typ_dg,
                                            ncoeffs, alpha_fv);
    free_array(u_interp);
  }
  else
  {
    troubled = evolve_Persson_array_trouble(VarA(node,iu), u_scale,
                                            node->pt_typ, ncoeffs, alpha);
  }
  return troubled;
}

/* DO NOT USE:
   compute Persson trouble based on all coeffs in node */
int evolve_Persson_trouble(tNode *node, int iu, double u_scale,
                           double alpha, double alpha_fv)
{
  int *ncoeffs = node->n;
  return evolve_Persson_trouble_ncoeffs(node, iu, u_scale, ncoeffs,
                                        alpha, alpha_fv);
}

/* set trouble score ts based on whether node is troubled, and dg or fv */
int trouble_score(tNode *node, int troubled)
{
  int fv = node->dat->info->use_fv;

  if(troubled)
  {
    if(fv) return  0; /* keep fv */
    else   return  1; /* switch to fv */
  }
  else
  {
    if(fv) return -1; /* switch to dg */
    else   return  0; /* keep dg */
  }
}

/* print debug info for the node named "nname" */
void trouble_print_if_name(tNode *node, const char *nname,
                           int trbl, const char *text)
{
  if(nodename_is(node, nname))
  {
    int fv = node->dat->info->use_fv;
    PRF;printf(": %s fv=%d, %s trbl=%d => score=%d\n", nname, fv, text,
               trbl, trouble_score(node, trbl));
  }
}

/* reset evo_troubled flag on all nodes */
int trouble_reset_evo_troubled_mesh(tMesh *mesh)
{
  formylnodes_noomp(mesh)
  {
    tNode *node = MyLnode;
    node->dat->info->evo_troubled = 0;
  }
  return 0;
}

/* return flag that tells us whether there was any trouble during an evo
   step i.e. during any RK substep */
int evolve_evosteps_troubled(tNode *node)
{
  return node->dat->info->evo_troubled;
}
