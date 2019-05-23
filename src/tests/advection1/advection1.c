/* advection1.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "advection1.h"

#define PR 1


/* func pointer for numerical flux */
void (*advection1_numflux)(tDGinfo *d);

/* frequently used pars */
double advection1_direction[3];

/* func to init frequently used pars */
int advection1_init_global_pars(tMesh *mesh)
{
  char *advdir = Gets(Par("advection1_direction"));
  sscanf(advdir, "%lg %lg %lg",
         &(advection1_direction[0]),
         &(advection1_direction[1]),
         &(advection1_direction[2]));
  printf("setting: advection1_direction[0] = %g\n", advection1_direction[0]);
  printf("setting: advection1_direction[1] = %g\n", advection1_direction[1]);
  printf("setting: advection1_direction[2] = %g\n", advection1_direction[2]);
  return 0;
}


/* flux in direction norm */
void advection1_flux1d(tMesh *mesh, int ncons, double *f, double norm[3],
                       double *u)
{
  double *nd = advection1_direction;

  /* flux times norm */
  f[0] = (norm[0]*nd[0] + norm[1]*nd[1] + norm[2]*nd[2]) * u[0];
}

/* eigenvalue in direction norm */
void advection1_eigenval1d(tMesh *mesh, int ncons, double *lam, double norm[3])
{
  double *nd = advection1_direction;

  /* eigenvalue */
  lam[0] = (norm[0]*nd[0] + norm[1]*nd[1] + norm[2]*nd[2]);
}

/* flux and its derivs for adv. eqn: f^i = n^i u */
void advection1_f_df(tMesh *mesh, tVarList *vlu)
{
  int iu = vlu->index[0];
  int ifx  = Ind("advection1_fx");
  int ify = ifx+1;
  int ifz = ifx+2;
  int ifxx = Ind("advection1_fxx");
  int ifyx = ifxx+3;
  int ifzx = ifxx+6;

  /* compute derivs */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    double *u = Vard(node, iu);
    double *fx = Vard(node, ifx);
    double *fy = Vard(node, ifx+1);
    double *fz = Vard(node, ifx+2);
    int i;

    /* flux at each point */
    forpoints(node, i)
    {
      double u_i = u[i];
      double no[3] = { 1., 0., 0. };
      advection1_flux1d(mesh,1, &(fx[i]),no, &u_i);
      no[0] = 0;
      no[1] = 1.;
      advection1_flux1d(mesh,1, &(fy[i]),no, &u_i);
      no[1] = 0.;
      no[2] = 1.;
      advection1_flux1d(mesh,1, &(fz[i]),no, &u_i);
    }

    /* flux derivs */
    cart_partials(node, ifx, ifxx);
    cart_partials(node, ify, ifyx);
    cart_partials(node, ifz, ifzx);
  }
}

/* function that sets fluxes and eigenvals on both sides of a node surface.
   In: vlu. Out: ui,ua, fi,fa, lami,lama */
void advection1_fluxes_pt(tDGinfo *d)
{
  tVarList *vlu = d->vlu;
  tNode *node = d->node;
  tMesh *mesh = node->pat->mesh;
  int *n = node->n;
  int f = d->face;
  int dir = f/2;
  int ijk = Ind_n(d->i,d->j,d->k, n);
  int JK = Ind_n_norm(d->i,d->j,d->k, n, dir);
  int nvars = vlu->n;
  double norm[3];
  int l;

  forvl(vlu, l)
  {
    int vi = Vind(vlu,l);
    double *u = Vard_(node, vi);
    double *uaj = Varaj(node, vi, f);

    /* cons var inside node, and cons var on adjacent side */
    d->ui[l] = u[ijk];
    if(uaj)
      d->ua[l] = uaj[JK];
    else /* do something special on outer boundary */
      d->ua[l] = 0.;
  }

  /* get face normal at point ijk */
  node_normal_at_ijk(node, f, ijk, norm);

  /* eigenval in dir norm */
  advection1_eigenval1d(mesh,nvars, d->lami,norm);
  d->lama[0] = d->lami[0]; // eigenval is same on both sides

  /* get inner and adjacent fluxes fi, fa */
  advection1_flux1d(mesh,nvars, d->fi,norm, d->ui);
  advection1_flux1d(mesh,nvars, d->fa,norm, d->ua);
}


/* set a BC on patch boundary */
void advection1_u_BC(tMesh *mesh, tVarList *vlr, tVarList *vlu)
{
  int ir = vlr->index[0];
  int iu = vlu->index[0];
  int ix = Ind("x");
  int sin_profile    = Getv(Par("advection1_profile"),"sin");
  int square_profile = Getv(Par("advection1_profile"),"square");
  char *advdir = Gets(Par("advection1_direction"));
  double nx,ny,nz, nmag2;

  /* prop. dir.*/
  sscanf(advdir, "%lg %lg %lg", &nx, &ny, &nz);
  nmag2 = (nx*nx + ny*ny + nz*nz);

  /* compute boundary flux terms */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    tPat *pat = node->pat;
    int *n = node->n;
    double *r = Vard(node, ir);
    double *u = Vard(node, iu);
    double *x = Vard(node, ix);
    double *y = Vard(node, ix+1);
    double *z = Vard(node, ix+2);
    double t = node->time;
    double norm[3];
    int face, dir, p, i,j,k, ijk;

    /* go over each face */
    for(face=0; face<6; face++)
    {
      tBface *bfaces = pat->bfaces[face];
      dir = face/2;
      p = (face%2)*(n[dir] - 1);

      if(node->patface[face] && bfaces && bfaces->outerbound)
        forplaneN(dir, i,j,k, n, p)
        {
          ijk = Ind_n(i,j,k, n);
          node_normal_at_ijk(node, face, ijk, norm);

          /* if stuff is coming in */
          if(norm[0]*nx + norm[1]*ny + norm[2]*nz < 0.)
          {
            if(sin_profile)
              r[ijk] = -nmag2*cos(nx*x[ijk] + ny*y[ijk] + nz*z[ijk] - nmag2*t);
            if(square_profile)
              r[ijk] = -u[ijk];
          }
        }
    }
  }
}

/* RHS of: d_t u = - d_i f^i */
int advection1_rhs_u(tMesh *mesh, tVarList *vlr, tVarList *vlu)
{
  int ir = vlr->index[0];
  //int iu = vlu->index[0];
  int ifxx = Ind("advection1_fxx");

  TIMER_START;

  /* compute flux */
  advection1_f_df(mesh, vlu);

  /* RHS */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    double *r  = Vard(node, ir);
    double *fxx = Vard(node, ifxx);
    double *fyy = Vard(node, ifxx+4);
    double *fzz = Vard(node, ifxx+8);
    int i;

    /* RHS at each point */
    forpoints(node, i) r[i] = -(fxx[i] + fyy[i] + fzz[i]);
  }

  /* get flux terms on boundary */
  dg_add_surface_fluxes(mesh, vlr, vlu,
                        advection1_fluxes_pt, advection1_numflux);

  /* impose outer BC */
  advection1_u_BC(mesh, vlr, vlu);

  TIMER_STOP;
  return 0;
}


/* set profile in var with index iu */
int advection1_set_profile(tMesh *mesh, int iu)
{
  int ix =  Ind("x");
  int sin_profile    = Getv(Par("advection1_profile"),"sin");
  int square_profile = Getv(Par("advection1_profile"),"square");
  char *advdir = Gets(Par("advection1_direction"));
  double nx,ny,nz, nmag2;

  /* prop. dir.*/
  sscanf(advdir, "%lg %lg %lg", &nx, &ny, &nz);
  nmag2 = (nx*nx + ny*ny + nz*nz);

  //if(PR) PRFs("\n");

  /* profile */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    double *u = Vard(node, iu);
    double *x = Vard(node, ix);
    double *y = Vard(node, ix+1);
    double *z = Vard(node, ix+2);
    double t = mesh->time;
    int i;

    /* profile */
    if(sin_profile)
      forpoints(node, i)
      {
        u[i] = sin(nx*x[i] + ny*y[i] + nz*z[i] - nmag2*t);
      }
    if(square_profile)
      forpoints(node, i)
      {
        double inx, iny;
        if(x[i]>=(-0.7 + nx*t) && x[i]<=(-0.3 + nx*t)) inx = 1.;
        else                                           inx = 0.;
        if(y[i]>=(-0.2 + ny*t) && y[i]<=(+0.2 + ny*t)) iny = 1.;
        else                                           iny = 0.;
        u[i] = inx*iny;
      }
  }
  return 0;
}


/* initialize advection1 */
int advection1_init(tMesh *mesh)
{
  int iu  = Ind("advection1_u");
  int ifx = Ind("advection1_fx");
  int ifxx = Ind("advection1_fxx");
  int iue = Ind("advection1_u_err");
  tVarList *vlu = vlalloc(mesh);
  int numflux = Par("advection1_numflux");
  int limiter = Par("advection1_limiter");
  char *advdir = Gets(Par("advection1_direction"));
  double nx,ny,nz;

  /* prop. dir.*/
  sscanf(advdir, "%lg %lg %lg", &nx, &ny, &nz);

  PRF;printf(": dt = %g\n", mesh->dt);

  /* varlist */
  vlpush(vlu, iu);

  /* enable all needed vars */
  enablevar(mesh, iu);
  enablevar(mesh, ifx);
  enablevar(mesh, ifxx);
  enablevar(mesh, iue);

  /* set initial profile, e.g. at t=0: set u=sin(x) */
  advection1_set_profile(mesh, iu);

  /* register u and its RHS with evolve */
  if(Getv(limiter, "MRS"))
    evolve_register_subsys_u_rhs_src_lim(mesh, vlu, advection1_rhs_u, 0,
                                         limdata_MRS, limiter_MRS);
  else if(Getv(limiter, "minmodB"))
    evolve_register_subsys_u_rhs_src_lim(mesh, vlu, advection1_rhs_u, 0,
                                         limdata_c000_100_010_001,
                                         limiter_minmodB);
  else
    evolve_register_subsys_u_rhs_src_lim(mesh, vlu, advection1_rhs_u, 0, 0,0);

  evolve_print_evosys(mesh);

  /* choose numerical flux */
  if(Getv(numflux, "LLF"))
    advection1_numflux = numflux1d_LLF;
  else
    advection1_numflux = numflux1d_upwind;

  /* test */
  /*
  init_all_myln_myindc_for_vl(mesh, vlu, 3);
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    limdata_MRS(node, vlu);
  }
  request_all_myln_indc_exchange_for_vl(mesh, vlu);

  get_all_myln_indc_for_vl(mesh, vlu);
  formylnodes(mesh)
  {
    //tNode *node = MyLnode;
    //limiter_MRS(node, vlu);
    //limiter_minmodB(node, vlu);
  }

  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    printvar_innode(node, iu);
    printvar_indc(node, iu);
  }

  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    limiter_MRS(node, vlu);
  }

  free_all_myln_indc_for_vl(mesh, vlu);

  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    printvar_innode(node, iu);
    printvar_indc(node, iu);
  }
  */
  //exit(9);

  return 0;
} 

/* calculate errors in u */
int advection1_analyze(tMesh *mesh)
{
  int iu  = Ind("advection1_u");
  int iue = Ind("advection1_u_err");

  if(PR) PRFs("\n");

  /* set correct profile in advection1_u_err */
  advection1_set_profile(mesh, iue);

  /*  compute errors: u_err = u - u_correct */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    double *u  = Vard(node, iu);
    double *ue = Vard(node, iue);
    int i;

    forpoints(node, i)
      ue[i] = u[i]- ue[i];
  }
  return 0;
}


/* refine if diff between min and max is above threshold */
int advection1_refine(tMesh *mesh)
{
  int iu  = Ind("advection1_u");
  int ref_method = PARENT_n;
  double rft  = 0.5;
  double urft = 0.1;

  if(PR) PRFs("\n");

  /* set rflag based on refinement criterion */
  formylnodes(mesh)
  {
    tNode *parent;
    tNode *node = MyLnode;
    //double *u  = Vard(node, iu);
    tArray *uA = VarA(node, iu);
    int imin, imax;
    double min = min_array(uA, &imin);
    double max = max_array(uA, &imax);

    /* should we refine, unrefine, or do nothing */
    if(max-min>rft)       node->rflag = ref_method;
    else if(max-min<urft) node->rflag = -ref_method;
    else                  node->rflag = 0;

    /* find parent if unset rflag if it exists */
    parent = node->parent;
    if(parent) parent->rflag = 0;

    /* do not allow unrefine if we are already at the root node level */
    if(node->rflag<0 && !parent) node->rflag = 0;

    if(node->rflag>0 && node->l>=4) node->rflag = 0;
    //node->rflag = 0;

    if(0 && node->rflag<0)
    {
      char s[100];
      printf("%s min=%g max=%g rflag=%d\n",
      nodename(node,s,99), min, max, node->rflag);
    }
  }

  /* unrefine where needed */
  remove_nodes_if_rflag(mesh, ref_method);
  update_mesh_myln_node_nid(mesh);

  /* refine where needed */
  hrefine_nodes_if_rflag(mesh, ref_method);
  update_mesh_myln_node_nid(mesh);

  return 0;
}
