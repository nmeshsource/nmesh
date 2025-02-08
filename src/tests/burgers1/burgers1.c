/* burgers1.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "burgers1.h"

#define PR 1



/* global vars/pars for burgers1.c */
tburgers1 burgers1[1];


/* func to init global vars/pars */
int burgers1_init_global_pars(tMesh *mesh)
{
  int d;
  char *advdir;
  double nd[3];

  PRFs(":\n");

  advdir = Gets(Par("burgers1_direction"));
  /* prop. dir.*/
  sscanf(advdir, "%lg %lg %lg", &(nd[0]), &(nd[1]), &(nd[2]));
  for(d=0; d<3; d++) burgers1->direction[d] = nd[d];

  printf(" burgers1->direction = {");
  for(d=0; d<3; d++) printf(" %.16g", burgers1->direction[d]);
  printf(" }\n");

  /* choose numerical flux */
  burgers1->numflux = numflux1d_scalarGodunov;

  return 0;
}


/* flux in direction norm */
void burgers1_flux1d(tNode *node, int ncons, double *f, double norm[3],
                     double *u)
{
  double *nd = burgers1->direction;

  /* flux times norm */
  f[0] = (norm[0]*nd[0] + norm[1]*nd[1] + norm[2]*nd[2]) * 0.5*u[0]*u[0];
}

/* eigenvalue in direction norm */
void burgers1_eigenval1d(tNode *node, int ncons, double *lam, double norm[3],
                         double *u)
{
  double *nd = burgers1->direction;

  /* eigenvalue */
  lam[0] = (norm[0]*nd[0] + norm[1]*nd[1] + norm[2]*nd[2]) * u[0];
}

/* flux and its derivs for adv. eqn: f^i = n^i u */
void burgers1_f_df(tNode *node, tVarList *vlu)
{
  tMesh *mesh = vlu->mesh;
  int ifx = Ind("burgers1_fx");
  int idivf = Ind("burgers1_divf");
  //char *advdir = Gets(Par("burgers1_direction"));
  int iu = vlu->index[0];
  //double nx,ny,nz;
  double *u = Vard(node, iu);
  double *fx = Vard(node, ifx);
  double *fy = Vard(node, ifx+1);
  double *fz = Vard(node, ifx+2);
  int i;

  /* prop. dir.*/
  //sscanf(advdir, "%lg %lg %lg", &nx, &ny, &nz);

  /* flux at each point */
  forpoints(node, i)
  {
    double u_i = u[i];
    double no[3] = { 1., 0., 0. };
    burgers1_flux1d(node,1, &(fx[i]),no, &u_i);
    no[0] = 0;
    no[1] = 1.;
    burgers1_flux1d(node,1, &(fy[i]),no, &u_i);
    no[1] = 0.;
    no[2] = 1.;
    burgers1_flux1d(node,1, &(fz[i]),no, &u_i);
  }

   /* flux divergence */
   cart_di_Ui(node, ifx, idivf);
}


/* function that sets fluxes and eigenvals on both sides of a node surface.
   In: vlu. Out: ui,ua, fi,fa, lami,lama */
void burgers1_fluxes_pt(tDGinfo *d)
{
  tVarList *vlu = d->vlu;
  tNode *node = d->node;
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
      d->ua[l] = d->ui[l];
  }

  /* get face normal at point ijk */
  node_normal_at_ijk(node, f, ijk, norm);

  /* eigenval in dir norm */
  burgers1_eigenval1d(node,nvars, d->lami,norm, d->ui);
  burgers1_eigenval1d(node,nvars, d->lama,norm, d->ua);

  /* get inner and adjacent fluxes fi, fa */
  burgers1_flux1d(node,nvars, d->fi,norm, d->ui);
  burgers1_flux1d(node,nvars, d->fa,norm, d->ua);
}


/* set a BC on patch boundary */
void burgers1_u_BC(tNode *node, tVarList *vlr, tVarList *vlu)
{
  int ir = vlr->index[0];
  //int iu = vlu->index[0];
  //int ix = Ind("x");

  /* compute boundary flux terms */
  //tPat *pat = node->pat;
  int *n = node->n;
  double *r = Vard(node, ir);
  //double *x = Vard(node, ix);
  //double *y = Vard(node, ix+1);
  //double *z = Vard(node, ix+2);
  //double t = node->time;
  double norm[3];
  int face, dir, p, i,j,k, ijk;

  /* go over each face */
  for(face=0; face<6; face++)
  {
    //tBface *bfaces = pat->bfaces[face];
    dir = face/2;
    p = (face%2)*(n[dir] - 1);

    //if(node->patface[face] && bfaces && bfaces->boundary)
    if(Elm_on_BOUND(node,face))


      forplaneN(dir, i,j,k, n, p)
      {
        ijk = Ind_n(i,j,k, n);
        node_normal_at_ijk(node, face, ijk, norm);

        /* do nothing on boundary */
        r[ijk] += 0.;
      }
  }
}

/* RHS of: d_t u = - d_i f^i */
int burgers1_vol_rhs_u(tNode *node, tEvoVars *evv)
{
  tVarList *vlr = EvoVars_vlr(evv);
  tVarList *vlu = EvoVars_vlu(evv);
  tMesh *mesh = vlu->mesh;
  int ir = vlr->index[0];
  //int iu = vlu->index[0];
  int idivf = Ind("burgers1_divf");
  double *r  = Vard(node, ir);
  double *divf = Vard(node, idivf);
  int i;

  /* compute flux */
  burgers1_f_df(node, vlu);

  /* RHS at each point */
  forpoints(node, i) r[i] = -divf[i];

  return 0;
}

/* surface terms in RHS of: d_t u */
int burgers1_surf_rhs_u(tNode *node, tEvoVars *evv)
{
  tVarList *vlr = EvoVars_vlr(evv);
  tVarList *vlu = EvoVars_vlu(evv);
  /* add boundary flux terms */
  dg_add_surface_fluxes(node, vlr, vlu, NULL,
                        burgers1_fluxes_pt, burgers1->numflux);

  /* impose outer BC */
  burgers1_u_BC(node, vlr, vlu);

  return 0;
}


/* initialize burgers1 */
int burgers1_init(tMesh *mesh)
{
  int iu  = Ind("burgers1_u");
  int ifx = Ind("burgers1_fx");
  int idivf = Ind("burgers1_divf");
  int ix =  Ind("x");
  int iue = Ind("burgers1_u_err");
  tVarList *vlu = vlalloc(mesh);
  int limiter = Par("burgers1_limiter");
  double nx = burgers1->direction[0]; /* prop. dir.*/
  double ny = burgers1->direction[1];
  double nz = burgers1->direction[2];

  PRF;printf(": dt = %g\n", mesh->dt);

  /* varlist */
  vlpush(vlu, iu);

  /* enable all needed vars */
  enablevar(mesh, iu);
  enablevar(mesh, ifx);
  enablevar(mesh, idivf);
  enablevar(mesh, iue);

  /* at t=0: set u=sin(x) */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    double *u = Vard(node, iu);
    double *x = Vard(node, ix);
    double *y = Vard(node, ix+1);
    double *z = Vard(node, ix+2);
    int i;
    forpoints(node, i)
    {
      double d = (nx*x[i] + ny*y[i] + nz*z[i]);
      if(d<1.0) u[i] = 1.5;
      else      u[i] = 0.5;
    }
  }

  /* register u and its RHS with evolve */
  evolve_register_vl(vlu);
  evolve_SetFun(VOLRHS,  burgers1_vol_rhs_u,  vlu);
  evolve_SetFun(SURFRHS, burgers1_surf_rhs_u, vlu);
  if(Getv(limiter, "MRS"))
  {
    evolve_SetFun(LIMDATA, limdata_MRS_evv, vlu);
    evolve_SetFun(LIMITER, limiter_MRS_evv, vlu);
  }
  else if(Getv(limiter, "minmodB"))
  {
    evolve_SetFun(LIMDATA, limdata_c000_100_010_001_evv, vlu);
    evolve_SetFun(LIMITER, limiter_minmodB_evv, vlu);
  }
  evolve_print_evosys(mesh);

  return 0;
} 

/* calculate errors in u */
int burgers1_analyze(tMesh *mesh)
{
  int iu  = Ind("burgers1_u");
  int iue = Ind("burgers1_u_err");
  int ix =  Ind("x");
  char *advdir = Gets(Par("burgers1_direction"));
  double nx,ny,nz;
  //double nmag2;

  /* prop. dir.*/
  sscanf(advdir, "%lg %lg %lg", &nx, &ny, &nz);
  //nmag2 = (nx*nx + ny*ny + nz*nz);

  if(PR) PRFs("\n");

  /*  compute errors */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    double *u = Vard(node, iu);
    double *ue = Vard(node, iue);
    double *x = Vard(node, ix);
    double *y = Vard(node, ix+1);
    double *z = Vard(node, ix+2);
    double t = mesh->time;
    int i;

    forpoints(node, i)
    {
      double ua;
      double d = (nx*x[i] + ny*y[i] + nz*z[i]);

      if(d < 1.0 + t) ua = 1.5;
      else            ua = 0.5;

      ue[i] = fabs(u[i]- ua);
    }
  }
  return 0;
}
