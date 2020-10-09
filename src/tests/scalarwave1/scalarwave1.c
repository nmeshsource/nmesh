/* scalarwave1.c */
/* Wolfgang Tichy, 9/2019 */

#include "nmesh.h"
#include "scalarwave1.h"

#define PR 1


/* global vars/pars */
tscalarwave1 scalarwave1[1];

/* func to init global vars/pars */
int scalarwave1_init_global_pars(tMesh *mesh)
{
  char *vec = Gets(Par("scalarwave1_k"));

  sscanf(vec, "%lg %lg %lg",
         &(scalarwave1->k[0]), &(scalarwave1->k[1]), &(scalarwave1->k[2]));

  printf("setting: scalarwave1->k[0] = %g\n", scalarwave1->k[0]);
  printf("setting: scalarwave1->k[1] = %g\n", scalarwave1->k[1]);
  printf("setting: scalarwave1->k[2] = %g\n", scalarwave1->k[2]);

  scalarwave1->sin_profile    = Getv(Par("scalarwave1_profile"),"sin");
  scalarwave1->square_profile = Getv(Par("scalarwave1_profile"),"square");

  return 0;
}


/* fluxes in direction norm, u = (pi, cx,cy,cz, phi) */
void scalarwave1_flux1d(tMesh *mesh, int ncons, double *f, double norm[3],
                        double *u)
{
  /* fluxes times norm */
  f[0] = -(u[1]*norm[0] + u[2]*norm[1] + u[3]*norm[2]);
  f[1] = -u[0]*norm[0];
  f[2] = -u[0]*norm[1];
  f[3] = -u[0]*norm[2];
  f[4] = 0.;
}

/* eigenvalues in direction norm */
void scalarwave1_eigenval1d(tMesh *mesh, int ncons, double *lam, double norm[3])
{
  /* eigenvalues in some order */
  lam[0] = 1.0;
  lam[1] = -1.0;
  lam[2] = 0.0;
  lam[3] = 0.0;
  lam[4] = 0.0;
}

/* flux and its derivs for adv. eqn: f^i = n^i u */
void scalarwave1_f_divf(tMesh *mesh, tVarList *vlu)
{
  int iphi = Vind(vlu, 4);
  int ipi  = Vind(vlu, 0);
  int icx  = Vind(vlu, 1);
  int if_pix = Ind("scalarwave1_f_pix");
  int if_cxx = Ind("scalarwave1_f_cxx");
  int if_cyx = if_cxx + 3;
  int if_czx = if_cxx + 6;
  int idivf_pi = Ind("scalarwave1_divf_pi");
  int idivf_cx = Ind("scalarwave1_divf_cx");

  /* compute derivs */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    double *phi = Vard(node, iphi);
    double *pi = Vard(node, ipi);
    double *cx = Vard(node, icx);
    double *cy = Vard(node, icx+1);
    double *cz = Vard(node, icx+2);

    double *f_pi_x = Vard(node, if_pix);
    double *f_cx_x = Vard(node, if_cxx);
    double *f_cy_x = Vard(node, if_cyx);
    double *f_cz_x = Vard(node, if_czx);

    double *f_pi_y = Vard(node, if_pix+1);
    double *f_cx_y = Vard(node, if_cxx+1);
    double *f_cy_y = Vard(node, if_cyx+1);
    double *f_cz_y = Vard(node, if_czx+1);

    double *f_pi_z = Vard(node, if_pix+2);
    double *f_cx_z = Vard(node, if_cxx+2);
    double *f_cy_z = Vard(node, if_cyx+2);
    double *f_cz_z = Vard(node, if_czx+2);
    int i;

    /* flux at each point */
    forpoints(node, i)
    {
      double u[] = { pi[i], cx[i],cy[i],cz[i], phi[i] };
      double f[5];
      double no[3] = { 1., 0., 0. };

      scalarwave1_flux1d(mesh,5, f,no, u);
      f_pi_x[i] = f[0];
      f_cx_x[i] = f[1];
      f_cy_x[i] = f[2];
      f_cz_x[i] = f[3];

      no[0] = 0;
      no[1] = 1.;
      scalarwave1_flux1d(mesh,5, f,no, u);
      f_pi_y[i] = f[0];
      f_cx_y[i] = f[1];
      f_cy_y[i] = f[2];
      f_cz_y[i] = f[3];

      no[1] = 0.;
      no[2] = 1.;
      scalarwave1_flux1d(mesh,5, f,no, u);
      f_pi_z[i] = f[0];
      f_cx_z[i] = f[1];
      f_cy_z[i] = f[2];
      f_cz_z[i] = f[3];
    }

    /* flux derivs */
    cart_di_Ui(node, if_pix, idivf_pi);
    cart_di_Ui(node, if_cxx, idivf_cx);
    cart_di_Ui(node, if_cyx, idivf_cx+1);
    cart_di_Ui(node, if_czx, idivf_cx+2);
  }
}

/* function that sets fluxes and eigenvals on both sides of a node surface.
   In: vlu. Out: ui,ua, fi,fa, lami,lama */
void scalarwave1_fluxes_pt(tDGinfo *d)
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
  double *u;
  double *uaj = NULL;

  /* get face normal norm at point ijk */
  node_normal_at_ijk(node, f, ijk, norm);

  /* eigenval in dir norm */
  scalarwave1_eigenval1d(mesh,nvars, d->lami,norm);
  scalarwave1_eigenval1d(mesh,nvars, d->lama,norm);

  /* loop over evo vars in vlu */
  forvl(vlu, l)
  {
    int vi = Vind(vlu,l);

    u = Vard_(node, vi);
    uaj = Varaj(node, vi, f);

    /* cons var inside node, and cons var on adjacent side */
    d->ui[l] = u[ijk];
    if(uaj) /* if there is an adjacent domain */
      d->ua[l] = uaj[JK];
  }

  /* no adjacent u, i.e. we are on outer boundary */
  if(!uaj) /* we check the last uaj here. */
  {
    int ix = Ind("x");
    double *x = Vard_(node, ix);
    double *y = Vard_(node, ix+1);
    double *z = Vard_(node, ix+2);
    double t = node->time;
    tPat *pat = node->pat;
    tBface *bfaces = pat->bfaces[f];

    /* compute boundary flux terms, if on outer bound */
    if(node->patface[f] && bfaces && bfaces->boundary)
    {
      double xyz[] = { x[ijk],y[ijk],z[ijk] };
      double u5[5];

      /* get u5 = (pi, cx,cy,cz, phi) */
      scalarwave1_set_profile_pt(xyz,t, 5, u5);
      for(l=0; l<nvars; l++) d->ua[l] = u5[l];
    }
  }

  /* get inner and adjacent fluxes fi, fa */
  scalarwave1_flux1d(mesh,nvars, d->fi,norm, d->ui);
  scalarwave1_flux1d(mesh,nvars, d->fa,norm, d->ua);
}


/* RHS of: d_t u = - d_i f^i */
int scalarwave1_vol_rhs_u(tMesh *mesh, tVarList *vlr, tVarList *vlu)
{
  int ipi = Vind(vlu, 0);
  int irpi = Vind(vlr, 0);
  int ircx = Vind(vlr, 1);
  int irphi = Vind(vlr, 4);
  int idivf_pi = Ind("scalarwave1_divf_pi");
  int idivf_cx = Ind("scalarwave1_divf_cx");

  TIMER_START;

  /* compute flux */
  if(Getv(Par("scalarwave1_nummethod"), "dg"))
    scalarwave1_f_divf(mesh, vlu);
  else
    scalarwave1_divf_FV(mesh, vlu);

  /* RHS */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    double *pi = Vard(node, ipi);
    double *rpi = Vard(node, irpi);
    double *rcx = Vard(node, ircx);
    double *rcy = Vard(node, ircx+1);
    double *rcz = Vard(node, ircx+2);
    double *rphi = Vard(node, irphi);
    double *divf_pi = Vard(node, idivf_pi);
    double *divf_cx = Vard(node, idivf_cx);
    double *divf_cy = Vard(node, idivf_cx+1);
    double *divf_cz = Vard(node, idivf_cx+2);
    int i;

    /* RHS at each point */
    forpoints(node, i)
    {
      rphi[i] = pi[i];
      rpi[i] = -divf_pi[i];
      rcx[i] = -divf_cx[i];
      rcy[i] = -divf_cy[i];
      rcz[i] = -divf_cz[i];
    }
  }

  TIMER_STOP;
  return 0;
}

/* surface terms in RHS of: d_t u */
int scalarwave1_surf_rhs_u(tMesh *mesh, tVarList *vlr, tVarList *vlu)
{
  TIMER_START;

  /* get flux terms on surfaces */
  dg_add_surface_fluxes(mesh, vlr, vlu, NULL,
                        scalarwave1_fluxes_pt, scalarwave1->numflux);

  TIMER_STOP;
  return 0;
}

/* Upwind flux, by A. Adhikari */
void scalarwave1_numflux1d_upwind(tDGinfo *d)
{
  int nf = 4; /* only 4 modes besides phi */
  tNode *node = d->node;
  int f = d->face;
  //int dir = f/2;
  int ijk;
  int l,m,s;
  double lam_p[4][4];// = {{0.}};  //positive eigenvalue diagonal matrix
  double lam_n[4][4];// = {{0.}};  //negative eigenvalue diagonal matrix
  double int_p[4][4];// = {{0.}};  //intermidiate calculation matrix
  double int_n[4][4];// = {{0.}};  //intermidiate calculation matrix
  double sim_p[4][4];// = {{0.}};  //similarity trans. of lam_p
  double sim_n[4][4];// = {{0.}};  //similarity trans. of lam_n
  double fnumi[4];// = {0.}; //numflux contribution from element u
  double fnuma[4];// = {0.}; //numflux contribution from adjacent u
  double norm[3];
  double nx, ny, nz;

  /* standard DG point */
  if(d->info == 0)
  {
    int *n = node->n;
    ijk = Ind_n(d->i,d->j,d->k, n);
    //JK = Ind_n_norm(d->i,d->j,d->k, n, dir);
    /* get face normal norm at point ijk */
    node_normal_at_ijk(node, f, ijk, norm);
  }
  else /* interpret i,j,k as midpoint indices */
  {
    int *n = node->n;
    ijk = Ind_n(d->i,d->j,d->k, n);
    //JK = Ind_n_norm(d->i,d->j,d->k, n, dir);
    /* get normal norm at midpoint ijk */
    node_normal_at_midpt_ijk(node, f, ijk, norm);
  }

  /* set face normal */
  nx = norm[0];
  ny = norm[1];
  nz = norm[2];

  for(l=0; l<nf; l++)
  {
    for(m=0; m<nf; m++)
    {
      lam_p[l][m] = 0.;
      lam_n[l][m] = 0.;
      int_p[l][m] = 0.;
      int_n[l][m] = 0.;
      sim_p[l][m] = 0.;
      sim_n[l][m] = 0.;
    }
    fnumi[l] = 0.;
    fnuma[l] = 0.;
  }

  for(l=0; l<nf; l++)
  {
    if(d->lami[l] > 0.)  lam_p[l][l] = d->lami[l];
    else                 lam_n[l][l] = d->lami[l];
  }


  /* For plane parallel to XY plane */
  if(nz!=0)
  {
    /* Transformation matrix */
    double R_XY[4][4] = {{1.0, 1.0, 0.0, 0.0},
                         {-nx, nx, 0.0, nz},
                         {-ny, ny, nz, 0.0},
                         {-nz, nz, -ny, -nx}};

    /* Inverse Transformation matrix */
    double inR_XY[4][4] = {{1.0/2.0, -nx/2.0, -ny/2.0, -nz/2.0},
                           {1.0/2.0, nx/2.0, ny/2.0, nz/2.0},
                           {0.0, -(nx*ny)/nz, (nx*nx + nz*nz)/nz, -ny},
                           {0.0, (ny*ny + nz*nz)/nz, -(nx*ny)/nz, -nx}};

    for(l=0; l<nf; l++)
      for(m=0; m<nf; m++)
        for(s=0; s<nf; s++)
        {
          int_p[l][m] += lam_p[l][s]*inR_XY[s][m];
          int_n[l][m] += lam_n[l][s]*inR_XY[s][m];
        }

    for(l=0; l<nf; l++)
      for(m=0; m<nf; m++)
        for(s=0; s<nf; s++)
        {
          sim_p[l][m] += R_XY[l][s]*int_p[s][m];
          sim_n[l][m] += R_XY[l][s]*int_n[s][m];
        }

    for(l=0; l<nf; l++)
      for(m=0; m<nf; m++)
      {
        fnumi[l] += sim_p[l][m]*d->ui[m];
        fnuma[l] += sim_n[l][m]*d->ua[m];
      }

    for(l=0; l<nf; l++)
      d->fnum[l] = fnumi[l] + fnuma[l];
  }

  /* For plane parallel to YZ plane */
  else if(nx!=0)
  {
    /* Transformation matrix */
    double R_YZ[4][4] = {{1.0, 1.0, 0.0, 0.0},
                         {-nx, nx, -nz, -ny},
                         {-ny, ny, 0.0, nx},
                         {-nz, nz, nx, 0.0}};

    /* Inverse Transformation matrix */
    double inR_YZ[4][4] = {{1.0/2.0, -nx/2.0, -ny/2.0, -nz/2.0},
                           {1.0/2.0, nx/2.0, ny/2.0, nz/2.0},
                           {0.0, -nz, -(ny*nz)/nx, (nx*nx + ny*ny)/nx},
                           {0.0, -ny, (nx*nx +nz*nz)/nx, -(ny*nz)/nx}};
    for(l=0; l<nf; l++)
      for(m=0; m<nf; m++)
        for(s=0; s<nf; s++)
        {
          int_p[l][m] += lam_p[l][s]*inR_YZ[s][m];
          int_n[l][m] += lam_n[l][s]*inR_YZ[s][m];
        }

    for(l=0; l<nf; l++)
      for(m=0; m<nf; m++)
        for(s=0; s<nf; s++)
        {
          sim_p[l][m] += R_YZ[l][s]*int_p[s][m];
          sim_n[l][m] += R_YZ[l][s]*int_n[s][m];
        }

    for(l=0; l<nf; l++)
      for(m=0; m<nf; m++)
      {
        fnumi[l] += sim_p[l][m]*d->ui[m];
        fnuma[l] += sim_n[l][m]*d->ua[m];
      }

    for(l=0; l<nf; l++)
      d->fnum[l] = fnumi[l] + fnuma[l];
  }

  /* For plane parallel to ZX plane */
  else
  {
    /* Transformation matrix */
    double R_ZX[4][4] = {{1.0, 1.0, 0.0, 0.0},
                         {-nx, nx, 0.0, ny},
                         {-ny, ny, -nz, -nx},
                         {-nz, nz, ny, 0.0}};

    /* Inverse Transformation matrix */
    double inR_ZX[4][4] = {{1.0/2.0, -nx/2.0, -ny/2.0, -nz/2.0},
                           {1.0/2.0, nx/2.0, ny/2.0, nz/2.0},
                           {0.0, -(nx*nz)/ny, -nz, (nx*nx + ny*ny)/ny},
                           {0.0, (ny*ny + nz*nz)/ny, -nx, -(nx*nz)/ny}};
    for(l=0; l<nf; l++)
      for(m=0; m<nf; m++)
        for(s=0; s<nf; s++)
        {
          int_p[l][m] += lam_p[l][s]*inR_ZX[s][m];
          int_n[l][m] += lam_n[l][s]*inR_ZX[s][m];
        }

    for(l=0; l<nf; l++)
      for(m=0; m<nf; m++)
        for(s=0; s<nf; s++)
        {
          sim_p[l][m] += R_ZX[l][s]*int_p[s][m];
          sim_n[l][m] += R_ZX[l][s]*int_n[s][m];
        }

    for(l=0; l<nf; l++)
      for(m=0; m<nf; m++)
      {
        fnumi[l] += sim_p[l][m]*d->ui[m];
        fnuma[l] += sim_n[l][m]*d->ua[m];
      }

    for(l=0; l<nf; l++)
      d->fnum[l] = fnumi[l] + fnuma[l];
  }

  /* the flux for scalar field phi itself is zero */
  d->fnum[4] = 0.;
}


/* set profile in var with index iphi */
void scalarwave1_set_profile_pt(double xyz[3], double t, int nv, double *u)
{
  double *kv = scalarwave1->k;
  double kx = kv[0];
  double ky = kv[1];
  double kz = kv[2];
  double om = sqrt(kx*kx + ky*ky + kz*kz);

  /* profile */
  if(scalarwave1->sin_profile)
  {
    double arg = kx*xyz[0] + ky*xyz[1] + kz*xyz[2] - om*t;
    double s = sin(arg);
    double c = cos(arg);

    /* u = (pi, cx,cy,cz, phi) */
    u[4] = s;
    u[0] = -om*c;
    u[1] = kx*c;
    u[2] = ky*c;
    u[3] = kz*c;
  }
  if(scalarwave1->square_profile)
  {
    double inx, iny;
    if(xyz[0]>=(-0.7 + kx*t) && xyz[0]<=(-0.3 + kx*t)) inx = 1.;
    else                                               inx = 0.;
    if(xyz[1]>=(-0.2 + ky*t) && xyz[1]<=(+0.2 + ky*t)) iny = 1.;
    else                                               iny = 0.;

    errorexit("this is wrong!!!");
    /* u = (pi, cx,cy,cz, phi) */
    u[4] = inx*iny;
    u[0] = 0.;
    u[1] = 0.;
    u[2] = 0.;
    u[3] = 0.;
  }
}

/* set profile in var with index iphi */
int scalarwave1_set_profile(tVarList *vlu)
{
  tMesh *mesh = vlu->mesh;
  int ix =  Ind("x");
  int iphi = Vind(vlu, 4);
  int ipi  = Vind(vlu, 0);
  int icx  = Vind(vlu, 1);
  double t = mesh->time;

  /* profile */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    double *x = Vard(node, ix);
    double *y = Vard(node, ix+1);
    double *z = Vard(node, ix+2);
    double *phi = Vard(node, iphi);
    double *pi = Vard(node, ipi);
    double *cx = Vard(node, icx);
    double *cy = Vard(node, icx+1);
    double *cz = Vard(node, icx+2);
    double u5[5] = {0.};
    int i;

    forpoints(node, i)
    {
      double xyz[] = { x[i],y[i],z[i] };

      /* get u5 = (pi, cx,cy,cz, phi) */
      scalarwave1_set_profile_pt(xyz,t, 5, u5);
      pi[i]  = u5[0];
      cx[i]  = u5[1];
      cy[i]  = u5[2];
      cz[i]  = u5[3];
      phi[i] = u5[4];
    }
  }
  return 0;
}


/* initialize scalarwave1 */
int scalarwave1_init(tMesh *mesh)
{
  int iphi = Ind("scalarwave1_phi");
  int ipi  = Ind("scalarwave1_pi");
  int icx  = Ind("scalarwave1_cx");
  int if_pix = Ind("scalarwave1_f_pix");
  int if_cxx = Ind("scalarwave1_f_cxx");
  //int if_cyx = if_cxx + 3;
  //int if_czx = if_cxx + 6;
  int idivf_pi = Ind("scalarwave1_divf_pi");
  int idivf_cx = Ind("scalarwave1_divf_cx");
  int iphie = Ind("scalarwave1_err_phi");
  int ipie  = Ind("scalarwave1_err_pi");
  int icxe  = Ind("scalarwave1_err_cx");
  tVarList *vlu = vlalloc(mesh);
  int numflux = Par("scalarwave1_numflux");
  int limiter = Par("scalarwave1_limiter");

  PRF;printf(": dt = %g\n", mesh->dt);

  /* varlist */
  vlpush(vlu, ipi);
  vlpush(vlu, icx);
  vlpush(vlu, iphi);

  /* enable all needed vars */
  enablevarlist(vlu);
  enablevar(mesh, if_pix);
  enablevar(mesh, if_cxx);
  enablevar(mesh, idivf_pi);
  enablevar(mesh, idivf_cx);
  enablevar(mesh, iphie);
  enablevar(mesh, ipie);
  enablevar(mesh, icxe);

  /* set initial profile, e.g. at t=0: set u=sin(x) */
  scalarwave1_set_profile(vlu);

  /* register u and its RHS with evolve */
  if(Getv(limiter, "MRS"))
    evolve_register_subsys_u_rhs_lim(mesh, vlu, scalarwave1_vol_rhs_u,
                                     scalarwave1_surf_rhs_u,
                                     limdata_MRS, limiter_MRS);
  else if(Getv(limiter, "minmodB"))
    evolve_register_subsys_u_rhs_lim(mesh, vlu, scalarwave1_vol_rhs_u,
                                   scalarwave1_surf_rhs_u,
                                   limdata_c000_100_010_001, limiter_minmodB);
  else
    evolve_register_subsys_u_rhs_lim(mesh, vlu, scalarwave1_vol_rhs_u,
                                     scalarwave1_surf_rhs_u, 0,0);
  evolve_print_evosys(mesh);

  /* choose numerical flux */
  if(Getv(numflux, "LLF"))
    scalarwave1->numflux = numflux1d_LLF;
  else
    scalarwave1->numflux = scalarwave1_numflux1d_upwind;

  return 0;
}

/* calculate errors in u */
int scalarwave1_analyze(tMesh *mesh)
{
  int iphi = Ind("scalarwave1_phi");
  int ipi  = Ind("scalarwave1_pi");
  int icx  = Ind("scalarwave1_cx");
  int iphie = Ind("scalarwave1_err_phi");
  int ipie  = Ind("scalarwave1_err_pi");
  int icxe  = Ind("scalarwave1_err_cx");
  tVarList *vle = vlalloc(mesh);

  if(PR) PRFs("\n");

  vlpush(vle, ipie);
  vlpush(vle, icxe);
  vlpush(vle, iphie);

  /* set correct profile in scalarwave1_err... */
  scalarwave1_set_profile(vle);

  /*  compute errors: u_err = u - u_correct */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    double *phi = Vard(node, iphi);
    double *pi  = Vard(node, ipi);
    double *cx  = Vard(node, icx);
    double *cy  = Vard(node, icx+1);
    double *cz  = Vard(node, icx+2);
    double *phie = Vard(node, iphie);
    double *pie  = Vard(node, ipie);
    double *cxe  = Vard(node, icxe);
    double *cye  = Vard(node, icxe+1);
    double *cze  = Vard(node, icxe+2);
    int i;

    forpoints(node, i)
    {
      pie[i] = pi[i]- pie[i];
      cxe[i] = cx[i]- cxe[i];
      cye[i] = cy[i]- cye[i];
      cze[i] = cz[i]- cze[i];
      phie[i] = phi[i]- phie[i];
    }
  }

  vlfree(vle);
  return 0;
}


/***********************************************************************/
/* funcs needed for finite volume method in nmesh */
/***********************************************************************/

/* compute d_i f^i with finite vol. methods */
void scalarwave1_divf_FV(tMesh *mesh, tVarList *vlu)
{
  int nvars = vlu->n - 1; /* only 4 since phi never contributes to any flux */
  int ipi = Vind(vlu, 0);
  int icx = Vind(vlu, 1);
  int idivf_pi = Ind("scalarwave1_divf_pi");
  int idivf_cx = Ind("scalarwave1_divf_cx");
  int iXm_sqrtgdiagx = Ind("Xm_sqrtgdiagx");
  int iYm_sqrtgdiagx = Ind("Ym_sqrtgdiagx");
  int iZm_sqrtgdiagx = Ind("Zm_sqrtgdiagx");
  tVarList *vldivf = vlalloc(mesh);

  TIMER_START;

  /* make var list for div of fluxes */
  vlpush(vldivf, idivf_pi);
  vlpush(vldivf, idivf_cx);

  /* RHS */
  formylnodes(mesh)
  {
    tDGinfo *d = alloc_DGinfo(vlu, NULL); /* each thread gets its own d */
    tNode *node = MyLnode;
    double *pi = Vard(node, ipi);
    double *cx = Vard(node, icx);
    double *cy = Vard(node, icx+1);
    double *cz = Vard(node, icx+2);
    double *m_sqrtgdiag[3][3] =
      { { Vard(node, iXm_sqrtgdiagx), Vard(node, iXm_sqrtgdiagx+1),
                                           Vard(node, iXm_sqrtgdiagx+2) },
        { Vard(node, iYm_sqrtgdiagx), Vard(node, iYm_sqrtgdiagx+1),
                                           Vard(node, iYm_sqrtgdiagx+2) },
        { Vard(node, iZm_sqrtgdiagx), Vard(node, iZm_sqrtgdiagx+1),
                                           Vard(node, iZm_sqrtgdiagx+2) } };
    int *n = node->n;
    int maxn = max3(n[0],n[1],n[2]);
    double *Xbm = dmalloc(maxn);
    double *dXb = dmalloc(maxn);
    double (*uc)[maxn] = dtensor(nvars*maxn);    //array for the 4 pi,cx,cy,cz
    double (*fnumR)[maxn] = dtensor(nvars*maxn); //array for the 4 fluxes
    double um_p[5]; /* array for the 5 pi,cx,cy,cz, phi */
    double um_m[5]; /* array for the 5 pi,cx,cy,cz, phi */
    int dir;

    /* write node into d because numflux needs this */
    d->node = node;
    d->info = 1; // anything other than 0 triggers normals on midpoints

    /* force node into FV mode */
    node->dat->use_fv = 1;

    /* add fluxes in each direction to RHS */
    for(dir=0; dir<3; dir++)
    {
      double *msqrtgdiag = (double *) (m_sqrtgdiag[dir][0]);
      //double *Xb;
      int i,j,k;

      /* get grid points */
      // Xb = Xbpts(node, dir);
      /* we can shift grid points of the nodefaces with this:
           shift_Xb0_XbN_toward_Xbm0_XbmN(Xbm, n[dir], Xb);
         but then really changes the points in the node!
         BETTER make Xb a copy of Xbpts(node, dir). */
      /* get midpoints */
      set_nm_nodemidpoints_Xb_dir(node, n[dir]-1,0, dir, Xbm);
      set_nm_nodemidpoint_distsXb_dir(node, dir, Xbm, dXb);

      /* loop over plane */
      forplaneN(dir, i,j,k, n, 0)
      {
        int l; /* field index */
        int i1 = i1_norm(i,j,k, dir); /* 1st and 2nd index in plane */
        int i2 = i2_norm(i,j,k, dir);
        int i0;                       /* index orthogonal to plane */
        int ic,jc,kc, ccc;

        /* fill field arrays uc, i0 runs orth. to plane */
        for(i0=0; i0<n[dir]; i0++)
        {
          ijk_inplaneN(dir, ic,jc,kc, i1,i2, i0);
          ccc = Ind_n(ic,jc,kc, n);
          uc[0][i0] = pi[ccc];
          uc[1][i0] = cx[ccc];
          uc[2][i0] = cy[ccc];
          uc[3][i0] = cz[ccc];
        }

        /* loop over points in dir */
        for(i0=0; i0<n[dir]; i0++)
        {
          double normR[3];
          int im0, im0m1, im,jm,km, cccm1;
          double i0g0, i0lN, gd_ow_m, gd_ow_m1, wm, wm1;

          /* set points and their index */
          ijk_inplaneN(dir, ic,jc,kc, i1,i2, i0);
          ccc = Ind_n(ic,jc,kc, n);

          /* set index of right and left midpoint and some flags if we
             are at endpoints */
          if(i0>0) { i0g0=1; im0m1 = i0-1; }
          else     { i0g0=0; im0m1 = i0; /* save value */ }
          if(i0<n[dir]-1) { i0lN=1; im0 = i0; }
          else            { i0lN=0; im0 = i0-1; /* save value */ }

          /* set midpoint to left and its index */
          ijk_inplaneN(dir, im,jm,km, i1,i2,im0m1);
          cccm1 = Ind_n(im,jm,km, n);

          /* if im0 is within range */
          if(im0<n[dir]-1)
          {
            /* interpolate fields 0-3 towards the midpoints */
            for(l=0; l<nvars; l++)
            {
              /* reconstruct from both sides of midpoint at i0m */
              um_p[l] = rec1d_p_0(n[dir], uc[l], im0);
              um_m[l] = rec1d_m_0(n[dir], uc[l], im0);
            }

            /* get right pointing normal at midpoint ccc */
            node_normal_at_midpt_ijk(node, dir*2+1, ccc, normR);

            /* set fields */
            for(l=0; l<nvars; l++)
            {
              d->ui[l] = um_p[l];
              d->ua[l] = um_m[l];
            }

            /* eigenval for normR */
            scalarwave1_eigenval1d(mesh,nvars, d->lami,normR);
            scalarwave1_eigenval1d(mesh,nvars, d->lama,normR);

            /* compute fluxes from fields */
            scalarwave1_flux1d(mesh,nvars, d->fi,normR, d->ui);
            scalarwave1_flux1d(mesh,nvars, d->fa,normR, d->ua);

            /* compute numerical flux */
            scalarwave1->numflux(d);
            for(l=0; l<nvars; l++)
              fnumR[l][im0] = d->fnum[l];
            /* we could now also compute fnumL with normL=-normR, but I think
               this results in fnumL = -fnumR */
          }

          /* factors in flux terms on RHS at right and left midpoint */
          wm  = dXb[im0];
          wm1 = dXb[im0m1];
          gd_ow_m  = i0lN * msqrtgdiag[ccc]/wm;
          gd_ow_m1 = i0g0 * msqrtgdiag[cccm1]/wm1;

          /* get divf with FV method */
          for(l=0; l<nvars; l++)
          {
            int idivf = Vind(vldivf,l);
            double *divf = Vard_(node, idivf);
            double *fnum = fnumR[l];

            divf[ccc] = fnum[im0]*gd_ow_m - fnum[im0m1]*gd_ow_m1;
          }
        }
      } /* end plane loop */
    } /* end dir-loop*/

    /* release mem */
    free(fnumR);
    free(uc);
    free(dXb);
    free(Xbm);
    free_DGinfo(d);
  }

  vlfree(vldivf);

  TIMER_STOP;
}
