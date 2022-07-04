/* scalarwave1.c */
/* Wolfgang Tichy, 9/2019 */

#include "nmesh.h"
#include "scalarwave1.h"

#define PR 1


/* use DGglobals */
extern tDGglobals DGglobals[1];


/* global vars/pars for scalarwave1.c */
tscalarwave1 scalarwave1[1];

/* func to init global vars/pars */
int scalarwave1_init_global_pars(tMesh *mesh)
{
  int d;

  PRFs(":\n");

  sscanf(Gets(Par("scalarwave1_k")), "%lg %lg %lg",
         &(scalarwave1->k[0]), &(scalarwave1->k[1]), &(scalarwave1->k[2]));

  printf(" scalarwave1->k = {");
  for(d=0; d<3; d++) printf(" %.16g", scalarwave1->k[d]);
  printf(" }\n");

  scalarwave1->sin_profile    = Getv(Par("scalarwave1_profile"),"sin");
  scalarwave1->square_profile = Getv(Par("scalarwave1_profile"),"square");

  return 0;
}


/* fluxes in direction norm, u = (pi, cx,cy,cz, phi) */
void scalarwave1_flux1d(tNode *node, int ncons, double *f, double norm[3],
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
void scalarwave1_eigenval1d(tNode *node, int ncons, double *lam, double norm[3])
{
  /* eigenvalues in some order */
  lam[0] = 1.0;
  lam[1] = -1.0;
  lam[2] = 0.0;
  lam[3] = 0.0;
  lam[4] = 0.0;
}

/* flux and its derivs for adv. eqn: f^i = n^i u */
void scalarwave1_f_divf(tNode *node, tVarList *vlu)
{
  int iphi = Vind(vlu, 4);
  int ipi  = Vind(vlu, 0);
  int icx  = Vind(vlu, 1);
  int if_pix = scalarwave1->if_pix;
  int if_cxx = scalarwave1->if_cxx;
  int if_cyx = if_cxx + 3;
  int if_czx = if_cxx + 6;
  int idivf_pi = scalarwave1->idivf_pi;
  int idivf_cx = scalarwave1->idivf_cx;

  /* compute derivs */
  {
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

      scalarwave1_flux1d(node,5, f,no, u);
      f_pi_x[i] = f[0];
      f_cx_x[i] = f[1];
      f_cy_x[i] = f[2];
      f_cz_x[i] = f[3];

      no[0] = 0;
      no[1] = 1.;
      scalarwave1_flux1d(node,5, f,no, u);
      f_pi_y[i] = f[0];
      f_cx_y[i] = f[1];
      f_cy_y[i] = f[2];
      f_cz_y[i] = f[3];

      no[1] = 0.;
      no[2] = 1.;
      scalarwave1_flux1d(node,5, f,no, u);
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
  int set_aj = !( (d->info) & DGINFO_INNONLY );
  int nvars = vlu->n;
  double norm[3];
  int l;
  double *u;
  double *uaj = NULL;

  /* set normal based on d->info */
  node_normal_from_DGinfo(d, norm);

  /* eigenval in dir norm */
  scalarwave1_eigenval1d(node,nvars, d->lami,norm);
  if(set_aj) scalarwave1_eigenval1d(node,nvars, d->lama,norm);

  /* reset d->Ffac to default */
  d->Ffac = 1.;

  /* loop over evo vars in vlu */
  forvl(vlu, l)
  {
    int vi = Vind(vlu,l);

    u = Vard_(node, vi);
    if(set_aj) uaj = Varaj(node, vi, f);

    /* cons var inside node */
    d->ui[l] = u[ijk];

    /* cons var on adjacent side */
    if(uaj) /* if there is an adjacent domain */
      d->ua[l] = uaj[JK];
  }

  /* no adjacent u, i.e. we are on outer boundary */
  if(!uaj && set_aj) /* we check the last uaj here. */
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

      /* set flux factor so that we can zero fluxes at outer boundary */
      //d->Ffac = DGglobals->outerBC_flux_fac[dir];
    }
  }

  /* get inner and adjacent fluxes fi, fa */
  scalarwave1_flux1d(node,nvars, d->fi,norm, d->ui);
  if(set_aj) scalarwave1_flux1d(node,nvars, d->fa,norm, d->ua);
}


/* compute div of flux */
void scalarwave1_set_divf(tNode *node, tVarList *vlu)
{
  int use_fv = node->dat->info->use_fv;
  tMesh *mesh = vlu->mesh;
  tVarList *vldivf = vlalloc(mesh);
  vlpush(vldivf, scalarwave1->idivf_pi);
  vlpush(vldivf, scalarwave1->idivf_cx);

  if(use_fv)
  {
    /* linear interpolation to moved in point */
    if(DGglobals->fv_rec_mode >= FV_REC_WENO3if2away_1)
      rec1d_uface_to_uin_1(node, vlu, 1);
    /* WARNING: The interploation does not help for 2 FV neighbors!!!
                BUT maybe it would help if going from FV to DG???  */

    /* compute d_i f^i with finite vol. methods on one node */
    fv_divf(node, vldivf, vlu, vlu,NULL, scalarwave1_rec_u_f_lam,
            scalarwave1_fluxes_pt, scalarwave1->numflux);
  }
  else
  {
    scalarwave1_f_divf(node, vlu);
  }

  vlfree(vldivf);
}

/* RHS of: d_t u = - d_i f^i */
int scalarwave1_vol_rhs_u(tNode *node, tVarList *vlr, tVarList *vlu)
{
  int ipi = Vind(vlu, 0);
  int irpi = Vind(vlr, 0);
  int ircx = Vind(vlr, 1);
  int irphi = Vind(vlr, 4);
  int idivf_pi = scalarwave1->idivf_pi;
  int idivf_cx = scalarwave1->idivf_cx;

  /* set div of flux */
  scalarwave1_set_divf(node, vlu);

  /* RHS */
  {
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
    int fc[6];
    int i;

    /* RHS at each point */
    forpoints(node, i)
    {
      double interior = !(ind_on_nodeface(node, i, fc));
      rphi[i] = pi[i];
      rpi[i] = -divf_pi[i] * interior; //only add divf if we are in interior
      rcx[i] = -divf_cx[i] * interior;
      rcy[i] = -divf_cy[i] * interior;
      rcz[i] = -divf_cz[i] * interior;
    }
  }

  return 0;
}

/* surface terms in RHS of: d_t u */
int scalarwave1_surf_rhs_u(tNode *node, tVarList *vlr, tVarList *vlu)
{
  tMesh *mesh = vlu->mesh;
  int idivf_pi = Ind("scalarwave1_divf_pi");
  int idivf_cx = Ind("scalarwave1_divf_cx");
  tVarList *vldivf;

  /* make var list with div of fluxes */
  vldivf = vlalloc(mesh);
  vlpush(vldivf, idivf_pi);
  vlpush(vldivf, idivf_cx);

  /* extrapolate u back to face on fv nodes */
  if(DGglobals->fv_rec_mode >= FV_REC_WENO3if2away_1)
    rec1d_uface_to_uin_1(node, vlu, 0);

  /* add flux terms on surfaces to divf */
  dg_add_surface_fluxes_sign(node, 1., vldivf, vlu, NULL,
                             scalarwave1_fluxes_pt, scalarwave1->numflux);

  /* extrapolate divf to face on fv nodes */
  if(DGglobals->fv_divf_extrap_mode == FV_DIVF_EXTRAP1)
    rec1d_uface_to_uin_1(node, vldivf, 0);

  /* add divf on node surfaces */
  vladdto_onfaces_node(node, vlr, -1., vldivf);

  vlfree(vldivf);
  return 0;
}

/* Upwind flux, by A. Adhikari */
void scalarwave1_numflux1d_upwind(tDGinfo *d)
{
  int nf = 4; /* only 4 modes besides phi */
  //tNode *node = d->node;
  //int f = d->face;
  //int dir = f/2;
  //int ijk;
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

  /* set normal based on d->info */
  node_normal_from_DGinfo(d, norm);

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


/* set profile in array u */
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

/* set d/dt(profile) in array dudt */
void scalarwave1_set_dt_profile_pt(double xyz[3], double t, int nv,
                                   double *dudt)
{
  double *kv = scalarwave1->k;
  double kx = kv[0];
  double ky = kv[1];
  double kz = kv[2];
  double om = sqrt(kx*kx + ky*ky + kz*kz);
  double u[5];

  /* set u */
  scalarwave1_set_profile_pt(xyz, t, nv, u);

  /* profile */
  if(scalarwave1->sin_profile)
  {
    /* set du/dt = d/dt (pi, cx,cy,cz, phi) */
    dudt[0] = -om*om*u[4];
    dudt[1] = om*kx*u[4];
    dudt[2] = om*ky*u[4];
    dudt[3] = om*kz*u[4];
    dudt[4] = u[0];
  }
  if(scalarwave1->square_profile)
  {
    errorexit("this is wrong!!!");
  }
}

/* set profile in vlu and its RHS in vlr */
int scalarwave1_set_profile(tVarList *vlu, tVarList *vlr)
{
  tMesh *mesh = vlu->mesh;
  int ix =  Ind("x");
  int iphi = Vind(vlu, 4); /* get indices out of vlu */
  int ipi  = Vind(vlu, 0);
  int icx  = Vind(vlu, 1);
  double t = mesh->time;
  int iphir, ipir, icxr;

  /* get indices out of vlr if it is set */
  if(vlr && vlr->n == 5)
  {
    iphir = Vind(vlr, 4);
    ipir  = Vind(vlr, 0);
    icxr  = Vind(vlr, 1);
  }
  else
  {
    iphir = ipir = icxr = ipi; /* ipir = ipi means vlr is NULL */
  }

  /* profile */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    double *x = Vard(node, ix);
    double *y = Vard(node, ix+1);
    double *z = Vard(node, ix+2);
    double *phi = Vard(node, iphi); /* u = (pi, cx,cy,cz, phi) */
    double *pi = Vard(node, ipi);
    double *cx = Vard(node, icx);
    double *cy = Vard(node, icx+1);
    double *cz = Vard(node, icx+2);
    double *phir = Vard(node, iphir); /* RHS of u */
    double *pir = Vard(node, ipir);
    double *cxr = Vard(node, icxr);
    double *cyr = Vard(node, icxr+1);
    double *czr = Vard(node, icxr+2);
    double u5[5] = {0.};
    double du5dt[5] = {0.};
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

      /* get time deriv or RHS of u5 if ipir is set */
      if(ipir != ipi)
      {
        scalarwave1_set_dt_profile_pt(xyz,t, 5, du5dt);
        pir[i]  = du5dt[0];
        cxr[i]  = du5dt[1];
        cyr[i]  = du5dt[2];
        czr[i]  = du5dt[3];
        phir[i] = du5dt[4];
      }
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
  tVarList *vlu = vlalloc(mesh);
  int numflux = Par("scalarwave1_numflux");
  int limiter = Par("scalarwave1_limiter");

  PRF;printf(": dt = %g\n", mesh->dt);

  /* save some global var indices */
  scalarwave1->if_pix = if_pix;
  scalarwave1->if_cxx = if_cxx;
  scalarwave1->idivf_pi = idivf_pi;
  scalarwave1->idivf_cx = idivf_cx;

  /* varlist for evolution: vlu contains indices for (pi, cx,cy,cz, phi) */
  vlpush(vlu, ipi);
  vlpush(vlu, icx);
  vlpush(vlu, iphi);

  /* enable all needed vars */
  enablevarlist(vlu);
  enablevar(mesh, if_pix);
  enablevar(mesh, if_cxx);
  enablevar(mesh, idivf_pi);
  enablevar(mesh, idivf_cx);

  /* set initial profile, e.g. at t=0: set u=sin(x) */
  scalarwave1_set_profile(vlu, NULL);

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
  double time = mesh->time;
  int ipi  = Ind("scalarwave1_pi");
  int icx  = Ind("scalarwave1_cx");
  int iphi = Ind("scalarwave1_phi");
  tVarList *vlu = vlalloc(mesh);
  tVarList *vlr, *vle, *vler; /* RHS of u and errors in u and RHS of u */

  if(PR) PRFs("\n");

  /* u */
  vlpush(vlu, ipi);
  vlpush(vlu, icx);
  vlpush(vlu, iphi);

  /* RHS of u */
  vlr = AddDuplicateEnable(vlu, "_r", AUXVAR, 0);

  /* add vars for errors in vlu and vlr */
  vle  = AddDuplicateEnable(vlu, "_err", AUXVAR, 0);
  vler = AddDuplicateEnable(vlr, "_err", AUXVAR, 0);

  /* set correct profile in vle and correct RHS in vler */
  if(time > 0.)
  {
    mesh->time -= mesh->dt; // decr. mesh->time to t where RHS was evaluated
    scalarwave1_set_profile(vle, vler);
    mesh->time = time;      // restore mesh->time
  }
  scalarwave1_set_profile(vle, NULL);

  /*  compute errors: u_err = u - u_correct */
  vladd(vle,  1.,vlu, -1.,vle);  // vle  = vlu - vle
  vladd(vler, 1.,vlr, -1.,vler); // vler = vlr - vler

  vlfree(vler);
  vlfree(vle);
  vlfree(vlr);
  vlfree(vlu);
  return 0;
}


/***********************************************************************/
/* funcs needed for finite volume method in nmesh */
/***********************************************************************/

/* function that sets cons vars u, fluxes and eigenvals.
   In:  fv->nq,qc,npts,im,q_scale,rec1d_p,rec1d_m, d->node,info,i,j,k,face
   Out: fv->qm_p,qm_m, d->ui,ua,fi,fa,lami,lama */
void scalarwave1_rec_u_f_lam(tFVinfo *fv, tDGinfo *d)
{
  tNode *node = d->node;
  //int *n = node->n;
  int f = d->face;
  int right_face = f%2;
  //int ijk = Ind_n(d->i,d->j,d->k, n);
  int nvars = 5;
  double norm[3];
  int l;

  /* get normal at midpoint left or right of ijk */
  node_normal_from_DGinfo(d, norm);

  /* reconstruct at mid point */
  fv_rec1d_q_midpt(fv);

  /* in scalarwave1 q is u, so transfer qm_p,qm_m into ui,ua now */

  /* right face means, use right midpoint */
  if(right_face)
    for(l=0; l<fv->nq; l++)
    {
      d->ui[l] = fv->qm_p[l];
      d->ua[l] = fv->qm_m[l];
    }
  else
    for(l=0; l<fv->nq; l++)
    {
      d->ui[l] = fv->qm_m[l];
      d->ua[l] = fv->qm_p[l];
    }

  /* eigenval in dir norm */
  scalarwave1_eigenval1d(node,nvars, d->lami,norm);
  scalarwave1_eigenval1d(node,nvars, d->lama,norm);

  /* get inner and adjacent fluxes fi, fa */
  scalarwave1_flux1d(node,nvars, d->fi,norm, d->ui);
  scalarwave1_flux1d(node,nvars, d->fa,norm, d->ua);
}
