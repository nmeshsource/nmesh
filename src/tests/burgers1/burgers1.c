/* burgers1.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "burgers1.h"

#define PR 1


/* func pointer for numerical flux */
void (*burgers1_numflux)(tMesh *mesh, int nf, double *fnum,
                         double *uL, double *uR, double *fL, double *fR,
                         double *lamL, double *lamR);


/* flux in direction norm */
void burgers1_flux1d(tMesh *mesh, int ncons, double *f, double norm[3],
                     double *u)
{
  static int firstcall = 1;
  static double nd[3];

  if(firstcall)
  {
    char *advdir = Gets(Par("burgers1_direction"));
    /* prop. dir.*/
    sscanf(advdir, "%lg %lg %lg", &(nd[0]), &(nd[1]), &(nd[2]));
    firstcall = 0;
  }

  /* flux times norm */
  f[0] = (norm[0]*nd[0] + norm[1]*nd[1] + norm[2]*nd[2]) * 0.5*u[0]*u[0];
}

/* eigenvalue in direction norm */
void burgers1_eigenval1d(tMesh *mesh, int ncons, double *lam, double norm[3],
                         double *u)
{
  static int firstcall = 1;
  static double nd[3];

  if(firstcall)
  {
    char *advdir = Gets(Par("burgers1_direction"));
    /* prop. dir.*/
    sscanf(advdir, "%lg %lg %lg", &(nd[0]), &(nd[1]), &(nd[2]));
    firstcall = 0;
  }

  /* eigenvalue */
  lam[0] = (norm[0]*nd[0] + norm[1]*nd[1] + norm[2]*nd[2]) * u[0];
}

/* flux and its derivs for adv. eqn: f^i = n^i u */
void burgers1_f_df(tMesh *mesh, tVarList *vlu)
{
  int iu = vlu->index[0];
  int ifx = Ind("burgers1_fx");
  int ify = ifx+1;
  int ifz = ifx+2;
  int ifxx = Ind("burgers1_fxx");
  int ifyx = ifxx+3;
  int ifzx = ifxx+6;
  char *advdir = Gets(Par("burgers1_direction"));
  double nx,ny,nz;
  int myid;

  /* prop. dir.*/
  sscanf(advdir, "%lg %lg %lg", &nx, &ny, &nz);

  /* compute derivs */
  formylnodes(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
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
      burgers1_flux1d(mesh,1, &(fx[i]),no, &u_i);
      no[0] = 0;
      no[1] = 1.;
      burgers1_flux1d(mesh,1, &(fy[i]),no, &u_i);
      no[1] = 0.;
      no[2] = 1.;
      burgers1_flux1d(mesh,1, &(fz[i]),no, &u_i);
    }

   /* flux derivs */
   cart_partials(node, ifx, ifxx);
   cart_partials(node, ify, ifyx);
   cart_partials(node, ifz, ifzx);
  }
}


/* function that sets fluxes and eigenvals on both sides of a node surface.
   In: vlu. Out: ui,ua, fi,fa, lami,lama */
void burgers1_fluxes_pt(tNode *node, int face, int i, int j, int k,
                        tVarList *vlu, double *ui, double *ua,
                        double *fi,  double *fa, double *lami, double *lama)
{
  tMesh *mesh = node->pat->mesh;
  int *n = node->n;
  int dir = face/2;
  int ijk = Ind_n(i,j,k, n);
  int JK = Ind_n_norm(i,j,k, n, dir);
  int nvars = vlu->n;
  double norm[3];
  int l;

  forvl(vlu, l)
  {
    int vi = Vind(vlu,l);
    double *u = Vard_(node, vi);
    double *uaj = Varaj(node, vi, face);

    /* cons var inside node, and cons var on adjacent side */
    ui[l] = u[ijk];
    if(uaj)
      ua[l] = uaj[JK];
    else /* do something special on outer boundary */
      ua[l] = 0.;
  }

  /* get face normal at point ijk */
  node_normal_at_ijk(node, face, ijk, norm);

  /* eigenval in dir norm */
  burgers1_eigenval1d(mesh,nvars, lami,norm, ui);
  burgers1_eigenval1d(mesh,nvars, lama,norm, ua);

  /* get inner and adjacent fluxes fi, fa */
  burgers1_flux1d(mesh,nvars, fi,norm, ui);
  burgers1_flux1d(mesh,nvars, fa,norm, ua);
}


/* flux function f in Burgers eqn */
double burgers1_f(double u)
{
  return 0.5*u*u;
}

/* flux at interface for 1d Godunov method */
double F_interface(double ul, double ur)
{
  if(ul >= 0.0 && ur >= 0.0) return burgers1_f(ul);
  if(ul <= 0.0 && ur <= 0.0) return burgers1_f(ur);
  if(ul  < 0.0 && ur >= 0.0) return burgers1_f(0.0);
  else
  {
     double s2 = ul + ur;
     if(s2 > 0) return burgers1_f(ul);
     else       return burgers1_f(ur);
  }
}

/* use numerical flux FN^i to set F */
void burgers1_F(tMesh *mesh, tVarList *vlu)
{
  int iu = vlu->index[0];
  int ifx = Ind("burgers1_fx");
  int iF  = Ind("burgers1_F0");
  char *advdir = Gets(Par("burgers1_direction"));
  double nx,ny,nz;
  int myid;

  /* prop. dir.*/
  sscanf(advdir, "%lg %lg %lg", &nx, &ny, &nz);

  /* compute boundary flux terms */
  formylnodes(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
    int *n = node->n;
    double *u = Vard(node, iu);
    double *fx = Vard(node, ifx);
    double *fy = Vard(node, ifx+1);
    double *fz = Vard(node, ifx+2);
    double norm[3];
    double F_int, FNx,FNy,FNz;
    int face, dir, isP, p, i,j,k, ijk, JK;

    /* set F on each face */
    for(face=0; face<6; face++)
    {
      double *F = Vard(node, iF+face);
      double *uaj = Varaj(node, iu, face);
      dir = face/2;
      isP = face%2;
      p = isP*(n[dir] - 1);
      forplaneN(dir, i,j,k, n, p)
      {
        double ul, ur; //, F_int;
        ijk = Ind_n(i,j,k, n);
        JK = Ind_n_norm(i,j,k, n, dir);
        node_normal_at_ijk(node, face, ijk, norm);

        /* set ul, ur. Depends on if we are at upper (plus) or lower
           end of domain */
        if(isP)
        {
          ul = u[ijk];
          /* if there is an adjacent surface */
          if(uaj) ur = uaj[JK];
          else    ur = ul;
        }
        else
        {
          ur = u[ijk];
          /* if there is an adjacent surface */
          if(uaj) ul = uaj[JK];
          else    ul = ur;
        }
        /* set numerical flux */
        F_int  = F_interface(ul, ur);
        FNx = F_int * nx;
        FNy = F_int * ny;
        FNz = F_int * nz;

        /* project flux onto boundary normal norm[i] */
        F[JK] = (FNx - fx[ijk])*norm[0] +
                (FNy - fy[ijk])*norm[1] +
                (FNz - fz[ijk])*norm[2];
      }
    }
  }
}

/* set a BC on patch boundary */
void burgers1_u_BC(tMesh *mesh, tVarList *vlr, tVarList *vlu)
{
  int ir = vlr->index[0];
  //int iu = vlu->index[0];
  //int ix = Ind("x");
  //char *advdir = Gets(Par("burgers1_direction"));
  //double nx,ny,nz, nmag2;
  int myid;

  /* prop. dir.*/
  //sscanf(advdir, "%lg %lg %lg", &nx, &ny, &nz);
  //nmag2 = (nx*nx + ny*ny + nz*nz);

  /* compute boundary flux terms */
  formylnodes(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
    tPat *pat = node->pat;
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
      tBface *bfaces = pat->bfaces[face];
      dir = face/2;
      p = (face%2)*(n[dir] - 1);

      if(node->patface[face] && bfaces && bfaces->outerbound)
        forplaneN(dir, i,j,k, n, p)
        {
          ijk = Ind_n(i,j,k, n);
          node_normal_at_ijk(node, face, ijk, norm);

          /* do nothing on boundary */
          r[ijk] += 0.;
        }
    }
  }
}

/* RHS of: d_t u = - d_i f^i */
int burgers1_rhs_u(tMesh *mesh, tVarList *vlr, tVarList *vlu)
{
  int ir = vlr->index[0];
  //int iu = vlu->index[0];
  int ifxx = Ind("burgers1_fxx");
  int iF   = Ind("burgers1_F0");
  int iooJ = Ind("det_dXbdx");
  int isqrtdet2gamma0 = Ind("sqrtdet2gamma0");
  int myid;

  TIMER_START;

  /* compute flux */
  burgers1_f_df(mesh, vlu);

  /* RHS */
  formylnodes(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
    double *r  = Vard(node, ir);
    double *fxx = Vard(node, ifxx);
    double *fyy = Vard(node, ifxx+4);
    double *fzz = Vard(node, ifxx+8);
    int i;

    /* RHS at each point */
    forpoints(node, i) r[i] = -(fxx[i] + fyy[i] + fzz[i]);
  }

if(1)
{
  /* get surfaces so that we can compute fluxes */
  get_all_myln_surfaces(mesh);

  /* get flux terms on boundary */
  burgers1_F(mesh, vlu);

  /* add boundary flux terms */
  formylnodes(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
    int *n = node->n;
    double *r  = Vard(node, ir);
    double *ooJ = Vard(node, iooJ);
    int face;

    for(face=0; face<6; face++)
    {
      int dir = face/2;
      int p = (face%2)*(n[dir] - 1);
      //double sig = 2*(face%2) - 1;
      double *sqrtdet2gam = Vard(node, isqrtdet2gamma0+face);
      double *F = Vard(node, iF+face);
      double *w = Wquad(node, dir);
      int i,j,k, ijk, JK, i0;

      forplaneN(dir, i,j,k, n, p)
      {
        ijk = Ind_n(i,j,k, n);
        JK = Ind_n_norm(i,j,k, n, dir);
        i0 = i0_norm(i,j,k, dir);

        r[ijk] -= F[JK] * sqrtdet2gam[JK] * fabs(ooJ[ijk])/ w[i0];
      }
    }
  }
}
else
{
  dg_add_surface_fluxes(mesh, vlr, vlu,
                        burgers1_fluxes_pt, burgers1_numflux);
}
  /* impose outer BC */
  burgers1_u_BC(mesh, vlr, vlu);

  TIMER_STOP;
  return 0;
}


/* initialize burgers1 */
int burgers1_init(tMesh *mesh)
{
  int iu  = Ind("burgers1_u");
  int ifx = Ind("burgers1_fx");
  int ifxx = Ind("burgers1_fxx");
  int iF  = Ind("burgers1_F0");
  int ix =  Ind("x");
  int iue = Ind("burgers1_u_err");
  tVarList *vlu = vlalloc(mesh);
  char *advdir = Gets(Par("burgers1_direction"));
  double nx,ny,nz;
  int myid;

  /* prop. dir.*/
  sscanf(advdir, "%lg %lg %lg", &nx, &ny, &nz);

  PRF;printf(": dt = %g\n", mesh->dt);

  /* varlist */
  vlpush(vlu, iu);

  /* enable all needed vars */
  enablevar(mesh, iu);
  enablevar(mesh, ifx);
  enablevar(mesh, ifxx);
  enablevar(mesh, iF);
  enablevar(mesh, iue);

  /* at t=0: set u=sin(x) */
  formylnodes(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
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
  evolve_register_subsys_u_rhs_src_lim(mesh, vlu, burgers1_rhs_u, 0,
                                       limdata_MRS, limiter_MRS);
  evolve_print_evosys(mesh);

  /* choose numerical flux */
  burgers1_numflux = numflux1d_scalarGodunov;

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
  int myid;

  /* prop. dir.*/
  sscanf(advdir, "%lg %lg %lg", &nx, &ny, &nz);
  //nmag2 = (nx*nx + ny*ny + nz*nz);

  if(PR) PRFs("\n");

  /*  compute errors */
  formylnodes(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
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
