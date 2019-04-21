/* advection1.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "advection1.h"

#define PR 1


/* func pointer for numerical flux */
void (*advection1_numflux)(tNode *node, int face, int nf, double *fnum,
                           double *uL, double *uR, double *fL, double *fR,
                           double *lamL, double *lamR);


/* flux in direction norm */
void advection1_flux1d(tMesh *mesh, int ncons, double *f, double norm[3],
                       double *u)
{
  static int firstcall = 1;
  static double nd[3];

  if(firstcall)
  {
    char *advdir = Gets(Par("advection1_direction"));
    /* prop. dir.*/
    sscanf(advdir, "%lg %lg %lg", &(nd[0]), &(nd[1]), &(nd[2]));
    firstcall = 0;
  }

  /* flux times norm */
  f[0] = (norm[0]*nd[0] + norm[1]*nd[1] + norm[2]*nd[2]) * u[0];
}

/* eigenvalue in direction norm */
void advection1_eigenval1d(tMesh *mesh, int ncons, double *lam, double norm[3])
{
  static int firstcall = 1;
  static double nd[3];

  if(firstcall)
  {
    char *advdir = Gets(Par("advection1_direction"));
    /* prop. dir.*/
    sscanf(advdir, "%lg %lg %lg", &(nd[0]), &(nd[1]), &(nd[2]));
    firstcall = 0;
  }

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
  int myid;

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
void advection1_fluxes_pt(tNode *node, int face, int i, int j, int k,
                          tVarList *vlu,
                          double *ui, double *ua,
                          double *fi,  double *fa,
                          double *lami, double *lama)
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
  advection1_eigenval1d(mesh,nvars, lami,norm);
  lama[0] = lami[0]; // eigenval is same on both sides

  /* get inner and adjacent fluxes fi, fa */
  advection1_flux1d(mesh,nvars, fi,norm, ui);
  advection1_flux1d(mesh,nvars, fa,norm, ua);
}


/* set a BC on patch boundary */
void advection1_u_BC(tMesh *mesh, tVarList *vlr, tVarList *vlu)
{
  int ir = vlr->index[0];
  //int iu = vlu->index[0];
  int ix = Ind("x");
  char *advdir = Gets(Par("advection1_direction"));
  double nx,ny,nz, nmag2;
  int myid;

  /* prop. dir.*/
  sscanf(advdir, "%lg %lg %lg", &nx, &ny, &nz);
  nmag2 = (nx*nx + ny*ny + nz*nz);

  /* compute boundary flux terms */
  formylnodes(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
    tPat *pat = node->pat;
    int *n = node->n;
    double *r = Vard(node, ir);
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
            r[ijk] = -nmag2*cos(nx*x[ijk] + ny*y[ijk] + nz*z[ijk] - nmag2*t);
//printf("i,j,k: %d %d %d face%d nid%ld  ", i,j,k, face, node->nid);
//pr3v("norm",norm);printf("\n");
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
  int myid;

  TIMER_START;

  /* compute flux */
  advection1_f_df(mesh, vlu);

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

  /* get flux terms on boundary */
  dg_add_surface_fluxes(mesh, vlr, vlu,
                        advection1_fluxes_pt, advection1_numflux);

  /* impose outer BC */
  advection1_u_BC(mesh, vlr, vlu);

  TIMER_STOP;
  return 0;
}


/* initialize test */
int advection1_init(tMesh *mesh)
{
  int iu  = Ind("advection1_u");
  int ifx = Ind("advection1_fx");
  int ifxx = Ind("advection1_fxx");
  int iF  = Ind("advection1_F0");
  int ix =  Ind("x");
  int iue = Ind("advection1_u_err");
  tVarList *vlu = vlalloc(mesh);
  int numflux = Par("advection1_numflux");
  char *advdir = Gets(Par("advection1_direction"));
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
    forpoints(node, i) u[i] = sin(nx*x[i] + ny*y[i] + nz*z[i]);
  }

  /* register u and its RHS with evolve */
  evolve_register_subsys_u_rhs_src_lim(mesh, vlu, advection1_rhs_u, 0, 0,0);
  evolve_print_evosys(mesh);

  /* choose numerical flux */
  if(Getv(numflux, "LLF"))
    advection1_numflux = numflux1d_LLF;
  else
    advection1_numflux = numflux1d_upwind;

  return 0;
} 

/* calculate errors in u */
int advection1_analyze(tMesh *mesh)
{
  int iu  = Ind("advection1_u");
  int iue = Ind("advection1_u_err");
  int ix =  Ind("x");
  char *advdir = Gets(Par("advection1_direction"));
  double nx,ny,nz, nmag2;
  int myid;

  /* prop. dir.*/
  sscanf(advdir, "%lg %lg %lg", &nx, &ny, &nz);
  nmag2 = (nx*nx + ny*ny + nz*nz);

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
      double ua = sin(nx*x[i] + ny*y[i] + nz*z[i] - nmag2*t);
      ue[i] = fabs(u[i]- ua);
    }
  }
  return 0;
}
