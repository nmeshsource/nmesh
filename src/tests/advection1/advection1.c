/* advection1.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "advection1.h"

#define PR 1



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
  char *advdir = Gets(Par("advection1_direction"));
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
      fx[i] = nx*u_i;
      fy[i] = ny*u_i;
      fz[i] = nz*u_i;
    }

   /* flux derivs */
   cart_partials(node, ifx, ifxx);
   cart_partials(node, ify, ifyx);
   cart_partials(node, ifz, ifzx);
  }
}


/* use numerical flux FN^i to set F */
void advection1_F(tMesh *mesh, tVarList *vlu)
{
  int iu = vlu->index[0];
  int ifx = Ind("advection1_fx");
  int iF  = Ind("advection1_F0");
  char *advdir = Gets(Par("advection1_direction"));
  double nx,ny,nz;
  int myid;

  /* prop. dir.*/
  sscanf(advdir, "%lg %lg %lg", &nx, &ny, &nz);

  /* compute boundary flux terms */
  formylnodes(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
    int *n = node->n;
    double *fx = Vard(node, ifx);
    double *fy = Vard(node, ifx+1);
    double *fz = Vard(node, ifx+2);
    double FNx,FNy,FNz, norm[3];
    int face, dir, p, i,j,k, ijk, JK;

    /* set F on each face */
    for(face=0; face<6; face++)
    {
      double *F = Vard(node, iF+face);
      double *uaj = Varaj(node, iu, face);
      dir = face/2;
      p = (face%2)*(n[dir] - 1);
      forplaneN(dir, i,j,k, n, p)
      {
        ijk = Ind_n(i,j,k, n);
        JK = Ind_n_norm(i,j,k, n, dir);
        node_normal_at_ijk(node, face, ijk, norm);

        /* if stuff is coming in */
        if(norm[0]*nx + norm[1]*ny + norm[2]*nz < 0.)
        {
          /* if there is an adjacent surface */
          if(uaj)
          {
            FNx = uaj[JK] * nx;
            FNy = uaj[JK] * ny;
            FNz = uaj[JK] * nz;
//var_finite(node, iu);
if(0) if(!isfinite(uaj[JK]))
{
var_finite(node, iu);
printf("uaj = %p JK = %d uaj[JK] = %g\n", uaj, JK, uaj[JK]);
printf("face%d dir%d p%d %d %d %d: %d\n", face, dir, p, i,j,k, ijk);
printarray( VarAaj(node, iu, face) );
errorexit("nan uaj");
}
          }
          else
          {
            FNx = FNy = FNz = 0;
          }
        }
        else
        {
          FNx = fx[ijk];
          FNy = fy[ijk];
          FNz = fz[ijk];
        }

        F[JK] = (FNx - fx[ijk])*norm[0] +
                (FNy - fy[ijk])*norm[1] +
                (FNz - fz[ijk])*norm[2];
      }
    }
  }
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
  int iF   = Ind("advection1_F0");
  int iooJ = Ind("det_dXbdx");
  int isqrtdet2gamma0 = Ind("sqrtdet2gamma0");
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

  /* get surfaces so that we can compute fluxes */
  get_all_myln_surfaces(mesh);

  /* get flux terms on boundary */
  advection1_F(mesh, vlu);

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
