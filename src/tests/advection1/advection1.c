/* advection1.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "advection1.h"

#define PR 1



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
      double no[3] = { 0., 0., 0. };

      no[0] = 1.;
      advection1_flux1d(mesh,1, &(fx[i]),no, &u_i);
      no[0] = 0;

      no[1] = 1.;
      advection1_flux1d(mesh,1, &(fy[i]),no, &u_i);
      no[1] = 0.;

      no[2] = 1.;
      advection1_flux1d(mesh,1, &(fz[i]),no, &u_i);
      no[2] = 0.;
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
    double *u  = Vard(node, iu);
    double *fl[] = { Vard(node, ifx), Vard(node, ifx+1), Vard(node, ifx+2) };
    double fln, FN, lam, norm[3];
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

        /* flux times normal vector */
        fln = (norm[0]*fl[0][ijk] + norm[1]*fl[1][ijk] + norm[2]*fl[2][ijk]);

        /* eigenval in dir norm */
        advection1_eigenval1d(mesh,1, &lam,norm);

        /* if stuff is coming in */
        if(lam < 0.)
        {
          /* if there is an adjacent surface */
          if(uaj)
          {
            advection1_flux1d(mesh,1, &FN,norm, &(uaj[JK]));
          }
          else
          {
            FN = 0;
          }
        }
        else
        {
          FN = fln;
        }

        /* set physical fluxes on left and right */
/*
        cdir = 0;
        if(norm[cdir]>=0.)
        {
          uL = u[ijk];
          fL = fl[cdir][ijk];
          advection1_eigenval1d(mesh,1, &lamL,cdir);

          uR = uaj[JK];
          advection1_flux1d(mesh,1, &fR,cdir, &uR);
          advection1_eigenval1d(mesh,1, &lamR,cdir);
        }
        else
        {
          uR = u[ijk];
          fR = fl[cdir][ijk];
          advection1_eigenval1d(mesh,1, &lamR,cdir);

          uL = uaj[JK];
          advection1_flux1d(mesh,1, &fL,cdir, &uL);
          advection1_eigenval1d(mesh,1, &lamL,cdir);
        }

        numflux1d_LLF(mesh,1, &(FN[cdir]), &uL, &uR, &fL, &fR, &lamL, &lamR);


*/



        F[JK] = FN - fln;
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
