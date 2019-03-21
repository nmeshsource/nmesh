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

//printvar_innode(node, iu);

    /* set F on each face */
    for(face=0; face<6; face++)
    {
      double *F = Vard(node, iF+face);
      double *uaj = Varaj(node, iu, face);
tArray *aua  = alloc_array2d(9,9);
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

double nmag = (nx*nx + ny*ny + nz*nz);
int ix = Ind("x");
double *x = Vard(node, ix);
double *y = Vard(node, ix+1);
double *z = Vard(node, ix+2);
double t = mesh->time;
double ua = sin(nx*x[ijk] + ny*y[ijk] + nz*z[ijk] - nmag*t);
aua->d[JK] = ua;
//FNx = ua * nx;
//FNy = ua * ny;
//FNz = ua * nz;
if(0)// && i==0 && k==0)
{
printf("i,j,k: %d %d %d face%d nid%ld  ", i,j,k, face, node->nid);
pr3v("norm",norm);printf("\n");
//printvar_innode(node, iu);
//printvar_innode(node->fnb[face][0], iu);
//printvar_ajsurfdiff(node, iu);
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
if(face==3 && node->nid==0 && VarAaj(node, iu, face))
{
double rdif = Lp_norm_array_diff(VarAaj(node, iu, face), aua, 2.);
printf("rdif=%g\n", rdif);
if(0 && rdif>0.2)
{
array_diff(aua, VarAaj(node, iu, face),aua);
//printarray(VarAaj(node, iu, face));
printarray(aua);
rdif = Lp_norm_array(aua, 2);
printf("rdif2=%g\n", rdif);
}
}
free_array(aua);
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
  double nx,ny,nz, nmag;
  int myid;

  /* prop. dir.*/
  sscanf(advdir, "%lg %lg %lg", &nx, &ny, &nz);
  nmag = (nx*nx + ny*ny + nz*nz);

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
            r[ijk] = -nmag*cos(nx*x[ijk] + ny*y[ijk] + nz*z[ijk] - nmag*t);
//printf("i,j,k: %d %d %d face%d nid%ld  ", i,j,k, face, node->nid);
//pr3v("norm",norm);printf("\n");
          }
        }
    }
  }
}

/* RHS of: d_t u = - d_i f^i */
void advection1_rhs_u(tMesh *mesh, tVarList *vlr, tVarList *vlu)
{
  int ir = vlr->index[0];
  //int iu = vlu->index[0];
  int ifxx = Ind("advection1_fxx");
  int iF   = Ind("advection1_F0");
  int iooJ = Ind("det_dXbdx");
  tArray *a2J = alloc_array2d(3,2);   /* 3x2 for 2-Jacobian */
  tArray *a2gam = alloc_array2d(2,2); /* 2x2 for induced metric on surface */
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
    //double dXbdX[3];
    int face;

    ///* get dXb/dX needed for induced 2-metric 2gam */
    //dXbYbZb_dXYZ(node, dXbdX);

    for(face=0; face<6; face++)
    {
      int dir = face/2;
      int p = (face%2)*(n[dir] - 1);
      //double sig = 2*(face%2) - 1;
      double *F = Vard(node, iF+face);
      double *w = Wquad(node, dir);
      double sqrtdet2gam, det2dxdXb;
      int i,j,k, ijk, JK, i0;

      forplaneN(dir, i,j,k, n, p)
      {
        ijk = Ind_n(i,j,k, n);
        JK = Ind_n_norm(i,j,k, n, dir);
        i0 = i0_norm(i,j,k, dir);

        /* get 2d Jacobian of dx/dXb */
        array_2dxdXb(node, ijk, dir, a2J);
        /* compute 2-metric from a2J^T a2J */
        mm_array0_norestrict(a2J,a2J, a2gam); // for now assume gam = flat
        det2dxdXb = det_2_2array(a2gam);
        sqrtdet2gam = sqrt(det2dxdXb);
/*
if(face==1 && j==1 && k==2 && myid==0)
{
printf("nid%ld l%d i=%d ", node->nid, node->l, i);
printarray_matrix0(a2J);
printarray_matrix0(a2gam);


printnode(node);
double dXbdX[3];
dXbYbZb_dXYZ(node, dXbdX);
for(int l=0; l<3; l++) printf("WWWWW %g\n", dXbdX[l]);


double det_dXbYbZb_dXYZ = dXbdX[0] * dXbdX[1] * dXbdX[2];
printf("det_dXbYbZb_dXYZ=%g\n", det_dXbYbZb_dXYZ);

printf("det_dXbdx = ooJ[ijk]=%g\n", ooJ[ijk]);

double norm[3];
node_normal_at_ijk(node, face, ijk, norm);
printf("i,j,k: %d %d %d face%d norm[0]=%g\n", i,j,k, face, norm[0]);

exit(8);
//printf("i0=%d w[i0]=%g\n", i0, w[i0]);
//printf("sqrtdet2gam=%g\n", sqrtdet2gam);
//if(face==3 && node->nid==0) printf("ooJ[ijk]=%g\n", ooJ[ijk]);
}
*/
        r[ijk] -= F[JK] * sqrtdet2gam * fabs(ooJ[ijk])/ w[i0];
      }
    }
  }

  /* impose outer BC */
  advection1_u_BC(mesh, vlr, vlu);

  /* free arrays */
  free_array(a2gam);
  free_array(a2J);

  TIMER_STOP;
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
  evolve_register_subsys_u_rhs_src(mesh, vlu, advection1_rhs_u, 0);
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
  double nx,ny,nz, nmag;
  int myid;

  /* prop. dir.*/
  sscanf(advdir, "%lg %lg %lg", &nx, &ny, &nz);
  nmag = (nx*nx + ny*ny + nz*nz);

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
      double ua = sin(nx*x[i] + ny*y[i] + nz*z[i] - nmag*t);
      ue[i] = fabs(u[i]- ua);
    }
  }
  return 0;
}
