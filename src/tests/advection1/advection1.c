/* advection1.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "advection1.h"

#define PR 1


/* RHS of: d_t u = - n_i d_i u */
void advection1_rhs_u(tMesh *mesh, tVarList *vlr, tVarList *vlu)
{
  int ir = vlr->index[0];
  int iu = vlu->index[0];
  int iux = Ind("advection1_ux");
  char *advdir = Gets(Par("advection1_direction"));
  double nx,ny,nz;
  int myid;

  /* prop. dir.*/
  sscanf(advdir, "%lg %lg %lg", &nx, &ny, &nz);

  /* compute derivs */
  formylnodes(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
    cart_partials(node, iu, iux);
  }

  /* get surfaces so that we can compute fluxes */
  get_all_myln_surfaces(mesh);

  /* RHS */
  formylnodes(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
    double *r  = Vard(node, ir);
    double *ux = Vard(node, iux);
    double *uy = Vard(node, iux+1);
    double *uz = Vard(node, iux+2);
    int i;

    /* RHS at each point */
    forpoints(node, i) r[i] = -(nx*ux[i] + ny*uy[i] + nz*uz[i]);
  }
}


/* initialize test */
int advection1_init(tMesh *mesh)
{
  int iu  = Ind("advection1_u");
  int iux = Ind("advection1_ux");
  int ix =  Ind("x");
  int iue = Ind("advection1_u_err");
  tVarList *vlu = vlalloc(mesh);
  int myid;

  PRF;printf(": dt = %g\n", mesh->dt);

  /* varlist */
  vlpush(vlu, iu);

  /* enable all needed vars */
  enablevar(mesh, iu);
  enablevar(mesh, iux);
  enablevar(mesh, iue);

  /* at t=0: set u=sin(x) */
  formylnodes(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
    double *u = Vard(node, iu);
    double *x = Vard(node, ix);
    int i;
    forpoints(node, i) u[i] = sin(x[i]);
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
  int myid;

  if(PR) PRFs("\n");

  /*  compute errors */
  formylnodes(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
    double *u = Vard(node, iu);
    double *ue = Vard(node, iue);
    double *x = Vard(node, ix);
    double t = mesh->time;
    int i;

    forpoints(node, i)
    {
      double ua = sin(x[i]-t);
      ue[i] = fabs(u[i]/ua - 1.);
    }
  }
  return 0;
}
