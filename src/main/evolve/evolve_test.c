/* evolve.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "evolve.h"

#define PR 0


/* RHS of: d_t u = -u */
int evolve_test_rhs_u(tMesh *mesh, tVarList *vlr, tVarList *vlu)
{
  /* get surfaces so that we can compute fluxes */
  get_all_myln_surfaces(mesh);

  /* RHS */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    double *u = Vard(node, vlu->index[0]);
    double *r = Vard(node, vlr->index[0]);
    int i;

    /* RHS at each point */
    forpoints(node, i) r[i] = -u[i];
  }
  return 0;
}

/* RHS of: d_t v = s */
int evolve_test_rhs_v(tMesh *mesh, tVarList *vlr, tVarList *vlv)
{
  /* get surfaces so that we can compute fluxes */
  get_all_myln_surfaces(mesh);

  /* RHS */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    //double *v = Vard(node, vlv->index[0]);
    double *r = Vard(node, vlr->index[0]);
    double *s = Vard(node, Ind("evolve_test_s"));
    int i;

    /* RHS at each point */
    forpoints(node, i) r[i] = s[i];
  }
  return 0;
}

/* set source s for v: s = u */
int evolve_test_src_u(tMesh *mesh, tVarList *vlr, tVarList *vlu)
{
  /* get surfaces so that we can compute fluxes */
  get_all_myln_surfaces(mesh);

  /* Source */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    double *u = Vard(node, vlu->index[0]);
    double *s = Vard(node, Ind("evolve_test_s"));
    int i;

    forpoints(node, i) s[i] = u[i];
  }
  return 0;
}


/* initialize test */
int evolve_test_init(tMesh *mesh)
{
  int iu = Ind("evolve_test_u");
  int iv = Ind("evolve_test_v");
  int is = Ind("evolve_test_s");
  int iue = Ind("evolve_test_u_err");
  int ive = Ind("evolve_test_v_err");
  tVarList *vlu = vlalloc(mesh);
  tVarList *vlv = vlalloc(mesh);

  PRF;printf(": dt = %g\n", mesh->dt);

  /* two varlists */
  vlpush(vlu, iu);
  vlpush(vlv, iv);

  /* enable all needed vars */
  enablevar(mesh, iu);
  enablevar(mesh, iv);
  enablevar(mesh, is);
  enablevar(mesh, iue);
  enablevar(mesh, ive);

  /* at t=0: set u=1, leave v=0 */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    double *u = Vard(node, iu);
    int i;
    forpoints(node, i) u[i] = 1.;
  }

  /* register u,v and their RHS with evolve */
  evolve_register_subsys_u_rhs_src_lim(mesh, vlu,
                                       evolve_test_rhs_u, evolve_test_src_u,
                                       NULL, NULL);
  evolve_register_subsys_u_rhs_src_lim(mesh, vlv,
                                       evolve_test_rhs_v, NULL, NULL, NULL);
  evolve_print_evosys(mesh);
  return 0;
} 

/* calculate errors in u and v */
int evolve_test_analyze(tMesh *mesh)
{
  int iu = Ind("evolve_test_u");
  int iv = Ind("evolve_test_v");
  int iue = Ind("evolve_test_u_err");
  int ive = Ind("evolve_test_v_err");

  if(PR) PRFs("\n");

  /*  compute errors */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    double *u = Vard(node, iu);
    double *v = Vard(node, iv);
    double *ue = Vard(node, iue);
    double *ve = Vard(node, ive);
    double t = mesh->time;
    double ua = exp(-t);       /* analytic soln for u */
    double va = 1. - exp(-t);  /* analytic soln for v */
    double c;
    int i;

    if(node->dt>0.)
      c = pow(node->dt, Getd(Par("evolve_method_order")));
     else
      c = pow(mesh->dt, Getd(Par("evolve_method_order")));

    forpoints(node, i)
    {
      ue[i] = fabs(u[i] - ua)/c; //fabs(u[i]/ua - 1.)/c;
      ve[i] = fabs(v[i] - va)/c; //fabs(v[i]/va - 1.)/c;
      //printf("c=%g u[i]=%g ua=%g ue[i]=%g\n", c, u[i], ua, ue[i]);
    }
  }
  return 0;
}
