/* evolve.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "evolve.h"

#define PR 0


/* RHS of: d_t u = -u */
void evolve_test_rhs_u(tNode *node, tVarList *vlr, tVarList *vlu)
{
  //tMesh *mesh = node->pat->mesh;
  double *u = Vard(node, vlu->index[0]);
  double *r = Vard(node, vlr->index[0]);
  int i;

  /* get surfaces so that we can compute fluxes */
  get_all_surfaces(node);
  free_dat_reqs_after_Waitall_com_send(node);

  /* RHS at each point */
  forpoints(node, i) r[i] = -u[i];
}

/* RHS of: d_t v = s */
void evolve_test_rhs_v(tNode *node, tVarList *vlr, tVarList *vlv)
{
  tMesh *mesh = node->pat->mesh;
  //double *v = Vard(node, vlu->index[0]);
  double *r = Vard(node, vlr->index[0]);
  double *s = Vard(node, Ind("evolve_test_s"));
  int i;

  /* get surfaces so that we can compute fluxes */
  get_all_surfaces(node);
  free_dat_reqs_after_Waitall_com_send(node);

  /* RHS at each point */
  forpoints(node, i) r[i] = s[i];
}

/* set source s for v: s = u */
void evolve_test_src_u(tNode *node, tVarList *vlu)
{
  tMesh *mesh = node->pat->mesh;
  double *u = Vard(node, vlu->index[0]);
  double *s = Vard(node, Ind("evolve_test_s"));
  int i;
  forpoints(node, i) s[i] = u[i];
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
  int myid;

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
  formylnodes(mesh, myid)
  {
    tNode *node = GetMyNode(mesh, myid);
    double *u = Vard(node, iu);
    int i;
    forpoints(node, i) u[i] = 1.;
  }

  /* register u,v and their RHS with evolve */
  evolve_register_subsys_u_rhs_src(mesh, vlu,
                                   evolve_test_rhs_u, evolve_test_src_u);
  evolve_register_subsys_u_rhs_src(mesh, vlv,
                                   evolve_test_rhs_v, NULL);
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
  int myid;

  if(PR) PRFs("\n");

  /*  compute errors */
  formylnodes(mesh, myid)
  {
    tNode *node = GetMyNode(mesh, myid);
    double *u = Vard(node, iu);
    double *v = Vard(node, iv);
    double *ue = Vard(node, iue);
    double *ve = Vard(node, ive);
    double c = pow(node->dt, Getd(Par("evolve_method_order")));
    double t = mesh->time;
    double ua = exp(-t);       /* analytic soln for u */
    double va = 1. - exp(-t);  /* analytic soln for v */
    int i;

    forpoints(node, i)
    {
      ue[i] = fabs(u[i]/ua - 1.)/c; //fabs(u[i] - ua)/c;
      ve[i] = fabs(v[i]/va - 1.)/c; //fabs(v[i] - va)/c;
    }
  }
  return 0;
}
