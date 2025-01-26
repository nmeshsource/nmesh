/* advectionFDy.c */
/* Wolfgang Tichy, 12/2019 */

#include "nmesh.h"
#include "advectionFDy.h"

#define PR 1



/* frequently used pars */
tadvectionFDy advectionFDy[1];


/* func to init frequently used pars */
int advectionFDy_init_global_pars(tMesh *mesh)
{
  advectionFDy->sin_profile    = Getv(Par("advectionFDy_profile"),"sin");
  advectionFDy->square_profile = Getv(Par("advectionFDy_profile"),"square");

  return 0;
}


/* RHS of: d_t u = - d_y u */
int advectionFDy_vol_rhs_u(tNode *node, tEvoVars *evv)
{
  tVarList *vlr = EvoVars_vlr(evv);
  tVarList *vlu = EvoVars_vlu(evv);
  int ir = Vind(vlr, 0);
  int iu = Vind(vlu, 0);

  double *u  = Vard(node, iu);
  double *r  = Vard(node, ir);
  int *n = node->n;
  int dj = n[0];
  int i,j,k;
  //double hx = (node->bbox[1] - node->bbox[0])/n[0];
  double hy = (node->bbox[3] - node->bbox[2])/n[1];
  //double hz = (node->bbox[5] - node->bbox[4])/n[2];
  //double oohx2 = 1./(2.*hx);
  double oohy2 = 1./(2.*hy);
  //double oohz2 = 1./(2.*hz);

  /* RHS at each point in interior */
  for(k = 0; k < n[2]; k++)
  for(j = 2; j < n[1]; j++)
  for(i = 0; i < n[0]; i++)
  {
    int ccc = Ind_n(i,j,k, n);
    int cmc = ccc - dj;
    int cMc = cmc - dj;
    double uy = (3.*u[ccc] - 4.*u[cmc] + u[cMc])*oohy2;  /* FD deriv */

    r[ccc] = -uy;
  }

  return 0;
}

/* surface terms in RHS of: d_t u = -d_y u */
int advectionFDy_surf_rhs_u(tNode *node, tEvoVars *evv)
{
  tVarList *vlr = EvoVars_vlr(evv);
  tVarList *vlu = EvoVars_vlu(evv);
  tMesh *mesh = vlu->mesh;
  int ir = Vind(vlr, 0);
  int iu = Vind(vlu, 0);
  int ix = Ind("advectionFDy_x");

  double *u = Vard(node, iu);
  double *r = Vard(node, ir);
  double *x = Vard(node, ix);
  double *y = Vard(node, ix+1);
  double *z = Vard(node, ix+2);
  int *n = node->n;
  int dj = n[0];
  int i,j,k;
  //double hx = (node->bbox[1] - node->bbox[0])/n[0];
  double hy = (node->bbox[3] - node->bbox[2])/n[1];
  //double hz = (node->bbox[5] - node->bbox[4])/n[2];
  //double oohx2 = 1./(2.*hx);
  double oohy2 = 1./(2.*hy);
  //double oohz2 = 1./(2.*hz);
  tArray *Au_aj[] = { VarAaj(node, iu, 0), VarAaj(node, iu, 1),
                      VarAaj(node, iu, 2), VarAaj(node, iu, 3),
                      VarAaj(node, iu, 4), VarAaj(node, iu, 5) };

  /* RHS at each point at y-surface */
  for(k = 0; k < n[2]; k++)
  for(j = 0; j < 2;    j++)
  for(i = 0; i < n[0]; i++)
  {
    int ccc = Ind_n(i,j,k, n);
    int cmc = ccc - dj;
    double u_ccc = u[ccc];
    double u_cmc, u_cMc, uy;


    if(Au_aj[2]) /* there is a neighbor box */
    {
      double *uaj = Au_aj[2]->d;
      int    *naj = Au_aj[2]->n;

      if(j==1)
      {
        u_cmc = u[cmc];
        u_cMc = uaj[Ind_n(i,1,k, naj)];
      }
      else
      {
        u_cmc = uaj[Ind_n(i,1,k, naj)];
        u_cMc = uaj[Ind_n(i,0,k, naj)];
      }
    }
    else /* no neighbor */
    {
      /* impose outer BC */
      double t = mesh->time;
      double u1[1];
      double xyz[] = { x[ccc],y[ccc],z[ccc] };

      /* set boundary values for u_cmc, u_cMc */
      if(j==1)
      {
        xyz[1] -= 2.*hy;
        advectionFDy_set_profile_pt(xyz,t, 1, u1);
        u_cmc = u[cmc];
        u_cMc = u1[0];
      }
      else
      {
        xyz[1] -= hy;
        advectionFDy_set_profile_pt(xyz,t, 1, u1);
        u_cmc = u1[0];
        xyz[1] -= hy;
        advectionFDy_set_profile_pt(xyz,t, 1, u1);
        u_cMc = u1[0];
      }
    }

    /* FD deriv */
    uy = (3.*u_ccc - 4.*u_cmc + u_cMc)*oohy2;

    r[ccc] = -uy;
  }

  return 0;
}


/* set specifuc x,y,z coords */
int advectionFDy_set_coords(tMesh *mesh, int ix)
{
  /* coords */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    double *x = Vard(node, ix);
    double *y = Vard(node, ix+1);
    double *z = Vard(node, ix+2);
    int *n = node->n;
    int i,j,k;
    double hx = (node->bbox[1] - node->bbox[0])/n[0];
    double hy = (node->bbox[3] - node->bbox[2])/n[1];
    double hz = (node->bbox[5] - node->bbox[4])/n[2];

    for(k = 0; k < n[2]; k++)
    for(j = 0; j < n[1]; j++)
    for(i = 0; i < n[0]; i++)
    {
      int ccc = Ind_n(i,j,k, n);
      x[ccc] = node->bbox[0] + 0.5*hx + hx*i;
      y[ccc] = node->bbox[2] + 0.5*hy + hy*j;
      z[ccc] = node->bbox[4] + 0.5*hz + hz*k;
    }
  }
  return 0;
}


/* set profile in var with index iphi */
void advectionFDy_set_profile_pt(double xyz[3], double t, int nv, double *u)
{
  double nx = 0.;
  double ny = 1.;
  double nz = 0.;
  double nmag2 = (nx*nx + ny*ny + nz*nz);

  /* profile */
  if(advectionFDy->sin_profile)
  {
    u[0] = sin(nx*xyz[0] + ny*xyz[1] + nz*xyz[2] - nmag2*t);
  }
  if(advectionFDy->square_profile)
  {
    double inx, iny;
    if(xyz[0]>=(-0.7 + nx*t) && xyz[0]<=(-0.3 + nx*t)) inx = 1.;
    else                                               inx = 0.;
    if(xyz[1]>=(-0.2 + ny*t) && xyz[1]<=(+0.2 + ny*t)) iny = 1.;
    else                                               iny = 0.;
    u[0] = inx*iny;
  }
}

/* set profile in var with index iphi */
int advectionFDy_set_profile(tMesh *mesh, int iu)
{
  int ix = Ind("advectionFDy_x");
  double t = mesh->time;
  double u1[1] = {0.};

  /* profile */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    double *x = Vard(node, ix);
    double *y = Vard(node, ix+1);
    double *z = Vard(node, ix+2);
    double *u = Vard(node, iu);
    int i;

    forpoints(node, i)
    {
      double xyz[] = { x[i],y[i],z[i] };

      advectionFDy_set_profile_pt(xyz,t, 1, u1);
      u[i] = u1[0];
    }
  }
  return 0;
}


/* initialize advectionFDy */
int advectionFDy_init(tMesh *mesh)
{
  int iu  = Ind("advectionFDy_u");
  int iue = Ind("advectionFDy_u_err");
  int ix  = Ind("advectionFDy_x");
  tVarList *vlu = vlalloc(mesh);

  PRF;printf(": dt = %g\n", mesh->dt);
  //printmesh(mesh);

  /* set global pars */
  advectionFDy_init_global_pars(mesh);

  /* set surface to thickness 2 */
  MeshVarSetSurfInfo(mesh, iu, 2);

  /* varlist */
  vlpush(vlu, iu);

  /* enable all needed vars */
  enablevar(mesh, iu);
  enablevar(mesh, iue);
  enablevar(mesh, ix);
  enablevar(mesh, ix+1);
  enablevar(mesh, ix+2);

  /* init coords */
  advectionFDy_set_coords(mesh, ix);

  /* set initial profile, e.g. at t=0: set u=sin(x) */
  advectionFDy_set_profile(mesh, iu);

  /* register u and its RHS with evolve */
  evolve_register_subsys_u_rhs_lim(mesh, vlu, advectionFDy_vol_rhs_u,
                                   advectionFDy_surf_rhs_u, 0,0);
  evolve_print_evosys(mesh);

  return 0;
}

/* calculate errors in u */
int advectionFDy_analyze(tMesh *mesh)
{
  int iu  = Ind("advectionFDy_u");
  int iue = Ind("advectionFDy_u_err");

  if(PR) PRFs("\n");

  /* set correct profile in advectionFDy_u_err */
  advectionFDy_set_profile(mesh, iue);

  /*  compute errors: u_err = u - u_correct */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    double *u  = Vard(node, iu);
    double *ue = Vard(node, iue);
    int i;

    forpoints(node, i)
      ue[i] = u[i]- ue[i];
  }
  return 0;
}
