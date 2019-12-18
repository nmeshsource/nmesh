/* advectionGhoFDy.c */
/* Wolfgang Tichy, 12/2019 */

#include "nmesh.h"
#include "advectionGhoFDy.h"

#define PR 1



/* frequently used pars */
tadvectionGhoFDy advectionGhoFDy[1];


/* func to init frequently used pars */
int advectionGhoFDy_init_global_pars(tMesh *mesh)
{
  advectionGhoFDy->sin_profile    = Getv(Par("advectionGhoFDy_profile"),"sin");
  advectionGhoFDy->square_profile = Getv(Par("advectionGhoFDy_profile"),"square");

  return 0;
}


/* RHS of: d_t u = - d_y u */
int advectionGhoFDy_vol_rhs_u(tMesh *mesh, tVarList *vlr, tVarList *vlu)
{
  int ir = Vind(vlr, 0);
  int iu = Vind(vlu, 0);

  /* RHS */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
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
  }

  return 0;
}

/* surface terms in RHS of: d_t u = -d_y u */
int advectionGhoFDy_surf_rhs_u(tMesh *mesh, tVarList *vlr, tVarList *vlu)
{
  int ir = Vind(vlr, 0);
  int iu = Vind(vlu, 0);
  int ix = Ind("advectionGhoFDy_x");
  int ngho = Geti(Par("amr_nghosts"));
  int n2gho = 2*ngho;

  /* RHS */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    tPat *pat = node->pat;
    int f = 2; /* face at lowest y */
    tBface *bfaces = pat->bfaces[f];

    /* compute boundary terms, if on outer bound */
    if(node->patface[f] && bfaces && bfaces->boundary)
    {
      double *u = Vard(node, iu);
      double *r = Vard(node, ir);
      double *x = Vard(node, ix);
      double *y = Vard(node, ix+1);
      double *z = Vard(node, ix+2);
      int *n = node->n;
      int dj = n[0];
      int i,j,k;
      //int in0 = n[0] - n2gho;
      int in1 = n[1] - n2gho;
      //int in2 = n[2] - n2gho;
      //double hx = (node->bbox[1] - node->bbox[0])/in0;
      double hy = (node->bbox[3] - node->bbox[2])/in1;
      //double hz = (node->bbox[5] - node->bbox[4])/in2;
      //double oohx2 = 1./(2.*hx);
      double oohy2 = 1./(2.*hy);
      //double oohz2 = 1./(2.*hz);

      /* RHS at each ghost point at outer boundary */
      for(k = 0; k < n[2]; k++)
      for(j = 0; j < 2;    j++)
      for(i = 0; i < n[0]; i++)
      {
        int ccc = Ind_n(i,j,k, n);
        int cmc = ccc - dj;
        double u_ccc = u[ccc];
        double u_cmc, u_cMc, uy;
        double t = mesh->time;
        double u1[1];
        double xyz[] = { x[ccc],y[ccc],z[ccc] };

        /* no neighbor ==> impose outer BC */
        /* set boundary values for u_cmc, u_cMc */
        if(j==0)
        {
          xyz[1] -= hy;
          advectionGhoFDy_set_profile_pt(xyz,t, 1, u1);
          u_cmc = u1[0];
          xyz[1] -= hy;
          advectionGhoFDy_set_profile_pt(xyz,t, 1, u1);
          u_cMc = u1[0];
        }
        else
        {
          xyz[1] -= 2.*hy;
          advectionGhoFDy_set_profile_pt(xyz,t, 1, u1);
          u_cmc = u[cmc];
          u_cMc = u1[0];
        }

        /* FD deriv */
        uy = (3.*u_ccc - 4.*u_cmc + u_cMc)*oohy2;

        r[ccc] = -uy;
      }
    }
  }

  return 0;
}


/* set specific x,y,z coords */
int advectionGhoFDy_set_coords(tMesh *mesh, int ix)
{
  int ngho = Geti(Par("amr_nghosts"));
  int n2gho = 2*ngho;

  /* coords */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    double *x = Vard(node, ix);
    double *y = Vard(node, ix+1);
    double *z = Vard(node, ix+2);
    int *n = node->n;
    int i,j,k;
    int in0 = n[0] - n2gho;
    int in1 = n[1] - n2gho;
    int in2 = n[2] - n2gho;
    double hx = (node->bbox[1] - node->bbox[0])/in0;
    double hy = (node->bbox[3] - node->bbox[2])/in1;
    double hz = (node->bbox[5] - node->bbox[4])/in2;

    for(k = 0; k < n[2]; k++)
    for(j = 0; j < n[1]; j++)
    for(i = 0; i < n[0]; i++)
    {
      int ccc = Ind_n(i,j,k, n);
      x[ccc] = node->bbox[0] + 0.5*hx + hx*(i - ngho);
      y[ccc] = node->bbox[2] + 0.5*hy + hy*(j - ngho);
      z[ccc] = node->bbox[4] + 0.5*hz + hz*(k - ngho);
    }
  }
  return 0;
}


/* set profile in var with index iphi */
void advectionGhoFDy_set_profile_pt(double xyz[3], double t, int nv, double *u)
{
  double nx = 0.;
  double ny = 1.;
  double nz = 0.;
  double nmag2 = (nx*nx + ny*ny + nz*nz);

  /* profile */
  if(advectionGhoFDy->sin_profile)
  {
    u[0] = sin(nx*xyz[0] + ny*xyz[1] + nz*xyz[2] - nmag2*t);
  }
  if(advectionGhoFDy->square_profile)
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
int advectionGhoFDy_set_profile(tMesh *mesh, int iu)
{
  int ix = Ind("advectionGhoFDy_x");
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

      advectionGhoFDy_set_profile_pt(xyz,t, 1, u1);
      u[i] = u1[0];
    }
  }
  return 0;
}


/* initialize advectionGhoFDy */
int advectionGhoFDy_init(tMesh *mesh)
{
  int iu  = Ind("advectionGhoFDy_u");
  int iue = Ind("advectionGhoFDy_u_err");
  int ix  = Ind("advectionGhoFDy_x");
  tVarList *vlu = vlalloc(mesh);

  PRF;printf(": dt = %g\n", mesh->dt);
  //printmesh(mesh);

  /* set global pars */
  advectionGhoFDy_init_global_pars(mesh);

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
  advectionGhoFDy_set_coords(mesh, ix);

  /* set initial profile, e.g. at t=0: set u=sin(x) */
  advectionGhoFDy_set_profile(mesh, iu);

  /* register u and its RHS with evolve */
  evolve_register_subsys_u_rhs_lim(mesh, vlu, advectionGhoFDy_vol_rhs_u,
                                   advectionGhoFDy_surf_rhs_u, 0,0);
  evolve_print_evosys(mesh);

  return 0;
}

/* calculate errors in u */
int advectionGhoFDy_analyze(tMesh *mesh)
{
  int iu  = Ind("advectionGhoFDy_u");
  int iue = Ind("advectionGhoFDy_u_err");

  if(PR) PRFs("\n");

  /* set correct profile in advectionGhoFDy_u_err */
  advectionGhoFDy_set_profile(mesh, iue);

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
