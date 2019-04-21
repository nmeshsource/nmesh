/* advection1.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "advection1.h"

#define PR 1


/* func pointer for numerical flux */
void (*adevection1_numflux)(tMesh *mesh, int nf, double *fnum,
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


/* function that sets flux array lists and eigenval array lists on both
   sides of the node surfaces */
void advection1_fluxarrays(tNode *node, tVarList *vlu,
                           pArrList *Afn[6],  pArrList *Afnaj[6],
                           pArrList *Alam[6], pArrList *Alamaj[6])
{
  tMesh *mesh = node->pat->mesh;
  int *n = node->n;
  int nvars = vlu->n;
  double *ui   = dmalloc(nvars);
  double *fi   = dmalloc(nvars);
  double *lami = dmalloc(nvars);
  double *ua   = dmalloc(nvars);
  double *fa   = dmalloc(nvars);
  double *lama = dmalloc(nvars);
  int ifx = Ind("advection1_fx");
  int face;

  /* set fluxes on both sides of each face */
  for(face=0; face<6; face++)
  {
    int dir = face/2;
    int p = (face%2)*(n[dir] - 1);
    int i,j,k;

    forplaneN(dir, i,j,k, n, p)
    {
      int ijk = Ind_n(i,j,k, n);
      int JK = Ind_n_norm(i,j,k, n, dir);
      double norm[3];
      int l;

      /* outward pointing normal */
      node_normal_at_ijk(node, face, ijk, norm);

      /* loop over vars and set ui, ua, as well as fi */
      forvl(vlu, l)
      {
        double *u = Vard_(node, l);
        double *uaj = Varaj_(node, l, face);
        double *f[] = { Vard_(node, ifx),
                        Vard_(node, ifx+1), Vard_(node, ifx+2) };

        /* cons var and normal flux inside node */
        ui[l] = u[ijk];
        fi[l] = ( norm[0]*f[0][ijk] + norm[1]*f[1][ijk] + norm[2]*f[2][ijk] );

        /* cons var on adjacent side */
        ua[l] = uaj[JK];
      }

      /* eigenval in dir norm */
      advection1_eigenval1d(mesh,nvars, lami,norm);
      lama[0] = lami[0]; // eigenval is same on both sides

      /* get adjcent flux fa. fi was set earlier */
      advection1_flux1d(mesh,nvars, fa,norm, ua);

      forvl(vlu, l)
      {
        tArray *afn    = ListEntry(Afn[face], l);
        tArray *afnaj  = ListEntry(Afnaj[face], l);
        tArray *alam   = ListEntry(Alam[face], l);
        tArray *alamaj = ListEntry(Alamaj[face], l);
        double *fn    = Arrd_(afn);
        double *fnaj  = Arrd_(afnaj);
        double *lam   = Arrd_(alam);
        double *lamaj = Arrd_(alamaj);

        /* set flux arrays */
        fn[JK]   = fi[l];
        fnaj[JK] = fa[l];

        /* set eigenval arrays */
        lam[JK]   = lami[l];
        lamaj[JK] = lama[l];
      }

    } /* end plane loop */
  }

  free(lama);
  free(fa);
  free(ua);
  free(lami);
  free(fi);
  free(ui);
}


/* function that sets fluxes and eigenvals on both sides of a node surface.
   In: ui, ua, norm. Out: fi,fa, lami,lama */
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
    double fln, FN, norm[3];
    double uL, uR, fL, fR, lamL, lamR;
    int face, dir, plus, p, i,j,k, ijk, JK;

    /* set F on each face */
    for(face=0; face<6; face++)
    {
      double *F = Vard(node, iF+face);
      double *uaj = Varaj(node, iu, face);
      dir = face/2;
      plus = (face%2);
      p = plus*(n[dir] - 1);
      forplaneN(dir, i,j,k, n, p)
      {
        ijk = Ind_n(i,j,k, n);
        JK = Ind_n_norm(i,j,k, n, dir);
        node_normal_at_ijk(node, face, ijk, norm);

        /* flux times normal vector */
        fln = (norm[0]*fl[0][ijk] + norm[1]*fl[1][ijk] + norm[2]*fl[2][ijk]);

if(1)
{
        /* eigenval in dir norm */
        advection1_eigenval1d(mesh,1, &lamL,norm);

        /* if stuff is coming in */
        if(lamL < 0.)
        {
          /* if there is an adjacent surface */
          if(uaj)
          {
            advection1_flux1d(mesh,1, &FN,norm, &(uaj[JK]));
          }
          else
          {
            FN = 0.;
          }
        }
        else
        {
          FN = fln;
        }
}
else
{
        /* set physical fluxes and eigenvalues on left and right */
        uL = u[ijk];
        fL = fln;
        advection1_eigenval1d(mesh,1, &lamL,norm);
        /* if there is a neighbor and an adjacent surface */
        if(uaj)
        {
          uR = uaj[JK];
          advection1_flux1d(mesh,1, &fR,norm, &uR);
          advection1_eigenval1d(mesh,1, &lamR,norm);

          // use LLF
          numflux1d_LLF(mesh,1, &FN, &uL, &uR, &fL, &fR, &lamL, &lamR);
        }
        else /* do something special on outer boundary */
        {
          if(lamL < 0.)
            FN = 0.;
          else
            FN = fln;
        }
}

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

if(0)
{
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
}
else
{
  dg_add_surface_fluxes(mesh, vlr, vlu,
                        advection1_fluxes_pt, adevection1_numflux);
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

  /* choose numerical flux */
  adevection1_numflux = numflux1d_LLF;

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
