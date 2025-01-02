/* dissipation.c */
/* Wolfgang Tichy, May 2021 */


#include "nmesh.h"
#include "dg.h"


/***********************************************************************/
/* funcs needed for dissipation */
/***********************************************************************/

/* Add Kreiss-Oliger dissipation terms to vlr. vlr can be varlist for RHS.
   Here we use the same Kreiss-Oliger 4th order derivative operator as in
   BAM ( https://arxiv.org/pdf/gr-qc/0610128.pdf ):
   Note:
   D_{+}^2 D_{-}^2 u -> (u[i-2] - 4u[i-1] + 6u[i] - 4u[i+1] + u[i+2])/h^4
   BAM's d4stencil2nd in dissipation4 thus has { 1, -4,  6, -4, 1}
   The diss. term in Eq. (63) for rho=1 and r=2 is:
     \sigma (-h)^3 D_{+}^2 D_{-}^2 / 16 u
     = -(\sigma/16)*(u[i-2] - 4u[i-1] + 6u[i] - 4u[i+1] + u[i+2])/h
   So our dissfac below is given by dissfac = \sigma.
   The diss term converges to 0 at O(h^2) in the grid spacing h. It doesn't
   really fit into the DG or FV scheme, but since we are already doing
   several non-standard things, we can try this too.
   In:
     node
     vlu contains evolved fields
     dissfac is dissipation factor
   Out:
     vlr is the varlist to which we add dissipation terms */
void dissipation_add_KO4(tNode *node, tVarList *vlr, tVarList *vlu,
                         double dissfac)
{
  int *n = node->n;
  double *bb = node->bbox;
  int maxn = max3(n[0],n[1],n[2]);
  double *uc = dtensor(maxn);
  int dir;

  /* add dissipation in each direction to RHS */
  for(dir=0; dir<3; dir++)
  {
    double ooh = (n[dir]-1)/(bb[2*dir+1] - bb[2*dir]);// 1/dist betw. points
    double facoh = dissfac * ooh * 0.0625; /* dissfac/h * 1/16 */
    int i,j,k;

    /* do nothing if we have less than 5 grid points */
    if(n[dir]<5) continue;

    /* loop over plane */
    forplaneN(dir, i,j,k, n, 0)
    {
      int i1 = i1_norm(i,j,k, dir); /* 1st and 2nd index in plane */
      int i2 = i2_norm(i,j,k, dir);
      int i0;                       /* index orthogonal to plane */
      int ic,jc,kc, ccc;
      int l;                        /* field index */

      /* loop over fields */
      forvl(vlu, l)
      {
        double *ul = Vard(node, Vind(vlu, l)); /* field data pointer */
        double *rl = Vard(node, Vind(vlr, l)); /* RHS data pointer */

        /* fill field arrays uc, i0 runs orth. to plane */
        for(i0=0; i0<n[dir]; i0++)
        {
          /* set points and their index */
          ijk_inplaneN(dir, ic,jc,kc, i1,i2, i0);
          ccc = Ind_n(ic,jc,kc, n);
          /* set uc */
          uc[i0] = ul[ccc];
        }

        /* loop over inner points in dir */
        for(i0=2; i0<n[dir]-2; i0++)
        {
          ijk_inplaneN(dir, ic,jc,kc, i1,i2, i0);
          ccc = Ind_n(ic,jc,kc, n);

          /* add dissipation term to RHS */
          rl[ccc] += -facoh*( 6.*uc[i0] - 4.*(uc[i0-1] + uc[i0+1])
                                        +     uc[i0-2] + uc[i0+2] );
        }
      } /* end loop over fields */
    } /* end plane loop */
  } /* end dir-loop*/

  /* release mem */
  free(uc);
}

/* Add Kreiss-Oliger dissipation terms to vlr, which can be the RHS vlr */
void dissipation_add_KO4_mesh(tMesh *mesh, tVarList *vlr, tVarList *vlu,
                              double dissfac)
{
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    dissipation_add_KO4(node, vlr, vlu, dissfac);
  }
}



/* stencil weights for nth order diss operators on uniform grids */
double sw2[]  = {1, -2, 1}; //2nd order deriv stencil
double sw4[]  = {1, -4,  6, -4, 1}; //4th order deriv stencil
double sw6[]  = {1, -6, 15,-20, 15, -6, 1}; //6th order deriv stencil
double sw8[]  = {1, -8, 28, -56, 70, -56, 28, -8, 1};
double sw10[] = {1, -10, 45, -120, 210, -252, 210, -120, 45, -10, 1};
double sw12[] = {1, -12, 66, -220, 495, -792, 924, -792, 495, -220, 66, -12, 1};


/* Add Kreiss-Oliger dissipation terms to vlr. vlr can be varlist for RHS.
   Here we use the same Kreiss-Oliger nth order derivative operator as in
   BAM ( https://arxiv.org/pdf/gr-qc/0610128.pdf ).
   In:
     node
     vlu contains evolved fields
     dissfac is dissipation factor
     order is the order of the derivative operator we want (r=order/2)
   Out:
     vlr is the varlist to which we add dissipation terms */
void dissipation_add_KO_order(tNode *node, tVarList *vlr, tVarList *vlu,
                              double dissfac, int order)
{
  int *n = node->n;
  double *bb = node->bbox;
  int maxn = max3(n[0],n[1],n[2]);
  double *uc = dtensor(maxn);
  int dir;
  double *sw;                   /* stencil weights */
  double sgn = (-1. + (order%4))/(1 << order); /* (overall sign) / 2^order */
  int srad = order/2;           /* stencil radius */

  /* choose stencil weights */
  switch(order)
  {
  case 2:
    sw = sw2;
    break;
  case 4:
    sw = sw4;
    break;
  case 6:
    sw = sw6;
    break;
  case 8:
    sw = sw8;
    break;
  case 10:
    sw = sw10;
    break;
  case 12:
    sw = sw12;
    break;
  default:
    errorexit("order must be 2,4,6,8,10,12");
  }

  /* add dissipation in each direction to RHS */
  for(dir=0; dir<3; dir++)
  {
    double ooh = (n[dir]-1)/(bb[2*dir+1] - bb[2*dir]);// 1/dist betw. points
    double facoh = sgn * dissfac * ooh; /* (-1)^(1+order/2) dissfac/h */
    int i,j,k;

    /* do nothing if we have less than order+1 grid points */
    if(n[dir]<=order) continue;

    /* loop over plane */
    forplaneN(dir, i,j,k, n, 0)
    {
      int i1 = i1_norm(i,j,k, dir); /* 1st and 2nd index in plane */
      int i2 = i2_norm(i,j,k, dir);
      int i0;                       /* index orthogonal to plane */
      int ic,jc,kc, ccc;
      int l;                        /* field index */

      /* loop over fields */
      forvl(vlu, l)
      {
        double *ul = Vard(node, Vind(vlu, l)); /* field data pointer */
        double *rl = Vard(node, Vind(vlr, l)); /* RHS data pointer */

        /* fill field arrays uc, i0 runs orth. to plane */
        for(i0=0; i0<n[dir]; i0++)
        {
          /* set points and their index */
          ijk_inplaneN(dir, ic,jc,kc, i1,i2, i0);
          ccc = Ind_n(ic,jc,kc, n);
          /* set uc */
          uc[i0] = ul[ccc];
        }

        /* loop over inner points in dir */
        for(i0=srad; i0<n[dir]-srad; i0++)
        {
          double dis;
          int s;

          ijk_inplaneN(dir, ic,jc,kc, i1,i2, i0);
          ccc = Ind_n(ic,jc,kc, n);

          /* Set dissipation term. Here we assume the stencil is symmetric
             (i.e. same weight for uc[i0-s] and uc[i0+s]) and thus use only
             the top half of sw. */
          dis = sw[srad]*uc[i0];
          for(s=1; s<=srad; s++) dis += sw[s+srad]*(uc[i0-s] + uc[i0+s]);
          dis *= facoh;

          /* add dissipation term to RHS */
          rl[ccc] += dis;
        }
      } /* end loop over fields */
    } /* end plane loop */
  } /* end dir-loop*/

  /* release mem */
  free(uc);
}


/* Add dissipation terms to vlr. vlr can be varlist for RHS.
   In the interior we use the same Kreiss-Oliger nth order derivative
   operator as in dissipation_add_KO_order. But the order is dropped
   successively to 2 when we approach the boundary.
   So it can do dissipation everywhere, except on the boundary.
   E.g. if order=6 and these are the grid points (marked by o):
o       o       o       o       o       o       o       o  . . .
^       ^       ^       ^_______^_______^_______^_______^__
|       |       |       6th order diss. of size O(h^5) in all of interior
|       |       4th order diss. of size O(h^3)
|       2nd order diss. of size O(h)
no diss.
   In:
     node
     vlu contains evolved fields
     dissfac is dissipation factor
     order is the order of the derivative operator we want (r=order/2) in
           the interior
     cf    describes by which factor we change diss.fac. near boundary
           E.g. if order=8: cf[] = { 0, .7, .8, .9 } changes it by 0.7 one
           point away, by 0.8 two points away, and 0.9 three points away from
           boundary. At the boundary it is always unchanged, so cf[0] is
           always ignored, but cf needs to have order/2 entries!
   Out:
     vlr is the varlist to which we add dissipation terms */
void dissipation_add_taperedKO_order_cf(tNode *node, tVarList *vlr,
                                        tVarList *vlu, double dissfac,
                                        int order, double *cf)
{
  int *n = node->n;
  double *bb = node->bbox;
  int maxn = max3(n[0],n[1],n[2]);
  double *uc = dtensor(maxn);
  int dir;
  double *sw[] = {sw2, sw4, sw6, sw8, sw10, sw12};  /* stencil weights */
  int isw;                /* index into sw of weights we want */
  int srad = order/2;     /* stencil radius */
  int sr;
  double sgn_bou[srad+1]; /* signs near boundary, and last in interior */
  double fac_bou[srad+1]; /* factors near boundary, and last in interior */
  double facoh_bou[srad]; /* factors/h near boundary */

  if(order%2 || order<2 || order>12)
    errorexit("order must be 2,4,6,8,10,12");

  /* set signs/2^order near boundary, and in interior (last entry) */
  /* overall sign = (-1)^(1+order/2), we also devide by 2^ord */
  for(sr=0; sr<srad+1; sr++) sgn_bou[sr] = (-1. + 2*(sr%2))/(1 << (2*sr));

  /* set fac_bou, i.e. fac. near boundary, and in interior (last entry) */
  for(sr=0; sr<srad; sr++) fac_bou[sr] = sgn_bou[sr] * dissfac * cf[sr];
  fac_bou[srad] = sgn_bou[srad] * dissfac;

  /* add dissipation in each direction to RHS */
  for(dir=0; dir<3; dir++)
  {
    int ndir = n[dir];
    double ooh = (ndir-1)/(bb[2*dir+1] - bb[2*dir]);// 1/dist betw. points
    int ord;      /* order we actually use */
    double facoh; /* (-1)^(1+ord/2)/2^ord * dissfac/h */
    int i,j,k;
//if(dir==1 || dir==2) continue;

    /* do nothing if we have too few grid points */
    if(ndir<3) continue;

    /* reduce ord if we have less than order+1 grid points */
    if(ndir<=order) ord = ((ndir-1)/2) * 2;
    else            ord = order;

    /* reset stencil radius, and stencil weight index */
    srad = ord/2;
    isw = srad - 1;

    /* set interior fac. */
    facoh = fac_bou[srad] * ooh; /* (-1)^(1+ord/2)/2^ord * dissfac/h */

    /* set facoh_bou, i.e. fac. near boundary */
    for(sr=0; sr<srad; sr++) facoh_bou[sr] = fac_bou[sr] * ooh;

    /* loop over plane */
    forplaneN(dir, i,j,k, n, 0)
    {
      int i1 = i1_norm(i,j,k, dir); /* 1st and 2nd index in plane */
      int i2 = i2_norm(i,j,k, dir);
      int i0;                       /* index orthogonal to plane */
      int ic,jc,kc, ccc;
      int l;                        /* field index */
      int ib;                       /* index used near boundary */

      /* loop over fields */
      forvl(vlu, l)
      {
        double *ul = Vard(node, Vind(vlu, l)); /* field data pointer */
        double *rl = Vard(node, Vind(vlr, l)); /* RHS data pointer */

        /* fill field arrays uc, i0 runs orth. to plane */
        for(i0=0; i0<ndir; i0++)
        {
          /* set points and their index */
          ijk_inplaneN(dir, ic,jc,kc, i1,i2, i0);
          ccc = Ind_n(ic,jc,kc, n);
          /* set uc */
          uc[i0] = ul[ccc];
        }

        /* loop over inner points in dir */
        for(i0=srad; i0<ndir-srad; i0++)
        {
          double dis;
          int s;

          ijk_inplaneN(dir, ic,jc,kc, i1,i2, i0);
          ccc = Ind_n(ic,jc,kc, n);

          /* Set dissipation term. Here we assume the stencil is symmetric
             (i.e. same weight for uc[i0-s] and uc[i0+s]) and thus use only
             the top half of sw[isw]. */
          dis = sw[isw][srad]*uc[i0];
          for(s=1; s<=srad; s++) dis += sw[isw][s+srad]*(uc[i0-s] + uc[i0+s]);
          dis *= facoh;

          /* add dissipation term to RHS */
          rl[ccc] += dis;
          //printf("ndir=%d i0=%d  dis=%g\n", ndir, i0, dis);
        }
        //if(ndir==5) exit(9);

        /* loop over points near boundary in dir */
        for(ib=1; ib<srad; ib++)
        {
          int is = ib-1;   /* current stencil index */
          double dis;
          int cnt, s;

          /* set current stencil radius */
          sr = ib;

          /* left boundary is at i0 = ib */
          /* right boundary is at i0 = ndir-1 - ib */
          /* => loop over left and right: */
          for(cnt=0, i0=ib; cnt<2; cnt++, i0+=ndir-1 - 2*ib)
          {
            ijk_inplaneN(dir, ic,jc,kc, i1,i2, i0);
            ccc = Ind_n(ic,jc,kc, n);

            /* Set dissipation term. Here we assume the stencil is symmetric
               (i.e. same weight for uc[i0-s] and uc[i0+s]) and thus use only
               the top half of sw[isw]. */
            dis = sw[is][sr]*uc[i0];
            for(s=1; s<=sr; s++) dis += sw[is][s+sr]*(uc[i0-s] + uc[i0+s]);
            dis *= facoh_bou[sr];

            /* add dissipation term to RHS */
            rl[ccc] += dis;
            //printf("ndir=%d i0=%d  dis=%g\n", ndir, i0, dis);
          }
          //if(ndir==5) exit(9);
        }

      } /* end loop over fields */
    } /* end plane loop */
  } /* end dir-loop*/

  /* release mem */
  free(uc);
}

/* calc diss terms for taperedKO from uc and write them into rc */
void diss_taperedKO(int srad, double *sw[], int ndir, const double *uc,
                    double facoh, const double facoh_bou[srad], double *rc)
{
  int isw = srad - 1;
  int i0, ib;

  /* we do nothing at the boundary below */
  rc[0] = rc[ndir-1] = 0;

  /* loop over inner points */
  for(i0=srad; i0<ndir-srad; i0++)
  {
    double dis;
    int s;

    /* Set dissipation term. Here we assume the stencil is symmetric
       (i.e. same weight for uc[i0-s] and uc[i0+s]) and thus use only
       the top half of sw[isw]. */
    dis = sw[isw][srad]*uc[i0];
    for(s=1; s<=srad; s++) dis += sw[isw][s+srad]*(uc[i0-s] + uc[i0+s]);
    dis *= facoh;

    /* save dissipation term for RHS */
    rc[i0] = dis;
    //rc[i0] += dis;
    //printf("ndir=%d i0=%d  dis=%g\n", ndir, i0, dis);
  }
  //if(ndir==5) exit(9);

  /* loop over points near boundary */
  for(ib=1; ib<srad; ib++)
  {
    int is = ib-1;   /* current stencil index */
    double dis;
    int sr, cnt, s;

    /* set current stencil radius */
    sr = ib;

    /* left boundary is at i0 = ib */
    /* right boundary is at i0 = ndir-1 - ib */
    /* => loop over left and right: */
    for(cnt=0, i0=ib; cnt<2; cnt++, i0+=ndir-1 - 2*ib)
    {
      /* Set dissipation term. Here we assume the stencil is symmetric
         (i.e. same weight for uc[i0-s] and uc[i0+s]) and thus use only
         the top half of sw[is]. */
      dis = sw[is][sr]*uc[i0];
      for(s=1; s<=sr; s++) dis += sw[is][s+sr]*(uc[i0-s] + uc[i0+s]);
      dis *= facoh_bou[sr];

      /* save dissipation term for RHS */
      rc[i0] = dis;
      //rc[i0] += dis;
      //printf("ndir=%d i0=%d  dis=%g\n", ndir, i0, dis);
    }
    //if(ndir==5) exit(9);
  }
}
/* new tapered KO diss */
void dissipation_add_taperedKO_order_cf__new(tNode *node, tVarList *vlr,
                                        tVarList *vlu, double dissfac,
                                        int order, double *cf)
{
  int *n = node->n;
  double *bb = node->bbox;
  int maxn = max3(n[0],n[1],n[2]);
  double *uc = dtensor(maxn);
  double *rc = dtensor(maxn);
  int dir;
  double *sw[] = {sw2, sw4, sw6, sw8, sw10, sw12};  /* stencil weights */
  int srad = order/2;     /* stencil radius */
  int sr;
  double sgn_bou[srad+1]; /* signs near boundary, and last in interior */
  double fac_bou[srad+1]; /* factors near boundary, and last in interior */
  double facoh_bou[srad]; /* factors/h near boundary */

  if(order%2 || order<2 || order>12)
    errorexit("order must be 2,4,6,8,10,12");

  /* set signs/2^order near boundary, and in interior (last entry) */
  /* overall sign = (-1)^(1+order/2), we also devide by 2^ord */
  for(sr=0; sr<srad+1; sr++) sgn_bou[sr] = (-1. + 2*(sr%2))/(1 << (2*sr));

  /* set fac_bou, i.e. fac. near boundary, and in interior (last entry) */
  for(sr=0; sr<srad; sr++) fac_bou[sr] = sgn_bou[sr] * dissfac * cf[sr];
  fac_bou[srad] = sgn_bou[srad] * dissfac;

  /* add dissipation in each direction to RHS */
  for(dir=0; dir<3; dir++)
  {
    int ndir = n[dir];
    double ooh = (ndir-1)/(bb[2*dir+1] - bb[2*dir]);// 1/dist betw. points
    int ord;      /* order we actually use */
    double facoh; /* (-1)^(1+ord/2)/2^ord * dissfac/h */
    int i,j,k;
//if(dir==1 || dir==2) continue;

    /* do nothing if we have too few grid points */
    if(ndir<3) continue;

    /* reduce ord if we have less than order+1 grid points */
    if(ndir<=order) ord = ((ndir-1)/2) * 2;
    else            ord = order;

    /* reset stencil radius, and stencil weight index */
    srad = ord/2;

    /* set interior fac. */
    facoh = fac_bou[srad] * ooh; /* (-1)^(1+ord/2)/2^ord * dissfac/h */

    /* set facoh_bou, i.e. fac. near boundary */
    for(sr=0; sr<srad; sr++) facoh_bou[sr] = fac_bou[sr] * ooh;

    /* loop over plane */
    forplaneN(dir, i,j,k, n, 0)
    {
      int i1 = i1_norm(i,j,k, dir); /* 1st and 2nd index in plane */
      int i2 = i2_norm(i,j,k, dir);
      int i0;                       /* index orthogonal to plane */
      int ic,jc,kc, ccc;
      int l;                        /* field index */

      /* loop over fields */
      forvl(vlu, l)
      {
        double *ul = Vard(node, Vind(vlu, l)); /* field data pointer */
        double *rl = Vard(node, Vind(vlr, l)); /* RHS data pointer */

        /* fill field arrays uc, i0 runs orth. to plane */
        for(i0=0; i0<ndir; i0++)
        {
          /* set points and their index */
          ijk_inplaneN(dir, ic,jc,kc, i1,i2, i0);
          ccc = Ind_n(ic,jc,kc, n);
          /* set uc */
          uc[i0] = ul[ccc];
        }

        /* use uc to calc dissipation terms rc */
        diss_taperedKO(srad, sw, ndir, uc, facoh, facoh_bou, rc);

        /* add dissipation terms to RHS */
        for(i0=0; i0<ndir; i0++)
        {
          /* set points and their index */
          ijk_inplaneN(dir, ic,jc,kc, i1,i2, i0);
          ccc = Ind_n(ic,jc,kc, n);
          /* add rc to RHS */
          rl[ccc] += rc[i0];
        }
      } /* end loop over fields */
    } /* end plane loop */
  } /* end dir-loop*/




  /* release mem */
  free(rc);
  free(uc);
}



/* same as dissipation_add_taperedKO_order_cf, but use same diss. factor for
   all diss. orders near boundary and in interior */
void dissipation_add_taperedKO_order(tNode *node, tVarList *vlr, tVarList *vlu,
                                     double dissfac, int order)
{
  double cf[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1};
  dissipation_add_taperedKO_order_cf(node, vlr, vlu, dissfac, order, cf);
}
void dissipation_add_taperedKO_order__new(tNode *node, tVarList *vlr, tVarList *vlu,
                                     double dissfac, int order)
{
  double cf[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1};
  dissipation_add_taperedKO_order_cf__new(node, vlr, vlu, dissfac, order, cf);
}

/* use dissipation_add_taperedKO_order_cf with some cf zeroed near boundary,
   so that no diss below min_order is used */
void dissipation_add_taperedKO_order_min(tNode *node,
                                         tVarList *vlr, tVarList *vlu,
                                         double dissfac, int order,
                                         int min_order)
{
  double cf[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1};
  int i;
  for(i=0; i<min_order/2; i++) cf[i]=0.;
  dissipation_add_taperedKO_order_cf(node, vlr, vlu, dissfac, order, cf);
}
