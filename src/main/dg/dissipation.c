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
   So our dissfac below is given by dissfac = \sigma/16.
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
    double dissfacoh = dissfac * ooh; /* dissfac/h */
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
          rl[ccc] += -dissfacoh*( 6.*uc[i0] - 4.*(uc[i0-1] + uc[i0+1])
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
  double sgn = -1. + (order%4); /* overall sign */
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
    double dissfacoh = sgn * dissfac * ooh; /* (-1)^(1+order/2) dissfac/h */
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
          dis *= dissfacoh;

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
   successively to 2 when we approach the boundary
   In:
     node
     vlu contains evolved fields
     dissfac is dissipation factor
     order is the order of the derivative operator we want (r=order/2)
   Out:
     vlr is the varlist to which we add dissipation terms */
void dissipation_add_taperedKO_order(tNode *node, tVarList *vlr, tVarList *vlu,
                                     double dissfac, int order)
{
  int *n = node->n;
  double *bb = node->bbox;
  int maxn = max3(n[0],n[1],n[2]);
  double *uc = dtensor(maxn);
  int dir;
  double *sw[] = {sw2, sw4, sw6, sw8, sw10, sw12};  /* stencil weights */
  int isw;                      /* index into sw of weights we want */
  double sgn = -1. + (order%4); /* overall sign */
  int srad = order/2;           /* stencil radius */

  /* choose stencil weight index */
  isw = srad - 1;
  if(isw<0 || isw> 5) errorexit("order must be 2,4,6,8,10,12");

  /* add dissipation in each direction to RHS */
  for(dir=0; dir<3; dir++)
  {
    double ooh = (n[dir]-1)/(bb[2*dir+1] - bb[2*dir]);// 1/dist betw. points
    double dissfacoh = dissfac * ooh; /* dissfac/h */
    double sdissfacoh = sgn * dissfacoh; /* (-1)^(1+order/2) dissfac/h */
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
      int ib;                       /* index used near boundary */

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
             the top half of sw[isw]. */
          dis = sw[isw][srad]*uc[i0];
          for(s=1; s<=srad; s++) dis += sw[isw][s+srad]*(uc[i0-s] + uc[i0+s]);
          dis *= sdissfacoh;

          /* add dissipation term to RHS */
          rl[ccc] += dis;
        }

        /* loop over points near boundary in dir */
        for(ib=1; ib<srad; ib++)
        {
          int is = ib-1;   /* current stencil index */
          int sr = ib;     /* current stencil radius */
          int sign = -1. + 2*(sr%2); /* overall sign */
          double dis;
          int s;

          /* left boundary is at i0 = ib */
          /* right boundary is at i0 = n[dir]-1 - ib */
          /* => loop over left and right: */
          for(i0 = ib; i0 < n[dir]-1; i0 += n[dir]-1 - 2*ib)
          {
            ijk_inplaneN(dir, ic,jc,kc, i1,i2, i0);
            ccc = Ind_n(ic,jc,kc, n);

            /* Set dissipation term. Here we assume the stencil is symmetric
               (i.e. same weight for uc[i0-s] and uc[i0+s]) and thus use only
               the top half of sw[isw]. */
            dis = sw[is][sr]*uc[i0];
            for(s=1; s<=sr; s++) dis += sw[is][s+sr]*(uc[i0-s] + uc[i0+s]);
            dis *= sign * dissfacoh; // FIXME
            errorexit("fix wrong sign in dissfacoh");

            /* add dissipation term to RHS */
            rl[ccc] += dis;
          }
        }

      } /* end loop over fields */
    } /* end plane loop */
  } /* end dir-loop*/

  /* release mem */
  free(uc);
}
