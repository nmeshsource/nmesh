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
