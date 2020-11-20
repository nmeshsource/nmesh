/* WENO.c */
/* Wolfgang Tichy, Nov. 2020 */

/* WENO for non-uniform grid in 1D */

#include "nmesh.h"
#include "dg.h"



/* Here we use Shu-WENO-notes.pdf Eq. 2.20:

  Eq. 2.10:
  v^{-}_{i + 1/2} := \sum_{j=0}^{k-1} c^{(k)}_{rj} \bar{v}_{i-r+j}
  v^{+}_{i - 1/2} := \sum_{j=0}^{k-1} \tilde{c}^{(k)}_{rj} \bar{v}_{i-r+j}
  Note text after 2.10: \tilde{c}^{(k)}_{rj} = c^{(k)}_{r-1,j}

  Eq. 2.11:
  v_{i + 1/2} := \sum_{j=0}^{k-1} c^{(k)}_{rj} \bar{v}_{i-r+j}


  This is how we get the c^{(k)}_{rj} :

  n := k+1

  w_m := 1/( \prod_{l=0, l\neq m}^{n-1} (xp_m - xp_l) )

  w_m = w[m] can be obtained from
  void Lagrange_winterp(int n, const double *x, double *w_interp) :
     xp[l] = x_{i - r - 1/2 + l}
     Lagrange_winterp(k+1, xp, w);


  L_{ml}(x) := \prod_{q=0, q\neq m,l}^{n-1} (x - xp_q)

  L_{ml}(x) can be obtained from
    double Lagrange_prod2(int l, int m, double x, int np, const double *x_p) :
    L_ml_x = Lagrange_prod2(l,m, x, k+1, xp);

  c^{(k)}_{rj} := \Delta x_{i-r+j} *
                  \sum_{m=j+1}^k w_m \sum_{l=0, l\neq m}^k L_{ml}(x)

  where
    x = x_{i + 1/2}
    \Delta x_{i-r+j} = x_{i - r + j + 1/2} - x_{i - r + j - 1/2}


  For WENO we need:
  -----------------
  v_{i + 1/2} = \sum_{r=0}^{k-1} d_r v_{i + 1/2}^{(r)}
  where from Eq. 2.51:
    v_{i + 1/2}^{(r)} = \sum_{j=0}^{k-1} c^{(k)}_{rj} \bar{v}_{i-r+j}

  The d_r are picked s.t.
  v_{i + 1/2} = v(x_{i + 1/2}) + O(\Delta x)^{2k-1}   <-- Eq. 2.54

  This is achieved by noting that
  v_{i + 1/2}
   = \sum_{j=0}^{k} c^{(k+1)}_{r'j} \bar{v}_{i-r'+j} + O(\Delta x)^{2k-1}
   where we choose r' later

  Thus we need to demand that:
  \sum_{r=0}^{k-1} d_r v_{i + 1/2}^{(r)}
   = \sum_{j=0}^{k} c^{(k+1)}_{r'j} \bar{v}_{i-r'+j}

  From this we can get the d_r. E.g. for k=2:
  d_0 ( c^{(2)}_{00} \bar{v}_{i}   + c^{(2)}_{01} \bar{v}_{i+1} )
  d_1 ( c^{(2)}_{10} \bar{v}_{i-1} + c^{(2)}_{11} \bar{v}_{i} )
    = c^{(3)}_{10} \bar{v}_{i-1+0} + c^{(3)}_{11} \bar{v}_{i-1+1} +
      c^{(3)}_{12} \bar{v}_{i-1+2}
  [ here we set r'=1 so that the RHS has
    \bar{v}_{i-1}, \bar{v}_{i}, \bar{v}_{i+1} like thr LHS ]

  Looking at the terms \bar{v}_{i+1} and \bar{v}_{i-1} we see that we need:
  d_0 = c^{(3)}_{12} / c^{(2)}_{01}
  d_1 = c^{(3)}_{10} / c^{(2)}_{10}

*/


/* Compute the coeffs c^{(k)}_{rj} = c _{rj} as in Eq. 2.20
   of Shu-WENO-notes.pdf
   In: np, pt,  i,k,r,j
   the array pt contains all the np midpoints
   pt[...]  is [ x_{-1/2}, x_{1/2}, ..., x_{-1/2 + np-1} ] */
double Shu_WENO_c_k_rj(int np, const double *pt, int i, int k, int r, int j)
{
  double x;
  int m, l;
  const double *xp;
  double *w = dmalloc(k+1);
  double Deltax, sum, s2, coeff;
  int pr=0;

  /* set xp */
  if( i-r < 0 ) errorexit("i-r out of range");
  xp = pt + (i-r);
  if(i-r + k > np-1) errorexit("i-r + k out of range");

  /* get \Delta x_{i-r+j} */
  if(j < 0 || j>=k) errorexit("j out of range");
  Deltax = xp[j+1] - xp[j];
  if(pr) printf("np=%d: i=%d k=%d r=%d j=%d: Deltax=%g\n", np, i,k,r,j, Deltax);

  /* get weights w */
  Lagrange_winterp(k+1, xp, w);
  if(pr) for(m=0; m<=k; m++) printf("  w[m]=%g\n", w[m]);

  /* set x = x_{i+1/2} */
  if( i+1 < 0 || i+1 > np-1 ) errorexit("i+1 out of range");
  x = pt[i+1];
  if(pr) printf("  x=%g\n", x);

  /* sum to get c^{(k)}_{rj} */
  sum = 0.;
  for(m=j+1; m<=k; m++)
  {
    s2 = 0.;
    for(l=0; l<=k; l++)
    {
      if(l!=m)
      {
        double L_ml_x = Lagrange_prod2(l,m, x, k+1, xp);
        if(pr) printf("    L_%d%d_x=%g\n", m,l, L_ml_x);
        s2 += L_ml_x;
      }
    }
    sum += w[m] * s2;
  }
  /* set coeff c^{(k)}_{rj} */
  coeff = sum * Deltax;

  free(w);
  return coeff;
}

/* weights to interpolate to midpoint at i+1/2 with index i
   Here we interpolate in the positive direction (p) from the left of the
   midpoint to the midpoint to obtain v_{i+1/2}.
   We use v at the points i-1, i, i+1 */
void Shu_p_WENO3_weights(int np, const double *pt, int i, double d[2])
{
  /* we use:
  v_{i+1/2} =
    d_0 ( c^{(2)}_{00} \bar{v}_{i}   + c^{(2)}_{01} \bar{v}_{i+1} )
    d_1 ( c^{(2)}_{10} \bar{v}_{i-1} + c^{(2)}_{11} \bar{v}_{i} )
      = c^{(3)}_{10} \bar{v}_{i-1+0} + c^{(3)}_{11} \bar{v}_{i-1+1} +
        c^{(3)}_{12} \bar{v}_{i-1+2} */
  double c2_01 = Shu_WENO_c_k_rj(np,pt, i, 2, 0,1);
  double c2_10 = Shu_WENO_c_k_rj(np,pt, i, 2, 1,0);
  double c3_12 = Shu_WENO_c_k_rj(np,pt, i, 3, 1,2);
  double c3_10 = Shu_WENO_c_k_rj(np,pt, i, 3, 1,0);

  d[0] = c3_12/c2_01;
  d[1] = c3_10/c2_10;
}

/* weights to interpolate to midpoint at i+1/2 with index i
   Here we interpolate in the negative direction (m) from the right of the
   midpoint to the midpoint to obtain v_{i+1/2}.
   We use v at the points i, i+1, i+2 */
void Shu_m_WENO3_weights(int np, const double *pt, int i, double d[2])
{
  /* we use:
  v_{i+1/2} =
    d_0 ( c^{(2)}_{00} \bar{v}_{i}   + c^{(2)}_{01} \bar{v}_{i+1} )
    d_1 ( c^{(2)}_{-10} \bar{v}_{i+1} + c^{(2)}_{-11} \bar{v}_{i+2} )
      = c^{(3)}_{00} \bar{v}_{i+0} + c^{(3)}_{01} \bar{v}_{i+1} +
        c^{(3)}_{02} \bar{v}_{i+2} */
  double c2_00  = Shu_WENO_c_k_rj(np,pt, i, 2, 0,0);
  double c2_m11 = Shu_WENO_c_k_rj(np,pt, i, 2, -1,1);
  double c3_00 = Shu_WENO_c_k_rj(np,pt, i, 3, 0,0);
  double c3_02 = Shu_WENO_c_k_rj(np,pt, i, 3, 0,2);

  d[0] = c3_00/c2_00;
  d[1] = c3_02/c2_m11;
}






/* print the c_k_rj and the resulting WENO3 weights */
int pr_c_k_rj(tMesh *mesh)
{
  double pt[] = {0., 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.8 };
  int np = sizeof(pt)/sizeof(double);
  double c, d[2];
  int i, k, r,j;

  PRF;printf(": np=%d cell boundary points\n", np);

  i=4;
  printf("i=%d:\n", i);

  k=2;
  /* Shu_WENO_c_k_rj(int np, const double *pt, int i, int k, int r, int j) */
  for(r=-1; r<=k-1; r++)
    for(j=0; j<k; j++)
    {
      c = Shu_WENO_c_k_rj(np,pt, i, k, r, j);
      printf("  c%d_%d%d = %g\n", k, r,j, c);
    }

  printf("i=%d:\n", i);

  k=3;
  /* Shu_WENO_c_k_rj(int np, const double *pt, int i, int k, int r, int j) */
  for(r=-1; r<=k-1; r++)
    for(j=0; j<k; j++)
    {
      c = Shu_WENO_c_k_rj(np,pt, i, k, r, j);
      printf("  6*c%d_%d%d = %g\n", k, r,j, 6.*c);
    }


  i=0;
  printf("i=%d:\n", i);
  Shu_m_WENO3_weights(np,pt, i, d);
  printf("m: d0=%g d1=%g   d0/d1=%g\n", d[0], d[1], d[0]/d[1]);

  i=1;
  printf("i=%d:\n", i);
  Shu_m_WENO3_weights(np,pt, i, d);
  printf("m: d0=%g d1=%g   d0/d1=%g\n", d[0], d[1], d[0]/d[1]);
  Shu_p_WENO3_weights(np,pt, i, d);
  printf("p: d0=%g d1=%g   d0/d1=%g\n", d[0], d[1], d[0]/d[1]);

  i=2;
  printf("i=%d:\n", i);
  Shu_m_WENO3_weights(np,pt, i, d);
  printf("m: d0=%g d1=%g   d0/d1=%g\n", d[0], d[1], d[0]/d[1]);
  Shu_p_WENO3_weights(np,pt, i, d);
  printf("p: d0=%g d1=%g   d0/d1=%g\n", d[0], d[1], d[0]/d[1]);

  i=np-4;
  printf("i=%d:\n", i);
  Shu_p_WENO3_weights(np,pt, i, d);
  printf("p: d0=%g d1=%g   d0/d1=%g\n", d[0], d[1], d[0]/d[1]);
  Shu_m_WENO3_weights(np,pt, i, d);
  printf("m: d0=%g d1=%g   d0/d1=%g\n", d[0], d[1], d[0]/d[1]);

  i=np-3;
  printf("i=%d:\n", i);
  Shu_p_WENO3_weights(np,pt, i, d);
  printf("p: d0=%g d1=%g   d0/d1=%g\n", d[0], d[1], d[0]/d[1]);

exit(99);
  return 0;
}
