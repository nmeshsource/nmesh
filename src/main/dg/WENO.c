/* WENO.c */
/* Wolfgang Tichy, Nov. 2020 */

/* WENO for non-uniform grid in 1D */

#include "nmesh.h"
#include "dg.h"


/* global vars to keep some weights */
tWENOweights *WENOweights;



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
   We use v at the points i-1, i, i+1
   Note: the 1st point pt[0] is a face point and thus at i=-1 */
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
   We use v at the points i, i+1, i+2
   Note: the last point pt[np-1] is a face point and thus at i=np-2 */
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
void pr_Shu_c_k_rj_AND_d(int np, double pt[])
{
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
}


/* WENO comes form:
   XD. Liu, S. Osher, T. Chan
   "Weighted Essentially Non-Oscillatory Schemes",
   J. of Comput. Physics 115, p. 200 (1994)
   [ LiuOsherChan_WENO_JComputPhys115p200_1994.pdf ]

  In Sec. 3.3 they introduce stencils S_{j+k} around point x_{j+1/2}
  e.g. for r=2: k=0,1 (since k=0, ..., r-1) we get:
    S_{j}   = (x_{j-3/2}, x_{j-1/2}, x_{j+1/2})       <-- k=0
    S_{j+1} = (x_{j-1/2}, x_{j+1/2}, x_{j+3/2})       <-- k=1

  In Eq. 3.6 the interpolation polynomial is
  R_j(x) = \sum_{k=0}^{r-1} \alpha^j_k p'_{j+k}(x) / (\sum_l \alpha^j_l)

  in Eq. 3.9 they set
  \alpha^j_k = C^j_k / (\epsilon + IS_{j+k})^r

  Their WENO idea is to choose the C^j_k such that the accuracy of
  R_j(x) is improved in smoothe regions.

  After Eq. 3.11a they introduce
  a^j_k(x) = \sum_{s=0}^r \prod_{l=0, l\neq s}^r (x - x_{j + k - l + 1/2})

  The truncation error of R_j(x) is given on Eq. 3.11c. It is samller if
  \sum_{k=0}^{r-1} C^j_k a^j_k = 0
  THIS is the main idea choose C^j_k such that this is zero!

  If r=2
  a^j_k(x) = (x - x_{j+k-1/2})(x - x_{j+k-3/2}) +
             (x - x_{j+k+1/2})(x - x_{j+k-3/2}) +
             (x - x_{j+k+1/2})(x - x_{j+k-1/2})

  a^j_0(x) = (x - x_{j-1/2})(x - x_{j-3/2}) +
             (x - x_{j+1/2})(x - x_{j-3/2}) +
             (x - x_{j+1/2})(x - x_{j-1/2})

  a^j_1(x) = (x - x_{j+1/2})(x - x_{j-1/2}) +
             (x - x_{j+3/2})(x - x_{j-1/2}) +
             (x - x_{j+3/2})(x - x_{j+1/2})

  So at x = x_{j+1/2}
  a^j_0(x_{j+1/2}) = (x_{j+1/2} - x_{j-1/2})(x_{j+1/2} - x_{j-3/2}) +
                     (x_{j+1/2} - x_{j+1/2})(x_{j+1/2} - x_{j-3/2}) +
                     (x_{j+1/2} - x_{j+1/2})(x_{j+1/2} - x_{j-1/2})
                   = (x_{j+1/2} - x_{j-1/2})(x_{j+1/2} - x_{j-3/2})

  a^j_1(x_{j+1/2}) = (x_{j+1/2} - x_{j+1/2})(x_{j+1/2} - x_{j-1/2}) +
                     (x_{j+1/2} - x_{j+3/2})(x_{j+1/2} - x_{j-1/2}) +
                     (x_{j+1/2} - x_{j+3/2})(x_{j+1/2} - x_{j+1/2})
                   = (x_{j+1/2} - x_{j+3/2})(x_{j+1/2} - x_{j-1/2})

  So C^j_0 a^j_0(x_{j+1/2}) + C^j_1 a^j_1(x_{j+1/2}) = 0 gives
  C^j_1 = -[a^j_0(x_{j+1/2}) / a^j_1(x_{j+1/2})] C^j_0
        = -[ (x_{j+1/2} - x_{j-3/2}) / (x_{j+1/2} - x_{j+3/2}) ] C^j_0
        = +[ (x_{j+1/2} - x_{j-3/2}) / (x_{j+3/2} - x_{j+1/2}) ] C^j_0

  These C^j_1 , C^j_0 are the ideal weights [ C^j_0 + C^j_1 = 1 ]
  for interpolation from the left of x_{j+1/2} towards it (p).


  Do the same at x = x_{j-1/2}
  a^j_0(x_{j+1/2}) = (x_{j-1/2} - x_{j-1/2})(x_{j-1/2} - x_{j-3/2}) +
                     (x_{j-1/2} - x_{j+1/2})(x_{j-1/2} - x_{j-3/2}) +
                     (x_{j-1/2} - x_{j+1/2})(x_{j-1/2} - x_{j-1/2})
                   = (x_{j-1/2} - x_{j+1/2})(x_{j-1/2} - x_{j-3/2})

  a^j_1(x_{j+1/2}) = (x_{j-1/2} - x_{j+1/2})(x_{j-1/2} - x_{j-1/2}) +
                     (x_{j-1/2} - x_{j+3/2})(x_{j-1/2} - x_{j-1/2}) +
                     (x_{j-1/2} - x_{j+3/2})(x_{j-1/2} - x_{j+1/2})
                   = (x_{j-1/2} - x_{j+3/2})(x_{j-1/2} - x_{j+1/2})

  So C^j_0 a^j_0(x_{j-1/2}) + C^j_1 a^j_1(x_{j-1/2}) = 0 gives
  C^j_1 = -[a^j_0(x_{j-1/2}) / a^j_1(x_{j-1/2})] C^j_0
        = -[ (x_{j-1/2} - x_{j-3/2}) / (x_{j-1/2} - x_{j+3/2}) ] C^j_0
        = +[ (x_{j-1/2} - x_{j-3/2}) / (x_{j+3/2} - x_{j-1/2}) ] C^j_0

  These C^j_1 , C^j_0 are the ideal weights [ C^j_0 + C^j_1 = 1 ]
  for interpolation from the right of x_{j-1/2} towards it (m).   */


/* weights to interpolate to midpoint at i+1/2 with index i
   Here we interpolate in the positive direction (p) from the left of the
   midpoint to the midpoint to obtain v_{i+1/2}.
   We use v at the points i-1, i, i+1 */
double LiuOsherChan_p_WENO3_weight_ratio(int np, const double *pt, int i)
{
  /* we use:
     C^j_1 = +[ (x_{j+1/2} - x_{j-3/2}) / (x_{j+3/2} - x_{j+1/2}) ] C^j_0 */
  int j = i;
  double numer = pt[j+1] - pt[j-1];
  double denom = pt[j+2] - pt[j+1];

  return numer/denom; /* C^j_1 /  C^j_0 */
}

/* weights to interpolate to midpoint at i+1/2 with index i
   Here we interpolate in the negative direction (m) from the right of the
   midpoint to the midpoint to obtain v_{i+1/2}.
   We use v at the points i, i+1, i+2 */
double LiuOsherChan_m_WENO3_weight_ratio(int np, const double *pt, int i)
{
  /* this:
     C^j_1 = +[ (x_{j-1/2} - x_{j-3/2}) / (x_{j+3/2} - x_{j-1/2}) ] C^j_0
     interpolates to x_{j-1/2}. We need to shift this over by 1 to the right
     to get the result at x_{j+1/2} */
  int j = i+1;
  double denom = pt[j]   - pt[j-1];
  double numer = pt[j+2] - pt[j];

  return numer/denom; /* C^j_0 /  C^j_1 */
}



/* Set WENO3 weights for interpolation from left to right (p).
   We pass in the np face- and mid-points.
   E.g.:
       | x |   o   |   o   |   o   |   o   |   o   | x |
   pt[0]  pt[1]   pt[2]   ...                          pt[np-1]

   This func computes the np-2 sets of weights at the np-2 mid-points
   i=0      i=1      ... i=(np-2)-1
   pt[1]    pt[2]    ... pt[np-2]    */
void set_mid_p_WENO3weights(int np, const double *pt, tWENO3weight **W3)
{
  int pr=0;
  int i;

  if(pr)
  {
    printf("\n");
    PRFs(":\n");
    printf("midpoints (including face points):\n");
    for(i=-1; i<(np)-1; i++) printf("pt[%d]=%g  ", i, pt[i+1]);
    printf("\n");
  }

  for(i=1; i<(np-2); i++)
  {
    int k, r, j;
    double ratio;

    k=2;
    //printf("At midpoint with index i=%d:\n", i);
    for(r=1; r>=0; r--)
    {
      for(j=0; j<k; j++)
      {
        double c = Shu_WENO_c_k_rj(np,pt, i, k, r, j);
        W3[i]->lw[1-r][j] = c;
      }
    }
    ratio = LiuOsherChan_p_WENO3_weight_ratio(np,pt, i);
    W3[i]->optw[0] = 1.;
    W3[i]->optw[1] = ratio;
  }
}
/* Set WENO3 weights for interpolation from right to left (m). */
void set_mid_m_WENO3weights(int np, const double *pt, tWENO3weight **W3)
{
  int pr=0;
  int i;

  if(pr)
  {
    printf("\n");
    PRFs(":\n");
    printf("midpoints (including face points):\n");
    for(i=-1; i<(np)-1; i++) printf("pt[%d]=%g  ", i, pt[i+1]);
    printf("\n");
  }

  for(i=0; i<(np-2)-1; i++)
  {
    int k, r, j;
    double ratio;

    k=2;
    //printf("At midpoint with index i=%d:\n", i);
    for(r=0; r>=-1; r--)
    {
      for(j=0; j<k; j++)
      {
        double c = Shu_WENO_c_k_rj(np,pt, i, k, r, j);
        W3[i]->lw[-r][j] = c;
      }
    }
    ratio = LiuOsherChan_m_WENO3_weight_ratio(np,pt, i);
    W3[i]->optw[1] = 1.;
    W3[i]->optw[0] = ratio;
  }
}
/* print WENO3 weights in W3 */
void pr_WENO3weight(tWENO3weight *W3)
{
  int l;
  for(l=0; l<2; l++)
    printf("lw[%d][]=%g %g  ", l, W3->lw[l][0], W3->lw[l][1]);
  printf(":  optw[]=%g %g\n", W3->optw[0], W3->optw[1]);
}

/* print WENO3 weights in W3 */
void pr_WENO3weights(int nm, tWENO3weight **W3)
{
  int i;
  for(i=0; i<nm; i++)
  {
    printf("%d: ", i);
    pr_WENO3weight(W3[i]);
  }
}


/********************************************************************/
/* funcs to generate tWENOweights object */
/********************************************************************/

/* allocate array with WENO3 weights */
tWENO3weight **WENO3weights_alloc(int nm)
{
  int i;
  tWENO3weight **W3 = calloc(nm, sizeof(W3[0]));
  for(i=0; i<nm; i++)
    W3[i] = calloc(1, sizeof(W3[0][0]));
  return W3;
}
/* free WENO3 weights */
void WENO3weights_free(int nm, tWENO3weight **W3)
{
  int i;
  for(i=0; i<nm; i++) free(W3[i]);
  free(W3);
}

/* allocate the WENO weight object */
tWENOweights *WENOweights_alloc(int n)
{
  tWENOweights *w = calloc(1, sizeof(w[0]));
  w->n = n;
  w->p_WENO3 = WENO3weights_alloc(n);
  w->m_WENO3 = WENO3weights_alloc(n);
  return w;
}
/* free WENO weight object */
void WENOweights_free(tWENOweights *w)
{
  WENO3weights_free(w->n, w->p_WENO3);
  WENO3weights_free(w->n, w->m_WENO3);
  free(w);
}

/* set tWENOweights object */
void WENOweights_set(tWENOweights *w, int np, const double pt[])
{
  set_mid_p_WENO3weights(np,pt, w->p_WENO3);
  set_mid_m_WENO3weights(np,pt, w->m_WENO3);
}

void WENOweights_print(tWENOweights *w)
{
  int n = w->n;
  PRFs(":\n");
  printf(" ->n = %d:\n", n);
  printf("p_WENO3:\n");
  pr_WENO3weights(n, w->p_WENO3);
  printf("m_WENO3:\n");
  pr_WENO3weights(n, w->m_WENO3);
}

/* set global WENOweights var */
void WENOweights_global_init(void)
{
  /*    | x |   o   |   o   |   o   |   o   | x |
    i:      0       1       2       3       4        <-- np-2 diff. i's */
  double pt[] = { 0., 0.5, 1.5, 2.5, 3.5, 4.5, 5. }; // grid we use for FV
  int np = sizeof(pt)/sizeof(double);

  WENOweights = WENOweights_alloc(np-2);
  WENOweights_set(WENOweights, np,pt);
}
void WENOweights_global_free(void)
{
  WENOweights_free(WENOweights);
}
/* same but with nmesh interface */
int WENOweights_init_globals(tMesh *mesh)
{
  WENOweights_global_init();
  return 0;
}
int WENOweights_free_globals(tMesh *mesh)
{
  WENOweights_global_free();
  return 0;
}


/* get access to weights at some points */
tWENO3weight *WENOweights_global_p_WENO3_at_(int i)
{
  tWENO3weight *W3 = WENOweights->p_WENO3[i];
  return W3;
}
tWENO3weight *WENOweights_global_m_WENO3_at_(int i)
{
  tWENO3weight *W3 = WENOweights->m_WENO3[i];
  return W3;
}
tWENO3weight *WENOweights_global_p_WENO3_at_last_minus_(int l)
{
  int n = WENOweights->n;
  tWENO3weight *W3 = WENOweights->p_WENO3[n-1-l];
  return W3;
}
tWENO3weight *WENOweights_global_m_WENO3_at_last_minus_(int l)
{
  int n = WENOweights->n;
  tWENO3weight *W3 = WENOweights->m_WENO3[n-1-l];
  return W3;
}


/*********************************************************/
/* print stuff */
/*********************************************************/

/* consider this grid:
    | x |   o   |   o   |   o   |   o   |   o   | x |
i:      0       1       2       3       4       5       <-- np-2 diff. i's
   the bars are the cell boundaries and their coords are put into the
   pt array below */
void print_p_WENO3_weights(void)
{
  double pt[] = { 0., 0.5, 1.5, 2.5, 3.5, 4.5, 5. }; // grid example
  int np = sizeof(pt)/sizeof(double);
  int i;

  printf("\n");
  PRFs(":\n");
  printf("midpoints (including face points):\n");
  for(i=-1; i<(np)-1; i++) printf("pt[%d]=%g  ", i, pt[i+1]);
  printf("\n");

  for(i=1; i<(np-2); i++)
  {
    int k, r, j;
    double d[2];

    k=2;

    printf("At midpoint with index i=%d:\n", i);
    /* Shu_WENO_c_k_rj(int np, const double *pt, int i, int k, int r, int j) */
    printf("  weights for 2 linear reconstructions:\n");
    for(r=1; r>=0; r--)
    {
      printf("    reconstruction%d: ", 2-r);
      for(j=0; j<k; j++)
      {
        double c = Shu_WENO_c_k_rj(np,pt, i, k, r, j);
        printf("c%d_%d%d=%g ", k, r,j, c);
      }
      printf("\n");
    }

    Shu_p_WENO3_weights(np,pt, i, d);
    printf("  p_WENO3: ideal WENO3 weights: d[0]=%g d[1]=%g\n", d[0], d[1]);
    printf("  p_WENO3: ideal WENO3 weight ratio=%g\n",
           LiuOsherChan_p_WENO3_weight_ratio(np,pt, i));
  }

  tWENO3weight **W3 = WENO3weights_alloc(np-2);
  set_mid_p_WENO3weights(np,pt, W3);
  pr_WENO3weights(np-2, W3);
  WENO3weights_free(np-2, W3);
}
void print_m_WENO3_weights(void)
{
  double pt[] = { 0., 0.5, 1.5, 2.5, 3.5, 4.5, 5. };
  int np = sizeof(pt)/sizeof(double);
  int i;

  printf("\n");
  PRFs(":\n");
  printf("midpoints (including face points):\n");
  for(i=-1; i<(np)-1; i++) printf("pt[%d]=%g  ", i, pt[i+1]);
  printf("\n");

  for(i=0; i<(np-2)-1; i++)
  {
    int k, r, j;
    double d[2];

    k=2;

    printf("At midpoint with index i=%d:\n", i);
    /* Shu_WENO_c_k_rj(int np, const double *pt, int i, int k, int r, int j) */
    printf("  weights for 2 linear reconstructions:\n");
    for(r=0; r>=-1; r--)
    {
      printf("    reconstruction%d: ", 1-r);
      for(j=0; j<k; j++)
      {
        double c = Shu_WENO_c_k_rj(np,pt, i, k, r, j);
        printf("c%d_%d%d=%g ", k, r,j, c);
      }
      printf("\n");
    }

    Shu_m_WENO3_weights(np,pt, i, d);
    printf("  m_WENO3: ideal WENO3 weights: d[0]=%g d[1]=%g\n", d[0], d[1]);
    printf("  m_WENO3: ideal WENO3 weight ratio=%g\n",
           LiuOsherChan_m_WENO3_weight_ratio(np,pt, i));
  }

  tWENO3weight **W3 = WENO3weights_alloc(np-2);
  set_mid_m_WENO3weights(np,pt, W3);
  pr_WENO3weights(np-2, W3);
  WENO3weights_free(np-2, W3);
}



/* print the resulting WENO3 weights */
int pr_weight_ratios(tMesh *mesh)
{
  double pt[] = { -0.05, 0., 0.05, 0.15, 0.25, 0.35, 0.45,
                  0.55, 0.65, 0.75, 0.8 };
  int np = sizeof(pt)/sizeof(double);
  int i;
  tWENO3weight *W3;

  /* first Shu things */
  pr_Shu_c_k_rj_AND_d(np, pt);

  PRF;printf(": np=%d cell boundary points\n", np);

  i=0;
  printf("i=%d:\n", i);
  printf("m: ratio=%g\n", LiuOsherChan_m_WENO3_weight_ratio(np,pt, i));

  i=1;
  printf("i=%d:\n", i);
  printf("m: ratio=%g\n", LiuOsherChan_m_WENO3_weight_ratio(np,pt, i));
  printf("p: ratio=%g\n", LiuOsherChan_p_WENO3_weight_ratio(np,pt, i));

  i=2;
  printf("i=%d:\n", i);
  printf("m: ratio=%g\n", LiuOsherChan_m_WENO3_weight_ratio(np,pt, i));
  printf("p: ratio=%g\n", LiuOsherChan_p_WENO3_weight_ratio(np,pt, i));

  i=np-4;
  printf("i=%d:\n", i);
  printf("p: ratio=%g\n", LiuOsherChan_p_WENO3_weight_ratio(np,pt, i));
  printf("m: ratio=%g\n", LiuOsherChan_m_WENO3_weight_ratio(np,pt, i));

  i=np-3;
  printf("i=%d:\n", i);
  printf("p: ratio=%g\n", LiuOsherChan_p_WENO3_weight_ratio(np,pt, i));

  print_p_WENO3_weights();
  print_m_WENO3_weights();

  WENOweights_global_init();
  WENOweights_print(WENOweights);

  printf("some weights:\n");

  printf("m: 0: ");
  W3 = WENOweights_global_m_WENO3_at_(0);
  pr_WENO3weight(W3);
  printf("m: 1: ");
  W3 = WENOweights_global_m_WENO3_at_(1);
  pr_WENO3weight(W3);

  printf("p: 1: ");
  W3 = WENOweights_global_p_WENO3_at_(1);
  pr_WENO3weight(W3);
  printf("p: 2: ");
  W3 = WENOweights_global_p_WENO3_at_(2);
  pr_WENO3weight(W3);

  printf("p: -1: ");
  W3 = WENOweights_global_p_WENO3_at_last_minus_(1);
  pr_WENO3weight(W3);
  printf("p: -0: ");
  W3 = WENOweights_global_p_WENO3_at_last_minus_(0);
  pr_WENO3weight(W3);

  printf("m: -2: ");
  W3 = WENOweights_global_m_WENO3_at_last_minus_(2);
  pr_WENO3weight(W3);
  printf("m: -1: ");
  W3 = WENOweights_global_m_WENO3_at_last_minus_(1);
  pr_WENO3weight(W3);

  WENOweights_global_free();
  exit(99);
  return 0;
}
