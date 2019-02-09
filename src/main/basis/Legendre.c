/* WT 2/2018 
   some functions to deal with Legendre polynomials, especially at the
   Legendre Gauss-Lobatto (LGL) points.
*/

#include "nmesh.h"
#include "basis.h"



/* associated Legendre polynomial (not really needed here) */
/* use recursion recommended in numrec to get P_l^m for positive m */
double assocLegendreP(int l, int m, double x)
{
  double fact, Pkm, Pmm, Pmp1_m, somx2;
  int i, k;

  if(m < 0 || m > l || fabs(x) > 1.0) errorexit("Bad arguments");

  /* build P_m^m */
  Pmm = 1.0;
  if(m > 0)
  {
    somx2 = sqrt((1.0-x)*(1.0+x));
    fact=1.0;
    for(i=1; i<=m; i++)
    {
      Pmm *= -fact*somx2;
      fact += 2.0;
    }
  }
  /* are we done? */
  if(l == m) return Pmm;

  /* set P_{m+1}^m = x*(2*m+1) P_m^m */
  Pmp1_m = x*(2*m+1)*Pmm;
  /* is this all we need? */
  if(l == (m+1)) return Pmp1_m;

  /* use recursion: 
     P_l^m = ( x*(2*l-1)*P_{l-1}^m - (l+m-1)*P_{l-2}^m )/(l-m)   */
  for(k=m+2; k<=l; k++) /* k runs through all l we have */
  {
    Pkm    = (x*(2*k-1)*Pmp1_m - (k+m-1)*Pmm)/(k-m);
    Pmm    = Pmp1_m;
    Pmp1_m = Pkm;
  } /* now k=l */
  return Pkm;
}


/* use recursion recommended in numrec to get P_l */
double LegendreP(int l, double x)
{
  double Pk, P0, P1;
  int k;

  if(fabs(x) > 1.0) errorexit("Bad arguments");

  /* build P_m^m */
  P0 = 1.0;
  if(l == 0) return P0;

  /* set P_{m+1} = x P_m */
  P1 = x; // P1 = x*P0;

  /* is this all we need? */
  if(l == 1) return P1;

  /* use recursion: 
     P_l = ( x*(2*l-1)*P_{l-1} - (l-1)*P_{l-2} )/l   */
  for(k=2; k<=l; k++) /* k runs through all l we have */
  {
    Pk = (x*(2*k-1)*P1 - (k-1)*P0)/k;
    P0 = P1;
    P1 = Pk;
  } /* now k=l */
  return Pk;
}




/* compute the values of the Legendre polynomial P=P_l(x), dP=dP_l(x)/dx,
   Q(x) = P_{l+1}(x) - P_{l-1}(x), and of dQ=Q'(x) */
void Legendre_P_dP_Q_dQ(int l, double x,
                       double *P, double *dP, double *Q, double *dQ)
{
  int k;
  double P0 = 1., P1 = x,  dP0 = 0., dP1 = 1.;
  double Pk, dPk;

  if(l==0)
  {
    *P  = 1.;
    *dP = 0.;
    return;
  }

  if(l==1)
  {
    *P  = x;
    *dP = 1.;
    return;
  }

  /* use recursion relation for P_l and its deriv */
  for(k = 2; k <= l+1; k++)
  {
    Pk  = (P1 * x * (2*k-1) - P0 * (k-1))/k;
    dPk = P1 * (2*k-1) + dP0;
    if(k>=l+1) break;
    P0  = P1;
    P1  = Pk;
    dP0 = dP1;
    dP1 = dPk;
  }

  *P  = P1;
  *dP = dP1;
  *Q  = Pk - P0;
  *dQ = dPk - dP0;
}


/* integral from x to 1 of Legendre polynomial P_k.
   This is coming from
   \frac{d}{dx}[(1-x^{2})\frac{dP_k}{dx}] + k(k+1) P_k(x) = 0 */
double int_LegendreP_x_1(int k, double x)
{
  double I, P,dP, q,e;

  if(k == 0)
  {
    I = 1. - x;
  }
  else
  {
    Legendre_P_dP_Q_dQ(k,x, &P,&dP, &q,&e);
    I = dP * (1. - x*x) / (k*(k+1));
  }
  return I;
}

/* integral from a to b of Legendre polynomial P_k */
double int_LegendreP_a_b(int k, double a, double b)
{
  return int_LegendreP_x_1(k, a) - int_LegendreP_x_1(k, b);
}


/* ************************************************************************ */
/* various functions needed for Legendre Gauss-Lobatto points or nodes      */
/* ************************************************************************ */

/* NOTE: there is more in Legendre.c of titan:Archives/DG_test*/
