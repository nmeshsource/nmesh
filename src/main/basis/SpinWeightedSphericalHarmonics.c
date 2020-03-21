/* SpinWeightedSphericalHarmonics.c */
/* Wolfgang Tichy 2/2010 in sgrid */

#include "nmesh.h"
#include "basis.h"


/*
We use the same spin weighted spherical harmonics as:
-arXiv:0709.0093v2 [gr-qc] at http://arxiv.org/abs/0709.0093
-arXiv:gr-qc/0610128
-astro-ph/0508514v3
NOTE: The spin weight n of astro-ph/0508514v3 is denoted by s here. i.e. s=n.
BUT arXiv:0709.0093v2 and arXiv:gr-qc/0610128 seem to use s=-n.
In bam s=-n. So for spin weight n=-2 bam has s=2, while sgrid has s=-2.
*/

/* factorial */
double fact(double n)
{
  if(n<=1)
    return 1.0;
  else
  {
    n=n*fact(n-1);
    return n;
  }
}

/* Wigner d function, coded by WT */
double Wigner_d_function_WT(int l, int m, int s, double theta)
{
  double Wigd = 0;
  int k1 = s > m  ? 0 : m-s;
  int k2 = m < -s ? l+m : l-s;
  double costhetao2 = cos(theta*0.5);
  double sinthetao2 = sin(theta*0.5);
  int k;

  if(k1%2==0)  /* k1 even */
  {
    for(k = k1; k <= k2; k+=2)
      Wigd += pow(costhetao2, 2*l+m-s-2*k) * pow(sinthetao2, 2*k+s-m) /
              ( fact(l+m-k) * fact(l-s-k) * fact(k) * fact(k+s-m) );
    for(k = k1+1; k <= k2; k+=2)
      Wigd -= pow(costhetao2, 2*l+m-s-2*k) * pow(sinthetao2, 2*k+s-m) /
              ( fact(l+m-k) * fact(l-s-k) * fact(k) * fact(k+s-m) );
  }
  else
  {
    for(k = k1; k <= k2; k+=2)
      Wigd -= pow(costhetao2, 2*l+m-s-2*k) * pow(sinthetao2, 2*k+s-m) /
              ( fact(l+m-k) * fact(l-s-k) * fact(k) * fact(k+s-m) );
    for(k = k1+1; k <= k2; k+=2)
      Wigd += pow(costhetao2, 2*l+m-s-2*k) * pow(sinthetao2, 2*k+s-m) /
              ( fact(l+m-k) * fact(l-s-k) * fact(k) * fact(k+s-m) );
  }
  return (sqrt(fact(l+m) * fact(l-m) * fact(l+s) * fact(l-s)) * Wigd);
}

/* real part of spin-weighted spherical harmonic sYlm */
double Re_sYlm(int l, int m, int s, double theta, double phi)
{
  double c = sqrt( (2.0*l+1)/(4*PI) );
  if(s%2 != 0) c=-c; /* multiply by (-1)^s */
  return c*Wigner_d_function_WT(l, m, -s, theta) * cos(m*phi);
}

/* imaginary part of spin-weighted spherical harmonic sYlm */
double Im_sYlm(int l, int m, int s, double theta, double phi)
{
  double c = sqrt( (2.0*l+1)/(4*PI) );
  if(s%2 != 0) c=-c; /* multiply by (-1)^s */
  return c*Wigner_d_function_WT(l, m, -s, theta) * sin(m*phi);
}
