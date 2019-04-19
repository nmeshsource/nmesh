/* numflux1d.c */
/* Wolfgang Tichy, April 2019 */


#include "nmesh.h"
#include "dg.h"


/* several numerical fluxes in 1d */



/* flux at interface of 1d scalar Godunov method for Burgers eqn */
double numflux1d_scalarGodunov(int nf, double *f,
                               double *uL, double *uR, double *fL, double *fR,
                               double *lamL, double *lamR)
{
  int i;

  for(i=0; i<nf; i++)
  {
    double ul = uL[i];
    double ur = uR[i];
    double fi;

    if(ul >= 0. && ur >= 0.)      fi = fL[i];
    else if(ul <= 0. && ur <= 0.) fi = fR[i];
    else if(ul  < 0. && ur >= 0.) fi = 0.;
    else
    {
      double s2 = ul + ur;
      if(s2 > 0.) fi = fL[i];
      else        fi = fR[i];
    }
    f[i] = fi;
  }
}


/* LLF flux for nf fields. Numerical flux is written in f[nf] */
void numflux1d_LLF(int nf, double *f,
                   double *uL, double *uR, double *fL, double *fR,
                   double *lamL, double *lamR)
{
  double amax;
  double apl = 0.;
  double ami = 0.;
  int i;

  /* max speeds */
  for(i=0; i<nf; i++)
  {
    apl = max2(apl, max2(lamR[i], lamL[i]));
    ami = min2(ami, min2(lamR[i], lamL[i]));
  }

  apl = fabs(apl);
  ami = fabs(ami);
  amax = max2(ami,apl);

  /* build fluxes */
  for(i=0; i<nf; i++) 
    f[i] = 0.5*( fR[i] +  fL[i] - amax*(uR[i] -  uL[i]) );
}


#define TINY 1e-32

/* HLL flux */
void numflux1d_HLL(int nf, double *f,
                   double *uL, double *uR, double *fL, double *fR, 
                   double *lamL, double *lamR)
{
  double oda, apm;
  double apl = 0.;
  double ami = 0.;
  int i;

  /* max speeds */
  for(i=0; i<nf; i++)
  {
    apl = max2(apl, max2(lamR[i], lamL[i]));
    ami = min2(ami, min2(lamR[i], lamL[i]));
  }

  apl = fabs(apl);
  ami = fabs(ami);
  
  oda = 1./( apl + ami + TINY);
  apm = apl * ami;

  /* build fluxes */
  for(i=0; i<nf; i++) 
    f[i] = oda * ( apl * fL[i] + ami * fR[i] - apm * ( uR[i] -  uL[i]) );
}
