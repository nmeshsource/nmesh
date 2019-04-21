/* numflux1d.c */
/* Wolfgang Tichy, April 2019 */


#include "nmesh.h"
#include "dg.h"


/* several numerical fluxes in 1d */



/* flux at interface of 1d scalar Godunov method for Burgers eqn */
void numflux1d_scalarGodunov(tDGinfo *d)
{
  int isP = d->face%2;
  int nf = d->vlu->n;
  int i;

  for(i=0; i<nf; i++)
  {
    double ul,ur, fl,fr, fn;

    if(isP)
    {
      ul = d->ui[i];
      ur = d->ua[i];
      fl = d->fi[i];
      fr = d->fa[i];
    }
    else
    {
      ul = d->ua[i];
      ur = d->ui[i];
      fl = d->fa[i];
      fr = d->fi[i];
    }

    if(ul >= 0. && ur >= 0.)      fn = fl;
    else if(ul <= 0. && ur <= 0.) fn = fr;
    else if(ul  < 0. && ur >= 0.) fn = 0.;
    else
    {
      double s2 = ul + ur;
      if(s2 > 0.) fn = fl;
      else        fn = fr;
    }
    d->fnum[i] = fn;
  }
}


/* Upwind flux: Numerical flux is written in fnum[i] */
void numflux1d_upwind(tDGinfo *d)
{
  int nf = d->vlu->n;
  int i;

  /* build fluxes */
  for(i=0; i<nf; i++)
  {
    if(d->lami[i] < 0.) d->fnum[i] = d->fa[i];
    else                d->fnum[i] = d->fi[i];
  }
}


/* LLF flux for nf fields. Numerical flux is written in fnum[i] */
void numflux1d_LLF(tDGinfo *d)
{
  int nf = d->vlu->n;
  double amax;
  double apl = 0.;
  double ami = 0.;
  int i;

  /* max speeds */
  for(i=0; i<nf; i++)
  {
    apl = max2(apl, max2(d->lama[i], d->lami[i]));
    ami = min2(ami, min2(d->lama[i], d->lami[i]));
  }

  apl = fabs(apl);
  ami = fabs(ami);
  amax = max2(ami,apl);

  /* build fluxes */
  for(i=0; i<nf; i++) 
    d->fnum[i] = 0.5*( d->fa[i] +  d->fi[i] - amax*(d->ua[i] - d->ui[i]) );
}


#define TINY 1e-32

/* HLL flux */
void numflux1d_HLL(tDGinfo *d)
{
  int nf = d->vlu->n;
  double oda, apm;
  double apl = 0.;
  double ami = 0.;
  int i;

  /* max speeds */
  for(i=0; i<nf; i++)
  {
    apl = max2(apl, max2(d->lama[i], d->lami[i]));
    ami = min2(ami, min2(d->lama[i], d->lami[i]));
  }

  apl = fabs(apl);
  ami = fabs(ami);
  
  oda = 1./( apl + ami + TINY);
  apm = apl * ami;

  /* build fluxes */
  for(i=0; i<nf; i++) 
    d->fnum[i] = oda * ( apl * d->fi[i] + ami * d->fa[i] -
                         apm * (d->ua[i] - d->ui[i]) );
}
