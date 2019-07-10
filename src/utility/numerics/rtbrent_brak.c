/* rtbrent_brak.c,  Wolfgang Tichy 6/2019 */

#include <stdio.h>
#include <math.h>
#define DBL_EPSILON 1e-14 /* machine double floating point precision */

/* find root using Brent's method with bracketing in [x1,x2]:

   It is similar to sgrid's zbrent_itsP
   It allows for par and maxits
   args and returns error code. The actual root is returned in x0.
   returns iter>=0   if ok,     iter = number of iterations done
   you need to check if iter>maxits !!!
   returns <0        if error!!! */
int rtbrent_brak(double *x0, double (*func)(double,void *par),
                 double x1, double x2, void *par, int maxits, double xacc,
                 int pr)
{
  double a = x1;
  double b = x2;
  double fa = (*func)(a,par);
  double fb = (*func)(b,par);
  double c = b;
  double fc = fb;
  double d,e, p,q,r,s, tol, xm;
  int iter;

 /* check bracket for NANs */
  if( (!isfinite(x1)) || (!isfinite(x2)) )
  {
    if(pr) printf("rtbrent_brak: Bracket is not finite!  "
                  "x1=%g x2=%g\n", x1,x2);
    return -2*maxits-2;
  }

  /* check if *x0 is within [x1,x2] */
  if((fa < 0.0 && fb < 0.0) || (fa > 0.0 && fb > 0.0))
  {
    if(pr) printf("rtbrent_brak: Root is not bracketed!\n");
    return -2*maxits-3;
  }

  /* iterate */
  for(iter=1; iter<=maxits; iter++)
  {
    if((fb < 0.0 && fc < 0.0) || (fb > 0.0 && fc > 0.0))
    {
      c  = a;
      fc = fa;
      e  = d = b - a;
    }
    if(fabs(fc) < fabs(fb))
    {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }

    /* check if c-b is very small, i.e. for convergence */
    tol = 2.*DBL_EPSILON*fabs(b) + 0.5*xacc;
    xm = 0.5*(c - b);
    if(fabs(xm) <= tol || fb == 0.0)
    {
      *x0 = b;
      return iter;
    }

    /* decide if we use interpolation or bisection */
    if(fabs(e) >= tol && fabs(fa) > fabs(fb))
    {
      /* use inverse cubic interpolation */
      s = fb/fa;
      if(a == c)
      {
        p = 2.0*xm*s;
        q = 1.0-s;
      }
      else
      {
        q = fa/fc;
        r = fb/fc;
        p = s * (2.*xm*q*(q-r) - (b-a)*(r-1.));
        q = (q-1.) * (r-1.) * (s-1.);
      }
      if(p > 0.0) q = -q;
      p = fabs(p);

      /* as in GSL's roots/brent.c, check if interpolation works */
      if( 2.0*p < fmin(3.0*xm*q - fabs(tol*q), fabs(e*q)) )
      {
        e = d;       /* interpolation is ok */
        d = p/q;
      }
      else
      {
        e = d = xm;  /* interpolation failed, fall back to bisection */
      }
    }
    else
    {
      /* use bisection */
      e = d = xm;
    }

    /* copy b into a */
    a = b;
    fa = fb;

    /* update b */
    if(fabs(d) > tol)
      b += d;
    else
      b += (xm >= 0 ? +tol : -tol);
    fb = (*func)(b,par);
  }

  /* if we get here, it did not convergence */
  if(pr) printf("rtbrent_brak: Maximum number of iterations exceeded!  "
                "iter=%d > maxits=%d\n", iter, maxits);
  *x0 = b;
  return iter;
}
