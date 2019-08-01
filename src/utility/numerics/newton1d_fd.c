/* newton1d_fd.c,  Wolfgang Tichy 6/2019 */

#include <stdio.h>
#include <math.h>


/* funcs needed */
int rtbrent_brak(double *x0, double (*func)(double,void *par),
                 double x1, double x2, void *par, int maxits, double xacc,
                 int pr);


#define EPS  1e-7   /* should be approx square root of machine precision */
#define HMIN 1e-10  /* min h we use in fin. diff. computation of derivs */

/* fintite diff. deriv of func */
double fd_deriv(double x0, double f0,
                double (*func)(double x, void *par), void *par)
{
  double h, f, x, df;
  int k;

  /* fd resolution */
  h = EPS*fabs(x0);

  /* make sure h is never below HMIN */
  if(h < HMIN) h = HMIN;

  /* increase h at most by 8 orders of mag. */
  for(k=0; k<8; k++, h=h*10)
  {
    x = x0 + h;
    h = x - x0; /* better value for h */
    f = (*func)(x, par);

    df = (f - f0)/h;
    if(df != 0.0) break;
  }
  return df;
}


/* Newton-Raphson within a region [x1,x2]:

   It allows for par and maxits args and returns error code.
   The actual root is returned in x0.
   returns j>=0   if ok,        j = number of iterations done
   returns <0     if error!!!
   *x0 is used as starting guess */
int newton1d_fd_region(double *x0, double (*func)(double x, void *par),
                       double x1, double x2, void *par,
                       int maxits, double xacc, int pr)
{
  int j, badstep;
  double f, df;
  double xrt, dx,dxold;
  double tmp;

  /* check bracket */
  if( (!isfinite(x1)) || (!isfinite(x2)) )
  {
    if(pr) printf("newton1d_fd_region: Region is not finite!  "
                  "x1=%g x2=%g\n", x1,x2);
    return -2*maxits-2;
  }

  /* check if *x0 is within [x1,x2] */
  if( (x1<=*x0 && *x0<=x2) || (x2<=*x0 && *x0<=x1) )
  {
    xrt=*x0;
  }
  else /* we also get this if *x0 is NAN */
  {
    if(pr) printf("newton1d_fd_region: x0 is not inside region  "
                  "x1=%g x2=%g x0=%g\n", x1, x2, *x0);
    return -2*maxits-3;
  }

  /* get func values */
  f  = (*func)(*x0, par);

  /* early returns */
  if(f == 0.0) { return 0; }

  /* set dx ... */
  dxold = fabs(x2-x1);
  dx    = dxold;

  /* iterate */
  for(j=1; j<=maxits; j++)
  {
    /* eval func at new x=xrt */
    f = (*func)(xrt, par);
    df = fd_deriv(xrt, f, func, par);

    /* catch special cases */
    if(!isfinite(f)) { return -j; }
    if(f == 0.0) { *x0=xrt; return j; }

    /* test if Newton step could be possible */
    if( !isfinite(df) ) badstep = 1;
    else                badstep = 0;

    if(badstep==0) /* try Newton step */
    {
      dxold = dx;
      dx  = f/df;
      tmp = xrt;
      xrt -= dx;
      /* check if Newton step brings us outside interval */
      if( ((xrt-x2)*(xrt-x1) > 0.0) )
      {
        badstep = 1;
        xrt = tmp;
        dx = dxold;
      }
    }

    /* update x0 */
    *x0 = xrt;
    if(badstep) { return -j; }

    /* check accuracy goal */
    if(fabs(dx) < xacc) { return j; }
  }

  if(pr) printf("newton1d_fd_region: Maximum number of iterations exceeded!  "
                "j=%d > maxits=%d\n", j, maxits);
  return j;
}


/* find 2 roots within a region [x1,x2]:

   It allows for par and maxits args and returns error code.
   The actual roots are returned in x0[0] and x0[1].
   returns 2 if 2 roots are found
   returns 1 if [x1,x2] is a valid bracket, then x0[0]=x0[1] are this root
   returns <0     if error!!!
   x0[0] is used as starting guess */
int find_2roots_region(double x0[2],
                       double (*func)(double x, void *par),
                       double x1, double x2, void *par,
                       int maxits, double xacc, int pr)
{
  int j, ret;
  double f, fh,fl, df;
  double xh,xl, xrt, dx, xmid;

  /* check bracket */
  if( (!isfinite(x1)) || (!isfinite(x2)) )
  {
    if(pr) printf("find_2roots_region: Region is not finite!  "
                  "x1=%g x2=%g\n", x1,x2);
    return -2*maxits-2;
  }

  /* check if x0[0] is within [x1,x2] */
  if( (x1<=x0[0] && x0[0]<=x2) || (x2<=x0[0] && x0[0]<=x1) )
  {
    xrt=x0[0];
  }
  else /* we also get this if x0[0] is NAN */
  {
    if(pr) printf("find_2roots_region: x0 is not inside region  "
                  "x1=%g x2=%g x0=%g\n", x1, x2, x0[0]);
    return -2*maxits-3;
  }

  /* get func values */
  fl = (*func)(x1, par);
  fh = (*func)(x2, par);

  /* check if root is bracketed and if func values are not NAN */
  j=0;
  if(fl <= 0.0 && fh >= 0.0) j=1; /* if fl or fh is NAN j stays 0 */
  if(fl >= 0.0 && fh <= 0.0) j=1; /* if fl or fh is NAN j stays 0 */
  if(j)
  {
    /* root is bracketed */
    ret = rtbrent_brak(&(x0[0]), func, x1,x2, par, maxits, xacc, pr);
    if(ret<0) return ret; /* error */

    /* assume both roots are the same */
    x0[1] = x0[0];
    return 1; /* only 1 root found */
  }

  /* iterate to find 1 root */
  for(j=1; j<=maxits; j++)
  {
    ret = newton1d_fd_region(&xrt, func, x1,x2, par, maxits, xacc, pr);
    if(ret>=0)
    {
      x0[0] = xrt;
      break;
    }
    else
    {
      xrt = 0.333333333333333333*(x1 + xrt + x2);
      if(xrt<x1 || xrt>x2)
        xrt = 0.5*(x1 + x2);
    }
  }

  /* fail if there is not at least one root */
  if(ret<0)
  {
    if(pr) printf("find_2roots_region: newton1d_fd_region failed ret=%d: "
                  "x1=%g x2=%g xrt=%g\n", ret, x1, x2, xrt);
    return -maxits; /* no root found */
  }

  /* we found one root x0[0], now make bracket to find the other */
  f = (*func)(x0[0], par);

  if(f*fl<0)
  {
    xmid = x0[0];
  }
  else
  {
    df = fd_deriv(x0[0], f, func, par);

    /* 3 times Newton step should overshoot and change sign */
    dx  = 3.*(f + fl*EPS)/df;
    xmid = x0[0] - dx;
  }

  if(0) printf("find_2roots_region: xmid=%.16g\n", xmid);
  if(0) printf("find_2roots_region: x1=%.16g x2=%.16g\n", x1,x2);

  /* brackets are [x1,xmid] and [xmid,x2] */
  if(x0[0] <= xmid) { xl=xmid; xh=x2; }
  else              { xl=x1;   xh=xmid; }

  if(0) printf("find_2roots_region: xl=%.16g xh=%.16g\n", xl,xh);

  /* use root finder with bracket */
  ret = rtbrent_brak(&(x0[1]), func, xl,xh, par, maxits, xacc, pr);
  if(ret<0) return ret; /* error */

  return 2; /* found 2 roots */
}
