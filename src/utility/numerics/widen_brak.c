#include <stdio.h>
#include <math.h>

/* Widen bracket [x1,x2] until it contains a zero of func(x,par), using
   factor Fac > 0.
   It is similar to sgrid's zbrac_P (but Fac does not need to be 1.6)
   returns number of widenings it   if ok
   returns <0                       if failure */
int widen_brak_f(double (*func)(double,void *par),
                 double *x1, double *x2, void *par, double Fac, int ntries,
                 int pr)
{
  int it;
  double f1,f2;

  if(*x1 == *x2)
  {
    if(pr) printf("widen_brak_f: x1 and x2 must be unequal\n");
    return -ntries-2;
  }

  f1 = (*func)(*x1, par);
  f2 = (*func)(*x2, par);

  for(it=0; it<ntries; it++)
  {
    if(pr==2) printf("widen_brak_f: it%d x1=%g x2=%g\n", it, *x1, *x2);

    if(f1*f2 < 0.0)
      return it; /* bracket is good now */

    if(fabs(f1) < fabs(f2))
    {
      *x1 -= Fac*(*x2-*x1); /* widen on left */
      f1 = (*func)(*x1, par);
    }
    else
    {
      *x2 += Fac*(*x2-*x1); /* widen on right */
      f2 = (*func)(*x2, par);
    }
  }
  if(pr)
    printf("widen_brak_f: Could not find bracket after %d tries.\n", ntries);
  return -it-1;
}

/* widen bracket with standard factor of 1.6 */
int widen_brak(double (*func)(double,void *par),
               double *x1, double *x2, void *par, int ntries, int pr)
{
  return widen_brak_f(func, x1,x2, par, 1.6, ntries, pr);
}

/* Note: see also widen_brak_f_fdf in rtbrent_brak.c */
