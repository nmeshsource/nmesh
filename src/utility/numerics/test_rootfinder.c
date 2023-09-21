/* test root finders */

#include <math.h>
#define finit(x1) isfinite(x1)

#include "rtbrent_brak.c"
#include "newton1d_brak.c"
#include "newton1d_fd.c"

double f(double x, void *par)
{
  double f = x*(x+1)*(x-2)*(x-3);

  printf("x=%.16g f=%.16g\n", x, f);
  return f;
}

void fdf(double x, void *par, double *func, double *dfunc)
{
  *func  = f(x, par);
  *dfunc = fd_deriv(x,*func, f, par);
}


int main()
{
  int maxits = 30;
  double tol = 1e-6;
  int pr = 1;
  int ret;
  double x1,x2;
  double x0[2];

  x1=2.5; x2=4; x0[0]=3.9;
  ret = newton1d_brak(x0, fdf, x1,x2, 0, maxits,tol,pr);
  printf("newton1d_brak ret=%d: x0[0]=%.16g\n\n", ret, x0[0]);

  x1=2.5; x2=4; x0[0]=3.9;
  ret = rtbrent_brak(x0, f, x1,x2, 0, maxits,tol,pr);
  printf("rtbrent_brak ret=%d: x0[0]=%.16g\n\n", ret, x0[0]);

  x1=1.5; x2=4; x0[0]=3.9;
  ret = newton1d_fd_region(x0, f, x1,x2, 0, maxits,tol,pr);
  printf("newton1d_fd_region ret=%d: x0[0]=%.16g\n\n", ret, x0[0]);

  x1=1.5; x2=4; x0[0]=3.9;
  ret = find_2roots_region(x0, f, x1,x2, 0, maxits,tol,pr);
  printf("find_2roots_region ret=%d: x0[0]=%.16g x0[1]=%.16g\n\n",
         ret, x0[0], x0[1]);

  x1=-0.5; x2=2.5; x0[0]=1.6;
  ret = find_2roots_region(x0, f, x1,x2, 0, maxits,tol,pr);
  printf("find_2roots_region ret=%d: x0[0]=%.16g x0[1]=%.16g\n\n",
         ret, x0[0], x0[1]);

  x1=-5; x2=5; x0[0]=1.5;
  ret = find_2roots_region(x0, f, x1,x2, 0, maxits,tol,pr);
  printf("find_2roots_region ret=%d: x0[0]=%.16g x0[1]=%.16g\n\n",
         ret, x0[0], x0[1]);

  x1=3.01; x2=7; x0[0]=6.5;
  ret = find_2roots_region(x0, f, x1,x2, 0, maxits,tol,pr);
  printf("find_2roots_region ret=%d: x0[0]=%.16g x0[1]=%.16g\n\n",
         ret, x0[0], x0[1]);

  x1=1.5; x2=4; x0[0]=3.9;
  ret = rtbrent_brak(x0, f, x1,x2, 0, maxits,tol,pr);
  printf("rtbrent_brak ret=%d: x0[0]=%.16g\n\n", ret, x0[0]);


  x1=2.0; x2=2.0; x0[0]=3.9;
  ret = rtbrent_brak(x0, f, x1,x2, 0, maxits,tol,pr);
  printf("rtbrent_brak ret=%d: x0[0]=%.16g\n\n", ret, x0[0]);

  x1=2.0; x2=2.5; x0[0]=3.9;
  ret = rtbrent_brak(x0, f, x1,x2, 0, maxits,tol,pr);
  printf("rtbrent_brak ret=%d: x0[0]=%.16g\n\n", ret, x0[0]);

  x1=1.5; x2=2.0; x0[0]=3.9;
  ret = rtbrent_brak(x0, f, x1,x2, 0, maxits,tol,pr);
  printf("rtbrent_brak ret=%d: x0[0]=%.16g\n\n", ret, x0[0]);


  x1=0.0; x2=0.0; x0[0]=3.9;
  ret = rtbrent_brak(x0, f, x1,x2, 0, maxits,tol,pr);
  printf("rtbrent_brak ret=%d: x0[0]=%.16g\n\n", ret, x0[0]);

  x1=0.0; x2=0.5; x0[0]=3.9;
  ret = rtbrent_brak(x0, f, x1,x2, 0, maxits,tol,pr);
  printf("rtbrent_brak ret=%d: x0[0]=%.16g\n\n", ret, x0[0]);

  x1=-0.5; x2=0.0; x0[0]=3.9;
  ret = rtbrent_brak(x0, f, x1,x2, 0, maxits,tol,pr);
  printf("rtbrent_brak ret=%d: x0[0]=%.16g\n\n", ret, x0[0]);
}
