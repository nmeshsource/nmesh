/* Fourier.c */
/* Wolfgang Tichy 4/2004 in sgrid's explicit_Four_trafos.c */

#include "nmesh.h"
#include "basis.h"



/* Note here we have X = [a,b] and Xb=[-1,1]
   Def; L = b-a
   X = (L/2)*(Xb+1) + a   so that: Xb=-1 => X=a,  Xb=1 => X=b
   Xb = 2(X-a)/L - 1
*/

/* compute Four coeffs of deriv cder[0...N-1] from Four coeffs c[0...N-1] */
void Fourier_deriv(int N, const double c[], double cder[], double L)
{
  int j;
  double PI2_con;

  PI2_con = 2.0*PI/L;

  for(j=1; 2*j<N; j++)
  {
    cder[2*j-1] =  c[2*j] * PI2_con*j;
    cder[2*j]   = -c[2*j-1] * PI2_con*j;
  }
  cder[0] = 0.0;
  if( N%2 == 0 ) cder[N-1] = 0.0;

  //else  /* this is what I had in the old explicit_Four_trafos.c version */
  //{ cder[N-2] = 0.0;  cder[N-1] = 0.0; }
}


/* compute Four coeffs of integral cint[0...N-1] from Cheb coeffs c[0...N-1] */
void Fourier_int(int N, const double c[], double cint[], double L)
{
  int j;
  double PI2_con;
  double *u = (double*) calloc(N+1, sizeof(double));

  PI2_con = 2.0*PI/L;

  /* get terms coming from integrating c[0] */
/*
  for(j=1; 2*j<N; j++)
  {
    cint[2*j-1] = -0.5*L*c[0]/((double) N);
    cint[2*j]   = -0.5*L*c[0]/((double) N); // WRONG!!!!!
    // integrate the func 1 and find its coeffs instead!!!
    // multiply this by c[0]
    if(N!=4) errorexit("Fourier_int is wrong");
  }
  if( N%2 == 0 ) cint[N-1] = -0.5*L*c[0]/((double) N);
  cint[0] = 0.5*n*L*c[0]/((double) N);
*/
  for(j=0; j<N; j++) u[j]=j*L/N;
  Fourier_coeffs(N, cint, u); /* get coeffs of the integral of 1 */
  for(j=0; j<N; j++) cint[j] *= c[0]/N;

  /* add terms coming from integrating everything but the c[0] term */
  for(j=1; 2*j<N; j++)
  {
    cint[2*j-1] += -c[2*j] / (PI2_con*j);
    cint[2*j]   +=  c[2*j-1] / (PI2_con*j);
  }
  if( N%2 == 0 ) cint[N-1] += 0.0;

  /* free temp mem u */
  free(u);
}


/* compute Four coeffs c[0...N-1] from function u
   at x_k = k/N, k=0,...,N-1
NOTE: Fourier_coeffs returns c[] that are N times of those of four_coeffs_alt */
void Fourier_coeffs(int N, double c[], const double u[])
{
  int k, j;
  double Re_c_j, Im_c_j, PI2oN;

  PI2oN=2.0*PI/N;

  /* first j=0 , i.e. c_0 */
  c[0] = 0.0;
  for(k=0;k<N;k++)  c[0] += u[k];

  for(j=1; j<=N/2; j++)
  {
    Re_c_j = Im_c_j = 0.0;
    for(k=0;k<N;k++)
    {
      Re_c_j += cos(j*PI2oN*k)*u[k];
      Im_c_j += sin(j*PI2oN*k)*u[k];
    }
    c[2*j-1] = Re_c_j; /* real part of c_j */
    if(2*j<N)
      c[2*j] = Im_c_j; /* imaginary part of c_j */
  }
  /* we should use a FFT for everything above this line */
  /* see four_coeffs_FFTW3 in sgrid/src/utility/Spectral/FFTs_for_sgrid.c */
}


/* find function u from Four coeffs c[0...N-1], computed with Fourier_coeffs */
void Fourier_eval(int N, const double c[], double u[])
{
  int k, j;
  double sum, Re_c_k, Im_c_k, PI2oN;

  PI2oN=2.0*PI/N;

  for(j=0; j<N; j++)
  {
    sum = 0.0;
    for(k=1;k<=N/2;k++)
    {
      Re_c_k   = c[2*k-1];
      if(2*k<N)
        Im_c_k = c[2*k];
      else
        Im_c_k = 0.0;
      sum += cos(k*PI2oN*j)*Re_c_k + sin(k*PI2oN*j)*Im_c_k;
    }
    if( N%2 == 0 )
      u[j] = ( c[0] + 2.0*sum - c[N-1]*cos(PI*j) )/N;
    else
      u[j] = ( c[0] + 2.0*sum )/N;
  }
  /* we should use a FFT for eveything above this line */
  /* see four_eval_FFTW3 in sgrid/src/utility/Spectral/FFTs_for_sgrid.c */
}


/* find value of Fourier basis function B_k at X (in [a,b]) */
double Fourier_basisfunc(int N, int k, double X, double L)
{
  double K = 2.0*PI/L;
  int j = k/2 + k%2;

  if(k==0)              return 1.0/N;
  if(k==N-1 && N%2==0)  return cos(j*K*X)*1.0/N;;
  if(k%2!=0)            return cos(j*K*X)*2.0/N;
  return sin(j*K*X)*2.0/N;
}
