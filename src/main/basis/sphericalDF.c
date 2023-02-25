/* sphericalDF.c */
/* Wolfgang Tichy 2/2023 */

#include "nmesh.h"
#include "basis.h"


/* sphericalDF arrays have data with
   d[Ind_n(i,j,k,n)] is such that
   i <--> theta
   j <--> phi
   k <--> 3rd dim. e.g. for integral against several sph. harmonics */



/* set theta_i, phi_j for sphericalDF */
void sphericalDF_theta_phi(int i, int j, const int *n,
                           double *theta_i, double *phi_j)
{
  int N = n[0];
  double thm = 2*PI*i/N;
  *theta_i = thm + PI/((1+N%2)*N);
  *phi_j   = 2*PI*j/n[1];
}


/* Integrate over a 2-sphere where we have double covering, because
   0 <= theta < 2pi
   In:  3d array auijk = u(theta_i, phi_j, k)
        Here theta_i = 2*PI*i/n0 + PI/((1+n0%2)*n0)  n0 = n_theta
             phi_j   = 2*PI*j/n1                     n1 = n_phi
   Out: 1d array aUk with 2d integral uijk over theta and phi for each k
   Steps we take to  integrate over a 2d array:
   -do phi-integral for every theta -> 1d array(theta)
   -integrate 1d-array(theta), for this we need ana matrix At for Fourier */
void sphericalDF_2dIntegral(tArray *auijk, tArray *aUk)
{
  int *un = Arrn(auijk);
  int n0 = un[0], n1 = un[1], n2 = un[2];
  tArray *aUik = alloc_array2d(n0, n2);
  tArray *aCik = alloc_array2d(n0, n2);
  tArray *At = alloc_array2d(n0,n0);
  double *uijk = Arrd(auijk);
  double *Uik  = Arrd(aUik);
  double *Cik  = Arrd(aCik);
  double *Uk   = Arrd(aUk);
  double PI2 = 2.0*PI;
  double L;
  int k;

  /* integrate over phi, and write result into 2d array aUik */
  L = PI2; /* assume phi goes from 0 to 2Pi */
  for(k=0; k<n2; k++)
  {
    int i, j;
    for(i=0; i<n0; i++)
    {
      double c0;
      /* sum uijk over j: sum -> Uik */
      c0 = 0.;
      for(j=0; j<n1; j++) c0 += uijk[Ind_n(i,j,k,un)];
      Uik[i + n2*k] = c0 * L/n1;
    }
  }

  /* make analysis matrix At for Fourier in theta-dir */
  set_TrafoArray(At, Fourier_coeffs);

  /* get Fourier coeffs aCik using At */
  mm_array0(At, aUik, aCik);

  /* use coeffs in aCik to find integral over theta for every k */
  L = PI2; /* assume theta goes from 0 to 2Pi */
  for(k=0; k<n2; k++)
  {
    double sum;
    int n;
    int N = n0;
    /* double theta = thm + PI/((1+N%2)*N); // thm = 2*PI*i/N0 */
    double d = 1.0/(2.0*(1+N%2)*N);
    double Re_c_n, Im_c_n;

    /* sum over all theta-integrated terms */
    sum = (1.0/PI) * Cik[0 + n2*k];
    sum += 0.5*( sin(PI2*d) * Cik[1 + n2*k]
                +cos(PI2*d) * Cik[2 + n2*k] );
    for(n=2;n<N/2;n+=2)
    {
      Re_c_n = Cik[2*n-1 + n2*k];  /* c[2*n-1]; */
      Im_c_n = Cik[2*n   + n2*k];  /* c[2*n];   */
      sum += 2.0*( cos(PI2*n*d)/((1-n*n)*PI) * Re_c_n
                  +sin(PI2*n*d)/((n*n-1)*PI) * Im_c_n );
    }
    if( N%4 == 0 )
      sum += cos(PI2*(N/2)*d)/((1-N*N/4)*PI) * Cik[N-1 + n2*k];

    /* adjust sum for L and N to obtain integral over theta */
    sum *= L/N;

    /* write integral into Uk */
    Uk[k] = sum;
  }
  free_array(At);
  free_array(aCik);
  free_array(aUik);
}
