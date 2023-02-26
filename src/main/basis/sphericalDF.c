/* sphericalDF.c */
/* Wolfgang Tichy 2/2023 */

#include "nmesh.h"
#include "basis.h"


/* sphericalDF arrays have data with
   d[Ind_n(i,j,k,n)] is such that
   i <--> theta \in [0,2pi)
   j <--> phi   \in [0,2pi)
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
   But we still compute:
   I_k = \int_0^{pi) dtheta \int_0^{2pi) dphi u(theta,phi,k) |sin(theta)|
   In:  3d array auijk = u(theta_i, phi_j, k)
        Here theta_i = 2*PI*i/n0 + PI/((1+n0%2)*n0)  n0 = n_theta
             phi_j   = 2*PI*j/n1                     n1 = n_phi
   Out: 1d array aUk with 2d integral I_k for each k
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
      Uik[i + n0*k] = c0 * L/n1;
    }
  }

  /* make analysis matrix At for Fourier in theta-dir */
  set_TrafoArray(At, Fourier_coeffs);

  /* get Fourier coeffs aCik using At */
  mm_array0(At, aUik, aCik);

  ////printf("auijk"); printarray(auijk);
  //printf("At"); printarray(At);
  //printf("aUik"); printarray(aUik);
  //printf("aCik"); printarray(aCik);

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

    /* Note:
       if we expand in terms of coeffs we get:
       u(theta) = aa_0
                 +2\sum_{n=1}^{N-1}[ aa_n cos(n theta) + bb_n sin(n theta) ]
       we need I = \int_0^{\pi} u(theta) sin(theta) dtheta
       Integrate[Sin[x], {x, 0,Pi}]
        = 2
       Simplify[ Integrate[Sin[x]*Sin[x], {x, 0,Pi}] ]
        = Pi/2
       Simplify[ Integrate[Sin[n*x]*Sin[x], {x, 0,Pi}] ]
        = 0
       Simplify[ Integrate[Cos[n*x]*Sin[x], {x, 0,Pi}] ]
        = (1 + (-1)^n)/(1 - n^2)
       ==>
       I = aa_0 * 2 + 2bb_1 * Pi/2 + 2\sum_{even n} aa_n * 2/(1 - n^2)
       ---------------------------------------------------------------
       BUT we really expand in thm not theta = thm + PI/((1+N%2)*N)
       Thus we get the coeffs
       c_n := a_n + i b_n
       and
       u(thm) = c_0 + 2\sum_{n=1}^{N-1} c_n e^{-i n thm}
              = c_0 + 2\sum_{n=1}^{N-1} c_n e^{-i n theta} e^{i 2Pi d n}
       where d = 1./(2.*(1+N%2)*N)
       Thus
       u(theta) = cc_0 + 2\sum_{n=1}^{N-1} cc_n e^{-i n theta}
       with cc_n = e^{i 2Pi d n} c_n =: aa_n + i bb_n
       ==> aa_n = cos(2Pi d n) Re_c_n - sin(2Pi d n) Im_c_n
           bb_n = sin(2Pi d n) Re_c_n + cos(2Pi d n) Im_c_n   */

    /* sum over all theta-integrated terms */
    sum = (1.0/PI) * Cik[0 + n0*k];
    if(N >= 3)
      sum += 0.5*( sin(PI2*d) * Cik[1 + n0*k]
                  +cos(PI2*d) * Cik[2 + n0*k] );
    for(n=2;n<N/2;n+=2)
    {
      Re_c_n = Cik[2*n-1 + n0*k];  /* c[2*n-1]; */
      Im_c_n = Cik[2*n   + n0*k];  /* c[2*n];   */
      sum += 2.0*( cos(PI2*n*d)/((1-n*n)*PI) * Re_c_n
                  +sin(PI2*n*d)/((n*n-1)*PI) * Im_c_n );
    }
    if( N%4 == 0 )
      sum += cos(PI2*(N/2)*d)/((1-N*N/4)*PI) * Cik[N-1 + n0*k];

    /* adjust sum for L and N to obtain integral over theta */
    sum *= L/N;

    /* write integral into Uk */
    Uk[k] = sum;
  }
  free_array(At);
  free_array(aCik);
  free_array(aUik);
}


/* copy sphericalDF array into double covered regions */
void sphericalDF_copy_to_doubleCoveredPoints(tArray *AsDF)
{
  double *arr = Arrd(AsDF);
  int *n = Arrn(AsDF);
  int n0 = n[0], n1 = n[1], n2 = n[2];
  int k;

  /* check if we can copy data into double covered regions */
  if( n0%2 || n1%2 )
    errorexit("n[0] and n[1] must be even!");

  /* copy arr into double covered regions */
  for(k = 0; k < n2; k++)
  {
    int i,j;
    for(j = 0;    j < n1/2; j++)
      for(i = n0/2; i < n0; i++)
        arr[Ind_n(i,j,k ,n)] = arr[Ind_n(n0-i-1,j+n1/2,k ,n)];

    for(j = n1/2; j < n1; j++)
      for(i = n0/2; i < n0; i++)
        arr[Ind_n(i,j,k ,n)] = arr[Ind_n(n0-i-1,j-n1/2,k ,n)];
  }
}


void sphericalDF_test(void)
{
  int n[3];
  tArray *aF;
  tArray *aI;
  double *f;
  int i,j,k;

  n[0] = n[1] = 6;
  n[2] = 6;
  aF = alloc_array(n);
  aI = alloc_array1d(n[2]);
  f = Arrd(aF);

  for(k=0; k<n[2]; k++)
  {
    for(j=0; j<n[1]; j++)
    for(i=0; i<n[0]; i++)
    {
      double th,ph;
      sphericalDF_theta_phi(i,j, n, &th, &ph);
      //thm=i*2*PI/n[0];

      switch(k)
      {
      case 0:
        f[Ind_n(i,j,k,n)] = (i+1) + 10*(j+1);
        break;
      case 1:
        f[Ind_n(i,j,k,n)] = 1./(4*PI);
        break;
      case 2:
        f[Ind_n(i,j,k,n)] = cos(ph)*cos(ph)/(4*PI);
        break;
      case 3:
        f[Ind_n(i,j,k,n)] = sin(th)*sin(th)*cos(ph)*cos(ph)/(PI);
        break;
      case 4:
        f[Ind_n(i,j,k,n)] = cos(th)/(4*PI);
        break;
      case 5:
        f[Ind_n(i,j,k,n)] = cos(th)*cos(th)/(2*PI);
        break;
      }
    }
  }

  printf("f");printarray(aF);
  sphericalDF_copy_to_doubleCoveredPoints(aF);
  printf("f");printarray(aF);

  sphericalDF_2dIntegral(aF, aI);
  printf("I");printarray(aI);

  free_array(aI);
  free_array(aF);
}
