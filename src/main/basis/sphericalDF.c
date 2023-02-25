


/* sphericalDF arrays have data with
   d[Ind_n(i,j,k,n)] is such that
   i <--> theta
   j <--> phi
   k <--> 3rd dim. e.g. for integral against several sph. harmonics
*/


/* we need arrays for
   -data
   -Re and Im of Y_l^m(theta,phi) for several l,m
   -transpose of ana/syn matrices At, St for Fourier in theta_i
   -theta, phi coords
  don't need
   -integr. matrix for phi-Fourier
     -> instead use sum over phi_j top get c_0
 */


/*
sphericalDF arrays uijk have data with
   d[Ind_n(i,j,k,n)] is such that
   i <--> theta
   j <--> phi
   k <--> 3rd dim. e.g. for integral against several sph. harmonics

1st implement  sphericalDF2dIntegral  over a 3d-array
I.e. we integrate over a 2d array
+do phi-integral for every theta -> 1d array(theta)
+integrate 1d-array(theta), for this we need ana matrix At for Fourier
*/

/* In:  3d array uijk = u(theta_i, phi_j, k)
        Here theta_i = 2*PI*i/N0 + PI/((1+N0%2)*N0)  N0 = n_theta
             phi_j   = 2*PI*j/N1                     N1 = n_phi
   Out: 1d array Uk with 2d integral uijk over theta and phi for each k */
void sphericalDF_2dIntegral(tArray *uijk, tArray *Uk)
{
  int *n = Arrn(uijk);
  double *u = Arrd(uijk);
  double *U = Arrd(Uk);
  tArray *Cik = alloc_array2d(n[0], n[2]);
  double *C = Arrd(Cik);
  /* make ana matrix At for Fourier in theta-dir */
  // call func

  /* integrate over phi, and write result into 2d arr Uik */
  for(k=0; k<n[2]; k++)
  {
    for(i=0; i<n[0]; i++)
    {
      double c0;
      /* sum uijk over j: sum -> Uik */
      c0 = 0.;
      for(j=0; j<n[1]; j++)  c0 += u[Ind_n(i,j,k,n)];
      C[i + n[2]*k] = c0;
    }

    /* get Fourier coeffs Cik using At, maybe can use Uik to store Cik */
    mm_array0_norestrict(At, Cik, ACik);

    // use coeffs as in spec_sphericalDF2dIntegral_at_radial_index_i
    // to find integral, write into Uk
  }
  free_array(Cik);
}
