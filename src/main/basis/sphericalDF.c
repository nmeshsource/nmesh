


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
  /* make ana matrix At for Fourier in theta-dir */
  // call func

  /* integrate over phi, and write result into 2d arr Uik */
  for(k
  {
    for(i
    {
      /* sum uijk over j: sum -> Uik */
    }

    /* get Fourier coeffs Cik using At, maybe can use Uik to store Cik */

    // use coeffs as in spec_sphericalDF2dIntegral_at_radial_index_i
    // to find integral, write into Uk
  }
}
