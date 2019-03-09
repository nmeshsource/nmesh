/* FSurf_CubedSpheres.c */
/* Wolfgang Tichy, Oct 2018 */
/* functions to to get sigma on surface of cubed sphere */

#include "nmesh.h"
#include "coordinates.h"


/* global vars in this file */
int lmax;        /* max l in Ylm expansion */

/* funcs in this file */
int FSurf_CubSph_sigma01_func(tPat *pat, int si, double AB[2], double *sig);



/* return value of surface function sigma01 */
int FSurf_CubSph_sigma01_func(tPat *pat, int si, double AB[2], double *sig)
{
  tArray *Co = pat->CI->Fcoef[si];
  int N = Co->N;
  int nYs = N/2; /* number of Ylm's we use */
  int lmax = (sqrt(8*nYs + 1) - 3)/2;
  /* We need (lmax*(lmax+1))/2 + lmax+1  complex numbers at each point A,B
     to store the table.
     So when is (lmax*(lmax+1))/2 + lmax+1 = n0?
     set L = lmax ==> L^2/2 + 3L/2 + 1 = n0  <==> L^2 + 3 L + 2 - 2*n0 = 0
     so: 2L = -3 +- sqrt(9 - 4*(2 - 2*n0)) = -3 +- sqrt(8*n0 + 1) 
     L = (sqrt(8*n0 + 1) - 3)/2  */
  int ijk, l,m;
  double sm;
  double *co = Co->d;
  double fv, Theta,Phi;
  double *ReYtab = alloc_Plm_Tab(lmax);
  double *ImYtab = alloc_Plm_Tab(lmax);

  /* get Theta,Phi from A,B */
  ThetaPhi_of_AB_CubSph(pat, AB[0],AB[1], &Theta,&Phi);

  /* make tables of Ylm at Theta,Phi */
  set_YlmTabs(lmax, Theta,Phi, ReYtab, ImYtab);

  /* get func val fv at Theta,Phi */
  fv = 0.;
  /* loop over positive m, here Ylmm=Y_l^{-m}, sm = (-1)^m */
  ijk=0;
  for(l=0; l<=lmax; l++)
    for(sm=1., m=0;  m<=l;  m++, sm=-sm)
    {
      double Rclm, Iclm, Re_Ylm, Im_Ylm;
      double Rclmm, Iclmm, Re_Ylmm, Im_Ylmm;
      /* get real and imag part of coeffs in co */
      Rclm = co[ijk++];
      Iclm = co[ijk++];

      /* get Ylm at Theta,Phi */
      Ylm_from_Tabs(lmax, ReYtab, ImYtab, l,m, &Re_Ylm,&Im_Ylm);

      /* There is a choice of sign here: define the inner product by
         (f,g) = int f^* g
         and define
         psi_mode = (Ylm, psi)  */
      /* fv = \sum_{l,m} c_l^m Y_l^m
         fv = \sum_{l,m} (  Re_c_lm Re_Ylm -   Im_c_lm Im_Ylm +
                          i Re_c_lm Im_Ylm + i Im_c_lm Re_Ylm   )  */
      /* fv = \sum_{l,m} c_l^m Y_l^m
            = \sum_l [ c_l^0 Y_l^0 +
                      \sum_{m=1}^l ( c_l^m Y_l^m + c_l^{-m} Y_l^{-m} ) ]
         Note: Y_l^{-m} = (-1)^m (Y_l^m)^*  <--always
               c_l^{-m} = (-1)^m (c_l^m)^*  <--if fv is real */
      Re_Ylmm =  sm*Re_Ylm;
      Im_Ylmm = -sm*Im_Ylm;
      Rclmm =  sm*Rclm; /* assuming fv is real */
      Iclmm = -sm*Iclm;

      if(m==0)
        fv += Rclm*Re_Ylm;
      else /* assuming fv is real */
        fv += Rclm*Re_Ylm - Iclm*Im_Ylm + Rclmm*Re_Ylmm - Iclmm*Im_Ylmm;
    }
  free(ImYtab);
  free(ReYtab);
  *sig = fv;
  return 0;
}

/* compute values of surface function derivs */
int FSurf_CubSph_sigma01_derivs(tPat *pat, int si, double AB[2], 
                                double dsig[2])
{
  tArray *Co = pat->CI->Fcoef[si];
  int N = Co->N;
  int nYs = N/2; /* number of Ylm's we use */
  int lmax = (sqrt(8*nYs + 1) - 3)/2;
  int ijk, l,m;
  double sm;
  double *co = Co->d;
  double ft, fp; /* ft = sin(Theta) dsigma/dTheta, fp = dsigma/dPhi */
  double fth;    /* fth= dsigma/dTheta */
  double Theta,Phi, dThetadA,dThetadB, dPhidA,dPhidB;
  double *ReYtab = alloc_Plm_Tab(lmax);
  double *ImYtab = alloc_Plm_Tab(lmax);
  double *csdth = calloc(nYs*2, sizeof(double));
  double *cdphi = calloc(nYs*2, sizeof(double));
  double A=AB[0], B=AB[1];

  /* regularize case where A=B=0 <==> Theta=0:
     Note for Theta=0 we cannot use sin(Theta) d/dTheta Ylm
     to find d/dTheta Ylm. Also dPhidA blows up!!!
     So for now just add epsilon to B. */
  if(A==0. && B==0.) B = 1e-10;

  /* get Theta,Phi and their derivs from A,B */
  ThetaPhi_dThetaPhidAB_of_AB_CubSph(pat, A,B, &Theta,&Phi,
                                     &dThetadA,&dThetadB, &dPhidA,&dPhidB);
  /* make tables of Ylm at Theta,Phi */
  set_YlmTabs(lmax, Theta,Phi, ReYtab, ImYtab);

  /* get coeffs of derivs */
  SphHarm_sin_theta_dtheta_forRealFunc(co, csdth, lmax);
  SphHarm_dphi_forRealFunc(co, cdphi, lmax);

  /* get func vals ft,fp at Theta,Phi */
  fp = ft = 0.;
  ijk=0;
  /* loop over positive m, here Ylmm=Y_l^{-m}, sm = (-1)^m */
  for(l=0; l<=lmax; l++)
    for(sm=1., m=0;  m<=l;  m++, sm=-sm)
    {
      double Re_Ylm, Im_Ylm, Re_Ylmm, Im_Ylmm;
      double Rcsdth, Icsdth, Rcsdthm, Icsdthm;
      double Rcdphi, Icdphi, Rcdphim, Icdphim;

      /* get real and imag part of coeffs in arrays */
      Rcsdth = csdth[ijk];
      Rcdphi = cdphi[ijk++];
      Icsdth = csdth[ijk];
      Icdphi = cdphi[ijk++];
/*
if(!isfinite(Rcsdth+Rcdphi+Icsdth+Icdphi))
{
printf("l=%d m=%d ijk=%d Rcsdth=%g Rcdphi=%g Icsdth=%g Icdphi=%g\n",
l,m,ijk, Rcsdth,Rcdphi,Icsdth,Icdphi);
errorexit("NAN!");
}
*/
      /* get Ylm at Theta,Phi */
      Ylm_from_Tabs(lmax, ReYtab, ImYtab, l,m, &Re_Ylm,&Im_Ylm);

      /* There is a choice of sign here: define the inner product by
         (f,g) = int f^* g
         and define
         psi_mode = (Ylm, psi)  */
      /* fv = \sum_{l,m} c_l^m Y_l^m
         fv = \sum_{l,m} (  Re_c_lm Re_Ylm -   Im_c_lm Im_Ylm +
                          i Re_c_lm Im_Ylm + i Im_c_lm Re_Ylm   )  */
      /* fv = \sum_{l,m} c_l^m Y_l^m
            = \sum_l [ c_l^0 Y_l^0 +
                      \sum_{m=1}^l ( c_l^m Y_l^m + c_l^{-m} Y_l^{-m} ) ]
         Note: Y_l^{-m} = (-1)^m (Y_l^m)^*  <--always
               c_l^{-m} = (-1)^m (c_l^m)^*  <--if fv is real */
      Re_Ylmm =  sm*Re_Ylm;
      Im_Ylmm = -sm*Im_Ylm;
      Rcsdthm =  sm*Rcsdth; /* assuming func is real */
      Icsdthm = -sm*Icsdth;
      Rcdphim =  sm*Rcdphi; /* assuming func is real */
      Icdphim = -sm*Icdphi;

      if(m==0)
      {
        ft += Rcsdth*Re_Ylm;
        fp += Rcdphi*Re_Ylm;
      }
      else /* assuming fv is real */
      {
        ft += Rcsdth*Re_Ylm - Icsdth*Im_Ylm + Rcsdthm*Re_Ylmm - Icsdthm*Im_Ylmm;
        fp += Rcdphi*Re_Ylm - Icdphi*Im_Ylm + Rcdphim*Re_Ylmm - Icdphim*Im_Ylmm;
      }
    }
  /* get derivs from ft,fp, and dThetadA,dThetadB, dPhidA,dPhidB */
  fth = ft/sin(Theta);
  dsig[0] = fth*dThetadA + fp*dPhidA;
  dsig[1] = fth*dThetadB + fp*dPhidB;
  free(cdphi);
  free(csdth);
  free(ImYtab);
  free(ReYtab);
/*
if(!isfinite(dsig[0]) || !isfinite(dsig[1]))
{
printf("fth=%g fp=%g dThetadA=%g dThetadB=%g dPhidA=%g dPhidB=%g\n",
fth,fp, dThetadA,dThetadB, dPhidA,dPhidB);
errorexit("NAN!");
}
*/
  return 0;
}


/* set var pat->CI->iSurf and its derivs from FSurf_CubSph_sigma01_func */
int FSurf_CubSph_set_sigma01vars_from_sigma01_func(tNode *node, int si)
{
  tPat *pat = node->pat;
  tMesh *mesh = pat->mesh;
  tCoordInfo *CI = pat->CI;
  int iX = Ind("X");
  double *Yp = Vard(node, iX+1);
  double *Zp = Vard(node, iX+2);
  int isigma    = CI->iSurf[si];
  int isigma_dA = CI->idSurfdX[si][1];
  int isigma_dB = CI->idSurfdX[si][2];
  double *sigma    = Vard(node, isigma);
  double *sigma_dA = Vard(node, isigma_dA);
  double *sigma_dB = Vard(node, isigma_dB);
  double AB[2], dsig[2], sig[1];
  int *n  = node->n;
  int *ns = Varn(node, isigma);
  int i,j,k, ijk, sjk;

  forplane0(i,j,k, n, (n[0]-1)*(si==1))
  {
    /* get A,B at point ijk */
    ijk = Ind_n(i,j,k, n);
    AB[0] = Yp[ijk];
    AB[1] = Zp[ijk];

    /* set sigma01 var and derivs */
    sjk = Ind_n(si,j,k, ns);
    CI->FSurf[si](pat,si, AB, sig);
    sigma[sjk] = sig[0];
    CI->dFSurfdC[si](pat,si, AB, dsig);
    sigma_dA[sjk] = dsig[0];
    sigma_dB[sjk] = dsig[1];
  }

  return 0;
}


/* set CI->FSurf func pointer and alloc room for coeffs in CI->Fcoef */
int FSurf_CubSph_init6pats(tMesh *mesh, int pi_dom0)
{
  tPat *pat = mesh->pat[pi_dom0];
  int type = pat->CI->type;
  int dom  = pat->CI->dom;
  int i, si, si0, si1;
  int CubedSphere_sigma01_lmax = Par("CubedSphere_sigma01_lmax");

  if(dom!=0) return -1; /* do nothing if this is not dom0 */

  /* set lmax we use */
  lmax = Geti(CubedSphere_sigma01_lmax);
  if(lmax<1) errorexit("lmax<1 is suspicious!");
  /* We need (lmax*(lmax+1))/2 + lmax+1  complex numbers at each point A,B
     to store the table.
     So when is (lmax*(lmax+1))/2 + lmax+1 = n0?
     set L = lmax ==> L^2/2 + 3L/2 + 1 = n0  <==> L^2 + 3 L + 2 - 2*n0 = 0
     so: 2L = -3 +- sqrt(9 - 4*(2 - 2*n0)) = -3 +- sqrt(8*n0 + 1) 
     L = (sqrt(8*n0 + 1) - 3)/2  */

  /* figure out range of si */
  switch(type)
  {
  case outerCubedSphere:
    si0 = si1 = 1;
    break;
  case innerCubedSphere:
    si0 = si1 = 0;
    break;
  case CubedShell:
    si0 = 0;
    si1 = 1;
    break;
  default:
    si0 = +2; /* do not loop over si */
    si1 = -1;
  }

  /* loop if si0<=si1 */
  for(si=si0; si<=si1; si++)
  {
    int nYs = (lmax*(lmax+1))/2 + lmax+1; /* number of Ylm's we use */
    int nc[] = { nYs*2, 1, 1 };           /* number of coeffs we need */

    /* alloc memory for coeffs in one of the six patches */
    if(pat->CI->Fcoef[si])
      pat->CI->Fcoef[si] = redimension_array(pat->CI->Fcoef[si], nc);
    else
      pat->CI->Fcoef[si] = alloc_array(nc);

    /* loop over 6 patches */
    for(i=0; i<6; i++)
    {
      tPat *pati = mesh->pat[pi_dom0 + i];
      tCoordInfo *CI = pati->CI;

      /* set surface functions */
      CI->FSurf[si]    = FSurf_CubSph_sigma01_func;
      CI->dFSurfdC[si] = FSurf_CubSph_sigma01_derivs;
    }
  }

  return 0;
}
