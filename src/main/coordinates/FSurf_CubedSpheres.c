/* FSurf_CubedSpheres.c */
/* Wolfgang Tichy, Oct 2018 */
/* functions to to get sigma on surface of cubed sphere */

#include "nmesh.h"
#include "coordinates.h"


/* return coeffs for face0 of 1st patch */
tArray *FSurf_CubSph_sigma01_Fcoef(tPat *pat, int si)
{
  tMesh *mesh = pat->mesh;
  int p0 = pat->pg0;
  int np = pat->npg;
  int pi = (si ? p0+np : p0); /* if si!=0 go to next pat group */

  if(pi >= mesh->npats) return NULL; /* there are no coeffs */

  return mesh->pat[pi]->CI->Fcoef[0];
}

/* return value of surface function sigma01 */
int FSurf_CubSph_sigma01(tPat *pat, int si, double AB[2], double *sig)
{
  tArray *Co = FSurf_CubSph_sigma01_Fcoef(pat, si);

  if(Co) /* if we have coeffs we now compute sig from them */
  {

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

    errorexit("this function needs to be tested!!!");

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
  }
  else /* if there are no coeffs we return a constant sig */
  {
    *sig = pat->CI->s[si];
  }
  return 0;
}

/* compute values of surface function derivs */
int FSurf_CubSph_dsigma01(tPat *pat, int si, double AB[2], double dsig[2])
{
  tArray *Co = FSurf_CubSph_sigma01_Fcoef(pat, si);

  if(Co) /* if we have coeffs we now compute dsig from them */
  {
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

    errorexit("this function needs to be tested!!!");

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
if(!finit(Rcsdth+Rcdphi+Icsdth+Icdphi))
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
if(!finit(dsig[0]) || !finit(dsig[1]))
{
printf("fth=%g fp=%g dThetadA=%g dThetadB=%g dPhidA=%g dPhidB=%g\n",
fth,fp, dThetadA,dThetadB, dPhidA,dPhidB);
errorexit("NAN!");
}
*/
  }
  else /* if there are no coeffs sig is constant */
  {
    dsig[0] = dsig[1] = 0.;
  }
  return 0;
}


/* set var pat->CI->iSurf and its derivs from FSurf_CubSph_sigma01 */
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
  int i,j,k, ijk, pl, sjk;

  errorexit("this function may not be needed and it needs to be tested!!!");
  errorexit("I think coordinates_init_node sets the sigma01 vars already...");

  pl = (n[0]-1)*(si==1);
  forplane0(i,j,k, n, pl)
  {
    /* get A,B at point ijk */
    ijk = Ind_n(i,j,k, n);
    AB[0] = Yp[ijk];
    AB[1] = Zp[ijk];

    /* set sigma01 var and derivs */
    sjk = Ind_n(0,j,k, ns);
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
  tPat *pat0 = mesh->pat[pi_dom0];
  int npg = pat0->npg; /* pick patch group */
  int pg0 = pat0->pg0;
  tPat *pat = mesh->pat[pg0];
  int type = pat->CI->type;
  //int dom  = pat->CI->dom;
  int i, si, si0, si1;
  int CubedSphere_sigma01_lmax = Par("CubedSphere_sigma01_lmax");
  int lmax;

  errorexit("this function needs to be tested!!!");

  if(npg==0) errorexit("we need a patgroup, i.e. npg>0");
  if(pi_dom0!=pg0) return -1; /* do nothing if this is not dom0 */

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
    //FIXME: alloc mem only for pat0 in group, and only for si=0
    if(pat->CI->Fcoef[si])
      redimension_array(pat->CI->Fcoef[si], nc);
    else
      pat->CI->Fcoef[si] = alloc_array(nc);

    /* loop over 6 patches */
    for(i=0; i<npg; i++)
    {
      tPat *pati = mesh->pat[pg0 + i];
      tCoordInfo *CI = pati->CI;

      /* set surface functions */
      CI->FSurf[si]    = FSurf_CubSph_sigma01;
      CI->dFSurfdC[si] = FSurf_CubSph_dsigma01;
    }
  }

  return 0;
}


/***********************************************************************/
/* funcs to compute Ylm coeffs */
/***********************************************************************/

/* put the Ylm into two C-arrays Re_Ylmp, Im_Ylmp */
/* We can use array_GLquadrature2X to compute surface integrals.
   For each mode coeff we need, we could add one radial point.
   E.g. make vars for Re and Im where at
   i=0     we put rY_0^{0} and iY_0^{0}
   i=1     we put rY_1^{0} and iY_1^{0}
   i=2     we put rY_1^{1} and iY_1^{1}
   i=3     we put rY_2^{0} and iY_2^{0}
   i=4     we put rY_2^{1} and iY_2^{1}
   i=5     we put rY_2^{2} and iY_2^{2}
   I.e. we use only positive m, because Y_l^{-m} = (-1)^m (Y_l^m)^* .
   Then I can use array_GLquadrature2X over these vars to compute the all
   coeffs. */
/* NOTE: Re_Ylmp,Im_Ylmp have size N0*n1*n2, where N0=n0*S0, with S0>=1 */
int FSurf_CubSph_set_ReIm_Ylm(tNode *node, int S0,
                              double *Re_Ylmp, double *Im_Ylmp, int lmax)
{
  tMesh *mesh = Elm_mesh(node);
  tPat *pat = node->pat;
  int *nn = node->n;
  int n0=nn[0];
  int n1=nn[1];
  int n2=nn[2];
  int Ng=n0*n1*n2;
  double *Yp = Vard(node, Ind("Y"));
  double *Zp = Vard(node, Ind("Z"));
  double *ReYtab = alloc_Plm_Tab(lmax);
  double *ImYtab = alloc_Plm_Tab(lmax);
  int nYs = (lmax*(lmax+1))/2 + lmax+1; /* number of Ylm's we use */
  int N0 = n0*S0; /* i-range of Re_Ylmp,Im_Ylmp arrays */
  int l,m, i,j,k, ijk, Ijk;

  if(nYs>N0) errorexit("size of Re_Ylmp,Im_Ylmp arrays is too small");

  /* loop over l and m and set Ylm over surface. Put each Ylm at a
     different radial coord for each l,m */
  //printf("setting Ylm in box%d\n", box->b);
  i=0;
  for(k=0; k<n2; k++)
  for(j=0; j<n1; j++)
  {
    double A,B, theta,phi, Re_Ylm,Im_Ylm;

    /* point at i=0 */
    i=0;
    ijk=Ind_n(i,j,k, nn);

    /* get A,B at point ijk */
    A = Yp[ijk];
    B = Zp[ijk];

    /* get theta,phi from A,B */
    ThetaPhi_of_AB_CubSph(pat, A,B, &theta,&phi);

    /* make tables of Ylm at Theta,Phi */
    set_YlmTabs(lmax, theta,phi, ReYtab, ImYtab);

    /* set all Ylm with positive m, since Y_l^{-m} = (-1)^m (Y_l^m)^* */
    for(l=0; l<=lmax; l++)
    for(m=0; m<=l; m++)
    {
      /* get Ylm at theta,phi */
      Ylm_from_Tabs(lmax, ReYtab, ImYtab, l,m, &Re_Ylm,&Im_Ylm);

      /* set spherical harmonic Ylm at point ijk.
         NOTE: Re_Ylmp and Im_Ylmp may not be on the grid and thus have a
         different range for the index i */
      ijk=Ind_n(i%n0,j,k, nn);
      Ijk=ijk + Ng*(i/n0);
      Re_Ylmp[Ijk] = Re_Ylm;
      Im_Ylmp[Ijk] = Im_Ylm;
      i++;
    }
  }

  free(ImYtab);
  free(ReYtab);
  /* return total number of coeffs up to l=lmax */
  return i;
}


/* Compute integrals of (Ylm^* var) that have real part at varindex Re_vind and
   imag. part at Im_vind. Do integrals over surface with index i=s in X-dir.
   Add integrals to array aInteg.
   FSurf_CubSph_add_Ylm_integrals will be called inside a loop over nodes.
   Make sure aInteg is zeroed before we loop over these nodes.
   If var has zero imag. part set Im_vind=-1. */
   //NOTE: we probably only need s=0
int FSurf_CubSph_add_Ylm_integrals(tNode *node, int s, int Re_vind, int Im_vind,
                                   int lmax, tArray *aInteg)
{
  tMesh *mesh = Elm_mesh(node);
  tPat *pat = node->pat;
  int *nn = node->n;
  int l,m, i,j,k, ijk, Ijk, seg;
  int n0=nn[0];
  int n1=nn[1];
  int n2=nn[2];
  int Ng=n0*n1*n2;
  int nYs = (lmax*(lmax+1))/2 + lmax+1; /* number of Ylm's we use */
  int N0;                               /* i-range of Re_Ylmp,Im_Ylmp arrays */
  int S0;                               /* num. of segments: S0 = N0/n0 */
  //int offset; /* offset used to write into var Integ_ind */
  double *Re_varp = Vard(node, Re_vind);
  double *Im_varp;
  double *Re_Ylmp;
  double *Im_Ylmp;
  double *Re_Integp;
  double *Im_Integp;
  //double *Integ = Vard(node, Integ_ind);
  double *Integ = Arrd(aInteg);
  double *Yp = Vard(node, Ind("Y"));
  double *Zp = Vard(node, Ind("Z"));
  tArray *Re_Integ; // array needed for array_GLquadrature2X
  tArray *Im_Integ; // array needed for array_GLquadrature2X

  //printf("lmax=%d\n", lmax);

  //printf("VarName(Re_vind)=%s Re_vind=%d Re_varp[s]=%g\n",
  //VarName(Re_vind),Re_vind, Re_varp[s]);
  //quick_Vars_output(box, VarName(Re_vind), 7,7);

  /* do we have imag. part in our var? */
  if(Im_vind>0) Im_varp = Vard(node, Im_vind);
  else          Im_varp = NULL;

  /* make room for all the Ylm's */
  S0 = nYs/n0;     /* number of segments if nYs is divisible by n0 */
  if(nYs%n0) S0++; /* if there was a remainder, increase S0 */
  //PRF;printf(": b=%d lmax=%d nYs=%d n0=%d S0=%d\n", box->b, lmax, nYs, n0, S0);
  N0 = S0*n0;
  Re_Ylmp = calloc(N0*n1*n2, sizeof(double));
  Im_Ylmp = calloc(N0*n1*n2, sizeof(double));
  Re_Integp = calloc(N0*n1*n2, sizeof(double));
  Im_Integp = calloc(N0*n1*n2, sizeof(double));
  if(Re_Ylmp==NULL || Im_Ylmp==NULL || Re_Integp==NULL || Im_Integp==NULL)
    errorexit("out of memory for Re_Ylmp, Im_Ylmp, ...");

  /* precompute the Ylm */
  FSurf_CubSph_set_ReIm_Ylm(node, S0, Re_Ylmp, Im_Ylmp, lmax);

  /* set integrands */
  for(k=0; k<n2; k++)
  for(j=0; j<n1; j++)
  {
    i=0;
    for(l=0; l<=lmax; l++)
    for(m=0; m<=l; m++, i++) /* here we set only integrands for m>=0 */
    {
      double R,I, RYlm,IYlm, A,B;
      double Theta,Phi, dThetadA,dThetadB, dPhidA,dPhidB, Jac, fac;

      /* get A,B and Re, Im part of data at surface where i=s */
      ijk=Ind_n(s,j,k, nn);
      A = Yp[ijk];
      B = Zp[ijk];
      R = Re_varp[ijk];
      if(Im_vind<=0) I = 0.; /* imag. part is zero */
      else           I = Im_varp[ijk];

      /* We need to compute \int d\phi d\theta \sin(theta) (Y_l^m)^* var .
         Later we actually compute  \int dA dB (Integ).
         Now \int d\phi d\theta \sin(theta) = \int dA dB Jac \sin(theta)
         So we need to multiply by the Jacobian  */
      /* get Theta, Phi and their derivs */
      ThetaPhi_dThetaPhidAB_of_AB_CubSph(pat, A,B, &Theta,&Phi,
                                         &dThetadA,&dThetadB, &dPhidA,&dPhidB);
      Jac = fabs(dThetadA*dPhidB - dThetadB*dPhidA); /* Jacobian */
      fac = Jac * sin(Theta);
      R = R * fac;
      I = I * fac;

      /* get spherical harmonic Ylm */
      ijk=Ind_n(i%n0,j,k, nn);
      Ijk=ijk + Ng*(i/n0);
      RYlm = Re_Ylmp[Ijk];
      IYlm = Im_Ylmp[Ijk];

      /* There is a choice of sign here: define the inner product by
         (f,g) = int f^* g
         and define
         psi_Integ = (Y, psi)  */
      Re_Integp[Ijk] = RYlm * R + IYlm * I;
      Im_Integp[Ijk] = RYlm * I - IYlm * R;
      //printf("b%ds%d Jac=%g R=%g RYlm=%g IYlm=%g @ %g %g\n",
      //box->b,s, Jac, R, RYlm, IYlm, A,B);
      //quick_Array_output(box, Re_Integp, "Re_Integp", 8,8);
    }
  }

  /* setup arrays for integration */
  Re_Integ = alloc_empty_array_with_segs(nn, 0,1);
  Im_Integ = alloc_empty_array_with_segs(nn, 0,1);

  /* integrate over surfaces */
  /* If we have more than one segment (S0>1) we need array_GLquadrature2X
     calls for each segment! */
  for(seg=0; seg<S0; seg++)
  {
    int os = Ng*seg;

    //spec_2dIntegral(box, 1, Re_Integp+os, Re_Integp+os);
    //spec_2dIntegral(box, 1, Im_Integp+os, Im_Integp+os);

    /* point array to Re_Integp+os data and set correct dims */
    point_array_d_to_data(Re_Integ, Re_Integp+os, 1);
    point_array_d_to_data(Im_Integ, Im_Integp+os, 1);
    /* set array-sizes so that redimension_array will not realloc */
    Re_Integ->size = Im_Integ->size = n0*n1*n2 * sizeof(double);
    redimension_array(Re_Integ, nn);
    redimension_array(Im_Integ, nn);

    /* do the integrals */
    array_GLquadrature2X(node, 0, Re_Integ, Re_Integ);
    array_GLquadrature2X(node, 0, Im_Integ, Im_Integ);
  }
  //quick_Array_output(box, Re_Integp, "Re_Integp", 9,9);
  free_array(Re_Integ);
  free_array(Im_Integ);

  ///* Put Integs into var with index Integ_ind */
  //if(nYs>=Ng/4) errorexit("decrease lmax!");
  //offset = ((Ng/2)*s)/(n0-1);  /* offset for Integs in Integ */
  //ijk = offset;
  /* Add Integs to Integ */
  if(nYs>=ArrN(aInteg)) errorexit("decrease lmax!");
  //offset = ((Ng/2)*s)/(n0-1);  /* offset for Integs in Integ */
  ijk = 0;
  i = 0;
  for(l=0; l<=lmax; l++)
  for(m=0; m<=l; m++)
  {
    /* set Re and Im part of Integ */
    Ijk = (i%n0) + Ng*(i/n0);        //FIXME: is this correct ????

    //Integ[ijk++] = Re_Integp[Ijk];
    //Integ[ijk++] = Im_Integp[Ijk];
    GEN_Pragma(omp atomic update)
    { Integ[ijk] += Re_Integp[Ijk]; }
    ijk++;

    GEN_Pragma(omp atomic update)
    { Integ[ijk] += Im_Integp[Ijk]; }
    ijk++;

    i++;
  }

  free(Im_Integp);
  free(Re_Integp);
  free(Im_Ylmp);
  free(Re_Ylmp);
  return 0;
}

/* take integrals over all nodes at cub.sph. patch faces and accumulate
   them in aco, such that they become the coeffs in the Ylm expansion */
int FSurf_CubSph_set_Ylm_coefArray(tMesh *mesh, int s, int pi_dom0,
                                   int Re_vind, int Im_vind, int lmax,
                                   tArray *aco)
{
  tPat *pat0 = mesh->pat[pi_dom0];
  int npg = pat0->npg; /* pick patch group */
  int pg0 = pat0->pg0;
  int f = s ? 1 : 0; //pick face

  /* zero aco */
  memset(Arrd_(aco), 0, ArrN(aco));

  /* add integral pieces from all nodes on face0 on this rank to aco */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    tPat *pat = node->pat;
    int p = pat->p;

    if(p<pg0 || p>=pg0+npg) continue;
    if(!Node_patface(node, f)) continue;

    FSurf_CubSph_add_Ylm_integrals(node, s, Re_vind,Im_vind, lmax, aco);
  }

  /* use AllReduce to add up aco from different ranks */
  MCK( nMPI_Allreduce(nMPI_IN_PLACE, Arrd(aco), ArrN(aco),
                      nMPI_DOUBLE, nMPI_SUM) );
  return 0;
}


/* Set CI->Fcoef[0] in 1st patch of cub.sph. group. We get these coeffs
   by integrating CubedSphere_sigma0_def against Ylm's. */
int FSurf_CubSph_set_CI_Fcoef0(tMesh *mesh, int pi_dom0)
{
  tPat *pat0 = mesh->pat[pi_dom0];
  int npg = pat0->npg; /* pick patch group */
  int pg0 = pat0->pg0;
  tPat *pat = mesh->pat[pg0];
  int type = pat->CI->type;
  //int dom  = pat->CI->dom;
  int innerSphere;
  int isigma0 = Ind("CubedSphere_sigma0_def");
  int CubedSphere_sigma01_lmax = Par("CubedSphere_sigma01_lmax");
  int lmax;

  errorexit("this function needs to be tested!!!");

  if(npg==0) errorexit("we need a patgroup, i.e. npg>0");
  if(pi_dom0!=pg0) return -1; /* do nothing if this is not dom0 */

  /* set lmax we use */
  lmax = Geti(CubedSphere_sigma01_lmax);
  if(lmax<1) errorexit("lmax<1 is suspicious!");
  /* We need (lmax*(lmax+1))/2 + lmax+1  complex numbers at each point A,B
     to store the table.
     So when is (lmax*(lmax+1))/2 + lmax+1 = n0?
     set L = lmax ==> L^2/2 + 3L/2 + 1 = n0  <==> L^2 + 3 L + 2 - 2*n0 = 0
     so: 2L = -3 +- sqrt(9 - 4*(2 - 2*n0)) = -3 +- sqrt(8*n0 + 1)
     L = (sqrt(8*n0 + 1) - 3)/2  */

  /* figure out if face0 is spherical */
  switch(type)
  {
  case innerCubedSphere:
  case CubedShell:
    innerSphere = 1;
    break;
  default:
    innerSphere = 0;
  }

  if(innerSphere)
    FSurf_CubSph_set_Ylm_coefArray(mesh, 0, pg0, isigma0,-1, lmax,
                                   pat->CI->Fcoef[0]);
  return 0;
}
