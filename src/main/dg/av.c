/* av.c */
/* Wolfgang Tichy, June 2026 */


#include "nmesh.h"
#include "dg.h"


/***********************************************************************/
/* funcs needed for artificial viscosity in nmesh */
/***********************************************************************/


/*************************************************************************/
/* Funcs to set mu */
/*************************************************************************/

/* get coeff falloff tau for var u from Persson ratio */
double av_tau_from_Persson(tElm *elm, int iu, double filter_alp,
                           double filter_s, double filter_dn, double f_unfilt)
{
  //tMesh *mesh = elm->pat->mesh;
  int n_unfilt[3], n_max;
  double se, ser, tau;

  /* reduce n to take into account filter for GRHD_D on dg grid */
  unfiltered_range_of_expfilter1(elm->n, filter_alp,filter_s,filter_dn,
                                 f_unfilt, n_unfilt);

  /* use mod. Persson indicator from evolve */
  se = evolve_Pmod_indicator_ncoeffs(elm, iu, -2., n_unfilt);

  /* assume coeffs fall off as c_n = k_n n^{-tau}
     se = log10(c_m^2 / ( \sum_{n=0}^{m} c_n^2 )),  where m = n_max
     Since c_n = k_n (n)^{-tau}
     se ~ log10(k_m^2 (m)^{-2tau} / const) = log10(k_m^2/const) - 2tau log10(m)
        ~ -2 tau log10(m) */
  /* scale out log10(m) in se */
  n_max = max3(n_unfilt[0], n_unfilt[1], n_unfilt[2]);
  if(n_max>1) ser = se / log10(n_max);
  else        ser = se;

  /* estimate coeff falloff */
  tau = -0.5 * ser;

  if(0)
  {
    printeploc(elm->eploc);
    printf(" f_unfilt=%g n_max=%d ser=%g tau=%g\n", f_unfilt, n_max, ser, tau);
  }
  return tau;
}

/* compute mu from tau */
double av_mu_from_tau(tElm *elm, double tau,
                      double nL, double nH, double cmax)
{
  //tMesh *mesh = elm->pat->mesh;
  /* In Atteneder's thesis, he writes nL=1, nH=3 comes from
     [186] J. Yu and J. S. Hesthaven, “A study of several artificial
      viscosity models within the discontinuous galerkin framework”,
      Communications in Computational Physics 27, 1309 (2020),
      doi:10.4208/cicp.OA-2019-0118.  */
  double h = space_diagonal0_length(elm)/sqrt(3.);
  double N = max3(elm->n[0], elm->n[1], elm->n[2]) - 1;
  double mu_max;

  if(N<1.) N = 1.;

  mu_max = cmax * h/N;

  //TEST: use max mu for testing
  //return mu_max;

  if(tau < nL) return mu_max;
  if(tau > nH) return 0.;
  return mu_max * ( 1. - (tau - nL)/(nH - nL) );
}

/* set mu var in one elm */
void av_mu_elm(tElm *elm, tVarList *vlu, double cmax, int imu,
               double filter_alp,  double filter_s , double filter_dn,
               double f_unfilt)
{
  int vli, ijk;
  double *av_mu = Vard(elm, imu);
  double nL   = 1.;
  double nH   = 3.;
  double mu;

  {
    double taumin = DBL_MAX;

    /* find min tau for vars in vlu */
    forvl(vlu, vli)
    {
      int iu = Vind(vlu, vli);
      double tau = av_tau_from_Persson(elm, iu, filter_alp, filter_s,
                                       filter_dn, f_unfilt);
      if(tau < taumin) taumin = tau;
    }

    /* set mu */
    mu = av_mu_from_tau(elm, taumin, nL,nH, cmax);
  }

  /* make var av_mu const for now */
  forpoints(elm,ijk) av_mu[ijk]  = mu;
}


/*************************************************************************/
/* C^{\infty}_0 funcs from https://arxiv.org/abs/1810.02152 that can
   give mu a particular shape */
/*************************************************************************/

/* Gegenbauer viscosity */
double av_GegenbauerViscosity(double xb, double lam)
{
  return pow((1 - xb*xb), lam);
}

/* super Gaussian viscosity */
double av_SuperGaussianViscosity(double xb, double lam)
{
  double alp = 36; // = -ln(eps) with eps representing machine precision
  return exp(-alp * pow(xb, 2.*lam));
}

/* Gevrey viscosity */
double av_GevreyViscosity(double xb, double lam)
{
  double x2 = xb*xb;
  if(x2 < 1.) return exp(x2 / (lam*(x2 - 1.)));
  else        return 0.;
}
