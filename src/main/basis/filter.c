/* filter.c */
/* Wolfgang Tichy, 8/2019 */

#include "nmesh.h"
#include "basis.h"


#define PR 0


/***********************************************************************/
/* functions to filter */
/***********************************************************************/

/* apply exponential filter as in https://arxiv.org/abs/1804.02003 */
void expfilter_coeff_array(tArray *ca, double alp[3], double s[3])
{
  int *n = ca->n;
  double N0 = n[0] - 1;
  double N1 = n[1] - 1;
  double N2 = n[2] - 1;
  double *c = Arrd_(ca);
  double alp0 = alp[0];
  double alp1 = alp[1];
  double alp2 = alp[2];
  double s0 = s[0];
  double s1 = s[1];
  double s2 = s[2];
  double f0, f1, f2;
  int i,j,k, ijk;

  /* make sure we don't divide by zero */
  if(N0==0.) N0 = 1.;
  if(N1==0.) N1 = 1.;
  if(N2==0.) N2 = 1.;

  /* loop over array and let c -> c * f0*f1*f2 */
  for(k = 0; k < n[2]; k++)
  {
    f2 = exp( -alp2 * pow(k/N2, s2) );
    for(j = 0; j < n[1]; j++)
    {
      f1 = exp( -alp1 * pow(j/N1, s1) );
      for(i = 0; i < n[0]; i++)
      {
        f0 = exp( -alp0 * pow(i/N0, s0) );
        ijk = Ind_n(i,j,k, n);
        c[ijk] = c[ijk] * f0 * f1 * f2;
      }
    }
  }
}

/* get coeffs ca of array ua, filter and then reset ua */
void expfilter_array(tNode *node, tArray *ua, double alp[3], double s[3])
{
  DECL_STACK_ARRAY(ca, ua->n);

  basis_array_analysis3(node, ua, ca);
  expfilter_coeff_array(ca, alp, s);
  basis_array_synthesis3(node, ua, ca);
}

/* filter var with index ui */
int expfilter_var(tNode *node, int ui, double alp[3], double s[3])
{
  tArray *ua;
  tDat *dat = node->dat;
  if(dat)
  {
    ua = dat->v[ui];
    expfilter_array(node, ua, alp, s);
    return 1;
  }
  else
    return 0;
}


/***********************************************************************/
/* functions to check coeff falloff */
/***********************************************************************/

/* is ca falling off */
int expfalloff_coeff_array(tArray *ca, double alp[3], double s[3])
{
  int *n = ca->n;
  double N0 = n[0] - 1;
  double N1 = n[1] - 1;
  double N2 = n[2] - 1;
  double *c = Arrd_(ca);
  double absc;
  double alp0 = alp[0];
  double alp1 = alp[1];
  double alp2 = alp[2];
  double s0 = s[0];
  double s1 = s[1];
  double s2 = s[2];
  double f0, f1, f2;
  double norm;
  int i,j,k, ijk;

  /* make sure we don't divide by zero */
  if(N0==0.) N0 = 1.;
  if(N1==0.) N1 = 1.;
  if(N2==0.) N2 = 1.;

  /* L1 norm of coeffs */
  norm = Lp_norm_array(ca, 1.);

  /* loop over array and check if coeffs are falling off exponentially */
  for(k = 0; k < n[2]; k++)
  {
    f2 = exp( -alp2 * pow(k/N2, s2) );
    for(j = 0; j < n[1]; j++)
    {
      f1 = exp( -alp1 * pow(j/N1, s1) );
      for(i = 0; i < n[0]; i++)
      {
        f0 = exp( -alp0 * pow(i/N0, s0) );
        ijk = Ind_n(i,j,k, n);
        absc = fabs(c[ijk]);

        /* if absc is too big we do not have proper exp. falloff */
        if( (absc > norm * f0) || (absc > norm * f1) || (absc > norm * f2) )
          return 0;
      }
    }
  }
  return 1;
}

/* get coeffs ca of array ua, and check for falloff */
int expfalloff_array(tNode *node, tArray *ua, double alp[3], double s[3])
{
  DECL_STACK_ARRAY(ca, ua->n);
  basis_array_analysis3(node, ua, ca);
  return expfalloff_coeff_array(ca, alp, s);
}

/* check var for exp falloff */
int expfalloff_var(tNode *node, int ui, double alp[3], double s[3])
{
  tArray *ua;
  tDat *dat = node->dat;
  if(dat)
  {
    ua = dat->v[ui];
    return expfalloff_array(node, ua, alp, s);
  }
  return 1; /* return 1 if we do not have var ui on this proc */
}
