/* limiter.c */
/* Wolfgang Tichy April, 2019 */


#include "nmesh.h"
#include "limiter.h"


/* the limiters in here use the struct tINDIC within tDat to exchange data
   between neighbors. The exchange is done in evolve, which calls functions
   from main/amr/indicators.c to do the exchange. */


/* set data for MRS limiter */
int limdata_MRS(tNode *node, tVarList *vl)
{
  int np;
  int vli;
  tDat *dat;
  int nvals = 2; /* we need 2 data values per var */

  if(!node || !vl) return nvals;

  dat = node->dat;
  if(!dat) return nvals;

  np = node->np;

  forvl(vl, vli)
  {
    int iq = Vind(vl, vli);
    double *q = Vard(node, iq);
    int im;

    /* find min and max in q in this node, and write it into indicator */
    dat->ic[iq]->myindc->d[0] = min_in_1d_array(q, np, &im);
    dat->ic[iq]->myindc->d[1] = max_in_1d_array(q, np, &im);
  }
  return nvals;
}

/* funcs in MRS limiter */
double phiy(double y)
{
  if(y<0.) return 0.;      // is this correct??? Read papers!
  return min2(y/1.1, 1.);
}

double theta_Mm(double Mi, double qbar, double qMi)
{
  double r;
  double num = Mi - qbar;
  double den = qMi - qbar;
  if(den != 0.) r = num/den;
  else          return 1.;  // what should this be 0 or 1???
  return phiy(r);
}

/* limiter: limit u using data in dat->ic */
int limiter_MRS(tNode *node, tVarList *vl)
{
  double *bb = node->bbox;
  tDat *dat;
  int vli, f, ni, ijk;
  double qbar, qMi, qmi, Mi, mi, theta_Mi, theta_mi, theta_i;
  double alpha, h, alpha_h;

  dat = node->dat;
  if(!dat) return 0;

  /* set alpha_h */
  alpha = 5.0;
  h = max3(bb[1]-bb[0], bb[3]-bb[2], bb[5]-bb[4]);
  alpha_h = alpha * pow(h, 1.5);

  forvl(vl, vli)
  {
    int iq = Vind(vl, vli);
    double *q = Vard(node, iq);

    /* find node average qbar */
    qbar = var_nodeaverage(node, iq);
//printf("qbar=%g theta_i=%g\n", qbar, theta_i);
//exit(88);

    /* get min, max on node */
    qmi = dat->ic[iq]->myindc->d[0];
    qMi = dat->ic[iq]->myindc->d[1];

    /* find min and max of q in neighbors */
    Mi = -1e300;
    mi = 1e300;
    for(f=0; f<6; f++)
      for(ni=0; ni<node->nfnb[f]; ni++)
      {
        int ma = dat->ic[iq]->nbindc[f][ni]->d[0];
        int Ma = dat->ic[iq]->nbindc[f][ni]->d[1];
        if(ma < mi) mi = ma;
        if(Ma > Mi) Mi = Ma;
      }

    mi = min2(qbar - alpha_h, mi);
    Mi = max2(qbar + alpha_h, Mi);

    /* set thetas */
    theta_Mi = theta_Mm(Mi, qbar, qMi);
    theta_mi = theta_Mm(mi, qbar, qmi);
    theta_i = min3(1., theta_mi, theta_Mi);

if(!isfinite(qbar) || !isfinite(theta_i))
{
printf("qbar=%g theta_i=%g\n", qbar, theta_i);
abort();
exit(8);
}
    /* now limit q */
    forpoints(node, ijk)
      q[ijk] = qbar + theta_i*(q[ijk] - qbar);
  }

  return 0;
}
