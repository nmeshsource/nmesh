/* limiter.c */
/* Wolfgang Tichy April, 2019 */


#include "nmesh.h"
#include "limiter.h"


/* the limiters in here use the struct tINDIC within tDat to exchange data
   between neighbors. The exchange is done in evolve, which calls functions
   from main/amr/indicators.c to do the exchange. */


/* set data for MRS limiter. This one uses the cons vars in vl */
int limdata_MRS(tNode *node, tVarList *vl)
{
  int np;
  int vli;
  tDat *dat;
  int nvals = 3; /* we need 3 data values per var: min, max, node-average */

  if(!node || !vl) return nvals;

  dat = node->dat;
  if(!dat) return nvals;

  np = node->np;

  forvl(vl, vli)
  {
    int iq = Vind(vl, vli);
    double *q = Vard(node, iq);
    int im;

    /* also save node average */
    dat->ic[iq]->myindc->d[0] = var_nodeaverage(node, iq);

    /* find min and max in q in this node, and write it into indicator */
    dat->ic[iq]->myindc->d[1] = min_in_1d_array(q, np, &im);
    dat->ic[iq]->myindc->d[2] = max_in_1d_array(q, np, &im);
  }
  return nvals;
}

/* funcs in MRS limiter */
double MRS_phiy(double y)
{
  if(y<0.) return 0.;      // is this correct??? Read papers!
  return min2(y/1.1, 1.);
}

double MRS_theta_Mml(double Mi, double wbar, double wMi)
{
  double r;
  double num = Mi - wbar;
  double den = wMi - wbar;
  if(den != 0.) r = num/den;
  else          return 1.;  // what should this be 0 or 1???
  return MRS_phiy(r);
}

/* limiter: limit u using data in dat->ic */
int limiter_MRS(tNode *node, tVarList *vl)
{
  static int firstcall = 1;
  static int limiter_alpha;
  tMesh *mesh = node->pat->mesh;
  double *bb = node->bbox;
  tDat *dat;
  int vli, f, ni, ijk;
  double alpha, h, alpha_h, theta_i;
  double theta_Mi, theta_mi;

  dat = node->dat;
  if(!dat) return 0;

  if(firstcall)
  {
    limiter_alpha = Par("limiter_alpha");
    firstcall = 0;
  }

  /* set alpha_h from alpha (smaller alpha makes MRS more agressive) */
  alpha = Getd(limiter_alpha); //0.1; //5.0;
  h = max3(bb[1]-bb[0], bb[3]-bb[2], bb[5]-bb[4]);
  alpha_h = alpha * pow(h, 1.5);

  /* set thetas */
  theta_Mi = theta_mi = 1e300;
  forvl(vl, vli)
  {
    int iq = Vind(vl, vli);
    double wbar, wMi, wmi, Mi, mi, thM, thm;

    /* get min, max and av on node */
    wmi = dat->ic[iq]->myindc->d[1];
    wMi = dat->ic[iq]->myindc->d[2];
    wbar = dat->ic[iq]->myindc->d[0];

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

    mi = min2(wbar - alpha_h, mi);
    Mi = max2(wbar + alpha_h, Mi);

    /* find min thetas among all vl */
    thM = MRS_theta_Mml(Mi, wbar, wMi);
    thm = MRS_theta_Mml(mi, wbar, wmi);
    if(thM < theta_Mi) theta_Mi = thM;
    if(thm < theta_mi) theta_mi = thm;
  }

  /* set the theta_i we use for limiting vars in vl */
  theta_i = min3(1., theta_mi, theta_Mi);

  /* limit all cons vars q in vl */
  forvl(vl, vli)
  {
    int iq = Vind(vl, vli);
    double *q = Vard(node, iq);
    double qbar;

    /* find node average qbar */
    qbar = var_nodeaverage(node, iq);
//printf("qbar=%g theta_i=%g\n", qbar, theta_i);
//exit(88);
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
