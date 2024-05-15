/* limiter.c */
/* Wolfgang Tichy April, 2019 */


#include "nmesh.h"
#include "limiter.h"


/* global struct with frequently used limiter pars */
tlimiter limiter[1];


/* func to init frequently used pars */
int limiter_init_global_par_indices(tMesh *mesh)
{
  limiter->alpha = Par("limiter_alpha");
  limiter->beta  = Par("limiter_beta");
  limiter->scaleBound  = Par("limiter_scaleBound");
  return 0;
}


/* the limiters in here use the struct tINDIC within tDat to exchange data
   between neighbors. The exchange is done in evolve, which calls functions
   from main/amr/indicators.c to do the exchange. */

/*********************************************************************/
/* funcs for MRS limiter arXiv:1507.03024v1 [math.NA] */
/*********************************************************************/

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

    /* save node average in 1st indicator */
    dat->ic[iq]->myindc->d[0] = var_nodeaverage(node, iq);

    /* find min and max in q in this node, and write it into indicator */
    dat->ic[iq]->myindc->d[1] = min_in_1d_array(q, np, &im);
    dat->ic[iq]->myindc->d[2] = max_in_1d_array(q, np, &im);
  }
  return 0; /* signal that all is ok */
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

/* MRS limiter: limit u using data in dat->ic */
int limiter_MRS(tNode *node, tVarList *vl)
{
  tMesh *mesh = node->pat->mesh;
  double *bb = node->bbox;
  tDat *dat;
  int vli, f, ni, ijk;
  double alpha, h, alpha_h, theta_i;
  double theta_Mi, theta_mi;
  int nlim = 0;

  dat = node->dat;
  if(!dat) return 0;

  /* set alpha_h from alpha (smaller alpha makes MRS more agressive) */
  alpha = Getd(limiter->alpha); //0.1; //5.0;
  h = max3(bb[1]-bb[0], bb[3]-bb[2], bb[5]-bb[4]);
  alpha_h = alpha * pow(h, 1.5);

  /* set thetas */
  theta_Mi = theta_mi = DBL_MAX;
  forvl(vl, vli)
  {
    int iq = Vind(vl, vli);
    double wbar, wMi, wmi, Mi, mi, thM, thm;

    /* get min, max and av on node */
    wmi = dat->ic[iq]->myindc->d[1];
    wMi = dat->ic[iq]->myindc->d[2];
    wbar = dat->ic[iq]->myindc->d[0];

    /* find min and max of q in neighbors */
    Mi = -DBL_MAX;
    mi = DBL_MAX;
    for(f=0; f<6; f++)
      for(ni=0; ni<node->nfnb[f]; ni++)
      {
        double ma = dat->ic[iq]->nbindc[f][ni]->d[1];
        double Ma = dat->ic[iq]->nbindc[f][ni]->d[2];
        if(ma < mi) mi = ma;
        if(Ma > Mi) Mi = Ma;
      }

    //printf("nid%ld: alpha_h=%g  mi=%g Mi=%g\n", Node_eid(node), alpha_h, mi,Mi);
    mi = min2(wbar - alpha_h, mi);
    Mi = max2(wbar + alpha_h, Mi);

    /* find min thetas among all vl */
    thM = MRS_theta_Mml(Mi, wbar, wMi);
    thm = MRS_theta_Mml(mi, wbar, wmi);
    if(thM < theta_Mi) theta_Mi = thM;
    if(thm < theta_mi) theta_mi = thm;
    //printf("  wbar=%g  wmi=%g wMi=%g  mi=%g Mi=%g\n", wbar, wmi,wMi, mi,Mi);
  }

  /* set the theta_i we use for limiting vars in vl */
  theta_i = min3(1., theta_mi, theta_Mi);
  //printf("  theta_mi=%g theta_Mi=%g  theta_i=%g\n", theta_mi,theta_Mi, theta_i);

  /* limit all cons vars q in vl */
  forvl(vl, vli)
  {
    int iq = Vind(vl, vli);
    double *q = Vard(node, iq);
    double qbar;

    /* find node average qbar */
    qbar = var_nodeaverage(node, iq);

    if(!finit(qbar) || !finit(theta_i))
    {
      int *n = node->n;
      pr_nodename(node);
      printf(": n=%d %d %d: vli=%d iq=%d %s\n",
             n[0],n[1],n[2], vli, iq, VarName(iq));
      printf("qbar=%g theta_i=%g\n", qbar, theta_i);
      errorexit("qbar or theta_i is not finite!");
    }
    /* now limit q */
    forpoints(node, ijk)
      q[ijk] = qbar + theta_i*(q[ijk] - qbar);

    /* set non-zero nlim if limiting occured */
    if(theta_i < 1.) nlim++;
  }

  if(nlim) return 1; /* troubled, could decide to set some other bit */
  else     return 0; /* not troubled */
}


/*********************************************************************/
/* funcs for minmodB limiter arXiv:1506.06140v2 [astro-ph.CO] */
/*********************************************************************/

/* set data for average, and Xb-slopes. We actually save the
   coeffs c000, c100, c010, c001 in the basis expansion.
   This one uses the cons vars in vl */
int limdata_c000_100_010_001(tNode *node, tVarList *vl)
{
  tArray *Ac;
  double *c;
  int i100, i010, i001;
  int *n;
  int vli;
  tDat *dat;
  int nvals = 4; /* we need 4 data values per var: node-av + 3 slopes */

  if(!node || !vl) return nvals;

  dat = node->dat;
  if(!dat) return nvals;

  /* array for coeffs */
  n = node->n;
  Ac = alloc_array(n);
  c = Arrd(Ac);

  if(n[0]>1) i100 = Ind_n(1,0,0,n);
  else       i100 = 0;
  if(n[1]>1) i010 = Ind_n(0,1,0,n);
  else       i010 = 0;
  if(n[2]>1) i001 = Ind_n(0,0,1,n);
  else       i001 = 0;

  forvl(vl, vli)
  {
    int iq = Vind(vl, vli);
    tArray *Aq = VarA(node, iq);

    basis_array_analysis3(node, Aq, Ac);

    /* save c000 which has average info */
    dat->ic[iq]->myindc->d[0] = c[0];

    /* save c100,c010,c001 which have slope info */
    dat->ic[iq]->myindc->d[1] = c[i100] * (i100>0);
    dat->ic[iq]->myindc->d[2] = c[i010] * (i010>0);
    dat->ic[iq]->myindc->d[3] = c[i001] * (i001>0);
  }

  free_array(Ac);
  return 0; /* signal that all is ok */
}


/* return min. modulus times sign, or 0 if signs differ */
double minmod3(double a, double b, double c)
{
  if(a*b>=0 && b*c>=0)
    return (a>0. ? 1. : -1.) * min3(fabs(a),fabs(b),fabs(c));
  else
    return 0.;
}

/* bounded minmod */
double minmod3B(double a, double b, double c,  double aBound)
{
  if(fabs(a) <= aBound) return a;
  else                  return minmod3(a,b,c);
}

/* minmodB slope limiter as in : limit u using data in dat->ic */
int limiter_minmodB(tNode *node, tVarList *vl)
{
  tMesh *mesh = node->pat->mesh;
  int *n = node->n;
  double *bb = node->bbox;
  tDat *dat;
  int vli, f, ni;
  int scaleBound = Getb(limiter->scaleBound);
  double alpha, beta, bos3, h, Mt_h, bound;
  const double sqrt3 = sqrt(3.);
  int i100, i010, i001;
  int nlim;
  tArray *Ac;

  dat = node->dat;
  if(!dat) return 0;

  /* alloc array for coeffs */
  Ac = alloc_array(n);

  /* set pars (smaller alpha makes minmodB more agressive):
     alpha corresponds to \tilde{M} of 1506.06140v2, up to a possible
     factor of sqrt(2) because the Legendre poly \bar{P} in 1506.06140v2
     is nornmalized to sqrt(2), while nmesh's basis_normLegendreP is
     normalized to 1. \tilde{M} := M * h, where h is the size of a node,
     and h=\Delta x^K in 1506.06140v2 */
  alpha = Getd(limiter->alpha);
  h     = max3(bb[1]-bb[0], bb[3]-bb[2], bb[5]-bb[4]);
  Mt_h  = alpha * h;
  beta = Getd(limiter->beta);
  bos3 = beta/sqrt3;

  /* locations of coeffs */
  if(n[0]>1) i100 = Ind_n(1,0,0,n);
  else       i100 = 0;
  if(n[1]>1) i010 = Ind_n(0,1,0,n);
  else       i010 = 0;
  if(n[2]>1) i001 = Ind_n(0,0,1,n);
  else       i001 = 0;

  /* default return value */
  nlim = 0;

  /* set weights */
  forvl(vl, vli)
  {
    int iq = Vind(vl, vli);
    double w0[6], w0L1;
    double w0c, w100, w010, w001;
    double wl100, wl010, wl001;

    /* get coeffs with av. and slope info */
    w0c = dat->ic[iq]->myindc->d[0];
    w100 = dat->ic[iq]->myindc->d[1];
    w010 = dat->ic[iq]->myindc->d[2];
    w001 = dat->ic[iq]->myindc->d[3];

    /* get coeffs with av info of neighbors */
    w0L1 = 0.;
    for(f=0; f<6; f++)
    {
      int nnb = node->nfnb[f];
      double av = 0.;

      /* set average nb val on each side */
      for(ni=0; ni<nnb; ni++) av += dat->ic[iq]->nbindc[f][ni]->d[0];
      if(nnb>0) av = av/nnb;
      w0[f] = av;
      w0L1 += fabs(av) * 0.166666666666667;
    }

    /* set bound in bounded minmodB */
    if(scaleBound) bound = Mt_h * w0L1;
    else           bound = Mt_h;

    /* set limited weights, where bos3 = beta/sqrt3 */
    wl100 = minmod3B(w100, bos3*(w0c-w0[0]), bos3*(w0[1]-w0c), bound);
    wl010 = minmod3B(w010, bos3*(w0c-w0[2]), bos3*(w0[3]-w0c), bound);
    wl001 = minmod3B(w001, bos3*(w0c-w0[4]), bos3*(w0[5]-w0c), bound);

    /* limit q */
    if( !dequal(wl100,w100) || !dequal(wl010,w010) || !dequal(wl001,w001) )
    {
      tArray *Aq = VarA(node, iq);
      double *c  = Arrd(Ac);
      int k;

      /* zero coeffs */
      forarray(Ac, k) Arrd_(Ac)[k] = 0.;

      /* set non-zero coeffs */
      c[0] = w0c;
      if(i100>0) c[i100] = wl100;
      if(i010>0) c[i010] = wl010;
      if(i001>0) c[i001] = wl001;

      /* set q from coeffs */
      basis_array_synthesis3(node, Aq, Ac);

      /* make nlim non-zero if limiting occurs */
      nlim++;
    }
  }

  free_array(Ac);

  if(nlim) return 1; /* troubled, could decide to set some other bit */
  else     return 0; /* not troubled */
}
