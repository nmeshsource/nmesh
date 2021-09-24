/* nmesh_limiter.h */
/* Wolfgang Tichy, April 2005 */


/* structure that holds global limiter vars and pars */
typedef struct tLIMITER {

  int alpha;              /* Par("limiter_alpha") */
  int beta;               /* Par("limiter_beta") */
  int scaleBound;         /* Par("limiter_scaleBound") */
} tlimiter;


/* limiter.c */
int limdata_MRS(tNode *node, tVarList *vl);
int limiter_MRS(tNode *node, tVarList *vl);
int limdata_c000_100_010_001(tNode *node, tVarList *vl);
int limiter_minmodB(tNode *node, tVarList *vl);
double MRS_phiy(double y);
double MRS_theta_Mml(double Mi, double wbar, double wMi);
double minmod3(double a, double b, double c);
double minmod3B(double a, double b, double c,  double aBound);
