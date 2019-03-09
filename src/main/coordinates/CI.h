/* CI.h */
/* Wolfgang Tichy, 3/2019 */



struct tPAT; /* incomplete struct, tPAT is completed in main/amr/nmesh_amr.h */

/*************************************************************************/
/* struct that contains info that helps with coord trafo */
typedef struct tCOORDINFO {
  int iSurf[6];       /* index of surface at nodepoints, computed from FSurf[] */
  int idSurfdX[6][3]; /* index of derivs at nodepoints, dSurfdX[i][0]=dFSurf[i]/dX */
  double s[6];        /* some values, e.g. value of surface var in case it is const */
  double xc[3];       /* xc[0..2] = (x,y,z) of coord center for this box */
  int dom;    /* domain index, e.g. 0-5 to indicate cubed sphere wedge */
  int type;   /* coordinate type, e.g. outerCubedSphere */
  int (*FSurf[6])(struct tPAT *pat, int f, double C[2], double *F); /* 6 funcs that set surface val, e.g. FSurf[0] -> F=sigma */
  int (*dFSurfdC[6])(struct tPAT *pat, int f, double C[2], double dFdC[2]); /* set derivs of FSurf, dFSurfdC[4] -> dFdC[0]=dFSurf[4]/dC0 */
  struct tARRAY *Fcoef[6]; /* coeffs that FSurf might need */
} tCoordInfo;
