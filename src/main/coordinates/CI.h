/* CI.h */
/* Wolfgang Tichy, 3/2019 */



/*************************************************************************/
/* struct that contains info that helps with coord trafo */
typedef struct tCOORDINFO {
  int iFS[6]; /* Index of var that defines a surface, e.g. FS[0]=sigma of
                 AnsorgNS. Both the vars iSurf[] and FSurf[] are computed from
                 iFS[] by functions in Coordinates. The user should only ever
                 set the var with index iFS[] */
  int iSurf[6];       /* index of surface at grid points, computed from iFS[] */
  int idSurfdX[6][3]; /* index of derivs at grid points, dSurfdX[i][0]=dFS[i]/dX */
  double s[6];        /* some values, e.g. value of surface var in case it is const */
  double xc[3];       /* xc[0..2] = (x,y,z) of coord center for this box */
  int dom;  /* domain index, e.g. 0-5 to indicate cubed sphere wedge */
  int type; /* coordinate type, e.g. outerCubedSphere */
//  int useF; /* if useF=1 we use the funcs below to get values between grid
//               points and to initialize vals inside iSurf and idSurfdX */
  double (*FSurf[6])(struct tNODE *node, int f, double C[2]); /* 6 funcs that return surface val, e.g. FSurf[0]=sigma */
  double (*dFSurfdX[6][4])(struct tNODE *node, int f, double C[2]); /* funcs that return derivs of FSurf, dFSurfdX[0][1]=dSurf[0]/dX */
  struct tARRAY *Fcoef[6]; /* coeffs that FSurf might need */
} tCoordInfo;
