/* utilities.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "coordinates.h"


/* global vars */
extern tcoordinates coordinates[1];


/* find normal vector (n[0],n[1],n[2]) of patch face f at X,Y,Z */
double patch_normal_at_XYZ(tPat *pat, int f, const double X[3], double n[3])
{
  int dir = f/2;         /* get direction */
  int sig = 2*(f%2) - 1; /* get sign for outward direction */
  double smag;
  int j;

  if(f>5 || f<0) errorexit("f must be 0,1,2,3,4,5");

  if(pat->dXYZ_dxyz) /* not Cartesian */
  {
    double x[3], dXYZdxyz[3][3];

    /* set dXYZdxyz */
    pat->dXYZ_dxyz(pat, NULL,-1, X, x, dXYZdxyz);
    for(j=0; j<3; j++) n[j] = dXYZdxyz[dir][j];
  }
  else  /* use Cartesian normal */
  {
    n[0] = n[1] = n[2] = 0.0;
    n[dir] = 1.0;
  }

  /* normalize and set sign from sig */
  smag = sig * sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  if(smag == 0.0) smag = sig;
  n[0] /= smag;
  n[1] /= smag;
  n[2] /= smag;
  return smag;
}

/* find normal vector (n[0],n[1],n[2]) of node face f at point ijk */
double node_normal_at_ijk(tNode *node, int f, int ijk, double n[3])
{
  tPat *pat = node->pat;
  tMesh *mesh = pat->mesh;
  tDat *dat = node->dat;
  int dir = f/2;     /* get direction */
  int sig = 2*(f%2) - 1; /* get sign for outward direction */
  int idXd = Ind("dXdx");
  double smag;

  if(f>5 || f<0) errorexit("f must be 0,1,2,3,4,5");

  if(!dat) return 0.;

  /* do we need to init. coords? */
  if(!(dat->coords_set)) coordinates_init_node(node);

  /* get normal from derivs */
  {
    double *pdXd[3][3]
            = { {Vard(node,idXd),   Vard(node,idXd+1), Vard(node,idXd+2)},
                {Vard(node,idXd+3), Vard(node,idXd+4), Vard(node,idXd+5)},
                {Vard(node,idXd+6), Vard(node,idXd+7), Vard(node,idXd+8)} };
    n[0] = pdXd[dir][0][ijk];
    n[1] = pdXd[dir][1][ijk];
    n[2] = pdXd[dir][2][ijk];
  }

  /* normalize and set sign from sig */
  smag = sig * sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  if(smag == 0.0) smag = sig;
  n[0] /= smag;
  n[1] /= smag;
  n[2] /= smag;
  return smag;
}

/* find normal vector (n[0],n[1],n[2]) on face f at midpoint ijk */
double node_normal_at_midpt_ijk(tNode *node, int f, int ijk, double n[3])
{
  tPat *pat = node->pat;
  tMesh *mesh = pat->mesh;
  tDat *dat = node->dat;
  int dir = f/2;     /* get direction */
  int sig = 2*(f%2) - 1; /* get sign for outward direction */
  int idXd;
  double smag;

  if(f>5 || f<0) errorexit("f must be 0,1,2,3,4,5");

  if(!dat) return 0.;

  /* get correct var index */
  if(dir==0)        idXd = Ind("Xm_dXdx");
  else if (dir==1)  idXd = Ind("Ym_dXdx");
  else              idXd = Ind("Zm_dXdx");

  /* do we need to init. coords? */
  if(!(dat->coords_set)) coordinates_init_node(node);

  /* get normal from derivs */
  {
    double *pdXd[3][3]
            = { {Vard(node,idXd),   Vard(node,idXd+1), Vard(node,idXd+2)},
                {Vard(node,idXd+3), Vard(node,idXd+4), Vard(node,idXd+5)},
                {Vard(node,idXd+6), Vard(node,idXd+7), Vard(node,idXd+8)} };
    n[0] = pdXd[dir][0][ijk];
    n[1] = pdXd[dir][1][ijk];
    n[2] = pdXd[dir][2][ijk];
  }

  /* normalize and set sign from sig */
  smag = sig * sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  if(smag == 0.0) smag = sig;
  n[0] /= smag;
  n[1] /= smag;
  n[2] /= smag;
  return smag;
}


/* get transpose of a 3x3 matrix. */
void transp3Dmat_from_3Dmat(CONST double M[3][3], double transpM[3][3])
{
  transpM[0][0] = M[0][0];
  transpM[0][1] = M[1][0];
  transpM[0][2] = M[2][0];
  transpM[1][0] = M[0][1];
  transpM[1][1] = M[1][1];
  transpM[1][2] = M[2][1];
  transpM[2][0] = M[0][2];
  transpM[2][1] = M[1][2];
  transpM[2][2] = M[2][2];
}

/* determinant of 3x3 matrix */
double det_3Dmatrix(CONST double M[3][3])
{
  /* M = {{m11,m12,m13},{m21,m22,m23},{m31,m32,m33}}
     Det[M] = m11*m22*m33 - m13*m22*m31 + m12*m23*m31 +
              m13*m21*m32 - m11*m23*m32 - m12*m21*m33   */
  LDOUBLE DetM =M[0][0]*M[1][1]*M[2][2] -
                M[0][2]*M[1][1]*M[2][0] +
                M[0][1]*M[1][2]*M[2][0] +
                M[0][2]*M[1][0]*M[2][1] -
                M[0][0]*M[1][2]*M[2][1] -
                M[0][1]*M[1][0]*M[2][2];
  return DetM;
}

/* 3d matrix inverse and return det. */
double inv3Dmat_from_3Dmat(CONST double M[3][3], double invM[3][3])
{
  /* M = {{m11,m12,m13},{m21,m22,m23},{m31,m32,m33}}
     Det[M] = m11*m22*m33 - m13*m22*m31 + m12*m23*m31 +
              m13*m21*m32 - m11*m23*m32 - m12*m21*m33
     Inverse[M]*Det[M] =
       {{-(m23*m32) + m22*m33, m13*m32 - m12*m33, -(m13*m22) + m12*m23},
        {m23*m31 - m21*m33, -(m13*m31) + m11*m33, m13*m21 - m11*m23},
        {-(m22*m31) + m21*m32, m12*m31 - m11*m32, -(m12*m21) + m11*m22}}  */
  LDOUBLE DetM = det_3Dmatrix(M);
  LDOUBLE detM, sum;
  if(DetM==0.0) detM = dequaleps*dequaleps*dequaleps;
  else          detM = DetM;
  sum        = (-(M[1][2]*M[2][1]) + M[1][1]*M[2][2])/detM;
  invM[0][0] = sum;
  sum        = (  M[0][2]*M[2][1]  - M[0][1]*M[2][2])/detM;
  invM[0][1] = sum;
  sum        = (-(M[0][2]*M[1][1]) + M[0][1]*M[1][2])/detM;
  invM[0][2] = sum;
  sum        = (  M[1][2]*M[2][0]  - M[1][0]*M[2][2])/detM;
  invM[1][0] = sum;
  sum        = (-(M[0][2]*M[2][0]) + M[0][0]*M[2][2])/detM;
  invM[1][1] = sum;
  sum        = (  M[0][2]*M[1][0]  - M[0][0]*M[1][2])/detM;
  invM[1][2] = sum;
  sum        = (-(M[1][1]*M[2][0]) + M[1][0]*M[2][1])/detM;
  invM[2][0] = sum;
  sum        = (  M[0][1]*M[2][0]  - M[0][0]*M[2][1])/detM;
  invM[2][1] = sum;
  sum        = (-(M[0][1]*M[1][0]) + M[0][0]*M[1][1])/detM;
  invM[2][2] = sum;

  return DetM;
}

/* print 3x3 matrix */
void print3Dmat(CONST double M[3][3])
{
  int i,j;
  for(i=0; i<3; i++)
  {
    for(j=0; j<3; j++) printf("%15g", M[i][j]);
    printf("\n");
  }
}

/* det of symm. matrix */
double det_3Dsymmmat(double M11, double M12, double M13,
                     double M22, double M23, double M33)
{
  double MMinv11, MMinv12, MMinv13;

  MMinv11 = M22*M33 - M23*M23;
  MMinv12 = M13*M23 - M12*M33;
  MMinv13 = M12*M23 - M13*M22;
  return M11*MMinv11 + M12*MMinv12 + M13*MMinv13;
}

/* inverse of symm. matrix, returns det */
double inv2Dmat_from_2Dsymmmat(double  s11, double  s12, double  s22,
                               double *i11, double *i12, double *i22)
{
/*
var('s11 s12 s22')
Smat = matrix([[s11,s12],[s12,s22]])
Sinvdet = Smat.I * Smat.det()
# idet11 = Sinvdet[0,0].simplify_full()
# ...
invS = matrix([[i11,i12],[i12,i22]])
*/
  double detS, idet11,idet12,idet22;

  idet11 = s22;
  idet12 = -s12;
  idet22 = s11;

  detS = s11*idet11 + s12*idet12;

  *i11 = idet11/detS;
  *i12 = idet12/detS;
  *i22 = idet22/detS;

  return detS;
}

/* inverse of symm. matrix, returns det */
double inv3Dmat_from_3Dsymmmat(double M11, double M12, double M13,
                               double M22, double M23, double M33,
                               double *i11, double *i12, double *i13,
                               double *i22, double *i23, double *i33)
{
  double detM, MMinv11, MMinv12, MMinv13, MMinv22, MMinv23, MMinv33;

  MMinv11 = M22*M33 - M23*M23;
  MMinv12 = M13*M23 - M12*M33;
  MMinv13 = M12*M23 - M13*M22;
  MMinv22 = M11*M33 - M13*M13;
  MMinv23 = M12*M13 - M11*M23;
  MMinv33 = M11*M22 - M12*M12;
  detM = M11*MMinv11 + M12*MMinv12 + M13*MMinv13;
  *i11 = MMinv11/detM;
  *i12 = MMinv12/detM;
  *i13 = MMinv13/detM;
  *i22 = MMinv22/detM;
  *i23 = MMinv23/detM;
  *i33 = MMinv33/detM;
  return detM;
}

/* det of symm. matrix */
double det_4Dsymmmat(double s11,double s12,double s13,double s14,
                     double s22, double s23, double s24,
                     double s33, double s34, double s44)
{
/*
var('s11 s12 s13 s14 s22 s23 s24 s33 s34 s44')
Smat = matrix([[s11,s12,s13,s14],[s12,s22,s23,s24],[s13,s23,s33,s34],
[s14,s24,s34,s44]])
Sinvdet = Smat.I * Smat.det()
# idet11 = Sinvdet[0,0].simplify_full()
# ...
invS = matrix([[i11,i12,i13,i14],[i12,i22,i23,i24],[i13,i23,i33,i34],
[i14,i24,i34,i44]])
*/
  double detS, idet11,idet12,idet13,idet14;

  idet11 = -s24*s24*s33 + 2*s23*s24*s34 - s22*s34*s34 - (s23*s23 - s22*s33)*s44;
  idet12 = s14*s24*s33 + s12*s34*s34 - (s14*s23 + s13*s24)*s34 + (s13*s23 - s12*s33)*s44;
  idet13 = -s14*s23*s24 + s13*s24*s24 + (s14*s22 - s12*s24)*s34 - (s13*s22 - s12*s23)*s44;
  idet14 = s14*s23*s23 - s13*s23*s24 - (s14*s22 - s12*s24)*s33 + (s13*s22 - s12*s23)*s34;

  detS = s11*idet11 + s12*idet12 + s13*idet13 + s14*idet14;

  return detS;
}

/* inverse of symm. matrix, returns det */
double inv4Dmat_from_4Dsymmmat(double s11,double s12,double s13,double s14,
                               double s22, double s23, double s24,
                               double s33, double s34, double s44,
                               double *i11,double *i12,double *i13,double *i14,
                               double *i22, double *i23, double *i24,
                               double *i33, double *i34, double *i44)
{
/*
var('s11 s12 s13 s14 s22 s23 s24 s33 s34 s44')
Smat = matrix([[s11,s12,s13,s14],[s12,s22,s23,s24],[s13,s23,s33,s34],
[s14,s24,s34,s44]])
Sinvdet = Smat.I * Smat.det()
# idet11 = Sinvdet[0,0].simplify_full()
# ...
invS = matrix([[i11,i12,i13,i14],[i12,i22,i23,i24],[i13,i23,i33,i34],
[i14,i24,i34,i44]])
*/
  double detS, idet11,idet12,idet13,idet14;
  double       idet22,idet23,idet24, idet33,idet34, idet44;

  idet11 = -s24*s24*s33 + 2*s23*s24*s34 - s22*s34*s34 - (s23*s23 - s22*s33)*s44;
  idet12 = s14*s24*s33 + s12*s34*s34 - (s14*s23 + s13*s24)*s34 + (s13*s23 - s12*s33)*s44;
  idet13 = -s14*s23*s24 + s13*s24*s24 + (s14*s22 - s12*s24)*s34 - (s13*s22 - s12*s23)*s44;
  idet14 = s14*s23*s23 - s13*s23*s24 - (s14*s22 - s12*s24)*s33 + (s13*s22 - s12*s23)*s34;
  idet22 = -s14*s14*s33 + 2*s13*s14*s34 - s11*s34*s34 - (s13*s13 - s11*s33)*s44;
  idet23 = s14*s14*s23 - s13*s14*s24 - (s12*s14 - s11*s24)*s34 + (s12*s13 - s11*s23)*s44;
  idet24 = -s13*s14*s23 + s13*s13*s24 + (s12*s14 - s11*s24)*s33 - (s12*s13 - s11*s23)*s34;
  idet33 = -s14*s14*s22 + 2*s12*s14*s24 - s11*s24*s24 - (s12*s12 - s11*s22)*s44;
  idet34 = s13*s14*s22 - s12*s14*s23 - (s12*s13 - s11*s23)*s24 + (s12*s12 - s11*s22)*s34;
  idet44 = -s13*s13*s22 + 2*s12*s13*s23 - s11*s23*s23 - (s12*s12 - s11*s22)*s33;

  detS = s11*idet11 + s12*idet12 + s13*idet13 + s14*idet14;

  *i11 = idet11/detS;
  *i12 = idet12/detS;
  *i13 = idet13/detS;
  *i14 = idet14/detS;
  *i22 = idet22/detS;
  *i23 = idet23/detS;
  *i24 = idet24/detS;
  *i33 = idet33/detS;
  *i34 = idet34/detS;
  *i44 = idet44/detS;

  return detS;
}

/* store inverse of symm. S in I */
double invmat_from_symmmat(int n, const double *S, double *I)
{
  switch(n)
  {
  case 1:
    I[0] = 1./S[0];
    return S[0];
  case 2:
    return inv2Dmat_from_2Dsymmmat(S[0], S[1], S[2],  I, I+1, I+2);
  case 3:
    return inv3Dmat_from_3Dsymmmat(S[0], S[1], S[2], S[3], S[4], S[5],
                                   I,    I+1,  I+2,  I+3,  I+4,  I+5);
  case 4:
    return inv4Dmat_from_4Dsymmmat(S[0], S[1], S[2], S[3],
                                         S[4], S[5], S[6],  S[7],S[8],  S[9],
                                   I, I+1, I+2, I+3,
                                      I+4, I+5, I+6,  I+7,I+8,  I+9);
  default:
    errorexit("n must be 1,2,3,4");
  }
}

/* return memory index(a,b) for a symm. n*n matrix that is stored as a 1d
   C-array, e.g. M_{ab} could be stored as [Mxx,Mxy,Mxz,Myy,Myz,Mzz]
   or [Mtt,Mtx,Mty,Mtz,Mxx,Mxy,Mxz,Myy,Myz,Mzz]. Here a,b \in [0,n-1]. */
int index_symmmat(int n, int a, int b)
{
  int c,d, oc;

  /* swap indices a,b s.t. new indices c,d have c<=d */
  if(a>b) { c=b; d=a; }
  else    { c=a; d=b; }

  /* return memory index of M_{cd}:  */
  /* nc = n-c;
     oc = ((n+1)*n)/2 - ((nc+1)*nc)/2 - c // Note: \sum_{k=1}^N k = ((N+1)*N)/2
        = ( (n+1)*n - (nc+1)*nc )/2 - c = ( n*n + n - nc*nc - nc - 2*c )/2
        = ( 2*n*c - c*c + c - 2*c)/2 = ( c*(2*n - c - 1) )/2 */
  oc = ( c*(2*n - c - 1) )/2;
  //PRF;printf(": n=%d %d%d -> %d%d: oc+d=%d\n", n,a,b,c,d, oc + d);
  return oc + d; /* M_{ab} = M_{cd} */
}

/* return element M_{ab} of a symm. n*n matrix that is stored as a 1d
   C-array, e.g. Mab could be stored as [Mxx,Mxy,Mxz,Myy,Myz,Mzz]
   or [Mtt,Mtx,Mty,Mtz,Mxx,Mxy,Mxz,Myy,Myz,Mzz]. Here a,b \in [0,n-1]. */
double matel_symmmat(int n, const double *M, int a, int b)
{
  return M[index_symmmat(n, a,b)];
}

/* compute Mv_a = M_{ab} v^b for for a symm. n*n matrix that is stored as a
   1d C-array, e.g. M_{ab} could be stored as [Mxx,Mxy,Mxz,Myy,Myz,Mzz]
   or [Mtt,Mtx,Mty,Mtz,Mxx,Mxy,Mxz,Myy,Myz,Mzz]. Here a,b \in [0,n-1]. */
void symmmat_times_vec(int n, const double *M, const double *v, double *Mv)
{
  int a,b;
  for(a=0; a<n; a++)
  {
    double sum=0.;
    for(b=0; b<n; b++) sum += M[index_symmmat(n, a,b)] * v[b];
    Mv[a] = sum;
  }
}

/* compute u^a v_a */
double vec_times_vec(int n, const double *u, const double *v)
{
  int a;
  double sum=0.;
  for(a=0; a<n; a++) sum += u[a] * v[a];
  return sum;
}

/* compute g_{ab} v^a v^b */
double mag2_vector_metric(int n, const double *g, const double *v)
{
  double gv[n], prod;
  symmmat_times_vec(n, g, v, gv);
  prod = vec_times_vec(n, v, gv);
  return prod;
}

/* compute BM = B_{ab} M^{ab} */
double BM_symmmat(int n, const double *B, const double *M)
{
  int a,b;
  double sum=0.;

  for(a=0; a<n; a++)
  for(b=0; b<n; b++)
  {
    sum += B[index_symmmat(n, a,b)] * M[index_symmmat(n, a,b)];
  }
  return sum;
}

/* compute BMB_{ab} = B_a^c B_b^d M_{cd} */
void BMB_symmmat(int n, const double *B, const double *M, double *BMB)
{
  int a,b, c,d;

  for(a=0; a<n; a++)
  for(b=a; b<n; b++)
  {
    double sum=0.;
    for(c=0; c<n; c++)
    for(d=0; d<n; d++)
      sum += B[index_symmmat(n, a,c)]*B[index_symmmat(n, b,d)]*
             M[index_symmmat(n, c,d)];
    BMB[index_symmmat(n, a,b)] = sum;
  }
}

/* compute BBBM_{abc} = B_a^d B_b^e B_c^f M_{def} for one fixed a,
   where M_{def} = M_{dfe}  ==>  BBBM_{abc} = BBBM_{acb} */
void BBBMa_symm_bc(int n, const double *B, const double *M, double *BBBMa,
                   int a)
{
  int b,c, d,e,f;
  int ns = ( (n+1)*n )/2; /* number of elems in symm matrix */

  for(b=0; b<n; b++)
  for(c=b; c<n; c++)
  {
    double sum=0.;
    for(d=0; d<n; d++)
    for(e=0; e<n; e++)
    for(f=0; f<n; f++)
      sum += B[index_symmmat(n, a,d)]*B[index_symmmat(n, b,e)]*
             B[index_symmmat(n, c,f)]*M[ns*d + index_symmmat(n, e,f)];
    BBBMa[index_symmmat(n, b,c)] = sum;
  }
}
/* compute BBBM_{abc} = B_a^d B_b^e B_c^f M_{def},
   where M_{def} = M_{dfe}  ==>  BBBM_{abc} = BBBM_{acb} */
void BBBM_symm_bc(int n, const double *B, const double *M, double *BBBM)
{
  int a;
  int ns = ( (n+1)*n )/2; /* number of elems in symm matrix */

  for(a=0; a<n; a++)
  {
    BBBMa_symm_bc(n, B, M, BBBM + ns*a, a);
  }
}

/* compute BBBBM_{abcd} = B_a^e B_b^f B_c^g B_d^h M_{efgh} for fixed a,b
   where M_{efgh} = M_{efhg} and M_{efgh} = M_{fegh}
   ==>  BBBBM_{abcd} = BBBM_{abdc} and BBBBM_{abcd} = BBBBM_{bacd} */
void BBBBMab_symm_ab_cd(int n, const double *B, const double *M,
                        double *BBBBMab, int a, int b)
{
  int c,d, e,f,g,h;
  int ns = ( (n+1)*n )/2; /* number of elems in symm matrix */

  for(c=0; c<n; c++)
  for(d=c; d<n; d++)
  {
    double sum=0.;
    for(e=0; e<n; e++)
    for(f=0; f<n; f++)
    for(g=0; g<n; g++)
    for(h=0; h<n; h++)
      sum += B[index_symmmat(n, a,e)]*B[index_symmmat(n, b,f)]*
             B[index_symmmat(n, c,g)]*B[index_symmmat(n, d,h)]*
             M[ns*index_symmmat(n, e,f) + index_symmmat(n, g,h)];
    BBBBMab[index_symmmat(n, c,d)] = sum;
  }
}
/* compute BBBBM_{abcd} = B_a^e B_b^f B_c^g B_d^h M_{efgh} where
   M_{efgh} = M_{efhg} and M_{efgh} = M_{fegh}
   ==>  BBBBM_{abcd} = BBBM_{abdc} and BBBBM_{abcd} = BBBBM_{bacd} */
void BBBBM_symm_ab_cd(int n, const double *B, const double *M, double *BBBBM)
{
  int a,b;
  int ns = ( (n+1)*n )/2; /* number of elems in symm matrix */

  for(a=0; a<n; a++)
  for(b=a; b<n; b++)
  {
    BBBBMab_symm_ab_cd(n, B, M, BBBBM + ns*index_symmmat(n, a,b), a,b);
  }
}


/* invert a 2*2*1 symmetric array in place and return det */
double invert2x2x1symm_array(tArray *a)
{
  double i11,i12,i22;
  double *d = Arrd_(a);
  double det = inv2Dmat_from_2Dsymmmat(d[0],d[1],d[3], &i11,&i12,&i22);

  d[0] = i11;
  d[1] = i12;

  d[2] = i12;
  d[3] = i22;

  //PRF;printf(": det=%g\n", det);
  return det;
}

/* invert a 3*3*1 symmetric array in place and return det */
double invert3x3x1symm_array(tArray *a)
{
  double i11,i12,i13, i22,i23, i33;
  double *d = Arrd_(a);
  double det = inv3Dmat_from_3Dsymmmat(d[0],d[1],d[2], d[4],d[5], d[8],
                                       &i11,&i12,&i13, &i22,&i23, &i33);
  d[0] = i11;
  d[1] = i12;
  d[2] = i13;

  d[3] = i12;
  d[4] = i22;
  d[5] = i23;

  d[6] = i13;
  d[7] = i23;
  d[8] = i33;

  //PRF;printf(": det=%g\n", det);
  return det;
}

/* invert a 4*4*1 symmetric array in place and return det */
double invert4x4x1symm_array(tArray *a)
{
  double i11,i12,i13,i14, i22,i23,i24, i33,i34, i44;
  double *d = Arrd_(a);
  double det = inv4Dmat_from_4Dsymmmat(d[0],d[1],d[2],d[3],
                                       d[5],d[6],d[7], d[10],d[11], d[15],
                                       &i11,&i12,&i13,&i14,
                                       &i22, &i23, &i24, &i33, &i34, &i44);
  d[0]  = i11;
  d[1]  = i12;
  d[2]  = i13;
  d[3]  = i14;
  d[5]  = i22;
  d[6]  = i23;
  d[7]  = i24;
  d[10] = i33;
  d[11] = i34;
  d[15] = i44;

  d[4]  = i12;
  d[8]  = i13;
  d[9]  = i23;
  d[12] = i14;
  d[13] = i24;
  d[14] = i34;

  //PRF;printf(": det=%g\n", det);
  return det;
}

/* invert a N*N*1 symmetric array in place and return det */
double invertNxNx1symm_array(tArray *a)
{
  double det;

  /* note we expect that: a->N = N*N */
  switch(a->N)
  {
  case 1:
    det = a->d[0];
    a->d[0] = 1./det;
    return det;
  case 4:
    return invert2x2x1symm_array(a);
  case 9:
    return invert3x3x1symm_array(a);
  case 16:
    return invert4x4x1symm_array(a);
  default:
    errorexit("implement N>4 cases");
  }
}

/* compute V_i from V^i using metric g_{ij}, or V^i from V_i using
   metric g^{ij}, returns V^2 = V_i V^i */
double symmmat3D_times_vec(double gxx, double gxy, double gxz,
                           double gyy, double gyz, double gzz,
                           double Vx, double Vy, double Vz,
                           double *gVx, double *gVy, double *gVz)
{
  *gVx = gxx*Vx + gxy*Vy + gxz*Vz;
  *gVy = gxy*Vx + gyy*Vy + gyz*Vz;
  *gVz = gxz*Vx + gyz*Vy + gzz*Vz;
  return Vx*(*gVx) + Vy*(*gVy) + Vz*(*gVz);
}

/* compute V^2 = V_i V^i from metric g_{ij} and V^i, or from g^{ij} and V_i */
double mag2_vector_3Dmetric(double gxx, double gxy, double gxz,
                            double gyy, double gyz, double gzz,
                            double Vx, double Vy, double Vz)
{
  double gVx, gVy, gVz;
  gVx = gxx*Vx + gxy*Vy + gxz*Vz;
  gVy = gxy*Vx + gyy*Vy + gyz*Vz;
  gVz = gxz*Vx + gyz*Vy + gzz*Vz;
  return Vx*gVx + Vy*gVy + Vz*gVz;
}

/* find pat size L of pat */
double find_pat_size(tPat *pat)
{
  int i,j,n;
  double dx,dy,dz;
  double sum, L;
  double *bb = pat->bbox;
  double X[8][3], x[8][3];

  /* max and min X */
  X[0][0] = bb[0];  X[0][1] = bb[2];  X[0][2] = bb[4];
  X[1][0] = bb[1];  X[1][1] = bb[2];  X[1][2] = bb[4];
  X[2][0] = bb[0];  X[2][1] = bb[3];  X[2][2] = bb[4];
  X[3][0] = bb[1];  X[3][1] = bb[3];  X[3][2] = bb[4];
  X[4][0] = bb[0];  X[4][1] = bb[2];  X[4][2] = bb[5];
  X[5][0] = bb[1];  X[5][1] = bb[2];  X[5][2] = bb[5];
  X[6][0] = bb[0];  X[6][1] = bb[3];  X[6][2] = bb[5];
  X[7][0] = bb[1];  X[7][1] = bb[3];  X[7][2] = bb[5];

  /* x at max and min X */
  for(i=0; i<8; i++)
    set_xyz(pat, NULL,-1, X[i], x[i]);

  sum=0.0;
  n=0;
  for(i=1; i<8; i++)
  for(j=0; j<i; j++)
  {
    dx = x[i][0] - x[j][0];
    dy = x[i][1] - x[j][1];
    dz = x[i][2] - x[j][2];

    sum += dx*dx + dy*dy + dz*dz;
    n++;
  }

  L = sqrt(sum/n);
  return L;
}

/* find pat size L of smallest pat */
double smallest_pat_size(tMesh *mesh)
{
  int i;
  double L=DBL_MAX;
  forpatches(mesh, i)
  {
    tPat *pat = mesh->pat[i];
    double Li = find_pat_size(pat);
    if(Li<L) L=Li;
  }
  return L;
}


/* put dX^i/dx^j at point ind into array 3x3 aJ */
void array_dXdx(tNode *node, int ind, tArray *aJ)
{
  tMesh *mesh = node->pat->mesh;
  int idXdx = Ind("dXdx");
  double val;
  int i,j;

  for(i=0; i<3; i++)
  for(j=0; j<3; j++)
  {
    val = Vard_(node, idXdx + 3*i + j)[ind];
    Arrd_(aJ)[i + 3*j] = val; /* aJ is stored in column major form */
  }
}

/* put dx^i/dX^j at point ind into 3x3 array aJ */
void array_dxdX(tNode *node, int ind, tArray *aJ)
{
  tMesh *mesh = node->pat->mesh;
  int idXdx = Ind("dXdx");
  double M[3][3], invM[3][3];
  int i,j;

  for(i=0; i<3; i++)
  for(j=0; j<3; j++)
    M[i][j] = Vard_(node, idXdx + 3*i + j)[ind];

  inv3Dmat_from_3Dmat(M, invM);

  for(i=0; i<3; i++)
  for(j=0; j<3; j++)
    Arrd_(aJ)[i + 3*j] = invM[i][j];
}

/* put dx^i/dX^A and dx^i/dX^B at point ind into 3x2 array aJ,
   A,B are e.g. 0,2 if norm=1 */
void array_2dxdX(tNode *node, int ind, int norm, tArray *aJ)
{
  tMesh *mesh = node->pat->mesh;
  int idXdx = Ind("dXdx");
  double M[3][3], invM[3][3];
  int i,j, J;

  for(i=0; i<3; i++)
  for(j=0; j<3; j++)
    M[i][j] = Vard_(node, idXdx + 3*i + j)[ind];

  inv3Dmat_from_3Dmat(M, invM);

  for(J=0, j=0; j<3; j++)
    if(j!=norm)
    {
      for(i=0; i<3; i++) Arrd_(aJ)[i + 3*J] = invM[i][j];
      J++;
    }
}

/* put dx^i/dXb^A and dx^i/dXb^B at point ind into 3x2 array aJ,
   A,B are e.g. 0,2 if norm=1 */
void array_2dxdXb(tNode *node, int ind, int norm, tArray *aJ)
{
  tMesh *mesh = node->pat->mesh;
  int idXdx = Ind("dXdx");
  double M[3][3], invM[3][3];
  double dXbdX[3];
  int i,j, J;

  /* get dXb/dX */
  dXbYbZb_dXYZ(node, dXbdX);

  /* set M = dX^i/dx^j */
  for(i=0; i<3; i++)
  for(j=0; j<3; j++)
    M[i][j] = Vard_(node, idXdx + 3*i + j)[ind];

  /* set invM = dx^i/dX^j */
  inv3Dmat_from_3Dmat(M, invM);
/*
print3Dmat(M);

//printnode(node);
//printvar_innode(node, Ind("Y"));
//printvar_innode(node, Ind("x"));
//printvar_innode(node, Ind("y"));
//printvar_innode(node, Ind("z"));
tArray *A = alloc_array2d(3,3);
point_array_d_to_data(A, M, 1);
printarray(A);

print3Dmat(invM);
point_array_d_to_data(A, invM, 1);
printarray(A);
for(i=0; i<3; i++) printf("invM[i][1]=%g\n" , invM[i][1]);
for(i=0; i<3; i++) printf("invM[i][2]=%g\n" , invM[i][2]);
*/
  /* save dx^i/dXb^j = dx^i/dX^j / (dXb^j/dX^j) in array aJ */
  for(J=0, j=0; j<3; j++)
    if(j!=norm)
    {
      for(i=0; i<3; i++) Arrd_(aJ)[i + 3*J] = invM[i][j]/dXbdX[j];
      J++;
    }
}

/* put 3x3 matrix M in array aM */
void array_3Dmat(CONST double M[3][3], tArray *aM)
{
  int i,j;

  for(i=0; i<3; i++)
  for(j=0; j<3; j++)
    Arrd_(aM)[i + 3*j] = M[i][j];
}

/* get determinant of a 2x2 array */
double det_2_2_array(tArray *aM)
{
  double *d = Arrd_(aM);
  return d[0]*d[3] - d[1]*d[2];
}

/* divide or multiply var u by Jacobian J,
   if Jpow=1: u -> J * u,
   if Jpow=-1: u -> J^{-1} * u
   if Jpow=0: do nothing */
void var_to_var_times_JtoPower(tNode *node, int ui, int Jpower)
{
  int iooJ;
  double *ooJ, *u;
  int i;

  if(Jpower==0) return;

  iooJ = coordinates->idet_dXbdx;
  ooJ = Vard_(node, iooJ); /* contains 1/J */
  u   = Vard_(node, ui);

  switch(Jpower)
  {
  case 1:
    forvari(node, ui, i) u[i] = u[i] / ooJ[i];
    break;
  case -1:
    forvari(node, ui, i) u[i] = u[i] * ooJ[i];
    break;
  default:
    errorexit("Jpower must be 1, -1 or 0");
  }
}


/* find nb. point index in dir, pm = +/- 1 */
int nb_point_in_dir(int ijk, int *n, int pm, int dir)
{
  if(n[dir]>1)
  {
    switch(dir)
    {
    case 0:
      return ijk + pm;
    case 1:
      return ijk + pm*n[0];
    case 2:
      return ijk + pm*n[0]*n[1];
    default:
      errorexit("dir must be 0,1,2");
    }
  }
  else
    return ijk; /* we return ijk if there is no nb. */
}

/* find new hmin given point ccc and nb point nbr */
double hmin_new_ccc_nbr(tNode *node, int ccc, int nbr, double hmin_old)
{
  errorexit("works only if node has dat and var x is enabled");
  if(nbr!=ccc)
  {
    tPat *pat = node->pat;
    tMesh *mesh = pat->mesh;
    int ix = pat->xyz_of_XYZ ? Ind("x") : Ind("X");
    double *x = Vard(node, ix);
    double *y = Vard(node, ix+1);
    double *z = Vard(node, ix+2);
    double r, dx,dy,dz;

    dx = x[nbr] - x[ccc];
    dy = y[nbr] - y[ccc];
    dz = z[nbr] - z[ccc];
    r = sqrt(dx*dx + dy*dy + dz*dz);
    if(r < hmin_old) return r;
  }
  return hmin_old;
}

/* trivial Cartesian distance between two points in x-coords */
double Cart_distance_x0_x1(tNode *node, double x0[3], double x1[3])
{
  int d;
  double dist2;

  dist2 = 0.;
  for(d=0; d<3; d++)
  {
    double temp = x1[d] - x0[d];
    dist2 += temp*temp;
  }
  return sqrt(dist2);
}

/* find Cartesian distance between two points in X-coords */
double Cart_distance_X0_X1(tNode *node, double X0[3], double X1[3])
{
  tPat *pat = node->pat;
  int d;
  double c0[3], c1[3], dx[3];
  double *x0, *x1;
  double dist2;

  /* Note: calling xyz_of_XYZ may be slow. We could just also pass in
     the point indices and then use the vars x,y,z. */

  /* set x0, x1 from X0, X1 */
  if(pat->xyz_of_XYZ)
  {
    x0 = &(c0[0]);
    x1 = &(c1[0]);
    pat->xyz_of_XYZ(pat,node,-1, X0, x0);
    pat->xyz_of_XYZ(pat,node,-1, X1, x1);
  }
  else /* assume Cartesian */
  {
    x0 = &(X0[0]);
    x1 = &(X1[0]);
  }

  /* get distance from distance^2 */
  dist2 = 0.;
  for(d=0; d<3; d++)
  {
    dx[d] = x1[d] - x0[d];
    dist2 += dx[d]*dx[d];
  }
  //PRFs(": ");pr_nodename(node);
  //printf(" pat->xyz_of_XYZ=%p ", pat->xyz_of_XYZ);
  //pr3v("x0", x0); pr3v("x1", x1);
  //printf("-> dist=%g\n", sqrt(dist2));

  return sqrt(dist2);
}

/* Find distance to closest point, but return hmin_old if no point closer
   than hmin_old found. If a point closer than hmin_old exists, write its
   index into *ijk1, and the index of i,j,k into *ijk0,
   otherwise do not change *ijk0 and *ijk1 */
double distance_to_closest_point(tNode *node, int i, int j, int k,
                                 double hmin_old, int *ijk0, int *ijk1)
{
  int *n = node->n;
  tArray *XbA[3];
  double Xb0[3], Xb1[3], X0[3], X1[3];
  int ii, jj, kk; /* increments to i,j,k to get to nb. points */
  double dist;
  double hmin = hmin_old;

  node_Xb3(node, XbA); /* set node point arrays */

  Xb0[0] = XbA[0]->d[i];
  Xb0[1] = XbA[1]->d[j];
  Xb0[2] = XbA[2]->d[k];
  XYZ_of_XbYbZb(node, Xb0, X0);

  for(kk=-1; kk<=1; kk++)
  {
    int kn = k+kk;
    if(kn<0 || kn>=n[2]) continue;

    for(jj=-1; jj<=1; jj++)
    {
      int jn = j+jj;
      if(jn<0 || jn>=n[1]) continue;

      for(ii=-1; ii<=1; ii++)
      {
        int in = i+ii;
        if(in<0 || in>=n[0]) continue;

        if(in==i && jn==j && kn==k) continue;

        if(in!=i && jn!=j) continue;
        if(in!=i && kn!=k) continue;
        if(jn!=j && kn!=k) continue;

        Xb1[0] = XbA[0]->d[in];
        Xb1[1] = XbA[1]->d[jn];
        Xb1[2] = XbA[2]->d[kn];
        XYZ_of_XbYbZb(node, Xb1, X1);

        dist = Cart_distance_X0_X1(node, X0,X1);
        if(dist<hmin)
        {
          hmin = dist;
          *ijk0 = Ind_n(i,j,k, n);
          *ijk1 = Ind_n(in,jn,kn, n);
          //PRFs(": ");pr_nodename(node);
          //pr3v(" X0", X0); pr3v("X1", X1);
          //printf("pts=%d,%d -> dist=%g\n", *ijk0,*ijk1, dist);
        }
      }
    }
  }
  return hmin;
}


/* find smallest Cartesian grid spacing in all 8 nnode corners */
double find_hmin(tNode *node, int *ijk0, int *ijk1)
{
  int *n = node->n;
  double X0[] = { node->bbox[0], node->bbox[2], node->bbox[4] };
  double X1[] = { node->bbox[1], node->bbox[3], node->bbox[5] };
  double hmin;

  /* set hmin from node diagonal */
  hmin = Cart_distance_X0_X1(node, X0,X1)/sqrt(3.);

  /* first corner */
  hmin = distance_to_closest_point(node, 0,0,0, hmin, ijk0, ijk1);

  /* next corner */
  hmin = distance_to_closest_point(node, n[0]-1,0,0, hmin, ijk0, ijk1);

  /* next corner */
  hmin = distance_to_closest_point(node, 0,n[1]-1,0, hmin, ijk0, ijk1);

  /* next corner */
  hmin = distance_to_closest_point(node, n[0]-1,n[1]-1,0, hmin, ijk0, ijk1);

  /* next corner */
  hmin = distance_to_closest_point(node, 0,0,n[2]-1, hmin, ijk0, ijk1);

  /* next corner */
  hmin = distance_to_closest_point(node, n[0]-1,0,n[2]-1, hmin, ijk0, ijk1);

  /* next corner */
  hmin = distance_to_closest_point(node, 0,n[1]-1,n[2]-1, hmin, ijk0, ijk1);

  /* last corner */
  hmin = distance_to_closest_point(node, n[0]-1,n[1]-1,n[2]-1, hmin,
                                   ijk0, ijk1);

  //PRF;printf(": pts=%d,%d -> hmin=%g\n", *ijk0,*ijk1, hmin);
  return hmin;
}

/* Write opposite corner into ijk2 and find distance dist02 to it.
   ijk2 is found by going from ijk0 through ijk1 until we get to a corner.
   return 1 if ijk2 was found and set
   return 0 if ijk0 or ijk1 are negative which is an error */
int distance_to_opposite_corner_ijk2(tNode *node, int *ijk0, int *ijk1,
                                     int *ijk2, double *dist02)
{
  int *n = node->n;

  if(*ijk0>=0 && *ijk1>=0)
  {
    int i0,j0,k0, i1,j1,k1, di,dj,dk;
    int i2,j2,k2;
    double X0[3], X2[3];
    /* get i,j,k vals of both points */
    k0 = kOfInd_n(*ijk0, n);
    j0 = jOfInd_n_k(*ijk0, n, k0);
    i0 = iOfInd_n_jk(*ijk0, n, j0,k0);
    k1 = kOfInd_n(*ijk1, n);
    j1 = jOfInd_n_k(*ijk1, n, k1);
    i1 = iOfInd_n_jk(*ijk1, n, j1,k1);
    di = (i1 - i0);
    dj = (j1 - j0);
    dk = (k1 - k0);
    /* construct index i2,j2,k2 of furthest point along di,dj,dk */
    i2 = i0 + di*(n[0]-1);
    j2 = j0 + dj*(n[1]-1);
    k2 = k0 + dk*(n[2]-1);
    *ijk2 = Ind_n(i2,j2,k2, n);
    /* get X-coords of points at i0,j0,k0 and i2,j2,k2 */
    XYZ_of_ijk(node, i0,j0,k0, X0);
    XYZ_of_ijk(node, i2,j2,k2, X2);
    /* set distance between X0 and X2 */
    *dist02 = Cart_distance_X0_X1(node, X0,X2);
    return 1;
  }
  return 0;
}

/* Set the max possible 1D timestep dtmax in all 3 directions separately.
   This func usually takes ijk0 from find_hmin as input.
   The CFL fac dist/(2*n[0] - 1) is from page 16 of:
   Galerkin methods. In Discontinuous Galerkin Methods, pages 3-50. Springer,
   2000. Bernardo Cockburn and Chi-Wang Shu. Runge-Kutta discontinuous
   Galerkin methods for convection-dominated problems. Journal of Scientific
   Computing, 16(3):173–261, 2001.
   nmesh-extra/resources/DG-Methods/Runge-Kutta_DG-Methods.pdf */
int set_dtmax3_from_corner_ijk0(tNode *node, int *ijk0, double dtmax[3])
{
  int *n = node->n;

  if(*ijk0>=0)
  {
    double Xb_ijk0[3], X_ijk0[3];
    double Xb_d[3];
    int dir;

    /* set Xb and X of ijk0 */
    XbYbZb_of_ind(node, *ijk0, Xb_ijk0);
    XYZ_of_XbYbZb(node, Xb_ijk0, X_ijk0);

    /* get Xb_d[dir] by moving a bit in for each dir */
    for(dir=0; dir<3; dir++)
    {
      /* we assume Xb_ijk0 is on node face */
      if(Xb_ijk0[dir]>0.) Xb_d[dir] = +0.99;
      else                Xb_d[dir] = -0.99;
    }

    /* set dt_max from dist */
    for(dir=0; dir<3; dir++)
    {
      double Xb_in[3];
      double X[3], dist, dx_o_dXb, L;
      int d;

      /* set Xb-coords of points that are a bit in */
      for(d=0; d<3; d++) Xb_in[d] = Xb_ijk0[d];
      Xb_in[dir] = Xb_d[dir];

      /* get dist between point a bit in and X_ijk0 */
      XYZ_of_XbYbZb(node, Xb_in, X); // set X coords of point a bit in
      dist  = Cart_distance_X0_X1(node, X_ijk0, X);

      dx_o_dXb = dist/0.01; //ratio between dist in x and Xb coords
      L = dx_o_dXb*2.;      //length of node assuming local dx_o_dXb is const

      if(node->pt_typ[dir]==P_UNIFORM)
        dtmax[dir] = L / (n[dir] - 1);
      else
        dtmax[dir] = L / (2*n[dir] - 1);
    }
    return 1;
  }
  return 0;
}

/* check if we should reduce dt because hmin is only half of what
   find_hmin finds if we use fin.vol. (FV) */
int hmin_is_in_uniform_direction(tNode *node, int *ijk0, int *ijk1)
{
  int *n = node->n;

  if(*ijk0>=0 && *ijk1>=0)
  {
    int i0,j0,k0, i1,j1,k1, di,dj,dk;
    int dir = -1;
    /* get i,j,k vals of both points */
    k0 = kOfInd_n(*ijk0, n);
    j0 = jOfInd_n_k(*ijk0, n, k0);
    i0 = iOfInd_n_jk(*ijk0, n, j0,k0);
    k1 = kOfInd_n(*ijk1, n);
    j1 = jOfInd_n_k(*ijk1, n, k1);
    i1 = iOfInd_n_jk(*ijk1, n, j1,k1);
    di = abs(i1 - i0);
    dj = abs(j1 - j0);
    dk = abs(k1 - k0);
    if(di==1 && dj==0 && dk==0) dir = 0;
    if(di==0 && dj==1 && dk==0) dir = 1;
    if(di==0 && dj==0 && dk==1) dir = 2;
    if(dir>=0)
    {
      if(node->pt_typ[dir] == P_UNIFORM)
        return 1; /* FV cell near face has only half the usual length */
    }
    else
    {
      return 0;
    }
  }
  return 0;
}


/* change node->dt and mesh->dt based on the smallest grid spacing hmin,
   return this dt  */
double adapt_node_dt_and_mesh_dt(tNode *node, double dtfac,
                                 double uniform_dtfac)
{
  tMesh *mesh = node->pat->mesh;
  double dtm, hmin;
  int ijk0[] = {-1}, ijk1[] = {-1};

  if(mesh->dt < node->dt || node->dt <= 0.) node->dt = mesh->dt;
  hmin = find_hmin(node, ijk0,ijk1);
  dtm = dtfac * hmin;
  /* effectively hmin may be reduced by half if we use fin.vol. (FV) */
  if(hmin_is_in_uniform_direction(node, ijk0,ijk1))
    dtm = uniform_dtfac * hmin;

  if(dtm < node->dt || node->dt <= 0.)
  {
    node->dt = dtm;
    mesh->dt = dtm;
    PRFs(": ");pr_nodename(node);
    printf(" pts = %d,%d: hmin = %g\n", ijk0[0],ijk1[0], hmin);
    PRF;printf(": setting mesh->dt = %g\n", mesh->dt);
  }
  return dtm;
}


/* put approx. node center in X-coords into X */
void set_nodecenter_XYZ(tNode *node, double X[3])
{
  double *bb = node->bbox;
  int d;

  for(d=0; d<3; d++)
    X[d] = 0.5*(bb[2*d + 1] + bb[2*d]);
}

/* put approx. node center in Cart. coords into x */
void set_nodecenter_xyz(tNode *node, double x[3])
{
  double X[3];

  /* first get X of center */
  set_nodecenter_XYZ(node, X);

  /* now get xc from X */
  set_xyz(NULL, node, -1, X, x);
}
