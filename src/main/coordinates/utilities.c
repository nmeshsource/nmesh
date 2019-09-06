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

/* find normal vector (n[0],n[1],n[2]) of pat face f at point ijk */
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

/* get transpose of a 3x3 matrix. */
void transp3Dmat_from_3Dmat(double M[3][3], double transpM[3][3])
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
double det_3Dmatrix(double M[3][3])
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
double inv3Dmat_from_3Dmat(double M[3][3], double invM[3][3])
{
  /* M = {{m11,m12,m13},{m21,m22,m23},{m31,m32,m33}}
     Det[M] = m11*m22*m33 - m13*m22*m31 + m12*m23*m31 +
              m13*m21*m32 - m11*m23*m32 - m12*m21*m33
     Inverse[M]*Det[M] =
       {{-(m23*m32) + m22*m33, m13*m32 - m12*m33, -(m13*m22) + m12*m23},
        {m23*m31 - m21*m33, -(m13*m31) + m11*m33, m13*m21 - m11*m23},
        {-(m22*m31) + m21*m32, m12*m31 - m11*m32, -(m12*m21) + m11*m22}}  */
  LDOUBLE DetM = det_3Dmatrix(M);
  LDOUBLE sum;
  if(DetM==0.0) DetM=dequaleps*dequaleps*dequaleps;
  sum        = (-(M[1][2]*M[2][1]) + M[1][1]*M[2][2])/DetM;
  invM[0][0] = sum;
  sum        = (  M[0][2]*M[2][1]  - M[0][1]*M[2][2])/DetM;
  invM[0][1] = sum;
  sum        = (-(M[0][2]*M[1][1]) + M[0][1]*M[1][2])/DetM;
  invM[0][2] = sum;
  sum        = (  M[1][2]*M[2][0]  - M[1][0]*M[2][2])/DetM;
  invM[1][0] = sum;
  sum        = (-(M[0][2]*M[2][0]) + M[0][0]*M[2][2])/DetM;
  invM[1][1] = sum;
  sum        = (  M[0][2]*M[1][0]  - M[0][0]*M[1][2])/DetM;
  invM[1][2] = sum;
  sum        = (-(M[1][1]*M[2][0]) + M[1][0]*M[2][1])/DetM;
  invM[2][0] = sum;
  sum        = (  M[0][1]*M[2][0]  - M[0][0]*M[2][1])/DetM;
  invM[2][1] = sum;
  sum        = (-(M[0][1]*M[1][0]) + M[0][0]*M[1][1])/DetM;
  invM[2][2] = sum;

  return DetM;
}

/* print 3x3 matrix */
void print3Dmat(double M[3][3])
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

  return det;
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
  double L=1e300;
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
void array_3Dmat(double M[3][3], tArray *aM)
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
