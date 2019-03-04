/* utilities.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "coordinates.h"



/* find normal vector (n[0],n[1],n[2]) of patch face f at X,Y,Z */
void patface_normal_at_XYZ(tPat *pat, int f, const double X[3],
                           double n[3])
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
    n[1] = n[2] = n[3] = 0.0;
    n[dir] = 1.0;
  }

  /* normalize and set sign from sig */
  smag = sig * sqrt(n[1]*n[1] + n[2]*n[2] + n[3]*n[3]);
  if(smag == 0.0) smag = sig;
  n[1] /= smag;
  n[2] /= smag;
  n[3] /= smag;
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
    pat->xyz_of_XYZ(pat, NULL,-1, X[i], x[i]);

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
