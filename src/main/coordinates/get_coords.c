/* get_coords.c */
/* Wolfgang Tichy, 2/2019 */

#include "nmesh.h"
#include "coordinates.h"

#define PR 1



/* get Xb from i,j,k */
void XbYbZb_of_ijk(tNode *node, int i, int j, int k, double Xb[3])
{
  tArray *A[] = { node->Xb[0], node->Xb[1], node->Xb[2] };
  int m[] = { i,j,k };
  int dir;
  
  for(dir=0; dir<3; dir++) Xb[dir] = A[dir]->d[m[dir]];
}
/* get Xb from index ind */ 
void XbYbZb_of_ind(tNode *node, int ind, double Xb[3])
{
  int *n = node->n;
  int k = kOfInd_n(ind, n);
  int j = jOfInd_n_k(ind, n, k);
  int i = iOfInd_n_jk(ind, n, j,k);

  XbYbZb_of_ijk(node, i,j,k, Xb);
}

/* get X from Xb */
void XYZ_of_XbYbZb(tNode *node, const double Xb[3], double X[3])
{
  double *nbb = node->bbox;
  int dir;

  for(dir=0; dir<3; dir++)
  {
    int f = dir*2;
    X[dir] = 0.5*( (nbb[f+1] - nbb[f]) * Xb[dir] + (nbb[f+1] + nbb[f]) );
  } 
}

/* get Xb from X */
void XbYbZb_of_XYZ(tNode *node, double Xb[3], const double X[3])
{
  double *nbb = node->bbox;
  int dir;

  for(dir=0; dir<3; dir++)
  {
    int f = dir*2;
    Xb[dir] = ( 2.*X[dir] - (nbb[f+1] + nbb[f]) )/(nbb[f+1] - nbb[f]);
  } 
}
