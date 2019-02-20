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

/* find i,j,k closest to Xb0 */
void nearest_ijk_of_XbYbZb(tNode *node, int ijk[3], const double Xb0[3])
{
  int *n = node->n;
  double *Xb[] = { node->Xb[0]->d, node->Xb[1]->d, node->Xb[2]->d };
  int dir, i,j,k;
  double a,b;

  for(dir=0; dir<3; dir++)
  {
    i = 0;
    k = n[dir]-1;
    a = Xb[dir][i];
    b = Xb[dir][k];
    if(dgreater(Xb0[dir], b))
    {
      ijk[dir] = -k-1;
      continue;
    }
    if(dless(Xb0[dir], a))
    {
      ijk[dir] = -1;
      continue;
    }

    while(k-i>1)
    {
      j=(i+k)/2;

      a = Xb[dir][i] - Xb0[dir];
      b = Xb[dir][j] - Xb0[dir];
      if(a*b<=0.0) k=j;
      else         i=j;
    }
    if( fabs(Xb[dir][i] - Xb0[dir]) < fabs(Xb[dir][k] - Xb0[dir]) )
      ijk[dir] = i;
    else
      ijk[dir] = k;
  }
}

/* find i,j,k closest to X0 */
void nearest_ijk_of_XYZ(tNode *node, int ijk[3], const double X0[3])
{
  double Xb[3];
  XbYbZb_of_XYZ(node, Xb, X0);
  nearest_ijk_of_XbYbZb(node, ijk, Xb);
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
