/* get_coords.c */
/* Wolfgang Tichy, 2/2019 */

#include "nmesh.h"
#include "coordinates.h"



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

/* get one X from Xb in one direction */
void X_of_Xb_indir(tNode *node, int dir, double Xb, double *X)
{
  double *nbb = node->bbox;
  int f = dir*2;
  *X = 0.5*( (nbb[f+1] - nbb[f]) * Xb + (nbb[f+1] + nbb[f]) );
}

/* get X from Xb for entire arrays */
void array_XYZ_of_XbYbZb(tNode *node, tArray *aXb[3], tArray *aX[3])
{
  int dir, k;
  for(dir=0; dir<3; dir++)
  {
    forarray(aXb[dir], k)
    {
      double Xb, X;
      Xb = aXb[dir]->d[k];
      X_of_Xb_indir(node, dir, Xb, &X);
      aX[dir]->d[k] = X;
    }
  }
}

/* get X from i,j,k */
void XYZ_of_ijk(tNode *node, int i, int j, int k, double X[3])
{
  XbYbZb_of_ijk(node, i,j,k, X);
  XYZ_of_XbYbZb(node, X, X);
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

/* get one Xb from X in one direction */
void Xb_of_X_indir(tNode *node, int dir, double *Xb, const double X)
{
  double *nbb = node->bbox;
  int f = dir*2;
  *Xb = ( 2.*X - (nbb[f+1] + nbb[f]) )/(nbb[f+1] - nbb[f]);
}

/* get Xb from X for entire arrays */
void array_XbYbZb_of_XYZ(tNode *node, tArray *aXb[3], tArray *aX[3])
{
  int dir, k;
  for(dir=0; dir<3; dir++)
  {
    forarray(aXb[dir], k)
    {
      double Xb, X;
      X = aX[dir]->d[k];
      Xb_of_X_indir(node, dir, &Xb, X);
      aXb[dir]->d[k] = Xb;
    }
  }
}

/* is XYZ inside a node, also sets X to boundary value if it is very close
   to the boundary */
int XYZ_is_in_node(tNode *node, double X[3])
{
  double *nbb = node->bbox;
  int dir;

  /* return 0 if X is outside node */
  for(dir=0; dir<3; dir++)
  {
    int f = dir*2;
    if(dless(X[dir],nbb[f]))      return 0;
    if(dgreater(X[dir],nbb[f+1])) return 0;
  }

  /* set X to boundary value if it is very close to boundary */
  for(dir=0; dir<3; dir++)
  {
    int f = dir*2;
    if(X[dir] < nbb[f])   X[dir] = nbb[f];
    if(X[dir] > nbb[f+1]) X[dir] = nbb[f+1];
  }
  return 1;
}

/* write all points in aXP[0..1] within node into aX[0..2]
   NOTE: this re-dimensions aX[0..2] !!! */
void array_get_XYZ_in_node(tNode *node, tArray *aXP[3], tArray *aX[3])
{
  int dir, k, ai;

  ai = 0;
  forarray(aXP[0], k)
  {
    double X[] = { aXP[0]->d[k], aXP[1]->d[k], aXP[2]->d[k] };

    if(XYZ_is_in_node(node, X))
    {
      for(dir=0; dir<3; dir++) aX[dir]->d[ai] = X[dir];
      ai++;
    }
  }

  /* redimension array aX to have correct number of points */
  for(dir=0; dir<3; dir++) redim_array(aX[dir], ai,1,1);
}
