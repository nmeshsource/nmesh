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

/* get dX/dXb */
void dXYZ_dXbYbZb(tNode *node, double dXdXb[3])
{
  double *nbb = node->bbox;
  int dir;

  for(dir=0; dir<3; dir++)
  {
    int f = dir*2;
    dXdXb[dir] = 0.5*( nbb[f+1] - nbb[f] );
  } 
}
/* get dX/dXb */
void dXbYbZb_dXYZ(tNode *node, double dXbdX[3])
{
  double *nbb = node->bbox;
  int dir;

  for(dir=0; dir<3; dir++)
  {
    int f = dir*2;
    dXbdX[dir] = 2./( nbb[f+1] - nbb[f] );
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
    int Nm = min2(aXb[dir]->N, aX[dir]->N);
    for(k=0; k<Nm; k++)
    {
      double Xb, X;
      Xb = aXb[dir]->d[k];
      X_of_Xb_indir(node, dir, Xb, &X);
      aX[dir]->d[k] = X;
    }
  }
}

/* get X from Xb for arrays in plane perp to dir */
void array_Xplane_of_Xb(tNode *node, int dir, tArray *aCb[2], tArray *aC[2])
{
  int d, d3, k, Nm;
  for(d=0; d<2; d++)
  {
    switch(dir)
    {
    case 0:
      d3 = d+1;
      break;
    case 1:
      d3 = d*2;
      break;
    case 2:
      d3 = d;
      break;
    default:
      errorexit("dir must be 0,1,2");
    }
    Nm = min2(aCb[d]->N, aC[d]->N);
    for(k=0; k<Nm; k++)
    {
      double Xb, X;
      Xb = aCb[d]->d[k];
      X_of_Xb_indir(node, d3, Xb, &X);
      aC[d]->d[k] = X;
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
    int Nm = min2(aXb[dir]->N, aX[dir]->N);
    for(k=0; k<Nm; k++)
    {
      double Xb, X;
      X = aX[dir]->d[k];
      Xb_of_X_indir(node, dir, &Xb, X);
      aXb[dir]->d[k] = Xb;
    }
  }
}

/* get Xb from X for arrays in plane perp to dir */
void array_Xbplane_of_X(tNode *node, int dir, tArray *aCb[2], tArray *aC[2])
{
  int d, d3, k, Nm;
  for(d=0; d<2; d++)
  {
    switch(dir)
    {
    case 0:
      d3 = d+1;
      break;
    case 1:
      d3 = d*2;
      break;
    case 2:
      d3 = d;
      break;
    default:
      errorexit("dir must be 0,1,2");
    }
    Nm = min2(aCb[d]->N, aC[d]->N);
    for(k=0; k<Nm; k++)
    {
      double Xb, X;
      X = aC[d]->d[k];
      Xb_of_X_indir(node, d3, &Xb, X);
      aCb[d]->d[k] = Xb;
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

/* Mark all points in aXP[0..2] within node by writing their index into
   aI. If a point is not in the node we write -1 into aI. I.e. we return
   a mask of points in node in Ip. */
void array_find_XYZ_in_node(tNode *node, tArray *aXP[3], tArray *aI)
{
  int k;
  forarray(aXP[0], k)
  {
    double X[] = { aXP[0]->d[k], aXP[1]->d[k], aXP[2]->d[k] };
    if(XYZ_is_in_node(node, X)) aI->i[k] = k;
    else                        aI->i[k] = -1;
  }
}


/* is a point C (in XYZ coords) in a plane (normal to dir) inside a node,
   also sets C to boundary value if it is very close to the boundary */
int Xplane_is_in_node(tNode *node, int dir, double C[2])
{
  double *nbb = node->bbox;
  int d, f;

  /* return 0 if C is outside node */
  for(d=0; d<2; d++)
  {
    switch(dir)
    {
    case 0:
      f = 2*d + 2;
      break;
    case 1:
      f = 4*d;
      break;
    case 2:
      f = 2*d;
      break;
    default:
      errorexit("dir has to be 0,1,2");
    }
    if(dless(C[d],nbb[f]))      return 0;
    if(dgreater(C[d],nbb[f+1])) return 0;
  }

  /* set X to boundary value if it is very close to boundary */
  for(d=0; d<2; d++)
  {
    switch(dir)
    {
    case 0:
      f = 2*d + 2;
      break;
    case 1:
      f = 4*d;
      break;
    case 2:
      f = 2*d;
    }
    if(C[d] < nbb[f])   C[d] = nbb[f];
    if(C[d] > nbb[f+1]) C[d] = nbb[f+1];
  }
  return 1;
}

/* Mark all points in aCP[0..1] in plane p normal to dir within node by
   writing their index into aI. If a point is not in the node we write -1
   into aI. I.e. we return a mask of points in node in Ip. */
void array_find_Xplane_in_node(tNode *node,int dir, tArray *aCP[2], tArray *aI)
{
  int k;
  forarray(aCP[0], k)
  {
    
    double C[]  = { aCP[0]->d[k], aCP[1]->d[k] };
    if(Xplane_is_in_node(node, dir, C)) aI->i[k] = k;
    else                                aI->i[k] = -1;
  }
}
