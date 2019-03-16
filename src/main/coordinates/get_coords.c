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
    /* set a, b */
    Xb_of_X_indir(node, dir, &a, node->bbox[2*dir]);
    Xb_of_X_indir(node, dir, &b, node->bbox[2*dir + 1]);
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


/* find index of nb node in the node->fnb[f] list that contains
   the surface point C on face o_f of patch o_pat,
   return index of nb in fnb, or -1 if not found */
int fnb_containing_point(tNode *node, int f,
                         tPat *o_pat, int o_f, double C[2])
{
  int odir = o_f/2;
  int od1 = Dir1_norm(odir);
  int od2 = Dir2_norm(odir);
  int nfnb = node->nfnb[f];
  int i;

  for(i=0; i<nfnb; i++)
  {
    tNode *nb = node->fnb[f][i];
    double nbrct[] = { nb->bbox[2*od1], nb->bbox[2*od1+1],
                       nb->bbox[2*od2], nb->bbox[2*od2+1] };

    if(nb->pat != o_pat) continue;

    if(C_in_brct(nbrct , C)) return i;
  }
  return -1;
}


/* Mark all points in aC[0..1] that are in node->fnb[f][ni] by writing
   their index into aI. If a point is not in the node we write -1
   into aI. I.e. we return a mask of points in node in Ip.
   Note that aC and aoC contain the same points but in different coords.
   aC in X coords of the node and aoC in X coords of the neighbor. */
void brct_mark_points_in_fnb_f_ni(tNode *node, int f, int ni, tArray *aC[2],
                                  tArray *aoC[2], tArray *aI)
{
  tPat *pat = node->pat;
  tNode *nb = node->fnb[f][ni];
  tBface *obface;
  tPat *opat;
  int odir, od1, od2;
  double nbrct[4];
  int k;

  forarray(aC[0], k)
  {
    double C[]  = {  aC[0]->d[k],  aC[1]->d[k] };
    double oC[] = { aoC[0]->d[k], aoC[1]->d[k] };

    /* find bface on other side with C */
    obface = nbbface_of_bface_containing_point(nb, pat, f, C);
    if(!obface)
    {
      aI->i[k] = -1;
      continue;
    }

    /* get patch on other side */
    opat = obface->pat;
    if(!opat)
    {
      aI->i[k] = -1;
      continue;
    }

    /* get bounding rectangle of the node nb */
    odir = obface->f/2;
    od1 = Dir1_norm(odir);
    od2 = Dir2_norm(odir);
    nbrct[0] = nb->bbox[2*od1];
    nbrct[1] = nb->bbox[2*od1+1];
    nbrct[2] = nb->bbox[2*od2];
    nbrct[3] = nb->bbox[2*od2+1];

    if(0) //(f==1 && node->nid==1)
    {
      char s[100], snb[100];
      printf("ni=%d: %s  %s ", ni, nodename(node,s,99), nodename(nb,snb,99));
      prbbox(nbrct,2);printf("\n");
      printf("nb->pat=%d opat=%d  oC=%g %g  in=%d\n",
      nb->pat->p, opat->p, oC[0],oC[1], C_in_brct(nbrct , oC));
      //printnode(node);
      //printbfaces_on_f(node->pat, f);
      //printnode(nb);
      //printbfaces_on_f(nb->pat, nb_f);
    }

    /* check if C from opat is within bounding rectangle of nb */
    if( (nb->pat == opat) && (C_in_brct(nbrct , oC)) )
      aI->i[k] = k;
    else
      aI->i[k] = -1;
  }
}

/* Mark all points in aC[0..1] of node that are in node nb by writing
   their index into aI. If a point is not in the node we write -1
   into aI. I.e. we return a mask of points in node in Ip.
   Note that aC and aoC contain the same points but in different coords.
   aC in X coords of the node and aoC in X coords of the neighbor nb. */
void mark_points_in_nb_f(tNode *node, int f, tArray *aC[2],
                         tNode *nb, int nb_f, tArray *aoC[2], tArray *aI)
{
  //tPat *pat = node->pat;
  //tPat *opat = nb->pat;
  int odir, od1, od2;
  double nbrct[4];
  int k;

  /* get bounding rectangle of the node nb */
  odir = nb_f/2;
  od1 = Dir1_norm(odir);
  od2 = Dir2_norm(odir);
  nbrct[0] = nb->bbox[2*od1];
  nbrct[1] = nb->bbox[2*od1+1];
  nbrct[2] = nb->bbox[2*od2];
  nbrct[3] = nb->bbox[2*od2+1];

  forarray(aC[0], k)
  {
    //double C[]  = {  aC[0]->d[k],  aC[1]->d[k] };
    double oC[] = { aoC[0]->d[k], aoC[1]->d[k] };

    /* check if oC from opat is within bounding rectangle of nb */
    if( C_in_brct(nbrct, oC) )
      aI->i[k] = k;
    else
      aI->i[k] = -1;
  }
}


/* return patch index if x is inside this patch, if not return -1 */
int p_XYZ_of_xyz(tPat *pat, double X[3], const double x[3])
{
  int d, stat=0;

  /* get X */
  if(pat->XYZ_of_xyz)
    stat = pat->XYZ_of_xyz(pat, 0,-1, X, x);
  else
    for(d=0; d<3; d++) X[d] = x[d];

  //pr3v("|x",x);
  //pr3v("X",X);

  if(stat) return -1;

  for(d=0; d<3; d++)
    if(dless(X[d],pat->bbox[2*d]) || dless(pat->bbox[2*d+1],X[d]))
      return -1;

  //prbbox(pat->bbox,3);

  /* round X to inside box */
  for(d=0; d<3; d++)
  {
    if(X[d] < pat->bbox[2*d])   X[d] = pat->bbox[2*d];
    if(X[d] > pat->bbox[2*d+1]) X[d] = pat->bbox[2*d+1];
  }
  return pat->p;
}

/* go over pat list and find the one that contains x */
int p_XYZ_of_xyz_inpatlist(tMesh *mesh, intList *pl,
                           double X[3], const double x[3])
{
  tPat *pat;
  int i, p=-1;

  for(i=0; i<pl->n; i++)
  {
    p = pl->e[i];
    pat = mesh->pat[p];
    p = p_XYZ_of_xyz(pat, X, x);
    if(p>0) break;
  }
  return p;
}

/* return node location if x is inside this patch, if not return -1 */
long l_XYZ_of_xyz(tNode *node, int ind, double X[3], const double x[3])
{
  tPat *pat = node->pat;
  int d, stat=0;
  long loc = node_location(node); /* get node location */

  /* get X */
  if(pat->XYZ_of_xyz)
    stat = pat->XYZ_of_xyz(pat, node,ind, X, x);
  else
    for(d=0; d<3; d++) X[d] = x[d];

  if(stat) return -1;

  for(d=0; d<3; d++)
    if(dless(X[d],node->bbox[2*d]) || dless(node->bbox[2*d+1],X[d]))
      return -1;

  /* round X to inside box */
  for(d=0; d<3; d++)
  {
    if(X[d] < node->bbox[2*d])   X[d] = node->bbox[2*d];
    if(X[d] > node->bbox[2*d+1]) X[d] = node->bbox[2*d+1];
  }

  return loc;
}


/* find the faces a point X is on within tol, face[2]=1 if X is on face2  */
int XYZ_on_face_tol(tPat *pat, int *face, const double X[3], double tol)
{
  double *bb=pat->bbox;
  int f;
  int nf;

  /* find all faces we are on */
  for(nf=0, f=0; f<6; f++)
  {
    int d=f/2;
    if(dequal_tol(X[d], bb[f], tol)) { face[f]=1; nf++; }
    else                             { face[f]=0; }
  }
  return nf; /* number of faces point is on */
}

/* find the faces a point X is on, face[2]=1 if X is on face2  */
int XYZ_on_face(tPat *pat, int *face, const double X[3])
{
  return XYZ_on_face_tol(pat, face, X, 1e-10);
}


/* set x at X */
int set_xyz(tPat *pat, tNode *node, int ind, const double X[3], double x[3])
{
  if(!pat) pat = node ? node->pat : NULL;

  /* now set x, dXb/dx */
  if(pat->xyz_of_XYZ)
  {
    return pat->xyz_of_XYZ(pat, node, ind, X, x);
  }
  else /* assume X,Y,Z are Cartesian*/
  {
    int d;

    for(d=0; d<3; d++)
      x[d] = X[d];
    return 0;
  }
}

/* set X at x */
int set_XYZ(tPat *pat, tNode *node, int ind, double X[3], const double x[3])
{
  if(!pat) pat = node ? node->pat : NULL;

  /* now set x, dXb/dx */
  if(pat->XYZ_of_xyz)
  {
    return pat->XYZ_of_xyz(pat, node, ind, X, x);
  }
  else /* assume X,Y,Z are Cartesian*/
  {
    int d;

    for(d=0; d<3; d++)
      X[d] = x[d];
    return 0;
  }
}

/* set x and dXYZ/dxyz at X */
int set_xyz_dXYZdxyz(tPat *pat, tNode *node, int ind,
                     const double X[3], double x[3], double dXYZdxyz[3][3])
{
  /* now set x, dXdx, det(dXb/dx) */
  if(pat->dXYZ_dxyz)
  {
    return pat->dXYZ_dxyz(pat, node, ind, X, x, dXYZdxyz);
  }
  else /* assume X,Y,Z are Cartesian*/
  {
    int d, e;

    for(d=0; d<3; d++)
    {
      x[d] = X[d];
      for(e=0; e<3; e++)
        dXYZdxyz[d][e] = (d==e);
    }
    return 0;
  }
}

/* get the bounding rectangle of a nodeface, norm = face/2 */
void brct_nodeface(tNode *node, int norm, double brct[4])
{
  int i;
  switch(norm)
  {
  case 0:
    for(i=0; i<4; i++) brct[i] = node->bbox[i+2];
    break;
  case 1:
    for(i=0; i<2; i++) brct[i] = node->bbox[i];
    for(i=2; i<4; i++) brct[i] = node->bbox[i+2];
    break;
  case 2:
    for(i=0; i<4; i++) brct[i] = node->bbox[i];
    break;
  default:
    errorexit("norm must be 0,1,2");
  }
}

/* get the bounding rectangle of a patchface, norm = face/2 */
void brct_patface(tPat *pat, int norm, double brct[4])
{
  int i;
  switch(norm)
  {
  case 0:
    for(i=0; i<4; i++) brct[i] = pat->bbox[i+2];
    break;
  case 1:
    for(i=0; i<2; i++) brct[i] = pat->bbox[i];
    for(i=2; i<4; i++) brct[i] = pat->bbox[i+2];
    break;
  case 2:
    for(i=0; i<4; i++) brct[i] = pat->bbox[i];
    break;
  default:
    errorexit("norm must be 0,1,2");
  }
}

/* make bounding rectangle large enough to fit the point X[3],
   if expand=0, include just this one point */
void expand_brct_to_include_X(double brct[4], int norm,
                              const double X[3], int expand)
{
  int d;
  double C[2]; /* point coords in face */

  switch(norm)
  {
  case 0:
    C[0] = X[1];  C[1] = X[2];
    break;
  case 1:
    C[0] = X[0];  C[1] = X[2];
    break;
  case 2:
    C[0] = X[0];  C[1] = X[1];
    break;
  default:
    errorexit("norm must be 0,1,2");
  }
  if(expand)
  {
    /* expand rectangle */
    for(d=0; d<2; d++)
    {
      if(C[d] < brct[2*d])   brct[2*d]   = C[d];
      if(C[d] > brct[2*d+1]) brct[2*d+1] = C[d];
    }
  }
  else
  {
    brct[1] = brct[0] = C[0];
    brct[3] = brct[2] = C[1];
   }
}

/* put intersection of 2 bounding rectangles into brct, if intersection
   is empty return 0 */
int intersection_brct1_brct2(const double brct1[4], const double brct2[4],
                             double brct[4])
{
  int d, isec=1;

  for(d=0; d<2; d++)
  {
    if((brct1[2*d] >= brct2[2*d]) && (brct1[2*d] <= brct2[2*d+1]))
    {
      brct[2*d] = brct1[2*d];
      if(brct1[2*d+1] < brct2[2*d+1])
        brct[2*d+1] = brct1[2*d+1];
      else
        brct[2*d+1] = brct2[2*d+1];
    }
    else if((brct2[2*d] >= brct1[2*d]) && (brct2[2*d] <= brct1[2*d+1]))
    {
      brct[2*d] = brct2[2*d];
      if(brct2[2*d+1] < brct1[2*d+1])
        brct[2*d+1] = brct2[2*d+1];
      else
        brct[2*d+1] = brct1[2*d+1];
    }
    else
    {
      isec = 0;
      break;
    }
  }

  /* test if the intersection is empty */
  if(isec)
    for(d=0; d<2; d++)
      if(dequal(brct[2*d], brct[2*d+1]))
      {
        isec = 0;
        break;
      }

  return isec;
}

/* transform XYZ from patch1 to XYZ of other patch2 */
int XYZpat2_of_XYZpat1(tPat *pat1, const double X1[3],
                        tPat *pat2, double X2[3])
{
  double x[3];
  int p2;

  set_xyz(pat1, NULL, -1, X1, x);
  p2 = p_XYZ_of_xyz(pat2, X2, x);
  return p2;
}

/* get Coords on face from X */
void C_from_X_on_face(const double X[3], int face, double C[2])
{
  switch(face/2)
  {
  case 0:
    C[0] = X[1];
    C[1] = X[2];
    break;
  case 1:
    C[0] = X[0];
    C[1] = X[2];
    break;
  case 2:
    C[0] = X[0];
    C[1] = X[1];
    break;
  default:
    errorexit("face/2 must be 0,1,2");
  }
}

/* get X from Coords on face */
void X_from_C_on_face(tPat *pat, int face, const double C[2], double X[3])
{
  switch(face/2)
  {
  case 0:
    X[0] = pat->bbox[face];
    X[1] = C[0];
    X[2] = C[1];
    break;
  case 1:
    X[0] = C[0];
    X[1] = pat->bbox[face];
    X[2] = C[1];
    break;
  case 2:
    X[0] = C[0];
    X[1] = C[1];
    X[2] = pat->bbox[face];
    break;
  default:
    errorexit("face/2 must be 0,1,2");
  }
}

/* transform point C1 on face f1 of pat1 to point C2 on face2 of pat2 */
int Cpat2_of_Cpat1(tPat *pat1, int f1, const double C1[2],
                    tPat *pat2, int f2, double C2[2])
{
  double X1[3], X2[3];
  int p2;

  X_from_C_on_face(pat1,f1, C1, X1);
  p2 = XYZpat2_of_XYZpat1(pat1, X1, pat2, X2);
  /* if(p2<0) X2 is not within pat2 or NAN */

  C_from_X_on_face(X2,f2, C2);

  return p2;
}

/* transform brct from face f1 of patch1 to brct of face f2 on patch2 */
int brctpat2_of_brctpat1__old(tPat *pat1, int f1, const double brct1[4],
                         tPat *pat2, int f2, double brct2[4])
{
  double C1[2], C2[2], sw;
  int p2_1, p2_2, problem = 0;

  C1[0] = brct1[0];
  C1[1] = brct1[2];
  p2_1 = Cpat2_of_Cpat1(pat1,f1,C1,   pat2,f2,C2);
  /* set lower bounds of brct2 */
  brct2[0] = C2[0];
  brct2[2] = C2[1];


  C1[0] = brct1[1];
  C1[1] = brct1[3];
  p2_2 = Cpat2_of_Cpat1(pat1,f1,C1, pat2,f2,C2);
  /* set upper bounds of brct2 */
  brct2[1] = C2[0];
  brct2[3] = C2[1];

  /* check if C2 is not NAN */
  if(p2_1<0 || p2_2<0)
  {
    int j;

    /* check if all is ok, or if there is a NAN */
    for(j=0; j<4; j++) if(isnan(brct2[j])) problem = 1;
  }

  /* swap brct2 entries to have min in brct2[0], brct2[2] */
  if(brct2[1] < brct2[0])
  {
    sw = brct2[1];
    brct2[1] = brct2[0];
    brct2[0] = sw;
  }
  if(brct2[3] < brct2[2])
  {
    sw = brct2[3];
    brct2[3] = brct2[2];
    brct2[2] = sw;
  }
  return problem;
}

/* transform brct from face f1 of patch1 to brct of face f2 on patch2 */
int brctpat2_of_brctpat1(tPat *pat1, int f1, const double brct1[4],
                         tPat *pat2, int f2, double brct2[4])
{
  double C1[4][2]; /* 4 points, 3rd one is (C1[2][0],C1[2][1]) */
  double C2[4][2]; /* 4 points, 3rd one is (C2[2][0],C2[2][1]) */
  double f[4];
  int i, im[4], p2[4], problem = 0;

  C1[0][0] = brct1[0];
  C1[0][1] = brct1[2];
  C1[1][0] = brct1[1];
  C1[1][1] = brct1[2];
  C1[2][0] = brct1[0];
  C1[2][1] = brct1[3];
  C1[3][0] = brct1[1];
  C1[3][1] = brct1[3];

  for(i=0; i<4; i++)
    p2[i] = Cpat2_of_Cpat1(pat1,f1,C1[i],   pat2,f2,C2[i]);

  /* check if C2 is not NAN */
  if(p2[0]<0 || p2[1]<0 || p2[2]<0 || p2[3]<0)
  {
    int j;

    /* check if all is ok, or if there is a NAN */
    for(i=0; i<4; i++)
      for(j=0; j<2; j++)
        if(isnan(C2[i][j])) problem = 1;
  }

  /* find min and max in coords */
  for(i=0; i<4; i++) f[i] = C2[i][0];
  brct2[0] = min_in_1d_array(f,4, im);

  for(i=0; i<4; i++) f[i] = C2[i][0];
  brct2[1] = max_in_1d_array(f,4, im);

  for(i=0; i<4; i++) f[i] = C2[i][1];
  brct2[2] = min_in_1d_array(f,4, im);

  for(i=0; i<4; i++) f[i] = C2[i][1];
  brct2[3] = max_in_1d_array(f,4, im);

  return problem;
}

/* check if a point is in brct,
   also sets C to boundary value if it is very close to the boundary */
int C_in_brct(const double brct[4], double C[2])
{
  int d;

  /* return 0 if C is clearly outside */
  for(d=0; d<2; d++)
  {
    if(dless(C[d],    brct[2*d]))   return 0;
    if(dgreater(C[d], brct[2*d+1])) return 0;
  }

  /* set C to boundary value if it is very close to boundary */
  for(d=0; d<2; d++)
  {
    if(C[d] < brct[2*d])   C[d] = brct[2*d];
    if(C[d] > brct[2*d+1]) C[d] = brct[2*d+1];
  }
  return 1;
}


/* convert X-coords on face f of node to X-coords of neighboring node,
   and write them into nbC */
void array_nbXface_of_Xface(tNode *node, int f,
                            tNode *nb, int nb_f, tArray *nbC[2])
{
  tPat *pat = node->pat;
  tMesh *mesh = pat->mesh;
  int *n = node->n;
  int dir = f/2;
  int pl = (n[dir]-1)*(f%2);
  //int d1 = Dir1_norm(dir);
  //int d2 = Dir2_norm(dir);
  int iX = Ind("X");
  int ix = Ind("x");
  double *pX[] = { Vard(node,iX), Vard(node,iX+1), Vard(node,iX+2) };
  double *px[] = { Vard(node,ix), Vard(node,ix+1), Vard(node,ix+2) };
  double *oC[] = { Arrd(nbC[0]), Arrd(nbC[1]) };
  int i,j,k, p;
  tPat *opat = nb->pat;
  int odir = nb_f/2;
  int od1 = Dir1_norm(odir);
  int od2 = Dir2_norm(odir);
  double nbrct[4];

  forplaneN(dir, i,j,k, n, pl)
  {
    int ijk = Ind_n(i,j,k, n);
    int ind = Ind_n_norm(i,j,k, n, dir);
    double x[] = { px[0][ijk], px[1][ijk], px[2][ijk] };
    double X[] = { pX[0][ijk], pX[1][ijk], pX[2][ijk] };
    //double C[] = { X[d1], X[d2] };
    double oX[3];

    p = p_XYZ_of_xyz(opat, oX, x);
    if(p<0)
    {
      /* we couldn't find, the nb point oX, but maybe the two coords on the
         face are ok? */
      double oCl[] = { oX[od1], oX[od2] };

      brct_nodeface(nb, odir, nbrct);
      if(!C_in_brct(nbrct, oCl))
      {
        /* just use something outside the nb's bound. rect., to signal to
           e.g. mark_points_in_nb_f that this is not be used */
        oX[0] = nb->bbox[1] * 2.;
        oX[1] = nb->bbox[3] * 2.;
        oX[2] = nb->bbox[5] * 2.;
      }
    }
    oC[0][ind] = oX[od1];
    oC[1][ind] = oX[od2];

    if(0) // if(p<0)
    {
      printnode(node);
      printf("f=%d  ", f);
      pr3v("X", X);
      pr3v("x", x);
      printf("\n");
      printnode(nb);
      //printbface(obface);
      //printf("obface->f=%d  ", obface->f);
      pr3v("oX", oX);
      printf("\nopat:  ");
      printpatch(opat);
      errorexit("x should be in opat!!!");
    }
  }
}


/* convert X-coords on face f of node to X-coords of neighboring node,
   write them into nbC, and mark points outside nb with -1 in nb I */
void array_find_nbXface_of_Xface(tNode *node, int f, tNode *nb, int nb_f,
                                 tArray *nbC[2], tArray *nbI)
{
  tPat *pat = node->pat;
  tMesh *mesh = pat->mesh;
  int *n = node->n;
  int dir = f/2;
  int pl = (n[dir]-1)*(f%2);
  int ix = Ind("x");
  double *px[] = { Vard(node,ix), Vard(node,ix+1), Vard(node,ix+2) };
  double *oC[] = { Arrd(nbC[0]), Arrd(nbC[1]) };
  int *oI = Arri(nbI);
  int i,j,k;
  long loc;
  int odir = nb_f/2;
  int od1 = Dir1_norm(odir);
  int od2 = Dir2_norm(odir);

  forplaneN(dir, i,j,k, n, pl)
  {
    int ijk = Ind_n(i,j,k, n);
    int ind = Ind_n_norm(i,j,k, n, dir);
    double x[] = { px[0][ijk], px[1][ijk], px[2][ijk] };
    double oX[3];

    /* find point x in nb */
    loc = l_XYZ_of_xyz(nb,-1, oX, x);
    if(loc>=0) /* point was found inside nb */
    {
      oC[0][ind] = oX[od1];
      oC[1][ind] = oX[od2];
      oI[ind] = ind;
    }
    else /* point is not inside opat */
      oI[ind] = -1;
  }
}
