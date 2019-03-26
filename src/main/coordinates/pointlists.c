/* pointlists.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"



/* allocate and make a point index list along a line in direc. dir,
   needs to be freed by caller */
intList *pointindexList_line(tNode *node, int dir, int axis[])
{
  intList *plist = alloc_intList(); /* list that contains other point indices */
  int *n = node->n;
  int i,j,k;
  int imin, jmin, kmin;
  int imax, jmax, kmax;

  imin = jmin = kmin=0;
  imax = n[0] - 1;
  jmax = n[1] - 1;
  kmax = n[2] - 1;

  switch(dir)
  {
  case 0:
    if(axis[1]>jmax) jmin = jmax;
    else             jmin = jmax = axis[1];
    if(axis[2]>kmax) kmin = kmax;
    else             kmin = kmax = axis[2];
    break;
  case 1:
    if(axis[0]>imax) imin = imax;
    else             imin = imax = axis[0];
    if(axis[2]>kmax) kmin = kmax;
    else             kmin = kmax = axis[2];
    break;
  case 2:
    if(axis[0]>imax) imin = imax;
    else             imin = imax = axis[0];
    if(axis[1]>jmax) jmin = jmax;
    else             jmin = jmax = axis[1];
    break;
  default:
    errorexit("dir has to be 0,1,2");
  }

  /* go over line and add points to list */
  for(k=kmin; k<=kmax; k++)
    for(j=jmin; j<=jmax; j++)
      for(i=imin; i<=imax; i++)
      {
        int ind = Ind_n(i,j,k, n);
        push_intList(plist, ind);
      }

  return plist;
}

/* allocate and make a point index list in a single plane,
   needs to be freed by caller */
intList *pointindexList_plane(tNode *node, int normal, int plane[])
{
  intList *plist = alloc_intList(); /* list that contains other point indices */
  int *n = node->n;
  int i,j,k;
  int imin, jmin, kmin;
  int imax, jmax, kmax;

  imin = jmin = kmin=0;
  imax = n[0] - 1;
  jmax = n[1] - 1;
  kmax = n[2] - 1;

  switch(normal)
  {
  case 0:
    if(plane[0]>imax) imin = imax;
    else              imin = imax = plane[0];
    break;
  case 1:
    if(plane[1]>jmax) jmin = jmax;
    else              jmin = jmax = plane[1];
    break;
  case 2:
    if(plane[2]>kmax) kmin = kmax;
    else              kmin = kmax = plane[2];
    break;
  default:
    errorexit("dir has to be 0,1,2");
  }

  /* go over plane, with normal */
  for(k=kmin; k<=kmax; k++)
    for(j=jmin; j<=jmax; j++)
      for(i=imin; i<=imax; i++)
      {
        int ind = Ind_n(i,j,k, n);
        push_intList(plist, ind);
      }

  return plist;
}

/* allocate and make a point index list for a entire node,
   needs to be freed by caller */
intList *pointindexList_node(tNode *node)
{
  intList *plist = alloc_intList(); /* list that contains other point indices */
  int ind;

  forpoints(node,ind) push_intList(plist, ind);

  return plist;
}
