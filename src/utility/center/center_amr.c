/* center.c */
/* Wolfgang Tichy 6/2024 */

#include "nmesh.h"
#include "center.h"

/* global pars for center */
tcenter center[1];



/* return refinement level that we want for elm
   We look at concentric spheres differing by a factor of 2 in radius. */
int center_amr_l(tElm *elm)
{
  tMesh *mesh = Elm_mesh(elm);
  double radius = 12.;
  int l_max = 4;
  int Nmax = 2;
  double cx[3][3]; /* 2 star centers: e.g. cx[1][3] = z-coord of center1 */
  double r;
  int N, dir, l, l_ref;

  /* get 2 star centers */
  for(N=1; N<=Nmax; N++)
    for(dir=0; dir<3; dir++)
      cx[N][dir] = Getd(center->cx[1][dir]);

  /* set refinement level l_ref that we want for this elm */
  l_ref = 0;
  r = radius * pow(2., l_max-1); /* radius separating level 0 and 1 */
  for(l=0; l<l_max; l++)
  {
    for(N=1; N<=Nmax; N++)
      if(elmpoints_any_in_sphere(elm, cx[N], r))
        l_ref = l+1;

    if(l_ref<=l) break;
    r = r*0.5; /* shrink r for next level l */
  }
  return l_ref;
}
