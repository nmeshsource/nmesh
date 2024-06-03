/* center.c */
/* Wolfgang Tichy 6/2024 */

#include "nmesh.h"
#include "center.h"


/* global pars for center */
extern tcenter center[1];

/* struct for more pars */
typedef struct {
  int lmax;
  int radius;
} tcenter_amr_pars;



/* return refinement level that we want for elm
   We look at concentric spheres differing by a factor of 2 in radius. */
int center_amr_l(tElm *elm, void *pars)
{
  tMesh *mesh = Elm_mesh(elm);
  tcenter_amr_pars *center_amr_pars = pars;
  int l_max     = Getd(center_amr_pars->lmax);
  double radius = Getd(center_amr_pars->radius);
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

/* h-refine until each elm has ref.-level given by func center_amr_l */
int center_amr(tMesh *mesh)
{
  double dt = Getd(Par("center_amr_time"));
  tcenter_amr_pars pars = { .lmax   = Par("center_amr_lmax"),
                            .radius = Par("center_amr_radius")};

  if(dt >= 0. && TimeIsAt_di_dt(mesh, -1, dt))
    hadapt_to_desired_l(mesh, center_amr_l, &pars);

  return 0;
}
