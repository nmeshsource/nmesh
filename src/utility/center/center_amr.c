/* center.c */
/* Wolfgang Tichy 6/2024 */

#include "nmesh.h"
#include "center.h"


/* global pars for center */
extern tcenter center[1];

/* struct for more pars */
typedef struct {
  int lmax;         /* max ref.-level */
  double radius[3]; /* radius[1] is radius of max ref.-level around center1 */
  int N_first;      /* id of first center */
  int N_last;       /* id of last center */
} tcenter_amr_pars;



/* return refinement level that we want for elm
   We look at concentric spheres differing by a factor of 2 in radius. */
int center_amr_l(tElm *elm, void *pars)
{
  tMesh *mesh = Elm_mesh(elm);
  tcenter_amr_pars *center_amr_pars = pars;
  int l_max      = center_amr_pars->lmax;
  double *radius = center_amr_pars->radius;
  int N_first    = center_amr_pars->N_first;
  int N_last     = center_amr_pars->N_last;
  double cx[3][3]; /* 2 star centers: e.g. cx[1][3] = z-coord of center1 */
  double r[3];
  int N, dir, l, l_ref;

  /* get 2 star centers */
  for(N=N_first; N<=N_last; N++)
    for(dir=0; dir<3; dir++)
      cx[N][dir] = Getd(center->cx[N][dir]);
  /*
  PRF;
  for(N=N_first; N<=N_last; N++)
    for(dir=0; dir<3; dir++)
      printf(" %g", cx[N][dir]);
  printf("\n");
  */

  /* radius separating level 0 and 1 */
  for(N=N_first; N<=N_last; N++) r[N] = radius[N] * pow(2., l_max-1);

  /* set refinement level l_ref that we want for this elm */
  l_ref = 0;
  for(l=0; l<l_max; l++)
  {
    //for(N=N_first; N<=N_last; N++) printf("l%d r[%d]=%g\n", l, N, r[N]);
    for(N=N_first; N<=N_last; N++)
      if(elmpoints_any_in_sphere(elm, cx[N], r[N]))
      {
        //PRFs(": ");printeploc(elm->eploc);
        //printf(" has point in r[%d]\n", N);
        l_ref = l+1;
        break;
      }

    if(l_ref<=l) break;

    /* shrink r for next level l */
    for(N=N_first; N<=N_last; N++) r[N] *= 0.5;
  }
  //PRFs(": ");printeploc(elm->eploc);
  //printf(" l_ref=%d\n", l_ref);
  return l_ref;
}

/* h-refine until each elm has ref.-level given by func center_amr_l */
int center_amr(tMesh *mesh)
{
  double dt = Getd(Par("center_amr_time"));
  if(dt >= 0.)
  {
    tcenter_amr_pars pars = {.lmax    = Geti(Par("center1_amr_lmax")),
                             .radius  = {0., Getd(Par("center1_amr_radius")),
                                             Getd(Par("center2_amr_radius")) },
                             .N_first = 1,
                             .N_last  = 2};
    if(!Getv(Par("center2_amr_lmax"), "center1_amr_lmax"))
      errorexit("currently we need center2_amr_lmax = center1_amr_lmax");

    if(TimeIsAt_di_dt(mesh, -1, dt))
    {
      struct timespec tp0[1];
      struct timespec tp1[1];
      getRealTime(tp0);
      hadapt_to_desired_l(mesh, center_amr_l, &pars);
      getRealTime(tp1);
      PRF;
      printf(": hadapt_to_desired_l took %gs\n", getTimeDiffIn_s(tp1, tp0));
    }
  }
  return 0;
}
