/* setup_CubedSpheres.c */
/* Wolfgang Tichy, Nov 2017 */
/* create various arrangements of cubed spheres */

#include "nmesh.h"
#include "coordinates.h"


/************************************************************************/
/* make patches for 2 stars as in DNSdata  */
/************************************************************************/

/* make this:
    __________
   |\    9   /|
   | \  __  / |
   |  \/3 \/  |
   |6 /\__/\  |
   | |0|__|1| |
   |  \/2 \/ 7|
   |  /\__/\  |
   | /   8  \ |
   |/________\|

   returns the index of the pat right after the last converted pat */
int arrange_12CubSph_into_empty_cube(tMesh *mesh, int N, double *xc,
                                     double din, double dmid, double dout)
{
  int pl;
  int i;
  double Din[6];
  double Dmid[6];
  double Dout[6];

  /* set distances from center to din, dmid, dout for all domains */
  for(i=0; i<6; i++)
  {
    Din[i] = din;
    Dmid[i]= dmid;
    Dout[i]= dout;
  }
  /* convert the 12 pats */
  pl = add_N_CubedSphere_pats(mesh, N, outerCubedSphere,0,1, xc, Din,Dmid);
  pl = add_N_CubedSphere_pats(mesh, N, innerCubedSphere,0,1, xc, Dmid,Dout);
  return pl;
}

/* same as arrange_12CubSph_into_empty_cube, but add one cube at the center */
int arrange_1pat12CubSph_into_full_cube(tMesh *mesh, int N, double *xc,
                                        double din, double dmid, double dout)
{
  int pl;
  pl = add_1cube_pat(mesh, xc, din);
  pl = arrange_12CubSph_into_empty_cube(mesh,N, xc, din,dmid,dout);
  return pl;
}


/* make two big touching cubes like this
    __________ __________
   |\        /|\        /|
   | \  __  / | \  __  / |
   |  \/  \/  |  \/  \/  |
   |  /\__/\  |  /\__/\  |
   | | |__| | | | |__| | |
   |  \/  \/  |  \/  \/  |
   |  /\__/\  |  /\__/\  |
   | /      \ | /      \ |
   |/________\|/________\|
*/
int two_full_cubes_touching_at_x0(tMesh *mesh, int N, double dc,
                                  double din1, double dmid1,
                                  double din2, double dmid2)
{
  int pl;
  double xc[3];

  /* put centers on x-axis */
  xc[1] = xc[2] = 0.0;

  /* full cube 1 is centered at xc[1]=dc and has width 2dc */
  xc[0] = dc;
  pl = arrange_1pat12CubSph_into_full_cube(mesh,N, xc, din1,dmid1, dc);

  /* full cube 2 is centered at xc[1]=-dc and has width 2dc */
  xc[0] = -dc;
  pl = arrange_1pat12CubSph_into_full_cube(mesh,N, xc, din2,dmid2, dc);

  return pl;
}


/* surround two big touching cubes with cubed spheres that have
   A,B that are in an extended range
         _______
      __/       \__
     /             \
    /__     3     __\
   /   -- _____ --   \     e.g. dom0/1 have A = [-0.8, 0.8]
  |   0  |  |  |   1  |         dom2/3 have A = [-1.25,1.25]
  |      |__|__|      |
   \ __--       --__ /      r0 is radius of outer sphere
    \       2       /       2*dc is sidelength of one cube
     \__         __/
        \_______/
*/
int sphere_around_two_full_cubes_touching_at_x0(tMesh *mesh, int N,
        double dc, double din1, double dmid1, double din2, double dmid2,
        double r0)
{
  int pl;
  double xc[3], Din[6], Dout[6];
  int i;

  /* make the 2 full cubes */
  pl = two_full_cubes_touching_at_x0(mesh,N, dc, din1,dmid1, din2,dmid2);

  /* set distances to make 6 more cubed spheres around these 2 full cubes */
  for(i=0; i<6; i++)
  {
    if(i<2) Din[i] = 2.0*dc;
    else    Din[i] = dc;
    Dout[i] = r0;
  }
  xc[0] = xc[1] = xc[2] = 0.0;
  pl = add_N_CubedSphere_pats(mesh, N, outerCubedSphere,0,0, xc, Din,Dout);
  return pl;
}


/* put 6 stretchedCubedShell's around the sphere from
   sphere_around_two_full_cubes_touching_at_x0
                   ___________
             _____/    r1     \______
          __/                        \_
         /                             \
       _- \-                         -/
      /     \-       _______       -/
     |        \-  __/  r0   \__  -/   ...
    |           \/             \/
   /            /__           __\
  |            /   -- _____ --   \    r0 is radius of inner sphere (sphere0)
  |           |      |  |  |      |   r1 is radius of outer sphere (sphere1)
  |           |      |__|__|      |
  |            \ __--       --__ /
   \            \               /
    |          _/\__         __/
     |       _/     \_______/
      \    _/
       -_ /
         \      ...
*/
int two_spheres_around_two_full_cubes(tMesh *mesh, int N,
        double dc, double din1, double dmid1, double din2, double dmid2,
        double r0, double r1)
{
  int pl;
  double xc[3], Din[6], Dout[6];
  int i;

  /* make the 2 full cubes and sphere0 around them */
  pl = sphere_around_two_full_cubes_touching_at_x0(mesh,N, dc,
                                                   din1,dmid1, din2,dmid2, r0);
  /* set distances to make 6 more stretched cubed shells around sphere0 */
  for(i=0; i<6; i++)
  {
    Din[i]  = r0;
    Dout[i] = r1;
  }
  xc[0] = xc[1] = xc[2] = 0.0;
  pl = add_N_CubedSphere_pats(mesh, N, CubedShell,1,1, xc, Din,Dout);
  return pl;
}

/************************************************************************/
/* make patches that surround a central box */
/************************************************************************/

/* surround a box with cubed spheres that have
   A,B that are in an extended range
         _______
      __/       \__
     /             \
    /__     3     __\
   /   -- _____ --   \     e.g. dom0/1 have A = [-0.8, 0.8]
  |   0  |     |   1  |         dom2/3 have A = [-1.25,1.25]
  |      |_____|      |
   \ __--       --__ /      r0 is radius of outer sphere
    \       2       /       2*dc[i] is sidelength of box in dir i
     \__         __/
        \_______/
*/
int sphere_around_empty_box_at_xc(tMesh *mesh, int N,
                                  double xc[3], double dc[3], double r0)
{
  int pl;
  double Din[6], Dout[6];
  int f;

  /* set distances to make 6 cubed spheres around the box */
  for(f=0; f<6; f++)
  {
    Din[f]  = dc[f/2];
    Dout[f] = r0;
  }
  pl = add_N_CubedSphere_pats(mesh, N, outerCubedSphere,0,0, xc, Din,Dout);
  return pl;
}

/* same as sphere_around_empty_box_at_xc, but add one box at the center */
int sphere_around_full_box_at_xc(tMesh *mesh, int N,
                                 double xc[3], double dc[3], double r0)
{
  int pl;
  pl = add_1box_pat(mesh, xc, dc);
  pl = sphere_around_empty_box_at_xc(mesh,N, xc, dc, r0);
  return pl;
}

/* Surround a box with cubed sphere domains that have A,B that in extended
   range. Subdivide each domain.

         ________         e.g. dom0/1 have A = [-0.8, 0.8]
      __/   |    \__           dom2/6 have A = [-1.25,-0.625]
     /   \ 7| 8 /   \          dom3/7 have A = [-0.625,0]
    /__ 6 \ |  / 9 __\
   /   -- ______ --   \   nAB_x[i] is number of subdivisions of A or B
  |  0   |      |   1  |  in the x^i-direction. nBmax = max(nAB_x)
  |      |______|      |
   \ __-- / | \  --__ /   nlam_AB[j][k] is number of subdivisions of lam for
    \  2 /  |  \ 5   /    piece (j,k) in A,B.
     \__/ 3 | 4 \ __/
        \________/        r0 is radius of outer sphere
                          2*dc[i] is sidelength of box in dir i
*/
int sphere_nABnlam_around_empty_box_at_xc(tMesh *mesh, int N, double xc[3],
                                          double dc[3], double r0,
                                          int nAB_x[3], int nBmax,
                                          int (*nlam_AB)[nBmax])
{
  int pl;
  double Din[6], Dout[6];
  int f;

  /* set distances to make 6 cubed spheres around the box */
  for(f=0; f<6; f++)
  {
    Din[f]  = dc[f/2];
    Dout[f] = r0;
  }
  pl = add_N_CubedSphere_doms(mesh, N, outerCubedSphere,0,0, xc, Din,Dout,
                              nAB_x, nBmax, nlam_AB);
  return pl;
}
/* like sphere_nABnlam_around_empty_box_at_xc, but set all nlam_AB[j][k]=1 */
int sphere_nAB_around_empty_box_at_xc(tMesh *mesh, int N,
                                      double xc[3], double dc[3], double r0,
                                      int nAB_x[3])
{
  int nABmax = max3(nAB_x[0], nAB_x[1], nAB_x[2]);
  int (*nlam_AB)[nABmax] = malloc(nABmax*nABmax*sizeof(int));
  int j,k, pl;

  /* set all nlam_AB to 1 */
  for(k=0; k<nABmax; k++)
    for(j=0; j<nABmax; j++)
      nlam_AB[j][k] = 1;

  pl = sphere_nABnlam_around_empty_box_at_xc(mesh, N, xc, dc, r0,
                                             nAB_x, nABmax,nlam_AB);
  free(nlam_AB);
  return pl;
}



/* put 6 stretchedCubedShell's around the sphere from
   sphere_around_full_box_at_xc
                   ___________
             _____/    r1     \______
          __/                        \_
         /                             \
       _- \-                         -/
      /     \-       _______       -/
     |        \-  __/  r0   \__  -/   ...
    |           \/             \/
   /            /__           __\
  |            /   -- _____ --   \    r0 is radius of inner sphere (sphere0)
  |           |      |     |      |   r1 is radius of outer sphere (sphere1)
  |           |      |_____|      |
  |            \ __--       --__ /    if stretch=0 use xyz_of_lamAB_CubSph
   \            \               /     if stretch=1 use xyz_of_rhoAB_CubSph
    |          _/\__         __/      in outermost CubedShell
     |       _/     \_______/
      \    _/
       -_ /
         \      ...
*/
int two_spheres_around_box_at_xc(tMesh *mesh, int N,
                                 double xc[3], double dc[3],
                                 double r0, double r1, int stretch)
{
  int pl;
  double Din[6], Dout[6];
  int i;

  /* make the full box and sphere0 around them */
  pl = sphere_around_full_box_at_xc(mesh,N, xc, dc, r0);

  /* set distances to make 6 more stretched cubed shells around sphere0 */
  for(i=0; i<6; i++)
  {
    Din[i]  = r0;
    Dout[i] = r1;
  }
  pl = add_N_CubedSphere_pats(mesh, N, CubedShell,stretch,0, xc, Din,Dout);
  return pl;
}

/* same as two_spheres_around_box_at_xc, but put no box at center */
int two_spheres_around_empty_box_at_xc(tMesh *mesh, int N, double xc[3],
                                       double dc[3], double r0, double r1,
                                       int stretch)
{
  int pl;
  double Din[6], Dout[6];
  int i;

  /* make sphere0 around empty box */
  pl = sphere_around_empty_box_at_xc(mesh,N, xc, dc, r0);

  /* set distances to make 6 more stretched cubed shells around sphere0 */
  for(i=0; i<6; i++)
  {
    Din[i]  = r0;
    Dout[i] = r1;
  }
  pl = add_N_CubedSphere_pats(mesh, N, CubedShell,stretch,0, xc, Din,Dout);
  return pl;
}

/* put 6 more stretchedCubedShell's around the sphere from
   two_spheres_around_box_at_xc
                   ___________
             _____/    r2     \______
          __/        _______         \_
         /      ____/  r1   \____      \
       _- \-   /                 \   -/
      /     \-       _______       -/
     |     /  \-  __/  r0   \__  -/   ...
    |     /     \/             \/
   /     /      /__           __\
  |     |      /   -- _____ --   \    r0 is radius of inner sphere (sphere0)
  |     |     |      |     |      |   r1 is radius of middle sphere (sphere1)
  |     |     |      |_____|      |   r2 is radius of outer sphere (sphere2)
  |     |      \ __--       --__ /
   \     \      \               /
    |     \    _/\__         __/
     |     \ _/     \_______/         if stretch=0 use xyz_of_lamAB_CubSph
      \    _/                         if stretch=1 use xyz_of_rhoAB_CubSph
       -_ /                           in outermost CubedShell
         \      ...
*/
int three_spheres_around_box_at_xc(tMesh *mesh, int N,
                                   double xc[3], double dc[3],
                                   double r0, double r1, double r2,
                                   int stretch)
{
  int pl;
  double Din[6], Dout[6];
  int stretch1, i;

  /* make the full box and sphere0 around them */
  stretch1 = 0;
  pl = two_spheres_around_box_at_xc(mesh,N, xc, dc, r0, r1, stretch1);

  /* set distances to make 6 more stretched cubed shells around sphere1 */
  for(i=0; i<6; i++)
  {
    Din[i]  = r1;
    Dout[i] = r2;
  }
  pl = add_N_CubedSphere_pats(mesh, N, CubedShell,stretch,0, xc, Din,Dout);
  return pl;
}

/* same as three_spheres_around_box_at_xc, but put no box at center */
int three_spheres_around_empty_box_at_xc(tMesh *mesh, int N,
                                         double xc[3], double dc[3],
                                         double r0, double r1, double r2,
                                         int stretch)
{
  int pl;
  double Din[6], Dout[6];
  int stretch1, i;

  /* make the full box and sphere0 around them */
  stretch1 = 0;
  pl = two_spheres_around_empty_box_at_xc(mesh,N, xc, dc, r0, r1, stretch1);

  /* set distances to make 6 more stretched cubed shells around sphere1 */
  for(i=0; i<6; i++)
  {
    Din[i]  = r1;
    Dout[i] = r2;
  }
  pl = add_N_CubedSphere_pats(mesh, N, CubedShell,stretch,0, xc, Din,Dout);
  return pl;
}

/************************************************************************/
/* make patches that surround an inner sphere */
/************************************************************************/

/* surround an inner excised sphere with cubed spheres to form a shell
         _______
      __/       \__
     /_     3     _\
    /  \_       _/  \       rin is radius of inner excised sphere
   /     \_---_/     \      rout is radius of outer sphere
  |   0  /     \   1  |
  |      \_   _/      |
   \    _/ --- \_    /
    \__/         \_ /
     \__    2    __/
        \_______/
*/
int CubedSphere_shell_at_xc(tMesh *mesh, int N, double xc[3],
                            double rin, double rout)
{
  int pl;
  double Din[6], Dout[6];
  int f;

  /* set distances to make 6 cubed spheres around the empty inner sphere */
  for(f=0; f<6; f++)
  {
    Din[f]  = rin;
    Dout[f] = rout;
  }
  pl = add_N_CubedSphere_pats(mesh, N, CubedShell,0,0, xc, Din,Dout);
  return pl;
}

/* put 6 stretchedCubedShell's around the sphere from
   CubedSphere_shell_at_xc
                   ___________
             _____/    r2     \______
          __/                        \_
         /                             \
       _- \-                         -/
      /     \-       _______       -/
     |        \-  __/  r1   \__  -/   ...
    |           \/_           _\/
   /            /  \_       _/  \
  |            /     \_---_/     \    first shell extends from r0 to r1
  |           |      /     \      |   second shell extends from r1 to r2
  |           |      \_   _/      |
  |            \    _/ --- \_    /    if stretch=0 use xyz_of_lamAB_CubSph
   \            \__/         \_ /     if stretch=1 use xyz_of_rhoAB_CubSph
    |          _/\__         __/      in outermost CubedShell
     |       _/     \_______/
      \    _/
       -_ /
         \      ...
*/
int two_CubedSphere_shells_at_xc(tMesh *mesh, int N,
                                 double xc[3], double r0,
                                 double r1, double r2, int stretch)
{
  int pl;
  double Din[6], Dout[6];
  int f;

  /* make inner shell */
  pl = CubedSphere_shell_at_xc(mesh,N, xc, r0, r1);

  /* set distances to make 6 cubed spheres around the inner shell */
  for(f=0; f<6; f++)
  {
    Din[f]  = r1;
    Dout[f] = r2;
  }
  pl = add_N_CubedSphere_pats(mesh, N, CubedShell,stretch,0, xc, Din,Dout);

  return pl;
}


/************************************************************************/
/* make test patches */
/************************************************************************/

/* one wedge that touches two others, all are domain1
              _____
        ___--- p1  \
   __--    \   dom1 \          e.g. p0 has A = [-1, 1]
  |   p0    |________|         p1/2 A = [0, 1] / A = [-1, 0]
  |__ dom1  |  p2    |
     -- ___/   dom1 /          r0 is radius of inner sphere
           ---____ /           r1 is radius of outer sphere
*/
int two_wegdes_touching_1_wedge(tMesh *mesh, double dc, double r0, double r1)
{
  double xc[] = { 0., 0., 0. };
  double ABrct0[] = { -1., 1., -1., 1. };
  double ABrct1[] = {  0., 1., -1., 1. };
  double ABrct2[] = { -1., 0., -1., 1. };
  int pl;

  add_1_CubedSphere_pat(mesh, 1, outerCubedSphere,0,0, xc, dc,r0, ABrct0);
  add_1_CubedSphere_pat(mesh, 1, CubedShell,0,0, xc, r0,r1, ABrct1);
  pl = add_1_CubedSphere_pat(mesh, 1, CubedShell,0,0, xc, r0,r1, ABrct2);

  return pl;
}

/* one wedge that touches two others, where different domain types are
   involved  ___
            /   \__
           /       \_           e.g.:
          /   p2     \          p0 has A = [0, 2]
         __   dom3    |         p1 has A = [0, 1]
        /  \_       _/ \        p2 has A = [1, 2]
       /     \    _/    \_
      /       \__/        \     r0 is radius of inner sphere
     /   p0    |    p1     \    r1 is radius of outer sphere
    |    dom1   |   dom1    |
    |___________|___________|
*/
int two_diff_wegdes_touching_1_wedge(tMesh *mesh, double dc,
                                     double r0, double r1)
{
  double xc[] = { 0., 0., 0. };
  double ABrct0[] = { 0., 2., -1., 1. };
  double ABrct1[] = { 0., 1., -1., 1. };
  double ABrct2[] = { 0.5, 1., -1., 1. };
  int pl;

  add_1_CubedSphere_pat(mesh, 1, outerCubedSphere,0,0, xc, dc,r0, ABrct0);
  add_1_CubedSphere_pat(mesh, 1, CubedShell,0,0, xc, r0,r1, ABrct1);
  pl = add_1_CubedSphere_pat(mesh, 3, CubedShell,0,0, xc, r0,r1, ABrct2);

  return pl;
}


/************************************************************************/
/* basic functions to make one box or a cubed sphere shell */
/************************************************************************/

/* primitive default forconstant FSurf funcs */
int FSurf_is_CI_s(tPat *pat, int si, double C[2], double *F)
{
  *F = pat->CI->s[si];
  return 0;
}
int dFSurfdC_is_zero(tPat *pat, int si, double C[2], double dF[2])
{
  dF[0] = dF[1] = 0.;
  return 0;
}




/* add 6 pats with some kind of cubed spheres */
/* call this before any pats exist already */
/* type can be "PyramidFrustum", "innerCubedSphere", "outerCubedSphere",
   "CubedShell"
   xc[0..2] = (x,y,z) of coord center for the 6 cubed spheres
   Din[0...5] inner distance from center for cubed sph. domain 0-5
   Dout[0...5] outer distance from center for cubed sph. domain 0-5 */
/* It returns the index of the last added pat */
int add_6CubedSphere_pats(tMesh *mesh, int type, int stretch, int SigFunc,
                          double *xc, double *Din, double *Dout)
{
  return add_N_CubedSphere_pats(mesh, 6, type, stretch, SigFunc,
                                xc, Din,Dout);
}
/* same as above but add only N of the cubed sphere pats */
int add_N_CubedSphere_pats(tMesh *mesh, int N,
                           int type, int stretch, int SigFunc,
                           double *xc, double *Din, double *Dout)
{
  int i, ret;

  if(N<1 || N>6) errorexit("N must be 1,2,3,4,5,6");

  /* make the N domains */
  for(i=0; i<N; i++)
  {
    double Amin,Amax, Bmin,Bmax;
    double ABrct[4];

    /* set min/max in A-, B-directions */
    set_AB_min_max_from_Din(i, Din, &Amin,&Amax, &Bmin,&Bmax);
    ABrct[0] = Amin;
    ABrct[1] = Amax;
    ABrct[2] = Bmin;
    ABrct[3] = Bmax;

    /* add 1 Cubed Sphere */
    ret = add_1_CubedSphere_pat(mesh, i,type, stretch,SigFunc,
                                xc,Din[i],Dout[i], ABrct);
  }

  return ret; /* return pat index of last added pat */
}

/* given (Amin,Amax) compute angles on patch
   In: plus, Amin,Amax.  Out: alphamin,alphamax
   plus=0 means on left  (e.g. dom 0,2,4)
   plus=1 means on right (e.g. dom 1,3,5) */
void angle_of_Arange_or_Brange(int plus, double Amin, double Amax,
                             double *alphamin, double *alphamax)
{
  switch(plus)
  {
  case 0:
    *alphamin = Arg_plus(-1., -Amin);
    *alphamax = Arg_plus(-1., -Amax);
    break;
  case 1:
    *alphamin = Arg(1., Amin);
    *alphamax = Arg(1., Amax);
    break;
  }
}

/* get A at index j*/
double A_or_B_of_index(int nA, int j, int mode,
                       double Amin,double Amax,
                       double alphamin, double alphamax)
{
  double dA, dalpha;

  /* use exact end values for A or B */
  if(j==0)  return Amin;
  if(j==nA) return Amax;

  switch(mode)
  {
  case 0: /* linear in A */
      dA = (Amax - Amin)/nA;
    return (Amin + dA*j);
    break;
  case 1: /* linear in angle */
    dalpha = (alphamax - alphamin)/nA;
    return tan(alphamin + dalpha*j);
  default:
    errorexit("mode must be: 0,1");
  }
}

/* Add cubed sphere doms 0 to N-1:
   +Divide A or B region of a dom into nAB_x pieces. Here nAB_x[i] tells us
    how many pieces we want in the x^i-direction.
   +nBmax = max3(nAB_x[0], nAB_x[1], nAB_x[2])
   +Also divide lam into nlam_AB[A_piece][B_piece] pieces,
    where A/B_piece=0,...,nA/B */
int add_N_CubedSphere_doms(tMesh *mesh, int N,
                           int type, int stretch, int SigFunc,
                           double xc[3], double Din[6], double Dout[6],
//                           int AB_div_mode,
                           int nAB_x[3], int nBmax, int (*nlam_AB)[nBmax])
{
  int AB_div_mode = 1;
  int f, ret;

  if(N<1 || N>6) errorexit("N must be 1,2,3,4,5,6");

  /* make the N domains */
  for(f=0; f<N; f++)
  {
    double Amin,Amax, Bmin,Bmax;
    double alphamin, alphamax;
    double betamin,  betamax;
    int nA, nB, j, k;
    int dir = f/2;
    int pls = f%2;

    /* set min/max in A-, B-directions */
    set_AB_min_max_from_Din(f, Din, &Amin,&Amax, &Bmin,&Bmax);

    /* find number of pieces in A, B */
    switch(dir)
    {
    case 0:
      nA = nAB_x[1];
      nB = nAB_x[2];
      break;
    case 1:
      nA = nAB_x[0];
      nB = nAB_x[2];
      break;
    case 2:
      nA = nAB_x[1];
      nB = nAB_x[0];
      break;
    default:
      errorexit("dir must 0,1,2");
    }

    /* find angles from A and B extrema */
    angle_of_Arange_or_Brange(pls, Amin,Amax, &alphamin,&alphamax);
    angle_of_Arange_or_Brange(pls, Bmin,Bmax, &betamin,&betamax);

    /* make nlam*nA*nB patches */
    for(k=0; k<nB; k++)
    {
      double B0 = A_or_B_of_index(nB, k, AB_div_mode, Bmin,Bmax, betamin,betamax);
      double B1 = A_or_B_of_index(nB, k+1, AB_div_mode, Bmin,Bmax, betamin,betamax);
      for(j=0; j<nA; j++)
      {
        double A0 = A_or_B_of_index(nA, j, AB_div_mode, Amin,Amax, alphamin,alphamax);
        double A1 = A_or_B_of_index(nA, j+1, AB_div_mode, Amin,Amax, alphamin,alphamax);
        int nlam = nlam_AB[j][k];
        double dlam = 1./nlam;
        int i;
        for(i=0; i<nlam; i++)
        {
          double lam0 = dlam*i;
          double lam1 = dlam*(i+1);
          double bbox[] = {lam0,lam1, A0,A1, B0,B1};

          /* add 1 Cubed Sphere */
          ret = add_1_CubedSphere_pat_bbox(mesh, f,type, stretch,SigFunc,
                                           xc,Din[f],Dout[f], bbox);
        }
      }
    } /* end k-loop */
  }

  return ret; /* return pat index of last added pat */
}


/* find Amax,Amin, Bmax,Bmin in a domain using distances from center */
void set_AB_min_max_from_Din(int dom, double *Din,
                             double *Amin, double *Amax,
                             double *Bmin, double *Bmax)
{
  switch(dom)
  {
    case 0:
      *Amin = -Din[3]/Din[dom];
      *Amax =  Din[2]/Din[dom];
      *Bmin = -Din[5]/Din[dom];
      *Bmax =  Din[4]/Din[dom];
      break;

    case 1:
      *Amin = -Din[2]/Din[dom];
      *Amax =  Din[3]/Din[dom];
      *Bmin = -Din[4]/Din[dom];
      *Bmax =  Din[5]/Din[dom];
      break;

    case 2:
      *Amin = -Din[1]/Din[dom];
      *Amax =  Din[0]/Din[dom];
      *Bmin = -Din[5]/Din[dom];
      *Bmax =  Din[4]/Din[dom];
      break;

    case 3:
      *Amin = -Din[0]/Din[dom];
      *Amax =  Din[1]/Din[dom];
      *Bmin = -Din[4]/Din[dom];
      *Bmax =  Din[5]/Din[dom];
      break;

    case 4:
      *Amin = -Din[3]/Din[dom];
      *Amax =  Din[2]/Din[dom];
      *Bmin = -Din[1]/Din[dom];
      *Bmax =  Din[0]/Din[dom];
      break;

    case 5:
      *Amin = -Din[2]/Din[dom];
      *Amax =  Din[3]/Din[dom];
      *Bmin = -Din[0]/Din[dom];
      *Bmax =  Din[1]/Din[dom];
      break;

   default:
      errorexit("set_AB_min_max_from_D: unknown dom");
  }
}

/* add 1 cubed sphere with specific bounds and bbox */
int add_1_CubedSphere_pat_bbox(tMesh *mesh, int dom, int type,
                               int stretch, int SigFunc, double *xc,
                               double Din, double Dout, double bbox[6])
{
  int amr_n0 = Geti(Par("amr_n0"));
  int amr_n1 = Geti(Par("amr_n1"));
  int amr_n2 = Geti(Par("amr_n2"));
  int n[] = { amr_n0, amr_n1, amr_n2 };
  tPat *pat;
  int d;

  /* make new patch */
  pat = add_patch(mesh, bbox, NULL, n, 0);
  //PRFs(": ");printpatch(pat);

  /* set center */
  for(d=0; d<3; d++)  pat->CI->xc[d] = xc[d];

  /* set inner and outer distance from center */
  pat->CI->s[0] = Din;
  pat->CI->s[1] = Dout;

  /* set domain index and type */
  pat->CI->dom = dom;
  pat->CI->type= type;

  /* use info in pat and stretch, SigFunc to set rest of CI */
  set_1_CubedSphere_pat(pat, stretch, SigFunc);

  return pat->p; /* return pat index of last added pat */
}
/* add 1 cubed sphere with specific bounds, and lam in [0,1] */
int add_1_CubedSphere_pat(tMesh *mesh, int dom, int type,
                          int stretch, int SigFunc, double *xc,
                          double Din, double Dout, double ABrct[4])
{
  /* set min/max in each direction */
  double bbox[] = { 0.,1.,  ABrct[0],ABrct[1], ABrct[2],ABrct[3] };

  return add_1_CubedSphere_pat_bbox(mesh, dom,type, stretch,SigFunc,
                                    xc, Din,Dout, bbox);
}

/* set 1 cubed sphere pat from patch info, as well as stretch, SigFunc */
int set_1_CubedSphere_pat(tPat *pat, int stretch, int SigFunc)
{
  tMesh *mesh = pat->mesh;
  int isigma0    = Ind("CubedSphere_sigma0");
  int isigma0_dA = Ind("CubedSphere_dsigma0_dA");
  int isigma0_dB = Ind("CubedSphere_dsigma0_dB");
  int isigma1    = Ind("CubedSphere_sigma1");
  int isigma1_dA = Ind("CubedSphere_dsigma1_dA");
  int isigma1_dB = Ind("CubedSphere_dsigma1_dB");

  /* set coord trafos */
  switch(stretch)
  {
  case 0: /* stretch = 0, i.e. no stretching */
    pat->xyz_of_XYZ = xyz_of_lamAB_CubSph;
    pat->XYZ_of_xyz = lamAB_of_xyz_CubSph;
    pat->dXYZ_dxyz  = dlamAB_dxyz_CubSph;
    break;
  case 1:
    pat->xyz_of_XYZ = xyz_of_rhoAB_CubSph;
    pat->XYZ_of_xyz = rhoAB_of_xyz_CubSph;
    pat->dXYZ_dxyz  = drhoAB_dxyz_CubSph;
    break;
  case 2:
    pat->xyz_of_XYZ = xyz_of_rh2AB_CubSph;
    pat->XYZ_of_xyz = rh2AB_of_xyz_CubSph;
    pat->dXYZ_dxyz  = drh2AB_dxyz_CubSph;
    break;
  default:
    errorexiti("stretch=%d is unkonwn", stretch);
  }

  /* set pat remaining parts of CI struct: */
  /* set sigma vars and iSurf, idSurfdX for them */
  if( (stretch==0) && (SigFunc) )
  {
    switch(pat->CI->type)
    {
    case innerCubedSphere:
      /* now set coord. info structure */
      //pat->CI->iFS[0] = isigdef;
      pat->CI->iSurf[0] = isigma0;
      pat->CI->idSurfdX[0][1] = isigma0_dA;
      pat->CI->idSurfdX[0][2] = isigma0_dB;
      /* default surface functions */
      pat->CI->FSurf[0] = FSurf_is_CI_s;
      pat->CI->dFSurfdC[0] = dFSurfdC_is_zero;
      break;

    case outerCubedSphere:
      /* now set coord. info structure */
      //pat->CI->iFS[1] = isigdef;
      pat->CI->iSurf[1] = isigma1;
      pat->CI->idSurfdX[1][1] = isigma1_dA;
      pat->CI->idSurfdX[1][2] = isigma1_dB;
      /* default surface functions */
      pat->CI->FSurf[1] = FSurf_is_CI_s;
      pat->CI->dFSurfdC[1] = dFSurfdC_is_zero;
      break;

    case CubedShell:
      /* now set coord. info structure */
      //pat->CI->iFS[0] = isigdef;
      pat->CI->iSurf[0] = isigma0;
      pat->CI->idSurfdX[0][1] = isigma0_dA;
      pat->CI->idSurfdX[0][2] = isigma0_dB;
      //pat->CI->iFS[1] = isigdef;
      pat->CI->iSurf[1] = isigma1;
      pat->CI->idSurfdX[1][1] = isigma1_dA;
      pat->CI->idSurfdX[1][2] = isigma1_dB;
      /* default surface functions */
      pat->CI->FSurf[0] = FSurf_is_CI_s;
      pat->CI->FSurf[1] = FSurf_is_CI_s;
      pat->CI->dFSurfdC[0] = dFSurfdC_is_zero;
      pat->CI->dFSurfdC[1] = dFSurfdC_is_zero;
      break;

    default:
      errorexit("not sure what to do...");

    /* compute sigma derivs */
//    if(pat->v[isigdef]!=NULL)
//    {
//      // /* deform sigma for testing */
//      // if(1) deform_CubedSphere_sigma01(pat, isigdef, 0.2, -0.1);
//      /* enable isigma and its derivs, and set them  */
//      enablevar_inpatch(pat, isigma);
//      enablevar_inpatch(pat, isigma_dA);
//      enablevar_inpatch(pat, isigma_dB);
//        /* wait until dom=5 and then set all 6 pats from dom0 to dom5 */
//      if(pat->CI->dom==5)
//        init_6CubedSphereBoxes_from_CI_iFS(mesh, p0);
//    }
    }
  }

  /* set label for this patch */
  pat->CI->label = coordinates_get_label(pat);

  /* save stretch coeffs in pat->CI->co[] */
  switch(pat->CI->label)
  {
  case CubedSphere_Stretch2:
    pat->CI->co[0] = Getd(Par("amr_Stretch_w"));
    pat->CI->co[1] = Getd(Par("amr_Stretch_A"));
    break;
  }

  //PRF;printf(" pat->p=%d xyz_of_XYZ=%p\n", pat->p, pat->xyz_of_XYZ);
  return pat->p; /* return index of this patch */
}
