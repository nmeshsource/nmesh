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
int arrange_12CubSph_into_empty_cube(tMesh *mesh, double *xc,
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
  pl = add_6CubedSphere_pats(mesh, outerCubedSphere,0,1, xc, Din,Dmid);
  pl = add_6CubedSphere_pats(mesh, innerCubedSphere,0,1, xc, Dmid,Dout);
  return pl;
}

/* same as arrange_12CubSph_into_empty_cube, but add one cube at the center */
int arrange_1pat12CubSph_into_full_cube(tMesh *mesh, double *xc,
                                        double din, double dmid, double dout)
{
  int pl;
  pl = add_1cube_pat(mesh, xc, din);
  pl = arrange_12CubSph_into_empty_cube(mesh, xc, din,dmid,dout);
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
int two_full_cubes_touching_at_x0(tMesh *mesh, double dc,
                                  double din1, double dmid1,
                                  double din2, double dmid2)
{
  int pl;
  double xc[3];

  /* put centers on x-axis */
  xc[1] = xc[2] = 0.0;

  /* full cube 1 is centered at xc[1]=dc and has width 2dc */
  xc[0] = dc;
  pl = arrange_1pat12CubSph_into_full_cube(mesh, xc, din1,dmid1, dc);

  /* full cube 2 is centered at xc[1]=-dc and has width 2dc */
  xc[0] = -dc;
  pl = arrange_1pat12CubSph_into_full_cube(mesh, xc, din2,dmid2, dc);

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
int sphere_around_two_full_cubes_touching_at_x0(tMesh *mesh,
        double dc, double din1, double dmid1, double din2, double dmid2,
        double r0)
{
  int pl;
  double xc[3], Din[6], Dout[6];
  int i;

  /* make the 2 full cubes */
  pl = two_full_cubes_touching_at_x0(mesh, dc, din1,dmid1, din2,dmid2);

  /* set distances to make 6 more cubed spheres around these 2 full cubes */
  for(i=0; i<6; i++)
  {
    if(i<2) Din[i] = 2.0*dc;
    else    Din[i] = dc;
    Dout[i] = r0;
  }
  xc[0] = xc[1] = xc[2] = 0.0;
  pl = add_6CubedSphere_pats(mesh, outerCubedSphere,0,0, xc, Din,Dout);
  return pl;
}


/* put 6 stretchedCubedShell's around the sphere from
   sphere_around_two_full_cubes_touching_at_x0
                   ___________
             _____/           \______
          __/                        \_
         /                             \
       _- \-                         -/
      /     \-       _______       -/
     |        \-  __/       \__  -/   ...
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
int two_spheres_around_two_full_cubes(tMesh *mesh,
        double dc, double din1, double dmid1, double din2, double dmid2,
        double r0, double r1)
{
  int pl;
  double xc[3], Din[6], Dout[6];
  int i;

  /* make the 2 full cubes and sphere0 around them */
  pl = sphere_around_two_full_cubes_touching_at_x0(mesh, dc,
                                                   din1,dmid1, din2,dmid2, r0);
  /* set distances to make 6 more stretched cubed shells around sphere0 */
  for(i=0; i<6; i++)
  {
    Din[i]  = r0;
    Dout[i] = r1;
  }
  xc[0] = xc[1] = xc[2] = 0.0;
  pl = add_6CubedSphere_pats(mesh, CubedShell,1,1, xc, Din,Dout);
  return pl;
}

/************************************************************************/
/* make patches the surround a central box */
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
int sphere_around_empty_box_at_x0(tMesh *mesh, double dc[3], double r0)
{
  int pl;
  double xc[3], Din[6], Dout[6];
  int f;

  /* set distances to make 6 cubed spheres around the box */
  for(f=0; f<6; f++)
  {
    Din[f]  = dc[f/2];
    Dout[f] = r0;
  }
  xc[0] = xc[1] = xc[2] = 0.0;
  pl = add_6CubedSphere_pats(mesh, outerCubedSphere,0,0, xc, Din,Dout);
  return pl;
}

/* same as sphere_around_empty_box_at_x0, but add one box at the center */
int sphere_around_full_box_at_x0(tMesh *mesh, double dc[3], double r0)
{
  int pl;
  double xc[] = { 0., 0., 0.};
  pl = add_1box_pat(mesh, xc, dc);
  pl = sphere_around_empty_box_at_x0(mesh, dc, r0);
  return pl;
}

/* put 6 stretchedCubedShell's around the sphere from
   sphere_around_full_box_at_x0
                   ___________
             _____/           \______
          __/                        \_
         /                             \
       _- \-                         -/
      /     \-       _______       -/
     |        \-  __/       \__  -/   ...
    |           \/             \/
   /            /__           __\
  |            /   -- _____ --   \    r0 is radius of inner sphere (sphere0)
  |           |      |     |      |   r1 is radius of outer sphere (sphere1)
  |           |      |_____|      |
  |            \ __--       --__ /
   \            \               /
    |          _/\__         __/
     |       _/     \_______/
      \    _/
       -_ /
         \      ...
*/
int two_spheres_around_box_at_x0(tMesh *mesh, double dc[3],
                                 double r0, double r1)
{
  int pl;
  double xc[3], Din[6], Dout[6];
  int i;

  /* make the 2 full cubes and sphere0 around them */
  pl = sphere_around_full_box_at_x0(mesh, dc, r0);

  /* set distances to make 6 more stretched cubed shells around sphere0 */
  for(i=0; i<6; i++)
  {
    Din[i]  = r0;
    Dout[i] = r1;
  }
  xc[0] = xc[1] = xc[2] = 0.0;
  pl = add_6CubedSphere_pats(mesh, CubedShell,1,1, xc, Din,Dout);
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
     -- ___/   dom1 /           r0 is radius of inner sphere
           ---____ /            r1 is radius of outer sphere
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


/* add 1 pat as a box centered at xc[i], returns the index of the pat */
int add_1box_pat(tMesh *mesh, double xc[3], double dout[3])
{
  int amr_n = Geti(Par("amr_n"));
  int n1max = Geti(Par("amr_nmax"));
  int n[3] = { amr_n,amr_n,amr_n };
  double bbox[6];
  tPat *pat;
  int d;

  /* set min/max in each direction */
  for(d=0; d<3; d++)
  {
    bbox[2*d]   = xc[d] - dout[d];
    bbox[2*d+1] = xc[d] + dout[d];
  }

  /* make new patch */
  pat = add_patch(mesh, bbox, n, n1max);

  return pat->p; /* return pat index */
}

/* add 1 pat as a cube centered at xc[i], returns the index of the pat */
int add_1cube_pat(tMesh *mesh, double *xc, double rout)
{
  double dout[] = { rout, rout, rout };
  return add_1box_pat(mesh, xc, dout);
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
  int amr_n = Geti(Par("amr_n"));
  int n1max = Geti(Par("amr_nmax"));
  int n[3] = { amr_n,amr_n,amr_n };
  double bbox[6];
  tPat *pat;
  int isigma0    = Ind("CubedSphere_sigma0");
  int isigma0_dA = Ind("CubedSphere_dsigma0_dA");
  int isigma0_dB = Ind("CubedSphere_dsigma0_dB");
  int isigma1    = Ind("CubedSphere_sigma1");
  int isigma1_dA = Ind("CubedSphere_dsigma1_dA");
  int isigma1_dB = Ind("CubedSphere_dsigma1_dB");
  //int isigdef   = Ind("CubedSphere_sigma01_def");
  int i;

  if(N<1 || N>6) errorexit("N must be 1,2,3,4,5,6");

  /* set pat CI struct */
  for(i=0; i<N; i++)
  {
    int d;
    double Amin,Amax, Bmin,Bmax;

    /* set min/max in each direction */
    set_AB_min_max_from_Din(i, Din, &Amin,&Amax, &Bmin,&Bmax);
    bbox[0] = 0.;
    bbox[1] = 1.;
    bbox[2] = Amin;
    bbox[3] = Amax;
    bbox[4] = Bmin;
    bbox[5] = Bmax;

    /* make new patch */
    pat = add_patch(mesh, bbox, n, n1max);

    /* set coords trafos */
    if(stretch==0)
    {
      pat->xyz_of_XYZ = xyz_of_lamAB_CubSph;
      pat->XYZ_of_xyz = lamAB_of_xyz_CubSph;
      pat->dXYZ_dxyz  = dlamAB_dxyz_CubSph;
    }
    else
    {
      pat->xyz_of_XYZ = xyz_of_rhoAB_CubSph;
      pat->XYZ_of_xyz = rhoAB_of_xyz_CubSph;
      pat->dXYZ_dxyz  = drhoAB_dxyz_CubSph;
    }

    /* set center */
    for(d=0; d<3; d++)  pat->CI->xc[d] = xc[d];

    /* set inner and outer distance from center */
    pat->CI->s[0] = Din[i];
    pat->CI->s[1] = Dout[i];

    /* set domain index and type */
    pat->CI->dom = i;
    pat->CI->type= type;

    /* set sigma vars and iSurf, idSurfdX for them */
    if( (stretch==0) && (SigFunc) )
    {
      switch(type)
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
      }
      /* compute sigma derivs */
//      if(pat->v[isigdef]!=NULL)
//      {
//        // /* deform sigma for testing */
//        // if(1) deform_CubedSphere_sigma01(pat, isigdef, 0.2, -0.1);
//        /* enable isigma and its derivs, and set them  */
//        enablevar_inpatch(pat, isigma);
//        enablevar_inpatch(pat, isigma_dA);
//        enablevar_inpatch(pat, isigma_dB);
//          /* wait until dom=5 and then set all 6 pats from dom0 to dom5 */
//        if(pat->CI->dom==5)
//          init_6CubedSphereBoxes_from_CI_iFS(mesh, p0);
//      }
    }
  }

  return pat->p; /* return pat index of last added pat */
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

/* add 1 cubed sphere with specific bounds */
int add_1_CubedSphere_pat(tMesh *mesh, int dom, int type,
                          int stretch, int SigFunc, double *xc,
                          double Din, double Dout, double ABrct[4])
{
  int amr_n = Geti(Par("amr_n"));
  int n1max = Geti(Par("amr_nmax"));
  int n[3] = { amr_n,amr_n,amr_n };
  double bbox[6];
  tPat *pat;
  int isigma0    = Ind("CubedSphere_sigma0");
  int isigma0_dA = Ind("CubedSphere_dsigma0_dA");
  int isigma0_dB = Ind("CubedSphere_dsigma0_dB");
  int isigma1    = Ind("CubedSphere_sigma1");
  int isigma1_dA = Ind("CubedSphere_dsigma1_dA");
  int isigma1_dB = Ind("CubedSphere_dsigma1_dB");
  int d;

  /* set min/max in each direction */
  bbox[0] = 0.;
  bbox[1] = 1.;
  bbox[2] = ABrct[0];
  bbox[3] = ABrct[1];;
  bbox[4] = ABrct[2];
  bbox[5] = ABrct[3];

  /* make new patch */
  pat = add_patch(mesh, bbox, n, n1max);

  /* set pat CI struct */
  /* set coords trafos */
  if(stretch==0)
  {
    pat->xyz_of_XYZ = xyz_of_lamAB_CubSph;
    pat->XYZ_of_xyz = lamAB_of_xyz_CubSph;
    pat->dXYZ_dxyz  = dlamAB_dxyz_CubSph;
  }
  else
  {
    pat->xyz_of_XYZ = xyz_of_rhoAB_CubSph;
    pat->XYZ_of_xyz = rhoAB_of_xyz_CubSph;
    pat->dXYZ_dxyz  = drhoAB_dxyz_CubSph;
  }

  /* set center */
  for(d=0; d<3; d++)  pat->CI->xc[d] = xc[d];

  /* set inner and outer distance from center */
  pat->CI->s[0] = Din;
  pat->CI->s[1] = Dout;

  /* set domain index and type */
  pat->CI->dom = dom;
  pat->CI->type= type;

  /* set sigma vars and iSurf, idSurfdX for them */
  if( (stretch==0) && (SigFunc) )
  {
    switch(type)
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
    }
  }
  return pat->p; /* return pat index of last added pat */
}
