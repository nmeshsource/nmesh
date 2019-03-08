/* setup_CubedSpheres.c */
/* Wolfgang Tichy, Nov 2017 */
/* create various arrangements of cubed spheres */

#include "nmesh.h"
#include "coordinates.h"



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
int arrange_12CubSph_into_empty_cube(tMesh *mesh, int b0, double *xc,
                                     double din, double dmid, double dout)
{
  int bl=b0;
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
  bl = convert_6pats_to_CubedSphere(mesh, bl, outerCubedSphere,0,
                                    xc, Din,Dmid);
  bl = convert_6pats_to_CubedSphere(mesh, bl, innerCubedSphere,0,
                                    xc, Dmid,Dout);
  return bl;
}

/* same as arrange_12CubSph_into_empty_cube, but add one cube at the center */
int arrange_1pat12CubSph_into_full_cube(tMesh *mesh, int b0, double *xc,
                                        double din, double dmid, double dout)
{
  int bl=b0;
  bl = convert_1pat_to_cube(mesh, bl, xc, din);
  bl = arrange_12CubSph_into_empty_cube(mesh, bl, xc, din,dmid,dout);
  return bl;
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
int two_full_cubes_touching_at_x0(tMesh *mesh, int b0, double dc,
                                  double din1, double dmid1,
                                  double din2, double dmid2)
{
  int bl=b0;
  double xc[4];

  /* put centers on x-axis */
  xc[2] = xc[3] = 0.0;

  /* full cube 1 is centered at xc[1]=dc and has width 2dc */
  xc[1] = dc;
  bl = arrange_1pat12CubSph_into_full_cube(mesh, bl, xc, din1,dmid1, dc);

  /* full cube 2 is centered at xc[1]=-dc and has width 2dc */
  xc[1] = -dc;
  bl = arrange_1pat12CubSph_into_full_cube(mesh, bl, xc, din2,dmid2, dc);

  return bl;
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
int sphere_around_two_full_cubes_touching_at_x0(tMesh *mesh, int b0,
        double dc, double din1, double dmid1, double din2, double dmid2,
        double r0)
{
  int bl=b0;
  double xc[4], Din[6], Dout[6];
  int i;

  /* make the 2 full cubes */
  bl = two_full_cubes_touching_at_x0(mesh, b0, dc, din1,dmid1, din2,dmid2);

  /* set distances to make 6 more cubed spheres around these 2 full cubes */
  for(i=0; i<6; i++)
  {
    if(i<2) Din[i] = 2.0*dc;
    else    Din[i] = dc;
    Dout[i] = r0;
  }  
  xc[1] = xc[2] = xc[3] = 0.0;
  bl = convert_6pats_to_CubedSphere(mesh, bl, outerCubedSphere,0,
                                     xc, Din,Dout);
  return bl;
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
int two_spheres_around_two_full_cubes(tMesh *mesh, int b0,
        double dc, double din1, double dmid1, double din2, double dmid2,
        double r0, double r1)
{
  int bl=b0;
  double xc[4], Din[6], Dout[6];
  int i;

  /* make the 2 full cubes and sphere0 around them */
  bl = sphere_around_two_full_cubes_touching_at_x0(mesh, b0, dc,
                                                   din1,dmid1, din2,dmid2, r0);
  /* set distances to make 6 more stretched cubed shells around sphere0 */
  for(i=0; i<6; i++)
  {
    Din[i]  = r0;
    Dout[i] = r1;
  }  
  xc[1] = xc[2] = xc[3] = 0.0;
  bl = convert_6pats_to_CubedSphere(mesh, bl, CubedShell,1, xc, Din,Dout);
  return bl;
}



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


/* convert 1 pat to a cube centered at xc[i],
   returns the index of the pat right after the last converted pat */
int convert_1pat_to_cube(tMesh *mesh, int b0, double *xc, double dout)
{
  tPat *pat = mesh->pat[b0];
  int d;

  if(b0+1 > mesh->npats)
  {
    printf("npats is too small\n");
    return b0;
  }

  /* set min/max in each direction */
  for(d=0; d<3; d++)
  {
    pat->bbox[2*d]   = xc[d] - dout;
    pat->bbox[2*d+1] = xc[d] + dout;
  }

  ///* erase all bface info in this pat */
  //remove_all_bfaces(pat);

  return b0+1; /* return pat index right after last added pat */
}


/* convert 6 pats starting with p0 to some kind of cubed spheres */
/* call this after all pats exist already, so at POST_GRID */
/* type can be "PyramidFrustum", "innerCubedSphere", "outerCubedSphere",
   "CubedShell"
   xc[1..3] = (x,y,z) of coord center for the 6 cubed spheres
   Din[0...5] inner distance from center for cubed sph. domain 0-5
   Dout[0...5] outer distance from center for cubed sph. domain 0-5 */
/* It returns the index of the pat right after the last converted pat */
int convert_6pats_to_CubedSphere(tMesh *mesh, int p0, int type, int stretch,
                                  double *xc, double *Din, double *Dout)
{
  int isigma    = Ind("CubedSphere_sigma01");
  int isigma_dA = Ind("CubedSphere_dsigma01_dA");
  int isigma_dB = Ind("CubedSphere_dsigma01_dB");
  int isigdef   = Ind("CubedSphere_sigma01_def");
  int i;

  if(p0+6 > mesh->npats)
    printf("npats is too small\n");

  /* set pat CI struct */
  for(i=0; i<6; i++)
  {
    tPat *pat;
    int d;
    double Amin,Amax, Bmin,Bmax;

    /* break i-loop if we do not have enough pats */
    if(p0+i >= mesh->npats) break;
    pat = mesh->pat[p0+i];

    /* set min/max in each direction */
    set_AB_min_max_from_Din(i, Din, &Amin,&Amax, &Bmin,&Bmax);
    pat->bbox[0] = 0.;
    pat->bbox[1] = 1.;
    pat->bbox[2] = Amin;
    pat->bbox[3] = Amax;
    pat->bbox[4] = Bmin;
    pat->bbox[5] = Bmax;

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

    /* default surface functions */
    pat->CI->FSurf[0] = FSurf_is_CI_s;
    pat->CI->FSurf[1] = FSurf_is_CI_s;
    pat->CI->dFSurfdC[0] = dFSurfdC_is_zero;
    pat->CI->dFSurfdC[1] = dFSurfdC_is_zero;

    /* set sigma vars and iSurf, idSurfdX for them */
    if( (stretch==0) )
    {
      switch(type)
      {
        case innerCubedSphere:
          /* now set coord. info structure */
          pat->CI->iFS[0] = isigdef;
          pat->CI->iSurf[0] = isigma;
          pat->CI->idSurfdX[0][1] = isigma_dA;
          pat->CI->idSurfdX[0][2] = isigma_dB;
          break;
      
        case outerCubedSphere:
          /* now set coord. info structure */
          pat->CI->iFS[1] = isigdef;
          pat->CI->iSurf[1] = isigma;
          pat->CI->idSurfdX[1][1] = isigma_dA;
          pat->CI->idSurfdX[1][2] = isigma_dB;
          break;
      
        case CubedShell:
          /* now set coord. info structure */
          pat->CI->iFS[0] = isigdef;
          pat->CI->iSurf[0] = isigma;
          pat->CI->idSurfdX[0][1] = isigma_dA;
          pat->CI->idSurfdX[0][2] = isigma_dB;
          pat->CI->iFS[1] = isigdef;
          pat->CI->iSurf[1] = isigma;
          pat->CI->idSurfdX[1][1] = isigma_dA;
          pat->CI->idSurfdX[1][2] = isigma_dB;
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

    ///* erase all bface info in this pat */
    //remove_all_bfaces(pat);
  }

  return p0+i; /* return pat index right after last added pat */
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
