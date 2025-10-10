/* penalty.c */
/* Wolfgang Tichy, Oct 2025 */


#include "nmesh.h"
#include "dg.h"


/* get global amr and DG vars */
extern tAMR amr[1];
extern tDGglobals DGglobals[1];


/***************************************************************************/
/* some functions to add penalty terms on boundary */
/***************************************************************************/

/* Add penalty terms with a choice of sign (sign=+1 or sign=-1) to vlr.
   We have the flag use_fv, to decide whether fv mode is active at
   all inside this function. */
int dg_add_penalty_sign_fvflag(tElm *elm, double sign,
                               tVarList *vlr, tVarList *vlu, tVarList *vls,
                               void (*u_fnum)(tDGinfo *d),
                               int use_fv)
{
  tDGinfo *dgi = alloc_DGinfo(vlu, vls);
  tMesh *mesh = vlu->mesh;
  int skip_fv = DGglobals->fv_divf_adds_surface_fluxes;
  int add_surface_fluxes = 1; /* by default we want to set fluxes here */
  int useLGL_Wq = 0;
  double distXb[6] = {0};

  /* we now call get_all_myln_surfaces in evolve_setrhs_mesh
     so we do not need to do it here */
  ///* get surfaces so that we can compute fluxes */
  //get_all_myln_surfaces(mesh);

  /* special fv cases */
  if(use_fv)
  {
    /* do nothing if fv_divf has already taken care of surface fluxes */
    if(skip_fv)
    {
      add_surface_fluxes = 0;
    }
    else
    {
      /* find distance from faces to nearest midpoint, and don't use LGL Wq */
      set_nodemidpoints_to_face_distXb(elm, distXb);
      useLGL_Wq = 0;
    }
  }

  /* add boundary flux terms */
  if(add_surface_fluxes)
  {
    int *n = elm->n;
    int face;
    int not_fv = !use_fv;
    //int have_Xofx = ( elm->pat->XYZ_of_xyz ? 1 : 0 );
    double mod0 = not_fv;     /* set to 1 if we don't use fin. vol. */
    double mod1 = 1. - mod0;  /* set to 1 if we use fin. vol. */

    /* set DG node info */
    dgi->node = elm;

    for(face=0; face<6; face++)
    {
      int dir = face/2;
      int N = n[dir] - 1;
      int p = (face%2)*N;
      double twooLX = 2./(elm->bbox[2*dir+1] - elm->bbox[2*dir]);
      double *WQ = Wquad(elm, dir);
      double Wq = (useLGL_Wq) ? 2./(N*(N+1.)) : WQ[p];
      double Wqmod = fabs(distXb[face]);
      double oow = 1./(Wq*mod0 + Wqmod*mod1);
      int i,j,k;

      /* do nothing if dir is not active */
      if(!Getb(amr->dir_active[dir])) continue;

      /* set DG face info */
      dgi->face = face;

      forplaneN(dir, i,j,k, n, p)
      {
        int ijk = Ind_n(i,j,k, n);
        double pen = twooLX * oow; /* same as bamps' bndy->penalty */
        double Ffac;
        int l;

        /* set DG i,j,k info */
        dgi->i = i;
        dgi->j = j;
        dgi->k = k;
        /* we do not use:  dgi->info = use_fv;
           because our regular grid points on the faces already are
           considered midpoints in fin. vol. approach */

        /* set vars on both sides, set correction to rhs in fnum */
        u_fnum(dgi);
        /* Note: in u_fnum we need to use normal vectors v_i that are just
                 v_i = dX/dx^i
           So v_i is not normalized with the flat metric (see GH_penalty_pt).
           Then we must not scale twooLX with dg_scale_penalty_bamps(dgi),
           since our radial normal vectors already contain a factor of
           dlam/dr = 1/(CI->s[1] - CI->s[0]) */

        /* get Ffac, this can be set in u_f_lam or numflux */
        Ffac = dgi->Ffac; /* usually 1, set to 0 to turn off surface fluxes */

        /* get F from dgi and add penalty terms to vlr */
        forvl(vlr, l)
        {
          int ir = Vind(vlr,l);
          double *r = Vard_(elm, ir);
          double F;

          //if(i==2 && j==0 && k==0 &&  elm->eploc->eid==8  &&  l==10)
          //{
          //  printf("0:Ffac=%g sign=%g oow=%g pen=%g: r[ijk]=%.16g\n",
          //         Ffac, sign, oow, pen, r[ijk]);
          //}

          F = (dgi->fnum[l]) * Ffac;
          r[ijk] += F * sign * pen;

          //if(i==2 && j==0 && k==0 &&  elm->eploc->eid==8  &&  l==10)
          //{
          //  printf("1:F=%g F*sign*pen=%g: r[ijk]=%.16g\n",
          //         F, F*sign*pen, r[ijk]);
          //}
        }
      }
    } /* end loop over faces */
  }

  free_DGinfo(dgi);
  return 0;
}


/* use smaller penalty in radial direction as in bamps */
double dg_scale_penalty_bamps(tDGinfo *dgi)
{
  tElm *elm = dgi->node;
  tMesh *mesh = Elm_mesh(elm);
  int pi = elm->eploc->p;
  tPat *pat = mesh->pat[pi];
  tCoordInfo *CI = pat->CI;
  double scale = 1.;

  switch(CI->type)
  {
  case innerCubedSphere:
  case outerCubedSphere:
  case CubedShell:
    if(dgi->face <= 1) scale *= 1./(CI->s[1] - CI->s[0]);
  }

  //printf("scale=%g\n", scale);
  return scale;
}
