/* fv.c */
/* Wolfgang Tichy, July 2022 */


#include "nmesh.h"
#include "dg.h"


/*********************************************************************/
/* these fluxes are equal to the normal n_x, n_y, or nz */
/*********************************************************************/

/* use normal in x,y,z dir as test flux */
void testflux_normi(tDGinfo *d, int i)
{
  double norm[3];
  int l;

  node_normal_from_DGinfo(d, norm);
  forvl(d->vlu, l)
    d->fi[l] = d->fa[l] = d->fnum[l] = norm[i];
}

/* use testflux_normi in x-dir */
void test_flux_normx(tDGinfo *d)
{
  testflux_normi(d, 0);
}
void test_rec_flux_normx(tFVinfo *f, tDGinfo *d)
{
  testflux_normi(d, 0);
}

/* use testflux_normi in y-dir */
void test_flux_normy(tDGinfo *d)
{
  testflux_normi(d, 1);
}
void test_rec_flux_normy(tFVinfo *f, tDGinfo *d)
{
  testflux_normi(d, 1);
}

/* use testflux_normi in z-dir */
void test_flux_normz(tDGinfo *d)
{
  testflux_normi(d, 2);
}
void test_rec_flux_normz(tFVinfo *f, tDGinfo *d)
{
  testflux_normi(d, 2);
}


/*********************************************************************/
/* use test fluxes */
/*********************************************************************/

/* Compute div(flux) for the 3 test fluxes above. This must converge to 0
   with 2nd order in the interior. At the boundary it converges with 1st
   order if fv_divf_extrap=no, and 2nd order if fv_divf_extrap=dnfn_extrap1.
   The 3 fluxes are simply the normals in the 3 X-coord dirs. */
int fv_test_fv_divf(tNode *node, tVarList *vldivf)
{
  /* func pointers */
  void (*test_rec_flux[3])(tFVinfo *f, tDGinfo *d) =
       {test_rec_flux_normx, test_rec_flux_normy, test_rec_flux_normz};
  void (*test_flux[3])(tDGinfo *d) =
       {test_flux_normx, test_flux_normy, test_flux_normz};
  tMesh *mesh = node->pat->mesh;
  tVarList *vldivf_i = vlalloc(mesh);
  int i;

  /* loop over the three test fluxes */
  for(i=0; i<3; i++)
  {
    tVarList *vlu = vldivf_i; /* since testflux_normi does not need vlu */

    vlpushone(vldivf_i, Vind(vldivf, i));

    fv_divf(node, vldivf_i, vlu,vlu,NULL,
            test_rec_flux[i], test_flux[i], test_flux[i]);
    dg_add_surface_fluxes_sign(node, 1., vldivf, vlu, NULL,
                               test_flux[i], test_flux[i]);

    vldropn(vldivf_i, 0);
  }
  vlfree(vldivf_i);
  return 0;
}


/*************************************************************************/
/* run all the tests */
/*************************************************************************/

/* run the tests */
int fv_tests(tMesh *mesh)
{
  tVarList *vldivf = vlalloc(mesh);

  enablevar(mesh, Ind("fv_test_divfx"));
  vlpush(vldivf, Ind("fv_test_divfx"));

  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    fv_test_fv_divf(node, vldivf);
  }
  vlfree(vldivf);
  return 0;
}
