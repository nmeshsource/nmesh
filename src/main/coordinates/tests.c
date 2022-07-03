/* tests.c */
/* Wolfgang Tichy, 7/2022 */

#include "nmesh.h"
#include "coordinates.h"


/* global vars */
extern tcoordinates coordinates[1];




/* Compute \partial_{\bar{i}} (J \sqrt{g^{\bar{i}\bar{i}}} n^{\bar{i}}_i)
   for each i. Save in var divb_J_sqrtgdiag_nx.
   This should converge to zero, if sqrtgdiag is obtained by transforming
   the flat metric. */
int coordinates_set_divb_J_sqrtgdiag_n(tNode *node)
{
  tPat *pat = node->pat;
  tMesh *mesh = pat->mesh;
  //tDat *dat = node->dat;

  //int midpoint_data = coordinates->midpoint_data;
  //int sqrtdet2g_o_det3gamma = coordinates->sqrtdet2g_o_det3gamma;
  //int idXdx = coordinates->idXdx;
  int itmp1 = coordinates->itmp1;
  int itmp2 = itmp1 + 1;
  int idivb_J_sqrtgdiag_nx = Ind("divb_J_sqrtgdiag_nx");
  int idet_dXbdx = coordinates->idet_dXbdx;
  //int isqrtdet2g_o_det3gamma0 = coordinates->isqrtdet2g_o_det3gamma0;
  int isqrtgdiagx = coordinates->isqrtgdiagx;
  //int iXm_det_dXbdx = coordinates->iXm_det_dXbdx;
  //int iYm_det_dXbdx = coordinates->iYm_det_dXbdx;
  //int iZm_det_dXbdx = coordinates->iZm_det_dXbdx;
  //int iXm_sqrtgdiagx = coordinates->iXm_sqrtgdiagx;
  //int iYm_sqrtgdiagy = coordinates->iYm_sqrtgdiagy;
  //int iZm_sqrtgdiagz = coordinates->iZm_sqrtgdiagz;

  tVarList *vltmp2 = vlalloc(mesh);
  tVarList *vldivb_J_sqrtgdiag_n = vlalloc(mesh);
  int ib, i, ijk;

  vlpush(vltmp2, itmp2);

  /* loop over i of n^{\bar{i}_i */
  for(i=0; i<3; i++)
  {
    vlpushone(vldivb_J_sqrtgdiag_n, idivb_J_sqrtgdiag_nx+i);
    vlsetconstant_node(node, vldivb_J_sqrtgdiag_n, 0.);

    /* loop over ib=\bar{i} of n^{\bar{i}_i */
    for(ib=0; ib<3; ib++)
    {
      double *tmp1 = Vard(node, itmp1);
      double *sqrtgdiag = Vard(node, isqrtgdiagx+ib);
      double *ooJ = Vard(node, idet_dXbdx);

      /* set tmp1 = \sqrt{g^{\bar{i}\bar{i}}} J * n^{\bar{i}}_i
         for each \bar{i}=ib */
      forpoints(node,ijk)
      {
        int f = 2*ib + 1; /* results in right pointing normal */
        double norm[3];
        node_normal_at_ijk(node, f, ijk, norm);
        tmp1[ijk] = (sqrtgdiag[ijk] / ooJ[ijk]) * norm[i];
      }

      /* tmp2 = \partial_{\bar{i}} tmp1 */
      basis_var_deriv1(node, ib, itmp1, itmp2, NULL);
      /* divb_J_sqrtgdiag_n += tmp2 */
      vladdto_node(node, vldivb_J_sqrtgdiag_n, 1., vltmp2);
    }
    vldropn(vldivb_J_sqrtgdiag_n, 1);
  }

  vlfree(vldivb_J_sqrtgdiag_n);
  vlfree(vltmp2);
  return 0;
}

/* run the tests */
int coordinates_tests(tMesh *mesh)
{
  enablevar(mesh, Ind("divb_J_sqrtgdiag_nx"));

  formylnodes(mesh)
  {
    tNode *node = MyLnode;

    coordinates_set_divb_J_sqrtgdiag_n(node);
  }
  return 0;
}
