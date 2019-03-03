/* derivs.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "coordinates.h"




/* compute Cart. derivs, put du/dx^m into vars with index dui */
int cart_partials(tNode *node, int ui, int dui)
{
  tPat *pat = node->pat;
  tMesh *mesh = pat->mesh;
  tDat *dat = node->dat;
  double *du[] = { Vard(node,dui), Vard(node,dui+1), Vard(node,dui+2) };
  double dXbdX[3];
  int ret, ind, m,i;

  if(!dat) return 0;

  /* do we need to init. coords? */
  if(!(dat->coords_set)) coordinates_init_node(node);

  /* take derivs with respect to Xb: du/dXb */
  ret = basis_var_derivs(node, ui, dui);

  /* get dXb/dX */
  dXbYbZb_dXYZ(node, dXbdX);

  /* scale: du/dX = dXb/dX du/dXb = */
  forpoints(node,ind)
    for(m=0; m<3; m++) du[m][ind] *= dXbdX[m];

  /* transform to Cartesian coords */
  if(pat->dXYZ_dxyz)
  {
    int idXd = Ind("dXdx");
    double *dXdx[3][3]
              = { {Vard(node,idXd),   Vard(node,idXd+1), Vard(node,idXd+2)},
                  {Vard(node,idXd+3), Vard(node,idXd+4), Vard(node,idXd+5)},
                  {Vard(node,idXd+6), Vard(node,idXd+7), Vard(node,idXd+8)} };
    /* compute Cartesian derivs at all points */
    forpoints(node,ind)
    {
      double dv[3];

      /* Transform derivs to Cartesian coords */
      for(m=0; m<3; m++)
      {
        dv[m] = 0.;
        for(i=0; i<3; i++) dv[m] += dXdx[i][m][ind] * du[i][ind];
      }
      /* copy dv into du */
      for(m=0; m<3; m++) du[m][ind] = dv[m];
    }
  }
  return ret;
}
