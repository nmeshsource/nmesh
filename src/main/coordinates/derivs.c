/* derivs.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "coordinates.h"


#define Vd(node,i) GetVarDpointer(node,i)


/* compute Cart. derivs, put du/dx^m into vars with index dui */
int cart_partials(tNode *node, int ui, int dui)
{
  tPat *pat = node->pat;
  tMesh *mesh = pat->mesh;
  double *du[] = { Vd(node,dui), Vd(node,dui+1), Vd(node,dui+2) };
  double dXbdX[3];
  int ret, ind, m,i;

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
              = { {Vd(node,idXd),   Vd(node,idXd+1), Vd(node,idXd+2)},
                  {Vd(node,idXd+3), Vd(node,idXd+4), Vd(node,idXd+5)},
                  {Vd(node,idXd+6), Vd(node,idXd+7), Vd(node,idXd+8)} };
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
