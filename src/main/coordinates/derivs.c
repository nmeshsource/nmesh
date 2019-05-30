/* derivs.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "coordinates.h"




/* compute Cart. derivs, put du/dx^m into vars with index dui[0..2] */
int cart_partials(tNode *node, int ui, int dui[3])
{
  tPat *pat = node->pat;
  tMesh *mesh = pat->mesh;
  tDat *dat = node->dat;
  double *du[] = { Vard(node,dui[0]), Vard(node,dui[1]), Vard(node,dui[2]) };
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


/* compute first derivs U_{,i} of a scalar U in a node */
void cart_partials_U(tNode *node, int U, int dUi)
{
  int dU[] = { dUi, dUi+1, dUi+2 };

  /* compute partial derivs of all components in node */
  cart_partials(node, U, dU);
}


/* compute first derivs V_{i,j} of a vector V_{i} in a node */
void cart_partials_Vi(tNode *node, int Vi, int dVij)
{
  /* compute partial derivs of all components in node */
  cart_partials_U(node, Vi,   dVij);
  cart_partials_U(node, Vi+1, dVij+3);
  cart_partials_U(node, Vi+2, dVij+6);
}

/* compute first derivs S_{ij,k} of a symmetric tensor S_{ij} in a node */
void cart_partials_Sij(tNode *node, int Sij, int dSijk)
{
  int n;

  /* compute partial derivs of all components in node */
  for(n=0; n<6; n++)
    cart_partials_U(node, Sij + n, dSijk + 3*n);
}

/* compute first derivs T_{ij,k} of a general tensor T_{ij} in a node */
void cart_partials_Tij(tNode *node, int Tij, int dTijk)
{
  int n;

  /* compute partial derivs of all components in node */
  for(n=0; n<6; n++)
    cart_partials_U(node, Tij + n,   dTijk + 3*n);
}

/* compute 1st and 2nd order Cart. derivs of scalar U */
void cart_partials2_U(tNode *node, int U, int dUi, int ddUij)
{
  /* 1st derivs */
  cart_partials_U(node, U, dUi);

  /* 2nd derivs */
  cart_partials_Vi(node, dUi, ddUij);
}

/* compute 1st and 2nd derivs V_{i,jk} of a vector V_{i} in a node */
void cart_partials2_Vi(tNode *node, int Vi, int dVij, int ddVijk)
{
  /* 1st derivs */
  cart_partials_Vi(node, Vi, dVij);

  /* 2nd derivs */
  cart_partials_Tij(node, dVij, ddVijk);
}

/* compute 1st and 2nd derivs S_{ij,kl} of a symmetric tensor S_{ij} */
void cart_partials2_Sij(tNode *node, int Sij, int dSijk, int ddSijkl)
{
  int n;

  /* 1st derivs */
  cart_partials_Sij(node, Sij, dSijk);

  /* 2nd derivs */
  cart_partials_U(node, dSijk,   ddSijkl);
  cart_partials_U(node, dSijk,   ddSijkl);
  for(n=0; n<18; n++)
    cart_partials_U(node, dSijk + n, ddSijkl + 3*n);


  
  cart_partials_Sij(node, dSijk,   ddSijkl);
  cart_partials_Sij(node, dSijk+3, ddSijkl);
  cart_partials_Sij(node, dSijk+6, ddSijkl);









  /* compute partial derivs of all components in node */
  cart_partials_U(node, Sij,   dSijk);
  cart_partials_U(node, Sij+1, dSijk+3);
  cart_partials_U(node, Sij+2, dSijk+6);
  cart_partials_U(node, Sij+3, dSijk+9);
  cart_partials_U(node, Sij+4, dSijk+12);
  cart_partials_U(node, Sij+5, dSijk+15);
}
