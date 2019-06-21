/* derivs.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "coordinates.h"




/* compute Cart. derivs of u in array au, put du/dx^m into arrays dau[0..2] */
int array_cart_partials(tNode *node, tArray *au, tArray *dau[3])
{
  tPat *pat = node->pat;
  tMesh *mesh = pat->mesh;
  tDat *dat = node->dat;
  double *du[] = { Arrd(dau[0]), Arrd(dau[1]), Arrd(dau[2]) };
  double dXbdX[3];
  int ind, m,i;

  if(!dat) return 0;

  /* do we need to init. coords? */
  if(!(dat->coords_set)) coordinates_init_node(node);

  /* take derivs with respect to Xb: du/dXb */
  basis_array_derivs(node, au, dau);

  /* get dXb/dX */
  dXbYbZb_dXYZ(node, dXbdX);

  /* scale: du/dX = dXb/dX du/dXb */
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
  return 1;
}


/* compute Cart. derivs, put du/dx^m into vars with index dui[0..2] */
int cart_partials(tNode *node, int ui, int dui[3])
{
  tDat *dat = node->dat;
  tArray *au;
  tArray *dau[3];

  if(!dat) return 0;

  au     = dat->v[ui];
  dau[0] = dat->v[dui[0]];
  dau[1] = dat->v[dui[1]];
  dau[2] = dat->v[dui[2]];

  return array_cart_partials(node, au, dau);
}

/***********************************************************************/
/* 2 variants of cart_partials to get the 1st derivs of a scaler */
/***********************************************************************/

/* compute first derivs U_{,i} of a scalar U in a node */
void cart_3partials(tNode *node, int U, int dUx, int dUy, int dUz)
{
  int dU[] = { dUx, dUy, dUz };
  cart_partials(node, U, dU);
}

/* compute first derivs U_{,i} of a scalar U in a node */
void cart_partials_U(tNode *node, int U, int dUx)
{
  int dU[] = { dUx, dUx+1, dUx+2 };
  cart_partials(node, U, dU);
}

/***********************************************************************/
/* 1st derivs of vectors and tensors */
/***********************************************************************/

/* compute first derivs U_{i,j} of a vector U_{i} in a node */
void cart_partials_Ui(tNode *node, int Ux, int dUxx)
{
  /* compute partial derivs of all components in node */
  cart_partials_U(node, Ux,   dUxx);
  cart_partials_U(node, Ux+1, dUxx+3);
  cart_partials_U(node, Ux+2, dUxx+6);
}

/* compute first derivs S_{ij,k} of a symmetric tensor S_{ij} in a node */
void cart_partials_Sij(tNode *node, int Sxx, int dSxxx)
{
  int n;

  /* compute partial derivs of all components in node */
  for(n=0; n<6; n++)
    cart_partials_U(node, Sxx + n, dSxxx + 3*n);
}

/* compute first derivs U_{ij,k} of a general tensor U_{ij} in a node */
void cart_partials_Uij(tNode *node, int Uxx, int dUxxx)
{
  int n;

  /* compute partial derivs of all components in node */
  for(n=0; n<9; n++)
    cart_partials_U(node, Uxx + n, dUxxx + 3*n);
}

/***********************************************************************/
/* 2nd derivs of vectors and tensors */
/***********************************************************************/

/* compute 1st and 2nd order Cart. derivs of scalar U */
void cart_partials2_U(tNode *node, int U, int dUx, int ddUxx)
{
  /* 1st derivs */
  cart_partials_U(node, U, dUx);

  /* 2nd derivs */
  cart_3partials(node, dUx,   ddUxx,  ddUxx+1,ddUxx+2);
  cart_3partials(node, dUx+1, ddUxx+1,ddUxx+3,ddUxx+4);
  cart_3partials(node, dUx+2, ddUxx+2,ddUxx+4,ddUxx+5);
}

/* compute 1st and 2nd derivs U_{i,jk} of a vector U_{i} in a node */
void cart_partials2_Ui(tNode *node, int Ux, int dUxx, int ddUxxx)
{
  int n;

  /* 1st derivs */
  cart_partials_Ui(node, Ux, dUxx);

  /* 2nd derivs */
  for(n=0; n<3; n++) /* n=0: dUxj, n=1: dUyj, n=2: dUzj */
  {
    cart_3partials(node, dUxx+3*n ,  ddUxxx+6*n ,  ddUxxx+6*n+1, ddUxxx+6*n+2);
    cart_3partials(node, dUxx+3*n+1, ddUxxx+6*n+1, ddUxxx+6*n+3, ddUxxx+6*n+4);
    cart_3partials(node, dUxx+3*n+2, ddUxxx+6*n+2, ddUxxx+6*n+4, ddUxxx+6*n+5);
  }
}

/* compute 1st and 2nd derivs S_{ij,kl} of a symmetric tensor S_{ij} */
void cart_partials2_Sij(tNode *node, int Sxx, int dSxxx, int ddSxxxx)
{
  int n;

  /* 1st derivs */
  cart_partials_Sij(node, Sxx, dSxxx);

  /* 2nd derivs */
  for(n=0; n<6; n++) /* n=0: dSxxj, n=2: dSxzj, n=3: dSxyj */
  {
    cart_3partials(node, dSxxx+3*n ,  ddSxxxx+6*n ,  ddSxxxx+6*n+1, ddSxxxx+6*n+2);
    cart_3partials(node, dSxxx+3*n+1, ddSxxxx+6*n+1, ddSxxxx+6*n+3, ddSxxxx+6*n+4);
    cart_3partials(node, dSxxx+3*n+2, ddSxxxx+6*n+2, ddSxxxx+6*n+4, ddSxxxx+6*n+5);
  }
}

/***********************************************************************/
/* compute just one Cart. deriv or the divergence */
/***********************************************************************/

/* compute Cart. deriv in direction dir, put deriv into var with index dui */
int array_cart_1partial(tNode *node, int dir, tArray *u, tArray *du)
{
  tDat *dat = node->dat;
  tArray *dau[3];
  int d, ret;

  if(!dat) return 0;

  /* use du and 2 additional arrays to hold all 3 derivs of u */
  for(d=0; d<dir; d++)
  {
    if(d==dir) dau[d] = du;
    else       dau[d] = alloc_array(node->n);
  }

  /* get all 3 derivs of u */
  ret = array_cart_partials(node, u, dau);

  /* free 2 additional arrays */
  for(d=0; d<dir; d++)
  {
    if(d!=dir) free_array(dau[d]);
  }

  return ret;
}

/* compute Cart. deriv in direction dir, put deriv into var with index dui */
int cart_1partial(tNode *node, int dir, int ui, int dui)
{
  tDat *dat = node->dat;
  tArray *au;
  tArray *dau;

  if(!dat) return 0;

  /* arrays with u and du */
  au  = dat->v[ui];
  dau = dat->v[dui];

  return array_cart_1partial(node, dir, au, dau);
}

/* compute Cart. divergence d_i U^i of vector U^i with index Ux,
   put it into var divUi */
int cart_div_Ui(tNode *node, int Ux, int divUi)
{
  tDat *dat = node->dat;
  tArray *aU;
  tArray *divaU;
  tArray *daU[3];
  tArray *daU0;
  double *divU = Vard(node, divUi);
  double *dU;
  int i;

  if(!dat) return 0;

  /* 4 arrays: divU, and 3 temp. daU */
  divaU = dat->v[divUi];
  daU0   = alloc_array(node->n);
  daU[1] = alloc_array(node->n);
  daU[2] = alloc_array(node->n);

  /* set x-deriv in divaU */
  aU = dat->v[Ux];
  daU[0] = divaU;
  array_cart_partials(node, aU, daU);

  /* add y-deriv */
  aU     = dat->v[Ux+1];
  daU[0] = daU0;
  array_cart_partials(node, aU, daU);
  dU = Arrd(daU[1]);
  forpoints(node,i) divU[i] += dU[i];

  /* add z-deriv */
  aU     = dat->v[Ux+2];
  //daU[0] = daU0;
  array_cart_partials(node, aU, daU);
  dU = Arrd(daU[2]);
  forpoints(node,i) divU[i] += dU[i];

  free(daU0);
  free(daU[1]);
  free(daU[2]);

  return 1;
}
