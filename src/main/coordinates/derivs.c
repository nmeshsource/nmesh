/* derivs.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "coordinates.h"


/* global vars */
extern tcoordinates coordinates[1];


/***********************************************************************/
/* compute one deriv of a scalar with respect to one X^i coord */
/***********************************************************************/

/* compute coord deriv of u in array au, put du/dX^dir into array dau */
void coordinate_array_deriv1(tNode *node, int dir, tArray *au, tArray *dau,
                             tDerivOpt *opt)
{
  double *du = Arrd(dau);
  double dXbdX[3];
  int ind;

  /* take derivs with respect to Xb: du/dXb */
  basis_array_deriv1(node, dir, au, dau, opt);

  /* get dXb/dX */
  dXbYbZb_dXYZ(node, dXbdX);

  /* scale: du/dX = dXb/dX du/dXb */
  forpoints(node,ind)
    du[ind] *= dXbdX[dir];
}

/* compute coord deriv, put du/dX^dir into var with index dui */
int coordinate_deriv1(tNode *node, int dir, int ui, int dui, tDerivOpt *opt)
{
  tDat *dat = node->dat;
  tArray *au;
  tArray *dau;

  if(!dat) return 0;

  au  = dat->v[ui];
  dau = dat->v[dui];

  coordinate_array_deriv1(node, dir, au, dau, opt);
  return 1;
}

/* Transform Cartesian vector components u^i to patch coord components U^I.
   We can do it in place if aU = au. */
int coordinate_dXdx_times_vector_array(tNode *node, tArray *au[3],
                                       tArray *aU[3])
{
  tPat *pat = node->pat;
  //tMesh *mesh = pat->mesh;
  tDat *dat = node->dat;
  double *u[] = { Arrd(au[0]), Arrd(au[1]), Arrd(au[2]) };
  double *U[] = { Arrd(aU[0]), Arrd(aU[1]), Arrd(aU[2]) };
  int ind, m,i;

  if(!dat) return 0;

  /* do we need to init. coords? */
  if(!(dat->coords_set)) coordinates_init_node(node);

  /* transform from Cartesian to patch coords */
  if(pat->dXYZ_dxyz)
  {
    int idXd = coordinates->idXdx;
    double *dXdx[3][3]
              = { {Vard(node,idXd),   Vard(node,idXd+1), Vard(node,idXd+2)},
                  {Vard(node,idXd+3), Vard(node,idXd+4), Vard(node,idXd+5)},
                  {Vard(node,idXd+6), Vard(node,idXd+7), Vard(node,idXd+8)} };
    /* compute Cartesian derivs at all points */
    forpoints(node,ind)
    {
      double V[3];

      /* Transform derivs from Cartesian to patch coords */
      for(m=0; m<3; m++)
      {
        V[m] = 0.;
        for(i=0; i<3; i++) V[m] += dXdx[m][i][ind] * u[i][ind];
      }
      /* copy V into U */
      for(m=0; m<3; m++) U[m][ind] = V[m];
    }
  }
  /* X^I is Cartesian, so copy components of au into aU if needed */
  else if(aU != au)
  {
    for(m=0; m<3; m++)
      copy_array_data(au[m], aU[m]);
  }
  return 1;
}

/* Transform from coord derivs daU to Cartesian derivs dau.
   We can do it in place if dau = daU. */
int coordinate_dXdx_times_1form_array(tNode *node, tArray *daU[3],
                                      tArray *dau[3])
{
  tPat *pat = node->pat;
  //tMesh *mesh = pat->mesh;
  tDat *dat = node->dat;
  double *dU[] = { Arrd(daU[0]), Arrd(daU[1]), Arrd(daU[2]) };
  double *du[] = { Arrd(dau[0]), Arrd(dau[1]), Arrd(dau[2]) };
  int ind, m,i;

  if(!dat) return 0;

  /* do we need to init. coords? */
  if(!(dat->coords_set)) coordinates_init_node(node);

  /* transform to Cartesian coords */
  if(pat->dXYZ_dxyz)
  {
    int idXd = coordinates->idXdx;
    double *dXdx[3][3]
              = { {Vard(node,idXd),   Vard(node,idXd+1), Vard(node,idXd+2)},
                  {Vard(node,idXd+3), Vard(node,idXd+4), Vard(node,idXd+5)},
                  {Vard(node,idXd+6), Vard(node,idXd+7), Vard(node,idXd+8)} };
    /* compute Cartesian derivs at all points */
    forpoints(node,ind)
    {
      double dw[3];

      /* Transform derivs to Cartesian coords */
      for(m=0; m<3; m++)
      {
        dw[m] = 0.;
        for(i=0; i<3; i++) dw[m] += dXdx[i][m][ind] * dU[i][ind];
      }
      /* copy dw into du */
      for(m=0; m<3; m++) du[m][ind] = dw[m];
    }
  }
  /* daU is Cartesian, so copy components of daU into dau if needed */
  else if(dau != daU)
  {
    for(m=0; m<3; m++)
      copy_array_data(daU[m], dau[m]);
  }
  return 1;
}

/***********************************************************************/
/* compute all 3 Cartesian 1st derivs of a scalar (using arrays) */
/***********************************************************************/

/* compute Cart. derivs of u in array au, put du/dx^m into arrays dau[0..2] */
int array_cart_partials(tNode *node, tArray *au, tArray *dau[3],
                        tDerivOpt *opt)
{
  tDat *dat = node->dat;
  double *du[] = { Arrd(dau[0]), Arrd(dau[1]), Arrd(dau[2]) };
  double dXbdX[3];
  int ind, m;

  if(!dat) return 0;

  /* take derivs with respect to Xb: du/dXb */
  basis_array_derivs(node, au, dau, opt);

  /* get dXb/dX */
  dXbYbZb_dXYZ(node, dXbdX);

  /* scale: du/dX = dXb/dX du/dXb */
  forpoints(node,ind)
    for(m=0; m<3; m++) du[m][ind] *= dXbdX[m];

  /* this is slightly slower than the 3 steps above: */
  //for(m=0; m<3; m++)
  //  coordinate_array_deriv1(node,m, au, dau[m], opt);

  /* transform to Cartesian coords */
  coordinate_dXdx_times_1form_array(node, dau, dau);
  return 1;
}


/* compute Cart. derivs, put du/dx^m into vars with index dui[0..2] */
int cart_partials(tNode *node, int ui, int dui[3], tDerivOpt *opt)
{
  tDat *dat = node->dat;
  tArray *au;
  tArray *dau[3];

  if(!dat) return 0;

  au     = dat->v[ui];
  dau[0] = dat->v[dui[0]];
  dau[1] = dat->v[dui[1]];
  dau[2] = dat->v[dui[2]];

  return array_cart_partials(node, au, dau, opt);
}

/***********************************************************************/
/* 2 variants of cart_partials to get the 1st derivs of a scalar */
/***********************************************************************/

/* compute first derivs U_{,i} of a scalar U in a node */
void cart_3partials(tNode *node, int U, int dUx, int dUy, int dUz,
                    tDerivOpt *opt)
{
  int dU[] = { dUx, dUy, dUz };
  cart_partials(node, U, dU, opt);
}

/* compute first derivs U_{,i} of a scalar U in a node */
void cart_partials_diScalar(tNode *node, int U, int dUx, tDerivOpt *opt)
{
  int dU[] = { dUx, dUx+1, dUx+2 };
  cart_partials(node, U, dU, opt);
}

/***********************************************************************/
/* 1st derivs of some special 3-vectors and 3-tensors */
/***********************************************************************/

/* compute first derivs U_{i,j} of a vector U_{i} in a node */
void cart_partials_dUi_dj(tNode *node, int Ux, int dUxx)
{
  /* compute partial derivs of all components in node */
  cart_partials_diScalar(node, Ux,   dUxx,   NULL);
  cart_partials_diScalar(node, Ux+1, dUxx+3, NULL);
  cart_partials_diScalar(node, Ux+2, dUxx+6, NULL);
}

/* compute first derivs S_{ij,k} of a symmetric tensor S_{ij} in a node */
void cart_partials_dSij_dk(tNode *node, int Sxx, int dSxxx)
{
  int n;

  /* compute partial derivs of all components in node */
  for(n=0; n<6; n++)
    cart_partials_diScalar(node, Sxx + n, dSxxx + 3*n, NULL);
}

/* compute first derivs U_{ij,k} of a general tensor U_{ij} in a node */
void cart_partials_dUij_dk(tNode *node, int Uxx, int dUxxx)
{
  int n;

  /* compute partial derivs of all components in node */
  for(n=0; n<9; n++)
    cart_partials_diScalar(node, Uxx + n, dUxxx + 3*n, NULL);
}

/* compute first derivs T_{ijk,l} of a tensor with T_{ijk} = T_{ikj}
   or with T_{ijk} = T_{jik} (in both cases we have 18 comps in T) */
void cart_partials_dTijk_dl(tNode *node, int Txxx, int dTxxxx)
{
  int n;

  /* compute partial derivs of all components in node */
  for(n=0; n<18; n++)
    cart_partials_diScalar(node, Txxx + n, dTxxxx + 3*n, NULL);
}


/* compute first derivs d_i U_{j} of a vector U_{j} in a node */
void cart_partials_diUj(tNode *node, int Ux, int dUxx)
{
  int n;

  /* compute partial derivs of all components in node */
  for(n=0; n<3; n++)
    cart_3partials(node, Ux+n, dUxx+n, dUxx+3+n, dUxx+6+n, NULL);
}

/* compute first derivs d_i S_{jk} of a symmetric tensor S_{jk} in a node */
void cart_partials_diSjk(tNode *node, int Sxx, int dSxxx)
{
  int n;

  /* compute partial derivs of all components in node */
  for(n=0; n<6; n++)
    cart_3partials(node, Sxx+n, dSxxx+n, dSxxx+6+n, dSxxx+12+n, NULL);
}

/* compute first derivs d_i U_{jk} of a general tensor U_{jk} in a node */
void cart_partials_diUjk(tNode *node, int Uxx, int dUxxx)
{
  int n;

  /* compute partial derivs of all components in node */
  for(n=0; n<9; n++)
    cart_3partials(node, Uxx+n, dUxxx+n, dUxxx+9+n, dUxxx+18+n, NULL);
}

/* compute first derivs d_i T_{jkl} of a tensor with T_{jkl} = T_{jlk}
   or with T_{ijk} = T_{jik} (in both cases we have 18 comps in T) */
void cart_partials_diTjkl(tNode *node, int Txxx, int dTxxxx)
{
  int n;

  /* compute partial derivs of all components in node */
  for(n=0; n<18; n++)
    cart_3partials(node, Txxx+n, dTxxxx+n, dTxxxx+18+n, dTxxxx+36+n, NULL);
}


/***********************************************************************/
/* 1st derivs of some special 4-vectors and 4-tensors */
/***********************************************************************/

/* compute first derivs U_{a,i} of a vector U_{a} in a node */
void cart_partials_dUa_di(tNode *node, int Ut, int dUtx)
{
  /* compute partial derivs of all components in node */
  cart_partials_diScalar(node, Ut,   dUtx,   NULL);
  cart_partials_diScalar(node, Ut+1, dUtx+3, NULL);
  cart_partials_diScalar(node, Ut+2, dUtx+6, NULL);
  cart_partials_diScalar(node, Ut+3, dUtx+9, NULL);
}

/* compute first derivs S_{ab,k} of a symmetric tensor S_{ab} in a node */
void cart_partials_dSab_di(tNode *node, int Stt, int dSttx)
{
  int n;

  /* compute partial derivs of all components in node */
  for(n=0; n<10; n++)
    cart_partials_diScalar(node, Stt + n, dSttx + 3*n, NULL);
}

/* compute first derivs U_{ab,k} of a general tensor U_{ab} in a node */
void cart_partials_dUab_di(tNode *node, int Utt, int dUttx)
{
  int n;

  /* compute partial derivs of all components in node */
  for(n=0; n<16; n++)
    cart_partials_diScalar(node, Utt + n, dUttx + 3*n, NULL);
}

/* compute first derivs d_i U_{a} of a vector U_{a} in a node */
void cart_partials_diUa(tNode *node, int Ut, int dUxt)
{
  int n;

  /* compute partial derivs of all components in node */
  for(n=0; n<4; n++)
    cart_3partials(node, Ut+n, dUxt+n, dUxt+4+n, dUxt+8+n, NULL);
}

/* compute first derivs d_i S_{ab} of a symmetric tensor S_{ab} in a node */
void cart_partials_diSab(tNode *node, int Stt, int dSxtt)
{
  int n;

  /* compute partial derivs of all components in node */
  for(n=0; n<10; n++)
    cart_3partials(node, Stt+n, dSxtt+n, dSxtt+10+n, dSxtt+20+n, NULL);
}

/* compute first derivs d_i U_{ab} of a general tensor U_{ab} in a node */
void cart_partials_diUab(tNode *node, int Utt, int dUxtt)
{
  int n;

  /* compute partial derivs of all components in node */
  for(n=0; n<16; n++)
    cart_3partials(node, Utt+n, dUxtt+n, dUxtt+16+n, dUxtt+32+n, NULL);
}

/* compute first derivs d_i T_{jab} of a tensor T_{iab} = T_{iba} */
void cart_partials_diTjab(tNode *node, int Txtt, int dTxxtt)
{
  int n;

  /* compute partial derivs of all components in node */
  for(n=0; n<30; n++)
    cart_3partials(node, Txtt+n, dTxxtt+n, dTxxtt+30+n, dTxxtt+60+n, NULL);
}


/***********************************************************************/
/* 1st derivs of general tensors */
/***********************************************************************/

/* write number of components of var T0 and dT0 into nT and ndT, and also
   check if dT0 has correct number of components for spatial derivs of T0 */
void cart_partials_SetAndCheck_nT_ndT(tNode *node, int T0, int dT0,
                                      int *nT, int *ndT)
{
  tMesh *mesh = node->pat->mesh;
  *nT = MeshVarNComponents(mesh, T0);
  *ndT = MeshVarNComponents(mesh, dT0);

  if( (*ndT) != 3*(*nT) )
  {
    char *T = MeshVarName(mesh, T0);
    char *dT = MeshVarName(mesh, dT0);
    char *Tindices = MeshVarTensorIndices(mesh, T0);
    char *dTindices = MeshVarTensorIndices(mesh, dT0);
    printf("%s %s (T0=%d) has nT=%d components\n",
           T, Tindices, T0, *nT);
    printf("%s %s (dT0=%d) has ndT=%d components\n",
           dT, dTindices, dT0, *ndT);
    errorexit("To store all 3 spatial derivs we need ndT = 3*nT.");
  }
}

/* Compute first derivs T_{... ,k} of an arbitrary tensor T. The derivative
   index is understood to be the last index of the resultant dT, i.e.:
   dT_{...k} = T_{... ,k} */
void cart_partials_dTensor_di(tNode *node, int T0, int dT0, tDerivOpt *opt)
{
  int nT, ndT, n;

  /* get and check number of components in T0 and dT0 */
  cart_partials_SetAndCheck_nT_ndT(node, T0,dT0, &nT, &ndT);

  /* compute partial derivs of all components of Tensor */
  for(n=0; n<nT; n++)
    cart_partials_diScalar(node, T0 + n, dT0 + 3*n, opt);
}

/* Compute first derivs d_i T of an arbitrary tensor T. The derivative
   index is understood to be the first index of the resultant dT, i.e.:
   dT_{i...} = d_i T_{...} */
void cart_partials_diTensor(tNode *node, int T0, int dT0, tDerivOpt *opt)
{
  int nT, ndT, n;

  /* get and check number of components in T0 and dT0 */
  cart_partials_SetAndCheck_nT_ndT(node, T0,dT0, &nT, &ndT);

  /* compute partial derivs of all components of Tensor */
  for(n=0; n<nT; n++)
    cart_3partials(node, T0 + n, dT0 + n, dT0 + nT + n, dT0 + 2*nT + n, opt);
}


/* write number of components of var T0 and dT0 into nT and ndT, and also
   check if dT0 has correct number of components for one coord deriv of T0 */
void coordinate_deriv1_SetAndCheck_nT_ndT(tNode *node, int dir,
                                          int T0, int dT0, int *nT, int *ndT)
{
  tMesh *mesh = node->pat->mesh;
  int ok;

  *nT = MeshVarNComponents(mesh, T0);
  *ndT = MeshVarNComponents(mesh, dT0);

  /* dT either has to have same number of comps or 3 times as many */
  ok = ( (*ndT) == (*nT) ) || ( (*ndT) == 3*(*nT) );

  if(!ok)
  {
    char *T = MeshVarName(mesh, T0);
    char *dT = MeshVarName(mesh, dT0);
    char *Tindices = MeshVarTensorIndices(mesh, T0);
    char *dTindices = MeshVarTensorIndices(mesh, dT0);
    printf("%s %s (T0=%d) has nT=%d components\n",
           T, Tindices, T0, *nT);
    printf("%s %s (dT0=%d) has ndT=%d components\n",
           dT, dTindices, dT0, *ndT);
    errorexit("To store one coordinate deriv we need ndT=nT, or ndT=3*nT.");
  }
}

/* Compute first coord deriv dT/dX of an arbitrary tensor T. Here, either
   the resultant dT has no derivative index, or the derivative index is
   understood to be the last index of the resultant dT, i.e.:
   dT_{...k} = T_{... ,k} */
void coordinate_deriv1_dTensor_dX(tNode *node, int dir, int T0, int dT0,
                                  tDerivOpt *opt)
{
  int nT, ndT, off, n;

  /* get and check number of components in T0 and dT0 */
  coordinate_deriv1_SetAndCheck_nT_ndT(node, dir, T0,dT0, &nT, &ndT);

  /* set offset par in dT where we put derivs */
  if(ndT==3*nT) off = 1;
  else          off = 0;

  /* compute partial derivs of all components of Tensor */
  for(n=0; n<nT; n++)
    coordinate_deriv1(node, dir, T0 + n, dT0 + (1+2*off)*n + off*dir, opt);
}

/* Compute first coord deriv d_X T of an arbitrary tensor T. Here, either
   the resultant dT has no derivative index, or the derivative index is
   understood to be the first index of the resultant dT, i.e.:
   dT_{i...} = d_Xi T_{...} */
void coordinate_deriv1_dXTensor(tNode *node, int dir, int T0, int dT0,
                                tDerivOpt *opt)
{
  int nT, ndT, pos, n;

  /* get and check number of components in T0 and dT0 */
  coordinate_deriv1_SetAndCheck_nT_ndT(node, dir, T0,dT0, &nT, &ndT);

  /* set position in dT where we put derivs */
  if(ndT==3*nT) pos = dir*nT;
  else          pos = 0;

  /* compute partial derivs of all components of Tensor */
  for(n=0; n<nT; n++)
    coordinate_deriv1(node, dir, T0 + n, dT0 + pos + n, opt);
}


/***********************************************************************/
/* 2nd derivs of vectors and tensors */
/***********************************************************************/

/* compute 1st and 2nd order Cart. derivs of scalar U */
void cart_partials_dU_di_dU_dij(tNode *node, int U, int dUx, int ddUxx)
{
  /* 1st derivs */
  cart_partials_diScalar(node, U, dUx, NULL);

  /* 2nd derivs */
  cart_3partials(node, dUx,   ddUxx,  ddUxx+1,ddUxx+2, NULL);
  cart_3partials(node, dUx+1, ddUxx+1,ddUxx+3,ddUxx+4, NULL);
  cart_3partials(node, dUx+2, ddUxx+2,ddUxx+4,ddUxx+5, NULL);
}

/* compute 1st and 2nd derivs U_{i,jk} of a vector U_{i} in a node */
void cart_partials_dUi_dj_dUi_djk(tNode *node, int Ux, int dUxx, int ddUxxx)
{
  int n;

  /* 1st derivs */
  cart_partials_dUi_dj(node, Ux, dUxx);

  /* 2nd derivs */
  for(n=0; n<3; n++) /* n=0: dUxj, n=1: dUyj, n=2: dUzj */
  {
    cart_3partials(node, dUxx+3*n ,  ddUxxx+6*n ,  ddUxxx+6*n+1, ddUxxx+6*n+2, NULL);
    cart_3partials(node, dUxx+3*n+1, ddUxxx+6*n+1, ddUxxx+6*n+3, ddUxxx+6*n+4, NULL);
    cart_3partials(node, dUxx+3*n+2, ddUxxx+6*n+2, ddUxxx+6*n+4, ddUxxx+6*n+5, NULL);
  }
}

/* compute 1st and 2nd derivs S_{ij,kl} of a symmetric tensor S_{ij} */
void cart_partials_dSij_dk_dSij_dkl(tNode *node, int Sxx,
                                    int dSxxx, int ddSxxxx)
{
  int n;

  /* 1st derivs */
  cart_partials_dSij_dk(node, Sxx, dSxxx);

  /* 2nd derivs */
  for(n=0; n<6; n++) /* n=0: dSxxj, n=2: dSxzj, n=3: dSxyj */
  {
    cart_3partials(node, dSxxx+3*n ,  ddSxxxx+6*n ,  ddSxxxx+6*n+1, ddSxxxx+6*n+2, NULL);
    cart_3partials(node, dSxxx+3*n+1, ddSxxxx+6*n+1, ddSxxxx+6*n+3, ddSxxxx+6*n+4, NULL);
    cart_3partials(node, dSxxx+3*n+2, ddSxxxx+6*n+2, ddSxxxx+6*n+4, ddSxxxx+6*n+5, NULL);
  }
}

/***********************************************************************/
/* 2nd derivs of general tensors from pre-existing 1st derivs */
/***********************************************************************/

/* write number of components of var dT0 and ddT0 into ndT and nddT, and also
   check if ddT0 has correct number of components for 2nd derivs of T0 */
void cart_partials_SetAndCheck_nT_ndT_nddT(tNode *node, int dT0, int ddT0,
                                           int *nT, int *ndT, int *nddT)
{
  tMesh *mesh = node->pat->mesh;
  *ndT = MeshVarNComponents(mesh, dT0);
  *nddT = MeshVarNComponents(mesh, ddT0);

  /* number of components of tensor T0 that resulted in tensor deriv dT0 */
  *nT = *ndT/3;

  if( (*nddT) != 6*(*nT) ) /* ddT is symm. in deriv indices */
  {
    char *dT = MeshVarName(mesh, dT0);
    char *ddT = MeshVarName(mesh, ddT0);
    char *dTindices = MeshVarTensorIndices(mesh, dT0);
    char *ddTindices = MeshVarTensorIndices(mesh, ddT0);
    printf("%s %s (dT0=%d) has ndT=%d components\n",
           dT, dTindices, dT0, *ndT);
    printf("Thus the tensor without derivs (T0) must have nT=%d components\n",
           *nT);
    printf("%s %s (ddT0=%d) has nddT=%d components\n",
           ddT, ddTindices, ddT0, *nddT);
    errorexit("To store all 6 2nd deriv components we need nddT = 6*nT.");
  }
}

/* Compute 2nd derivs T_{... ,ij} from 1st derivs T_{... ,i} of a general
   tensor T_{...}. The resulting T_{... ,ij} needs to be defined as
   symmetric in the last 2 indices. */
void cart_partials_ddTensor_dij(tNode *node,
                                int dT0, int ddT0, tDerivOpt *opt)
{
  int nT, ndT, nddT;
  int n;

  /* get and check number of components of dT and ddT and also set nT */
  cart_partials_SetAndCheck_nT_ndT_nddT(node, dT0,ddT0, &nT, &ndT, &nddT);

  /* 2nd derivs */
  for(n=0; n<nT; n++) /* ddT is symm in last 2 indices ==> fac 6 below */
  {
    cart_3partials(node, dT0+3*n ,  ddT0+6*n ,  ddT0+6*n+1, ddT0+6*n+2, opt);
    cart_3partials(node, dT0+3*n+1, ddT0+6*n+1, ddT0+6*n+3, ddT0+6*n+4, opt);
    cart_3partials(node, dT0+3*n+2, ddT0+6*n+2, ddT0+6*n+4, ddT0+6*n+5, opt);
  }
}

/* Compute 2nd derivs d_i d_j T_{...} from 1st derivs d_i T_{...} of a
   general tensor T_{...}.  The resulting d_i d_j T_{...} needs to be
   defined as symmetric in the first 2 indices. */
void cart_partials_didjTensor(tNode *node,
                              int dT0, int ddT0, tDerivOpt *opt)
{
  int nT, ndT, nddT;
  int n;

  /* get and check number of components of dT and ddT and also set nT */
  cart_partials_SetAndCheck_nT_ndT_nddT(node, dT0,ddT0, &nT, &ndT, &nddT);

  //printf("nT=%d ndT=%d nddT=%d\n", nT, ndT, nddT);
  //errorexit("2nd deriv part of this function is untested");

  /* 2nd derivs */
  for(n=0; n<nT; n++) //each 2nd deriv ddT has nT comps => steps of nT below
  {
    cart_3partials(node,dT0+     n, ddT0      +n,ddT0+  nT+n,ddT0+2*nT+n, opt);
    cart_3partials(node,dT0+  nT+n, ddT0+   nT+n,ddT0+3*nT+n,ddT0+4*nT+n, opt);
    cart_3partials(node,dT0+2*nT+n, ddT0+ 2*nT+n,ddT0+4*nT+n,ddT0+5*nT+n, opt);
  }
}

/***********************************************************************/
/* 1st and 2nd derivs of general tensors */
/***********************************************************************/

/* write number of components of var T0 and ddT0 into nT and nddT, and also
   check if ddT0 has correct number of components for 2nd derivs of T0 */
void cart_partials_SetAndCheck_nT_nddT(tNode *node, int T0, int ddT0,
                                       int *nT, int *nddT)
{
  tMesh *mesh = node->pat->mesh;
  *nT = MeshVarNComponents(mesh, T0);
  *nddT = MeshVarNComponents(mesh, ddT0);

  if( (*nddT) != 6*(*nT) ) /* ddT is symm. in deriv indices */
  {
    char *T = MeshVarName(mesh, T0);
    char *ddT = MeshVarName(mesh, ddT0);
    char *Tindices = MeshVarTensorIndices(mesh, T0);
    char *ddTindices = MeshVarTensorIndices(mesh, ddT0);
    printf("%s %s (T0=%d) has nT=%d components\n",
           T, Tindices, T0, *nT);
    printf("%s %s (ddT0=%d) has nddT=%d components\n",
           ddT, ddTindices, ddT0, *nddT);
    errorexit("To store all 6 2nd deriv components we need nddT = 6*nT.");
  }
}

/* Compute 1st and 2nd derivs T_{... ,i} and T_{... ,ij} of a general
   tensor T_{...}.  The resulting T_{... ,ij} needs to be defined as
   symmetric in the last 2 indices. */
void cart_partials_dTensor_di_ddTensor_dij(tNode *node, int T0,
                                           int dT0, int ddT0, tDerivOpt *opt)
{
  int nT, nddT;

  /* 1st derivs */
  cart_partials_dTensor_di(node, T0, dT0, opt);

  /* get and check number of components of T and ddT */
  cart_partials_SetAndCheck_nT_nddT(node, T0,ddT0, &nT, &nddT);

  /* 2nd derivs */
  cart_partials_ddTensor_dij(node, dT0, ddT0, opt);
}

/* Compute 1st and 2nd derivs d_i T_{...} and d_i d_j T_{...} of a general
   tensor T_{...}.  The resulting d_i d_j T_{...} needs to be defined as
   symmetric in the first 2 indices. */
void cart_partials_diTensor_didjTensor(tNode *node, int T0,
                                       int dT0, int ddT0, tDerivOpt *opt)
{
  int nT, nddT;

  /* 1st derivs */
  cart_partials_diTensor(node, T0, dT0, opt);

  /* get and check number of components of T and ddT */
  cart_partials_SetAndCheck_nT_nddT(node, T0,ddT0, &nT, &nddT);

  /* 2nd derivs */
  cart_partials_didjTensor(node, dT0, ddT0, opt);
}


/***********************************************************************/
/* compute just one Cart. deriv or the divergence */
/***********************************************************************/

/* compute Cart. deriv in direction dir, put deriv into var with index dui */
int array_cart_1partial(tNode *node, int dir, tArray *u, tArray *du,
                        tDerivOpt *opt)
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
  ret = array_cart_partials(node, u, dau, opt);

  /* free 2 additional arrays */
  for(d=0; d<dir; d++)
  {
    if(d!=dir) free_array(dau[d]);
  }

  return ret;
}

/* compute Cart. deriv in direction dir, put deriv into var with index dui */
int cart_1partial(tNode *node, int dir, int ui, int dui, tDerivOpt *opt)
{
  tDat *dat = node->dat;
  tArray *au;
  tArray *dau;

  if(!dat) return 0;

  /* arrays with u and du */
  au  = dat->v[ui];
  dau = dat->v[dui];

  return array_cart_1partial(node, dir, au, dau, opt);
}

/* compute Cart. divergence d_i U^i of vector U^i with index Ux,
   put it into var divUi */
/* NOTE: maybe one could get d_i U^i from [d_A (sqrt(f) U^A)]/sqrt(f)
   where f is the flat metric in X^A coords. This might need less memory! */
int cart_div3Vector(tNode *node, int Ux, int divUi, tDerivOpt *opt)
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
  //daU0   = alloc_array(node->n);
  //daU[1] = alloc_array(node->n);
  //daU[2] = alloc_array(node->n);
  daU0   = VarA(node, coordinates->itmp1);
  daU[1] = VarA(node, coordinates->itmp1+1);
  daU[2] = VarA(node, coordinates->itmp1+2);

  /* set x-deriv in divaU */
  aU = dat->v[Ux];
  daU[0] = divaU;
  array_cart_partials(node, aU, daU, opt);

  /* add y-deriv */
  aU     = dat->v[Ux+1];
  daU[0] = daU0;
  array_cart_partials(node, aU, daU, opt);
  dU = Arrd(daU[1]);
  forpoints(node,i) divU[i] += dU[i];

  /* add z-deriv */
  aU     = dat->v[Ux+2];
  //daU[0] = daU0;
  array_cart_partials(node, aU, daU, opt);
  dU = Arrd(daU[2]);
  forpoints(node,i) divU[i] += dU[i];

  //free_array(daU0);
  //free_array(daU[1]);
  //free_array(daU[2]);

  return 1;
}
