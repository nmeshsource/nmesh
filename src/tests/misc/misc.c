/* misc.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "misc.h"

#define PR 1


double test_func(double x, double y, double z)
{
  return pow(x,4) + pow(y,3) + pow(z,2);
}


/* try some things */
int misc_test(tMesh *mesh)
{
  tNode *nd;
  int ui = Ind("misc_u");
  int vi = Ind("misc_v");
  int myid, dir;
  double *Xb[3];
  double X[3];
  double f, interp;
  tArray *coef;

  PRF;printf(": Hmmm.\n");
  enablevar(mesh, ui);
  enablevar(mesh, vi);

  formylnodes(mesh, myid)
  {
    int ijk;
    tNode *node = GetMyNode(mesh, myid);
    tArray *ua = GetVarArray(node, ui);

    for(dir=0; dir<3; dir++) Xb[dir] = node->Xb[dir]->d;

    /* set particular pattern in u */
    forarray(ua, ijk)
    {
      int k = kOfInd_n(ijk, ua->n);
      int j = jOfInd_n_k(ijk, ua->n, k);
      int i = iOfInd_n_jk(ijk, ua->n, j,k);
      double x = Xb[0][i];
      double y = Xb[1][j];
      double z = Xb[2][k];

      ua->d[ijk] = test_func(x,y,z);
    }
  }

  /* print var in one node */
  prdivider('I');
  nd = GetMyNode(mesh, 0); /* my first node */
  printarray(nd->Xb[2]);
  printarray(nd->WL[2]);

  prdivider('I');
  printnode(nd);
  printvar_innode(nd, ui);

  /* interpolate */
  f = Lagrange_of_x(1, 0., 3, nd->Xb[2]->d, nd->WL[2]->d);
  X[0]=0.9;
  X[1]=0.8;
  X[2]=0.7;
  printf("f=%g\n", f);
  f = test_func(X[0],X[1],X[2]);
  interp = Lagrange_array_interpolate(nd, GetVarArray(nd, ui), X);
  printf("(%g,%g,%g) -> f=%g interp-f=%g\n", X[0],X[1],X[2], f ,interp-f);

  coef = alloc_array(GetVarArray(nd, ui)->n);
  basis_array_analysis3(nd,  GetVarArray(nd, ui), coef);
  basis_array_synthesis3(nd, GetVarArray(nd, vi), coef);
  //printvar_innode(nd, vi);
  //printarray(coef);
  interp = basis_array_interpolate(nd, coef, X);
  printf("(%g,%g,%g) -> f=%g interp-f=%g\n", X[0],X[1],X[2], f ,interp-f);

  X[0]=-0.134;
  X[1]=-0.457;
  X[2]=+0.666;
  f = test_func(X[0],X[1],X[2]);
  interp = Lagrange_array_interpolate(nd, GetVarArray(nd, ui), X);
  printf("(%g,%g,%g) -> f=%g interp-f=%g\n", X[0],X[1],X[2], f ,interp-f);

  coef = alloc_array(GetVarArray(nd, ui)->n);
  basis_array_analysis3(nd,  GetVarArray(nd, ui), coef);
  basis_array_synthesis3(nd, GetVarArray(nd, vi), coef);
  //printvar_innode(nd, vi);
  //printarray(coef);
  interp = basis_array_interpolate(nd, coef, X);
  printf("(%g,%g,%g) -> f=%g interp-f=%g\n", X[0],X[1],X[2], f ,interp-f);

  free_array(coef);
  return 0;
}

