/* misc.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "misc.h"

#define PR 1


double test_func(double x, double y, double z)
{
  return pow(x-0.1, 4) + pow(y+0.2, 3) + pow(z-0.3, 2);
}


/* try some things */
int misc_test(tMesh *mesh)
{
  tNode *nd;
  int ui = Ind("misc_u");
  int vi = Ind("misc_v");
  int myid, dir, p, k;
  double *Xb[3];
  double X[3], Cb[2];
  double f, interp;
  tArray *coef, *Xp[3], *Cp[2];

  prdivider(0);
  PRF;printf(": Starting misc. tests.\n");
  enablevar(mesh, ui);
  enablevar(mesh, vi);

  formylnodes(mesh, myid)
  {
    int ijk;
    tNode *node = GetMyNode(mesh, myid);
    tArray *ua = GetVarArray(node, ui);

    for(dir=0; dir<3; dir++) Xb[dir] = node->Xb[dir]->d;

    /* set u to func test_func at grid points */
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
  nd = GetMyNode(mesh, 0); /* my first node */
  printarray(nd->Xb[2]);
  printarray(nd->WL[2]);
  f = Lagrange_of_x(1, 0., 3, nd->Xb[2]->d, nd->WL[2]->d);
  printf("Lagrange_of_x at Zb=0: f=%g\n", f);
  printnode(nd);
  printvar_innode(nd, ui);

  /* interpolate in 2 ways */
  /* get coeffs for interp. using basis, i.e. Legendre poly */
  coef = alloc_array(GetVarArray(nd, ui)->n); /* space for coeffs */
  basis_array_analysis3(nd,  GetVarArray(nd, ui), coef);
  basis_array_synthesis3(nd, GetVarArray(nd, vi), coef);
  printvar_innode(nd, vi);
  printarray(coef);

  PRF;printf(": 3d interp. at 2 points with Lagrange and Legendre:\n");
  X[0]=0.9;
  X[1]=0.8;
  X[2]=0.7;
  f = test_func(X[0],X[1],X[2]);
  interp = Lagrange_array_interpolate(nd, GetVarArray(nd, ui), X);
  printf("(%g,%g,%g) -> f=%g interp-f=%g\n", X[0],X[1],X[2], f, interp-f);
  interp = basis_array_interpolate(nd, coef, X);
  printf("(%g,%g,%g) -> f=%g interp-f=%g\n", X[0],X[1],X[2], f, interp-f);

  X[0]=-0.134;
  X[1]=-0.457;
  X[2]=+0.666;
  f = test_func(X[0],X[1],X[2]);
  interp = Lagrange_array_interpolate(nd, GetVarArray(nd, ui), X);
  printf("(%g,%g,%g) -> f=%g interp-f=%g\n", X[0],X[1],X[2], f, interp-f);
  interp = basis_array_interpolate(nd, coef, X);
  printf("(%g,%g,%g) -> f=%g interp-f=%g\n", X[0],X[1],X[2], f, interp-f);

  PRF;printf(": 2d interp. in 3 dir. with Lagrange:\n");
  dir = 0;
  p = 1;
  Cb[0] = X[1];
  Cb[1] = X[2];
  f = test_func(Xb[0][p],X[1],X[2]);
  interp = Lagrange_array_interpolate2d(nd, GetVarArray(nd, ui), dir,p, Cb);
  printf("%d %d: (%g,%g) -> f=%g interp-f=%g\n", dir,p, X[1],X[2], f, interp-f);

  dir = 1;
  p = 1;
  Cb[0] = X[0];
  Cb[1] = X[2];
  f = test_func(X[0],Xb[1][p],X[2]);
  interp = Lagrange_array_interpolate2d(nd, GetVarArray(nd, ui), dir,p, Cb);
  printf("%d %d: (%g,%g) -> f=%g interp-f=%g\n", dir,p, X[0],X[2], f, interp-f);

  dir = 2;
  p = 1;
  Cb[0] = X[0];
  Cb[1] = X[1];
  f = test_func(X[0],X[1],Xb[2][p]);
  interp = Lagrange_array_interpolate2d(nd, GetVarArray(nd, ui), dir,p, Cb);
  printf("%d %d: (%g,%g) -> f=%g interp-f=%g\n", dir,p, X[0],X[1], f, interp-f);

  prdivider(0);
  PRF;printf(": 3d interp. at all points with Lagrange:\n");
  Xp[0] = alloc_array(GetVarArray(nd, ui)->n);
  Xp[1] = alloc_array(GetVarArray(nd, ui)->n);
  Xp[2] = alloc_array(GetVarArray(nd, ui)->n);
  forvari(nd,vi, k) GetVarDpointer(nd,vi)[k] = 666;
  printvar_innode(nd, vi);
  fill_3arrays_with_nodepoints(nd, Xp);
  Lagrange_interpolate_topoints(nd, GetVarArray(nd, ui), Xp,
                                GetVarArray(nd, vi));
  printf("u and v should now be the same:\n");
  printvar_innode(nd, vi);
  printvar_innode(nd, ui);

  prdivider(0);
  PRF;printf(": 2d interp. in plane with Lagrange:\n");
  Cp[0] = alloc_array1d(200);
  Cp[1] = alloc_array1d(200);
  forvari(nd,vi, k) GetVarDpointer(nd,vi)[k] = 666;
  //printvar_innode(nd, vi);
  dir = 0;
  p = 1;
  redim_array(Cp[0], 12,0,0);
  fill_2arrays_with_nodepoints(nd, dir,p, Cp);
  Lagrange_interpolate2d_topoints(nd, GetVarArray(nd, ui), dir,p,
                                  Cp, Xp[0]);
  insert_array_inplane(nd, GetVarArray(nd, vi), dir,p, Xp[0]);
  //printvar_innode(nd, vi);
  
  dir = 1;
  p = 1;
  redim_array(Cp[0], 15,0,0);
  fill_2arrays_with_nodepoints(nd, dir,p, Cp);
  Lagrange_interpolate2d_topoints(nd, GetVarArray(nd, ui), dir,p,
                                  Cp, Xp[0]);
  insert_array_inplane(nd, GetVarArray(nd, vi), dir,p, Xp[0]);
  //printvar_innode(nd, vi);

  dir = 2;
  p = 1;
  redim_array(Cp[0], 20,0,0);
  fill_2arrays_with_nodepoints(nd, dir,p, Cp);
  Lagrange_interpolate2d_topoints(nd, GetVarArray(nd, ui), dir,p,
                                  Cp, Xp[0]);
  insert_array_inplane(nd, GetVarArray(nd, vi), dir,p, Xp[0]);
  printf("u and v should now be the same in plane 1 of all 3 dirs:\n");
  printvar_innode(nd, vi);
  printvar_innode(nd, ui);


  free_array(Cp[0]);
  free_array(Cp[1]);
  free_array(Xp[0]);
  free_array(Xp[1]);
  free_array(Xp[2]);
  free_array(coef);
  return 0;
}

