/* misc.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "misc.h"

#define PR 1


int test_point_interpolation(tMesh *mesh);
int test_point_finders(tMesh *mesh);
int test_parent_child_interpolation(tMesh *mesh);
int test_ajsurf(tMesh *mesh);


double test_func(double x, double y, double z)
{
  return pow(x-0.1, 4) + pow(y+0.2, 3) + pow(z-0.3, 2);
}


/* try some things */
int misc_test(tMesh *mesh)
{
  int ui = Ind("misc_u");
  int vi = Ind("misc_v");
  int myid, dir;
  double *Xb[3];

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

  test_point_interpolation(mesh);
  test_point_finders(mesh);
  test_parent_child_interpolation(mesh);
  test_ajsurf(mesh);

  return 0;
}


/* interpolate to a few points */
int test_point_interpolation(tMesh *mesh)
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
  formylnodes(mesh, myid)
  {
    tNode *node = GetMyNode(mesh, myid);
    for(dir=0; dir<3; dir++) Xb[dir] = node->Xb[dir]->d;
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
  fill_2arrays_with_nodepoints(nd, dir, Cp);
  Lagrange_interpolate2d_topoints(nd, GetVarArray(nd, ui), dir,p,
                                  Cp, Xp[0]);
  insert_array_inplane(GetVarArray(nd, vi), dir,p, Xp[0]);
  //printvar_innode(nd, vi);

  dir = 1;
  p = 1;
  redim_array(Cp[0], 15,0,0);
  fill_2arrays_with_nodepoints(nd, dir, Cp);
  Lagrange_interpolate2d_topoints(nd, GetVarArray(nd, ui), dir,p,
                                  Cp, Xp[0]);
  insert_array_inplane(GetVarArray(nd, vi), dir,p, Xp[0]);
  //printvar_innode(nd, vi);

  dir = 2;
  p = 1;
  redim_array(Cp[0], 20,0,0);
  fill_2arrays_with_nodepoints(nd, dir, Cp);
  Lagrange_interpolate2d_topoints(nd, GetVarArray(nd, ui), dir,p,
                                  Cp, Xp[0]);
  insert_array_inplane(GetVarArray(nd, vi), dir,p, Xp[0]);
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



/* print misc_u - test_func */
int print_u_minus_f(tNode *node)
{
  tMesh *mesh = node->pat->mesh;
  int ui = Ind("misc_u");
  int vi = Ind("misc_v");
  tArray *ua = GetVarArray(node, ui);
  tArray *va = GetVarArray(node, vi);
  int dir;
  double *Xb[3];

  prdivider(0);
  PRFs(": v = u - f:\n");

  if(ua!=NULL && va!=NULL)
  {
    int ijk;

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

      va->d[ijk] = ua->d[ijk] - test_func(x,y,z);
    }
  }

  printvar_innode(node, vi);

  return 0;
}

/* test interpolation inside make8children_in_mesh_lns_myln and
   destroy8siblings_in_mesh_lns_myln */
int test_parent_child_interpolation(tMesh *mesh)
{
  tNode *nd;
  int ui = Ind("misc_u");
  int vi = Ind("misc_v");
  int nn[] = { 3,5,5 };
  tNlist *el = mesh->lns;
  tNlist *l2 = NULL;
  tDat *d0;

  prdivider(0);
  PRF;printf(": Starting misc. tests.\n");
  enablevar(mesh, ui);
  enablevar(mesh, vi);

  print_u_minus_f(mesh->lns->node);

  /* get 1st node */
  el = mesh->lns;
  nd = el->node;
  d0 = nd->dat;
  printnode(el->node);
  printvar_innode(nd, ui);
  if(d0) printf("1 nd %p %p %d\n", nd, d0, d0->nv);

  make8children_in_mesh_lns_myln(el, nn);
  //printmesh(mesh);
  el = mesh->lns;
  l2 = NULL;
  l2 = addnode_to_nodelist_after(l2, el->next->next->node);
  addnode_to_nodelist_after(l2, el->next->next->next->node);
  move_nodelist_to_rank(l2, (nMPI_size()>1));
  free_nodelist(l2);
  l2 = NULL;
  //printmesh(mesh);

  printf("2 nd %p %p\n", nd, nd->dat);
  //if(d0) printf("2b nd %p %p %d\n", nd, d0, d0->nv);
  printf("2c: test_func=%g\n", test_func(-4,-2,-1));
  printf("2c: test_func=%g\n", test_func(-1,-1,-1));
  el = mesh->lns;
  printnode(el->node);
  printvar_innode(nd->child[0], ui);

  destroy8siblings_in_mesh_lns_myln(el);
  printf("3 nd %p %p\n", nd, nd->dat);
  //if(d0) printf("4 nd %p %p %d\n", nd, d0, d0->nv);
  el = mesh->lns;
  printnode(el->node);
  printvar_innode(nd, ui);

  print_u_minus_f(el->node);

  return 0;
}


/* try some things */
int test_point_finders(tMesh *mesh)
{
  tNode *nd;
  int dir;
  int nn[] = { 3,5,5 };
  double XX0[] = { -4.-1e-12, -1.1, -0.01 };
  double XX1[] = { -4.-0.1,   -1.1, -0.01 };
  double XX2[] = { -2,        -1.1, -0.01 };
  tArray *Xc[3], *Xd[3];

  prdivider(0);
  PRF;printf(": Starting misc. tests.\n");

  /* 1st node */
  nd = mesh->lns->node;

  Xc[0] = alloc_array(nn);
  Xc[1] = alloc_array(nn);
  Xc[2] = alloc_array(nn);
  Xd[0] = alloc_array(nn);
  Xd[1] = alloc_array(nn);
  Xd[2] = alloc_array(nn);
  for(dir=0; dir<3; dir++)
  {
    Xc[dir]->d[0] = XX0[dir];
    Xc[dir]->d[1] = XX1[dir];
    Xc[dir]->d[2] = XX2[dir];
  }
  redim_array(Xc[0], 3,1,1);
  redim_array(Xc[1], 3,1,1);
  redim_array(Xc[2], 3,1,1);
  redim_array(Xd[0], 3,1,1);
  redim_array(Xd[1], 3,1,1);
  redim_array(Xd[2], 3,1,1);
  printarray(Xc[0]);
  printarray(Xc[1]);
  printarray(Xc[2]);

  printarray(Xd[0]);
  printarray(Xd[1]);
  printarray(Xd[2]);

  printnode(nd);
  //array_get_XYZ_in_node(nd, Xc, Xd);
  array_find_XYZ_in_node(nd, Xc, Xd[0]);

  //array_XbYbZb_of_XYZ(nd, Xd, Xc);
  printarray(Xc[0]);
  printarray(Xc[1]);
  printarray(Xc[2]);
  printf("in=%d: XX1 %.15e %.15e %.15e\n",
          XYZ_is_in_node(nd, XX1), XX1[0],XX1[1],XX1[2]);
  printarray_int(Xd[0]);
  printarray(Xd[1]);
  printarray(Xd[2]);
  //array_get_XYZ_in_node(nd, Xd, Xc);
  //printarray(Xc[0]);
  //printarray(Xc[1]);
  //printarray(Xc[2]);

  free_array(Xc[0]);
  free_array(Xc[1]);
  free_array(Xc[2]);
  free_array(Xd[0]);
  free_array(Xd[1]);
  free_array(Xd[2]);
  return 0;
}


/* exchange some surfaces for testing */
int test_ajsurf(tMesh *mesh)
{
  tNode *nd;
  int vi = Ind("misc_v");
  double *Xbd[3];
  int myid;

  PRF;printf(": Hmmm.\n");
  enablevar(mesh, vi);

  formylnodes(mesh, myid)
  {
    int ijk, dir;
    tNode *node = GetMyNode(mesh, myid);
    tArray *va = GetVarArray(node, vi);

    for(dir=0; dir<3; dir++) Xbd[dir] = node->Xb[dir]->d;

    /* set v to func test_func at grid points */
    forarray(va, ijk)
    {
      int k = kOfInd_n(ijk, va->n);
      int j = jOfInd_n_k(ijk, va->n, k);
      int i = iOfInd_n_jk(ijk, va->n, j,k);
      double Xb[] = { Xbd[0][i], Xbd[1][j], Xbd[2][k] };
      double X[3];

      XYZ_of_XbYbZb(node, Xb, X);

      va->d[ijk] = test_func(X[0],X[1],X[2]) + 0.000 * node->nid * node->nid;
    }
  }

  /* print var */
  formylnodes(mesh, myid)
  {
    tNode *node = GetMyNode(mesh, myid);
    printnode(node);
    printvar_innode(node, vi);
  }

  /* print var in one node again */
  nd = GetMyNode(mesh, 0); /* my first node */
  printnode(nd);
  printvar_innode(nd, vi);

  /* exchange surfaces */
  prdivider(0);
  PRF;printf(": request_all_myln_surfaces_exchange\n");
  init_all_myln_surfaces(mesh);
  set_all_myln_mysurf(mesh);
  request_all_myln_surfaces_exchange(mesh);

  /* Here we can do work. MPI is now busy sending buffers */

  /* now get the surfaces and wait for buffers if necessary */
  get_all_myln_surfaces(mesh);

  /* get_all_myln_surfaces sets ajsurf via interpolation */
  PRF;printf(": get_all_myln_surfaces has set ajsurf via interpolation\n");

  /* print var in one node yet again with surfaces */
  nd = GetMyNode(mesh, 0); /* my first node */
  printnode(nd);
  printvar_innode(nd, vi);

  /* print var in all nodes */
  prdivider(0);
  PRF;printf(": ajsurfdiff on all nodes:\n");
  formylnodes(mesh, myid)
  {
    tNode *node = GetMyNode(mesh, myid);
    printnode(node);
    printvar_ajsurfdiff(node, vi);
  }

  /* free redundant nbsurf stuff */
  free_all_myln_nbsurf_only(mesh);

  /* print var in all nodes again */
  prdivider(0);
  PRF;printf(": ajsurfdiff on all nodes after freeing nbsurf:\n");
  formylnodes(mesh, myid)
  {
    tNode *node = GetMyNode(mesh, myid);
    printnode(node);
    printvar_ajsurfdiff(node, vi);
  }

  /* after we have printed them, we no longer need the surfaces */
  free_all_myln_surfaces(mesh);
  return 0;
}
