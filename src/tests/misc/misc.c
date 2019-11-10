/* misc.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "misc.h"

#define PR 1




double test_func(double x, double y, double z)
{
  return pow(x-0.1, 4) + pow(y+0.2, 3) + pow(z-0.3, 2);
}


double test_func2(double x[3], int lmax[3], int n[3])
{
  int d, l;
  double sum[] = { 0.,0.,0. };

  for(d=0; d<3; d++)
    for(l=0; l<=lmax[d]; l++)
    {
      sum[d] += (d+1. + 0.001*(d+1)*l)*basis_normLegendreP(l, x[d], n[d]);
    }

  return sum[0] * sum[1] * sum[2];
  return sum[0] + sum[1] + sum[2];
}


/* try some things */
int misc_test(tMesh *mesh)
{
  int ui = Ind("misc_u");
  int vi = Ind("misc_v");

  prdivider(0);
  PRF;printf(": Starting misc. tests.\n");
  enablevar(mesh, ui);
  enablevar(mesh, vi);

  formylnodes(mesh)
  {
    int ijk;
    tNode *node = MyLnode;
    tArray *ua = VarA(node, ui);
    int dir;
    double *Xb[3];

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
  test_indc(mesh);
  test_node_av(mesh);
  test_ajsurf(mesh);
  test_0doutput(mesh);
  test_filter(mesh);

  return 0;
}


/* interpolate to a few points */
int test_point_interpolation(tMesh *mesh)
{
  tNode *nd;
  int n2;
  int ui = Ind("misc_u");
  int vi = Ind("misc_v");
  int dir, p, k;
  double *Xb[] = { NULL, NULL, NULL };
  double X[3], Cb[2];
  double f, interp;
  tArray *coef, *Xp[3], *Cp[2];

  prdivider(0);
  PRF;printf(": Starting misc. tests.\n");
  nd = Lnode_myid(mesh, mesh->myln->nm-1); /* my last node??? */
  for(dir=0; dir<3; dir++) Xb[dir] = nd->Xb[dir]->d;

  /* print var in one node */
  nd = Lnode_myid(mesh, 0); /* my first node */
  printarray(nd->Xb[2]);
  printarray(nd->WL[2]);
  n2 = nd->n[2];
  f = Lagrange_of_x(n2/2, 0., n2, nd->Xb[2]->d, nd->WL[2]->d);
  printf("Lagrange_of_x at Zb=0: f=%g\n", f);
  printnode(nd);
  printvar_innode(nd, ui);

  /* interpolate in 2 ways */
  /* get coeffs for interp. using basis, i.e. Legendre poly */
  coef = alloc_array(VarA(nd, ui)->n); /* space for coeffs */
  basis_array_analysis3(nd,  VarA(nd, ui), coef);
  basis_array_synthesis3(nd, VarA(nd, vi), coef);
  printvar_innode(nd, vi);
  printarray(coef);

  PRF;printf(": 3d interp. at 2 points with Lagrange and Legendre:\n");
  X[0]=0.9;
  X[1]=0.8;
  X[2]=0.7;
  f = test_func(X[0],X[1],X[2]);
  interp = Lagrange_array_interpolate(nd, VarA(nd, ui), X);
  printf("(%g,%g,%g) -> f=%g interp-f=%g\n", X[0],X[1],X[2], f, interp-f);
  interp = basis_array_interpolate(nd, coef, X);
  printf("(%g,%g,%g) -> f=%g interp-f=%g\n", X[0],X[1],X[2], f, interp-f);

  X[0]=-0.134;
  X[1]=-0.457;
  X[2]=+0.666;
  f = test_func(X[0],X[1],X[2]);
  interp = Lagrange_array_interpolate(nd, VarA(nd, ui), X);
  printf("(%g,%g,%g) -> f=%g interp-f=%g\n", X[0],X[1],X[2], f, interp-f);
  interp = basis_array_interpolate(nd, coef, X);
  printf("(%g,%g,%g) -> f=%g interp-f=%g\n", X[0],X[1],X[2], f, interp-f);

  PRF;printf(": 2d interp. in 3 dir. with Lagrange:\n");
  dir = 0;
  p = (n2>1);
  Cb[0] = X[1];
  Cb[1] = X[2];
  f = test_func(Xb[0][p],X[1],X[2]);
  interp = Lagrange_array_interpolate2d(nd, VarA(nd, ui), dir,p, Cb);
  printf("%d %d: (%g,%g) -> f=%g interp-f=%g\n", dir,p, X[1],X[2], f, interp-f);

  dir = 1;
  p = (n2>1);
  Cb[0] = X[0];
  Cb[1] = X[2];
  f = test_func(X[0],Xb[1][p],X[2]);
  interp = Lagrange_array_interpolate2d(nd, VarA(nd, ui), dir,p, Cb);
  printf("%d %d: (%g,%g) -> f=%g interp-f=%g\n", dir,p, X[0],X[2], f, interp-f);

  dir = 2;
  p = (n2>1);
  Cb[0] = X[0];
  Cb[1] = X[1];
  f = test_func(X[0],X[1],Xb[2][p]);
  interp = Lagrange_array_interpolate2d(nd, VarA(nd, ui), dir,p, Cb);
  printf("%d %d: (%g,%g) -> f=%g interp-f=%g\n", dir,p, X[0],X[1], f, interp-f);

  prdivider(0);
  PRF;printf(": 3d interp. at all points with Lagrange:\n");
  Xp[0] = alloc_array(VarA(nd, ui)->n);
  Xp[1] = alloc_array(VarA(nd, ui)->n);
  Xp[2] = alloc_array(VarA(nd, ui)->n);
  forvari(nd,vi, k) Vard(nd,vi)[k] = 666;
  printvar_innode(nd, vi);
  fill_3arrays_with_nodepoints(nd, Xp);
  Lagrange_interpolate_topoints(nd, VarA(nd, ui), Xp,
                                VarA(nd, vi));
  printf("u and v should now be the same:\n");
  printvar_innode(nd, vi);
  printvar_innode(nd, ui);

  prdivider(0);
  PRF;printf(": 2d interp. in plane with Lagrange:\n");
  Cp[0] = alloc_array1d(100*n2*n2);
  Cp[1] = alloc_array1d(100*n2*n2);
  forvari(nd,vi, k) Vard(nd,vi)[k] = 666;
  //printvar_innode(nd, vi);
  dir = 0;
  p = (n2>1);
  redim_array(Cp[0], 11*p+1,0,0);
  fill_2arrays_with_nodepoints(nd, dir, Cp);
  Lagrange_interpolate2d_topoints(nd, VarA(nd, ui), dir,p,
                                  Cp, Xp[0]);
  insert_array_inplane(VarA(nd, vi), dir,p, Xp[0]);
  //printvar_innode(nd, vi);

  dir = 1;
  p = (n2>1);
  redim_array(Cp[0], 14*p+1,0,0);
  fill_2arrays_with_nodepoints(nd, dir, Cp);
  Lagrange_interpolate2d_topoints(nd, VarA(nd, ui), dir,p,
                                  Cp, Xp[0]);
  insert_array_inplane(VarA(nd, vi), dir,p, Xp[0]);
  //printvar_innode(nd, vi);

  dir = 2;
  p = (n2>1);
  redim_array(Cp[0], 19*p+1,0,0);
  fill_2arrays_with_nodepoints(nd, dir, Cp);
  Lagrange_interpolate2d_topoints(nd, VarA(nd, ui), dir,p,
                                  Cp, Xp[0]);
  insert_array_inplane(VarA(nd, vi), dir,p, Xp[0]);
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
  tArray *ua = VarA(node, ui);
  tArray *va = VarA(node, vi);
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
  update_mesh_myln_node_nid(mesh);
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
  int iX = Ind("X");
  int ix = Ind("x");
  //double *Xbd[3];
  double sum, Sum;
  int test_func_lamAB = Getv(Par("misc_ajsurf_v_init"), "test_func(lam,A,B)");

  prdivider(0);
  PRF;printf(": Hmmm.\n");
  enablevar(mesh, vi);

  /* above we messed with all kinds of things,
     so make sure all coords are set again */
  coordinates_init(mesh);
//nd = Lnode_myid(mesh, 1);
//printvar_innode(nd, ix);

  formylnodes(mesh)
  {
    int ijk;
    //int dir;
    tNode *node = MyLnode;
    tArray *va = VarA(node, vi);

    //for(dir=0; dir<3; dir++) Xbd[dir] = node->Xb[dir]->d;

    /* set v to func test_func at grid points */
    forarray(va, ijk)
    {
      //int k = kOfInd_n(ijk, va->n);
      //int j = jOfInd_n_k(ijk, va->n, k);
      //int i = iOfInd_n_jk(ijk, va->n, j,k);
      //double Xb[] = { Xbd[0][i], Xbd[1][j], Xbd[2][k] };
      //double X[3];
      double X[] = { Vard(node, iX)[ijk],
                     Vard(node, iX+1)[ijk], Vard(node, iX+2)[ijk] };
      double x[] = { Vard(node, ix)[ijk],
                     Vard(node, ix+1)[ijk], Vard(node, ix+2)[ijk] };
      tPat *pat = node->pat;
      int p = pat->p;
      double dlam = (p != 0);
      //double dlam = (p > 5);
      double lam = X[0] + dlam;
      //double A = dom!=3 ? X[1] : 2.*(1.-X[1]) + 1.;
      //double A = dom/2 ? X[1] : -X[1];
      double A = X[1];
      double B = X[2];
      //int dom = pat->CI->dom;
      //double B = dom/2 ? X[2] : -X[2];

      //XYZ_of_XbYbZb(node, Xb, X);
      va->d[ijk] = test_func(x[0],x[1],x[2]) + 0.000 * node->nid * node->nid;
      //va->d[ijk] = test_func(X[0],X[1],X[2]) + 0.000 * node->nid * node->nid;
      if(test_func_lamAB) va->d[ijk] = test_func(lam,A,B);
      //va->d[ijk] = 1000*node->nid +100*X[0] +10*X[1] +X[2];
    }
  }

  /* print var */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    printnode(node);
    printvar_innode(node, vi);
  }

  /* print var in one node again */
  nd = Lnode_myid(mesh, 0); /* my first node */
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
//  nd = Lnode_myid(mesh, 7);
//  printnode(nd);
  //printvar_innode(nd, ix);
  //printvar_innode(nd, Ind("oC0_1"));
  //printvar_innode(nd, Ind("oC1_1"));
//  printvar_innode(nd, vi);

//  nd = Lnode_myid(mesh, 78);
//  printnode(nd);
  //printvar_innode(nd, ix);
  //printvar_innode(nd, Ind("oC0_1"));
  //printvar_innode(nd, Ind("oC1_1"));
//  printvar_innode(nd, vi);

  /* print var in all nodes */
  prdivider(0);
  PRF;printf(": ajsurfdiff on all nodes:\n");
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    printnode(node);
    printvar_ajsurfdiff(node, vi);
  }

  /* free redundant nbsurf stuff */
  free_all_myln_nbsurf_only(mesh);

  /* print var in all nodes again */
  prdivider(0);
  PRF;printf(": ajsurfdiff on all nodes after freeing nbsurf:\n");
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    printnode(node);
    printvar_ajsurfdiff(node, vi);
  }

  PRF;printf(": L2 norm of ajsurfdiff:\n");
  sum = 0.;
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    int f;
    double norm_n_f, norm_n_f2;
    char s[100];

    for(f=0; f<6; f++)
    {
      tSurface *sf = node->dat->s[f][vi];
      if(!sf) continue;

      norm_n_f = Lp_norm_array_reldiff(sf->ajsurf, sf->mysurf, 2);
      printf("  nid%ld %s f%d: %g\n", node->nid,
             nodename(node,s,99), f, norm_n_f);
      norm_n_f2 = pow(norm_n_f, 2);
      //#pragma omp atomic
      GEN_Pragma(omp atomic)
      sum += norm_n_f2;
    }
  }
  printf("on this proc: total %.15g\n", sqrt(sum));
  Sum = sum;
  nMPI_Allreduce(&sum, &Sum, 1, nMPI_DOUBLE, nMPI_SUM);
  printf("on all procs: total %.15g\n", sqrt(Sum));

  /* after we have printed them, we no longer need the surfaces */
  free_all_myln_surfaces(mesh);
  return 0;
}

/* exchange some indicators for testing */
int test_indc(tMesh *mesh)
{
  tNode *nd;
  int ui = Ind("misc_u");
  int vi = Ind("misc_v");
  tVarList *vl;

  /* varlist with ui and vi */
  vl = vlalloc(mesh);
  vlpush(vl, ui);
  vlpush(vl, vi);

  prdivider(0);
  PRF;printf(": indc.\n");
  enablevar(mesh, vi);

  /* above we messed with all kinds of things,
     so make sure all coords are set again */
  coordinates_init(mesh);

  /* print var */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    printnode(node);
    printvar_innode(node, vi);
  }

  /* print var in one node again */
  nd = Lnode_myid(mesh, 0); /* my first node */
  printnode(nd);
  printvar_innode(nd, vi);

  /* exchange indc */
  prdivider(0);
  PRF;printf(": request_all_myln_indc_exchange_for_vl\n");
  init_all_myln_myindc_for_vl(mesh, vl, 3);
  /* set indc to something ... */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    node->dat->ic[vi]->myindc->d[0] = node->nid; 
    node->dat->ic[vi]->myindc->d[1] = -node->nid; 
  }
  request_all_myln_indc_exchange_for_vl(mesh, vl);

  /* Here we can do work. MPI is now busy sending buffers */

  /* now get the indicators and wait for MPI buffers if necessary */
  get_all_myln_indc_for_vl(mesh, vl);

  /* nbindc should be set now */
  PRF;printf(": get_all_myln_indc_for_vl has set nbindc\n");

  /* print var in all nodes */
  prdivider(0);
  PRF;printf(": print the var indc on all nodes:\n");
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    printnode(node);
    printvar_indc(node, vi);
  }

  /* after we have printed them, we no longer need the indicators */
  free_all_myln_indc_for_vl(mesh, vl);

  vlfree(vl);
  return 0;
}

/* compute some node avrages */
int test_node_av(tMesh *mesh)
{
  int ui = Ind("misc_u");
  int vi = Ind("misc_v");

  prdivider(0);
  PRF;printf(": node average\n");
  enablevar(mesh, vi);

  /* above we messed with all kinds of things,
     so make sure all coords are set again */
  coordinates_init(mesh);

  /* print var and its average */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;

    prdivider(0);
    printnode(node);
    //printvar_innode(node, vi);
    basis_var_analysis3(node, vi, ui);
    printf("c_{000}*sqrt(2)^3 of vi / 8  = %g\n",
           Vard(node, ui)[0] * pow(sqrt(2.),3.)/8.);
    printf("var_GLquadrature3 of vi / 8  = %g\n",
           var_GLquadrature3(node, vi)/8.);
    printf("var_nodeaverage of vi        = %g\n", var_nodeaverage(node, vi));
  }

  return 0;
}

/* 0d output */
int test_0doutput(tMesh *mesh)
{
  int ui = Ind("misc_u");

  prdivider(0);
  PRF;printf(":\n");
  enablevar(mesh, ui);

  /* above we messed with all kinds of things,
     so make sure all coords are set again */
  coordinates_init(mesh);

  /* print var and its average */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    double *u = Vard(node, ui);
    int ijk;

    forpoints(node, ijk) u[ijk] = 1.0; //node_location(node);
  }

  return 0;
}

/* filter */
int test_filter(tMesh *mesh)
{
  //tNode *nd;
  int ui = Ind("misc_u");
  int vi = Ind("misc_v");
  int iX = Ind("X");
  int ix = Ind("x");
  double alp[3], s[3];
  int dir;

  prdivider(0);
  PRF;printf(": ...\n");
  enablevar(mesh, ui);
  enablevar(mesh, vi);

  /* above we messed with all kinds of things,
     so make sure all coords are set again */
  coordinates_init(mesh);

  formylnodes(mesh)
  {
    int ijk;
    //int dir;
    tNode *node = MyLnode;
    double *ud = Vard(node, ui);

    /* set v to func test_func at grid points */
    forpoints(node, ijk)
    {
      int lmax[] = { 7,7,7 };// { 3,3,3 };
      //double Xb[] = { Xbd[0][i], Xbd[1][j], Xbd[2][k] };
      double X[] = { Vard(node, iX)[ijk],
                     Vard(node, iX+1)[ijk], Vard(node, iX+2)[ijk] };
      double x[] = { Vard(node, ix)[ijk],
                     Vard(node, ix+1)[ijk], Vard(node, ix+2)[ijk] };
      //tPat *pat = node->pat;
      //int p = pat->p;
      //double dlam = (p != 0);
      //double dlam = (p > 5);
      //double lam = X[0] + dlam;
      //double A = dom!=3 ? X[1] : 2.*(1.-X[1]) + 1.;
      //double A = dom/2 ? X[1] : -X[1];
      //double A = X[1];
      //double B = X[2];
      //int dom = pat->CI->dom;
      //double B = dom/2 ? X[2] : -X[2];

      //XYZ_of_XbYbZb(node, Xb, X);
      ud[ijk] = test_func2(x, lmax, node->n) + 0.000 * node->nid * node->nid;
      ud[ijk] = test_func2(X, lmax, node->n);
      //ud[ijk] = 1000*node->nid +100*X[0] +10*X[1] +X[2];
    }
  }

  PRF;printf(": before filter\n");

  /* print coeffs of ui in vi */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    printnode(node);
    basis_var_analysis3(node, ui, vi);
    printvar_innode(node, vi);
  }

  /* filter var ui */
  //Setd(Par("basis_expfilter_JacobianPower"), 0.);
  for(dir=0; dir<3; dir++)
  {
    alp[dir] = 36.;
    s[dir] = 32.;
  }
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    expfilter_var(node, ui, alp, s);
  }

  PRF;printf(": after expfilter_var\n");

  /* print coeffs of ui in vi again */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    printnode(node);
    basis_var_analysis3(node, ui, vi);
    printvar_innode(node, vi);
  }

  return 0;
}
