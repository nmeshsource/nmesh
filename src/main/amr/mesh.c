/* mesh.c */
/* Wolfgang Tichy, 1/2019 */


#include "nmesh.h"
#include "amr.h"

#define PR 0




/* make an empty mesh, into which we an then initialize or into which
   we can e.g. copy the contents of an existing mesh */
tMesh *make_empty_mesh(int pr)
{
  tMesh *mesh;

  /* print info */
  if(pr) prdivider(0);

  mesh = alloc_mesh(0);

  if(pr) { PRFs(":\n"); printmesh(mesh); }

  /* return pointer to newly created mesh */
  return mesh;
}

/* add a patch to the mesh */
tPat *add_patch(tMesh *mesh, double bbox[6], int nroot[3], int nmax,
                int datrank)
{
  tNlist *nlist;
  tPat *pat;
  int p = mesh->npats;
  double dg;
  int i, ni, dir;

  /* check if we have enough space for diff. and other matrices */
  for(dir=0; dir<3; dir++)
    if(nroot[dir] > nmax)
      errorexiti("nmax=%d is too small. Maybe increase par amr_nmax.", nmax);

  /* make room for new patch in mesh and then add an empty patch */
  realloc_patlist_in_mesh(mesh, p + 1);
  pat = alloc_patch(mesh, p, nmax);
  mesh->pat[p] = pat;

  /* set bbox and bbdiag */
  for(i=0; i<6; i++) pat->bbox[i] = bbox[i];
  pat->bbdiag = 0.;
  for(i=0; i<3; i++)
  {
    dg = bbox[2*i+1] - bbox[2*i];
    pat->bbdiag += dg*dg;
  }
  pat->bbdiag = sqrt(pat->bbdiag);

  /* set diff, and other matrices */
  for(dir=0; dir<3; dir++)
  {
    for(ni=1; ni<=nmax; ni++)
    {
      double *Xb = pat->Xb[ni][dir]->d;
      double *Wq = pat->Wq[ni][dir]->d;
      double *WL = pat->WL[ni][dir]->d;
      double *DT = pat->Dt[ni][dir]->d;
      double *AT = pat->At[ni][dir]->d;
      double *ST = pat->St[ni][dir]->d;

      /* since we currently use Legendre Gauss-Lobatto in all 3 dirs
         we do not need e.g. pat->Xb[ni] 3 times */
      if(dir==0) /* set them only for dir0 */
      {
        /* get Legendre Gauss-Lobatto points and integration weights */
        LGL_x_wquad(ni, Xb, Wq);

        /* diff matrix DT for Lagrange interp. poly basis */
        Lagrange_winterp(ni, Xb, WL);
        Lagrange_DT(ni, Xb, WL, DT);

        /* get analysis & synthesis matrix for Legendre basis,
           could be useful for filtering, but not needed for interpolation */
        LGL_AT_ST_matrices(ni, Xb, Wq, AT, ST);
      }
      else /* point dir1 & 2 arrays to the same mem as dir1 arrays */
      {
        free_array(pat->Xb[ni][dir]);
        free_array(pat->Wq[ni][dir]);
        free_array(pat->WL[ni][dir]);
        free_array(pat->Dt[ni][dir]);
        free_array(pat->At[ni][dir]);
        free_array(pat->St[ni][dir]);
        pat->Xb[ni][dir] = pat->Xb[ni][0];
        pat->Wq[ni][dir] = pat->Wq[ni][0];
        pat->WL[ni][dir] = pat->WL[ni][0];
        pat->Dt[ni][dir] = pat->Dt[ni][0];
        pat->At[ni][dir] = pat->At[ni][0];
        pat->St[ni][dir] = pat->St[ni][0];
      }
    }
    /* set Legendre polys as basis since AT and ST are for Legendre basis */
    pat->basis[dir] = basis_normLegendreP;
  }

  /* setup root node */
  pat->rnode = make_root_node(pat, nroot, datrank);
  /* add root node to global mesh->lns list */
  nlist = alloc_nodelist(pat->rnode);
  append_nodelist_to_mesh_lns_myln(mesh, nlist);

  return pat;
}

/* remove all patches from mesh */
void remove_all_patches(tMesh *mesh)
{
  realloc_patlist_in_mesh(mesh, 0);
}



/* select mesh */
int amr_setup_mesh(tMesh *mesh)
{
  int mesh_type = Par("amr_mesh_type");
  int luni = Geti(Par("amr_luni"));
  int refp = Geti(Par("amr_refine_p"));
  int ret;

  if(Getv(mesh_type, "BoxMesh"))
    ret = setup_box_mesh(mesh);
  else if(Getv(mesh_type, "CubedSpheres"))
    ret = setup_CubedSphere_mesh(mesh);
  else if(Getv(mesh_type, "l2_mesh"))
    ret = setup_l2_mesh(mesh);
  else if(Getv(mesh_type, "3patchl2_mesh"))
    ret = setup_3patchl2_mesh(mesh);
  else
    ret = setup_test_mesh(mesh);

  /* load balance root nodes */
  simple_load_balance(mesh);
  //printmesh(mesh);

  /* refine mesh */
  hrefine_mesh_to_level_loadbalance(mesh, luni);
  if(refp >= 0) hrefine_pat(mesh, refp);
/*
hrefine_pat(mesh, 1);
Yo(1);printf("%ld %d\n", mesh->nln, mesh->myln->nm);
 simple_load_balance(mesh);
Yo(1);printf("%ld %d\n", mesh->nln, mesh->myln->nm);

hrefine_pat(mesh, 1);
Yo(2);printf("%ld %d\n", mesh->nln, mesh->myln->nm);
 simple_load_balance(mesh);
Yo(2);printf("%ld %d\n", mesh->nln, mesh->myln->nm);

hrefine_pat(mesh, 0);
Yo(3);printf("%ld %d\n", mesh->nln, mesh->myln->nm);
 simple_load_balance(mesh);
Yo(3);printf("%ld %d\n", mesh->nln, mesh->myln->nm);

hrefine_mesh_to_level(mesh, 3);
Yo(4);printf("%ld %d\n", mesh->nln, mesh->myln->nm);
 simple_load_balance(mesh);
Yo(4);printf("%ld %d\n", mesh->nln, mesh->myln->nm);

hcoarsen_pat(mesh, 0);
Yo(5);printf("%ld %d\n", mesh->nln, mesh->myln->nm);
 simple_load_balance(mesh);
Yo(5);printf("%ld %d\n", mesh->nln, mesh->myln->nm);

hcoarsen_mesh_to_level(mesh, 2);
Yo(6);printf("%ld %d\n", mesh->nln, mesh->myln->nm);
 simple_load_balance(mesh);
Yo(6);printf("%ld %d\n", mesh->nln, mesh->myln->nm);

hcoarsen_pat(mesh, 0);
Yo(7);printf("%ld %d\n", mesh->nln, mesh->myln->nm);
*/

/* // test:
simple_load_balance(mesh);
hcoarsen_mesh_to_level(mesh, 2);
simple_load_balance(mesh);
hcoarsen_pat(mesh, 0);
*/

  /* load balance full mesh */
  simple_load_balance(mesh);
  if(PR) { printmesh(mesh); }

  return ret;
}


/* init neighbor info of root nodes */
int amr_set_bfaces_and_rnode_nfaces_fnb(tMesh *mesh, int pr)
{
  /* setup all bfaces */
  amr_set_all_bfaces(mesh);
  if(pr) printallbfaces(mesh);

  /* now setup root node connections, i.e. setup neighbors of root nodes */
  update_all_rnode_nfaces_fnb(mesh);
  if(pr) printmesh(mesh);
  return 0;
}


/* setup mesh made out of boxes */
int setup_box_mesh(tMesh *mesh)
{
  int mesh_type = Par("amr_mesh_type");
  int npats = Geti(mesh_type);
  char *BoxMesh_xc = Gets(Par("amr_BoxMesh_xc"));
  double d = Getd(Par("amr_BoxMesh_dout"));
  double xc[]   = { 0., 0., 0. };
  double dout[] = { d, d, d };

  PRFs(":\n");

  sscanf(BoxMesh_xc, "%lg %lg %lg", &(xc[0]), &(xc[1]), &(xc[2]));
  //pr3v("xc", xc);printf("\n");

  mesh->dt = Getd(Par("dt"));
  mesh->time = 0.;
  mesh->iteration = 0;

  /* remove all patches to mesh, so we can just ad new pristine ones */
  remove_all_patches(mesh);

  if(Getv(mesh_type, "Line"))
  {
    /* put npats boxes in x-dir */
    add_Nbox_pats_indir(mesh, xc, dout, npats, 0);
  }
  else if(Getv(mesh_type, "Plane"))
  {
    /* arrange npats into one plane */
    int rnpats = pow(npats, 0.5);
    int N[] = { rnpats, rnpats, 1 };
    arrange_box_pats_inBox(mesh, xc, dout, N);
  }
  else
  {
    /* arrange npats into one big box */
    int crnpats = pow(npats, 0.3333333333334);
    int N[] = { crnpats, crnpats, crnpats };
    arrange_box_pats_inBox(mesh, xc, dout, N);
  }

  /* setup all bfaces and root node connections */
  amr_set_bfaces_and_rnode_nfaces_fnb(mesh, 1);

  return 0;
}

/* a mesh with a number of cubed spheres */
int setup_CubedSphere_mesh(tMesh *mesh)
{
  int mesh_type = Par("amr_mesh_type");
  int npats = Geti(mesh_type);
  double rf_surf1 = 0.25;
  double rf_surf2 = 0.25;
  double dc = 0.5;
  double csize = 0.375;
  double ssfac = 4.;
  double obfac = 6.; // 10000 in DNSdata
  double xc[3], rc[3];
  double ABrct[] = { -1.,1., -1.,1. };

  PRFs(":\n");

  mesh->dt = Getd(Par("dt"));
  mesh->time = 0.;
  mesh->iteration = 0;

  /* remove all patches to mesh, so we can just ad new pristine ones */
  remove_all_patches(mesh);

  /* setup cubed spheres */
  switch(npats)
  {
    case 1:
      xc[0] = xc[1] = xc[2] = 0.0;
      add_1_CubedSphere_pat(mesh, 0, outerCubedSphere,0,0,
                            xc, dc, ssfac*dc, ABrct);
      break;
    case 2:
      xc[0] = xc[1] = xc[2] = 0.0;
      add_1_CubedSphere_pat(mesh, 0, outerCubedSphere,0,0,
                            xc, dc, ssfac*dc, ABrct);
      add_1_CubedSphere_pat(mesh, 2, outerCubedSphere,0,0,
                            xc, dc, ssfac*dc, ABrct);
      break;
    case 3:
      if(Getv(mesh_type, "dom1"))
        two_wegdes_touching_1_wedge(mesh, 1.0, 2.0, 3.0);
      else
        two_diff_wegdes_touching_1_wedge(mesh, 1.0, 3.0, 5.0);
      break;
    case 6:
      rc[0] = rc[1] = rc[2] = dc;
      sphere_around_empty_box_at_x0(mesh, rc, ssfac*dc);
      break;
    case 7:
      rc[0] = rc[1] = rc[2] = dc;
      sphere_around_full_box_at_x0(mesh, rc, ssfac*dc);
      break;
    case 12:
      rc[1] = rc[2] = dc; //dc*0.5;
      rc[0] = dc;
      two_spheres_around_empty_box_at_x0(mesh, rc, ssfac*dc, ssfac*obfac*dc);
      break;
    case 13:
      xc[1] = xc[2] = 0.0;
      xc[0] = dc;
      arrange_1pat12CubSph_into_full_cube(mesh, xc,
                                          csize*rf_surf1, rf_surf1, dc);
      break;
    case 26:
      two_full_cubes_touching_at_x0(mesh, dc,
                                    csize*rf_surf1, rf_surf1,
                                    csize*rf_surf2, rf_surf2);
      break;
    case 32:
      sphere_around_two_full_cubes_touching_at_x0(mesh, dc,
                                                  csize*rf_surf1, rf_surf1,
                                                  csize*rf_surf2, rf_surf2,
                                                  ssfac*dc);
      break;
    case 38:
      two_spheres_around_two_full_cubes(mesh, dc,
                                        csize*rf_surf1, rf_surf1,
                                        csize*rf_surf2, rf_surf2,
                                        ssfac*dc, obfac*dc);
      break;
    default:
      errorexiti("amr_mesh_type = %d CubedSpheres  <--not implemented", npats);
  }

coordinates_init(mesh);
outputPatchPlanes_meshvar(mesh, "x", 0,0);
outputPatchPlanes_meshvar(mesh, "y", 0,0);
outputPatchPlanes_meshvar(mesh, "z", 0,0);
//exit(9);

  /* setup all bfaces and root node connections */
  amr_set_bfaces_and_rnode_nfaces_fnb(mesh, 1);

  return 0;
}


/* set up a mesh with 2 levels  */
int setup_l2_mesh(tMesh *mesh)
{
  int amr_n0 = Geti(Par("amr_n0"));
  int amr_n1 = Geti(Par("amr_n1"));
  int amr_n2 = Geti(Par("amr_n2"));
  int n[] = { amr_n0, amr_n1, amr_n2 };
  double bbox[6] = { -4,4, -2,2, -1,1 };
  int n1max = Geti(Par("amr_nmax"));
  tNlist *el, *en;

  PRFs(":\n");

  mesh->dt = Getd(Par("dt"));
  mesh->time = 0.;
  mesh->iteration = 0;

  remove_all_patches(mesh);
  add_patch(mesh, bbox, n, n1max, 0);

  make8children_in_mesh_lns_myln(mesh->lns, n);

  el = mesh->lns;
  for(en = el->next; el; en = el ? el->next : 0)
  {
    if(el->node->l < 2)
    {
      make8children_in_mesh_lns_myln(el, n);
      el = en;
    }
  }

  simple_load_balance(mesh);
  printmesh(mesh);

  return 0;
}

/* set up a mesh with 2 levels  */
int setup_3patchl2_mesh(tMesh *mesh)
{
  int amr_n0 = Geti(Par("amr_n0"));
  int amr_n1 = Geti(Par("amr_n1"));
  int amr_n2 = Geti(Par("amr_n2"));
  int n[] = { amr_n0, amr_n1, amr_n2 };
  double bbox0[6] = { -4,4, -2,2, -1,1 };
  double bbox1[6] = { -4,0,  2,4, -1,1 };
  double bbox2[6] = {  0,4,  2,4, -1,1 };
  int n1max = Geti(Par("amr_nmax"));
  tNlist *el, *en;

  PRFs(":\n");

  mesh->dt = Getd(Par("dt"));
  mesh->time = 0.;
  mesh->iteration = 0;

  remove_all_patches(mesh);
  add_patch(mesh, bbox0, n, n1max, 0);
  add_patch(mesh, bbox1, n, n1max, 0);
  add_patch(mesh, bbox2, n, n1max, 0);

  /* setup all bfaces and root node connections */
  amr_set_bfaces_and_rnode_nfaces_fnb(mesh, 1);

  /* 8 children in patch0 */
  make8children_in_mesh_lns_myln(mesh->lns, n);
  printmesh(mesh);

  /* 8 more in each patch */
  el = mesh->lns;
  for(en = el->next; el; en = el ? el->next : 0)
  {
    if(el->node->l < 2)
    {
      make8children_in_mesh_lns_myln(el, n);
      el = en;
    }
  }

  simple_load_balance(mesh);
  printmesh(mesh);

  return 0;
}



void test_array_thingies(tMesh *mesh)
{
  int n[3];
  int i,j,k;
  double  A[6] = { 1,2,
                   3,4,
                   5,6 };
  /*
  double B0[24] = { 1,2,3,4,5,6,7,8,9,10,11,12,
                    13,14,15,16,17,18,19,20,21,22,23,24 };
  */
  double B0[24], B1[24], B2[24];
  //double AB[36];
  int nB0[] = {2,3,4};
  int nB1[] = {3,2,4};
  int nB2[] = {4,3,2};

  for(i=0; i<12; i++) { B0[2*i] = 2*i+1; B0[2*i+1] = 2*i+2; }

  for(k=0; k<4; k++)
  for(j=0; j<3; j++)
  for(i=0; i<2; i++)
  B2[Ind_n(k,j,i, nB2)] = B1[Ind_n(j,i,k, nB1)] = B0[Ind_n(i,j,k, nB0)];

  n[0]=2; n[1]=2; n[2]=1;
  tArray *Aa = alloc_array(n);
  //Aa->d = A;
  point_array_d_to_data(Aa, A, 1);
  //Aa->n[0]=3; Aa->n[1]=2; Aa->n[1]=1;

  tArray *B0a = alloc_array(nB0);
  //B0a->d = B0;
  point_array_d_to_data(B0a, B0, 1);
  tArray *B1a = alloc_array(nB1);
  //B1a->d = B1;
  point_array_d_to_data(B1a, B1, 1);
  tArray *B2a = alloc_array(nB2);
  //B2a->d = B2;
  point_array_d_to_data(B2a, B2, 1);

  int nC0[] = {2,3,4};
  tArray *C0a = alloc_array(nC0);
  double *C0 = C0a->d;
  //printarray_matrix0(Aa);
  printarray(Aa);
  //printarray_matrix0(B0a);
  printarray(B0a);
  mm_array0(Aa,B0a, C0a);
  //printarray_matrix0(C0a);
  printarray(C0a);

  Yo(1);
  int nC1[] = {3,2,4};
  tArray *C1a = alloc_array(nC1);
  double *C1 = C1a->d;
  //printarray_matrix1(B1a);
  printarray(B1a);
  mm_array1(Aa,B1a, C1a);
  //printarray_matrix0(Ca1);
  set_const_array(C0a, 0.);
  for(k=0; k<4; k++)
  for(j=0; j<2; j++)
  for(i=0; i<3; i++)
  C0[Ind_n(j,i,k, nC0)] = C1[Ind_n(i,j,k, nC1)];
  printarray_matrix0(C0a);
  printarray(C1a);

  Yo(2);
  int nC2[] = {4,3,2};
  tArray *C2a = alloc_array(nC2);
  double *C2 = C2a->d;
  //printarray_matrix2(B2a);
  printarray(B2a);
  mm_array2(Aa,B2a, C2a);
  //printarray_matrix0(Ca1);
  set_const_array(C0a, 0.);
  for(k=0; k<2; k++)
  for(j=0; j<3; j++)
  for(i=0; i<4; i++)
  C0[Ind_n(k,j,i, nC0)] = C2[Ind_n(i,j,k, nC2)];
  printarray_matrix0(C0a);
  printarray(C2a);

  //free_array(Aa);
  //free_array(B0a);
  //free_array(B1a);
  //free_array(B2a);
  free(Aa);
  free(B0a);
  free(B1a);
  free(B2a);
  free_array(C0a);
  free_array(C1a);
  free_array(C2a);
}


/* a function just for testing */
int setup_test_mesh(tMesh *mesh)
{
  double bbox[6] = { -4,4, -2,2, -1,1 };
  int n1max = Geti(Par("amr_nmax"));
  int n[3] = { 5,4,3 };
  tNlist *el, *el2;
  int i;

  PRFs(":\n");

  mesh->dt = Getd(Par("dt"));
  mesh->time = 0.;
  mesh->iteration = 0;

  remove_all_patches(mesh);

//tNode *tnode = alloc_node();
//mesh->pat[0]->rnode = 0;
  //realloc_patlist_in_mesh(mesh, 1);
  add_patch(mesh, bbox, n, n1max, 0);

  enablevar(mesh, Ind("SurfExchange_u"));
  enablevar(mesh, Ind("SurfExchange_v"));
  enablevar(mesh, Ind("X"));
//  tNlist *nlist;
//  tNode *nd;
//  nd = mesh->pat[0]->rnode;
//  nlist = make8_child_nodes(nd, n);
//  replace1_in_mesh_lns_myln(mesh->lns, nlist);
  make8children_in_mesh_lns_myln(mesh->lns, n);

  //printnodelist(nlist);
  printmesh(mesh);

  el = mesh->lns;
  for(i=1; i<=1; i++) el = el->next;
//  nd = el->node;
//  nlist = make8_child_nodes(nd, n);
//  replace1_in_mesh_lns_myln(el, nlist);
  make8children_in_mesh_lns_myln(el, n);

  el = mesh->lns;
  for(i=1; i<=8+2; i++) el = el->next;
//  nd = el->node;
//  nlist = make8_child_nodes(nd, n);
//  replace1_in_mesh_lns_myln(el, nlist);
  el = make8children_in_mesh_lns_myln(el, n);

  //printnodelist(nlist);
  printmesh(mesh);
  printnodelist_and_neighbors(mesh->lns);

  destroy8siblings_in_mesh_lns_myln(el);
  printmesh(mesh);
  printnodelist_and_neighbors(mesh->lns);

  //test_array_thingies(mesh);
  //abort();

  //printarray(mesh->lns->next->node->St[1]);
  printarray_matrix0(mesh->lns->next->node->St[1]);
  printarray_matrix0(mesh->lns->next->node->Dt[1]);

  printarray(mesh->lns->next->node->Xb[1]);
  printarray(mesh->lns->next->node->Wq[1]);

//  el = mesh->lns;
//  printnodelist(el);

  el = alloc_nodelist(mesh->pat[0]->rnode->child[1]);
  printnode_and_neighbors(el->node);

Yo(1);
  printnodelist(el);


Yo(2);
  el2 = ldescendants_along_face(el, 0, &i);
Yo(3);
  printf("i=%d\n",i);
  printnodelist(el2);
  free_nodelist(el);
  free_nodelist(el2);

Yo(4);
//  printnodelist(el2);
//  tNlist *make_mesh_neighbor_list(tNode *node, int face)
  el2 = make_mesh_neighbor_list(mesh->pat[0]->rnode->child[1]->child[2], 0);
  printnodelist(el2);

  free_nodelist(el2);

Yo(5);
prdivider(2);
el = mesh->lns;
for(i=1; i<=8+2; i++) el = el->next;
printnode(el->node);

double *d = Vard(el->node, Ind("SurfExchange_u"));
if(d) d[3] = 3;
printvar_innode(el->node, Ind("SurfExchange_u"));

simple_load_balance(mesh);
printnode(el->node);
printvar_innode(el->node, Ind("SurfExchange_u"));
//printmesh(mesh);
prdivider('^');
//  fflush(stdout);
//  nMPI_barrier();

  return 0;
}
