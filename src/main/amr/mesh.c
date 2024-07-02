/* mesh.c */
/* Wolfgang Tichy, 1/2019 */


#include "nmesh.h"
#include "amr.h"

#define PR 0


/* use gridpoints from basis/gridpoints.c */
extern tGridPoints gridpoints[1];
extern tAMR amr[1];


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


/* choose grid points x_i and weights w_i for patch p */
int choose_patch_points(tMesh *mesh, int p)
{
  int ret;
  intList *pl = alloc_intList();

  /* push all ints from amr_uniform_p into pl */
  str_to_intList(Gets(Par("amr_uniform_p")), " ", pl);

  /* now check if patch p is mentioned in amr_uniform_p */
  if(in_intList(pl, p))
    ret = P_UNIFORM;
  else
    ret = P_LGL; /* default is Legendre Gauss-Lobatto */

  free_intList(pl);
  return ret;
}

/* set use_fv on some nodes */
int amr_set_use_fv_flag(tMesh *mesh)
{
  intList *pl = alloc_intList();

  /* push all ints from scalarwave1_fv_p into pl */
  str_to_intList(Gets(Par("amr_fv_p")), " ", pl);

  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    int p = node->pat->p;
    int use_fv = in_intList(pl, p);

    if(use_fv)
    {
      /* put node into FV mode */
      node->dat->info->use_fv = 1;
    }
  }
  free_intList(pl);
  return 0;
}

/* switch nodes with uniform spacing into fv mode */
int amr_use_fv_if_P_UNIFORM(tMesh *mesh)
{
  int ptUNI[] = { P_UNIFORM, P_UNIFORM, P_UNIFORM };
  int ptLGL[] = { P_LGL, P_LGL, P_LGL };

  /* switch on fv on all uniform nodes */
  refine_set_use_fv_if_pt_typ(mesh, ptUNI, 1);

  /* switch off fv on all LGL nodes */
  refine_set_use_fv_if_pt_typ(mesh, ptLGL, 0);

  return 0;
}


/* add a patch to the mesh */
tPat *add_patch_without_rnode(tMesh *mesh, double bbox[6])
{
  tPat *pat;
  int p = mesh->npats;
  double dg;
  int i;

  //PRFs(":\n");

  /* make room for new patch in mesh and then add an empty patch */
  realloc_patlist_in_mesh(mesh, p + 1);
  pat = alloc_patch(mesh, p);
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

  return pat;
}

/* add a patch to the mesh, and put a root node/elm in it */
tPat *add_patch(tMesh *mesh, double bbox[6],
                int *pt_typ_root, int nroot[3], int datrank)
{
  int nmax = gridpoints->nmax;
  tPat *pat;
  int p = mesh->npats;
  int dir;
  int pt_typ[3];

  //PRFs(":\n");

  /* check if we have enough space for diff. and other matrices */
  for(dir=0; dir<3; dir++)
    if(nroot[dir] > nmax)
      errorexiti("nmax=%d is too small. Maybe increase par amr_nmax.", nmax);

  /* make patch first, add root elm later */
  pat = add_patch_without_rnode(mesh, bbox);

  /* get points (e.g. Legendre Gauss-Lobatto) and integration weights */
  if(pt_typ_root)
  {
    for(dir=0; dir<3; dir++) pt_typ[dir] = pt_typ_root[dir];
  }
  else
  {
    int typ = choose_patch_points(mesh, p);
    for(dir=0; dir<3; dir++) pt_typ[dir] = typ;
  }

  /* setup root node element */
  make_and_add_root_elm(pat, nroot, pt_typ, datrank);

  //if(1) print_matrices_innode(pat->rnode);
  //exit(88);

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
  int hrefp = Par("amr_hrefine_p");
  int prefp = Par("amr_prefine_p");
  int sph_l = Geti(Par("amr_hrefine_sphere_levels"));
  double sph_r = Geti(Par("amr_hrefine_sphere_radius"));
  double x0[3] = {0.};
  double BoxMesh_rc[] = {0.,0.,0.};
  double CubedSphere_r = 0.;
  int ret = 0;

  /* remove all patches from mesh, so we can just add new pristine ones */
  PRFs(": calling remove_all_patches\n");
  remove_all_patches(mesh);

  /* call different funcs based on par amr_mesh_type */
  PRFs(": setup patches according to pars\n");
  if(Getv(mesh_type, "BoxMesh"))
    ret = setup_box_mesh(mesh, BoxMesh_rc);
  if(Getv(mesh_type, "CubedSpheres"))
    ret = setup_CubedSphere_mesh(mesh, BoxMesh_rc, &CubedSphere_r);
  if(Getv(mesh_type, "Shell"))
    ret = setup_Shell_mesh(mesh, CubedSphere_r);

  if( (Getv(mesh_type, "BoxMesh")) || (Getv(mesh_type, "CubedSpheres")) ||
      (Getv(mesh_type, "Shell")) )
  {
    /* setup all bfaces and root node connections */
    PRFs(": calling amr_set_bfaces_and_rnode_nbinfo_fnb\n");
    amr_set_bfaces_and_rnode_nbinfo_fnb(mesh, 0,1);
  }
  else  /* test meshes */
  {
    PRFs(": make some test mesh\n");
    if(Getv(mesh_type, "test_mesh"))
      ret = test_mesh(mesh);
    else if(Getv(mesh_type, "l2_mesh"))
      errorexit("add back:  ret = setup_l2_mesh(mesh);");
    else if(Getv(mesh_type, "3patchl2_mesh"))
      ret = setup_3patchl2_mesh(mesh);
    else
      ret = setup_elm_mesh1(mesh); // test case
    /* NOTE: all test meshes call amr_set_bfaces_and_rnode_nbinfo_fnb
             by themselves */
  }

  /* apply par amr_innerboundary_faces */
  if(GetLen(Par("amr_innerboundary_faces")) > 0)
  {
    char *flist = Gets(Par("amr_innerboundary_faces"));
    char *fl, *str, *sav;

    PRFs(": marking some faces as inner boundary\n");
    fl = strdup(flist);
    for(str=strtok_r(fl, " ", &sav); str!=NULL;
        str=strtok_r(NULL, " ", &sav))
    {
      int f = atoi(str);
      if(f>=0 && f<6)
        mark_bfacesonface_without_op_as_boundary(mesh, f, INNERBOUND);
    }
    free(fl);
  }

  /* after amr_set_bfaces_and_rnode_nbinfo_fnb and
     mark_bfacesonface_without_op_as_boundary it is good to see bfaces */
  PRFs(":\n");
  printallbfaces(mesh);

  ///* print info about all patches we have now */
  //forpatches(mesh,pind)
  //{
  //  tPat *pat = mesh->pat[pind];
  //  printpatch(pat);
  //  printCI(pat);
  //}

  /* load balance root nodes */
  PRFs(": calling simple_load_balance\n");
  simple_load_balance(mesh);
  //printmesh(mesh);

  /* h-refine mesh uniformly */
  PRFs(": calling hrefine_mesh_to_level_loadbalance\n");
  hrefine_mesh_to_level_loadbalance(mesh, luni);

  /* now h-refine the patches listed in amr_hrefine_p */
  if(GetLen(hrefp) > 0)
  {
    char *plist = Gets(hrefp);
    char *pl, *str, *sav;

    PRFs(": calling hrefine_pat for some patches\n");
    pl = strdup(plist);
    for(str=strtok_r(pl, " ", &sav); str!=NULL;
        str=strtok_r(NULL, " ", &sav))
    {
      int p = atoi(str);
      if(p>=0) hrefine_pat(mesh, p);
    }
    free(pl);
  }

  /* h-refine further in nested sphere regions */
  PRFs(": calling hrefine_sphere_loadbalance\n");
  hrefine_sphere_loadbalance(mesh, sph_r, x0, sph_l);

  /* now p-refine the patches listed in amr_prefine_p */
  if(GetLen(prefp) > 0)
  {
    int n[] = { Geti(Par("amr_prefine_n0")),
                Geti(Par("amr_prefine_n1")),
                Geti(Par("amr_prefine_n2")) };
    char *plist = Gets(prefp);
    char *pl, *str, *sav;

    PRFs(": calling prefine_pat for some patches\n");
    pl = strdup(plist);
    for(str=strtok_r(pl, " ", &sav); str!=NULL;
        str=strtok_r(NULL, " ", &sav))
    {
      int p = atoi(str);
      if(p>=0) prefine_pat(mesh, p, n);
    }
    free(pl);
  }

  /* load balance full mesh */
  PRFs(": calling simple_load_balance at the end\n");
  simple_load_balance(mesh);

  PRFs(": The resulting mesh is:\n");
  printmesh(mesh);
  //write_mylnodes(mesh, "mesh:", 0);

  return ret;
}


/* init neighbor info of root nodes */
int amr_set_bfaces_and_rnode_nbinfo_fnb(tMesh *mesh,
                                        int pr_bfaces, int pr_mesh)
{
  /* setup all bfaces */
  amr_set_all_bfaces(mesh);
  if(pr_bfaces)
  {
    PRFs(":\n");
    printallbfaces(mesh);
  }

  /* set mesh->myelms and eids */
  update_mesh_myelms_elm_eid_dt(mesh);

  /* now setup root node connections, i.e. setup neighbors of root nodes */
  amr_update_elm_nbinfo_if_nnbinfo_negative(mesh);
  amr_elm_nbinfo_to_elm_fnb(mesh);
  amr_elm_nbinfo_set_nnbinfo_mesh(mesh, 1); //make nnbinfo positive */
  amr_get_nbelm_elmheaders(mesh);

  if(pr_mesh)
  {
    PRFs(":\n");
    printmesh(mesh);
  }

  return 0;
}


/* setup mesh made out of boxes, and set distances from center in boxes_rc */
int setup_box_mesh(tMesh *mesh, double boxes_rc[3])
{
  int mesh_type = Par("amr_mesh_type");
  int BoxMesh_npatches = Geti(Par("amr_BoxMesh_npatches"));
  int npats = Geti(mesh_type);
  char *mesh_xc = Gets(Par("amr_mesh_xc"));
  double d = Getd(Par("amr_BoxMesh_dout"));
  double xc[]   = { 0., 0., 0. };
  double dout[] = { d, d, d };
  int dd;

  PRFs(":\n");

  /* re-set npats */
  if(npats<=0 || BoxMesh_npatches>0)
    npats = BoxMesh_npatches;

  /* do nothing if npats=0 */
  if(npats<=0) return 0;

  sscanf(mesh_xc, "%lg %lg %lg", &(xc[0]), &(xc[1]), &(xc[2]));
  //pr3v("xc", xc);printf("\n");

  mesh->dt = Getd(Par("dt"));
  mesh->time = 0.;
  mesh->iteration = 0;

  if(Getv(mesh_type, "Line"))
  {
    int dir;

    for(dd=0; dd<3; dd++) boxes_rc[dd] = dout[dd];

    /* set direction of line */
    if(Getv(mesh_type, "Dir2"))      dir = 2;
    else if(Getv(mesh_type, "Dir1")) dir = 1;
    else                             dir = 0;

    /* some rare case */
    if(Getv(mesh_type, "HalfEndPatches"))
    {
      double douto2[] = { d, d, d };
      double x[]   = { 0., 0., 0. };
      double c = 0.5 * (!(npats%2));
      int s = npats/2;

      /* end patches have only half the thickness */
      douto2[dir] = 0.5*d;

      /* left end patch */
      x[dir] = xc[dir] + 2.*d*(0 - s + c) + 0.5*d;
      add_1box_pat(mesh, x, douto2);
      boxes_rc[dir] = douto2[dir];

      /* middle patches */
      add_Nbox_pats_indir(mesh, xc, dout, npats-2, dir);
      boxes_rc[dir] += dout[dir]*(npats-2);

      /* right end patch */
      x[dir] = xc[dir] + 2.*d*(npats-1 - s + c) - 0.5*d;
      add_1box_pat(mesh, x, douto2);
      boxes_rc[dir] += douto2[dir];
    }
    else /* all npats patches are the same */
    {
      /* put npats boxes in dir */
      add_Nbox_pats_indir(mesh, xc, dout, npats, dir);
      boxes_rc[dir] *= npats;
    }
  }
  else if(Getv(mesh_type, "Plane"))
  {
    /* arrange npats into one plane */
    int rnpats = pow(npats, 0.5);
    int N[] = { rnpats, rnpats, 1 };
    arrange_box_pats_inBox(mesh, xc, dout, N);
    for(dd=0; dd<3; dd++) boxes_rc[dd] = dout[dd]*N[dd];
  }
  else
  {
    /* arrange npats into one big box */
    int crnpats = pow(npats, 0.3333333333334);
    int N[] = { crnpats, crnpats, crnpats };
    arrange_box_pats_inBox(mesh, xc, dout, N);
    for(dd=0; dd<3; dd++) boxes_rc[dd] = dout[dd]*N[dd];
  }

  for(dd=0; dd<3; dd++) printf("boxes_rc[%d]=%g\n", dd, boxes_rc[dd]);

  return 0;
}

/* a mesh with a number of cubed spheres, IN: mesh, BoxMesh_rc */
int setup_CubedSphere_mesh(tMesh *mesh, double BoxMesh_rc[3],
                           double *CubedSphere_r)
{
  int mesh_type = Par("amr_mesh_type");
  int npats = Geti(mesh_type);
  int CubedSphere_ndomains = Geti(Par("amr_CubedSphere_ndomains"));
  double rf_surf1 = 0.25;
  double rf_surf2 = 0.25;
  double BoxMesh_rc_max = max3(BoxMesh_rc[0],BoxMesh_rc[1],BoxMesh_rc[2]);
  double dc = Getd(Par("amr_CubedSphere_dc"));
  double csize = 0.375; //extent of inner cubes from center (must be below ~1/sqrt(3))
  double ssfac = Getd(Par("amr_CubedSphere_r0fac")); //DNSdata_OuterShellStart
  double obfac = Getd(Par("amr_CubedSphere_r1fac")); //DNSdata_OuterBoundary
  double r2fac = Getd(Par("amr_CubedSphere_r2fac"));
  /* stretch type in cubed spheres for outermost shell */
  int stretch = Geti(Par("amr_OuterShellStretch"));
  double rc[3];
  double ABrct[] = { -1.,1., -1.,1. };
  int nAB[] = { 1,1,1 };
  double xc[] = { 0., 0., 0. };
  int CubedSphere_domain_divide = Par("amr_CubedSphere_domain_divide");
  double AB_div_mode;
  char *domain_nAB = Gets(Par("amr_CubedSphere_domain_nAB"));
  char *mesh_xc = Gets(Par("amr_mesh_xc"));
  sscanf(mesh_xc, "%lg %lg %lg", &(xc[0]), &(xc[1]), &(xc[2]));
  sscanf(domain_nAB, "%d %d %d", &(nAB[0]), &(nAB[1]), &(nAB[2]));

  /* reset dc from BoxMesh_rc if BoxMesh_rc>0 */
  if(BoxMesh_rc_max>0.) dc = BoxMesh_rc_max;

  PRFs(":\n");

  mesh->dt = Getd(Par("dt"));
  mesh->time = 0.;
  mesh->iteration = 0;

  /* perform some par sanity checks */
  if(ssfac <= sqrt(3.))
    errorexit("ssfac=amr_CubedSphere_r0fac needs to be larger than sqrt(3)");
  if(npats>0 && CubedSphere_ndomains>0)
    errorexit("remove the number of patches from amr_mesh_type if "
              "amr_CubedSphere_ndomains>0");

  /* set AB_div_mode from par amr_CubedSphere_domain_divide */
  if(Getv(CubedSphere_domain_divide,"AB"))
    AB_div_mode = 1.;
  else if(Getv(CubedSphere_domain_divide,"angles"))
    AB_div_mode = 0.;
  else if(Getv(CubedSphere_domain_divide,"AB_angles"))
    AB_div_mode = Getd(CubedSphere_domain_divide);
  else
    errorexits("unknown amr_CubedSphere_domain_divide = %s",
               Gets(CubedSphere_domain_divide));

  /* setup cubed spheres based on npats read from amr_mesh_type */
  switch(npats)
  {
  case 0:
    /* do nothing if amr_mesh_type contains no number */
    break;
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
  case 5:
    rc[0] = rc[1] = rc[2] = dc;
    /* use only the 4 cubed spheres in the xy-plane: */
    sphere_around_full_box_at_xc(mesh,4, xc, rc, ssfac*dc);
    break;
  case 6:
    rc[0] = rc[1] = rc[2] = dc;
    sphere_around_empty_box_at_xc(mesh,6, xc, rc, ssfac*dc);
    break;
  case 7:
    rc[0] = rc[1] = rc[2] = dc;
    sphere_around_full_box_at_xc(mesh,6, xc, rc, ssfac*dc);
    break;
  case 12:
    rc[1] = rc[2] = dc; //dc*0.5;
    rc[0] = dc;
    two_spheres_around_empty_box_at_xc(mesh,6, xc,
                                       rc, ssfac*dc, obfac*dc, stretch);
    break;
  case 13:
    rc[1] = rc[2] = dc; //dc*0.5;
    rc[0] = dc;
    two_spheres_around_box_at_xc(mesh,6, xc, rc, ssfac*dc, obfac*dc, stretch);
    break;
  /* 13 patches but with 2 centers as in sgrid:
  case 13:
    xc[1] = xc[2] = 0.0;
    xc[0] = dc;
    arrange_1pat12CubSph_into_full_cube(mesh,6, xc,
                                        csize*rf_surf1, rf_surf1, dc);
    break; */
  case 19:
    rc[1] = rc[2] = rc[0] = dc;
    three_spheres_around_box_at_xc(mesh,6, xc, rc, ssfac*dc, obfac*dc,
                                   r2fac*dc, stretch);
    break;
  case 26:
    two_full_cubes_touching_at_x0(mesh,6, dc,
                                  csize*rf_surf1, rf_surf1,
                                  csize*rf_surf2, rf_surf2);
    break;
  case 32:
    sphere_around_two_full_cubes_touching_at_x0(mesh,6, dc,
                                                csize*rf_surf1, rf_surf1,
                                                csize*rf_surf2, rf_surf2,
                                                ssfac*dc);
    break;
  case 38:
    two_spheres_around_two_full_cubes(mesh,6, dc,
                                      csize*rf_surf1, rf_surf1,
                                      csize*rf_surf2, rf_surf2,
                                      ssfac*dc, obfac*dc);
    break;
  default:
    errorexiti("amr_mesh_type = %d CubedSpheres  <--not implemented", npats);
  }

  /* setup cubed spheres based on CubedSphere_ndomains */
  switch(CubedSphere_ndomains)
  {
  case 0:
    /* do nothing */
    break;
  case 6:
    //sphere_around_empty_box_at_xc(mesh,6, xc, BoxMesh_rc, ssfac*dc);
    sphere_nAB_around_empty_box_at_xc(mesh,6, xc, BoxMesh_rc, ssfac*dc,
                                      AB_div_mode, nAB);
    *CubedSphere_r = ssfac*dc;
    break;
  case 12:
    errorexit("we should use Shell in amr_mesh_type");
    two_spheres_around_empty_box_at_xc(mesh,6, xc, BoxMesh_rc,
                                       ssfac*dc, obfac*dc, stretch);
    break;
  case 18:
    errorexit("we should use Shell in amr_mesh_type");
    three_spheres_around_empty_box_at_xc(mesh,6, xc, BoxMesh_rc,
                                         ssfac*dc, obfac*dc,
                                         r2fac*dc, stretch);
    break;
  }

  /*
  coordinates_init(mesh);
  outputPatchPlanes_meshvar(mesh, "x", 0,0);
  outputPatchPlanes_meshvar(mesh, "y", 0,0);
  outputPatchPlanes_meshvar(mesh, "z", 0,0);
  //exit(9);
  */
  return 0;
}


/* a shell made out of a number of cubed spheres */
int setup_Shell_mesh(tMesh *mesh, double CubedSphere_r)
{
  double rin  = Getd(Par("amr_Shell_rin"));
  double r1   = Getd(Par("amr_Shell_r1"));
  double r2   = Getd(Par("amr_Shell_r2"));
  double rout = Getd(Par("amr_Shell_rout"));
  int nr = Geti(Par("amr_Shell_nr"));

  /* stretch type in cubed spheres for outermost shell */
  int stretch = Geti(Par("amr_OuterShellStretch"));
  double xc[] = { 0., 0., 0. };
  char *mesh_xc = Gets(Par("amr_mesh_xc"));
  sscanf(mesh_xc, "%lg %lg %lg", &(xc[0]), &(xc[1]), &(xc[2]));

  /* reset rin from CubedSphere_r if CubedSphere_r>0 */
  if(CubedSphere_r>0.) rin = CubedSphere_r;

  PRFs(":\n");

  mesh->dt = Getd(Par("dt"));
  mesh->time = 0.;
  mesh->iteration = 0;

  /* setup cubed spheres in form of a shell */
  if(nr == 1)
  {
    if(r1 <= rin)
    {
      CubedSphere_shell_at_xc(mesh,6, xc, rin, rout);
    }
    else
    {
      if(r2 <= r1)
        two_CubedSphere_shells_at_xc(mesh,6, xc, rin,r1,rout, stretch);
      else
        three_CubedSphere_shells_at_xc(mesh,6, xc, rin,r1,r2,rout, stretch);
    }
  }
  else /* divide radial dir between rin and rout into nr pieces */
  {
    CubedSphere_shells_at_xc(mesh,6, xc, nr, rin,rout);
  }

  return 0;
}


/***************************************************************************/
/* all the functions below are only for elm testing and can be
   removed later */
/***************************************************************************/

int setup_elm_mesh1(tMesh *mesh)
{
  int amr_n0 = Geti(Par("amr_n0"));
  int amr_n1 = Geti(Par("amr_n1"));
  int amr_n2 = Geti(Par("amr_n2"));
  int n[] = { amr_n0, amr_n1, amr_n2 };
  int pt_typ[] = { P_LGL, P_LGL, P_LGL };
  double bbox0[6] = { -4,4, -2,2, -1,1 };
  double bbox1[6] = { -4,0,  2,4, -1,1 };
  double bbox2[6] = {  0,4,  2,4, -1,1 };
  double bbox3[6] = {  4,8,  2,4, -1,1 };
  tRef ref[1];
  struct list_head *pos;
  struct list_head fnb_head;
  tElm *elm;
  int i;
  ref->method = GIVEN_n_P_UNIFORM; /* use uniform grid spacing */

  PRFs(":\n");

  mesh->dt = Getd(Par("dt"));
  mesh->time = 0.;
  mesh->iteration = 0;

  Setd(amr->sibl1to7_weight, 1.);

  add_patch(mesh, bbox0, pt_typ, n, 0);
  add_patch(mesh, bbox1, pt_typ, n, 0);
  add_patch(mesh, bbox2, pt_typ, n, 0);
  add_patch(mesh, bbox3, pt_typ, n, 0);

  /* setup all bfaces and root node connections */
  amr_set_bfaces_and_rnode_nbinfo_fnb(mesh, 1,1);

  /* we can already use mesh->myelm[2] */
  //enablevar(mesh, Ind("advection1_u"));
  if(nMPI_rank()==0)
    enablevarcomp_innode(mesh->myelm[2], Ind("advection1_u"));

  simple_load_balance(mesh);

  Yo(200);
  printmyelms(mesh);

Yo(33);
  /* refine!!! */
  formyelms(mesh)
  {
    elm = MyElm;
    if(elm->eploc->eid==2) elm->rflag = 1;
  }
  ref->method = PARENT_n;
  hrefine_elms_if_rflag(mesh, ref);
  printmyelms(mesh);

Yo(33.1);
  /* refine again !!! */
  formyelms(mesh)
  {
    elm = MyElm;
    if(elm->eploc->eid==7) elm->rflag = 1;
  }
  ref->method = PARENT_n;
  hrefine_elms_if_rflag(mesh, ref);
  printmyelms(mesh);

  simple_load_balance(mesh);
  printmyelms(mesh);

Yo(33.2);
  /* refine again !!! */
  formyelms(mesh)
  {
    elm = MyElm;
    if(elm->eploc->eid==7) elm->rflag = 1;
  }
  ref->method = PARENT_n;
  hrefine_elms_if_rflag(mesh, ref);
  printmyelms(mesh);

  simple_load_balance(mesh);
  printmyelms(mesh);



  printf("mesh->myelm:\n");
  printelmarray(mesh->nmyelm, mesh->myelm);
  printmyelms(mesh);

  printf("mesh->nbelm:\n");
  printelmarray(mesh->nnbelm, mesh->nbelm);
  printnbelms(mesh);


  /* print a var */
  /*
  list_for_each(pos, &mesh->myelm_head)
  {
    elm = list_entry(pos, tElm, list);
    printvar_innode(elm, Ind("advection1_u"));
  }
  */


/*
  printmesh(mesh);
//nMPI
MCK( nMPI_barrier() );
RunFun(FINALIZE);
finalize_all_and_exit(mesh, 0); //<--exit code 0
*/




  // try to find nb
  INIT_LIST_HEAD(&fnb_head);
  elm = mesh->myelm[6];
  if(nMPI_rank()==0)
    amr_make_fnb_list(elm, 1, mesh->nmyelm, mesh->myelm, &fnb_head);
  printf("%zu in fnb_head\n", list_count_nodes(&fnb_head));
  list_for_each(pos, &fnb_head)
  {
    elm = glist_entry(pos);
    printelm(elm);
  }
  glist_free_elems(&fnb_head);

  test_eploc();



  tArray *arr = alloc_array1d((18*sizeof(tEploc))/sizeof(double));
  ulong ui = 3;
  ulong nef = 1*(255<<24) + (3<<16) + 11;
  tEploc eploc[7];
  tEloc eloc0;
  tEloc eloc;

  eloc0.p=16;
  eloc0.l=7;
  eloc0.eid=16;
  for(i=0; i<NLOCS; i++) eloc0.loc[i]=5 + '0';
  //strcpy(eloc0.loc, "1234567");
  eloc_to_eploc(&eloc0, &eploc[0]);
  eloc_from_eploc(&eloc, &eploc[0]);
  printeloc_s(&eloc, "\n");

  redim_array(arr, 50,1,1);
  forarray(arr, i) arr->d[i] = 1;

  printf("nef=%lu\n", nef);
  memcpy_to_array_redim(arr, sizeof(tEploc), ui,
                        &nef, 8);
  printarray(arr);

  memcpy_to_array_redim(arr, sizeof(tEploc), ui,
                        &eploc[0], sizeof(tEploc));
  printarray(arr);
  printf("from arr: ");
  eloc_from_eploc(&eloc, &(arr->eploc[ui]));
  printeloc_s(&eloc, "\n");
  printf("arr");printarray_eploc(arr, 1);

  elm = mesh->myelm[6];
  for(i=1; i<6; i++) eploc[i] = eploc[0];
  for(i=1; i<6; i++) eploc[i].l = i;

  amr_elm_nbinfo_add_nbeploc(elm, 3, 3,&eploc[0]);
  printf("amr_elm_nbinfo3: ");print_amr_elm_nbinfo(elm, 3);printf("\n");
  amr_elm_nbinfo_add_nbeploc(elm, 3, 2,&eploc[3]);
  printf("amr_elm_nbinfo3: ");print_amr_elm_nbinfo(elm, 3);printf("\n");
  //amr_elm_nbinfo_add_nbeploc(elm, 3, 2,&eploc[1]);
  //printf("amr_elm_nbinfo3: ");print_amr_elm_nbinfo(elm, 3);printf("\n");
  free_array(arr);

  //disablevar_innode(elm, Ind("amr_elm_nbinfo0"));


  /* set flag to update all fnb, not needed because alloc_dat does this */
  formyelms(mesh)
  {
    tElm *Elm = MyElm;
    for(int f=0; f<6; f++)
    {
      Elm->dat->info->nnbinfo[f] = -1; //make nnbinfo negative
      //Elm->nfnb[f] = -1;
    }
  }
  amr_update_elm_nbinfo_if_nnbinfo_negative(mesh);
  amr_elm_nbinfo_to_elm_fnb(mesh);
  amr_elm_nbinfo_set_nnbinfo_mesh(mesh, 1); //make nnbinfo positive


  formyelms(mesh)
  {
    tElm *Elm = MyElm;
    for(int f=0; f<6; f++)
    {
      Elm->dat->info->nnbinfo[f] = -1; //make nnbinfo negative
      //Elm->nfnb[f] = -1;
    }
  }
  amr_update_elm_nbinfo_if_nnbinfo_negative(mesh);
  amr_elm_nbinfo_to_elm_fnb(mesh);
  amr_elm_nbinfo_set_nnbinfo_mesh(mesh, 1); //make nnbinfo positive


  printf("mesh->myelm:\n");
  printelmarray(mesh->nmyelm, mesh->myelm);
  printmyelms(mesh);

  amr_get_nbelm_elmheaders(mesh);

  printf("mesh->nbelm:\n");
  printelmarray(mesh->nnbelm, mesh->nbelm);
  printnbelms(mesh);


  tElm *e1 = mesh->myelm[6];
  tElm *e2 = mesh->myelm[5];
  tElm *ar[] = {e1, e2};
  printelm(e1);
  printelm(e2);
  printf(":::%d\n", eploc1_eploc2_agree_upto_l_max(e1->eploc, e2->eploc, 1));
  printf(":::%d\n", amr_elms_are_siblings(2, ar));

  ulong eidarr[99];
  tElm0 elm0[99] = {0};
  int narr;

  printf("amr_get_otherrank_elm0_for_eids:\n");
  if(nMPI_rank()==0)
  {
    eidarr[0] = 1;
    eidarr[1] = 2;
    eidarr[0] = 18;
    eidarr[1] = 17;
    eidarr[2] = 3;
    narr = 3;
    amr_get_elm0_for_eids(mesh, narr, eidarr, elm0);
    for(i=0; i<narr; i++) printelm0(&elm0[i],"\n");
  }
  else
  {
    eidarr[0] = 17;
    eidarr[1] = 18;
    eidarr[2] = 16;
    eidarr[0] = 2;
    eidarr[1] = 3;
    eidarr[2] = 1;
    eidarr[3] = 17;
    narr = 4;
    amr_get_elm0_for_eids(mesh, narr, eidarr, elm0);
    for(i=0; i<narr; i++) printelm0(&elm0[i],"\n");
  }


  Setd(amr->sibl1to7_weight, 0.);
  load_balance(mesh, 1);
  Setd(amr->sibl1to7_weight, 1.);

  printf("mesh->myelm:\n");
  printelmarray(mesh->nmyelm, mesh->myelm);
  printmyelms(mesh);

  printf("mesh->nbelm:\n");
  printelmarray(mesh->nnbelm, mesh->nbelm);
  printnbelms(mesh);



  formyelms(mesh)
  {
    elm = MyElm;
    //if(Elm_l(elm)>2) elm->rflag = -1;
    //if(Elm_eid(elm)==11) elm->rflag = 0;
    if(Elm_eid(elm)>=2 && Elm_eid(elm)<=16) elm->rflag = -1;
  }
  ref->method = PARENT_nO2_P1; // PARENT_n;
  remove_elms_if_rflag(mesh, ref);
  printmyelms(mesh);

  simple_load_balance(mesh);
  printmyelms(mesh);





//nMPI
MCK( nMPI_barrier() );
RunFun(FINALIZE);
finalize_all_and_exit(mesh, 0); //<--exit code 0
  return 0;
}




/***************************************************************************/
/* all the functions below were just for very early testing and could be
   removed now */
/***************************************************************************/

/* set up a mesh with 2 levels  */
int setup_3patchl2_mesh(tMesh *mesh)
{
  int amr_n0 = Geti(Par("amr_n0"));
  int amr_n1 = Geti(Par("amr_n1"));
  int amr_n2 = Geti(Par("amr_n2"));
  int n[] = { amr_n0, amr_n1, amr_n2 };
  int pt_typ[] = { P_LGL, P_LGL, P_LGL };
  double bbox0[6] = { -4,4, -2,2, -1,1 };
  double bbox1[6] = { -4,0,  2,4, -1,1 };
  double bbox2[6] = {  0,4,  2,4, -1,1 };

  PRFs(":\n");

  mesh->dt = Getd(Par("dt"));
  mesh->time = 0.;
  mesh->iteration = 0;

  add_patch(mesh, bbox0, pt_typ, n, 0);
  add_patch(mesh, bbox1, pt_typ, n, 0);
  add_patch(mesh, bbox2, pt_typ, n, 0);

  /* setup all bfaces and root node connections */
  amr_set_bfaces_and_rnode_nbinfo_fnb(mesh, 1,1);

  errorexit("add back the stuff below");
  /* 8 children in patch0 */
  hrefine_pat(mesh, 0);

  /* 8 more in each patch */
  hrefine_pat(mesh, 0);
  hrefine_pat(mesh, 1);
  hrefine_pat(mesh, 2);

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
