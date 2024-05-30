/* nmesh_amr.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "amr.h"



int nmesh_amr(tMesh *mesh)
{
  printf("Adding amr\n");

  /* functions */
  /* NOTE: amr_setup_mesh(mesh) is called directly from main */
  AddFun(POST_INITLIBS, amr_print_thread_info);
  //AddFun(POST_INITLIBS, timing_mm_speed); //this is too inaccurate (& slow)
  AddFun(POST_PARAMETERS, amr_init_global_pars);
  AddFun(POST_INITMESH, amr_set_use_fv_flag);
  AddFun(LOADBALANCING, load_balance_if_needed);

  /* variables */
  AddAuxVar("amr_elm_nbinfo", "@", "vars with tEploc arrays for nbinfo[6]");
  MeshVarSetType(mesh, Ind("amr_elm_nbinfo0"), LBTVAR); //transf. in loadbal.
  //AddAuxVar("Xb", "", "coord0 inside each node"); // we may not need these???
  //AddAuxVar("Yb", "", "coord1 inside each node");
  //AddAuxVar("Zb", "", "coord2 inside each node");

  /* parameters */
  AddPar("amr_mesh_type", "test_mesh",
	 "initial mesh [BoxMesh,CubedSpheres,Shell,test_mesh,l2_mesh]");
  AddPar("amr_n0", "5", "number of points in dir. 0 in one node");
  AddPar("amr_n1", "5", "number of points in dir. 1 in one node");
  AddPar("amr_n2", "5", "number of points in dir. 2 in one node");
  AddPar("amr_nmax", "55", "max number of points in all 3 dir. in one node");
  AddPar("amr_mesh_xc", "0 0 0", "center for boxes or spheres");
  AddPar("amr_BoxMesh_npatches", "0", "number of BoxMesh patches "
         "(0 means par is no set");
  AddPar("amr_BoxMesh_dout", "1", "box radius");
  AddPar("amr_CubedSphere_n0", "-1", "number of points in dir. 0");
  AddPar("amr_CubedSphere_n1", "-1", "number of points in dir. 1");
  AddPar("amr_CubedSphere_n2", "-1", "number of points in dir. 2");
  AddPar("amr_Shell_n0", "-1", "number of points in dir. 0");
  AddPar("amr_Shell_n1", "-1", "number of points in dir. 1");
  AddPar("amr_Shell_n2", "-1", "number of points in dir. 2");
  AddPar("amr_CubedSphere_ndomains", "0", "number of CubedSphere domains");
  AddPar("amr_CubedSphere_domain_nAB", "1 1 1", "number of domain divisions "
         "in x, y, and z-dirs");
  AddPar("amr_CubedSphere_domain_divide", "angles", "how we divide domain "
         "[AB,angles,0.8 AB_angles]");
  AddPar("amr_CubedSphere_dc", "0.5", "length scale for cubed sphere setup");
  AddPar("amr_CubedSphere_r0fac", "4", "size factor for outer radius r0 of "
         "cubed sphere patches: r0 = r0fac*dc");
  AddPar("amr_CubedSphere_r1fac", "6", "size factor for outer radius r1 of "
         "cubed sphere patches: r1 = r1fac*dc");
  AddPar("amr_CubedSphere_r2fac", "10", "size factor for outer radius r2 of "
         "cubed sphere patches: r2 = r2fac*dc");
  AddPar("amr_OuterShellStretch", "no", "if and how we stretch outer "
         "shell [0,1,2]");
  AddPar("amr_Stretch_w", "1", "width par for amr_OuterShellStretch");
  AddPar("amr_Stretch_A", "1", "amplitude par for amr_OuterShellStretch");
  AddPar("amr_Shell_rin", "0.5", "inner radius of shell");
  AddPar("amr_Shell_r1",  "0" ,  "intermediate radius1 [for 0 it's off]");
  AddPar("amr_Shell_rout", "1", "outer radius of shell");
  AddPar("amr_uniform_p", "", "list of patches with uniform grid spacing");
  AddPar("amr_fv_p", "", "patches where we use finite volume");
  AddPar("amr_fv_if_uniform", "no", "set use_fv in uniform elms [no,yes]");
  if(Getb(Par("amr_fv_if_uniform")))
    AddFun(PRE_COORDINATES, amr_use_fv_if_P_UNIFORM);

  /* pars that control how many dimensions we use */
  AddPar("amr_dir_active0", "yes", "whether dir. 0 is used [yes,no]");
  AddPar("amr_dir_active1", "yes", "whether dir. 1 is used [yes,no]");
  AddPar("amr_dir_active2", "yes", "whether dir. 2 is used [yes,no]");

  /* refinement related pars for initial mesh creation */
  AddPar("amr_luni", "0",  "level up to which each patch is refined initially");
  AddPar("amr_hrefine_p", "", "patch list that we h-refine 1 level further");
  AddPar("amr_prefine_p", "", "patch list that we p-refine");
  AddPar("amr_prefine_n0", "10", "n0 if amr_prefine_p is used");
  AddPar("amr_prefine_n1", "10", "n1 if amr_prefine_p is used");
  AddPar("amr_prefine_n2", "10", "n2 if amr_prefine_p is used");
  AddPar("amr_hrefine_sphere_levels", "0",
         "number of nested sphere h-refinement levels");
  AddPar("amr_hrefine_sphere_radius", "10",
         "radius of innermost sphere in nested sphere h-refinement");

  /* pars that determine how load is balanced */
  AddPar("amr_load_balance", "no", "how load is balanced [no,yes]");
  AddPar("amr_sibl1to7_weight", "0", "if 0 all 8 siblings are kept on the "
         "same MPI rank, else they can end up on different ranks [0,1]");

  /* MPI related pars */
  AddPar("amr_MPIexchange", "1", "type of MPI exchange we use: "
         "1: surfaces, 2: ghosts");
  /* ghost related pars */
  AddPar("amr_nghosts", "1", "number of ghost zones");
  AddPar("amr_N0", "-1", "total number of inner points in dir. 0 "
         "[any positive #, -1 means don't use it]");
  AddPar("amr_N1", "-1", "total number of inner points in dir. 1");
  AddPar("amr_N2", "-1", "total number of inner points in dir. 2");

  /* pars related to finding neighbors */
  AddPar("amr_nbsearch_n", "6", "num. of points n used to search for nbs in "
         "func find_elmfacepoints_in_nbface");

  /* pars related to interpolation during mesh refinement */
  AddPar("amr_Lagrange_interp_order", "n", "interp order [n,4,6,8,...]");
  AddPar("amr_WENO_interp_order", "6", "interpolation order [4,6]");

  /* bface pars */
  AddPar("bface_options", "face2_order3", "how we set some bface flags "
         "[none,face2_order0,face2_order1,face2_order2,face2_order3]");

  /* Old parameters that are now banned */
  BanPar("amr_BoxMesh_xc", "use amr_mesh_xc instead");
  BanPar("amr_refine_p", "use amr_hrefine_p instead");
  BanPar("amr_refine_sphere_levels", "use amr_hrefine_sphere_levels instead");
  BanPar("amr_refine_sphere_radius", "use amr_hrefine_sphere_radius instead");

  return 0;
}
