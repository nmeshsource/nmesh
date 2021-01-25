/* nmesh_amr.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "amr.h"



int nmesh_amr(tMesh *mesh)
{
  printf("Adding amr\n");

  /* functions */
  /* NOTE: amr_setup_mesh(mesh) is called directly from main */
  AddFun(FIRST, amr_print_thread_info);
  //AddFun(FIRST, timing_mm_speed); // this is too slow right now
  AddFun(POST_PARAMETERS, amr_init_global_pars);
  AddFun(INITIALDATA, amr_set_use_fv_flag);
  AddFun(LOADBALANCING, load_balance_if_needed);

  /* variables */
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
  AddPar("amr_BoxMesh_xc", "0 0 0", "center for boxes");
  AddPar("amr_BoxMesh_dout", "1", "box radius");
  AddPar("amr_CubedSphere_dc", "0.5", "length scale for cubed sphere setup");
  AddPar("amr_CubedSphere_r0fac", "4", "size factor for outer radius r0 of "
         "cubed sphere patches: r0 = r0fac*dc");
  AddPar("amr_CubedSphere_r1fac", "6", "size factor for outer radius r1 of "
         "cubed sphere patches: r1 = r1fac*dc");
  AddPar("amr_OuterShellStretch", "no", "if and how we stretch outer "
         "shell [0,1,2]");
  //AddPar("amr_Stretch_w", "1", "width par for amr_OuterShellStretch");
  //AddPar("amr_Stretch_A", "1", "amplitude par for amr_OuterShellStretch");
  AddPar("amr_Shell_rin", "0.5", "inner radius of shell");
  AddPar("amr_Shell_rout", "1", "outer radius of shell");
  AddPar("amr_uniform_p", "", "list of patches with uniform grid spacing");
  AddPar("amr_fv_p", "", "patches where we use finite volume");

  /* pars that control how many dimensions we use */
  AddPar("amr_dir_active0", "yes", "whether dir. 0 is used [yes,no]");
  AddPar("amr_dir_active1", "yes", "whether dir. 1 is used [yes,no]");
  AddPar("amr_dir_active2", "yes", "whether dir. 2 is used [yes,no]");

  /* refinement related pars for initial mesh creation */
  AddPar("amr_luni", "0",  "level up to which each patch is refined initially");
  AddPar("amr_refine_p", "", "patch list that we refine one level further");
  AddPar("amr_refine_sphere_levels", "0",
         "number of nested sphere refinement levels");
  AddPar("amr_refine_sphere_radius", "10",
         "radius of innermost sphere in nested sphere refinement");

  /* pars that determine how load is balanced */
  AddPar("amr_load_balance", "no", "[no,simple,timingbased]");

  /* MPI related pars */
  AddPar("amr_MPIexchange", "1", "type of MPI exchange we use: "
         "1: surfaces, 2: ghosts");
  /* ghost related pars */
  AddPar("amr_nghosts", "1", "number of ghost zones");
  AddPar("amr_N0", "-1", "total number of inner points in dir. 0 "
         "[any positive #, -1 means don't use it]");
  AddPar("amr_N1", "-1", "total number of inner points in dir. 1");
  AddPar("amr_N2", "-1", "total number of inner points in dir. 2");

  /* bface pars */
  AddPar("bface_options", "face2_order3", "how we set some bface flags "
         "[none,face2_order0,face2_order1,face2_order2,face2_order3]");

  return 0;
}
