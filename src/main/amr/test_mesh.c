/* test_mesh.c */
/* Wolfgang Tichy, 5/2023 */


#include "nmesh.h"
#include "amr.h"

#define PR 0



/* refine a randomly chosen fraction frac of nodes, do this N times */
void random_refine_frac_nodes_N_times(tMesh *mesh, double frac, int N)
{
  uint32_t ran=1; // seed for rand_32primitive
  tRef ref[1];
  int i, j;

  /* do N refines */
  for(i=0; i<N; i++)
  {
    long nnodes = mesh->myln->nm;

    /* set rflag Nr times */
    for(j=0; j<nnodes*frac+1; j++)
    {
      long myid = nnodes * rand_32primitive_u01(&ran);
      tNode *node = Lnode_myid(mesh, myid);
      node->rflag = 1;
    }
    ref->type = rand_32primitive_u01(&ran) * 1.99999;
    ref->method = GIVEN_n;
    ref->n[0] = rand_32primitive_u01(&ran) * 6;
    ref->n[1] = rand_32primitive_u01(&ran) * 6;
    ref->n[2] = rand_32primitive_u01(&ran) * 6;
    hp_refine_nodes_if_rflag(mesh, ref);
    update_mesh_myln_node_nid(mesh);
  }
}


/* make some patches, refine and unrefine several times and output
   the resulting meshes */
int test_mesh(tMesh *mesh)
{
  int mode=1;
  PRFs(":\n");

  mesh->dt = Getd(Par("dt"));
  mesh->time = 0.;
  mesh->iteration = 0;

  /* remove all patches from mesh, so we can just add new pristine ones */
  remove_all_patches(mesh);

  /* from cubed sph tests */
  two_diff_wegdes_touching_1_wedge(mesh, 1.0, 3.0, 5.0);

  /* setup all bfaces and root node connections */
  amr_set_bfaces_and_rnode_nfaces_fnb(mesh, 1);

  /* write initial mesh into files */
  write_mynodeelms(mesh, "initial patches:", mode);

  /* 1st random refine */
  random_refine_frac_nodes_N_times(mesh, 0.2, 6);
  write_mynodeelms(mesh, "ref1:", mode);

  return 0;
}
