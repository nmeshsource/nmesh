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
  int i;

  /* do N refines */
  for(i=0; i<N; i++)
  {
    /* set rflag on a fraction */
    formylnodes(mesh)
    {
      tNode *node = MyLnode;
      node->rflag = rand_32primitive_u01(&ran) + frac;
    }
    ref->type = rand_32primitive_u01(&ran) * 1.99999;
    ref->method = GIVEN_n;
    ref->n[0] = rand_32primitive_u01(&ran) * 6 + 1;
    ref->n[1] = rand_32primitive_u01(&ran) * 6 + 1;
    ref->n[2] = rand_32primitive_u01(&ran) * 6 + 1;
    //hp_refine_nodes_if_rflag(mesh, ref);
    if(ref->type==0) hrefine_nodes_if_rflag(mesh, ref);
    else             prefine_nodes_if_rflag(mesh, ref);
    update_mesh_myln_node_nid(mesh);
  }
}

/* unrefine a randomly chosen fraction frac of nodes, do this N times */
void random_remove_frac_nodes_N_times(tMesh *mesh, double frac, int N)
{
  uint32_t ran=1; // seed for rand_32primitive
  tRef ref[1];
  int i, j, rflag;

  /* do N unrefines */
  for(i=0; i<N; i++)
  {
    /* set rflag on a fraction */
    j=0;
    formylnodes(mesh)
    {
      tNode *node = MyLnode;
      if(j%8==0) rflag = -(rand_32primitive_u01(&ran) + frac);
      node->rflag = rflag * (Node_l(node)>0);
      j++;
    }
    //ref->type = rand_32primitive_u01(&ran) * 1.99999;
    ref->method = GIVEN_n;
    ref->n[0] = rand_32primitive_u01(&ran) * 6 + 1;
    ref->n[1] = rand_32primitive_u01(&ran) * 6 + 1;
    ref->n[2] = rand_32primitive_u01(&ran) * 6 + 1;
    remove_nodes_if_rflag(mesh, ref);
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
  update_mesh_myln_node_nid(mesh);
  simple_load_balance(mesh);
  update_mesh_myln_node_nid(mesh);

  /* write initial mesh into files */
  write_mylnodes(mesh, "initial patches:", mode);

  /* 1. random refine */
  random_refine_frac_nodes_N_times(mesh, 0.2, 12);
  simple_load_balance(mesh);
  update_mesh_myln_node_nid(mesh);
  write_mylnodes(mesh, "ref1:", mode);

  /* 2. random unrefine */
  random_remove_frac_nodes_N_times(mesh, 0.2, 12);
  simple_load_balance(mesh);
  update_mesh_myln_node_nid(mesh);
  write_mylnodes(mesh, "ref2:", mode);

  return 0;
}
