/* test_mesh.c */
/* Wolfgang Tichy, 5/2023 */


#include "nmesh.h"
#include "amr.h"

#define PR 0


extern tAMR amr[1];


/* refine a randomly chosen fraction frac of nodes, do this N times */
void random_refine_frac_nodes_N_times(tMesh *mesh, double frac, int N)
{
  int size = nMPI_size();
  int rank = nMPI_rank();
  uint32_t ran1=1; // seed1 for rand_32primitive
  uint32_t ran2=1; // seed2 for rand_32primitive
  tRef ref[1];
  int i, rk;

  /* do N refines */
  for(i=0; i<N; i++)
  {
    for(rk=0; rk<size; rk++)
    {
      if(rk==rank)
      {
        /* set rflag on a fraction */
        formylnodes_noomp(mesh)
        {
          tNode *node = MyLnode;
          /* ran2 depends on the number of nodes on this proc */
          node->rflag = rand_32primitive_u01(&ran2) + frac;
        }
      }
      /* broadcast ran2 to other ranks */
      MCK( nMPI_Bcast(&ran2,1, nMPI_DOUBLE, rk) );
    }
    /* All MPI procs must do the lines below in the same way, so ran1 must
       be the same no matter how many procs we have! */
    ref->type = rand_32primitive_u01(&ran1) * 1.99999;
    ref->method = GIVEN_n;
    ref->n[0] = rand_32primitive_u01(&ran1) * 6 + 1;
    ref->n[1] = rand_32primitive_u01(&ran1) * 6 + 1;
    ref->n[2] = rand_32primitive_u01(&ran1) * 6 + 1;
    //hp_refine_nodes_if_rflag(mesh, ref);
    if(ref->type==0) hrefine_nodes_if_rflag(mesh, ref);
    else             prefine_nodes_if_rflag(mesh, ref);
    update_mesh_myln_node_nid(mesh);
//printmyelms(mesh);
//abort();
    simple_load_balance(mesh);
    update_mesh_myln_node_nid(mesh);
  }
}

/* unrefine a randomly chosen fraction frac of nodes, do this N times */
void random_remove_frac_nodes_N_times(tMesh *mesh, double frac, int N)
{
  int size = nMPI_size();
  int rank = nMPI_rank();
  uint32_t ran1=1; // seed1 for rand_32primitive
  uint32_t ran2=1; // seed2 for rand_32primitive
  tRef ref[1];
  int i, rk, j, rflag;

  /* do N unrefines */
  for(i=0; i<N; i++)
  {
    //prdivider(3);
    //prdivider(3);
    //prdivider(3);
    //prdivider(3);
    //PRF;printf(": i=%d\n", i);
    j=0;
    for(rk=0; rk<size; rk++)
    {
      if(rk==rank)
      {
        /* set rflag on a fraction */
        formylnodes_noomp(mesh)
        {
          tNode *node = MyLnode;
          /* ran2 depends on the number of nodes on this proc */
          if(j%8==0) rflag = -(rand_32primitive_u01(&ran2) + frac);
          node->rflag = rflag * (Node_l(node)>0);
          j++;
        }
      }
      /* broadcast ran2 and j to other ranks */
      MCK( nMPI_Bcast(&ran2,1, nMPI_DOUBLE, rk) );
      MCK( nMPI_Bcast(&j,1, nMPI_INT, rk) );
      MCK( nMPI_Bcast(&rflag,1, nMPI_INT, rk) );
    }
    /* All MPI procs must do the lines below in the same way, so ran1 must
       be the same no matter how many procs we have! */
    //ref->type = rand_32primitive_u01(&ran1) * 1.99999;
    ref->method = GIVEN_n;
    ref->n[0] = rand_32primitive_u01(&ran1) * 6 + 1;
    ref->n[1] = rand_32primitive_u01(&ran1) * 6 + 1;
    ref->n[2] = rand_32primitive_u01(&ran1) * 6 + 1;
    remove_elms_if_rflag(mesh, ref);
    update_mesh_myln_node_nid(mesh);
    simple_load_balance(mesh);
    update_mesh_myln_node_nid(mesh);
  }
}


/* make some patches, refine and unrefine several times and output
   the resulting meshes */
int test_mesh(tMesh *mesh)
{
  int mode=0;
  tRef ref[1];
  int iX, vi;

  PRFs(":\n");

  AddEvoVar("test_mesh_var", "", "some test var");
  vi = Ind("test_mesh_var");
  iX = Ind("X");

  mesh->dt = Getd(Par("dt"));
  mesh->time = 0.;
  mesh->iteration = 0;

  /* give equal weight to all siblings */
  Setd(amr->sibl1to7_weight, 1.);

  /* from cubed sph tests */
  two_diff_wegdes_touching_1_wedge(mesh, 1.0, 3.0, 5.0);

  /* setup all bfaces and root node connections */
  amr_set_bfaces_and_rnode_nbinfo_fnb(mesh, 1);
  update_mesh_myln_node_nid(mesh);

  /* refine all once to have more nodes */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    node->rflag = 1;
  }
  ref->method = PARENT_n;
  hrefine_nodes_if_rflag(mesh, ref);
  update_mesh_myln_node_nid(mesh);
  simple_load_balance(mesh);
  update_mesh_myln_node_nid(mesh);

  /* write some data in test_mesh_var */
  enablevar(mesh, vi);
  formyelms(mesh)
  {
    tElm *elm = MyElm;
    int *n = elm->n;
    int i,j,k;
    double *X = Vard(elm,iX);
    double *Y = Vard(elm,iX+1);
    double *Z = Vard(elm,iX+2);
    double *v = Vard(elm,vi);

    forijk(i,j,k, n)
    {
      int ijk = Ind_n(i,j,k, n);
      double Xi = X[ijk];
      double Yi = Y[ijk];
      double Zi = Z[ijk];
      v[ijk] = Xi*Xi*Xi*Xi + 2*Yi*Yi*Yi + 3*Zi*Zi;
    }
  }

  /* write initial mesh into files */
  write_mylnodes(mesh, "initial patches:", mode);

  /* 1. random refine */
  random_refine_frac_nodes_N_times(mesh, 0.2, 8);
  write_mylnodes(mesh, "ref1:", mode);

  //does not work with old AMR:
  /* 2. random unrefine */
  random_remove_frac_nodes_N_times(mesh, 0.4, 9);
  write_mylnodes(mesh, "ref2:", mode);
  /* write_mylnodes crashes while printing a fnb that is messed
     up (e.g. pat=NULL). Probably remove_nodes_if_rflag has a bug!!! */

  return 0;
}
