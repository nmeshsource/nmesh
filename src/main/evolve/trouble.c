/* evolve.c */
/* Wolfgang Tichy, 5/2022 */

#include "nmesh.h"
#include "evolve.h"

#define PR 0

/* Determine and set trouble score in a node,
   i.e. set node->dat->info->trouble .
   Note evosys->f[TROUBLE] returns a score of 1, 0, or -1:
    1 trouble  =>  switch to more robust method (e.g. fv)
    0 node ok with current method  =>  do nothing
   -1 all ok   =>  switch to more accurate but delicate method (e.g. dg)
   examples: if shock in dg return 1
             if shock in fv return 0
             if no shock in dg return 0
             if no shock in fv return -1
   Here we also accumulate the 1 and -1 scores and return them */
int evolve_set_trouble_score(tNode *node)
{
  tMesh *mesh = node->pat->mesh;
  tEvoSys *evosys = mesh->evosys;
  pVLList *u_p = evosys->u_p;
  pVLList *u   = evosys->u;
  pVLList *r   = evosys->rhs;
  int max_trouble = 1073741824;
  int troubled = 0; /* default is to assume that we change nothing */
  int tr_max = -max_trouble; /* init to low value */
  int i;

  if(node->dat == NULL) errorexit("node->dat is NULL");

  if(PR) PRFs(":\n");

  /* check all evo systems for trouble and put max into tr_max */
  forList(u, i)
  {
    tVarList *vlu_p = ListEntry(u_p,i);
    tVarList *vlu   = ListEntry(u,i);
    tVarList *vlr   = ListEntry(r,i);

    /* run TROUBLE func */
    if(tr_max<=0) /* need to check only if there no trouble yet */
    {
      if(ListEntry(evosys->f[TROUBLE],i))
      {
        int tr = 0;
        tr = ListEntry(evosys->f[TROUBLE],i)(node, vlr, vlu, vlu_p);
        if(tr>tr_max) tr_max = tr;
      }
    }
  }

  /* set troubled flag to 1, 0, or -1 */
  if(tr_max>0) /* new trouble found, switch to fv */
    troubled = 1;
  else if(tr_max==0 || tr_max==-max_trouble) /* ok, keep as is */
    troubled = 0;
  else /* all is very good */
    troubled = -1;
  //pr_nodename(node);
  //printf(" tr_max=%d troubled=%d", tr_max, troubled);

  /* set node->dat->info->trouble */
  if(troubled>0) /* i.e. there is trouble now */
  {
    /* if there was no trouble, we set the trouble score to 1 */
    if(node->dat->info->trouble<=0) node->dat->info->trouble  = 1;
    /* if there was trouble before, we continuously increase the score */
    else                            node->dat->info->trouble += 1;
  }
  else if(troubled==0) /* is ok, keep node as is */
  {
    node->dat->info->trouble = 0;
  }
  else /* all is good, can switch back to dg */
  {
    /* if there was trouble, we set the trouble score to -1 */
    if(node->dat->info->trouble>=0) node->dat->info->trouble  = -1;
    /* if there was no trouble, we continuously lower the score */
    else                            node->dat->info->trouble -= 1;
  }

  /* make sure trouble score does not become too large or small */
  if(node->dat->info->trouble < -max_trouble)
    node->dat->info->trouble = -max_trouble;
  if(node->dat->info->trouble > max_trouble)
    node->dat->info->trouble = max_trouble;

  //printf(" ->trouble=%d\n", node->dat->info->trouble);

  return node->dat->info->trouble;
}

/* determine and set trouble score in each node,
   i.e. set node->dat->info->trouble */
int evolve_set_trouble_score_mesh(tMesh *mesh)
{
  tEvoSys *evosys = mesh->evosys;
  pVLList *u_p = evosys->u_p;
  int Max_trb, max_trb=INT_MIN;
  int Min_trb, min_trb=INT_MAX;
  if(PR) PRFs(":\n");

  /* collect min/max values all nodes and their of neighbors */
  evolve_collect_u_p_data_mesh(mesh, u_p);

  /* loop over all nodes, check for trouble, and node-info trouble score */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    evolve_set_trouble_score(node);
  }

  /* now find rank-local max of node->dat->info->trouble */
  formylnodes_noomp(mesh)
  {
    tNode *node = MyLnode;
    int trb = node->dat->info->trouble;
    if(trb>max_trb) max_trb = trb; /* max trouble */
    if(trb<min_trb) min_trb = trb; /* min trouble */
  }

  /* Max over all ranks */
  Max_trb = max_trb;
  nMPI_Allreduce(&max_trb, &Max_trb, 1, nMPI_INT, nMPI_MAX);

  /* Min over all ranks */
  Min_trb = min_trb;
  nMPI_Allreduce(&min_trb, &Min_trb, 1, nMPI_INT, nMPI_MIN);

  /* if there is bad trouble somewhere */
  if(Max_trb>0)
    return Max_trb; /* returns max. trouble value of all nodes */
  else
    return Min_trb; /* no trouble, return smallest to signal need for dg */
}


/* switch bewteen fv and dg based on node->dat->info->trouble flag */
void evolve_trouble_switch_dg_fv_mesh(tMesh *mesh)
{
  tRef ref[1]; /* for ref info */
  int ptUNI[] = { P_UNIFORM, P_UNIFORM, P_UNIFORM };
  int ptLGL[] = { P_LGL, P_LGL, P_LGL };

  if(PR) PRFs(":\n");

  /* free surfaces & indc since they will change now anyway */
  evolve_free_communication_structs(mesh);

  /* set ref to uniform */
  ref->method = PARENT_n_P_UNIFORM; /* use uniform grid with same n */

  /* loop over all nodes, and flag all troubled nodes for uniform grid */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;

    /* make sure time on node is correct */
    node->time = mesh->time;

    /* mark troubled nodes as to be refined */
    if(node->dat->info->trouble > 0)
      node->rflag = ref->method;
    else
      node->rflag = 0;
  }

  /* do p-refinement to desired n and point type */
  prefine_nodes_if_rflag(mesh, ref);
  ///* also p-refine neighbors to make sure the star surface stays within fv
  //   region when the radius changes a bit */
  //prefine_nodes_if_nb_uniform_in_any_dir(mesh, ref);

  /* set ref to LGL */
  ref->method = PARENT_n_P_LGL; /* use LGL grid with same n */

  /* loop over all nodes, and flag all non-troubled nodes for LGL grid */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;

    /* mark non-troubled nodes as to be refined */
    if(node->dat->info->trouble <= -NOTROUBLES)
      node->rflag = ref->method;
    else
      node->rflag = 0;
  }

  /* do p-refinement to desired n and point type */
  prefine_nodes_if_rflag(mesh, ref);
  ///* also p-refine neighbors to make sure the star surface stays within fv
  //   region when the radius changes a bit */
  //prefine_nodes_if_nb_uniform_in_any_dir(mesh, ref);

  /* clear rflag on all leaf nodes */
  refine_set_rflag_forall_nodes(mesh, 0);

  /* update, nids won't change but hmin and thus dt might */
  update_mesh_myln_node_nid(mesh);

  /* switch on fv on all uniform nodes */
  refine_set_use_fv_if_pt_typ(mesh, ptUNI, 1);

  /* switch off fv on all LGL nodes */
  refine_set_use_fv_if_pt_typ(mesh, ptLGL, 0);

  /* now that nodes are changed re-init surfaces & indc */
  evolve_init_communication_structs(mesh);

  /* balance load, now that some nodes use a different method */
  //FIXME: do something better than simple_load_balance
}


/* set u = u_p, and then switch to fv */
void evolve_prepare_do_over_mesh(tMesh *mesh)
{
  tEvoSys *evosys = mesh->evosys;
  double  t = mesh->time;
  double dt = mesh->dt;
  pVLList *u_p = evosys->u_p;
  pVLList *u   = evosys->u;

  if(PR) PRFs(":\n");

  /* go back to t-dt and set u = u_p */
  copy_pVLList(u, u_p, vlcopy,0);
  mesh->time = t-dt;

  /* loop over all nodes, and set time */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    node->time = mesh->time;
  }

  /* switch nodes based on trouble flag */
  evolve_trouble_switch_dg_fv_mesh(mesh);
}

/* switch back to dg */
void evolve_switch_nontroubled_nodes_mesh(tMesh *mesh)
{
  /* switch nodes based on trouble flag */
  evolve_trouble_switch_dg_fv_mesh(mesh);
}

/* Set myindc for u_p to get min/max of u_p needed for a RDMP trouble indicator.
   Here we use limdata_MRS to get min,max,average. */
void evolve_collect_u_p_data_mesh(tMesh *mesh, pVLList *u_p)
{
  tEvoSys *evosys = mesh->evosys;
  int i;

  if(PR) PRFs(":\n");

  /* loop over list of varlists */
  forList(u_p, i)
  {
    tVarList *vl = ListEntry(u_p,i);

    /* set limiter data in indicators (indc) if trouble indicator needs it */
    if(ListEntry(evosys->f[TROUBLE],i))
    {
      /* set data RDMP indicator needs in myindc arrays of each node */
      formylnodes(mesh)
      {
        tNode *node = MyLnode;
        limdata_MRS(node, vl); /* sets min,max,average */
      }

      /* initiate indc exchange */
      request_all_myln_indc_exchange_for_vl(mesh, vl);
      /* After this we could do work. MPI is now busy sending buffers */

      /* now get the indicators and wait for MPI buffers if necessary */
      get_all_myln_indc_for_vl(mesh, vl);
    }
  } /* end forList */
}
