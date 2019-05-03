/* refine.c */
/* Wolfgang Tichy, 5/2019 */


#include "nmesh.h"
#include "amr.h"

#define PR 0



/* refine nodes with nids in array, we assume long *nid is sorted in ascending
   order. We do not update nid's in here */
void refine_nodes_without_nid_update__old(tMesh *mesh, long nnodes, long *nid)
{
  tNlist *elem = mesh->lns;
  long i;

  for(i=0; i<nnodes; i++)
  {
    tNode *parent;
    tNlist *children;
    tNlist *lastchild;
    int *n;

    /* forward to node with nid[i] */
    for(; elem->node->nid != nid[i]; elem = elem->next) ;

    /* make children */
    parent = elem->node;
    n = elem->node->n; /* pick n */
    children = make8_child_nodes(parent, n);

    /* update mesh->lns if needed and add children to list */
    if(elem == mesh->lns) mesh->lns = first_nodelist(children);
    lastchild = replace1_in_nodelist(elem, children, 1);

    /* set elem to last child we added */
    elem = lastchild;
  }
}


/* refine nodes with nids in array, we assume long *nid is sorted in ascending
   order. We do not update nid's in here */
void refine_nodes_without_nid_update(tMesh *mesh, long nnodes, long *nid,
                                     int ref_method)
{
  tNlist **replace  = calloc(nnodes, sizeof(replace[0]));
  tNlist **children = calloc(nnodes, sizeof(children[0]));
  long i;

  if(!replace || !children) errorexit("no memory for replace, children");

  FORNODES_Pragma(omp parallel)
  {
    tNlist *elem = mesh->lns;

    FORNODES_Pragma(omp for)
    for(i=0; i<nnodes; i++)
    {
      tNode *parent;
      int nc[3], *n, d;

      /* forward to node with nid[i] */
      for(; elem->node->nid != nid[i]; elem = elem->next) ;

      /* make children */
      parent = elem->node;

      /* pick n */
      switch(ref_method)
      {
      case PARENT_n_O2:
        for(d=0; d<3; d++)
        {
          nc[d] = parent->n[d] / 2;
          if(nc[d]<1) nc[d] = 1;  /* do not allow n[d]<1 */
          n = nc;
        }
        break;
      case PARENT_n:
      default:
        n = parent->n;
      }
      children[i] = make8_child_nodes(parent, n);

      /* save elem that has to be replaced by children[i] later */
      replace[i] = elem;
    }
  }

  /* update mesh->lns using info in replace[i], children[i] */
  for(i=0; i<nnodes; i++)
  {
    tNlist *elem = replace[i];

    /* update mesh->lns if needed and add children to list */
    if(elem == mesh->lns) mesh->lns = first_nodelist(children[i]);
    replace1_in_nodelist(elem, children[i], 1); // can return last child
  }
  free(children);
  free(replace);
}

/* Refine all nodes on all MPI procs if indicated by func needs_refine. */
void refine_all_nodes_if_needed(tMesh *mesh, int (*needs_refine)(tNode *n),
                                int ref_method)
{
  int rank = nMPI_rank();
  int size = nMPI_size();
  int r, myid, done, flag;
  int nnodes = (mesh->myln->nncats)*(mesh->myln->nm);
  long *my_nid   = calloc(nnodes, sizeof(my_nid[0]));
  nMPI_Req *req  = calloc(size, sizeof(req[0]));
  int *nn        = calloc(size, sizeof(nn[0]));
  long **ref_nid = calloc(size, sizeof(ref_nid[0]));

  if(!my_nid || !req || !nn || !ref_nid)
    errorexit("no memory for my_nid, req, nn, ref_nid");

  /* record which nodes we want to refine */
  formylnodes(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
    node->refine = needs_refine(node);
  }

  /* save all nids where we need refinement in my_nid */
  nnodes = 0;
  formylnodes_noomp(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
    if(node->refine)
      my_nid[nnodes++] = node->nid;
  }
  nn[rank] = nnodes;

  /* broadcast number of nodes nn to all MPI jobs */
  for(r=0; r<size; r++)
  {
    nMPI_Bcast(&(nn[r]), 1, nMPI_INT, r);
    if(r==rank)
      ref_nid[r] = my_nid;
    else
      ref_nid[r] = calloc(nn[r], sizeof(ref_nid[r][0]));
  }

  /* broadcast my_nid to all MPI jobs */
  for(r=0; r<size; r++)
  {
    nMPI_Ibcast(&(ref_nid[r]), nn[r], nMPI_LONG, r, &(req[r]));
  }

  /* refine my_nid */
  refine_nodes_without_nid_update(mesh, nnodes, my_nid, PARENT_n_O2);
  nn[rank] = 0;
  done = 1;

  /* check for incoming broadcasts and then work on them */
  r = 0;
  while(done<size)
  {
    nMPI_Test(&(req[r]), &flag, nMPI_STATUS_IGNORE);
    if(flag && nn[r]>0)
    {
      /* work on ref_nid[r] */
      refine_nodes_without_nid_update(mesh, nn[r], ref_nid[r], PARENT_n_O2);
      nn[r] = 0;
      done++;
    }
    r++;
    if(r>=size) r = 0;
  }

  /* free ref_nid content */
  for(r=0; r<size; r++)
  {
    free(ref_nid[r]);
    ref_nid[r] = NULL;
  }
  /* my_nid was freed as one of ref_nid */

  /* free rest */
  free(ref_nid);
  free(nn);
  free(req);
}


/* remove nodes with nids in array nid0 and their 7 other siblings,
   we assume long *nid0 is sorted in ascending order.
   We do not update nid's in here */
void unrefine_nodes_without_nid_update(tMesh *mesh, long nnodes, long *nid0,
                                       int ref_method)
{
  tNlist **replace  = calloc(nnodes, sizeof(replace[0]));
  tNlist **children = calloc(nnodes, sizeof(children[0]));
  long i;

  if(!replace || !children) errorexit("no memory for replace, children");

  /* update mesh->lns by removing all  in nid0 and their siblings */
  for(i=0; i<nnodes; i++)
  {
    tNlist *elem;

    /* forward to node with nid0[i] */
    for(; elem->node->nid != nid0[i]; elem = elem->next) ;

    /* update mesh->lns if needed and add children to list */
    if(elem == mesh->lns) mesh->lns = first_nodelist(children[i]);
    replace1_in_nodelist(elem, children[i], 1); // can return last child
  }


  FORNODES_Pragma(omp parallel)
  {
    tNlist *elem = mesh->lns;

    FORNODES_Pragma(omp for)
    for(i=0; i<nnodes; i++)
    {
      tNode *parent;
      int nc[3], *n, d;

      /* forward to node with nid0[i] */
      for(; elem->node->nid != nid0[i]; elem = elem->next) ;

      /* destroy 8 siblings */
      parent = elem->node->parent;
      destroy_children(parent);

      /* save elem that has to be replaced by children[i] later */
      replace[i] = elem;
    }
  }

  free(children);
  free(replace);
}
