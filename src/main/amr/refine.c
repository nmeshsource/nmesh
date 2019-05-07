/* refine.c */
/* Wolfgang Tichy, 5/2019 */


#include "nmesh.h"
#include "amr.h"

#define PR 0



/* h-refine nodes with nids in array, we assume long *nid is sorted in ascending
   order. We do not update nid's in here */
void hrefine_nodes_without_nid_update__old(tMesh *mesh, long nnodes, long *nid)
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


/* h-refine nodes with nids in array, we assume nid[] is sorted in ascending
   order. We do not update nids in here */
void create_children_no_nid_update(tMesh *mesh, long nnodes, long *nid,
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
      for(; elem && elem->node->nid != nid[i]; elem = elem->next) ;
      //if(!elem) errorexiti("could not find nid[i]=%d", nid[i]);

      /* make children */
      parent = elem->node;

      /* pick n */
      switch(ref_method)
      {
      case PARENT_nO2:
        for(d=0; d<3; d++)
        {
          nc[d] = parent->n[d]/2;
          if(nc[d]<1) nc[d] = 1;  /* do not allow n[d]<1 */
          n = nc;
        }
        break;
      case PARENT_nO2_P1:
        for(d=0; d<3; d++)
        {
          nc[d] = parent->n[d]/2;
          if(nc[d]>1) nc[d] += 1; /* add 1, unless it is 1 already */
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

/* h-refine all nodes on all MPI procs if indicated by func needs_refine. */
void hrefine_nodes_if_rflag(tMesh *mesh, int ref_method)
{
  int rank = nMPI_rank();
  int size = nMPI_size();
  int r, myid, todo, done, flag;
  int nnodes = (mesh->myln->nncats)*(mesh->myln->nm);
  long *my_nid   = calloc(nnodes, sizeof(my_nid[0]));
  nMPI_Req *req  = calloc(size, sizeof(req[0]));
  int *nn        = calloc(size, sizeof(nn[0]));
  long **ref_nid = calloc(size, sizeof(ref_nid[0]));

  if(!my_nid || !req || !nn || !ref_nid)
    errorexit("no memory for my_nid, req, nn, ref_nid");

  ///* record which nodes we want to refine */
  //formylnodes(mesh, myid)
  //{
  //  tNode *node = MyNode(mesh, myid);
  //  node->rflag = needs_refine(node);
  //}

  /* save all nids where we need refinement in my_nid */
  nnodes = 0;
  formylnodes_noomp(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
    if(node->rflag)
      my_nid[nnodes++] = node->nid;
  }
  nn[rank] = nnodes;

  /* number of procs where we have to refine and how many we have done so far */
  todo = size;
  done = 0;

  /* broadcast number of nodes nn to all MPI jobs */
  for(r=0; r<size; r++)
  {
    nMPI_Bcast(&(nn[r]), 1, nMPI_INT, r);

    if(nn[r]<=0) todo--; /* there is nothing to do if nn[r]=0 */

    if(r==rank)
      ref_nid[r] = my_nid;
    else
      ref_nid[r] = calloc(nn[r], sizeof(ref_nid[r][0]));
  }

  /* broadcast ref_nid to all MPI jobs */
  for(r=0; r<size; r++)
  {
    if(nn[r]>0)
      nMPI_Ibcast(&(ref_nid[r][0]), nn[r], nMPI_LONG, r, &(req[r]));
  }

  /* refine my_nid */
  if(nnodes>0)
  {
    create_children_no_nid_update(mesh, nnodes, my_nid, ref_method);
    nn[rank] = 0;
    done++;
  }

  /* check for incoming broadcasts and then work on them */
  r = 0;
  while(done<todo)
  {
    if(nn[r]>0)
    {
      nMPI_Test(&(req[r]), &flag, nMPI_STATUS_IGNORE);
      if(flag)
      {
        /* work on ref_nid[r] */
        create_children_no_nid_update(mesh, nn[r], ref_nid[r], ref_method);
        nn[r] = 0;
        done++;
      }
    }
    r++;
    if(r>=size) r = 0;
  }

#if 0
r = rank^1;
create_children_no_nid_update(mesh, nn[r], ref_nid[r], ref_method);

//for(r=size-1; r>=0; r--)
////for(r=0; r<size; r++)
//create_children_no_nid_update(mesh, nn[r], ref_nid[r], ref_method);


r = rank;
//r=1;
create_children_no_nid_update(mesh, nn[r], ref_nid[r], ref_method);
r = rank^1;
//r=0;
create_children_no_nid_update(mesh, nn[r], ref_nid[r], ref_method);
#endif

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
   We assume long *nid0 is sorted in ascending order.
   nid0[2*i]   has the node id of sibling0, 
   nid0[2*i+1] has the number of siblings that need unrefinement
   We do not update node->nid in here, but we change nid0 and return
   what is left to process, i.e. all nid0 where we had less than 8 siblings
   in this MPI proc. */
long destroy_nodes_no_nid_update(tMesh *mesh, long nnodes, long *nid0)
{
  tNlist **elem_parent = calloc(nnodes, sizeof(elem_parent[0]));
  long i;
  long n_remain;

  if(!elem_parent) errorexit("no memory for elem_parent");

  /* update mesh->lns by removing all in nid0 and their siblings */
  n_remain = 0;
  for(i=0; i<nnodes; i++)
  {
    tNlist *elem;
    int update_lns;

    /* forward to node with nid0[2*i] */
    for(; elem->node->nid != nid0[2*i]; elem = elem->next) ;

    /* go to next node if not all 8 siblings need unrefinement,
       but save it for later in nid0 */
    if(nid0[2*i+1] < 8)
    {
      nid0[2*n_remain]   = nid0[2*i];
      nid0[2*n_remain+1] = nid0[2*i+1];
      n_remain++;
      continue;
    }

    /* update mesh->lns if needed and remove 8 siblings */
    if(elem == mesh->lns)
      update_lns = 1;
    else
      update_lns = 0;

    elem_parent[i] = remove8siblings_in_mesh_lns(elem);
    if(update_lns) mesh->lns = elem_parent[i];

    /* set elem to parent */
    elem = elem_parent[i];
  }

  /* Now mesh->lns is up to date, next destroy the children */
  FORNODES_Pragma(omp parallel for)
  for(i=0; i<nnodes; i++)
  {
    if(elem_parent[i])
    {
      tNode *parent = elem_parent[i]->node;
      /* destroy 8 children */
      destroy_children(parent);
    }
  }

  free(elem_parent);

  /* return number of nid0s not destroyed yet */
  return n_remain;
}

/* merge nid0-list in nid0b into nid0-list in nid0, the allocated size
   of nid0 has to be 2*(n+nb) */
long merge_nid0b_into_nid0(long n, long *nid0, long nb, long *nid0b)
{
  long i, j, nn=n;
  for(i=0; i<nn; i++)
  {
    long nid0_2i = nid0[2*i];

    /* find nid0_2i in nid0b */
    for(j=0; j<nb; j++)
      if(nid0b[2*j] == nid0_2i) break;

    if(j<nb) /* found nid0_2i in nid0b */
      nid0[2*i+1] += nid0b[2*j+1];
    else     /* append nid0b to nid0 */
      { nid0[2*nn] = nid0b[2*j];  nid0[2*nn+1] = nid0b[2*j+1];  nn++; }
  }
  /* return new length of nid0 */
  return nn;
}

/* compare nid0 numbers for qsort */
int nid0_compar(const void *x1, const void *x2)
{
  long *const *n1 = x1;
  long *const *n2 = x2;

  if(n1 < n2) return -1;
  if(n1 > n2) return 1;
  return 0;
}

/* Unrefine all nodes on all MPI procs if indicated by func needs_refine. */
void remove_nodes_if_rflag(tMesh *mesh, int ref_method)
{
  int rank = nMPI_rank();
  int size = nMPI_size();
  int r, myid, done, flag;
  int nnodes;
  int myn = (mesh->myln->nncats)*(mesh->myln->nm);
  long *my_unr  = calloc(2*myn, sizeof(my_unr[0]));
  nMPI_Req *req = calloc(size, sizeof(req[0]));
  int *nn       = calloc(size, sizeof(nn[0]));
  long **unref  = calloc(size, sizeof(unref[0]));
  tNode *sib0;

  if(!my_unr || !req || !nn || !unref)
    errorexit("no memory for my_unr, req, nn, unref");

  ///* record which nodes we want to remove */
  //formylnodes(mesh, myid)
  //{
  //  tNode *node = MyNode(mesh, myid);
  //  node->rflag = unrefine(node);
  //}

  /* save all sibling0 nids where we need refinement in my_unr */
  sib0 = NULL;
  nnodes = 0;
  formylnodes_noomp(mesh, myid)
  {
    tNode *node = MyNode(mesh, myid);
    tNode *sib[8];
    int ijk = node->ijk;
    int n[] = {2,2,2};
    int i,j,k, uref;
    long nid0;

    /* get sibling0 */
    sib[0] = node;
    k = kOfInd_n(ijk,n);
    j = jOfInd_n_k(ijk,n,k);
    i = iOfInd_n_jk(ijk,n,j,k);
    if(i) sib[0] = sib[0]->nb[0];
    if(j) sib[0] = sib[0]->nb[2];
    if(k) sib[0] = sib[0]->nb[4];
 
    /* check if we processed sib[0] already */
    if(sib[0] == sib0) continue;
    sib0 = sib[0]; /* save sib[0] as last one processed */

    /* other siblings */
    sib[1] = sib[0]->nb[1];
    sib[2] = sib[0]->nb[3];
    sib[4] = sib[0]->nb[5];
    sib[3] = sib[2]->nb[1];
    sib[5] = sib[4]->nb[1];
    sib[6] = sib[4]->nb[3];
    sib[7] = sib[3]->nb[5];

    /* check if all that we have on this proc needs to unrefined */
    uref = 0;
    for(ijk=0; ijk<8; ijk++)
    {
      if(sib[ijk]->dat && !(sib[ijk]->rflag))
        goto continue_with_next_node;
      uref++;
    }

    /* save nid0 and number of siblings that need a to be unrefined */
    nid0 = sib[0]->nid;
    my_unr[2*nnodes]   = nid0;
    my_unr[2*nnodes+1] = uref;
    nnodes++;

  continue_with_next_node:
    ;
  }
  nn[rank] = nnodes;

  /* broadcast number of nodes nn to all MPI jobs */
  for(r=0; r<size; r++)
  {
    nMPI_Bcast(&(nn[r]), 1, nMPI_INT, r);
    if(r==rank)
      unref[r] = my_unr;
    else
      unref[r] = calloc(2*nn[r], sizeof(unref[r][0]));
  }
  my_unr = NULL; /* we do not need my_unr anymore */

  /* broadcast unref to all MPI jobs */
  for(r=0; r<size; r++)
  {
    nMPI_Ibcast(&(unref[r][0]), 2*nn[r], nMPI_LONG, r, &(req[r]));
  }

  /* unrefine my own unref[rank] */
  nn[rank] = -destroy_nodes_no_nid_update(mesh, nn[rank], unref[rank]);
  done = 1;

  /* check for incoming broadcasts and then work on them */
  r = 0;
  while(done<size)
  {
    nMPI_Test(&(req[r]), &flag, nMPI_STATUS_IGNORE);
    if(flag && nn[r]>0)
    {
      /* work on unref[r] */
      nn[r] = -destroy_nodes_no_nid_update(mesh, nn[r], unref[r]);
      done++;
    }
    r++;
    if(r>=size) r = 0;
  }

  /* flip sign on nn[r] since we made it negative above */
  for(r=0; r<size; r++) nn[r] = -nn[r];


  /* merge the arrays with left over nids into one */
  for(r=0; r<size; r++)
  {
    if(r != rank)
    {
      int n;

      /* check if we need more room */
      n = nn[rank] + nn[r];
      if(n > myn)
      {
        unref[rank] = realloc(unref[rank], 2*n*sizeof(unref[rank][0]));
        myn = n;
      }

      /* now merge the arrays */
      n = merge_nid0b_into_nid0(nn[rank], unref[rank], nn[r], unref[r]);
      nn[rank] = n;
    }
  }

  /* sort new longer unref[rank] with left over nid0s */
  qsort(unref[rank], nn[rank], 2*sizeof(unref[rank][0]), nid0_compar);

  /* finally remove the stuff in unref[rank] */
  nn[rank] = destroy_nodes_no_nid_update(mesh, nn[rank], unref[rank]);
  if(nn[rank]) errorexit("nn[rank] is supposed to be 0 now!");

  /* free unref content */
  for(r=0; r<size; r++)
  {
    free(unref[r]);
    unref[r] = NULL;
  }
  /* old my_unr is freed as one of unref */

  /* free rest */
  free(unref);
  free(nn);
  free(req);
}



/* refine all nodes up to level l */
//void refine_mesh_to_level(tMesh *mesh, int l)
void refine_mesh_to_level__old(tMesh *mesh, int l)
{
  tNlist *el;

  for(el=mesh->lns; el; el = el->next)
  {
    while(el->node->l < l)
      el = make8children_in_mesh_lns_myln(el, el->node->n);
  }
}

/* refine patch number p in mesh */
//void refine_pat(tMesh *mesh, int p)
void refine_pat__old(tMesh *mesh, int p)
{
  tPat *pat = mesh->pat[p];
  tNlist *el, *en;

  el = mesh->lns;
  for(en = el->next; el; en = el ? el->next : 0)
  {
    if(el->node->pat == pat)
      make8children_in_mesh_lns_myln(el, el->node->n);
    el = en;
  }
}



/* refine all nodes up to level l */
void refine_mesh_to_level(tMesh *mesh, int l)
{
  int i, ref;

  for(i=0; i<l; i++)
  {
    tNlist *el;
    ref = 0;
    fornodelist(mesh->lns, el)
    {
      tNode *node = el->node;
      if(node->l < l)
      {
        node->rflag = 1; /* flag node for refinement */
        ref++;           /* count number of nodes that need refinement */
      }
      else
      {
        node->rflag = 0;
      }
    }

    if(ref)
    {
      hrefine_nodes_if_rflag(mesh, PARENT_n);
      update_mesh_myln_node_nid(mesh);
    }
    else
    {
      break;
    }
  }
}

/* refine patch number p in mesh */
void refine_pat(tMesh *mesh, int p)
{
  tPat *pat = mesh->pat[p];
  tNlist *el;

  fornodelist(mesh->lns, el)
  {
    tNode *node = el->node;
    if(node->pat == pat) node->rflag = 1;
    else                 node->rflag = 0;
  }
  hrefine_nodes_if_rflag(mesh, PARENT_n);
  update_mesh_myln_node_nid(mesh);
}
