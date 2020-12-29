/* refine.c */
/* Wolfgang Tichy, 5/2019 */


#include "nmesh.h"
#include "amr.h"

#define PR 0



/* Set n and pt_typ depending on ref->method. pnode is the parent node
   if we do h-refinement, and the current node if we do a p-refinement. */
void hp_refine_set_n_pt_typ(tNode *pnode, tRef *ref, int *n, int *pt_typ)
{
  int d;

  /* default point type for refined node(s) is taken from pnode */
  for(d=0; d<3; d++) pt_typ[d] = pnode->pt_typ[d];

  /* pick n */
  switch(ref->method)
  {
  case PARENT_nO2:
    for(d=0; d<3; d++)
    {
      n[d] = pnode->n[d]/2;
      if(n[d]<1) n[d] = 1;  /* do not allow n[d]<1 */
    }
    break;
  case PARENT_nO2_P1:
    for(d=0; d<3; d++)
    {
      n[d] = pnode->n[d]/2 + 1;
    }
    break;
  case PARENT_nO2_P1IFnG3:
    for(d=0; d<3; d++)
    {
      int pn = pnode->n[d];
      if(pn>3) n[d] = pn/2 + 1; /* add 1, unless pnode->n <= 3 */
      else     n[d] = pn/2;
      if(n[d]<1) n[d] = 1;  /* do not allow n[d]<1 */
    }
    break;
  case PARENT_nO2_P1MOD:
    for(d=0; d<3; d++)
    {
      int pn = pnode->n[d];
      if(pn>3) n[d] = pn/2 + 1;
      else     n[d] = pn - 1;
      if(n[d]<1) n[d] = 1;  /* do not allow n[d]<1 */
    }
    break;
  case GIVEN_n:
    for(d=0; d<3; d++) n[d] = ref->n[d];
    break;

  case PARENT_n_P_LGL:
    for(d=0; d<3; d++)
    {
      n[d] = pnode->n[d];
      pt_typ[d] = P_LGL;
    }
    break;
  case PARENT_n_P_UNIFORM:
    for(d=0; d<3; d++)
    {
      n[d] = pnode->n[d];
      pt_typ[d] = P_UNIFORM;
    }
    break;

  case GIVEN_n_P_LGL:
    for(d=0; d<3; d++)
    {
      n[d] = ref->n[d];
      pt_typ[d] = P_LGL;
    }
    break;
  case GIVEN_n_P_UNIFORM:
    for(d=0; d<3; d++)
    {
      n[d] = ref->n[d];
      pt_typ[d] = P_UNIFORM;
    }
    break;

  case PARENT_n:
  default:
    for(d=0; d<3; d++) n[d] = pnode->n[d];
  }
}


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
    int *n, *pt_typ;

    /* forward to node with nid[i] */
    for(; elem->node->nid != nid[i]; elem = elem->next) ;

    /* make children */
    parent = elem->node;
    pt_typ = elem->node->pt_typ; /* pick pt_typ */
    n = elem->node->n; /* pick n */
    children = make8_child_nodes(parent, pt_typ, n);

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
                                   tRef *ref)
{
  tNlist **replace, **children;
  long i;

  if(nnodes<=0) return;

  /* get mem. */
  replace  = calloc(nnodes, sizeof(replace[0]));
  children = calloc(nnodes, sizeof(children[0]));
  if(!replace || !children) errorexit("no memory for replace, children");

  NODELEVEL_Pragma(omp parallel)
  {
    tNlist *elem = mesh->lns;

    NODELEVEL_Pragma(omp for)
    for(i=0; i<nnodes; i++)
    {
      tNode *parent;
      int pt_typ[3], n[3];

      /* forward to node with nid[i] */
      //for(; elem && elem->node->nid != nid[i]; elem = elem->next) ;
      //if(!elem) errorexiti("could not find nid[i]=%d", nid[i]);
      for(; elem->node->nid != nid[i]; elem = elem->next) ;

      /* find parent */
      parent = elem->node;

      /* set n and pt_typ */
      hp_refine_set_n_pt_typ(parent, ref, n, pt_typ);

      /* make children */
      children[i] = make8_child_nodes(parent, pt_typ, n);

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

/* hrefine nids in the order in which we recv the MPI messages.
   NOTE: The order in which we recv is not certain => we do not have the same
   order for all MPI procs. => order of node->nfaces and thus node->fnb
   differs between different MPI procs!!! => nb index ni differs between
   MPI procs => request_surfaces_exchange_for_all_vars deadlocks because
   tags in there need a unique ni!!! */
void hrefine_nids_in_recv_order(tMesh *mesh, nMPI_Req *req,
                                int *nn, long **ref_nid,
                                int todo, tRef *ref)
{
  int rank = nMPI_rank(); /* my own rank */
  int size = nMPI_size();
  int r, done, flag;

  /* refine my own ref_nid[rank] */
  if(nn[rank]>0)
    create_children_no_nid_update(mesh, nn[rank], ref_nid[rank], ref);

  /* check for incoming broadcasts and then work on them */
  r = 0;
  done = 0;
  while(done<todo)
  {
    if(nn[r]>0)
    {
      nMPI_Test(&(req[r]), &flag, nMPI_STATUS_IGNORE);
      if(flag)
      {
        /* work on ref_nid[r] */
        if(r != rank) /* r=rank has been done already above */
        {
          if(ref->type == H_REFINE)
          {
            /* do h-refinement */
            create_children_no_nid_update(mesh, nn[r], ref_nid[r], ref);
          }
          else
          {
            /* do p-refinement */
            errorexit("implement p-refinement");
          }
        }
        nn[r] = 0;
        done++;
      }
    }
    r++;
    if(r>=size) r = 0;
  }
}

/* hrefine nids in the order of the MPI ranks. */
void hrefine_nids_in_rank_order(tMesh *mesh, nMPI_Req *req,
                                int *nn, long **ref_nid,
                                int todo, tRef *ref)
{
  int size = nMPI_size();
  int r;

  /* wait for incoming broadcasts and then work on them */
  for(r=0; r<size; r++)
  {
    if(nn[r]>0)
    {
      nMPI_Wait(&(req[r]), nMPI_STATUS_IGNORE);

      /* work on ref_nid[r] */
      if(ref->type == H_REFINE)
      {
        /* do h-refinement */
        create_children_no_nid_update(mesh, nn[r], ref_nid[r], ref);
        nn[r] = 0;
      }
      else
      {
        /* do p-refinement */
        errorexit("implement p-refinement");
      }
    }
  } /* end for loop */
}

/* h- or p-refine all nodes on all MPI procs if indicated by node->rflag.
   If hrefine=1 we actually create child nodes,
   if hrefine=0 we just change the number (and possibly the spacing)
   of points */
void hp_refine_nodes_if_rflag(tMesh *mesh, tRef *ref)
{
  int rank = nMPI_rank();
  int size = nMPI_size();
  int r, todo;
  int nnodes = (mesh->myln->nncats)*(mesh->myln->nm);
  long *my_nid   = calloc(nnodes+1, sizeof(my_nid[0]));
  nMPI_Req *req  = calloc(size, sizeof(req[0]));
  int *nn        = calloc(size, sizeof(nn[0]));
  long **ref_nid = calloc(size, sizeof(ref_nid[0]));

  if(!my_nid || !req || !nn || !ref_nid)
    errorexit("no memory for my_nid, req, nn, ref_nid");

  ///* record which nodes we want to refine */
  //formylnodes(mesh)
  //{
  //  tNode *node = MyLnode;
  //  node->rflag = needs_refine(node);
  //}

  /* save all nids where we need refinement in my_nid */
  nnodes = 0;
  formylnodes_noomp(mesh)
  {
    tNode *node = MyLnode;
    if(node->rflag > 0)
      my_nid[nnodes++] = node->nid;
  }
  nn[rank] = nnodes;

  /* number of procs where we have to refine and how many we have done so far */
  todo = size;

  /* broadcast number of nodes nn to all MPI jobs */
  for(r=0; r<size; r++)
  {
    nMPI_Bcast(&(nn[r]), 1, nMPI_INT, r);

    if(nn[r]<=0) todo--; /* there is nothing to do if nn[r]=0 */

    if(r==rank)
    {
      ref_nid[r] = my_nid;
    }
    else
    {
      if(nn[r]>0)
        ref_nid[r] = calloc(nn[r], sizeof(ref_nid[r][0]));
      else
        ref_nid[r] = NULL;
    }
  }

  /* broadcast ref_nid to all MPI jobs */
  for(r=0; r<size; r++)
  {
    if(nn[r]>0)
      nMPI_Ibcast(&(ref_nid[r][0]), nn[r], nMPI_LONG, r, &(req[r]));
  }

  /* check for incoming broadcasts and then work on them */
  if(ref->type == H_REFINE)
  {
    /* do h-refinement */
    //hrefine_nids_in_recv_order(mesh, req, nn, ref_nid, todo, ref);
    hrefine_nids_in_rank_order(mesh, req, nn, ref_nid, todo, ref);
  }
  else
  {
    /* do p-refinement */
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

/* h-refine all nodes on all MPI procs if indicated by node->rflag */
void hrefine_nodes_if_rflag(tMesh *mesh, tRef *ref)
{
  ref->type = H_REFINE;
  hp_refine_nodes_if_rflag(mesh, ref);
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

  //PRF;printnodelist(mesh->lns);
  //prlarray("nid0", 2*nnodes, nid0);

  /* update mesh->lns by removing all in nid0 and their siblings */
  n_remain = 0;
  for(i=0; i<nnodes; i++)
  {
    tNlist *elem = mesh->lns;
    int update_lns;

    /* forward to node with nid0[2*i] */
    //for(; elem && elem->node->nid != nid0[2*i]; elem = elem->next) ;
    //if(!elem) errorexiti("could not find nid0[2*i]=%d", nid0[2*i]);
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

    /* update mesh->lns if needed */
    if(elem == mesh->lns)
      update_lns = 1;
    else
      update_lns = 0;

    //printf("B:rm %ld: ", i);
    //printnodelistelement_and_neighbors_flag(elem,2);
    ////printnodelist(elem);

    /* remove 8 siblings */
    elem_parent[i] = remove8siblings_in_mesh_lns(elem);
    if(update_lns) mesh->lns = elem_parent[i];

    /* set elem to parent */
    elem = elem_parent[i];
    //printf("E:rm %ld: ", i);
    //printnodelistelement_and_neighbors_flag(elem,2);
    ////printnodelist(elem);
  }

  /* Now mesh->lns is up to date, next destroy the children */
  NODELEVEL_Pragma(omp parallel for)
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

  /* do nothing if nid0b is empty */
  if(nb<=0) return n;

  /* if nid0 is not empty, merge nid0b into nid0 */
  if(n>0)
  {
    for(j=0; j<nb; j++)
    {
      long nid0b_2j = nid0b[2*j];

      /* find nid0b_2j in nid0 */
      for(i=0; i<n; i++)
        if(nid0[2*i] == nid0b_2j) break;

      if(i<n)  /* found nid0b_2j in nid0 */
        nid0[2*i+1] += nid0b[2*j+1];
      else     /* append nid0b to nid0 */
        { nid0[2*nn] = nid0b[2*j];  nid0[2*nn+1] = nid0b[2*j+1];  nn++; }
    }
  }
  else /* copy nid0b into nid0 */
  {
    for(j=0; j<2*nb; j++) nid0[j] = nid0b[j];
    nn = nb;
  }

  /* return new length of nid0 */
  return nn;
}

/* compare nid0 numbers for qsort */
int nid0_compar(const void *x1, const void *x2)
{
  const long *n1, *n2;
  n1 = x1;
  n2 = x2;

  if(*n1 < *n2) return -1;
  if(*n1 > *n2) return 1;
  return 0;
}

/* destroy nids in the order in which we recv the MPI messages.
   NOTE: The order in which we recv is not certain => we do not have the same
   order for all MPI procs. => order of node->nfaces and thus node->fnb
   differs between different MPI procs!!! => nb index ni differs between
   MPI procs => request_surfaces_exchange_for_all_vars deadlocks because
   tags in there need a unique ni!!! */
void destroy_nids_in_recv_order(tMesh *mesh, nMPI_Req *req,
                                int *nn, long **unref,
                                int todo, tRef *ref)
{
  int rank = nMPI_rank(); /* my own rank */
  int size = nMPI_size();
  int r, done, flag, nn_rank;

  /* unrefine my own unref[rank] */
  nn_rank = nn[rank];
  if(nn[rank]>0)
    nn_rank = destroy_nodes_no_nid_update(mesh, nn[rank], unref[rank]);
  //PRFs(" 2: ");prlarray("unref[rank]", 2*nn[rank], unref[rank]);

  /* check for incoming broadcasts and then work on them */
  r = 0;
  done = 0;
  while(done<todo)
  {
    if(nn[r]>0)
    {
      nMPI_Test(&(req[r]), &flag, nMPI_STATUS_IGNORE);
      if(flag)
      {
        /* work on unref[r] */
        if(r != rank) /* r=rank has been done already above */
          nn[r] = destroy_nodes_no_nid_update(mesh, nn[r], unref[r]);
        else
          nn[r] = nn_rank;
        nn[r] = -nn[r]; /* make nn[r] negative */
        done++;
      }
    }
    r++;
    if(r>=size) r = 0;
  }

  /* flip sign on nn[r] since we made it negative above */
  for(r=0; r<size; r++) nn[r] = -nn[r];
  //PRFs(" 3: ");prlarray("unref[rank]", 2*nn[rank], unref[rank]);
}

/* destroy nids in the order of the MPI ranks. */
void destroy_nids_in_rank_order(tMesh *mesh, nMPI_Req *req,
                                int *nn, long **unref,
                                int todo, tRef *ref)
{
  int size = nMPI_size();
  int r;

  /* check for incoming broadcasts and then work on them */
  for(r=0; r<size; r++)
  {
    if(nn[r]>0)
    {
      nMPI_Wait(&(req[r]), nMPI_STATUS_IGNORE);

      /* work on unref[r] */
      nn[r] = destroy_nodes_no_nid_update(mesh, nn[r], unref[r]);
    }
  }
}

/* Unrefine all nodes on all MPI procs if indicated by node->rflag */
void remove_nodes_if_rflag(tMesh *mesh, tRef *ref)
{
  int rank = nMPI_rank();
  int size = nMPI_size();
  int r, todo;
  int nnodes;
  int myn = (mesh->myln->nncats)*(mesh->myln->nm);
  long *my_unr  = calloc(2*myn+2, sizeof(my_unr[0]));
  nMPI_Req *req = calloc(size, sizeof(req[0]));
  int *nn       = calloc(size, sizeof(nn[0]));
  long **unref  = calloc(size, sizeof(unref[0]));
  tNode *sib0;

  if(!my_unr || !req || !nn || !unref)
    errorexit("no memory for my_unr, req, nn, unref");

  //PRF;printnodelist(mesh->lns);
  //PRF;printmesh(mesh);
  //printNlistarray(mesh->myln->nm, mesh->myln->ln[0]);

  ///* record which nodes we want to remove */
  //formylnodes(mesh)
  //{
  //  tNode *node = MyLnode;
  //  node->rflag = unrefine(node);
  //}

  /* save all sibling0 nids where we need refinement in my_unr */
  sib0 = NULL;
  nnodes = 0;
  formylnodes_noomp(mesh)
  {
    tNode *node = MyLnode;
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

    if(sib0->l==0)
    {
      if(sib0->rflag < 0) errorexit("root node cannot be removed!");
      else continue;
    }

    /* other siblings */
    sib[1] = sib[0]->nb[1];
    sib[2] = sib[0]->nb[3];
    sib[4] = sib[0]->nb[5];
    sib[3] = sib[2]->nb[1];
    sib[5] = sib[4]->nb[1];
    sib[6] = sib[4]->nb[3];
    sib[7] = sib[3]->nb[5];

    /* check if all that we have on this proc needs to be unrefined */
    uref = 0;
    for(ijk=0; ijk<8; ijk++)
    {
      if(sib[ijk]->child[0]) /* cannot unrefine if there are any children */
        goto continue_with_next_node;

      if(sib[ijk]->dat)
      {
        if(sib[ijk]->rflag < 0) uref++;
        else                    goto continue_with_next_node;
      }
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

  /* number of procs where we have to uref and how many we have done so far */
  todo = size;

  /* broadcast number of nodes nn to all MPI jobs */
  for(r=0; r<size; r++)
  {
    nMPI_Bcast(&(nn[r]), 1, nMPI_INT, r);

    if(nn[r]<=0) todo--; /* there is nothing to do if nn[r]=0 */

    if(r==rank)
    {
      unref[r] = my_unr;
    }
    else
    {
      if(nn[r]>0)
        unref[r] = calloc(2*nn[r], sizeof(unref[r][0]));
      else
        unref[r] = NULL;
    }
  }
  my_unr = NULL; /* we do not need my_unr anymore */

  /* broadcast unref to all MPI jobs */
  for(r=0; r<size; r++)
  {
    if(nn[r]>0)
      nMPI_Ibcast(&(unref[r][0]), 2*nn[r], nMPI_LONG, r, &(req[r]));
  }

  //PRFs(" 1: ");prlarray("unref[rank]", 2*nn[rank], unref[rank]);

  /* check for incoming broadcasts and then work on them */
  //destroy_nids_in_recv_order(mesh, req, nn, unref, todo, ref);
  destroy_nids_in_rank_order(mesh, req, nn, unref, todo, ref);

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
  //PRFs(" 4: ");prlarray("unref[rank]", 2*nn[rank], unref[rank]);

  /* finally remove the stuff in unref[rank] */
  nn[rank] = destroy_nodes_no_nid_update(mesh, nn[rank], unref[rank]);
  if(nn[rank])
  {
    //printmesh(mesh);
    prlarray("unref[rank]", 2*nn[rank]+4, unref[rank]);
    errorexiti("nn[rank]=%d is supposed to be 0 now!", nn[rank]);
  }

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


/***************************************************************************/
/* functions we can call to refine in some particular way */
/***************************************************************************/

/* Refine all nodes that have neighbors whose level is greater by an
   amount dl than that of each node. This makes only sense if dl>=1 */
void hrefine_nodes_if_nb_finer_by_dl(tMesh *mesh, int dl, tRef *ref)
{
  /* go over mesh */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    int f, ni;

    /* mark nodes as to be refined if nb is more refined */
    for(f=0; f<6; f++)
      for(ni=0; ni<node->nfnb[f]; ni++)
      {
        tNode *nb = node->fnb[f][ni];
        if(nb->l - node->l >= dl) node->rflag = ref->method;
      }
  }
  hrefine_nodes_if_rflag(mesh, ref);
}

/* refine all nodes that have finer neighbors */
void hrefine_nodes_if_nb_finer(tMesh *mesh, tRef *ref)
{
  hrefine_nodes_if_nb_finer_by_dl(mesh, 1, ref);
}


/* refine all nodes up to level l */
//void refine_mesh_to_level(tMesh *mesh, int l)
void hrefine_mesh_to_level__old(tMesh *mesh, int l)
{
  tNlist *el;

  for(el=mesh->lns; el; el = el->next)
  {
    while(el->node->l < l)
      el = make8children_in_mesh_lns_myln(el, el->node->pt_typ, el->node->n);
  }
}

/* refine patch number p in mesh */
//void refine_pat(tMesh *mesh, int p)
void hrefine_pat__old(tMesh *mesh, int p)
{
  tPat *pat = mesh->pat[p];
  tNlist *el, *en;

  el = mesh->lns;
  for(en = el->next; el; en = el ? el->next : 0)
  {
    if(el->node->pat == pat)
      make8children_in_mesh_lns_myln(el, el->node->pt_typ, el->node->n);
    el = en;
  }
}



/* refine all nodes up to level l */
void hrefine_mesh_to_level(tMesh *mesh, int l)
{
  int i, ref;
  tRef rf[1];
  rf->method = PARENT_n;

  for(i=0; i<l; i++)
  {
    tNlist *el;
    ref = 0;
    fornodelist(mesh->lns, el)
    {
      tNode *node = el->node;
      if(node->l < l)
      {
        node->rflag = rf->method; /* flag node for refinement */
        ref++;                    /* count number of nodes that need refinement */
      }
      else
      {
        node->rflag = 0;
      }
    }

    if(ref)
    {
      hrefine_nodes_if_rflag(mesh, rf);
      update_mesh_myln_node_nid(mesh);
      if(PR)
      {
        PRF;printf(": On rank%d mesh is now:\n", nMPI_rank());
        printmesh(mesh);
      }
    }
    else
    {
      break;
    }
  }
}

/* like hrefine_mesh_to_level, but load balance after each level is created */
void hrefine_mesh_to_level_loadbalance(tMesh *mesh, int l)
{
  int i;
  for(i=0; i<=l; i++)
  {
    hrefine_mesh_to_level(mesh, i);
    simple_load_balance(mesh);
  }
}

/* coarsen all nodes up to level l */
void hcoarsen_mesh_to_level(tMesh *mesh, int l)
{
  int ref;
  tRef rf[1];
  rf->method = PARENT_n;

  do
  {
    tNlist *el;
    ref = 0;
    fornodelist(mesh->lns, el)
    {
      tNode *node = el->node;
      if(node->l > l)
      {
        node->rflag = -rf->method; /* flag node for unrefinement */
        ref++;                     /* count number of nodes that need refinement */
      }
      else
      {
        node->rflag = 0;
      }
    }

    if(ref)
    {
      remove_nodes_if_rflag(mesh, rf);
      update_mesh_myln_node_nid(mesh);
      if(PR)
      {
        PRF;printf(": On rank%d mesh is now:\n", nMPI_rank());
        printmesh(mesh);
      }
    }
  } while(ref);
}

/* refine patch number p in mesh */
void hrefine_pat(tMesh *mesh, int p)
{
  tPat *pat = mesh->pat[p];
  tNlist *el;
  tRef rf[1];
  rf->method = PARENT_n;

  fornodelist(mesh->lns, el)
  {
    tNode *node = el->node;
    if(node->pat == pat) node->rflag = rf->method;
    else                 node->rflag = 0;
  }
  hrefine_nodes_if_rflag(mesh, rf);
  update_mesh_myln_node_nid(mesh);
  if(PR)
  {
    PRF;printf(": On rank%d mesh is now:\n", nMPI_rank());
    printmesh(mesh);
  }
}

/* unrefine patch number p in mesh */
void hcoarsen_pat(tMesh *mesh, int p)
{
  tPat *pat = mesh->pat[p];
  tNlist *el;
  tRef rf[1];
  rf->method = PARENT_n;

  fornodelist(mesh->lns, el)
  {
    tNode *node = el->node;
    if(node->pat == pat) node->rflag = -rf->method;
    else                 node->rflag = 0;
  }
  remove_nodes_if_rflag(mesh, rf);
  update_mesh_myln_node_nid(mesh);
  if(PR)
  {
    PRF;printf(": On rank%d mesh is now:\n", nMPI_rank());
    printmesh(mesh);
  }
}


/* refine to better resolve nodes where limiting occured */
void hrefine_pcoarsen_nodes_if_nlim(tMesh *mesh)
{
  tRef ref[1];
  ref->method = PARENT_nO2_P1;

  /* look for nodes where limiting occured and refine them */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    /* flag refinement if we have more than 1 point */
    if(node->dat->info->nlim && node->np > 1)
      node->rflag = ref->method;
    else
      node->rflag = 0;

    //if(node->rflag)
    //{
    //  char s[100];
    //  printf("%s: nlim=%d rflag=%d\n", nodename(node,s,99),
    //         node->dat->info->nlim, node->rflag);
    //}
  }
  hrefine_nodes_if_rflag(mesh, ref);
  update_mesh_myln_node_nid(mesh);

  if(PR)
  {
    PRF;printf(": On rank%d mesh is now:\n", nMPI_rank());
    printmesh(mesh);
  }
  //exit(9);
}


/* remove previously refined nodes if no limiting occured */
void undo_hrefine_pcoarsen_nodes_if_zero_nlim(tMesh *mesh)
{
  tRef ref[1];
  ref->method = PARENT_nO2_P1;

  /* look for nodes where no limiting occured */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    tNode *parent = node->parent;
    int p_np, p_rflag;

    if(parent)
    {
      p_np = parent->np;
      p_rflag = parent->rflag;
    }
    else
    {
      p_np = p_rflag = 0;
    }

    /* flag h-unrefinement if parent has more points */
    if( (node->dat->info->nlim==0) && (node->np < p_np) && (p_rflag == ref->method) )
      node->rflag = -ref->method;
    else
      node->rflag = 0;

    //if(node->rflag)
    //{
    //  char s[100];
    //  printf("%s: nlim=%d rflag=%d\n", nodename(node,s,99),
    //         node->dat->info->nlim, node->rflag);
    //}
  }
  remove_nodes_if_rflag(mesh, ref);
  update_mesh_myln_node_nid(mesh);

  if(PR)
  {
    PRF;printf(": On rank%d mesh is now:\n", nMPI_rank());
    printmesh(mesh);
  }
}

/* use limiter data in node->dat->info->nlim to decide if and where we refine */
int resolve_shocks_using_nlim(tMesh *mesh)
{
  undo_hrefine_pcoarsen_nodes_if_zero_nlim(mesh);
  hrefine_pcoarsen_nodes_if_nlim(mesh);
  return 0;
}


/************************************************************************/
/* h-refine in nested spherical regions */
/************************************************************************/

/* refine within a sphere centered on (xc[0],xc[1],xc[2]) */
int hrefine_once_within_sphere(tMesh *mesh, double radius, double xc[3],
                               tRef *ref)
{
  int cnt = 0;

  /* go over mesh */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    double nc[3], r;

    /* find position of node center nc */
    set_nodecenter_xyz(node, nc);

    /* distance from xc to nc */
    r = Cart_distance_x0_x1(node, nc, xc);

    //pr3v("xc",xc);
    //pr3v("nc",nc);
    //printf(" radius=%g r=%g\n", radius, r);

    if(r<=radius)
    {
      node->rflag = ref->method;
      cnt++;
    }
  }
  hrefine_nodes_if_rflag(mesh, ref);
  update_mesh_myln_node_nid(mesh);
  return cnt;
}

/************************************************************************/
/* h-refine concentric nested spherical regions
         _______
      __/       \__
     /    _____    \     radius : radius of innermost sphere
    /    /     \    \    xc     : center of spheres
   /    /   _   \    \   levels : number of refinement levels we add
  |    |   |_|   |    |
   \    \       /    /   this should be controled by pars like:
    \    \_____/    /    amr_refine_sphere_radius : radius of innermost sphere
     \__         __/     amr_refine_sphere_levels : how many spheres we add
        \_______/
 */
void hrefine_sphere(tMesh *mesh, double radius, double xc[3], int levels)
{
  tRef ref[1];
  ref->method = PARENT_n;
  int i;
  double r;

  r = radius * pow(2., levels-1);
  for(i=0; i<levels; i++)
  {
    int cnt = hrefine_once_within_sphere(mesh, r, xc, ref);
    PRF;printf(": refined %d nodes within r=%g\n", cnt, r);
    r = r*0.5;
  }
}

/* use hrefine_sphere together with simple_load_balance */
void hrefine_sphere_loadbalance(tMesh *mesh, double radius, double xc[3],
                                int levels)
{
  int i;
  double r;

  r = radius * pow(2., levels-1);
  for(i=0; i<levels; i++)
  {
    hrefine_sphere(mesh, r, xc, 1);
    simple_load_balance(mesh);
    r = r*0.5;
  }
}
