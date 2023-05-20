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

  case PARENT_nO2_P_LGL:
    for(d=0; d<3; d++)
    {
      if(pnode->n[d]>1) n[d] = pnode->n[d]/2;
      else              n[d] = pnode->n[d];   /* do not allow n=0 */
      pt_typ[d] = P_LGL;
    }
    break;
  case PARENT_nO2_P_UNIFORM:
    for(d=0; d<3; d++)
    {
      if(pnode->n[d]>1) n[d] = pnode->n[d]/2;
      else              n[d] = pnode->n[d];   /* do not allow n=0 */
      pt_typ[d] = P_UNIFORM;
    }
    break;

  case PARENT_2n_P_LGL:
    for(d=0; d<3; d++)
    {
      if(pnode->n[d]>1) n[d] = pnode->n[d]*2;
      else              n[d] = pnode->n[d];   /* keep n=1 */
      pt_typ[d] = P_LGL;
    }
    break;
  case PARENT_2n_P_UNIFORM:
    for(d=0; d<3; d++)
    {
      if(pnode->n[d]>1) n[d] = pnode->n[d]*2;
      else              n[d] = pnode->n[d];   /* keep n=1 */
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


////////////////////////////////////////////////
//Rename these 2, but keep old names as defines:

/* h-refine all nodes on all MPI procs if indicated by node->rflag */
void hrefine_nodes_if_rflag(tMesh *mesh, tRef *ref)
{
  ref->type = H_REFINE;
  hp_refine_elms_if_rflag(mesh, ref);
}

/* p-refine all nodes on all MPI procs if indicated by node->rflag */
void prefine_nodes_if_rflag(tMesh *mesh, tRef *ref)
{
  ref->type = P_REFINE;
  hp_refine_elms_if_rflag(mesh, ref);
}

////////////////////////////////////////////////////


// ====
// NEW:
// ====
/* h- or p-refine all elms if indicated by elm->rflag.
   If we h-refine we actually create child elms,
   otherwise we just change the number (and possibly the spacing)
   of points */
void hp_refine_elms_if_rflag(tMesh *mesh, tRef *ref)
{
  struct list_head *pos, *sav;

  /* loop over list with elms */
  list_for_each_safe(pos, sav, &mesh->myelm_head)
  {
    tElm *elm = list_entry(pos, tElm, list);
    int pt_typ[3], n[3];

    /* refine only if rflag is set */
    if(elm->rflag > 0)
    {
      /* set n and pt_typ */
      hp_refine_set_n_pt_typ(elm, ref, n, pt_typ);

      if(ref->type == H_REFINE) /* replace elm with its children */
      {
        tElm *child0 = replace_parent_by_8children(elm, n, pt_typ);
        set_children_nbinfo_remove_parent(child0, elm);
      }
      else /* p-refine by changing n and pt_typ of elm */
      {
        update_node_n_pt_typ(elm, n, pt_typ);
      }
    }
  }
  /* FIXME: This does not update the list mesh->myelm!
            Only the linked list mesh->myelm_head is changed here */
}

/* set some nbinfo and then remove and free the old parent */
void set_children_nbinfo_remove_parent(tElm *child0, tElm *parent)
{
  /* #pragma omp critical (change_mesh_myelm_list) */
  /* NOTE: For some reason gcc's -fsanitize=thread throws a ?false? positive
           if I use a named critical section!
           So replace "GEN_Pragma(omp critical (change_mesh_myelm_list))"
           by "GEN_Pragma(omp critical)" when debugging races!!! */
  //GEN_Pragma(omp critical (change_mesh_myelm_list))
  GEN_Pragma(omp critical)
  {
    int nbs_on_other_rank;
    /* NOTE: The new children have all zero for nfnb, fnb, and nnbinfo=-1.
             Also, all their neighbors have now the wrong nfnb and nbinfo.
             Even worse, all its neighbors have fnb pointers pointing
             to the parent which will be removed!!!  */
    /* FIXME: Maybe go over parent's nbs and set whatever
              nb-info we can!!! */
    /* we just created the children (child0), now set local nb-info */
    //int amr_set_nbinfo_of_new_children(tElm *child0, tElm *parent)
    // use: connections_get_nbloc_InsidePat

    /* for now we just invalidate a lot and remove the parent */
    nbs_on_other_rank = amr_invalidate_nbinfo_of_all_nbs(parent);
    if(nbs_on_other_rank)
      amr_remove_mesh_nbelm(Elm_mesh(parent));

    /* free parent and all data on parent */
    free_elm(parent);
  }
}

/* Unrefine all nodes on all MPI procs if indicated by node->rflag */
void remove_nodes_if_rflag(tMesh *mesh, tRef *ref)
{
  int rank = nMPI_rank();
  ulong myeidlim = mesh->eidlim[rank];
  //struct list_head *pos;

  /* arrays with missing eids and elmheaders */
  ulong neids;
  ulong eidarr[8];
  tElm0 elm0[8];
  int num, uref;
  tElm0 *elmar[8];
  int i;

  num = uref = 0;
  formyelms(mesh)
  {
    tElm *elm  = MyElm;
    ulong myid = MyID;
    //ulong eid  = Elm_eid(elm);
    int ijk = elm_get_ijk(elm);

    if(ijk==0)
    {
      /* try to get 8 elms, return how many we have on this rank */
      num  = amr_get_8elms_at_myid(mesh, myid, elmar);

      /* check if the num elms in elmar are indeed siblings and count
         how many want to be refined */
      uref = 0;
      if(amr_elms_are_siblings(num, elmar))
        for(uref=0, i=0; i<num; i++)
          if(elmar[i]->rflag < 0) uref++;

      //FIXME: do we want this:
      ///* if not all want to be refined, erase their rflags */
      //if(uref<num)
      //  for(uref=0, i=0; i<num; i++)
      //    elmar[i]->rflag = 0;

      /* if not all want to be refined, continue with next elm */
      if(uref<num)
      {
        num = 0;  /* do not look beyond these */
        continue;
      }
    }
  }

  /* ask for all the eids that we need beyond my rank */
  neids = 8 - num;

  /* set the eids I need beyond what I have */
  for(i=0; i<neids; i++) eidarr[i] = myeidlim+i;

  /* get elmheader for all in eidarr */
  amr_get_elm0_for_eids(mesh, neids, eidarr, elm0);

  /* add elm0 from other ranks to elmar */
  for(i=0; i<neids; i++)
    elmar[num+i] = &(elm0[i]);

  uref = 0;
  if(amr_elms_are_siblings(8, elmar))
    for(i=0; i<8; i++)
      if(elmar[i]->rflag < 0) uref++;

  /* get the neids missing elms onto my rank */

  /* change myelm_head only if all 8 want to be refined */
  if(uref==8)
  {
    /* insert elm0 after current end of list */
    for(i=0; i<neids; i++)
    {
      tElm *elm = alloc_elm_init_pat(mesh, elm0[i].eploc->p); /* fresh elm */
      memcpy(elm, &(elm0[i]), sizeof(tElm0)); /* init elm from r_elms[i] */
      /* now add elm to the end of list in mesh */
      list_add_tail(&elm->list, &mesh->myelm_head);
    }
  }

  /////////////////////////////////////////////////////////////////////////
  printf("PROBLEM !!!\n"
         "How do I tell the other ranks to mark the elms that I want?\n"
         "NOTE: load_exchange_dat_after_moving_elms only works if they set"
         "their dat->info->desrank equal to my rank!!!");
  errorexit("I give up because this PROBLEM has no easy solution!");
  /////////////////////////////////////////////////////////////////////////

  /* free surfaces & indc since they will change now anyway */
  evolve_free_communication_structs(mesh);

  /* move dat to correct ranks now */
  load_exchange_dat_after_moving_elms(mesh);

  //alloc_and_set_mesh_myelm(mesh);
  //NOTE: update_mesh_myln_node_nid call causes an update of mesh->myelm

  //FIXME: adapt  update_mesh_myln_node_nid
  update_mesh_myln_node_nid(mesh);

  //FIXME: call function that set's up elm->fnb and such...
  //       maybe also update_mesh_myln_node_nid ???

  /* now that nodes are elsewhere re-init surfaces & indc */
  evolve_init_communication_structs(mesh);


  ////////////////////////////////////////////////////////////////////


  /* FIXME: This does not update the list mesh->myelm!
            Only the linked list mesh->myelm_head is changed here */

  errorexit("finish: void remove_elms_if_rflag(tMesh *mesh, tRef *ref)");
  //0th  get all 8 sibling children onto one MPI rank
  //     (we have most of this above)

  /* if not all 8 want to be refined we do nothing */
  if(uref<8) return;

  //1st Call:
  //tElm *replace_8localchildren_by_parent(tElm *child0, int n[3], int pt_typ[3],
  //                                       struct list_head *ch_head)

  //2nd Call:
  //void set_parent_nbinfo_remove_children(tElm *parent,
  //                                       struct list_head *ch_head)
}

/* set some nbinfo and then free the children in ch_head */
void set_parent_nbinfo_remove_children(tElm *parent,
                                       struct list_head *ch_head)
{
  /* #pragma omp critical (change_mesh_myelm_list) */
  /* NOTE: For some reason gcc's -fsanitize=thread throws a ?false? positive
           if I use a named critical section!
           So replace "GEN_Pragma(omp critical (change_mesh_myelm_list))"
           by "GEN_Pragma(omp critical)" when debugging races!!! */
  //GEN_Pragma(omp critical (change_mesh_myelm_list))
  GEN_Pragma(omp critical)
  {
    int nbs_on_other_rank;
    struct list_head *pos_ijk, *sav;
    /* NOTE: This new parent has all zero for nfnb, fnb, and nnbinfo=-1.
             Also, all its neighbors have now the wrong nfnb and nbinfo.
             Even worse, all its neighbors have fnb pointers pointing
             to the children which will be removed!!!  */
    /* FIXME: Maybe go over children's nbs and set whatever
              nb-info we can!!! */
    /* we just created the parent, now set local nb-info */
    //int amr_set_nbinfo_of_new_parent(tElm *parent, struct list_head *ch_head)
    // use: connections_get_nbloc_InsidePat

    /* for now we just invalidate a lot and remove the children */
    list_for_each(pos_ijk, ch_head)
    {
      tElm *ch_ijk = list_entry(pos_ijk, tElm, list);

      nbs_on_other_rank = amr_invalidate_nbinfo_of_all_nbs(ch_ijk);
      if(nbs_on_other_rank)
        amr_remove_mesh_nbelm(Elm_mesh(parent));
    }

    /* free children */
    list_for_each_safe(pos_ijk, sav, ch_head)
    {
      tElm *ch_ijk = list_entry(pos_ijk, tElm, list);
      list_del(&ch_ijk->list); //del from ch_head
      free_elm(ch_ijk);        //free mem of child ch_ijk
    }
  }
}


/***************************************************************************/
/* old stuff. Remove later: */
/***************************************************************************/

/* p-refine nodes with nids in array, we assume nid[] is sorted in ascending
   order. We do not update nids in here */
void prefine_nid_list(tMesh *mesh, long nnodes, long *nid, tRef *ref)
{
  long i;

  if(nnodes<=0) return;

  NODELEVEL_Pragma(omp parallel)
  {
    tNlist *elem = mesh->lns;

    NODELEVEL_Pragma(omp for)
    for(i=0; i<nnodes; i++)
    {
      tNode *node;
      int pt_typ[3], n[3];

      /* forward to node with nid[i] */
      //for(; elem && elem->node->nid != nid[i]; elem = elem->next) ;
      //if(!elem) errorexiti("could not find nid[i]=%d", nid[i]);
      for(; elem->node->eploc->eid != nid[i]; elem = elem->next) ;

      /* find node */
      node = elem->node;

      /* set n and pt_typ */
      hp_refine_set_n_pt_typ(node, ref, n, pt_typ);

      /* p-refine by changing n and pt_typ of node */
      update_node_n_pt_typ(node, n, pt_typ);
    }
  }
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



/***************************************************************************/
/* helper functions */
/***************************************************************************/

/* helper func to reset rflag in all nodes */
void refine_set_rflag_forall_nodes(tMesh *mesh, int rflag)
{
  /* go over mesh */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    node->rflag = rflag;
  }
}

/* helper func to set use_fv flag in all nodes that have rflag set */
void refine_set_use_fv_if_rflag(tMesh *mesh, int use_fv)
{
  /* go over mesh */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    if(node->rflag) node->dat->info->use_fv = use_fv;
  }
}

/* helper func to set use_fv flag in all nodes that have pt_typ */
void refine_set_use_fv_if_pt_typ(tMesh *mesh, int pt_typ[3], int use_fv)
{
  /* go over mesh */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    int d, equal;

    equal=1;
    for(d=0; d<3; d++)
      if(node->pt_typ[d] != pt_typ[d]) { equal=0; break; }
    if(equal) node->dat->info->use_fv = use_fv;
  }
}

/* Synchronize ref->method on all procs. This assumes that some procs have
   ref->method=REF_METH_DONOTHING, while others have one particular
   value that is higher. We use MPI_Allreduce to set them all to the max
   ref->method. */
int refine_synchronize_ref_method(tRef *ref)
{
  int Max_method = ref->method;

  nMPI_Allreduce(&(ref->method), &Max_method, 1, nMPI_INT, nMPI_MAX);
  ref->method = Max_method;
  return Max_method;
}


/***************************************************************************/
/* functions we can call to h-refine in some particular way */
/***************************************************************************/

/* h-refine all nodes that have neighbors whose level is greater by an
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
        if(Elm_l(nb) - Elm_l(node) >= dl) node->rflag = ref->method;
      }
  }
  hrefine_nodes_if_rflag(mesh, ref);
}

/* h-refine all nodes that have finer neighbors */
void hrefine_nodes_if_nb_finer(tMesh *mesh, tRef *ref)
{
  hrefine_nodes_if_nb_finer_by_dl(mesh, 1, ref);
}


/* h-refine all nodes up to level l */
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
      if(Elm_l(node) < l)
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
      if(Elm_l(node) > l)
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

/* h-refine patch number p in mesh */
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
     /    _____    \    radius : radius of innermost sphere
    /    /     \    \   xc     : center of spheres
   /    /   _   \    \  levels : number of refinement levels we add
  |    |   |_|   |    |
   \    \       /    /  this should be controled by pars like:
    \    \_____/    /   amr_hrefine_sphere_radius : radius of innermost sphere
     \__         __/    amr_hrefine_sphere_levels : how many spheres we add
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


/************************************************************************/
/* p-refine under certain conditions */
/************************************************************************/

/* p-refine all nodes that have neighbors that have uniform grid spacing
   in any direction */
void prefine_nodes_if_nb_uniform_in_any_dir(tMesh *mesh, tRef *ref)
{
  /* go over mesh */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    int f, ni;

    /* mark nodes as to be refined if nb has uniform grid spacing */
    for(f=0; f<6; f++)
      for(ni=0; ni<node->nfnb[f]; ni++)
      {
        tNode *nb = node->fnb[f][ni];

        if(nb->pt_typ[0]==P_UNIFORM ||
           nb->pt_typ[1]==P_UNIFORM ||
           nb->pt_typ[2]==P_UNIFORM)
          node->rflag = ref->method;
        /* FIXME: skip nodes that already have the same refinement as asked
                  for in ref */
      }
  }
  prefine_nodes_if_rflag(mesh, ref);
}

/* p-refine patch number p in mesh to have n[] points */
void prefine_pat(tMesh *mesh, int p, int n[3])
{
  tPat *pat = mesh->pat[p];
  tNlist *el;
  tRef rf[1];
  rf->method = GIVEN_n;
  rf->n[0] = n[0];
  rf->n[1] = n[1];
  rf->n[2] = n[2];

  fornodelist(mesh->lns, el)
  {
    tNode *node = el->node;
    if(node->pat == pat) node->rflag = rf->method;
    else                 node->rflag = 0;
  }
  prefine_nodes_if_rflag(mesh, rf);
  update_mesh_myln_node_nid(mesh);
  if(PR)
  {
    PRF;printf(": On rank%d mesh is now:\n", nMPI_rank());
    printmesh(mesh);
  }
}
