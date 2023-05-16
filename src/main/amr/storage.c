/* storage.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "amr.h"

#define PR 0


/* world comm from main */
extern nMPI_Comm main_comm;

/* use gridpoints from basis/gridpoints.c */
extern tGridPoints gridpoints[1];

/* amr stuff */
extern tAMR amr[1];


/**************************************************************************/
/* basic memory management */
/**************************************************************************/

/**********************************************************************/
/* storage for arrays */
/**********************************************************************/
/* allocate an empty array with ns segments, and extra space of size Ne */
tArray *alloc_empty_array_with_segs(int n[3], int Ne, int ns)
{
  tArray *array;
  int i;

  if(!n) return NULL;

  array = calloc(1, sizeof(tArray));
  if(!array) errorexit("out of memory");

  array->N = n[0] * n[1] * n[2];
  for(i=0; i<3; i++)  array->n[i] = n[i];

  array->Ne = Ne;
  array->ns = ns;
  array->size = 0;

  return array;
}

/* allocate an array with ns segments, and extra space of size Ne */
tArray *alloc_array_with_segs(int n[3], int Ne, int ns)
{
  tArray *array;
  int Nt, nmemb;
  int size1 = max3(sizeof(array->d[0]), sizeof(array->i[0]),
                   sizeof(array->l[0]));
  if(!n) return NULL;

  array = alloc_empty_array_with_segs(n, Ne, ns);

  /* memory for data */
  Nt = array->N + array->Ne;
  nmemb = Nt * ns;
  if(nmemb>0)
  {
    array->d = calloc(nmemb, size1);
    if(!array->d) errorexit("out of memory for array->d");
  }
  /* else array->d = NULL; //because array was alloced with calloc */
  array->size = nmemb * size1;

  return array;
}

/* allocate a 1d array */
tArray *alloc_array1d_with_segs(int N, int Ne, int ns)
{
  int n[] = { N,1,1 };
  return alloc_array_with_segs(n, Ne, ns);
}

/* allocate an array with extra space but just 1 segment */
tArray *alloc_array_extra(int n[3], int Ne)
{
  return alloc_array_with_segs(n, Ne, 1);
}

/* allocate a standard array (without extra space and just 1 segment) */
tArray *alloc_array(int n[3])
{
  return alloc_array_with_segs(n, 0, 1);
}
/* allocate a 1d array */
tArray *alloc_array1d(int N)
{
  int n[] = { N,1,1 };
  return alloc_array(n);
}
/* allocate a 2d array */
tArray *alloc_array2d(int n0, int n1)
{
  int n[] = { n0,n1,1 };
  return alloc_array(n);
}
/* allocate an empty 2d array */
tArray *alloc_empty_array2d(int n0, int n1)
{
  int n[] = { n0,n1,1 };
  return alloc_empty_array_with_segs(n, 0, 1);
}

/* get array that starts at segment si */
tArray *get_array_seg(tArray *array, int si)
{
  if(si>0)
  {
    int Nt = array->N + array->Ne;
    tArray *as = calloc(1, sizeof(tArray));
    *as = *array; /* shallow copy */
    as->si = si;
    as->d  = array->d + Nt * si;
    return as;
  }
  else /* if si=0 we just return the array without allocating a new array */
    return array;
}

/* replace pointer to data with something else */
void point_array_d_to_data(tArray *array, void *data, int nofree)
{
  /* free data only if this is segm. 0, and only if it should be freed */
  if( (array->si == 0) && (array->d_nofree == 0) )
    free(array->d);
  array->d = data;
  array->d_nofree = nofree;
}

/* re-dimension array, return values:
   0 if no new memory was allocated,
   1 if new memory was allocated, but array->d pointer is unchanged
   2 if new memory was allocated and array->d pointer has changed */
int redimension_array_with_segs(tArray *array, int n[3], int Ne, int ns)
{
  int size1 = max3(sizeof(array->d[0]), sizeof(array->i[0]),
                   sizeof(array->l[0]));
  int Nt;
  size_t size_new;
  int reallocd_d = 0;
  int i;

  for(i=0; i<3; i++) if(n[i]>=0) array->n[i] = n[i];
  array->N = array->n[0] * array->n[1] * array->n[2];
  array->Ne = Ne;
  array->ns = ns;

  Nt = array->N + array->Ne;
  size_new = Nt * ns * size1;
  if(size_new > array->size)
  {
    double *d = realloc(array->d, size_new);
    if(!d) errorexit("out of memory for array->d");
    if(d != array->d) reallocd_d = 2;
    else              reallocd_d = 1;
    array->d = d;
    array->size = size_new;
  }
  return reallocd_d;
}
int redimension_array(tArray *array, int n[3])
{
  return redimension_array_with_segs(array, n, array->Ne, array->ns);
}
int redim_array(tArray *array, int n0, int n1, int n2)
{
  int n[] = { n0,n1,n2 };
  return redimension_array(array, n);
}

/* free an array */
void free_array(tArray *array)
{
  if(!array) return;
  if( (array->si == 0) && (array->d_nofree == 0) )
    free(array->d); /* free data only if this segm. 0 and it should be freed */
  free(array);
}

/* free 3 arrays that might point to same mem */
void free_3_arrays(tArray *array[3])
{
  int free1, free2;

  free2 = free1 = 1;

  if(array[2] == array[1] || array[2] == array[0]) free2 = 0;
  if(array[1] == array[0]) free1 = 0;

  if(free2) free_array(array[2]);
  if(free1) free_array(array[1]);
  free_array(array[0]);
}

/* copy n bytes from a tArray (at bytestride-position pos) into dest */
void *memcpy_from_array(const tArray *ar, size_t bytestride, size_t pos,
                        void *dest, size_t n)
{
  return memcpy(dest, ar->c + bytestride*pos, n);
}

/* copy n bytes into a tArray at bytestride-position pos */
void *memcpy_to_array(tArray *ar, size_t bytestride, size_t pos,
                      const void *src, size_t n)
{
  return memcpy(ar->c + bytestride*pos, src, n);
}
/* same as memcpy_to_array but redim array ar */
void *memcpy_to_array_redim(tArray *ar, size_t bytestride, size_t pos,
                            const void *src, size_t n)
{
  size_t nbs = (n+bytestride-1)/bytestride; //num. bytestrides in n bytes
  size_t nbytes0 = bytestride*pos + n;  //min number of bytes needed in ar->c
  size_t nbytes = bytestride*(pos+nbs); //number of bytes used in redim
  size_t sd = sizeof(ar->d[0]);         //sizeof double
  size_t nd = (nbytes+sd-1)/sd;         //num. of doubles in nbytes

  //printf("bytestride=%zu pos=%zu nbytes=%zu sd=%zu nd=%zu\n",
  //        bytestride, pos, nbytes, sd, nd);
  redim_array(ar, nd,1,1);

  /* zero part of ar->c that memcpy will not write to */
  memset(ar->c + bytestride*pos + n, 0, nbytes-nbytes0);
  /* copy src into ar->c */
  return memcpy(ar->c + bytestride*pos, src, n);
}

/* return number of tEplocs in array */
int array_Neplocs(tArray *ar)
{
  size_t se = sizeof(ar->eploc[0]);  //sizeof eploc
  size_t sd = sizeof(ar->d[0]);      //sizeof double
  size_t nbytes = sd * ar->N;        //number of bytes in ar
  return (nbytes)/se;                //number of eplocs in nbytes
}

/****************************************************************************/
/* elm storage */
/****************************************************************************/

/* allocate one elm */
tElm *alloc_elm(tMesh *mesh)
{
  tElm *elm = calloc(1, sizeof(*elm));
  if(!elm) errorexit("out of memory");

  //FIXME: once elm has mesh, set it here
  //elm->mesh = mesh;

  elm->eploc->eid = EID_INVALID;  /* mark eid as not set */

  return elm;
}

/* also init pat: REMOVE this func once we do not need elm->pat anymore */
tElm *alloc_elm_init_pat(tMesh *mesh, int p)
{
  tElm *elm = alloc_elm(mesh);
  elm->pat = mesh->pat[p];
  return elm;
}

/* alloc a new elm using ELMHEADER info */
tElm *alloc_elm_of_elmheader(tMesh *mesh, tElm0 *elmheader)
{
  int p = elmheader->eploc->p;
  tElm *elm = alloc_elm(mesh);

  /* first copy entire elmheader into new elm */
  memcpy(elm, elmheader, sizeof(elmheader[0]));

  /* now set anything else we need */
  //FIXME: once elm has mesh, set it here
  //elm->mesh = mesh;

  elm->pat = mesh->pat[p];
  return elm;
}

/* free one elm and its dat */
void free_elm(tElm *elm)
{
  int face;

  if(!elm) return;

  /* free variable data */
  free_dat(elm->dat);

  /* free surface neigbhor list */
  for(face=0; face<6; face++) free(elm->fnb[face]);

  free(elm);
}


/* allocate and set mesh->myelm array from mesh->myelm_head */
long alloc_and_set_mesh_myelm(tMesh *mesh)
{
  long ei;
  struct list_head *pos;

  /* free whatever we had so far in myelm */
  free(mesh->myelm);
  mesh->myelm = NULL;

  /* get number of elms on this rank */
  mesh->nmyelm = list_count_nodes(&mesh->myelm_head);
  PRF;printf(": mesh->nmyelm=%ld\n", mesh->nmyelm);
  /* return if there are no elms */
  if(mesh->nmyelm == 0) return 0;

  /* alloc sufficient space */
  mesh->myelm = calloc(mesh->nmyelm, sizeof(mesh->myelm[0]));
  if(!mesh->myelm) errorexit("no memory for mesh->myelm");

  /* store entries of list mesh->myelm_head in array mesh->myelm */
  ei = 0;
  list_for_each(pos, &mesh->myelm_head)
  {
    mesh->myelm[ei] = list_entry(pos, tElm, list);
    ei++;
  }
  if(ei != mesh->nmyelm) errorexit("Whaaaat??????");
  return mesh->nmyelm;
}

/* free array mesh->myelm array */
void free_mesh_myelm(tMesh *mesh)
{
  free(mesh->myelm);
  mesh->myelm = NULL;
  mesh->nmyelm = 0;
}


/* make root node element and add it to mesh */
tElm *make_and_add_root_elm(tPat *pat, int n[3], int pt_typ[3], int datrank)
{
  /* make root element only on the MPI rank that owns it */
  if(nMPI_rank()==datrank)
  {
    tMesh *mesh = pat->mesh;
    tElm *elm = alloc_elm(mesh);
    tEploc *eploc = elm->eploc;
    int i;

    /* check for overflow in eploc->p */
    i = (1 << (sizeof(eploc->p)*8)) - 1;
    if(pat->p > i) errorexiti("cannot have more than %d patches", i);

    /* fill in info */
    eploc->p = pat->p;
    eploc->l = 0; /* root node */
    eploc->ploc[0] = 0;
    amr_set_elm_pat(mesh, elm);
    amr_set_elm_bbox(elm);

    /* save n and pt_typ for root node */
    for(i=0; i<3; i++)
    {
      elm->n[i] = n[i];
      elm->pt_typ[i] = pt_typ[i];
    }
    elm->np = n[0] * n[1] * n[2];
    elm->eploc->eid = EID_INVALID;    /* mark eid as not set */

    /* see where dat needs to be allocated */
    elm->datrank = datrank;
    elm->dat = alloc_dat(elm);

    /* first set fnb */
    //update_node_fnb_only(node);
    //FIXME ???

    /* add new root element to list mesh->myelm_head */
    list_add_tail(&elm->list, &mesh->myelm_head);

    return elm;
  }
  else
  {
    return NULL;
  }
}


/* make a child node element */
tElm *make_child_elm(tElm *parent, int n[3], int pt_typ[3], int ijk)
{
  tMesh *mesh = parent->pat->mesh;
  tElm *elm = alloc_elm(mesh);
  int d, vi,nvdb;

  /* transfer parent time info */
  elm->time = parent->time;
  elm->dt = parent->dt;  // FIXME: For now all elms have same dt

  /* mark eid as not set */
  elm->eploc->eid = EID_INVALID;

  /* fill in info */
  amr_set_child_eploc(parent->eploc, ijk, elm->eploc);
  amr_set_elm_pat(mesh, elm);
  amr_set_elm_bbox(elm);

  for(d=0; d<3; d++)
  {
    elm->n[d] = n[d];
    elm->pt_typ[d] = pt_typ[d]; /* save point type */
  }
  elm->np = n[0] * n[1] * n[2];

  nvdb = elm->pat->mesh->nvdb;

  /* default is same proc as parent */
  elm->datrank = parent->datrank;

  /* if parent has dat the child will have it too */
  if(parent->dat)
  {
    tArray *Xp[3];

    /* alloc dat for child */
    elm->dat = alloc_dat(elm);

    /* array memory to store points of elm */
    Xp[0] = alloc_array(n);
    Xp[1] = alloc_array(n);
    Xp[2] = alloc_array(n);
    fill_3arrays_with_nodepoints(elm, Xp);
    /* convert from Xb of elm to X to Xb of parent */
    array_XYZ_of_XbYbZb(elm, Xp, Xp);
    array_XbYbZb_of_XYZ(parent, Xp, Xp);

    /* use interpolation to get vars from parent to child elm */
    for(vi=0; vi<nvdb; vi++)
      if(parent->dat->v[vi])
      {
        int vt = MeshVarType(mesh, vi);
        /* enable same vars in this dat as in parent->dat */
        enablevarcomp_innode(elm, vi);

        /* fill elm->dat with interpolation data from parent */
        if( (vt==EVOVAR) || (vt==DATAVAR) ) /* exclude Aux. vars */
        {
          basis_interp_topoints(parent, parent->dat->v[vi],
                                Xp, elm->dat->v[vi], Lagrange_of_x);
        }
      } /* end: if parent has dat->v[vi] */
    free_array(Xp[2]);
    free_array(Xp[1]);
    free_array(Xp[0]);

    /* init coords in this new elm */
    coordinates_init_node(elm);
  }
  return elm;
}

/* make 8 children, insert them into mesh->myelm_head, and return child0 */
tElm *replace_parent_by_8children(tElm *parent, int n[3], int pt_typ[3])
{
  //tMesh *mesh = parent->pat->mesh;
  struct list_head elist;
  tElm *elm, *elm0;
  int ijk;

  INIT_LIST_HEAD(&elist);

  /* make children */
  for(ijk=0; ijk<8; ijk++)
  {
    elm = make_child_elm(parent, n, pt_typ, ijk);
    if(elm) list_add_tail(&elm->list, &elist);
    if(ijk==0) elm0 = elm; /* save first child */
  }

  /* #pragma omp critical (change_mesh_myelm_list) */
  /* NOTE: For some reason gcc's -fsanitize=thread throws a ?false? positive
           if I use a named critical section!
           So replace "GEN_Pragma(omp critical (change_mesh_myelm_list))"
           by "GEN_Pragma(omp critical)" when debugging races!!! */
  //GEN_Pragma(omp critical (change_mesh_myelm_list))
  GEN_Pragma(omp critical)
  {
    /* now replace parent by elist in mesh->myelm_head */
    list_splice(&elist, &parent->list);
    list_del(&parent->list);
  }

  /* free parent and all data on parent */
  free_elm(parent);

  //printf("Created:\n");
  //printnodes_in_list(nlist);

  return elm0;
}

// Equivalent of
// tNode *destroy_children(tNode *parent)
// FIXME: it also needs to be then used in refine.c
/* Remove 8 children and replace them by parent.
   We assume that all 8 have been moved to this rank before
   replace_8localchildren_by_parent is called! */
tElm *replace_8localchildren_by_parent(tElm *child0, int n[3], int pt_typ[3])
{
  tMesh *mesh = Elm_mesh(child0);
  tElm *parent = alloc_elm(mesh);
  int d;

  /* transfer child0 time info */
  parent->time = child0->time;
  parent->dt = child0->dt;  // FIXME: For now all elms have same dt

  /* fill in info */
  amr_set_parent_eploc(child0->eploc, parent->eploc);
  amr_set_elm_pat(mesh, parent);
  amr_set_elm_bbox(parent);

  /* mark eid as not set */
  parent->eploc->eid = EID_INVALID;

  for(d=0; d<3; d++)
  {
    parent->n[d] = n[d];
    parent->pt_typ[d] = pt_typ[d]; /* save point type */
  }
  parent->np = n[0] * n[1] * n[2];

  /* set parent's datrank to the same as child0 */
  parent->datrank = child0->datrank;
  if(child0->dat)
  {
    int nvdb = mesh->nvdb;
    int vi, ijk;
    tArray *Xp[3], *Ip[8], *Xc[8][3], *Res[8];
    struct list_head *pos_ijk;

    if(!parent->dat) parent->dat = alloc_dat(parent);

    /* enable same vars in parent->dat as in child0->dat */
    for(vi=0; vi<nvdb; vi++)
      if(child0->dat->v[vi]) enablevarcomp_innode(parent, vi);

    /* but disable nb-info vars */
    enablevar_innode(parent, amr->elm_nbinfo0);

    /* array memory to store points of parent in X coords */
    Xp[0] = alloc_array(parent->n);
    Xp[1] = alloc_array(parent->n);
    Xp[2] = alloc_array(parent->n);
    fill_3arrays_with_nodepoints(parent, Xp);
    /* convert from Xb of parent to X for parent,
       these X are spread over the 8 child nodes */
    array_XYZ_of_XbYbZb(parent, Xp, Xp);

    /* fill parent->dat with interpolation data from children */
    /* 1. set children coords within parent */
    pos_ijk = &child0->list;
    for(ijk=0; ijk<8; ijk++)
    {
      tElm *child = list_entry(pos_ijk, tElm, list);

      /* array memory to store points of child in Xb coords */
      Ip[ijk] = alloc_array(parent->n);
      Xc[ijk][0] = alloc_array(parent->n);
      Xc[ijk][1] = alloc_array(parent->n);
      Xc[ijk][2] = alloc_array(parent->n);
      Res[ijk] = alloc_array(parent->n);

      /* find points inside child node -> mask is returned in Ip */
      array_find_XYZ_in_node(child, Xp, Ip[ijk]);
      // NOTE: Parent points on the boundary of child are found in several
      // children (several ijk). Res[ijk] will contain the result
      // interpolated from child ijk. This is a problem if data in
      // children is not smooth. We need to average somehow.

      /* convert Xp to child's internal basis coords */
      array_XbYbZb_of_XYZ(child, Xc[ijk], Xp);

      /* pos of next child */
      pos_ijk = pos_ijk->next;
    }

    /* 2. use interpolation to get vars from child to parent */
    for(vi=0; vi<nvdb; vi++)
    {
      int vt = MeshVarType(mesh, vi);
      /* fill parent->dat with interpolation data from child */
      if( (vt==EVOVAR) || (vt==DATAVAR) ) /* exclude Aux. vars */
      {
        if(child0->dat->v[vi])
        {
          int k, cnt;

          /* interpolate for each child and save results in Res */
          pos_ijk = &child0->list;
          for(ijk=0; ijk<8; ijk++)
          {
            tElm *child = list_entry(pos_ijk, tElm, list);
            tArray *var = child->dat->v[vi];
            if(var)
              basis_interp_toIpoints(child, var, Xc[ijk],Ip[ijk], Res[ijk],
                                     Lagrange_of_x);
            /* pos of next child */
            pos_ijk = pos_ijk->next;
          }

          /* take average of results from different child nodes */
          forarray(parent->dat->v[vi], k)
          {
            cnt = 0;
            parent->dat->v[vi]->d[k] = 0.;
            for(ijk=0; ijk<8; ijk++)
            {
              if(Ip[ijk]->i[k] >= 0) /* if child has the point */
              {
                cnt++; /* count num of chidren who have this point */
                parent->dat->v[vi]->d[k] += Res[ijk]->d[k];
              }
            }
            /* average if there was more than one child with this point */
            if(cnt>1) parent->dat->v[vi]->d[k] /= cnt;
          }
        } /* end: if(child0->dat) */
      }
    }
    /* 3. free temp arrays for coords */
    for(ijk=0; ijk<8; ijk++)
    {
      free_array(Res[ijk]);
      free_array(Xc[ijk][2]);
      free_array(Xc[ijk][1]);
      free_array(Xc[ijk][0]);
      free_array(Ip[ijk]);
    }
    free_array(Xp[2]);
    free_array(Xp[1]);
    free_array(Xp[0]);

    /* init coords in parent */
    coordinates_init_node(parent);
  }

  /* NOTE: This new parent has all zero for nfnb, fnb, and nbinfo.
           Also, all its neighbors have now the wrong nfnb and nbinfo.
           Even worse, all its neighbors have fnb pointers pointing
           to the children which will be removed!!!  */

  /* NOTE: Should we go over children's nbs and set whatever
           nb-info we can??? */


  /* replace children by parent mesh->myelm_head list */
  /* #pragma omp critical (change_mesh_myelm_list) */
  /* NOTE: For some reason gcc's -fsanitize=thread throws a ?false? positive
           if I use a named critical section!
           So replace "GEN_Pragma(omp critical (change_mesh_myelm_list))"
           by "GEN_Pragma(omp critical)" when debugging races!!! */
  //GEN_Pragma(omp critical (change_mesh_myelm_list))
  GEN_Pragma(omp critical)
  {
    int ijk;
    /* now replace children by parent in mesh->myelm_head */
    /* first remove child 1-7 */
    for(ijk=1; ijk<8; ijk++)
    {
      struct list_head *pos_ijk = (child0->list).next; //pos after child0
      tElm *ch_ijk = list_entry(pos_ijk, tElm, list);  //child after child0
      list_del(&ch_ijk->list);
    }
    /* now replace child0 by parent in mesh->myelm_head */
    list_add(&parent->list, &child0->list);
    list_del(&child0->list);
  }

  /* NOTE: right now the children have only been removed from the
           mesh->myelm_head list but otherwise are still there! */
  return parent;
}



/**************************************************************************/
/* node storage */
/**************************************************************************/

/* allocate one node*/
tNode *alloc_node(int initcomm)
{
  tNode *node = calloc(1, sizeof(*node));
  if(!node) errorexit("out of memory");

  // /* set node MPI communicator node->comm */
  //if(initcomm)
  //  node->comm = nMPIvars_get_comm(initcomm-1);
  //else
  //  node->comm = main_comm;

  return node;
}

/* free one node only, leaves children hanging */
void free_this_node_only(tNode *node)
{
  tNode *parent = node->parent;
  int ijk, face;

  if(!node) return;

  /* free surface neigbhor list */
  for(face=0; face<6; face++) free(node->fnb[face]);

  /* remove nfaces */
  remove_all_nfaces(node);

  /* remove parent's pointer to it */
  ijk = node->ijk;
  if(parent) parent->child[ijk] = NULL;

  /* should we also remove pointer of neighbors to it??? */


  /* free variable data */
  free_dat(node->dat);

  free(node);
}

/* free one node and all its children */
void free_node(tNode *node)
{
  tNode *chld;
  int ijk;

  //PRF;printf(": nid=%d l=%d node=%p\n", get_node_nid(node), node->l, node);
  //printnode(node);
  if(!node) return;

  for(ijk=0; ijk<8; ijk++)
  {
    chld = node->child[ijk];
    if(chld) free_node(chld);
  }
  //PRF;printf(": nid=%d l=%d node=%p\n", get_node_nid(node), node->l, node);
  free_this_node_only(node);
}


/* make root node */
tNode *make_root_node(tPat *pat, int pt_typ[3], int n[3], int datrank)
{
  tNode *node = alloc_node(0);
  int i;

  /* fill in info */
  node->pat = pat;
  /* node->nb is left uninitialized here !!! */
  /* we assume that it has no neighbors */

  /* take bounding boxes from pat, and mark it as touching all patch faces */
  for(i=0; i<6; i++)
  {
    node->bbox[i] = pat->bbox[i];
    node->patface[i] = 1;
  }

  /* save n and pt_typ for root node */
  for(i=0; i<3; i++)
  {
    node->n[i] = n[i];
    node->pt_typ[i] = pt_typ[i];
  }
  node->np = n[0] * n[1] * n[2];
  node->l = 0;
  node->leaf = 1;    /* make this a leaf node */
  node->nid = -1;    /* mark nid as not set */

  /* see where dat needs to be allocated */
  node->datrank = datrank;
  if(nMPI_rank()==datrank)
    node->dat = alloc_dat(node);

  /* lock new node and its face nbs, to do so first set fnb */
  update_node_fnb_only(node);
  node_and_fnbs_lock(node);

  /* initialize surface neigbhor list in node and it neighbors */
  update_node_and_neighbors_nfaces_fnb(node);

  /* unlock new node and its face nbs */
  node_and_fnbs_unlock(node);

  return node;
}

/* make a child node */
tNode *make_child_node(tNode *parent, int pt_typ[3], int n[3], int ijk)
{
  tMesh *mesh = parent->pat->mesh;
  tNode *node = alloc_node(0);
  double mid[3];
  int i,j,k, d, vi,nvdb;
  int ns[] = {2,2,2};

  /* transfer parent time info */
  node->time = parent->time;
  node->dt = parent->dt;  // FIXME: For now all nodes have same dt

  /* mark nid as not set */
  node->nid = -1;

  /* register this child with the parent */
  parent->child[ijk] = node;
  parent->leaf = 0;  /* parent is now no longer a leaf node */

  /* node coords from node index ijk */
  k = kOfInd_n(ijk, ns);
  j = jOfInd_n_k(ijk, ns,k);
  i = iOfInd_n_jk(ijk, ns,j,k);

  /* mid point in parent node */
  for(d=0; d<3; d++)
    mid[d] = 0.5*(parent->bbox[2*d] + parent->bbox[2*d+1]);

  /* set new bounding boxes */
  /* at first take bbox from pat, patface from parent */
  for(d=0; d<6; d++)
  {
    node->bbox[d]    = parent->bbox[d];
    node->patface[d] = parent->patface[d];
  }
  if(i%2) { node->bbox[0] = mid[0]; node->patface[0] = 0; }
  else    { node->bbox[1] = mid[0]; node->patface[1] = 0; }
  if(j%2) { node->bbox[2] = mid[1]; node->patface[2] = 0; }
  else    { node->bbox[3] = mid[1]; node->patface[3] = 0; }
  if(k%2) { node->bbox[4] = mid[2]; node->patface[4] = 0; }
  else    { node->bbox[5] = mid[2]; node->patface[5] = 0; }

  /* fill in info */
  node->pat    = parent->pat;
  node->parent = parent;
  /* node->nb is left uninitialized here !!! */

  for(d=0; d<3; d++)
  {
    node->n[d] = n[d];
    node->pt_typ[d] = pt_typ[d]; /* save point type */
  }
  node->np = n[0] * n[1] * n[2];

  node->l = parent->l + 1;
  node->leaf = 1;    /* make this a leaf node */
  node->ijk = ijk;
  nvdb = node->pat->mesh->nvdb;

  /* default is same proc as parent */
  node->datrank = parent->datrank;

  /* if parent has dat the child will have it too */
  if(parent->dat)
  {
    tArray *Xp[3];

    /* alloc dat for child */
    node->dat = alloc_dat(node);

    /* array memory to store points of node */
    Xp[0] = alloc_array(n);
    Xp[1] = alloc_array(n);
    Xp[2] = alloc_array(n);
    fill_3arrays_with_nodepoints(node, Xp);
    /* convert from Xb of node to X to Xb of parent */
    array_XYZ_of_XbYbZb(node, Xp, Xp);
    array_XbYbZb_of_XYZ(parent, Xp, Xp);

    /* use interpolation to get vars from parent to child node */
    for(vi=0; vi<nvdb; vi++)
      if(parent->dat->v[vi])
      {
        int vt = MeshVarType(mesh, vi);
        /* enable same vars in this dat as in parent->dat */
        enablevarcomp_innode(node, vi);

        /* fill node->dat with interpolation data from parent */
        if( (vt==EVOVAR) || (vt==DATAVAR) ) /* exclude Aux. vars */
        {
          basis_interp_topoints(parent, parent->dat->v[vi],
                                Xp, node->dat->v[vi], Lagrange_of_x);
        }
      } /* end: if parent has dat->v[vi] */
    free_array(Xp[2]);
    free_array(Xp[1]);
    free_array(Xp[0]);

    /* init coords in this new node */
    coordinates_init_node(node);
  }
  return node;
}

/* make 8 childern and return them in a short list */
tNlist *make8_child_nodes(tNode *parent, int pt_typ[3], int n[3])
{
  tNlist *nlist = NULL;
  tNlist *elem = NULL;
  tNode parent_tmp[1];
  char *pt1, *pt2;
  tNode *node;
  tNode *narray[8];
  int ijk;

  /* shallow copy of parent,
     but exclude nb info parts to not throw off data race detectors */
  //parent_tmp[0] = parent[0];
  pt1 = (char *) parent;
  pt2 = (char *) &(parent->nb[0]);
  memcpy(&(parent_tmp[0]), parent, pt2-pt1);

  /* use parent_tmp to make children */
  for(ijk=0; ijk<8; ijk++)
  {
    node = make_child_node(parent_tmp, pt_typ, n, ijk);
    elem = addnode_to_nodelist_after(elem, node);
    if(ijk==0) nlist = elem; // save list head
    narray[ijk] = node; /* save nodes also in an array */
  }

  /*
  Below we use critical sections because the locking of nearest neighbors
  in node_and_fnbs_lock (while probably working) is not enough.
  update_node_and_neighbors_nfaces_fnb calls add_nfaces_outside_patch, which
  in turn calls leafdescendants_along_face. Yet leafdescendants_along_face
  is called with the root node! From there it descends to all nodes along
  the face of a patch. Many of these nodes may not be locked and are far
  away from parent. I.e. we access nodes that other threads are working on
  so that we get a race condition! So for now we let only one thread go
  here. However, a better solution would be to improve
  add_nfaces_outside_patch such that it starts from some known neighbors in
  the other patch and not from the root node. Then there is a chance to keep
  the descends local. So that locking nearest neighbors would be effective.
  */

  /* #pragma omp critical (make_or_destroy_nodes) */
  /* NOTE: For some reason gcc's -fsanitize=thread throws a ?false? positive
           if I use a named critical section!
           So replace "GEN_Pragma(omp critical (make_or_destroy_nodes))"
           by "GEN_Pragma(omp critical)" when debugging races!!! */
  //GEN_Pragma(omp critical (make_or_destroy_nodes))
  GEN_Pragma(omp critical)
  {
    ///* aquire lock for change of connections */
    //node_and_fnbs_lock(parent);

    /* now attach children to parent */
    // parent[0] = parent_tmp[0];
    parent->leaf = 0;  // parent is now no longer a leaf node
    for(ijk=0; ijk<8; ijk++)
    {
      parent->child[ijk] = parent_tmp->child[ijk];
      parent->child[ijk]->parent = parent;
    }

    /* fill in neighbor info, as far as these 8 are concerned */
    connect8_with_neighbors(narray, 1);

    /* update fnb of all in narray and their neighbors */
    for(ijk=0; ijk<8; ijk++)
      update_node_and_neighbors_nfaces_fnb(narray[ijk]);

    ///* release locks */
    //node_and_fnbs_unlock(parent);
  }

  /* free all data on parent */
  free_dat(parent->dat);
  parent->dat = NULL;

  //printf("Created:\n");
  //printnodes_in_list(nlist);

  return nlist;
}

/* update node->n (and node->pt_typ if pt_typ != NULL) on one node,
   { int *n, int *pt_typ are really int n[3], int pt_typ[3] },
   should be called for all 8 siblings */
void update_node_n_pt_typ__old(tNode *node, int *n, int *pt_typ)
{
  tMesh *mesh = node->pat->mesh;
  int nvdb = mesh->nvdb;
  tNode node_old[1];
  int d, vi;

  /* backup old node info */
  memcpy(node_old, node, sizeof(node_old[0]));

  /* update node info */
  if(n)
  {
    for(d=0; d<3; d++) node->n[d] = n[d];
    node->np = n[0] * n[1] * n[2];
  }
  if(pt_typ)
    for(d=0; d<3; d++) node->pt_typ[d] = pt_typ[d];

  /* if node has dat, we need to interpolate vars */
  if(node_old->dat)
  {
    tArray *Xp[3];

    /* alloc new dat for node */
    node->dat = alloc_dat(node);

    /* array memory to store points of node */
    Xp[0] = alloc_array(n);
    Xp[1] = alloc_array(n);
    Xp[2] = alloc_array(n);
    fill_3arrays_with_nodepoints(node, Xp);
    /*FIXME: I think instead of alloc and fill_3arrays_with_nodepoints, we
      could use: node_Xb3(node, Xp); to get new node points into Xp */

    /* use interpolation to get vars from old dat to new node->dat */
    for(vi=0; vi<nvdb; vi++)
      if(node_old->dat->v[vi])
      {
        int vt = MeshVarType(mesh, vi);
        /* enable same vars in new dat as in dat_old */
        enablevarcomp_innode(node, vi);

        /* fill node->dat with interpolation data from old dat */
        if( (vt==EVOVAR) || (vt==DATAVAR) ) /* exclude Aux. vars */
        {
          basis_interp_topoints(node_old, node_old->dat->v[vi],
                                Xp, node->dat->v[vi], Lagrange_of_x);
        }
      } /* end: if parent has dat->v[vi] */
    free_array(Xp[2]);
    free_array(Xp[1]);
    free_array(Xp[0]);
    free_dat(node_old->dat);

    /* init coords in this new node */
    coordinates_init_node(node);
  }
}

/* update node->n (and node->pt_typ if pt_typ != NULL) on one node,
   { int *n, int *pt_typ are really int n[3], int pt_typ[3] },
   should be called for all 8 siblings
   NOTE: After calling update_node_n_pt_typ_return_node_old we have to call:
           update_node_n_pt_typ_free_node_old(node_old);
         to free the copy node_old that is allocated and returned here. */
tNode *update_node_n_pt_typ_return_node_old(tNode *node, int *n, int *pt_typ)
{
  tMesh *mesh = node->pat->mesh;
  int nvdb = mesh->nvdb;
  int d, vi;
  tNode *node_old = calloc(1, sizeof(*node_old));
  if(!node_old) errorexit("no memory for node_old");

  /* backup old node info */
  memcpy(node_old, node, sizeof(node_old[0]));

  /* update node info */
  if(n)
  {
    for(d=0; d<3; d++) node->n[d] = n[d];
    node->np = n[0] * n[1] * n[2];
  }
  if(pt_typ)
    for(d=0; d<3; d++) node->pt_typ[d] = pt_typ[d];

  /* if node has dat, we need to interpolate vars */
  if(node_old->dat)
  {
    tArray *Xp[3];

    /* alloc new dat for node */
    node->dat = alloc_dat(node);

    /* array memory to store points of node */
    Xp[0] = alloc_array(n);
    Xp[1] = alloc_array(n);
    Xp[2] = alloc_array(n);
    fill_3arrays_with_nodepoints(node, Xp);
    /*FIXME: I think instead of alloc and fill_3arrays_with_nodepoints, we
      could use: node_Xb3(node, Xp); to get new node points into Xp */

    /* use interpolation to get vars from old dat to new node->dat */
    for(vi=0; vi<nvdb; vi++)
      if(node_old->dat->v[vi])
      {
        int vt = MeshVarType(mesh, vi);
        /* enable same vars in new dat as in dat_old */
        enablevarcomp_innode(node, vi);

        /* fill node->dat with interpolation data from old dat */
        if( (vt==EVOVAR) || (vt==DATAVAR) ) /* exclude Aux. vars */
        {
          basis_interp_topoints(node_old, node_old->dat->v[vi],
                                Xp, node->dat->v[vi], Lagrange_of_x);
        }
      } /* end: if parent has dat->v[vi] */
    free_array(Xp[2]);
    free_array(Xp[1]);
    free_array(Xp[0]);

    /* shallow copy of node_old->dat->info, to get e.g. load timers */
    node->dat->info[0] = node_old->dat->info[0];

    /* init coords in this new node */
    coordinates_init_node(node);
  }
  return node_old;
}

/* free the node_old returned by update_node_n_pt_typ_return_node_old */
void update_node_n_pt_typ_free_node_old(tNode *node, tNode *node_old)
{
  if(node_old->dat != node->dat) free_dat(node_old->dat);
  free(node_old);
}

/* Undo a call to update_node_n_pt_typ_return_node_old. This still does not
   free node_old! */
void update_node_n_pt_typ_restore_from_node_old(tNode *node, tNode *node_old)
{
  /* free any newly allocated dat in node */
  free_dat(node->dat);

  /* use info in node_old to restore node */
  memcpy(node, node_old, sizeof(node[0]));
}

/* update node->n (and node->pt_typ if pt_typ != NULL) on one node,
   { int *n, int *pt_typ are really int n[3], int pt_typ[3] },
   should be called for all 8 siblings */
void update_node_n_pt_typ(tNode *node, int *n, int *pt_typ)
{
  tNode *node_old;
  node_old = update_node_n_pt_typ_return_node_old(node, n, pt_typ);
  update_node_n_pt_typ_free_node_old(node, node_old);
}

/* update node->n on one node, should be called for all 8 siblings */
void update_node_n(tNode *node, int n[3])
{
  update_node_n_pt_typ(node, n, NULL);
}

/* update node->n (and possibly node->pt_typ) on all 8 siblings,
   { int *n, int *pt_typ are really int n[3], int pt_typ[3] },
   must be called by all MPI procs */
void update8_node_n_pt_typ(tNode *node, int *n, int *pt_typ)
{
  tNode *parent = node->parent;

  /* update all 8 siblings, unless this is root node */
  if(parent)
  {
    int ijk;
    for(ijk=0; ijk<8; ijk++)
    {
      tNode *sib = parent->child[ijk];
      update_node_n_pt_typ(sib, n, pt_typ);
    }
  }
  else
  {
    update_node_n_pt_typ(node, n, pt_typ);
  }
}

/* update node->n on all 8 siblings, must be called by all MPI procs */
void update8_node_n(tNode *node, int n[3])
{
  update8_node_n_pt_typ(node, n, NULL);
}

/* remove children */
tNode *destroy_children(tNode *parent)
{
  tNlist *el, *clist = NULL;
  tNode *narray[8];
  tNode *child0 = parent->child[0];
  int ijk;

  /* set parent time to the same as child0 */
  parent->time = child0->time;

 /* save children in an array and nodelist clist */
  el = NULL;
  for(ijk=0; ijk<8; ijk++)
  {
    narray[ijk] = parent->child[ijk];
    el = addnode_to_nodelist_after(el, parent->child[ijk]);
    if(ijk==0) clist = el; /* save begin of nodelist */
  }
  /* now move these 8 children to rank of child0.
     This is needed so that the interpolation below can work within
     one rank without MPI calls. */
  move_nodelist_to_rank(clist, child0->datrank);
  free_nodelist(clist); /* we do not need clist anymore */

  /* set parent's datrank to the same as child0 */
  parent->datrank = child0->datrank;
  if(child0->dat)
  {
    tMesh *mesh = parent->pat->mesh;
    int nvdb = mesh->nvdb;
    int vi;
    tArray *Xp[3], *Ip[8], *Xc[8][3], *Res[8];

    if(!parent->dat) parent->dat = alloc_dat(parent);
    /* enable same vars in parent->dat as in child0->dat */
    for(vi=0; vi<nvdb; vi++)
      if(child0->dat->v[vi]) enablevarcomp_innode(parent, vi);

    /* array memory to store points of parent in X coords */
    Xp[0] = alloc_array(parent->n);
    Xp[1] = alloc_array(parent->n);
    Xp[2] = alloc_array(parent->n);
    fill_3arrays_with_nodepoints(parent, Xp);
    /* convert from Xb of parent to X for parent,
       these X are spread over the 8 child nodes */
    array_XYZ_of_XbYbZb(parent, Xp, Xp);

    /* fill parent->dat with interpolation data from children */
    /* 1. set children coords within parent */
    for(ijk=0; ijk<8; ijk++)
    {
      tNode *child = narray[ijk];

      /* array memory to store points of child in Xb coords */
      Ip[ijk] = alloc_array(parent->n);
      Xc[ijk][0] = alloc_array(parent->n);
      Xc[ijk][1] = alloc_array(parent->n);
      Xc[ijk][2] = alloc_array(parent->n);
      Res[ijk] = alloc_array(parent->n);

      /* find points inside child node -> mask is returned in Ip */
      array_find_XYZ_in_node(child, Xp, Ip[ijk]);
      // NOTE: Parent points on the boundary of child are found in several
      // children (several ijk). Res[ijk] will contain the result
      // interpolated from child ijk. This is a problem if data in
      // children is not smooth. We need to average somehow.

      /* convert Xp to child's internal basis coords */
      array_XbYbZb_of_XYZ(child, Xc[ijk], Xp);
    }

    /* 2. use interpolation to get vars from child to parent */
    for(vi=0; vi<nvdb; vi++)
    {
      int vt = MeshVarType(mesh, vi);
      /* fill parent->dat with interpolation data from child */
      if( (vt==EVOVAR) || (vt==DATAVAR) ) /* exclude Aux. vars */
      {
        if(child0->dat->v[vi])
        {
          int k, cnt;

          /* interpolate for each child and save results in Res */
          for(ijk=0; ijk<8; ijk++)
          {
            tNode *child = narray[ijk];
            tArray *var = child->dat->v[vi];
            if(var)
              basis_interp_toIpoints(child, var, Xc[ijk],Ip[ijk], Res[ijk],
                                     Lagrange_of_x);
          }

          /* take average of results from different child nodes */
          forarray(parent->dat->v[vi], k)
          {
            cnt = 0;
            parent->dat->v[vi]->d[k] = 0.;
            for(ijk=0; ijk<8; ijk++)
            {
              if(Ip[ijk]->i[k] >= 0) /* if child has the point */
              {
                cnt++; /* count num of chidren who have this point */
                parent->dat->v[vi]->d[k] += Res[ijk]->d[k];
              }
            }
            /* average if there was more than one child with this point */
            if(cnt>1) parent->dat->v[vi]->d[k] /= cnt;
          }
        } /* end: if(child0->dat) */
      }
    }
    /* 3. free temp arrays for coords */
    for(ijk=0; ijk<8; ijk++)
    {
      free_array(Res[ijk]);
      free_array(Xc[ijk][2]);
      free_array(Xc[ijk][1]);
      free_array(Xc[ijk][0]);
      free_array(Ip[ijk]);
    }
    free_array(Xp[2]);
    free_array(Xp[1]);
    free_array(Xp[0]);

    /* init coords in parent */
    coordinates_init_node(parent);
  }

  /*
  The long comment before critical section in make8_child_nodes applies here
  as well. add_nfaces_outside_patch makes the leafdescendants_along_face
  start with the root node. Thus we access far away nodes that cannot be
  locked with parent_and_fnbs_lock.
  */

  /* #pragma omp critical (make_or_destroy_nodes) */
  /* NOTE: For some reason gcc's -fsanitize=thread throws a ?false? positive
           if I use a named critical section!
           So replace "GEN_Pragma(omp critical (make_or_destroy_nodes))"
           by "GEN_Pragma(omp critical)" when debugging races!!! */
  //GEN_Pragma(omp critical (make_or_destroy_nodes))
  GEN_Pragma(omp critical)
  {
    ///* obtain lock on face neighbors of narray in name of parent */
    //parent_and_fnbs_lock(narray, parent);

    /* update neighbor info */
    /* set neighbor info to NULL, as far as these 8 are concerned */
    connect8_with_neighbors(narray, 0);

    /* unlink child nodes, and remove all nfaces of narray */
    for(ijk=0; ijk<8; ijk++)
    {
      if(narray[ijk]->child[0])
        errorexit("cannot orphan child that itself has child[0]");
      parent->child[ijk] = NULL;
      remove_all_nfaces(narray[ijk]);
    }
    /* the narray is disconnected now, but still has lock on face neighbors */

    /* update fnb on parent and its neighbors */
    update_node_and_neighbors_nfaces_fnb(parent);

    /* parent is now a leaf node */
    parent->leaf = 1;

    /////* release lock on face neighbors */
    ////parent_and_fnbs_unlock(narray, parent);

    /* free child nodes */
    for(ijk=0; ijk<8; ijk++) free_node(narray[ijk]);

    ///* release lock on parent and its face neighbors */
    //node_and_fnbs_unlock(parent);
  }

  return parent;
}

/**************************************************************************/
/* patch and mesh storage */
/**************************************************************************/

/* allocate patch */
tPat *alloc_patch(tMesh *mesh, int p)
{
  tPat *pat;

  pat = calloc(1, sizeof(*pat));
  if(!pat) errorexit("out of memory");

  pat->mesh = mesh;
  pat->p = p;

  /* Bfaces */

  return pat;
}

/* free arrays in pat->CI*/
void free_pat_CI(tPat *pat)
{
  int d;

  for(d=0; d<6; d++) free_array(pat->CI->Fcoef[d]);
}

/* free pat, currently leaves mesh untouched */
void free_patch(tPat *pat)
{
  tMesh *mesh = pat->mesh;
  tNlist *elem;
  struct list_head *pos, *sav;

  if(!pat) return;

  //PRFs(":\n");

  /* free all in CI coordinfo, and also all bfaces  */
  free_pat_CI(pat);
  remove_all_bfaces(pat);


  /* remove all elms in this patch from mesh->myelm_head */
  list_for_each_safe(pos, sav, &mesh->myelm_head)
  {
    tElm *elm = list_entry(pos, tElm, list);
    if(elm->pat == pat)
    {
      list_del(&elm->list);
      free_elm(elm);
    }
  }


  /* remove all nodes in this patch from mesh->lns */
  /* 1. look at all except the head mesh->lns */
  for(elem=mesh->lns->next; elem; )
  {
    if(elem->node->pat == pat) elem = remove1_in_nodelist(elem, 1);
    else                       elem = elem->next;
  }
  /* 2. now examine the head */
  elem = mesh->lns;
  if(elem->node->pat == pat) elem = remove1_in_nodelist(elem, 1);
  /* and update the head */
  mesh->lns = first_nodelist(elem);

  /* update myln lists */
  update_mesh_myln_node_nid(mesh);

  /* free root node and all its children ... */
  free_node(pat->rnode);

  free(pat);
}



/* allocate mesh */
tMesh *alloc_mesh(int npats)
{
  tMesh *mesh;

  mesh = calloc(1, sizeof(*mesh));
  if(!mesh) errorexit("out of memory for mesh");

  mesh->eidlim = calloc(nMPI_size(), sizeof(mesh->eidlim[0]));

  realloc_patlist_in_mesh(mesh, npats);

  /* init list heads in mesh */
  INIT_LIST_HEAD(&mesh->myelm_head);

  /* init mesh mutex */
  MUTEX_INIT(mesh->mutex);

  return mesh;
}

/* make room for more patches */
void realloc_patlist_in_mesh(tMesh *mesh, int npats)
{
  int opats = mesh->npats;

//printf("npats=%d sss=%d\n", npats, sizeof(mesh->pat[0]));

  /* alloc list of pointers to patches */
  if(npats > opats)
  {
    mesh->pat = realloc(mesh->pat, npats*sizeof(mesh->pat[0]));
    if(!mesh->pat) errorexit("out of memory for mesh->pat");

    /* zero newly allocated part */
    memset(&(mesh->pat[opats]), 0, (npats-opats)*sizeof(mesh->pat[0]));
  }
  if(npats < opats)
  {
    int p;

    for(p=npats; p<opats; p++) free_patch(mesh->pat[p]);

    if(npats>0)
    {
      mesh->pat = realloc(mesh->pat, npats*sizeof(mesh->pat[0]));
      if(!mesh->pat) errorexit("shrinking mesh->pat failed");
    }
    else
    {
      free(mesh->pat);
      mesh->pat = NULL;
    }
  }
  mesh->npats = npats;
}

/* free mesh contents */
void free_mesh_patches_and_nodes(tMesh *mesh)
{
  int i;

  if(!mesh) return;

  if(PR) PRFs(":\n");

  /* free patches */
  for(i = 0; i < mesh->npats; i++)
    free_patch(mesh->pat[i]);

  /* now free patch list array */
  free(mesh->pat);

  /* node list in mesh */
  free_nodelist(mesh->lns);
  mesh->lns = NULL;
  realloc_myln_nncats(mesh->myln, 0);

  /* free mesh->myelm */
  free_mesh_myelm(mesh);

  /* set patch and node stuff to 0 */
  mesh->npats = 0;
  mesh->pat = NULL;
  mesh->lns = NULL;
  mesh->nln = 0;
  memset(&(mesh->myln[0]), 0, sizeof(mesh->myln[0]));
}

/* free mesh contents, except for MeshFuns that make up func skeleton */
void free_mesh_contents_exceptMeshFuns(tMesh *mesh)
{
  if(!mesh) return;

  if(PR) PRFs(":\n");

  /* free patches and nodes */
  free_mesh_patches_and_nodes(mesh);

  /* free vdb and pdb in mesh */
  free_mesh_vdb_contents(mesh);
  free(mesh->vdb);
  free_mesh_pdb_contents(mesh);
  free(mesh->pdb);

  /* free mesh mutex */
  MUTEX_DESTROY(mesh->mutex);

  /* free mesh->eidlim at the end */
  free(mesh->eidlim);
}

/* free all mesh contents */
void free_all_mesh_contents(tMesh *mesh)
{
  if(!mesh) return;

  if(PR) PRFs(":\n");

  free_mesh_contents_exceptMeshFuns(mesh);

  /* free skeleton in mesh */
  remove_all_MeshFuns(mesh);

  /* now set all in mesh back to 0 */
  memset(mesh, 0, sizeof(mesh[0]));
}

/* free mesh */
/*
void free_mesh(tMesh *mesh)
{
  free_mesh_contents(mesh);
  free(mesh);
}
*/

/**********************************************************************/
/* storage for lists of nodes */
/**********************************************************************/
/* allocate a node list with one node */
/* NOTE: we can also add to a nodelist that is NULL, so alloc_nodelist
   is not always needed */
tNlist *alloc_nodelist(tNode *node)
{
  tNlist *nlist;
  nlist = calloc(1, sizeof(*nlist));
  if(!nlist) errorexit("out of memory for nlist");
  nlist->node = node;
  return nlist;
}

/* add one node to nodelist after elem, and return new nodelist element
   that now contains the node */
tNlist *addnode_to_nodelist_after(tNlist *elem, tNode *node)
{
  tNlist *after = alloc_nodelist(node);
  return insertnodelist_into_nodelist_after(elem, after);
}
/* add one node to nodelist before elem, and return new nodelist element
   that now contains the node */
tNlist *addnode_to_nodelist_before(tNlist *elem, tNode *node)
{
  tNlist *before = alloc_nodelist(node);
  return insertnodelist_into_nodelist_before(elem, before);
}

/* make a copy */
tNlist *copy_of_nodelist(tNlist *elem)
{
  tNlist *dest = NULL;
  tNlist *src  = first_nodelist(elem);
  tNlist *el;

//printnodelist(src);
  for(el=src; el; el=el->next)
    dest = addnode_to_nodelist_after(dest, el->node);
//printnodelist(dest);

  return dest;
}

/* count num. of elem. in list */
int count_elements_nodelist(tNlist *list)
{
  tNlist *beg = first_nodelist(list);
  tNlist *el;
  int count=0;

  /* count elem. in list */
  for(el=beg; el; el=el->next) count++;

  return count;
}

/* replace element elem in nodelist by node */
tNlist *replacenode_in_nodelist(tNlist *elem, tNode *node)
{
  tNlist *repl = alloc_nodelist(node);
  return replace1_in_nodelist(elem, repl, 0);
}

/* insert nodelist "list" into another nodelist after elem,
   and return the end of "list" */
tNlist *insertnodelist_into_nodelist_after(tNlist *elem, tNlist *list)
{
  tNlist *elem2;
  tNlist *lend;
  tNlist *lbeg;

  /* if "list" is empty do nothing */
  if(!list) return elem;

  /* find end and beginning of tNlist *list */
  for(lend=list; lend->next; lend=lend->next) ;
  for(lbeg=list; lbeg->prev; lbeg=lbeg->prev) ;

  if(!elem) return lend;

  elem2 = elem->next;
  lend->next = elem2;
  lbeg->prev = elem;
  elem->next = lbeg;
  if(elem2) elem2->prev = lend;
  return lend;
}

/* insert nodelist "list" into another nodelist before elem,
   and return the first of "list" */
tNlist *insertnodelist_into_nodelist_before(tNlist *elem, tNlist *list)
{
  tNlist *elem2;
  tNlist *lend;
  tNlist *lbeg;

  /* if "list" is empty do nothing */
  if(!list) return elem;

  /* find end and beginning of tNlist *list */
  for(lend=list; lend->next; lend=lend->next) ;
  for(lbeg=list; lbeg->prev; lbeg=lbeg->prev) ;

  if(!elem) return lbeg;

  elem2 = elem->prev;
  lend->next = elem;
  lbeg->prev = elem2;
  elem->prev = lend;
  if(elem2) elem2->next = lbeg;
  return lbeg;
}

/* replace 1 element in a nodelist by a list and then free the element */
tNlist *replace1_in_nodelist(tNlist *elem, tNlist *list, int return_lend)
{
  tNlist *left;
  tNlist *right;
  tNlist *lend;
  tNlist *lbeg;

  /* when list is NULL do nothing */
  if(!list) return elem;

  /* find end and beginning of tNlist *list */
  for(lend=list; lend->next; lend=lend->next) ;
  for(lbeg=list; lbeg->prev; lbeg=lbeg->prev) ;

  if(!elem) return lbeg;

  left = elem->prev;
  right= elem->next;

  lend->next = right;
  lbeg->prev = left;
  if(right) right->prev = lend;
  if(left)  left->next = lbeg;

  free(elem);

  if(return_lend) return lend;
  else            return lbeg;
}

/* replace 1 element in a nodelist by a list, then free the element,
   return the very first element of the new list */
tNlist *first_replace1_in_nodelist(tNlist *elem, tNlist *list)
{
  tNlist *newlist = replace1_in_nodelist(elem, list, 0);
  return first_nodelist(newlist);
}

/* remove 1 element from nodelist, and
   return element after elem if return_next=1 */
tNlist *remove1_in_nodelist(tNlist *elem, int return_next)
{
  tNlist *left;
  tNlist *right;
  if(!elem) return 0;

  left = elem->prev;
  right= elem->next;
  if(right) right->prev = left;
  if(left)  left->next = right;

  free(elem);
  if(return_next) return right;
  else		  return left;
}

/* make a new list with all nodes that are children of the nodes in nlist */
tNlist *childnodelist_of_nodelist(tNlist *nlist)
{
  tNlist *elem, *cnlist, *clast;

  cnlist = clast = NULL; /* child node list is NULL at first */
  fornodelist(nlist, elem)
  {
    tNode *node = elem->node;
    tNode *child0 = node->child[0];
    if(child0)
    {
      int ijk;
      /* add all 8 children if child0 exists */
      for(ijk=0; ijk<8; ijk++)
      {
        clast = addnode_to_nodelist_after(clast, node->child[ijk]);
        if(!cnlist) cnlist = clast; /* save first entry of childlist */
      }
    }
  }
  return cnlist;
}

/* return 1st element in a nodelist */
tNlist *first_nodelist(tNlist *list)
{
  tNlist *lbeg;

  if(!list) return NULL;

  /* find beginning of tNlist *list */
  for(lbeg=list; lbeg->prev; lbeg=lbeg->prev) ;
  return lbeg;
}
/* return last element in a nodelist */
tNlist *last_nodelist(tNlist *list)
{
  tNlist *lend;

  if(!list) return NULL;

  /* find end of tNlist *list */
  for(lend=list; lend->next; lend=lend->next) ;
  return lend;
}


/* remove all from nodelist and free it */
void free_nodelist(tNlist *elem)
{
  tNlist *tmp;

  if(!elem) return;

  /* remove all after elem */
  for(tmp=elem->next; tmp; )
    tmp = remove1_in_nodelist(tmp, 1);

  /* remove all before elem */
  for(tmp=elem->prev; tmp; )
    tmp = remove1_in_nodelist(tmp, 0);

  /* remove elem */
  remove1_in_nodelist(elem, 0);
}

/* free all nodes in a list */
void free_nodesinlist(tNlist *elem)
{
  tNlist *tmp;

  /* free nodes in elem and all after in */
  for(tmp=elem; tmp; tmp=tmp->next)
    free_node(tmp->node);

  /* free nodes in all before elem */
  for(tmp=elem->prev; tmp; tmp=tmp->prev)
    free_node(tmp->node);
}

/**********************************************************************/
/* functions to allocate and free tMylnodes */
/**********************************************************************/
/* allocate or free categories in myln */
int realloc_myln_nncats(tMylnodes *myln, int nncats)
{
  int nncats_old = myln->nncats;
  int *ncat = myln->ncat;
  tNlist ***ln = myln->ln;
  int dc = nncats - nncats_old;
  int c;

  /* do nothing if old and new cat. numbers are the same */
  if(nncats_old == nncats) return nncats;

  /* free content of ln[c] for all c>=nncats */
  for(c=nncats; c<nncats_old; c++) free(ln[c]);

  if(nncats)
  {
    //Yo(1);
    myln->ncat = realloc(ncat, sizeof(ncat[0])*nncats);
    if(!myln->ncat) errorexit("no memory for myln->ncat");

    myln->ln = realloc(ln, sizeof(ln[0])*nncats);
    if(!myln->ln) errorexit("no memory for myln->ln");
  }
  else
  {
    //Yo(3);
    free(myln->ncat);
    free(myln->ln);
    myln->ncat = NULL;
    myln->ln = NULL;
    myln->nncats = 0;
    myln->nm = 0;
  }

  /* set new stuff to 0 */
  if(dc>0)
  {
    //Yo(4);
    memset(myln->ncat + nncats_old, 0, sizeof(ncat[0])*dc);
    memset(myln->ln  + nncats_old, 0, sizeof(ln[0])*dc);
  }
  myln->nncats = nncats;

  //if(nncats) printf("myln->ncat[0]=%d\n", myln->ncat[0]);
  return nncats;
}

/* make room for elements in category c, add cat. c if it's not existing yet */
int realloc_myln_ln_c(tMylnodes *myln, int c, int nelem)
{
  int ainc = 256;
  int chunks = 1 + nelem/ainc;
  tNlist ***ln = myln->ln;
  int nncats = myln->nncats;

  if(c>=nncats) realloc_myln_nncats(myln, c+1);

  myln->ln[c] = realloc(ln[c], sizeof(ln[0][0])*(chunks*ainc));
  if(!myln->ln[c]) errorexit("no memory for myln->ln[c]");

  return nelem;
}

/* add one element to myln in cat. c */
int addto_myln_ln_c(tMylnodes *myln, int c, tNlist *elem)
{
  int old_nelem, nelem;

  /* if myln is empty add one category */
  if(!myln->nncats) realloc_myln_nncats(myln, 1);

  /* make room for new el */
  old_nelem = myln->ncat[c];
  nelem = realloc_myln_ln_c(myln, c, old_nelem+1);
  /* reset nm in case now we have more */
  if(nelem > myln->nm) myln->nm = nelem;

  /* add elem at the end of ln[c] */
  myln->ln[c][old_nelem] = elem;
  myln->ncat[c] = nelem;
  //PRF;printf(": nelem=%d, myln->ncat[%d]=%d\n", nelem, c, myln->ncat[c]);
  return myln->ncat[c];
}

/* total number of nodes in all categories of myln */
int total_nnodes_in_myln(tMylnodes *myln)
{
  int c, sum=0;
  for(c=0; c < myln->nncats; c++)
    sum += myln->ncat[c];
  return sum;
}

/**********************************************************************/
/* functions to update elm->eploc->eid */
/**********************************************************************/
/* Update array of elms on this proc, set eids.
   Also update elm->dt and mesh->dt if auto_dt!=0 */
ulong update_elm_eid_dt(tMesh *mesh, double dt, int auto_dt,
                        double dtfac, double uniform_dtfac)
{
  int size = nMPI_size();
  int rank = nMPI_rank();
  int rk;
  struct list_head *pos;
  ulong eid;
  double dt_old = mesh->dt;
  double dt_new;
  if(auto_dt)
  {
    if(dt>0.) mesh->dt = dt;
    else      mesh->dt = DBL_MAX*0.1; /* reset mesh->dt to giant value */
  }

  /* set all my eids and Bcast my mesh->eidlim to each rank rk */
  eid = 0;
  for(rk=0; rk<size; rk++)
  {
    if(rk == rank)
    {
      /* set my dt info and all my eids */
      list_for_each(pos, &mesh->myelm_head)
      {
        tElm *elm = list_entry(pos, tElm, list);

        /* check if we need to change elm->dt and mesh->dt */
        if(auto_dt)
          adapt_node_dt_and_mesh_dt(elm, auto_dt, dtfac, uniform_dtfac);

        /* set my eid */
        elm->eploc->eid = eid++;
      }
      /* last elm->eploc->eid+1 is eidlim for me */
      mesh->eidlim[rk] = eid;
    }
    /* we use blocking MPI here */
    nMPI_Bcast(&(mesh->eidlim[rk]),1, nMPI_UNSIGNED_LONG, rk);
    /* This blocks until we get mesh->eidlim[rk] from rank rk. */

    /* update eid to start value for next rk iteration */
    eid = mesh->eidlim[rk];
  }

  //FIXME: should this be here???
  /* set elm array */
  alloc_and_set_mesh_myelm(mesh);

  /* if there are no nodes do not update dt mesh->dt */
  if(eid==0)
    mesh->dt = dt_old;

  /* now make sure we use the min dt of all ranks */
  dt_new = mesh->dt;
  nMPI_Allreduce(&dt_new, &(mesh->dt), 1, nMPI_DOUBLE, nMPI_MIN);

  if(mesh->dt != dt_old)
  { PRF;printf(": mesh->dt = %g\n", mesh->dt); }

  return eid;
}


/**********************************************************************/
/* functions to update the nodelist and node array in mesh */
/**********************************************************************/
/* Update array of leaf nodes on this proc, set nid.
   Also update node->dt and mesh->dt if auto_dt!=0 */
long update_mesh_myln_node_nid_dt(tMesh *mesh, double dt, int auto_dt,
                                  double dtfac, double uniform_dtfac)
{
  tNlist *elem;
  long nid = 0;
  //int lid = 0;

  /* delete mylns contents */
  realloc_myln_nncats(mesh->myln, 0);

  /* go over leaves if mesh->lns is not NULL */
  if(mesh->lns)
  {
    double dt_old = mesh->dt;
    if(auto_dt)
    {
      if(dt>0.) mesh->dt = dt;
      else      mesh->dt = DBL_MAX*0.1; /* reset mesh->dt to max value */
    }

    fornodelist(mesh->lns, elem)
    {
      tNode *node = elem->node;
      tNode *parent = node->parent;

      if(node->dat)
      {
        /* for now we put all leaves in cat. 0 */
        addto_myln_ln_c(mesh->myln, 0, elem);
        //PRF;printf(": myln->ncat[0]=%d %p\n", mesh->myln->ncat[0], elem);

        /* set lid and invalidate parent's lid */
        //node->lid = lid++;
        //if(parent) parent->lid = -lid;
      }
      /* set nid and invalidate parent's nid */
      Node_eid(node) = nid++;
      if(parent) parent->nid = -nid;
      //PRF;printf(": nmyln%ld nid%ld\n", nmyln,nid);

      /* set node MPI communicator */
      //i = nMPIvars_get_ncomms();
      //i = Node_eid(node) % i;
      //node->comm = nMPIvars_get_comm(i);
      // //PRF;printf(": i=%d node->comm=%d\n", i, node->comm);

      /* check if we need to change node->dt and mesh->dt */
      if(auto_dt)
        adapt_node_dt_and_mesh_dt(node, auto_dt, dtfac, uniform_dtfac);
    } /* end fornodelist */

    if(mesh->dt != dt_old) { PRF;printf(": mesh->dt = %g\n", mesh->dt); }
  }
  else /* mesh->lns is NULL, so free myln */
  {
    realloc_myln_nncats(mesh->myln, 0);
  }

  mesh->nln = nid;
  return nid;
}

/* update array of leaf nodes on this proc, set nid */
ulong update_mesh_myln_node_nid(tMesh *mesh)
{
  int Par_dt   = Par("dt");
  double dt    = Getd(Par_dt);
  /* auto_dt can be 0,1,2: */
  int auto_dt  = 1*Getv(Par_dt, "auto") + 2*Getv(Par_dt, "auto2");
  double dtfac = Getd(Par("dtfac"));
  double uniform_dtfac = Getd(Par("uniform_dtfac"));
  ulong ret;

  ret = update_elm_eid_dt(mesh, dt, auto_dt, dtfac, uniform_dtfac);

  //FIXME: remove this call
  ret = update_mesh_myln_node_nid_dt(mesh, dt, auto_dt,
                                     dtfac, uniform_dtfac);
  return ret;
}

/* return nid or -1 */
long get_node_nid(tNode *node)
{
  return node ? Node_eid(node) : -1;
}

/* return a local node id */
int calc_node_lid(tNode *node)
{
  tMesh *mesh = node->pat->mesh;
  long nnodes = mesh->nln;
  long size = nMPI_size();
  long npr2 = 2*nnodes/size + 1;
  long tmp = (Node_eid(node)) % npr2;
  int lid = tmp;

  return lid;
}

/* append a node list to mesh->lns and also update mesh->myln */
tNlist *append_nodelist_to_mesh_lns_myln(tMesh *mesh, tNlist *list)
{
  tNlist *lnl = NULL;
  if(mesh->lns)
  {
    lnl = last_nodelist(mesh->lns); /* last elem. in mesh->lns */
    lnl = insertnodelist_into_nodelist_after(lnl, list);
  }
  else
    mesh->lns = first_nodelist(list);

  /* update nids but leave node->dt alone */
  update_mesh_myln_node_nid_dt(mesh, -1., 0, 0.25, 0.125);

  return lnl;
}

/* replace elem in mesh->lns by nlist, return first of nlist */
tNlist *replace1_in_mesh_lns_myln(tNlist *elem, tNlist *nlist)
{
  tNlist *nlist_beg;
  tMesh *mesh = NULL;
  int update_lns;

  if(elem) mesh = elem->node->pat->mesh;
  else     errorexit("elem is NULL!!!");

  if(elem == mesh->lns) update_lns = 1;
  else                  update_lns = 0;

  nlist_beg = replace1_in_nodelist(elem, nlist, 0);
  if(update_lns) mesh->lns = nlist_beg;

  update_mesh_myln_node_nid(mesh);
  return nlist_beg;
}

/* replace current entry in leaf node list with its 8 new childern,
   return element with 0th child */
tNlist *make8children_in_mesh_lns_myln(tNlist *elem, int pt_typ[3], int n[3])
{
  tNode *parent;
  tNlist *children;
  tNlist *children0;

  TIMER_START;

  parent = elem->node;
  children = make8_child_nodes(parent, pt_typ, n);
  children0 = replace1_in_mesh_lns_myln(elem, children);

  TIMER_STOP;

  return children0;
}

/* replace siblings at element sib of mesh->lns by parent,
   node with parent is returned so we can use it later */
tNlist *remove8siblings_in_mesh_lns(tNlist *sib)
{
  tNode *parent, *node0;
  tNlist *elem, *elem0;
  tMesh *mesh = NULL;
  int update_lns;
  int ijk;

  if(sib==NULL) errorexit("sib is NULL!!!");
  mesh = sib->node->pat->mesh;
  parent = sib->node->parent;
  if(parent==NULL) errorexit("parent is NULL!!!");

  /* find sibling 0 */
  elem0=sib;
  for(ijk=sib->node->ijk; ijk>0; ijk--)
    elem0=elem0->prev;
  node0 = elem0->node;
  if(node0->parent != parent || node0->ijk != 0)
    errorexit("elem0 has wrong parent!");

  /* set elem to sibling 1 and remove the 7 after sibling 0 */
  elem=elem0->next;
  for(ijk=1; ijk<8; ijk++)
  {
    if(elem->node->parent != parent) errorexit("elem has wrong parent!");
    elem = remove1_in_nodelist(elem, 1);
  }

  /* do we need to update mesh->lns? */
  if(elem0 == mesh->lns) update_lns = 1;
  else                   update_lns = 0;

  /* replace sibling 0 by parent in list */
  elem0 = replacenode_in_nodelist(elem0, parent);

  /* reset mesh lists */
  if(update_lns) mesh->lns = elem0;

  return elem0;
}

/* replace siblings at element sib of mesh->lns by their parent, and then
   destroy the 8 siblings */
void destroy8siblings_in_mesh_lns_myln(tNlist *sib)
{
  tNlist *elem_parent;
  tNode *parent;

  TIMER_START;

  elem_parent = remove8siblings_in_mesh_lns(sib);
  parent = elem_parent->node;
  destroy_children(parent);
  update_mesh_myln_node_nid(parent->pat->mesh);

  TIMER_STOP;
}



/**********************************************************************/
/* storage for dat lists in the nodes */
/**********************************************************************/
/* allocate room for nv variables that can be enabled or disabled */
tDat *alloc_dat(tNode *node)
{
  int nv = node->pat->mesh->nvdb;
  tDat *dat;
  int f;

  dat = calloc(1, sizeof(tDat));
  if(!dat) errorexit("out of memory for dat");

  dat->node = node;
  dat->nv = nv;
  if(nv==0) return dat;

  dat->v = calloc(nv, sizeof(tArray *));
  if(!dat->v) errorexit("out of memory for dat->v");

  for(f=0; f<6; f++)
  {
    dat->s[f] = calloc(nv, sizeof(tSurface *));
    if(!dat->s[f]) errorexit("out of memory for dat->s[f]");

    /* create com[f] with double buffers,
       that will be freed by free_com in free_dat */
    dat->com[f] = alloc_com(sizeof(double), 1);
  }

  dat->ic = calloc(nv, sizeof(dat->ic[0]));
  if(!dat->ic) errorexit("out of memory for dat->ic");

  dat->icom = alloc_com(sizeof(double), 1);
  dat->gcom = alloc_com(sizeof(double), 1);

  return dat;
}
/* free dat and all arrays within it */
void free_dat(tDat *dat)
{
  int i,f;

  if(!dat) return;

  /* free contents */
  for(i=0; i<dat->nv; i++)
  {
    free_array(dat->v[i]);
    for(f=0; f<6; f++) free_surface(dat->s[f][i]);
    free_indc(dat->ic[i]);
  }
  free(dat->v);
  free(dat->ic);

  for(f=0; f<6; f++)
  {
    free(dat->s[f]);
    free_com(dat->com[f]);
  }

  free_com(dat->icom);
  free_com(dat->gcom);

  free(dat);
}

/* change dat->nv  to  dat->nv=nv_new */
void realloc_datvariables(tDat *dat, int nv_new)
{
  int i,f;

  if(nv_new<dat->nv) errorexit("implement var removal");

  if(nv_new)
  {
    dat->v = realloc(dat->v, nv_new*sizeof(tArray *));
    if(!dat->v) errorexit("out of memory for dat->v");

    dat->ic = realloc(dat->ic, nv_new*sizeof(dat->ic[0]));
    if(!dat->ic) errorexit("out of memory for dat->ic");
  }
  else
  {
    free(dat->v);
    free(dat->ic);
    return;
  }

  for(f=0; f<6; f++)
  {
    dat->s[f] = realloc(dat->s[f], nv_new*sizeof(tSurface *));
    if(!dat->s[f]) errorexit("out of memory for dat->s");
  }

  /* set newly added var pointers to NULL */
  for(i=dat->nv; i<nv_new; i++)
  {
    dat->v[i] = NULL;
    for(f=0; f<6; f++) dat->s[f][i] = NULL;
    dat->ic[i] = NULL;
  }
  dat->nv = nv_new;
}

/**********************************************************************/
/* storage for variable data base vdb in mesh */
/**********************************************************************/
void realloc_nodevariables(tNode *node, int nvdb_new)
{
  tDat *dat = node->dat;
  if(dat) realloc_datvariables(dat, nvdb_new);
}

void realloc_meshvariables(tMesh *mesh, int nvdb_new)
{
  int nvdb_old = mesh->nvdb;

  if(0) printf("realloc_meshvariables from %d to %d\n",
               mesh->nvdb, nvdb_new);

  /* realloc list on mesh struct */
  if(nvdb_new)
    mesh->vdb = realloc(mesh->vdb, sizeof(tVar)*(nvdb_new));
  else
    free(mesh->vdb);

  /* set newly added stuff to 0 */
  if(nvdb_new > nvdb_old)
  {
    tVar *newv0 = &(mesh->vdb[nvdb_old]);
    memset(newv0, 0, sizeof(tVar)*(nvdb_new-nvdb_old));
  }
  mesh->nvdb = nvdb_new;

  /* now make sure dat in nodes is also reallocated */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    realloc_nodevariables(node, nvdb_new);
  }
}

/**********************************************************************/
/* storage for variables in the nodes */
/**********************************************************************/

/* enable one component of variable i on one node */
void enablevarcomp_innode(tNode *node, int i)
{
  tDat *dat = node->dat;

  /* do nothing if no data is stored on this node */
  if(dat==NULL) return;

  if(i>=dat->nv) errorexiti("var comp %i does not exist", i);
  if(!dat->v[i])
  {
    tMesh *mesh = node->pat->mesh;
    int *ns = MeshVar_n_special(mesh, i);
    int Ne = MeshVar_Nextra(mesh, i);
    int dir, n[3];

    for(dir=0; dir<3; dir++)
    {
      int node_n = node->n[dir];

      /* use n_special of this var comp */
      switch(ns[dir])
      {
      case NOT_USED:
      case NODE_n:
        n[dir] = node_n;
        break;
      case NODE_nM1:
        n[dir] = node_n - 1;
        break;
      case NODE_nP1:
        n[dir] = node_n + 1;
        break;
      default:
        n[dir] = ns[dir];
      }
      if(n[dir]<=0) errorexiti("forbidden n_special in var comp %i", i);
    }
    dat->v[i] = alloc_array_extra(n, Ne);
    dat->nvenabled++;
    if(PR)
    {
      PRF;printf(": var_%d: %s\n", i, MeshVarName(node->pat->mesh, i));
    }
  }
}

/* disable one component of variable */
void disablevarcomp_innode(tNode *node, int i)
{
  tDat *dat = node->dat;

  /* do nothing if no data is stored on this node */
  if(dat==NULL) return;

  if(i>=dat->nv) errorexiti("var comp %i does not exist", i);
  if(dat->v[i])
  {
    free_array(dat->v[i]);
    dat->v[i] = NULL;
    dat->nvenabled--;
    if(PR)
    {
      PRF;printf(": var_%d: %s\n", i, MeshVarName(node->pat->mesh, i));
    }
  }
}

/* enable all components of a variable on one node */
void enablevar_innode(tNode *node, int i)
{
  tMesh *mesh = node->pat->mesh;
  int j, n = MeshVarNComponents(mesh, i);
  for(j=0; j<n; j++) enablevarcomp_innode(node, i+j);
}

/* disable all components of a variable on one node */
void disablevar_innode(tNode *node, int i)
{
  tMesh *mesh = node->pat->mesh;
  int j, n = MeshVarNComponents(mesh, i);
  for(j=0; j<n; j++) disablevarcomp_innode(node, i+j);
}

/* enable onr component of a variable on one pat */
void enablevarcomp_inpatch(tPat *pat, int i)
{
  tMesh *mesh = pat->mesh;

  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    if(node->pat == pat) enablevarcomp_innode(node, i);
  }
}

/* disable one component of a variable on one pat */
void disablevarcomp_inpatch(tPat *pat, int i)
{
  tMesh *mesh = pat->mesh;

  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    if(node->pat == pat) disablevarcomp_innode(node, i);
  }
}

/* enable all components of a variable on one pat */
void enablevar_inpatch(tPat *pat, int i)
{
  tMesh *mesh = pat->mesh;

  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    if(node->pat == pat) enablevar_innode(node, i);
  }
}

/* disable all components of a variable on one pat */
void disablevar_inpatch(tPat *pat, int i)
{
  tMesh *mesh = pat->mesh;

  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    if(node->pat == pat) disablevar_innode(node, i);
  }
}

/* enable all components of a variable on one mesh */
void enablevar(tMesh *mesh, int i)
{
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    enablevar_innode(node, i);
  }
}

/* disable all components of a variable on one mesh */
void disablevar(tMesh *mesh, int i)
{
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    disablevar_innode(node, i);
  }
}


/* enable variable list in a node */
void enablevarlist_innode(tNode *node, tVarList *vl)
{
  int i;
  if(vl) for(i=0; i<vl->n; i++) enablevarcomp_innode(node, vl->index[i]);
}

/* disable variable list in a node */
void disablevarlist_innode(tNode *node, tVarList *vl)
{
  int i;
  if(vl) for(i=0; i<vl->n; i++) disablevarcomp_innode(node, vl->index[i]);
}

/* enable variable list */
void enablevarlist(tVarList *vl)
{
  if(vl)
  {
    tMesh *mesh = vl->mesh;
    if(mesh)
    {
      formylnodes(mesh)
      {
        tNode *node = MyLnode;
        enablevarlist_innode(node, vl);
      }
    }
  }
}

/* disable variable list */
void disablevarlist(tVarList *vl)
{
  if(vl)
  {
    tMesh *mesh = vl->mesh;
    if(mesh)
    {
      formylnodes(mesh)
      {
        tNode *node = MyLnode;
        disablevarlist_innode(node, vl);
      }
    }
  }
}

// may need enablevarlist rather for patch or mesh ???




/* Bfaces stuff */
// /* add a bface to a pat, f denotes the pat face (0 to 5),
//    return index in bface list */
// int add_empty_bface(tPat *pat, int f)
// {
//   int fi = pat->nbfaces; /* add bface in this pos. in bface list */
//   void *ret;
//
//   /* increase size of bface list */
//   ret = realloc( pat->bface, (sizeof( *(pat->bface) ))*(fi+1) );
//   if(ret==NULL)  errorexit("add_bface: not enough memory for pat->bface");
//   pat->bface = ret;
//
//   /* mem for new bface */
//   ret = calloc( 1, sizeof( *(pat->bface[fi]) ) );
//   if(ret==NULL) errorexit("add_bface: not enough memory for pat->bface[n]");
//   /* add new bface */
//   pat->bface[fi] = ret;
//   pat->nbfaces = fi+1;
//
//   /* set some bface info */
//   pat->bface[fi]->mesh = pat->mesh;
//   pat->bface[fi]->b    = pat->b;
//   pat->bface[fi]->f    = f;
//   pat->bface[fi]->fi   = fi;
//   pat->bface[fi]->ob   = -1; /* other pat is not known yet */
//   pat->bface[fi]->ofi  = -1; /* fi in other pat is not known yet */
//   pat->bface[fi]->oXi  = -1; /* var indices of other coords not known yet */
//   pat->bface[fi]->oYi  = -1; /* var indices of other coords not known yet */
//   pat->bface[fi]->oZi  = -1; /* var indices of other coords not known yet */
//   return fi;
// }
//
// /* add a point ijk on face f to a bface with index fi, returns fi */
// /* if called with fi<0, it first calls add_empty_bface and returns the new fi */
// int add_point_to_bface_inpat(tPat *pat, int fi, int ijk, int f)
// {
//   tBface *bface;
//   /* make new bface when needed */
//   if(fi<0 || fi>=pat->nbfaces) fi = add_empty_bface(pat, f);
//   bface = pat->bface[fi];
//   /* make PointList id needed */
//   if(bface->fpts==NULL)
//     bface->fpts = AllocatePointList(pat->mesh);
//   /* add point ijk */
//   AddToPointList(bface->fpts, pat->b, ijk);
//   /* check if we have more than one face on this bface */
//   if(bface->f!=f) bface->f = -1;
//   return fi;
// }
//
// /* duplicate bface without pointlist fpts */
// tBface *duplicate_bface_without_fpts_for_mesh(tBface *bface0, tMesh *mesh)
// {
//   tBface *bface;
//   void *ret;
//   if(bface0==NULL) return NULL;
//
//   /* mem for new bface */
//   ret = calloc( 1, sizeof( *(bface) ) );
//   if(ret==NULL)
//     errorexit("duplicate_bface_without_fpts_for_mesh: not enough memory");
//   bface = ret;
//
//   /* make a shallow copy of the struct */
//   *bface = *bface0;
//   /* now set mesh pointer */
//   bface->mesh = mesh;
//   /* remove pointer to bface0->fpts */
//   bface->fpts = NULL;
//
//   return bface;
// }
//
// /* duplicate bface with pointlist fpts */
// tBface *duplicate_bface_for_mesh(tBface *bface0, tMesh *mesh)
// {
//   tBface *bface = duplicate_bface_without_fpts_for_mesh(bface0, mesh);
//   /* copy fpts */
//   bface->fpts = DuplicatePointList_for_mesh(bface0->fpts, mesh);
//   return bface;
// }
//
// /* free a bface */
// void free_bface(tBface *bface)
// {
//   if(bface!=NULL)
//   {
//     /* free the point lists of the faces */
//     FreePointList(bface->fpts);
//     free(bface);
//   }
// }
//
// /* free all bfaces in pat */
// void free_all_bfaces(tPat *pat)
// {
//   int i;
//   for(i=0; i<pat->nbfaces; i++)  free_bface(pat->bface[i]);
//   free(pat->bface);
//   pat->bface = NULL;
//   pat->nbfaces = 0;
// }
//
// /* remove a bface with index fi from a pat, return number of bfaces removed */
// int remove_bface(tPat *pat, int fi)
// {
//   void *ret;
//   int nbfaces = pat->nbfaces;
//   int i;
//
//   /* return 0 if bface does not exist */
//   if(fi<0 || fi>=nbfaces) return 0;
//
//   /* free the bface */
//   free_bface(pat->bface[fi]);
//
//   /* shift bfaces behind fi one position to the front */
//   for(i=fi; i<nbfaces-1; i++)
//   {
//     pat->bface[i] = pat->bface[i+1];
//     pat->bface[i]->fi = i;
//   }
//
//   /* reduce size of bface list */
//   nbfaces--;
//   ret = realloc( pat->bface, (sizeof( *(pat->bface) ))*(nbfaces) );
//   if(ret==NULL && nbfaces!=0)
//     errorexit("remove_bface: not enough memory for pat->bface");
//   pat->bface = ret;
//   pat->nbfaces = nbfaces;
//
//   return 1;
// }
//
// /* look for empty bfaces and remove them */
// int remove_bfaces_with_NULL_fpts(tPat *pat)
// {
//   int n, fi;
//   n=0;
//   for(fi=0; fi<pat->nbfaces; fi++)
//     if(pat->bface[fi]->fpts == NULL)
//     {
//       int r = remove_bface(pat, fi); /* this decreases pat->nbfaces */
//       fi=fi-r;  /* go back by one in fi so that for-loop covers fi again */
//       n++;      /* count number of bfaces removed */
//     }
//   return n;
// }


