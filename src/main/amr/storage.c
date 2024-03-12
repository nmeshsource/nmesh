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

/* alloc mem for range part of array */
void alloc_2darray_irange_of_j(tArray *array)
{
  int k;

  for(k=0; k<2; k++)
  {
    array->range[k] = calloc(array->n[1], sizeof(array->range[k][0]));
    if(!array->range[k]) errorexit("out of memory");
  }
}

/* free an array */
void free_array(tArray *array)
{
  int k;

  if(!array) return;

  if( (array->si == 0) && (array->d_nofree == 0) )
    free(array->d); /* free data only if this segm. 0 and it should be freed */

  for(k=0; k<2; k++)
    free(array->range[k]); /* free range */

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

/* redim array in eploc size increments */
int redim_array_Neplocs(tArray *ar, int Neplocs)
{
  size_t se = sizeof(ar->eploc[0]);  //sizeof eploc
  size_t nbytes = se*Neplocs;        //number of bytes needed
  size_t sd = sizeof(ar->d[0]);      //sizeof double
  size_t nd = (nbytes+sd-1)/sd;      //num. of doubles in nbytes
  return redim_array(ar, nd,1,1);
}


/****************************************************************************/
/* functions to select interp. scheme for mesh refinement */
/****************************************************************************/

/* Set order and scheme from elm info and amr pars. */
void amr_interp_get_order_scheme(tElm *elm, int *order, int *scheme)
{
  tMesh *mesh = Elm_mesh(elm);
  int *pt_typ = elm->pt_typ;
  int is_UNI, schm, npts;
  int d;

  /* check grid */
  is_UNI=1;
  for(d=0; d<3; d++)
    if(pt_typ[d]!=P_UNIFORM) { is_UNI=0; break; }

  /* pick scheme */
  if(is_UNI) schm = INTERP_WENO;
  else       schm = INTERP_LAGRANGE;

  /* BUT amr->force_interp_scheme can override scheme */
  if(amr->force_interp_scheme) schm = amr->force_interp_scheme;

  /* set npts */
  switch(schm)
  {
  case INTERP_LAGRANGE:
    npts = Geti(amr->Lagrange_interp_order);
    break;
  case INTERP_WENO:
    npts = Geti(amr->WENO_interp_order);
    break;
  default:
    errorexit("unknown scheme");
  }

  /* set order, scheme */
  *order  = npts;
  *scheme = schm;
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

/* alloc a new elm using tEploc info */
tElm *alloc_elm_of_eploc(tMesh *mesh, tEploc *eploc)
{
  tElm0 elm0[1] = {0};
  amr_init_elm0_from_eploc(mesh, eploc, elm0);
  return alloc_elm_of_elmheader(mesh, elm0);
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
ulong alloc_and_set_mesh_myelm(tMesh *mesh)
{
  ulong ei;
  struct list_head *pos;

  /* free whatever we had so far in myelm */
  free(mesh->myelm);
  mesh->myelm = NULL;

  /* get number of elms on this rank */
  mesh->nmyelm = list_count_nodes(&mesh->myelm_head);
  //PRF;printf(": mesh->nmyelm=%ld\n", mesh->nmyelm);
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
void make_and_add_root_elm(tPat *pat, int n[3], int pt_typ[3], int datrank)
{
  int rank = nMPI_rank();
  tMesh *mesh = pat->mesh;
  tElm0 elm0[1] = {0}; /* elm header with root node info */
  tEploc *eploc = elm0->eploc;
  int i;

  /* check for overflow in eploc->p */
  i = (1 << (sizeof(eploc->p)*8)) - 1;
  if(pat->p > i) errorexiti("cannot have more than %d patches", i);

  /* fill in info */
  eploc->p = pat->p;
  eploc->l = 0; /* root node */
  eploc->ploc[0] = 0;
  amr_set_elm0_bbox(mesh, elm0);

  /* save n and pt_typ for root node */
  for(i=0; i<3; i++)
  {
    elm0->n[i] = n[i];
    elm0->pt_typ[i] = pt_typ[i];
  }
  elm0->np = n[0] * n[1] * n[2];
  elm0->eploc->eid = EID_INVALID;    /* mark eid as not set */

  /* set where dat needs to be allocated */
  elm0->datrank = datrank;

  /* make root element only on the MPI rank that owns it */
  if(datrank==rank)
  {
    tElm *elm = alloc_elm_of_elmheader(mesh, elm0);

    /* dat needs to be allocated */
    elm->dat = alloc_dat(elm);

    /* add new root element to list mesh->myelm_head */
    list_add_tail(&elm->list, &mesh->myelm_head);
  }
  else if((datrank<0) && (rank==0))
  {
    /* if datrank<0 create an elm without dat on rank0 */
    tElm *elm = alloc_elm_of_elmheader(mesh, elm0);
    list_add_tail(&elm->list, &mesh->myelm_head);
  }

  ///* set rnode info */
  //pat->rnode[0] = elm0[0];
}


/* make a child node element */
tElm *make_child_elm(tElm *parent, int n[3], int pt_typ[3], int ijk)
{
  tMesh *mesh = parent->pat->mesh;
  tElm *elm = alloc_elm(mesh);
  int d, vi,nvdb, f;

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
    int order, scheme;

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
    amr_interp_get_order_scheme(parent, &order, &scheme);
    for(vi=0; vi<nvdb; vi++)
      if(parent->dat->v[vi])
      {
        int vt = MeshVarType(mesh, vi);
        /* enable same vars in this dat as in parent->dat */
        enablevarcomp_innode(elm, vi);

        /* fill elm->dat with interpolation data from parent */
        if( (vt==EVOVAR) || (vt==DATAVAR) ) /* exclude Aux. vars */
        {
          interp_topoints(parent, parent->dat->v[vi], Xp,
                          order, scheme, 1., elm->dat->v[vi]);
        }
      } /* end: if parent has dat->v[vi] */
    free_array(Xp[2]);
    free_array(Xp[1]);
    free_array(Xp[0]);

    /* init coords in this new elm */
    coordinates_init_node(elm);

    /* mark nbinfo as not set */
    disablevar_innode(elm, amr->elm_nbinfo0);
    for(f=0; f<6; f++) elm->dat->info->nnbinfo[f]=-1;
  }
  return elm;
}

/* make 8 children, insert them into mesh->myelm_head, and return child0 */
tElm *replace_parent_by_8children(tElm *parent, int n[3], int pt_typ[3])
{
  //tMesh *mesh = parent->pat->mesh;
  struct list_head clist;
  tElm *child, *child0=NULL; //set child0=NULL to avoid stupid gcc warning
  int ijk;

  INIT_LIST_HEAD(&clist);

  /* make children */
  for(ijk=0; ijk<8; ijk++)
  {
    child = make_child_elm(parent, n, pt_typ, ijk);
    list_add_tail(&child->list, &clist);
    if(ijk==0) child0 = child; /* save first child */
  }

  /* #pragma omp critical (change_mesh_myelm_list) */
  /* NOTE: For some reason gcc's -fsanitize=thread throws a ?false? positive
           if I use a named critical section!
           So replace "GEN_Pragma(omp critical (change_mesh_myelm_list))"
           by "GEN_Pragma(omp critical)" when debugging races!!! */
  //GEN_Pragma(omp critical (change_mesh_myelm_list))
  GEN_Pragma(omp critical)
  {
    /* NOTE: The new children have all zero for nfnb, fnb, and nnbinfo=-1.
             Also, all their neighbors have now the wrong nfnb and nbinfo.
             Even worse, all its neighbors have fnb pointers pointing
             to the parent which will be removed!!!  */
    /* now replace parent by clist in mesh->myelm_head */
    list_splice(&clist, &parent->list);
    list_del(&parent->list);
  }

  return child0;
}

/* use info in child0 to recreate the parent */
tElm *make_parent_elm(tElm *child0, int n[3], int pt_typ[3])
{
  tMesh *mesh = Elm_mesh(child0);
  tElm *parent = alloc_elm(mesh);
  int d;
  struct list_head *pos_ijk;
  int ijk;

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
    int vi, f;
    tArray *Xp[3], *Ip[8], *Xc[8][3], *Res[8];

    if(!parent->dat) parent->dat = alloc_dat(parent);

    /* enable same vars in parent->dat as in child0->dat */
    for(vi=0; vi<nvdb; vi++)
      if(child0->dat->v[vi]) enablevarcomp_innode(parent, vi);

    /* but disable nb-info vars, mark them as not set */
    disablevar_innode(parent, amr->elm_nbinfo0);
    for(f=0; f<6; f++) parent->dat->info->nnbinfo[f]=-1;

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
            {
              int order, scheme;
              amr_interp_get_order_scheme(child, &order, &scheme);
              interp_toIpoints(child, var, Xc[ijk],Ip[ijk],
                               order, scheme, 1., Res[ijk]);
            }
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

  return parent;
}

/* Remove 8 children and replace them by parent.
   We assume that all 8 have been moved to this rank before
   replace_8localchildren_by_parent is called!
   In: child0,n,pt_typ
   Out: ch_head <- list of removed children, (initialized by caller)
   Returns: parent */
tElm *replace_8localchildren_by_parent(tElm *child0, int n[3], int pt_typ[3],
                                       struct list_head *ch_head)
{
  tElm *parent;
  unsigned char l_ch0 = child0->eploc->l;
  struct list_head *pos_ijk;
  int ijk;

  /* sanity checks */
  pos_ijk = &child0->list;
  for(ijk=0; ijk<8; ijk++)
  {
    tElm *child = list_entry(pos_ijk, tElm, list);
    tEloc eloc[1];

    if(!child->dat)
      errorexit("all 8 children need to be on this proc");
    if(child->eploc->l != l_ch0)
      errorexit("all 8 starting with child0 must be on same level");

    eloc_from_eploc(eloc, child->eploc);
    if(eloc->loc[l_ch0-1] != '0' + ijk)
      errorexiti("this is not child%d", ijk);

    pos_ijk = pos_ijk->next;  /* pos of next child */
  }

  /* make new parent elm */
  parent = make_parent_elm(child0, n, pt_typ);

  /* replace children by parent mesh->myelm_head list */
  /* #pragma omp critical (change_mesh_myelm_list) */
  /* NOTE: For some reason gcc's -fsanitize=thread throws a ?false? positive
           if I use a named critical section!
           So replace "GEN_Pragma(omp critical (change_mesh_myelm_list))"
           by "GEN_Pragma(omp critical)" when debugging races!!! */
  //GEN_Pragma(omp critical (change_mesh_myelm_list))
  GEN_Pragma(omp critical)
  {
    /* NOTE: This new parent has all zero for nfnb, fnb, and nnbinfo=-1.
             Also, all its neighbors have now the wrong nfnb and nbinfo.
             Even worse, all its neighbors have fnb pointers pointing
             to the children which will be removed!!!  */

    /* FIXME: We should go over children's nbs and set whatever
              nb-info we can!!! */

    /* now replace children by parent in mesh->myelm_head */
    /* first insert parent before child0 */
    list_add_tail(&parent->list, &child0->list);
    /* remove children from mesh list, but add them to ch_head */
    for(ijk=0; ijk<8; ijk++)
    {
      pos_ijk = (parent->list).next; //pos after parent
      list_del(pos_ijk);
      list_add_tail(pos_ijk, ch_head);
    }
  }

  return parent;
}


/**************************************************************************/
/* update n or pt_typ in a node or elm */
/**************************************************************************/

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
    int order, scheme;

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
    amr_interp_get_order_scheme(node_old, &order, &scheme);
    for(vi=0; vi<nvdb; vi++)
      if(node_old->dat->v[vi])
      {
        int vt = MeshVarType(mesh, vi);
        /* enable same vars in new dat as in dat_old */
        enablevarcomp_innode(node, vi);

        /* fill node->dat with interpolation data from old dat */
        if( (vt==EVOVAR) || (vt==DATAVAR) ) /* exclude Aux. vars */
        {
          interp_topoints(node_old, node_old->dat->v[vi], Xp,
                          order, scheme, 1., node->dat->v[vi]);
        }
        /* copy nbinfo vars */
        if( (vi >= amr->elm_nbinfo0) && (vi < amr->elm_nbinfo0+6) )
        {
          tArray *nbinfo_old = node_old->dat->v[vi];
          tArray *nbinfo     = node->dat->v[vi];
          int Neplocs = array_Neplocs(nbinfo_old);
          redim_array_Neplocs(nbinfo, Neplocs);
          copy_array_data(nbinfo_old, nbinfo);
        }
      } /* end: if node_old has dat->v[vi] */
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
  struct list_head *pos, *sav;

  if(!pat) return;

  //PRFs(":\n");

  /* all nb info will be bad if a patch is removed, so delete it */
  amr_remove_mesh_nbelm(mesh, 0);

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

  /* update myln lists */
  update_mesh_myelms_elm_eid_dt(mesh);

  /* free all in CI coordinfo, and also all bfaces  */
  free_pat_CI(pat);
  remove_all_bfaces(pat);

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
  //MUTEX_INIT(mesh->mutex);

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

  ///* node list in mesh */
  //realloc_myln_nncats(mesh->myln, 0);

  /* free mesh->myelm */
  free_mesh_myelm(mesh);

  /* set patch and node stuff to 0 */
  mesh->npats = 0;
  mesh->pat = NULL;
  //memset(&(mesh->myln[0]), 0, sizeof(mesh->myln[0]));
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
  //MUTEX_DESTROY(mesh->mutex);

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
/* functions to allocate and free tMylnodes */
/**********************************************************************/
/* allocate or free categories in myln */
int realloc_myln_nncats(tMylnodes *myln, int nncats)
{
  int nncats_old = myln->nncats;
  int *ncat = myln->ncat;
  tElm ***ln = myln->ln;
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
  tElm ***ln = myln->ln;
  int nncats = myln->nncats;

  if(c>=nncats) realloc_myln_nncats(myln, c+1);

  myln->ln[c] = realloc(ln[c], sizeof(ln[0][0])*(chunks*ainc));
  if(!myln->ln[c]) errorexit("no memory for myln->ln[c]");

  return nelem;
}

/* add one element to myln in cat. c */
int addto_myln_ln_c(tMylnodes *myln, int c, tElm *elm)
{
  int old_nelm, nelm;

  /* if myln is empty add one category */
  if(!myln->nncats) realloc_myln_nncats(myln, 1);

  /* make room for new el */
  old_nelm = myln->ncat[c];
  nelm = realloc_myln_ln_c(myln, c, old_nelm+1);
  /* reset nm in case now we have more */
  if(nelm > myln->nm) myln->nm = nelm;

  /* add elem at the end of ln[c] */
  myln->ln[c][old_nelm] = elm;
  myln->ncat[c] = nelm;
  //PRF;printf(": nelm=%d, myln->ncat[%d]=%d\n", nelm, c, myln->ncat[c]);
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
/* functions to update elm->eploc->eid and mesh->eidlim */
/**********************************************************************/
/* Update elm eids and set mesh->eidlim.
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

        ///* save old eid */
        //elm->oid = elm->eploc->eid;

        /* set my new eid */
        elm->eploc->eid = eid++;
      }
      /* last elm->eploc->eid+1 is eidlim for me */
      mesh->eidlim[rk] = eid;
    }
    /* we use blocking MPI here */
    MCK( nMPI_Bcast(&(mesh->eidlim[rk]),1, nMPI_UNSIGNED_LONG, rk) );
    /* This blocks until we get mesh->eidlim[rk] from rank rk. */

    /* update eid to start value for next rk iteration */
    eid = mesh->eidlim[rk];
  }
  /* if there are no elms do not update dt mesh->dt */
  if(eid==0)
    mesh->dt = dt_old;

  //for(rk=0; rk<size; rk++)
  //{ PRF;printf(": rank%d mesh->eidlim[%d]=%lu\n", rank, rk, mesh->eidlim[rk]); }

  /* now make sure we use the min dt of all ranks */
  dt_new = mesh->dt;
  MCK( nMPI_Allreduce(&dt_new, &(mesh->dt), 1, nMPI_DOUBLE, nMPI_MIN) );

  if(mesh->dt != dt_old)
  { PRF;printf(": mesh->dt = %g\n", mesh->dt); }

  return eid;
}

/* update array of leaf nodes or elms on this proc, set eid and dt */
ulong update_mesh_myelms_elm_eid_dt(tMesh *mesh)
{
  int Par_dt   = Par("dt");
  double dt    = Getd(Par_dt);
  /* auto_dt can be 0,1,2: */
  int auto_dt  = 1*Getv(Par_dt, "auto") + 2*Getv(Par_dt, "auto2");
  double dtfac = Getd(Par("dtfac"));
  double uniform_dtfac = Getd(Par("uniform_dtfac"));
  ulong ret;

  /* set elm array */
  alloc_and_set_mesh_myelm(mesh);

  /* update elm eids and mesh->eidlim */
  ret = update_elm_eid_dt(mesh, dt, auto_dt, dtfac, uniform_dtfac);

  return ret;
}

/**********************************************************************/
/* functions to calculate node IDs */
/**********************************************************************/

/* return a local node id */
int calc_node_lid(tNode *node)
{
  tMesh *mesh = node->pat->mesh;
  ulong size = nMPI_size();
  ulong nnodes = mesh->eidlim[size-1];
  ulong npr2 = 2*nnodes/size + 1;
  ulong tmp = (Node_eid(node)) % npr2;
  int lid = tmp;

  return lid;
}

/* return a local elm id */
ulong calc_elm_lid(tNode *elm)
{
  tMesh *mesh = elm->pat->mesh;
  ulong size = nMPI_size();
  ulong nelms = mesh->eidlim[size-1];
  ulong npr2 = 2*nelms/size + 1;
  ulong tmp = (Node_eid(elm)) % npr2;
  ulong lid = tmp;

  return lid;
}

/* return another local elm id */
ulong calc_local_elm_id(tNode *elm)
{
  tMesh *mesh = elm->pat->mesh;
  ulong nmyelms = mesh->nmyelm;
  return (Node_eid(elm)) % nmyelms;
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
  for(f=0; f<6; f++) dat->info->nnbinfo[f]=-1; //amr_elm_nbinfo is invalid
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


