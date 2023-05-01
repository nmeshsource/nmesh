/* connections.c */
/* Wolfgang Tichy, 3/2023 */

#include "nmesh.h"
#include "amr.h"




/****************************************************************************/
/* functions that determine the order of element locations */
/****************************************************************************/

/* Function that orders locations described in in tEloc:
   return -1,0,1 if loc is before,at,after eloc
   This can be use in qsort. */
int loccmp(const void *loc, const void *eloc)
{
  const tEloc *lc = (const tEloc *) loc;
  const tEloc *el = (const tEloc *) eloc;
  int i;

  /* if not in same patch p move right or left in search */
  if(lc->p > el->p) return  1; /* after */
  if(lc->p < el->p) return -1; /* before */

  /* ok, if we get here, lc and el are in same patch */
  if(el->l >= lc->l)
  {
    for(i=0; i<lc->l; i++)
    {
      if(lc->loc[i] == el->loc[i]) continue;
      if(lc->loc[i] >  el->loc[i]) return  1;
      else                         return -1;
    }
    return 0; /* lc and el are equal up the first lc->l */
  }
  else
  {
    for(i=0; i<el->l; i++)
    {
      if(lc->loc[i] == el->loc[i]) continue;
      if(lc->loc[i] >  el->loc[i]) return  1;
      else                         return -1;
    }
    return 1; /* make binarysearch move to right */
  }
}

/* return -1,0,1 if loc is before,at,after elem location,
   this is used in binarysearch */
int lecmp(const void *loc, const void *elem, void *arg)
{
  const tEloc *lc = (const tEloc *) loc;
  //const tElm **elm_arr = elem; //(const tElm **) elem;
  const tElm *const*elm_arr = elem; //(const tElm **) elem;
  const tEloc *elc = elm_arr[0]->eloc;
  int cmp = loccmp(lc, elc);
  PRFs(": ");printeloc_s(lc, " ");printeloc_s(elc, " ");
  printf("--> cmp=%d\n", cmp);
  //printelm(elm_arr[0]);
  return cmp;
}

/****************************************************************************/
/* primitive functions that work on integers and strings */
/****************************************************************************/

/* find ijk from l,loc by reading last in loc */
int connections_get_ijk(int l, const char loc[LOCSMAX])
{
  if(l<1) return 0;
  return loc[l-1] - '0';
}

/* return 1 if ijk face touches parentnode-face */
int connections_ijk_is_at_parentface(int ijk, int face)
{
  int ns[] = {2,2,2};
  /* set node's i,j,k */
  int k = kOfInd_n(ijk, ns);
  int j = jOfInd_n_k(ijk, ns,k);
  int i = iOfInd_n_jk(ijk, ns,j,k);

  switch(face)
  {
  case 0:  return i^1;
  case 1:  return i;
  case 2:  return j^1;
  case 3:  return j;
  case 4:  return k^1;
  case 5:  return k;
  default: errorexit("face must be 0,1,2,3,4,5");
  }
  return -1; /* cannot get here */
}

/* return ijk of nb in same level, assuming ijk is not on patch face */
int connections_get_inner_nb_ijk(int ijk, int dir)
{
  int ns[] = {2,2,2};
  int i,j,k;  /* index of loc */

  //if(l<1) return 0;

  /* get i,j,k */
  k = kOfInd_n(ijk, ns);
  j = jOfInd_n_k(ijk, ns,k);
  i = iOfInd_n_jk(ijk, ns,j,k);

  /* i, j, or k to nb  value */
  switch(dir)
  {
  case 0:  return Ind_n(i^1,j,k, ns);
  case 1:  return Ind_n(i,j^1,k, ns);
  case 2:  return Ind_n(i,j,k^1, ns);
  default: errorexit("dir must be 0,1,2");
  }
  return -1;
}


/* Out: return value: number of faces l,loc is on
        patface[f] = 1 if l,loc is on patch face f */
int connections_loc_on_patchface(int l, const char loc[LOCSMAX],
                                 int patface[6])
{
  int ll, f, npatfaces;

  /* set result if l=0 */
  npatfaces = 6;
  for(f=0; f<6; f++) patface[f] = 1;

  if(l<1) return npatfaces; // <--- not needed

  for(ll=1; ll<=l; ll++)
  {
    int ijk = connections_get_ijk(ll, loc);

    for(f=0; f<6; f++)
    {
      if(!connections_ijk_is_at_parentface(ijk, f))
      {
        npatfaces--;
        patface[f] = 0;
      }
      if(npatfaces<=0) break;
    }
    if(npatfaces<=0) break;
  }
  return npatfaces;
}


/* get nbloc,nf of neighbor on face of elm with l,loc,face
   BUT this func works only if the face is not on patch face.
   In:  l, loc, face
   Out: nbloc, nb_f */
int connections_get_nbloc_SameLevel_InsidePat(int l, const char loc[LOCSMAX],
                                              int face,
                                              char nbloc[LOCSMAX], int *nb_f)
{
  int patface[6];
  int ijk, nb_ijk;

  PRF;printf(": l=%d loc=%s face=%d\n", l, loc, face);

  //nfaces = connections_loc_on_patchface(l,loc, patface);
  connections_loc_on_patchface(l,loc, patface);

  if(patface[face])
    errorexiti("face%d of loc is on patch surface", face);

  /* find ijk of node and ijk of nb */
  ijk = connections_get_ijk(l, loc);
  nb_ijk = connections_get_inner_nb_ijk(ijk, face/2);
  *nb_f  = face^1;
  printf("  ijk=%d nb_ijk=%d *nb_f=%d\n", ijk, nb_ijk, *nb_f);

  if(connections_ijk_is_at_parentface(ijk, face))
  {
    char pnbloc[LOCSMAX]; /* location of parent nb */
    PRF;printf(" at parentface\n");

    /* l-1,loc is parent */
    connections_get_nbloc_SameLevel_InsidePat(l-1,loc, face, pnbloc, nb_f);
    pnbloc[l] = 0; /* add string-end marker */
    strncpy(nbloc, pnbloc, LOCSMAX);
    nbloc[l-1] = nb_ijk;
    if(l<LOCSMAX) nbloc[l] = 0;
    return l;
  }
  else
  {
    strncpy(nbloc, loc, LOCSMAX);
    nbloc[l-1] = nb_ijk + '0';
    return l;
  }
  return l;
}

/* get nbloc,nf of neighbor on face of elm with l,loc,face
   In:  elm, face
   Out: nbeloc,nbface */
int connections_get_nb_eloc_face(const tElm *elm, int elmface,
                                 tEloc nbeloc[1], int *nbface)
{
  const tEloc *eloc = elm->eloc;
  int patface[6]; //, nfaces;
  int l;

  PRFs(": ");printeloc(eloc);printf(" elmface=%d\n", elmface);

  connections_loc_on_patchface(eloc->l, eloc->loc, patface);
  /* deal with more complicated case where elmface is on patch surface */
  if(patface[elmface])
  {
    tPat *pat = elm->pat;

    errorexit("deal with pat face");
    /*
    tBface *bfaces = pat->bfaces[face];
    // loop over bfaces
    forbfacesonface(pat, f, bface) ;
    // same as: for(bface=bfaces; bface; bface=bface->next) ;

    //if(bfaces && bfaces->boundary==OUTERBOUND)
    */

    tBface *bface;
    tNlist *nbl, *nblist1, *elem;
    tEloc nbeloc[1];
    int nc, nb_f;
    int nnb = 0;   /* number of nfaces added */


    /* loop over all bfaces on face and find nb */
    forbfacesonface(pat, elmface, bface)
    {
      tBface *obface = bface->obface;
      int touch;

      /* do nothing if no other patch face */
      if(!obface) continue;

      /* eloc and face of root elm in other patch */
      nbeloc->l      = 0;
      nbeloc->loc[0] = 0;
      nb_f = obface->f;

      //dummies:
      tNode *node = elm;
      tNode *nb   = obface->pat->rnode;



      /* so now we have a neighbor loc, but is it childless? */
      nc = count_children(nb);
      if(nc==0) /* neighbor has 0 children */
      {
        nblist1 = alloc_nodelist(nb);
      }
      else
      {
        if(nc!=8) errorexiti("nb has %d children, not 8!!!", nc);

        /* find nblist1 with all leaves on face nb_f */
        nblist1 = leafdescendants_along_face(nb, nb_f, NULL);
      }




      /* beginning of nblist1 */
      nbl = first_nodelist(nblist1);

      /* go over nbl and remove all who do not have common face points
         with the node */
      nblist1 = NULL;
      fornodelist(nbl, elem)
      {
      nbl_loop_start:

        /* get neigh. and check if node and nb have common points */
        nb = elem->node;
        touch = common_facepoints(elm,elmface, nb,nb_f);
        if(touch)
        {
          nblist1 = elem; /* save elem that touches our node */
          continue;
        }


        /* remove nb=elem->node from nbl */
        elem = remove1_in_nodelist(elem, 1); /* now elem has the next one */
        if(elem) goto nbl_loop_start;
        else     break;
      }

      /* rewind nblist1 so that the fornodelist loop below works */
      nblist1 = first_nodelist(nblist1);

      /* add all in nblist1 as nfaces */
      fornodelist(nblist1, elem)
      {
        nb = elem->node;
        add_nface(node, elmface, nb, nb_f);
        nnb++; /* count neighbors */
      }

      /* free node lists */
      free_nodelist(nblist1);
    }














    /*
    tPat *pat = elm->pat;
    tBface *bfaces = pat->bfaces[face];
    // loop over bfaces
    forbfacesonface(pat, f, bface) ;
    // same as: for(bface=bfaces; bface; bface=bface->next) ;

    //if(bfaces && bfaces->boundary==OUTERBOUND)
    */
    l=-9999;
  }
  else /* elmface is a refinement boundary in patch interior */
  {
    nbeloc->p = eloc->p;
    nbeloc->l = eloc->l;
    l = connections_get_nbloc_SameLevel_InsidePat(eloc->l, eloc->loc, elmface,
                                                  nbeloc->loc, nbface);
  }

  return l;
}






/****************************************************************************/
/* functions that work on eloc */
/****************************************************************************/

//void

/****************************************************************************/
/* functions that work on elm */
/****************************************************************************/

/* get ijk of elm */
int elm_get_ijk(tElm *elm)
{
  tEloc *eloc = elm->eloc;
  return connections_get_ijk(eloc->l, eloc->loc);
}


/****************************************************************************/
/* functions to initialize tElm */
/****************************************************************************/

/* find patch of elm and save it in elm->pat */
void amr_set_elm_pat(tMesh *mesh, tElm *elm)
{
  tEloc *eloc = elm->eloc;
  int p = eloc->p;
  //elm->pat = elm->mesh->pat[p];
  elm->pat = mesh->pat[p];
}

/* find bbox of elm and save it in elm->bbox */
void amr_set_elm_bbox(tElm *elm)
{
  tPat *pat = elm->pat;
  tEloc *eloc = elm->eloc;
  int l = eloc->l; /* get level number */
  char *loc = eloc->loc;
  double *bbox  = elm->bbox;
  double LX[3];
  int f, d, ll;

  /* copy patch bbox values, and put lengths into LX */
  for(f=0; f<6; f++) bbox[f] = pat->bbox[f];
  for(d=0; d<3; d++) LX[d]   = bbox[2*d+1] - bbox[2*d];

  /* cut pat->bbox values into 2 for each level */
  for(ll=1; ll<=l; ll++)
  {
    int ijk = connections_get_ijk(ll, loc);
    for(d=0; d<3; d++)
    {
      int bt = 1 << d; /* bit we check: bt = 1, 2, or 4 */
      LX[d] /= 2.;
      if(ijk & bt) bbox[2*d] = bbox[2*d+1] - LX[d];
      else         bbox[2*d+1] = bbox[2*d] + LX[d];
    }
  }
}

/* set eloc of child */
int amr_set_child_eloc(tEloc *parentloc, int ijk, tEloc *eloc)
{
  int l   = parentloc->l;
  if(l >= LOCSMAX-1)
    errorexit("parentloc is at limit ==> no further child possible!");
  eloc->p = parentloc->p;
  eloc->l = l + 1;
  strncpy(eloc->loc, parentloc->loc, LOCSMAX);
  eloc->loc[l]   = '0' + ijk;
  eloc->loc[l+1] = 0;
  return l+1;
}

/****************************************************************************/


//////////////////////////////////////////////////////////////////////////
// replaces l_XYZ_of_xyz
//FIXME: pick a good file for this func
//       maybe around l_XYZ_of_xyz in main/coordinates/get_coords.c
/* set X and return 1 if x is inside this elm, otherwise return 0 */
int elmXYZ_of_xyz(tElm *elm, int ind, double X[3], const double x[3])
{
  tPat *pat = elm->pat;
  int d, stat=0;

  /* get X */
  if(pat->XYZ_of_xyz)
    //stat = pat->XYZ_of_xyz(pat, (tNode *)elm,ind, X, x);
    stat = pat->XYZ_of_xyz(pat, elm,ind, X, x);
  else
    for(d=0; d<3; d++) X[d] = x[d];

  if(stat) return 0;

  for(d=0; d<3; d++)
    if(dless(X[d],elm->bbox[2*d]) || dless(elm->bbox[2*d+1],X[d]))
      return 0;

  /* round X to inside box */
  for(d=0; d<3; d++)
  {
    if(X[d] < elm->bbox[2*d])   X[d] = elm->bbox[2*d];
    if(X[d] > elm->bbox[2*d+1]) X[d] = elm->bbox[2*d+1];
  }

  return 1;
}
////////////////////////////////////////////////////////////////////////////



/*
in add_nfaces_outside_patch study:
                            =====
      nblist1 = leafdescendants_along_face(nb, nb_f, NULL);

      touch = common_facepoints(node,face, nb,nb_f);

in common_facepoints study:
  f1 = find_nodefacepoints_in_nbface(node,f, nb,nb_f);

*/


// equivalent to find_nodefacepoints_in_nbface:
/* find out if any elm points on face f are on face nb_f of elm nb */
int find_elmfacepoints_in_nbface(tElm *elm, int f, tElm *nb, int nb_f)
{
  double *bbox  = elm->bbox;
  int dir = f/2;
  int n[] = { 3,3,3 };        /* we use 3 points */
  double X0[3], LX[3], dX[3]; /* grid of points */
  int dd;
  int i,j,k, plane, ret0, ret;

  /* make a grid of points, that excludes endpoints */
  for(dd=0; dd<3; dd++)
  {
    X0[dd] = bbox[2*dd];
    LX[dd] = bbox[2*dd+1] - X0[dd];
    dX[dd] = LX[dd]/(n[dd]);
    X0[dd] += dX[dd] * 0.5;
  }

  /* loop over points */
  plane = (n[dir] - 1) * (f%2);
  forplaneN(dir, i,j,k, n, plane)
  {
    double X[3], x[3], oX[3];
    int nbface[6];

    /* point grid, that never includes edges */
    X[0] = X0[0] + dX[0] * i;
    X[1] = X0[1] + dX[1] * j;
    X[2] = X0[2] + dX[2] * k;

    /* pick one of X,Y,Z on boundary */
    X[dir] = bbox[f];

    /* get x,y,z of X,Y,Z and then oX,oY,oZ in nb */
    //set_xyz(NULL, (tNode *)elm,-1, X, x);
    set_xyz(NULL, elm,-1, X, x);
    ret0 = elmXYZ_of_xyz(nb,-1, oX, x);

    /* try another point, if this one is not in nb */
    if(!ret0) continue;

    /* check if this point is on nb face */
    ret = XYZ_on_face(nb->pat, nbface, oX);

    /* if this point is only in face nb_f of nb we are done */
    if(ret==1 && nbface[nb_f]) return 1;

    /* if this point is in several faces try another point */
    if(ret>1) continue;

    if(ret==0)
    {
      errorexiti("oX was supposed to be on 1 face, not %d faces!!!", ret);
    }
  }

  return 0;
}









/**/
tElm *amr_get_parent(tElm *elm)
{
return NULL;
}


void amr_get_fnb(tElm *elm, int patface, int *nfnb, tElm **fnb)
{
}


/* myelm contains all elms on this proc
   nbelm contains all elms that are neighbors on other procs
   both can be searched to find a specific nb of one elm in myelm. */

/* NOTE: myelm and nbelm need to be sorted (use qsort) for searching
   we also keep a linked list for myelm to easily remove or add elms */

/* we need: function to find nb in all myelm */
//... use wolfGIT/c/binarysearch.c
// 1. search using comparfunc that is equal even grandparents agree
// 2. search using comparfunc that is equal even parents agree
// 3. search using comparfunc that is equal if elm themselves agree
// 4. search using comparfunc that is equal if ...


/* we need: function to find nb in all nbelm */
//...



/* pick nb location on face f of elm (currently in same patch and thus
   on same level), and write its loc into nbeloc */
void amr_get_nbeloc_nbface(const tElm *elm, int elmface,
                           tEloc nbeloc[1], int *nbface)
{
  const tEloc *eloc = elm->eloc;

  //FIXME: call connections_get_nb_eloc_face

  nbeloc->p = eloc->p; // so far look in same pat
  nbeloc->l = eloc->l;
  connections_get_nbloc_SameLevel_InsidePat(eloc->l, eloc->loc, elmface,
                                            nbeloc->loc, nbface);
}

/* Look in elm-array arr (in [arr+off,arr+num-1]) to find the elm
   with loc s_eloc and face s_f.
   *s_eloc is a loc where we start searching
   *But we start searching 1st for s_eloc's ancestor on level l0.
   *Return list with elms on face s_f */
int amr_elms_on_eloc_face(long narr, const tElm **arr,
                          size_t off0, size_t num0,
                          tEloc s_eloc[1], int s_f, int l0,
                          struct list_head *f_elms_head)
{
  tElm **f_elm;
  size_t off, num;
  //tEloc *eloc = elm->eloc;
  tEloc f_eloc[1];
  tEloc cheloc[1];
  int l, lret;
  int lmax = s_eloc->l;;
  int mor, ijk;

  /* init */
  off = off0;
  num = num0;

  PRF;printf(": off=%zu num=%zu s_eloc=", off, num);
  printeloc_s(s_eloc, "\n  ");
  printelmarray(narr, arr);

  for(l = l0; l <= lmax; l++)
  {
Yo(l);
    /* search for ancestor of s_eloc of level l */
    f_eloc[0] = s_eloc[0];
    f_eloc->l = l;
    f_elm = binarysearch(f_eloc, arr, &off, &num, sizeof(*arr), lecmp, NULL);

//printf("******** This changes all the time:\n");
printf("off=%zu num=%zu  f_elm pos=%zu\n",
off, num, (size_t) ((const tElm **)f_elm - arr));
printf("got ");printelm(*f_elm);

    if(!f_elm) return -s_eloc->l - 1000; /* found nothing */

    /* is there only one f_elm? */
    mor=binarysearchmore(f_eloc, arr, narr, sizeof(*arr), f_elm, lecmp, NULL);
printf("mor=%d\n", mor);
    if(!mor)  /* if there is only one */
    {
      //add f_elm to list and then return
      glist_entry_add_tail(*f_elm, f_elms_head);
      return l;
    }
  }
  /* if we get here, nb at s_eloc has children */

Yo(555);
Yo(l);
  /* search on children one level higher */
  l = s_eloc->l + 1;
  cheloc[0] = s_eloc[0];
  cheloc->l = l;
  /* we set cheloc->loc below */

  /* default return value */
  lret = -l - 1000; /* found nothing */

  /* get the 4 children elocs on nb face s_f */
  for(ijk = 0; ijk<8; ijk++) /* loop over children */
  {
    if(connections_ijk_is_at_parentface(ijk, s_f))
    {
      int ret;
      /* child ijk */
      cheloc->loc[l-1] = ijk + '0';
      if(l<LOCSMAX) cheloc->loc[l] = 0;
      ret = amr_elms_on_eloc_face(narr, arr, off, num, cheloc, s_f, l,
                                  f_elms_head);
      if(ret >=0) lret = l; /* record success */
    }
  }

  /* finally signal failure or success with at least one nb child */
  printf("final lret=%d\n", lret);
  return lret;
}

/* Look in elm-array arr (in [arr+off,arr+num-1]) to find the nb of
   elm on face elmface. */
int amr_set_fnb_list(tElm *elm, int elmface, long narr, const tElm **arr,
                     struct list_head *fnb_head)
{
  tEloc *eloc = elm->eloc;
  tEloc nbeloc[1];
  int nb_f;

  PRFs(": ");printeloc(elm->eloc);printf(" f=%d", elmface);

  /* Before calling amr_elms_on_eloc_face, set nb_f and nbeloc by calling
     amr_get_nbeloc_nbface(elm, elmface, nbeloc, &nb_f); */
  amr_get_nbeloc_nbface(elm, elmface, nbeloc, &nb_f);
  printf(" -> nbeloc=");printeloc(nbeloc);printf(" nb_f=%d\n", nb_f);

  return
    amr_elms_on_eloc_face(narr, arr, 0, narr, nbeloc, nb_f, 0, fnb_head);
}


/****************************************************************************/
/* functions to exchange info between rank_i and rank_{i \pm 1} */
/****************************************************************************/

/* exchange first and last elm between neigboring ranks
  Out: myfl  <-- first and last elm on this rank
       rm1fl <-- first and last on rank-1
       rp1fl <-- first and last on rank+1 */
void elmfl_exchange_between_nbranks(tMesh *mesh, tElmfl myfl[1],
                                    tElmfl rm1fl[1], tElmfl rp1fl[1])
{
  int size = nMPI_size();
  int rank = nMPI_rank();
  nMPI_Req req[4];
  nMPI_Stat stat[4];
  int nreqs, tag;
  /* 1st and last entry in mesh->myelm_head list (can be NULL) */
  tElm *first = NULL;
  tElm *last  = NULL;

  if(!list_empty(&mesh->myelm_head))
    first = list_first_entry(&mesh->myelm_head, tElm, list);

  /* if first != NULL there are entries */
  if(first)
  {
    last = list_last_entry(&mesh->myelm_head, tElm, list);
    myfl->elm_fl[0] = first[0]; /* shallow copy */
    myfl->elm_fl[1] = last[0];  /* shallow copy */
    myfl->nelms = mesh->nmyelm;
  }
  else /* zero all of myfl */
  {
    memset(&myfl[0], 0, sizeof(myfl[0]));
  }

  /* MPI exchanges */
  nreqs = 0;
  tag = 1;
  if(rank > 0)       /* send info to rank-1 */
    nMPI_Isend((char *)myfl,sizeof(myfl[0]), nMPI_CHAR, rank-1, tag,
               WORLD, &req[nreqs++]);
  if(rank < size-1)  /* Receive info from rank+1 */
    nMPI_Irecv((char *)rp1fl,sizeof(rp1fl[0]), nMPI_CHAR, rank+1, tag,
               WORLD, &req[nreqs++]);

  tag = 2;
  if(rank < size-1)  /* send info to rank+1 */
    nMPI_Isend((char *)myfl,sizeof(myfl[0]), nMPI_CHAR, rank+1, tag,
               WORLD, &req[nreqs++]);
  if(rank > 0)       /* Receive info from rank-1 */
    nMPI_Irecv((char *)rm1fl,sizeof(rm1fl[0]), nMPI_CHAR, rank-1, tag,
               WORLD, &req[nreqs++]);

  /* wait until all MPI requests are done */
  nMPI_Waitall(nreqs, req, stat);

  /**/
  /*FIXME: update rflag by using nb rank info,
           if all 8 siblings are on different ranks,
           run MPI exchanges up to 5 times more */
}

/****************************************************************************/
/* functions to build elm lists */
/****************************************************************************/
//...
