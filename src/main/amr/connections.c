/* connections.c */
/* Wolfgang Tichy, 3/2023 */

#include "nmesh.h"
#include "amr.h"


// move this into .h file
int lecmp(const void *loc, const void *elem, void *arg);
int loccmp(const void *loc, const void *eloc);


/****************************************************************************/
/* primitive funxtion that work on integers and strings */
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
  return -1;
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

/* Out: nbloc */
int connections_get_nbloc_SameLevel_InsidePat(int l, const char loc[LOCSMAX],
                                              int face,
                                              char nbloc[LOCSMAX])
{
  int nfaces, patface[6];
  int ijk, nb_ijk;

  nfaces = connections_loc_on_patchface(l,loc, patface);

  if(patface[face])
  {
    errorexit("deal with pat face");
    /*
    tPat *pat = node->pat;
    tBface *bfaces = pat->bfaces[face];
    // loop over bfaces
    forbfacesonface(pat, f, bface) ;
    // same as: for(bface=bfaces; bface; bface=bface->next) ;

    //if(bfaces && bfaces->boundary==OUTERBOUND)
    */
  }


  /* find ijk of node and ijk of nb */
  ijk = connections_get_ijk(l, loc);
  nb_ijk = connections_get_inner_nb_ijk(ijk, face/2);

  if(connections_ijk_is_at_parentface(ijk, face))
  {
    char pnbloc[LOCSMAX]; /* location of parent nb */
    /* l-1,loc is parent */
    connections_get_nbloc_SameLevel_InsidePat(l-1,loc, face, pnbloc);
    pnbloc[l] = 0; /* add string-end marker */
    strncpy(nbloc, pnbloc, LOCSMAX);
    nbloc[l-1] = nb_ijk;
    if(l<LOCSMAX) nbloc[l] = 0;
    return l;
  }
  else
  {
    strncpy(nbloc, loc, LOCSMAX);
    nbloc[l-1] = nb_ijk;
    return l;
  }
  return l;
}


/****************************************************************************/
/* functions to initialize tElm */
/****************************************************************************/

/* find patch of elm and save it in elm->pat */
void amr_set_elm_pat(tMesh *mesh, tElm *elm)
{
  tEloc *eloc = elm->eloc;
  int p = eloc->p;
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
      if(ijk & bt) bbox[2*d] = bbox[2*d+1] - LX[2*d];
      else         bbox[2*d+1] = bbox[2*d] + LX[2*d];
    }
  }
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



/* pick nb location (with index inbu2) at 2 levels up from elm,
   and write loc into nbu2eloc */
void amr_get_nbu2loc(const tElm *elm, int f, int inbu2,
                     tEloc nbu2eloc[1])
{
  /* pick nb loc (with index inbu2) at 2 levels up */
  nbu2eloc->p = 000; //???
  nbu2eloc->l = elm->eloc->l+2; //???
  //nbu2eloc->loc = "12352"; //???
}


/* Look in elm array arr to find all nb of elm on face f with nb loc
   that is 2 levels up.*/
void amr_set_fnb(int narr, const tElm **arr, const tElm *elm, int f,
                 tElm *fnb[11111])
{
  int inbu2;
  tEloc nbu2eloc[1];

  //for inbu2:
  {
    int lf;
    int nbi=0; //???

    /* pick nb loc (with index inbu2) at 2 levels up */
    amr_get_nbu2loc(elm, f, inbu2, nbu2eloc);
    lf = amr_set1_fnb(narr,arr, elm, f, nbu2eloc, fnb[nbi]);
    //if lf ...
  }

}

/* Look in elm array arr to find the nb of elm on face f with nb loc
   nbu2eloc. nbu2eloc is a nb that is 2 levels up.*/
int amr_set1_fnb(int narr, const tElm **arr, const tElm *elm,
                 int f, tEloc nbu2eloc[1],
                 tElm *fnb[1])
{
  tElm *nbelm;
  size_t off, num;
  //tEloc *eloc = elm->eloc;
  tEloc nbfeloc[1];
  int mor;

  /* pick nb loc (with index inbu2) at 2 levels up */
  int inbu2=0; //???
  amr_get_nbu2loc(elm, f, inbu2, nbu2eloc);

  /* init */
  off = 0;
  num = narr;

  if(nbfeloc->l > 4)
  {
    /* search for grand-grand-grand parent of nbu2eloc (l-4) */
    nbfeloc[0] = nbu2eloc[0];
    nbfeloc->l = nbu2eloc->l - 4;
    nbelm = binarysearch(nbfeloc, arr, &off, &num, sizeof(*arr), lecmp, NULL);
    if(!nbelm) return -9999;
  }

  /* save nbelm, may need to alloc fnb ??? */
  fnb[0] = nbelm;
  mor=binarysearchmore(nbfeloc, arr, narr, sizeof(*arr), nbelm, lecmp, NULL);
  if(!mor) return -4; /* if there is only one */

  if(nbfeloc->l > 3)
  {
    /* search for grand-grand parent of nbu2eloc (l-3) */
    nbfeloc[0] = nbu2eloc[0];
    nbfeloc->l = nbu2eloc->l - 3;
    nbelm = binarysearch(nbfeloc, arr, &off, &num, sizeof(*arr), lecmp, NULL);
    if(!nbelm) return -4;
  }

  /* save nbelm, may need to alloc fnb ??? */
  fnb[0] = nbelm;
  mor=binarysearchmore(nbfeloc, arr, narr, sizeof(*arr), nbelm, lecmp, NULL);
  if(!mor) return -3; /* if there is only one */

  /* search for grand parent of nbu2eloc (l-2) */
  nbfeloc[0] = nbu2eloc[0];
  nbfeloc->l = nbu2eloc->l - 2;
  nbelm = binarysearch(nbfeloc, arr, &off, &num, sizeof(*arr), lecmp, NULL);
  if(!nbelm) return -3;

  /* save nbelm, may need to alloc fnb ??? */
  fnb[0] = nbelm;
  mor=binarysearchmore(nbfeloc, arr, narr, sizeof(*arr), nbelm, lecmp, NULL);
  if(!mor) return -2; /* if there is only one */

  /* search for parent of nbu2eloc (l-1) */
  nbfeloc[0] = nbu2eloc[0];
  nbfeloc->l = nbu2eloc->l - 1;
  nbelm = binarysearch(nbfeloc, arr, &off, &num, sizeof(*arr), lecmp, NULL);
  if(!nbelm) return -2;

  /* save nbelm, may need to alloc fnb ??? */
  fnb[0] = nbelm;
  mor=binarysearchmore(nbfeloc, arr, narr, sizeof(*arr), nbelm, lecmp, NULL);
  if(!mor) return -1; /* if there is only one */

  /* search for nbu2eloc (l) */
  nbfeloc[0] = nbu2eloc[0];
  nbfeloc->l = nbu2eloc->l;
  nbelm = binarysearch(nbfeloc, arr, &off, &num, sizeof(*arr), lecmp, NULL);
  if(!nbelm) return -1;

  /* save nbelm, may need to alloc fnb ??? */
  fnb[0] = nbelm;
  mor=binarysearchmore(nbfeloc, arr, narr, sizeof(*arr), nbelm, lecmp, NULL);
  if(!mor) return 0; /* if there is only one */
  else errorexit("2 levels up there should be only one nb");
}





/* return -1,0,1 if loc is before,at,after elem location */
int lecmp(const void *loc, const void *elem, void *arg)
{
  const tEloc *lc = (const tEloc *) loc;
  const tElm *elm = (const tElm *) elem;
  const tEloc *elc = elm->eloc;
  return loccmp(lc, elc);
}


/* return -1,0,1 if loc is before,at,after eloc */
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
