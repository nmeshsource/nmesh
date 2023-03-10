/* connections.c */
/* Wolfgang Tichy, 3/2023 */

#include "nmesh.h"
#include "amr.h"




/* find ijk from loc */
int connections_get_ijk(int l, const char loc[LOCSMAX])
{
  if(l<1) return 0;
  return loc[l-1] - '0';
}

/* return 1 if ijk touches node-face */
int connections_ijk_is_at_nodeface(int ijk, int face)
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

/* get ijk of nb in same level, assuming ijk is not on patch face */
int connections_get_inner_nb_ijk(int ijk, int dir)
{
  int ns[] = {2,2,2};
  int i,j,k;  /* index of loc */

  if(l<1) return 0;

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
int connections_loc_on_patch_face(int l, const char loc[LOCSMAX],
                                  int patface[6])
{
  int i, f, npatfaces;

  npatfaces = 6;
  for(f=0; f<6; f++) patface[f] = 1;

  if(l<1) return npatfaces;

  for(i=0; i<l; i++)
  {
    int ijk = connections_get_ijk(i, loc);

    for(f=0; f<6; f++)
    {
      if(!connections_ijk_is_at_nodeface(ijk, f))
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
  int ijk;

  nfaces = connections_loc_on_patch_patface(l,loc, patface);

  if(patface[face])
  {
    errorexit("deal with pat face");

    tPat *pat = node->pat;
    tBface *bfaces = pat->bfaces[face];
    // loop over bfaces
    forbfacesonface(pat, f, bface) ;
    // same as: for(bface=bfaces; bface; bface=bface->next) ;

    //if(bfaces && bfaces->boundary==OUTERBOUND)

  }


  /* find ijk of node and ijk of nb */
  ijk = connections_get_ijk(l, loc);
  nb_ijk = connections_get_inner_nb_ijk(ijk, face/2);

  if(connections_ijk_is_at_nodeface(ijk, face))
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


/* look in elm array arr to find nbs of elm on face f */
void amr_AllocAndSet_fnb(int narr, const tElm **arr, const tElm *elm, int f,
                         int *nfnb, tElm **fnb)
{
  tElm *nbelm;
  size_t off, num;
  tEloc *eloc = elm->eloc;
  tEloc nbeloc[1];
  tEloc nbfeloc[1];

  /* find nb at same level */
  nbeloc->p = 000; //???
  nbeloc->l = l; //???
  nbeloc->loc = "12352"; //???

  /* init */
  *nfnb = 0 ;
  off = 0;
  num = narr;

  /* search for grand parent of nbeloc (l-2) */
  nbfeloc[0] = nbeloc[0];
  nbfeloc->l = nbeloc->l - 2;
  nbelm = binarysearch(nbfeloc, arr, &off, &num, sizeof(*arr), loccmp, NULL);
  if(!nbelm) return;

  /* save nbelm, may need to alloc fnb ??? */
  fnb[*nfnb] = nbelm;
  *nfnb = 1;
  if(num<=1) return; /* if there is only one */

  /* search for parent of nbeloc (l-1) */
  nbfeloc[0] = nbeloc[0];
  nbfeloc->l = nbeloc->l - 1;
  nbelm = binarysearch(nbfeloc, arr, &off, &num, sizeof(*arr), loccmp, NULL);
  if(!nbelm) return;

  /* save nbelm, may need to alloc fnb ??? */
  fnb[*nfnb] = nbelm;
  *nfnb = 1;
  if(num<=1) return; /* if there is only one */

  /* search for nbeloc */
  nbfeloc[0] = nbeloc[0];
  nbfeloc->l = nbeloc->l - 1;
  nbelm = binarysearch(nbfeloc, arr, &off, &num, sizeof(*arr), loccmp, NULL);
  if(!nbelm) return;

  /* save nbelm, may need to alloc fnb ??? */
  fnb[*nfnb] = nbelm;
  *nfnb = 1;
  if(num<=1) return; /* if there is only one */

  /* if we get here, there are several neighbors */
  //...
  errorexit("deal with several neighbors");
}

/* return -1,0,1 if loc is before,at,after elem location */
int loccmp(const void *loc, const void *elem, void *arg)
{
  tEloc *lc = (tEloc *) loc;
  tElm *elm = (tElm *) elem;
  tEloc *el = elm->eloc;
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
      else                         return -1
    }
    return 0; /* lc and el are equal up the first lc->l */
  }
  else
  {
    for(i=0; i<el->l; i++)
    {
      if(lc->loc[i] == el->loc[i]) continue;
      if(lc->loc[i] >  el->loc[i]) return  1;
      else                         return -1
    }
    return 1; /* make binarysearch move to right */
  }
}
