/* connections.c */
/* Wolfgang Tichy, 3/2023 */

#include "nmesh.h"
#include "amr.h"


extern tnMPIvars nMPIvars[1];
extern tAMR amr[1];


/****************************************************************************/
/* functions that determine the order of element locations */
/****************************************************************************/

/* We need functions to find elm in myelm or nbelm.
   For this we use binarysearch.
   Note: myelm and nbelm need to be sorted (use qsort) for searching.
   We also keep a linked list for myelm to easily remove or add elms */


/* Function that orders locations described in in tEloc:
   return -1,0,1 if loc is before,at,after eloc
   This is used in our special binarysearch way with mode=0.
   It can also be used in qsort, or bsearch if mode=1.
   About mode:
   0: return  0 (equal) if el and lc agree
      return  0 (equal) even if el->l > lc->l, as long as first lc->l agree
      return  1 (lc>el) if el->l < lc->l, even if first el->l agree
   1: return  0 (equal) if el and lc agree
      return -1 (lc<el) if el->l > lc->l, even if first lc->l agree
      return  1 (lc>el) if el->l < lc->l, even if first el->l agree
   2: return  0 (equal) if el and lc agree
      return  0 (equal) even if el->l > lc->l, as long as first lc->l agree
      return  0 (equal) even if el->l < lc->l, as long first el->l agree
   3: return  0 (equal) if el and lc agree
      return -1 (lc<el) if el->l > lc->l, even if first lc->l agree
      return  0 (equal) even if el->l < lc->l, as long first el->l agree
   mode 1 is good for normal qsort or bsearch
   mode 0 is good for finding nbs even if they are more refined
   mode 3 is good for finding nbs even if they are less refined
   mode 2 is probably good for nothing...  */
int loccmp_mode(const void *loc, const void *eloc, int mode)
{
  const tEloc *lc = loc;
  const tEloc *el = eloc;
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
    if( (mode & 1) && (el->l > lc->l) )
      return -1; /* make el greater than lc, and thus move to left */
    else
      return 0;  /* lc and el are equal up the first lc->l */
  }
  else
  {
    for(i=0; i<el->l; i++)
    {
      if(lc->loc[i] == el->loc[i]) continue;
      if(lc->loc[i] >  el->loc[i]) return  1;
      else                         return -1;
    }
    if( (mode & 2) )
      return 0; /* lc and el are equal up the first el->l */
    else
      return 1; /* make binarysearch move to right */
  }
}

/* Function that orders locations described in in tEloc:
   return -1,0,1 if loc is before,at,after eloc
   This is used in our special binarysearch way. */
/*
int loccmp(const void *loc, const void *eloc)
{
  return loccmp_mode(loc, eloc, 0);
}
*/

/* return -1,0,1 if loc is before,at,after elem location,
   this is used in binarysearch */
int lecmp_0(const void *loc, const void *elem, void *arg)
{
  const tEloc *lc = loc;
  //const tElm0 **elm0_arr = elem;
  const tElm0 *const*elm0_arr = elem;
  const tEploc *pelc = elm0_arr[0]->eploc;
  tEloc elc[1];
  int cmp;

  /* NOTE: we could optimize this by caching an unpacked loc in each elm0: */
  eloc_from_eploc(elc, pelc);

  cmp = loccmp_mode(lc, elc, 0);
  //PRFs(": ");printeloc_s(lc, " ");printeloc_s(elc, " ");
  //printf("--> cmp=%d\n", cmp);
  ////printelm0(elm0_arr[0]);
  return cmp;
}

/* return -1,0,1 if key_elem location is before,at,after elem location,
   this is used in binarysearch */
int eecmp_0(const void *key_elem, const void *elem, void *arg)
{
  const tElm0 *const*kelm0 = key_elem;
  const tEploc *keploc = kelm0[0]->eploc;
  const tElm0 *const*elm0_arr = elem;
  const tEploc *eploc = elm0_arr[0]->eploc;
  tEloc klc[1], elc[1];
  int cmp;

  /* NOTE: we could optimize this by caching an unpacked loc in each elm0: */
  eloc_from_eploc(klc, keploc);
  eloc_from_eploc(elc, eploc);

  cmp = loccmp_mode(klc, elc, 0);
  //PRFs(": ");printeloc_s(klc, " ");printeloc_s(elc, " ");
  //printf("--> cmp=%d\n", cmp);
  ////printelm0(elm0_arr[0]);
  return cmp;
}


/* return -1,0,1 if loc is before,at,after elem location,
   this is used in qsort and bsearch, and returns 0 only if they are
   strictly equal */
int lecmp_q(const void *loc, const void *elem)
{
  const tEloc *lc = loc;
  //const tElm0 **elm0_arr = elem;
  const tElm0 *const*elm0_arr = elem;
  const tEploc *pelc = elm0_arr[0]->eploc;
  tEloc elc[1];
  int cmp;

  /* NOTE: we could optimize this by caching an unpacked loc in each elm0: */
  eloc_from_eploc(elc, pelc);

  cmp = loccmp_mode(lc, elc, 1);
  //PRFs(": ");printeloc_s(lc, " ");printeloc_s(elc, " ");
  //printf("--> cmp=%d\n", cmp);
  ////printelm0(elm0_arr[0]);
  return cmp;
}

/* return -1,0,1 if key_elem location is before,at,after elem location,
   this is used qsort and bsearch, and returns 0 only if they are
   strictly equal */
int eecmp_q(const void *key_elem, const void *elem)
{
  const tElm0 *const*kelm0 = key_elem;
  const tEploc *keploc = kelm0[0]->eploc;
  const tElm0 *const*elm0_arr = elem;
  const tEploc *eploc = elm0_arr[0]->eploc;
  tEloc klc[1], elc[1];
  int cmp;

  /* NOTE: we could optimize this by caching an unpacked loc in each elm0: */
  eloc_from_eploc(klc, keploc);
  eloc_from_eploc(elc, eploc);

  cmp = loccmp_mode(klc, elc, 1);
  //PRFs(": ");printeloc_s(klc, " ");printeloc_s(elc, " ");
  //printf("--> cmp=%d\n", cmp);
  ////printelm0(elm0_arr[0]);
  return cmp;
}

/****************************************************************************/
/* primitive functions that work on integers and strings */
/****************************************************************************/

/* find ijk from l,loc by reading last in loc */
int connections_get_ijk(int l, const char loc[NLOCS])
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


/* check if l,loc is on patchface f
   Out: return value: 1 if on patch face, 0 if not */
int connections_loc_on_patchface_f(int l, const char loc[NLOCS], int f)
{
  int ll;

  //PRF;printf(": l=%d loc=%s f=%d\n", l, loc, f);
  for(ll=1; ll<=l; ll++)
  {
    int ijk = connections_get_ijk(ll, loc);
    //printf("  %d in %s is ijk=%d\n", ll, loc, ijk);
    if(!connections_ijk_is_at_parentface(ijk, f))
      return 0;
  }
  return 1;
}

/* Out: return value: number of faces l,loc is on
        patface[f] = 1 if l,loc is on patch face f */
int connections_loc_on_patchface(int l, const char loc[NLOCS],
                                 int patface[6])
{
  int npatfaces, f;
  for(npatfaces=0, f=0; f<6; f++)
  {
    patface[f] = connections_loc_on_patchface_f(l, loc, f);
    npatfaces += patface[f];
  }
  //for(f=0; f<6; f++) printf("patface[%d]=%d ", f, patface[f]);
  //printf("npatfaces=%d\n", npatfaces);
  return npatfaces;
}

/* Out: return value: number of faces l,loc is on
        patface[f] = 1 if l,loc is on patch face f */
int connections_loc_on_patchface_WRONG(int l, const char loc[NLOCS],
                                 int patface[6])
{
  int ll, f, npatfaces;

  /* set result if l=0 */
  npatfaces = 6;
  for(f=0; f<6; f++) patface[f] = 1;

  if(l<1) return npatfaces; // <--- not needed

  printf("WWWWWWWWWWWWWWWWWWWWWWWW\n");
  PRF;printf(": l=%d loc=%s\n", l, loc);
  for(ll=1; ll<=l; ll++)
  {
    int ijk = connections_get_ijk(ll, loc);
    PRF;printf(": %d in %s is ijk=%d\n", ll, loc, ijk);

    for(f=0; f<6; f++)
    {
      if(!connections_ijk_is_at_parentface(ijk, f))
      {
        npatfaces--;
        patface[f] = 0;
        PRF;printf("f%d: %d=>\n", f, connections_ijk_is_at_parentface(ijk, f));
      }
      if(npatfaces<=0) break;
    }
    if(npatfaces<=0) break;
  }
  printf("npatfaces=%d\n", npatfaces);
  return npatfaces;
}


/* get nbloc,nf of neighbor on face of elm with l,loc,face
   BUT this func works only if the face is not on patch face.
   In:  l, loc, face
   Out: nbloc, nb_f */
int connections_get_nbloc_InsidePat(int l, const char loc[NLOCS], int face,
                                    char nbloc[NLOCS], int *nb_f)
{
  //int patface[6];
  int ijk, nb_ijk;
  int lret;

  //PRF;printf(": l=%d loc=%s face=%d\n", l, loc, face);

  //nfaces = connections_loc_on_patchface(l,loc, patface);

  if(connections_loc_on_patchface_f(l,loc, face))
    errorexiti("face%d of loc is on patch surface", face);

  /* find ijk of node and ijk of nb */
  ijk = connections_get_ijk(l, loc);
  nb_ijk = connections_get_inner_nb_ijk(ijk, face/2);
  *nb_f  = face^1;
  //printf("  ijk=%d nb_ijk=%d *nb_f=%d\n", ijk, nb_ijk, *nb_f);

  if(connections_ijk_is_at_parentface(ijk, face))
  {
    int pl;
    char pnbloc[NLOCS]; /* location of parent nb */
    int pnb_f;
    //PRF;printf(" at parentface\n");
    ////printf(" pnbloc=%s\n", pnbloc);

    /* l-1,loc is parent, write parent nb loc into pnbloc */
    pl = connections_get_nbloc_InsidePat(l-1,loc, face, pnbloc, &pnb_f);
    pnbloc[l-1] = 0; /* add string-end marker */
    //printf("   pnbloc=%s\n", pnbloc);
    strncpy(nbloc, pnbloc, NLOCS);
    //printf("    nbloc=%s?  pl=%d\n", nbloc, pl);
    lret = pl;
  }
  else
  {
    strncpy(nbloc, loc, NLOCS);
    lret = l;
  }
  nbloc[l-1] = nb_ijk + '0';
  if(l<NLOCS) nbloc[l] = 0; /* add string-end marker */
  //printf("  nbloc=%s\n", nbloc);
  return lret;
}


/****************************************************************************/
/* functions that may need to be moved into other files */
/****************************************************************************/

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


// equivalent to find_nodefacepoints_in_nbface:
/* find out if any elm points on face f are on face nb_f of elm nb,
   Returns: 1 if elm,f and nb,nb_f have points in common
            0 otherwise */
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
    ret = XYZ_on_elmface(nb, nbface, oX);
    //This used to be:  ret = XYZ_on_face(nb->pat, nbface, oX);
    //Which was very wrong, because we don't want to know it's on a patface

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

/* check if elm and nb has common points on faces f and nb_f
   since res in find_elmfacepoints_in_nbface is low, try it both ways,
   Returns: 1 if elm,f and nb,nb_f have points in common
            0 otherwise */
int elm_common_facepoints(tElm *elm, int f, tElm *nb, int nb_f)
{
  int f1, f2;

  f1 = find_elmfacepoints_in_nbface(elm,f, nb,nb_f);
  if(f1) return 1;

  f2 = find_elmfacepoints_in_nbface(nb,nb_f, elm,f);

  return f2 || f1;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////



/****************************************************************************/
/* functions that work on eloc */
/****************************************************************************/

/* unpack ploc into loc */
void connections_loc_from_ploc(char loc[NLOCS],
                               const unsigned char ploc[NPBYTES])
{
  int i, p1,p2, b1,b2, bi1,bi2;
  unsigned char ch, c1,c2;
  //p2=0;

  /* translate ploc to loc */
  for(i=0; i<NLOCS; i++)
  {
    int i3 = i*3;
    b1 = i3;        /* bit number of 1st bit in ploc */
    b2 = i3+2;      /* bit number of 3rd bit in ploc */

    p1 = b1/8;      /* position of 1st byte in ploc */

    bi1 = b1 % 8;   /* starting bit index we need in c1 */
    bi2 = b2 % 8;   /* starting bit index we need in c2 */

    /* if all is in c1 */
    c1 = ploc[p1];
    ch = c1>>bi1;   /* shift 1st bit into correct position */
    ch = ch & 7;    /* keep only lowest 3 bits */
    //printf("ch=%o\n",ch);
    if(bi2<2) /* not all is in c1 */
    {
      p2 = b2/8;       /* position of 2nd byte in ploc */
      c2 = ploc[p2];
      //printf("c2=%o\n",c2);
      c2 = c2<<(2-bi2); /* shift 3rd bit into correct position */
      //printf("c2=%o\n",c2);
      c2 = c2 & 7;     /* keep only lowest 3 bits */
      //printf("c2=%o\n",c2);
      ch = ch | c2;
    }
    loc[i] = ch + '0'; /* set loc */

    //PRF;printf(": i=%d  b1=%d b2=%d bi1=%d bi2=%d p1=%d p2=%d %c(%d)<-%o",
    //i, b1,b2, bi1,bi2, p1,p2, loc[i],loc[i]-'0', ploc[p1]);
    //if(bi2<2) printf(",%o", ploc[p2]);
    //printf("\n");
  }
  /* do not add string end marker in loc[eloc->l], this would kill
     stuff below eloc->l which we intend to keep, because sometimes we
     just decrease l in eloc */
}

/* pack loc into ploc */
void connections_loc_to_ploc(const char loc[NLOCS],
                             unsigned char ploc[NPBYTES])
{
  int i, p1,p2, b1,b2, bi1,bi2;
  unsigned char ch, c1,c2, pc;
  //p2=0;

  ploc[0] = 0;         /* init first char in ploc */
  //ploc[NPBYTES-1] = 0; /* init last char in ploc */

  /* translate ploc to loc */
  //ploc[0] = 0; // not needed
  for(i=0; i<NLOCS; i++)
  {
    int i3 = i*3;
    b1 = i3;        /* bit number of 1st bit in ploc */
    b2 = i3+2;      /* bit number of 3rd bit in ploc */

    p1 = b1/8;      /* position of 1st byte in ploc */

    bi1 = b1 % 8;   /* starting bit index we need in c1 */
    bi2 = b2 % 8;   /* starting bit index we need in c2 */

    ch = loc[i];
    //if(ch == 0) break; //Note: do not brak like this!!!
    //ch = ch - '0';
    ch = ch & 7; // this also subtracts '0'
    c1 = ch<<bi1;   /* shift 1st bit into correct position */
    pc = ploc[p1];
    pc = pc<<(8-bi1);   /* clear all left of bi1 */
    pc = pc>>(8-bi1);   /* still needed if ploc[0]=0 */
    ploc[p1] = pc | c1; /* then add the bits form c1 */

    if(bi2<2) /* not all is in one ploc */
    {
      c2 = ch>>(2-bi2); /* shift 3rd bit into correct position */
      p2 = b2/8;        /* position of 2nd byte in ploc */
      ploc[p2] = c2;    /* this also zeros all above the bit bi2 */
    }
    //PRF;printf(": i=%d  b1=%d b2=%d bi1=%d bi2=%d p1=%d p2=%d %c(%d)->%o",
    //i, b1,b2, bi1,bi2, p1,p2, loc[i],loc[i]-'0', ploc[p1]);
    //if(bi2<2) printf(",%o", ploc[p2]);
    //printf("\n");
  }
}

/* unpack eploc into eloc */
void eloc_from_eploc(tEloc eloc[1], const tEploc eploc[1])
{
  const unsigned char *ploc = eploc->ploc;
  char *loc  = eloc->loc;

  /* trivial copies */
  eloc->eid = eploc->eid;
  eloc->p = eploc->p;
  eloc->l = eploc->l;

  /* translate ploc to loc */
  connections_loc_from_ploc(loc, ploc);
}

/* pack eloc into eploc */
void eloc_to_eploc(const tEloc eloc[1], tEploc eploc[1])
{
  unsigned char *ploc = eploc->ploc;
  const char *loc  = eloc->loc;

  /* trivial copies */
  eploc->eid = eloc->eid;
  eploc->p = eloc->p;
  eploc->l = eloc->l;

  /* translate ploc to loc */
  connections_loc_to_ploc(loc, ploc);
}

/* test eloc_to_eploc and eloc_from_eploc */
void test_eploc(void)
{
  tEloc eloc[1];
  tEloc eloc2[1];
  tEloc eloc3[1];
  tEploc eploc[1];
  int i;

  for(i=0; i<NPBYTES; i++) eploc->ploc[i] = 255;

  for(i=0; i<NLOCS; i++) eloc->loc[i] = (i%8) +'0';
  eloc->loc[NLOCS-1]= 4+'0';
  eloc->p=9;
  eloc->l=NLOCS;
  eloc->eid=98;
  strcpy(eloc->loc, "1234567564321012345");
  eloc->loc[NLOCS-3] = '1';
  printeloc_s(eloc, "\n");

  //printf("eploc->ploc[NPBYTES-1]=%o\n", eploc->ploc[NPBYTES-1]);
  eloc_to_eploc(eloc, eploc);
  //eploc->ploc[NPBYTES-1]=128 + 32+  8+4 + 1;
  //printf("eploc->ploc[NPBYTES-1]=%o\n", eploc->ploc[NPBYTES-1]);
  eloc_from_eploc(eloc2, eploc);

  eloc_to_eploc(eloc2, eploc);
  eloc_from_eploc(eloc3, eploc);

  printeloc_s(eloc2, "\n");
  printeloc_s(eloc3, "\n");
}


/* check if a eplockey is in the array eploc
   return val: 1st position in peloc, where eplockey is found,
               -1 id not found */
int eploc_key_in_eplocarray(const tEploc *eplockey,
                            ulong n, const tEploc eploc[n])
{
  ulong i;

  //printf("n=%lu eplockey=",n);printeploc_s(eplockey, ":\n");
  for(i=0; i<n; i++)
  {
    int dif = memcmp(eplockey, &(eploc[i]), sizeof(eploc[0]));
    //printf("%lu:",i);printeploc(&eploc[i]);printf("dif=%i\n",dif);
    if(dif==0) return i;
  }
  return -1;
}


/* check if 2 elocs agree up to level l_max, return 0/1 if no/yes */
int eloc1_eloc2_agree_upto_l_max(const tEloc *eloc1,
                                 const tEloc *eloc2, int l_max)
{
  int p1 = eloc1->p;
  int l1 = eloc1->l;
  int p2 = eloc2->p;
  int l2 = eloc2->l;
  int i;

  if(p1!=p2)   return 0;
  if(l1<l_max) return 0;
  if(l2<l_max) return 0;

  for(i=0; i<l_max; i++)
  {
    if(eloc1->loc[i] != eloc2->loc[i]) return 0;
  }
  return 1;
}


/* check if 2 eplocs agree up to level l_max, return 0/1 if no/yes */
int eploc1_eploc2_agree_upto_l_max(const tEploc *eploc1,
                                   const tEploc *eploc2, int l_max)
{
  tEloc eloc1[1], eloc2[1];
  eloc_from_eploc(eloc1, eploc1);
  eloc_from_eploc(eloc2, eploc2);
  return eloc1_eloc2_agree_upto_l_max(eloc1, eloc2, l_max);
}


/* try to get the 8 elms in mesh->myelm_head, return how many we actually got
     this should really be:
     int amr_get_8elms_at_pos(tElm *elm_start, tElm *elm[8]);
     void *ptr_elm is really tElm *elm[8],
     but this way we can also pass in tElm0 *elm[8] */
int amr_get_8elms_at_elm_start(tElm *elm_start, void *ptr_elm)
{
  tElm **elm = ptr_elm; // like: tElm *elm[8];
  tMesh *mesh = Elm_mesh(elm_start);
  struct list_head *pos = (elm_start->list).prev;
  int cnt = 0;
  list_for_each_continue(pos, &mesh->myelm_head)
  {
    tElm *elm_i = list_entry(pos, tElm, list);
    elm[cnt++] = elm_i;
    if(cnt>=8) break;
  }
  return cnt;
}

/* try to get the 8 elms in mesh->myelm list, return how many we actually got
     this should really be:
     int amr_get_8elms_at_myid(tMesh *mesh, ulong myid, tElm *elm[8]);
     void *ptr_elm is really tElm *elm[8],
     but this way we can also pass in tElm0 *elm[8] */
int amr_get_8elms_at_myid(tMesh *mesh, ulong myid, void *ptr_elm)
{
  tElm **elm = ptr_elm; // like: tElm *elm[8];
  int cnt = 0;
  formyelms_s_n(mesh, myid, 8)
  {
    elm[cnt++] = MyElm;
  }
  return cnt;
}

/* check if all n in elm[8] are siblings
     this should really be:
     int amr_elms_are_siblings(int n, tElm0 *elm[8]);
     i.e. void *ptr_elm is really tElm0 *elm[8],
     but this way we can also pass in tElm *elm[8] */
int amr_elms_are_siblings(int n, void *ptr_elm)
{
  tElm0 **elm = ptr_elm; // like: tElm0 *elm[8];
  int p0 = Elm_p(elm[0]);
  int l0 = Elm_l(elm[0]);
  int are_sibs;
  int i;
  tEloc eloc0[1], eloci[1];

  eloc_from_eploc(eloc0, elm[0]->eploc);

  for(i=1; i<n; i++)
  {
    int pi = Elm_p(elm[i]);
    int li = Elm_l(elm[i]);
    if(pi!=p0) return 0;
    if(li!=l0) return 0;
    eloc_from_eploc(eloci, elm[i]->eploc);
    are_sibs = eloc1_eloc2_agree_upto_l_max(eloc0, eloci, l0-1);
    if(!are_sibs) return 0;
  }
  return 1;
}

/* write the info of elm's sibling sib_ijk into sib */
void amr_set_sibling_elm0(const tElm *elm, int sib_ijk, tElm *sib)
{
  int l = elm->eploc->l;
  tEloc s_eloc[1];

  /* copy header from elm into sib to get dt, time and such */
  memset(sib, 0, sizeof(tElm)); // 1st zero sib
  memcpy(sib, elm, sizeof(tElm0));

  /* set sib eploc */
  eloc_from_eploc(s_eloc, elm->eploc);
  if(l<1) errorexit("root node has no siblings");
  s_eloc->loc[l-1] = sib_ijk + '0';
  eloc_to_eploc(s_eloc, sib->eploc);

  /* set sib bbox */
  amr_set_elm_bbox(sib);

  /* set datrank to -1 since this is just a copy */
  sib->datrank = -1;
}

/****************************************************************************/
/* functions for elm-name and elm-location strings */
/****************************************************************************/

/* construct a unique string that describes node location in patch:
   e.g. 743 in octal is node on level 3 that has
    ijk=7 on l1,  ijk=4 on l2,  ijk=3 on l3   */
char *elm_location_str(tElm *elm, char *s, int slen)
{
  int i, l;
  tEloc eloc[1];
  eloc_from_eploc(eloc, elm->eploc);
  l = eloc->l;

  if(slen <= l) errorexit("slen is too small");

  for(i=0; i<l; i++)
    s[i] = eloc->loc[i];
  s[l] = 0;

  return s;
}

/* convert string from elm_location_str into a unsigned long int */
ulong elm_location__old(tElm *elm)
{
  ulong uloc;

  if(!elm) return 0;

  memcpy(&uloc, &(elm->eploc->ploc[0]), sizeof(uloc));

  /* remove MSB from uloc */
  uloc = uloc<<1;
  uloc = uloc>>1;

  return uloc;
}

/* use elm_location_str to make a unique node name that also contains the
   patch number */
char *elmname(tElm *elm, char *s, int slen)
{
  char loc[NLOCS+2];
  if(elm)
  {
    elm_location_str(elm, loc,NLOCS+2);
    snprintf(s,slen, "%d_%s", elm->eploc->p, loc);
  }
  else
  {
    snprintf(s,slen, "-");
  }
  return s;
}

/* get eploc of elm in a patch from string produced by elm_location_str */
void eloc_from_location_str(tEloc *eloc, int p, char *loc)
{
  int i;

  for(i=0; (i<NLOCS) && (loc[i]!=0); i++)
    eloc->loc[i] = loc[i];

  eloc->l = i;
  eloc->p = p;
}

/* get eploc of elm in a patch from string produced by elm_location_str */
void eploc_from_location_str(tEploc *eploc, int p, char *loc)
{
  tEloc eloc[1];
  eloc_from_location_str(eloc, p, loc);
  eloc_to_eploc(eloc, eploc);
}

/* get eloc of elm from full elmname */
void eloc_from_elmname(tEloc *eloc, char *name)
{
  int i, p;
  char *loc;
  int max = 99+NLOCS;

  /* find pos i of '_' */
  for(i=0; i<max; i++) if(name[i]=='_') break;

  /* get patch */
  p = atoi(name); /* atoi ignores '_' and all after it */
  //printf("name=%s => p=%d\n", name, p);

  /* get location str. */
  loc = name + i+1;
  eloc_from_location_str(eloc, p, loc);
}

/* get eploc of elm from full elmname */
void eploc_from_elmname(tEploc *eploc, char *name)
{
  tEloc eloc[1];
  eloc_from_elmname(eloc, name);
  eloc_to_eploc(eloc, eploc);
}

/* find an eloc in the elm list of all ranks,
   In: mesh, eloc   Out: eid
   Returns: the elm if it is on my rank, otherwise NULL */
tElm *amr_elm_eid_from_eloc(tMesh *mesh, tEloc *eloc, ulong *eid)
{
  int size = nMPI_size();
  int rank = nMPI_rank();
  int r, sum;
  char *found = checked_calloc(size, sizeof(found[0]));
  char *Found = checked_calloc(size, sizeof(found[0]));
  tElm **f_elm;

  f_elm = amr_elmarray_bsearch_eloc(mesh->nmyelm, mesh->myelm, eloc);
  if(f_elm) found[rank] = 1;
  Found[rank] = found[rank];
  nMPI_Allreduce(found, Found, size, nMPI_CHAR, nMPI_BOR);

  /* check result */
  for(sum=0, r=0; r<size; r++) sum += Found[r];
  if(sum>1) errorexit("at most one rank should have this eloc");

  /* get eid */
  *eid = EID_INVALID;
  if(f_elm) *eid = Elm_eid(f_elm[0]);

  /* get rank r that has eid and Bcast it */
  for(r=0; r<size; r++) if(Found[r]) break;
  nMPI_Bcast(eid,1, nMPI_UNSIGNED_LONG, r);

  free(Found);
  free(found);
  if(f_elm)  return f_elm[0];
  else       return NULL;
}

/* get elm,eid in the mesh from its full elmname,
   returns NULL if elm is on other rank */
tElm *elm_eid_from_elmname(tMesh *mesh, char *name, ulong *eid)
{
  tEloc eloc[1];
  eloc_from_elmname(eloc, name);
  return amr_elm_eid_from_eloc(mesh, eloc, eid);
}

/* get elm in the mesh from its full elmname
   WARNING: better use elm_eid_from_elmname */
tNode *elm_from_elmname(tMesh *mesh,  char *name)
{
  ulong eid;
  return elm_eid_from_elmname(mesh, name, &eid);
}

/* check if a elm has the name in string nname */
int elmname_is(tNode *elm, const char *nname)
{
  char myname[NLOCS+12];
  elmname(elm, myname,NLOCS+12); /* get name of elm into myname */
  if(strcmp(myname, nname)==0)   /* if myname=nname */
    return 1;
  else
    return 0;
}

/* write elocname into string s */
char *elocname(tEloc *eloc, char *s, int slen)
{
  int l = eloc->l;
  int l1, i;

  /* patch part */
  snprintf(s,slen, "%d_", eloc->p);
  l1 = strlen(s);
  //printf("|l=%d|1s=%s|", l,s);

  /* copy loc */
  for(i=0; (i<l) && (l1+i<slen-1); i++)
    s[l1+i] = eloc->loc[i];
  //printf("|2s=%s|", s);

  /* add str-end marker */
  if(l1+i<slen) s[l1+i] = 0;
  //printf("|l1+i=%d|3s=%s|", l1+i, s);

  return s;
}

/* check if a eloc has the name in string ename */
int elocname_is(tEloc *eloc, const char *ename)
{
  char myname[NLOCS+12];
  elocname(eloc, myname,NLOCS+12); /* get name of eloc into myname */
  if(strcmp(myname, ename)==0)     /* if myname=ename */
    return 1;
  else
    return 0;
}


//instead of elm_from_elmname and elm_from_location_str use:
/* get elmindex and elmrank for eid
   In: mesh, eid   Out: elmindex, elmrank
   Returns: elm if it is on my rank, otherwise NULL */
tElm *elm_from_eid(tMesh *mesh, ulong eid, ulong *elmindex, int *elmrank)
{
  amr_elmindex_and_elmrank_of_eid(mesh, eid, elmindex, elmrank);
  //PRF;printf(": eid=%lu\n", eid);
  //printf("*elmindex=%lu *elmrank=%d\n", *elmindex, *elmrank);

  if( nMPI_rank() == *elmrank ) return mesh->myelm[*elmindex];
  else                          return NULL;
}


/****************************************************************************/
/* functions that work on elm */
/****************************************************************************/

/* get ijk of elm */
int elm_get_ijk(tElm *elm)
{
  tEploc *eploc = elm->eploc;
  tEloc eloc[1];
  eloc_from_eploc(eloc, eploc);
  return connections_get_ijk(eloc->l, eloc->loc);
}

/* set patface of elm, return number of patfaces elm is on */
int elm_set_patface(tElm *elm, int patface[6])
{
  int npatfaces;
  tEploc *eploc = elm->eploc;
  tEloc eloc[1];
  eloc_from_eploc(eloc, eploc);
  npatfaces = connections_loc_on_patchface(eloc->l, eloc->loc, patface);
  return npatfaces;
}

/* returns: 1 if elm's face f in on patchface, 0 otherwise */
int elm_is_on_patface(tElm *elm, int f)
{
  tEploc *eploc = elm->eploc;
  tEloc eloc[1];
  eloc_from_eploc(eloc, eploc);
  return connections_loc_on_patchface_f(eloc->l, eloc->loc, f);
}


/* find node facenb in the node->fnb lists, returns 1  if it is there
   returns face and nb index in vars: face and ni */
int locate_facenb_in_fnbs(tNode *node, tNode *facenb, int *face, int *ni)
{
  int found = 0;
  int f;
  for(f=0; f<6; f++)
  {
    int nfnb = node->nfnb[f];
    int i;
    for(i=0; i<nfnb; i++)
    {
      if(node->fnb[f][i] == facenb)
      {
        found = 1;
        *face = f;
        *ni   = i;
        break;
      }
    }
    if(found) break;
  }
  return found;
}


/* get face of nb that touches elm,elmface
   Returns: nb-face nb_f if successful
            -1 otherwise */
int amr_get_nbface(tElm *elm, int elmface, tElm *nb)
{
  //int patface[6]; //, nfaces;
  tPat *pat = elm->pat;
  tBface *bface;
  int nb_f;
  const tEploc *eploc = elm->eploc;
  tEloc eloc[1];

  /* unpack eploc */
  eloc_from_eploc(eloc, eploc);

  /* check if elmface is on patch-face */
  //connections_loc_on_patchface(eloc->l, eloc->loc, patface);
  if(!connections_loc_on_patchface_f(eloc->l, eloc->loc, elmface))
  {
    /* simple case where elm and nb are in same patch */
    //PRF;printf(": nb is in same patch\n");
    nb_f  = elmface^1;
    return nb_f;
  }
  else
  {
    /* nb is outside elm's patch */
    //PRF;printf(": nb is outside patch\n");
    /* loop over all bfaces on face and find nb-face nb_f */
    forbfacesonface(pat, elmface, bface)
    {
      tBface *obface = bface->obface;
      int touch;

      /* do nothing if no other patch face */
      if(!obface) continue;

      /* check if we touch on the face of the other bface */
      nb_f = obface->f;
      touch = elm_common_facepoints(elm,elmface, nb,nb_f);
      if(touch) return nb_f;
    }
  }
  return -1;
}


/****************************************************************************/
/* functions to initialize tElm */
/****************************************************************************/

/* find patch of elm and save it in elm->pat */
void amr_set_elm_pat(tMesh *mesh, tElm *elm)
{
  tEploc *eploc = elm->eploc;
  int p = eploc->p;
  //elm->pat = elm->mesh->pat[p];
  elm->pat = mesh->pat[p];
}

/* find bbox of elmheader elm0 and save it in elm0->bbox */
void amr_set_elm0_bbox(tMesh* mesh, tElm0 *elm0)
{
  tEploc *eploc = elm0->eploc;
  int p = eploc->p; /* get patch number */
  int l = eploc->l; /* get level number */
  tPat *pat = mesh->pat[p];
  double *bbox  = elm0->bbox;
  double LX[3];
  int f, d, ll;
  tEloc eloc[1];
  char *loc;

  eloc_from_eploc(eloc, eploc);
  loc = eloc->loc;

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

/* find bbox of elm and save it in elm->bbox */
void amr_set_elm_bbox(tElm *elm)
{
  tMesh* mesh = elm->pat->mesh;
  union { tElm *elm; tElm0 *elm0; } e2e0;
  e2e0.elm = elm;
  amr_set_elm0_bbox(mesh, e2e0.elm0);
}


/* set eploc of child */
int amr_set_child_eploc(tEploc *parenteploc, int ijk, tEploc *eploc)
{
  tEloc eloc[1];
  int l = parenteploc->l;
  if(l >= NLOCS-1)
    errorexit("parenteploc is at limit ==> no further child possible!");
  eloc_from_eploc(eloc, parenteploc);
  eloc->p = parenteploc->p;
  eloc->l = l + 1;
  eloc->eid = EID_INVALID;
  eloc->loc[l]   = '0' + ijk;
  eloc->loc[l+1] = 0;
  eloc_to_eploc(eloc, eploc);
  return l+1;
}

/* set eploc of child */
int amr_set_parent_eploc(tEploc *eploc, tEploc *parenteploc)
{
  int l = eploc->l;
  if(l < 1)
    errorexit("eploc is root node, there is no parent");

  /* copy everything */
  parenteploc[0] = eploc[0];
  /* now change l and eid */
  parenteploc->l = l-1;
  parenteploc->eid = EID_INVALID;

  return l+1;
}

/****************************************************************************/
/* functions to initialize tElm0 */
/****************************************************************************/

/* compare two ulong numbers */
int cmp_ulong(const void *key, const void *ar, void *arg)
{
  ulong eid = *((const ulong *) key);
  ulong eidlim = *((const ulong *) ar);
  /* return eid - eidlim; // is not overflow safe */
  if(eid > eidlim) return 1;
  if(eid < eidlim) return -1;
  return 0;
}

/* return the rank that an elm with eid is on */
int amr_rank_of_eid(tMesh* mesh, ulong eid)
{
  int size = nMPI_size();
  ulong *eidlim = mesh->eidlim;
  const ulong *li;
  size_t off, num, nn;

  /* if eid is on rank rk: eidlim[rk-1] <= eid < eidlim[rk] */
  if(eid >= eidlim[size-1])
    errorexit("eid is outside the bounds of mesh->eidlim");

  if(eid < eidlim[0]) return 0; //rank0

  /* search in eidlim for eid */
  off = 0;
  num = size;
  li=bisectionsearch(&eid, eidlim, &off, &num, sizeof(eidlim[0]),
                     cmp_ulong, NULL);
  if(!li) errorexit("mesh->eidlim seems wrong");

  /* in case there are duplicates in eidlim do a linear search
     starting at off */
  //printf("off=%zu\n", off);
  for(nn=0; nn<8; nn++)
  {
    if(eidlim[off+nn] > eid) break;
  }
  //printf("nn=%zu\n", nn);
  //printf("off+1=%zu\n", off+1);

  return off+nn;
}

/* return the rank that an elm with eploc is on */
/*
int amr_rank_of_elm_eploc(tMesh* mesh, tEploc *eploc)
{
  ulong eid = eploc->eid;

  return amr_rank_of_eid(mesh, eid);
}
*/

/* return the number of elms on a rank */
ulong amr_nelms_on_rank(tMesh *mesh, int rank)
{
  ulong *eidlim = mesh->eidlim;
  if(rank>0) return eidlim[rank] - eidlim[rank-1];
  else       return eidlim[0];
}

/* return the eid of the first elm on a rank rank */
ulong amr_1st_eid_on_rank(tMesh *mesh, int rank)
{
  ulong *eidlim = mesh->eidlim;
  if(rank>0) return eidlim[rank-1];
  else       return 0;
}

/* find the index into mesh->myelm (on the rank it is on) that an elm
   with eid is on */
void amr_elmindex_and_elmrank_of_eid(tMesh* mesh, ulong eid,
                                     ulong *elmindex, int *elmrank)
{
  ulong *eidlim = mesh->eidlim;

  *elmrank = amr_rank_of_eid(mesh, eid);

  if(*elmrank==0) *elmindex = eid;
  else            *elmindex = eid - eidlim[*elmrank-1];
}

/* return the index into mesh->myelm (on the rank it is) that an elm
   with eploc is on */
/*
void amr_elmindex_and_elmrank_of_elm_eploc(tMesh* mesh, tEploc *eploc,
                                           ulong *elmindex, int *elmrank)
{
  ulong *eidlim = mesh->eidlim;
  ulong eid = eploc->eid;

  *elmrank = amr_rank_of_eid(mesh, eid);

  if(*elmrank==0) *elmindex = eid;
  else            *elmindex = eid - eidlim[*elmrank-1];
}
*/

/* init elm0 data as far as possible from eploc */
void amr_init_elm0_from_eploc(tMesh* mesh, tEploc *eploc, tElm0 *elm0)
{
  /* set eploc in elm0 */
  elm0->eploc[0] = eploc[0];

  /* set bbox and datrank, all else in elm0 remains unchanged */
  amr_set_elm0_bbox(mesh, elm0);
  elm0->datrank = amr_rank_of_eid(mesh, eploc->eid);
}


/****************************************************************************/
/* functions that sort and search elm lists contained in C-arrays */
/****************************************************************************/

/* sort elmarray with qsort */
void amr_elmarray_qsort(ulong narr, tElm **arr)
{
  qsort(arr, narr, sizeof(arr[0]), eecmp_q);
}

/* find elm with (eloc->p, eloc->l, eloc->loc) in elmarray arr
   (eloc->eid is not looked at) */
tElm **amr_elmarray_bsearch_eloc(ulong narr, tElm **arr, tEloc *eloc)
{
  tElm **f_elm;
  //PRFs(": ");printelmarray(narr, arr);
  f_elm = bsearch(eloc, arr, narr, sizeof(arr[0]), lecmp_q);
  //if(f_elm) { printf("found ");printeloc_s(eloc,"\n"); }
  //else      { printf("could not find ");printeloc_s(eloc,"\n"); }
  return f_elm;
}

/* find elm in elmarray arr */
tElm **amr_elmarray_bsearch(ulong narr, tElm **arr, tElm *elm)
{
  tElm **f_elm;
  tEloc eloc[1];
  eloc_from_eploc(eloc, elm->eploc); //could optimize if elm also has eloc
  //PRFs(": ");printelmarray(narr, arr);
  f_elm = bsearch(eloc, arr, narr, sizeof(arr[0]), lecmp_q);
  //if(f_elm) { printf("found ");printeloc_s(eloc,"\n"); }
  //else      { printf("could not find ");printeloc_s(eloc,"\n"); }
  return f_elm;
}

/* insert elm into sorted elmarray, modifies narrp, arrp */
void amr_elmarray_add_sort(ulong *narrp, tElm ***arrp, tElm *elm)
{
  ulong narr = *narrp;
  ulong narrp1 = narr+1; //len of extended newarr
  tElm **arr = *arrp;    //arr is elmarray
  tElm **newarr = realloc(arr, (narrp1)*sizeof(arr[0])); //make arr 1 bigger
  if(!newarr) errorexit("out of memory for newarr");
  newarr[narr] = elm; //add elm in last position
  amr_elmarray_qsort(narrp1, newarr);
  *narrp = narrp1; //increase narrp
  *arrp  = newarr; //point arrp to realloced mem.
}

tElm **amr_elmarray_linsearch(ulong narr, tElm **arr, tElm *elm)
{
  tElm **f_elm;
  ulong i;
  tEloc eloc[1];
  eloc_from_eploc(eloc, elm->eploc); //could optimize if elm also has eloc
  //PRFs(": ");printelmarray(narr, arr);
  f_elm = NULL;
  for(i=0; i<narr; i++)
  {
    if(lecmp_q(eloc, &(arr[i]))==0)
    {
      f_elm = &(arr[i]);
      break;
    }
  }
  //if(f_elm) { printf("found ");printeloc_s(eloc,"\n"); }
  //else      { printf("could not find ");printeloc_s(eloc,"\n"); }
  return f_elm;
}

/****************************************************************************/
/* functions that sort and search elm0 C-arrays */
/****************************************************************************/

/* return -1,0,1 if loc is before,at,after elem location,
   this is used in bsearch */
int le0cmp_q(const void *loc, const void *elem0)
{
  const tEloc *lc = loc;
  const tElm0 *elm0_arr = elem0;
  const tEploc *pelc = elm0_arr->eploc;
  tEloc elc[1];
  int cmp;

  /* NOTE: we could optimize this by caching an unpacked loc in each elm0: */
  eloc_from_eploc(elc, pelc);

  cmp = loccmp_mode(lc, elc, 1);
  //PRFs(": ");printeloc_s(lc, " ");printeloc_s(elc, " ");
  //printf("--> cmp=%d\n", cmp);
  ////printelm0(elm0_arr[0]);
  return cmp;
}

/* return -1,0,1 if key_elem location is before,at,after elem location,
   this is used in qsort */
int e0e0cmp_q(const void *key_elem0, const void *elem0)
{
  const tElm0 *kelm0 = key_elem0;
  const tEploc *keploc = kelm0->eploc;
  const tElm0 *elm0_arr = elem0;
  const tEploc *eploc = elm0_arr->eploc;
  tEloc klc[1], elc[1];
  int cmp;

  /* NOTE: we could optimize this by caching an unpacked loc in each elm0: */
  eloc_from_eploc(klc, keploc);
  eloc_from_eploc(elc, eploc);

  cmp = loccmp_mode(klc, elc, 1);
  //PRFs(": ");printeloc_s(klc, " ");printeloc_s(elc, " ");
  //printf("--> cmp=%d\n", cmp);
  ////printelm0(elm0_arr[0]);
  return cmp;
}

/* sort elm0 C-array with qsort */
void amr_elm0array_qsort(ulong narr, tElm0 *arr)
{
  qsort(arr, narr, sizeof(arr[0]), e0e0cmp_q);
}

/* find elm0 in a C-array of tElm0 */
tElm0 *amr_elm0array_bsearch(ulong narr, tElm0 *arr, tElm0 *elm0)
{
  tElm0 *f_elm0;
  tEloc eloc[1];
  eloc_from_eploc(eloc, elm0->eploc); //could optimize if elm0 also has eloc
  //PRF;for(int i=0; i<narr; i++) printeploc_s(arr[i].eploc, " ");
  //printf("\nbsearch: ");printeloc_s(eloc," ");
  f_elm0 = bsearch(eloc, arr, narr, sizeof(arr[0]), le0cmp_q);
  //if(f_elm0) printf("found\n");
  //else       printf("not found\n");
  return f_elm0;
}

/* insert elm0 into sorted elm0array, modifies narrp, arrp */
void amr_elm0array_add_sort(ulong *narrp, tElm0 **arrp, tElm0 *elm0)
{
  ulong narr = *narrp;
  ulong narrp1 = narr+1; //len of extended newarr
  tElm0 *arr = *arrp;    //arr is elm0array
  tElm0 *newarr = realloc(arr, (narrp1)*sizeof(arr[0])); //make arr 1 bigger
  if(!newarr) errorexit("out of memory for newarr");
  newarr[narr] = elm0[0]; //add content of elm0 in last position
  //PRFs(":1: ");
  //for(int i=0; i<narrp1; i++) printeploc_s(newarr[i].eploc, " ");
  //printf("\n");
  amr_elm0array_qsort(narrp1, newarr);
  //PRFs(":2: ");
  //for(int i=0; i<narrp1; i++) printeploc_s(newarr[i].eploc, " ");
  //printf("\n");
  *narrp = narrp1; //increase narrp
  *arrp  = newarr; //point arrp to realloced mem.
}

/* Add elm0 into sorted elm0array if it is not there yet.
   This modifies narrp, arrp. */
void amr_elm0array_unionadd_sort(ulong *narrp, tElm0 **arrp, tElm0 *elm0)
{
  tElm0 *f_elm0 = amr_elm0array_bsearch(*narrp, *arrp, elm0);
  if(!f_elm0)
    amr_elm0array_add_sort(narrp, arrp, elm0);
}




/****************************************************************************/
/* functions that use tEloc or tEploc to find and send elms */
/****************************************************************************/

/* Look in elm-array arr (in [arr+off,arr+num-1]) to find the elm
   with loc s_eloc and face s_f.
   *s_eloc is the loc caller wants to search for
   *But we start searching 1st for s_eloc's ancestor on level l0.
   *If s_eloc has children we increase the level number l
   Out: list f_elms_head of elms with loc s_eloc and face s_f
   Returns: On success: level number of descendant(s) of s_eloc
   Returns: On failure: an int below -999
   Note: elems in f_elms_head have to be freed later! */
int amr_make_elms_on_eloc_face_list(long narr, tElm **arr,
                                    size_t off0, size_t num0,
                                    tEloc s_eloc[1], int s_f, int l0,
                                    struct list_head *f_elms_head)
{
  tElm *const*f_elm;
  int f_elm_l;
  size_t off, num;
  //tEloc *eloc = elm->eloc;
  tEloc f_eloc[1];
  tEloc cheloc[1];
  int l, lret;
  int lmax = s_eloc->l;
  //unsigned mor, atB;
  int ijk;

  /* We start searching 1st for s_eloc's ancestor on level l0,
     thus l0 cannot be bigger than lmax s_eloc->l. */
  if(l0>lmax) errorexit("l0>lmax");

  /* init */
  off = off0;
  num = num0;

  //PRF;printf(": off=%zu num=%zu l0=%d s_eloc=", off, num, l0);
  //printeloc(s_eloc);
  //printf(" s_f=%d\n  ", s_f);
  ////printelmarray(narr, arr);

  for(l = l0; l <= lmax; l++)
  {
    /* search for ancestor of s_eloc of level l */
    f_eloc[0] = s_eloc[0];
    f_eloc->l = l;
    f_elm = binarysearch(f_eloc, arr, &off, &num, sizeof(*arr), lecmp_0, NULL);

    //printf("off=%zu num=%zu  f_elm=%p <-pos=%zu\n",
    //       off, num, f_elm, (size_t) ((tElm **)f_elm - arr));
    //if(f_elm) { printf("got ");printelm(*f_elm); }

    if(!f_elm) return -s_eloc->l - 1000; /* found nothing */

    /* if we find an exact match we are done */
    f_elm_l = Elm_l(*f_elm);
    if(f_elm_l == l)
    {
      /* looks like we found a unique elm, so add f_elm to list and then
         return with success */
      glist_entry_add_tail(*f_elm, f_elms_head);
      return l;
    }

    ///* if f_elm is longer than what we were looking for,
    //   we need to try again with a larger f_eloc->l = l */
    //if(f_elm_l > l) continue;


    ///* is there only one f_elm? */
    //binarysearchmore(f_eloc, arr, narr, sizeof(*arr), f_elm, lecmp_0, NULL,
    //                 &mor, &atB);
    //printf("mor=%d\n", mor);
    //if(!mor)  /* if there is only one */
    //{
    //  f_elm_l = Elm_l(*f_elm);
    //  /* Note: If we are at arr-bound (atB) there could be more than one.
    //           If f_elm_l > s_eloc->l and we are at arr-bound there are
    //           likely more on other ranks and we need to further verify the
    //           one we found */
    //  if( (f_elm_l <= s_eloc->l) || (!atB) )
    //  {
    //    /* looks like we found a unique elm, so add f_elm to list and then
    //       return with success */
    //    glist_entry_add_tail(*f_elm, f_elms_head);
    //    return l;
    //  }
    //}
  }
  /* if we get here, nb at s_eloc has children */

  /* search on children one level higher */
  l = s_eloc->l + 1;
  cheloc[0] = s_eloc[0];
  cheloc->l = l;
  /* we set cheloc->loc below */

  /* default return value */
  lret = -l - 1000; /* found nothing */

  /* get the 4 children elocs on nb face s_f */
  for(ijk = 0; ijk<8; ijk++) /* loop over all children */
  {
    if(connections_ijk_is_at_parentface(ijk, s_f)) /* only 4 are relevant */
    {
      int ret;
      /* child ijk */
      cheloc->loc[l-1] = ijk + '0';
      if(l<NLOCS) cheloc->loc[l] = 0;
      ret = amr_make_elms_on_eloc_face_list(narr, arr, off, num, cheloc, s_f,
                                            l, f_elms_head);
      if(ret >=0) lret = l; /* record success */
    }
  }

  /* finally signal failure or success with at least one nb child */
  //printf("final lret=%d\n", lret);
  return lret;
}


/* find a neighbors on patchface of elm,elmface in elmarray narr,arr
   In:  elm,elmface, narr,arr
   Out: will append to list fnb_head
   Note: elems in fnb_head have to be freed later! */
int amr_make_patchface_fnb_list(tElm *elm, int elmface,
                                long narr, tElm **arr,
                                struct list_head *fnb_head)
{
  //int patface[6]; //, nfaces;
  int l  = - 9999; /* found nothing */
  int l2 = l;
  tPat *pat = elm->pat;
  tBface *bface;
  tEloc nbeloc[1];
  int nb_f;
  struct list_head *pos, *sav;
  const tEploc *eploc = elm->eploc;
  tEloc eloc[1];

  /* unpack eploc */
  eloc_from_eploc(eloc, eploc);
  //PRFs(": ");printeloc(eloc);printf(" elmface=%d\n", elmface);

  /* sanity check */
  if(!connections_loc_on_patchface_f(eloc->l, eloc->loc, elmface))
    errorexit("call this for patch faces only!!!");

  //printbfaces_on_f(pat, elmface);

  /* loop over all bfaces on face and find nb */
  forbfacesonface(pat, elmface, bface)
  {
    tBface *obface = bface->obface;

    //printf("bface\n");
    //printbface(bface);
    //printf("obface\n");
    //printbface(obface);

    /* do nothing if no other patch face */
    if(!obface) continue;

    /* eloc and face of root elm in other patch */
    nbeloc->l      = 0;
    nbeloc->loc[0] = 0;
    nbeloc->p      = obface->pat->p;
    nb_f = obface->f;

    /* set pos to last entry to start of list_for_each_safe_continue below */
    pos = fnb_head->prev;

    //printeloc(nbeloc);printf(" nb_f=%d\n", nb_f);
    /* add all elms in arr on face nb_f to list fnb_head */
    l2=amr_make_elms_on_eloc_face_list(narr, arr, 0, narr, nbeloc, nb_f,
                                       0, fnb_head);

    /* Go over newly added part of fnb_head list and remove all that have no
       face points in common with the elm */
    list_for_each_safe_continue(pos, sav, fnb_head)
    {
      /* get neigh. and check if elm and nb have common points */
      tGlist *elem = list_entry(pos, tGlist, list);
      tElm *nb = elem->entry;
      int touch = elm_common_facepoints(elm,elmface, nb,nb_f);

      if(!touch) /* remove elem with nb from list fnb_head */
        glist_elem_del(elem);
      else
        l = fmax(l,l2);
    }
  } /* end forbfacesonface */

  //PRFs(": fnb list\n");
  //list_for_each(pos, fnb_head)
  //{
  //  tElm *fnb = glist_entry(pos);
  //  tEploc *fnbeploc = fnb->eploc;
  //  printeploc_s(fnbeploc, " ");
  //}
  //printf("\n");

  return l;
}



/* add a list of nbeplocs to the var amr_elm_nbinfo on one face */
void amr_elm_nbinfo_add_nbeploc(tElm *elm, int face,
                                int nnb, const tEploc nbeploc[nnb])
{
  int i_nbinfo = amr->elm_nbinfo0 + face;
  tArray *nbinfo;

  /* -if nbinfo is enabled we assume it has some nbs already
      and that we now want to add nbeploc
     -if it is not enabled we assume it has been cleared before and we
      just set it to nbeploc */
  nbinfo = VarA(elm, i_nbinfo);
  //if(elmname_is(elm, "2_7") && face==2)
  //{
  //  PRF;printeploc(elm->eploc);
  //  printf(": 2_7 face%d\n", face);
  //  printf("nbinfo1: ");print_amr_elm_nbinfo(elm,face);printf("\n");
  //  //printf("nbinfo2: ");printarray_eploc(nbinfo, 1);//printf("\n");
  //  printf("nnb=%d nbeploc: ", nnb);
  //  for(int ni=0; ni<nnb; ni++) printeploc_s(&(nbeploc[ni]), " ");
  //  printf("\n");
  //}
  if(!nbinfo) /* no nb info yet at all */
  {
    /* switch on var */
    enablevarcomp_innode(elm, i_nbinfo);
    nbinfo = VarA(elm, i_nbinfo);

    /* write all of nbeploc into var amr_elm_nbinfo */
    memcpy_to_array_redim(nbinfo, sizeof(tEploc), 0,
                          nbeploc, sizeof(nbeploc[0])*nnb);
  }
  else /* combine info in amr_elm_nbinfo and info in nnb,nbeploc */
  {
    int nbinfo_nnb = array_Neplocs(nbinfo); //num. of nbs we already have
    int ni;

    /* loop over nbs to check if they had been added earlier already,
       this should never happen! */
    for(ni=0; ni<nnb; ni++)
    {
      int ret = eploc_key_in_eplocarray(&(nbeploc[ni]),
                                        nbinfo_nnb, nbinfo->eploc);
      //printf("ret=%d\n", ret);
      if(ret>=0)
        errorexiti("nbeploc[%d] is already in nbinfo->eploc", ret);
    }

    /* add all of nbeploc to var amr_elm_nbinfo */
    memcpy_to_array_redim(nbinfo, sizeof(tEploc), nbinfo_nnb,
                          nbeploc, sizeof(nbeploc[0])*nnb);
  }
}

/* set trivial nb info between siblings */
void amr_set_intersibling_nbinfo_nnbinfo(tElm *sib0)
{
  tElm *sib[8];
  int ijk;
  struct list_head *pos_ijk = &sib0->list; /* pos of 1st sib */

  /* record 8 sibs in sib[8] array */
  for(ijk=0; ijk<8; ijk++)
  {
    sib[ijk] = list_entry(pos_ijk, tElm, list);
    pos_ijk = pos_ijk->next;  /* pos of next sib */
  }

  /* set inter-sib nbinfo for each sib */
  for(ijk=0; ijk<8; ijk++)
  {
    tElm *sib_ijk = sib[ijk];
    int f;
    for(f=0; f<6; f++)
    {
      tElm *nb;
      if(!connections_ijk_is_at_parentface(ijk, f))
      {
        int nb_ijk = connections_get_inner_nb_ijk(ijk, f/2);
        nb =  sib[nb_ijk];

        /* set nbinfo and nnbinfo */
        disablevarcomp_innode(sib_ijk, amr->elm_nbinfo0+f);
        amr_elm_nbinfo_add_nbeploc(sib_ijk, f, 1, nb->eploc);
        sib_ijk->dat->info->nnbinfo[f] = 1;

        /* now also set fnb pointers */
        free(sib_ijk->fnb[f]);
        sib_ijk->fnb[f] = checked_calloc(1, sizeof(sib_ijk->fnb[f][0]));
        sib_ijk->fnb[f][0] = nb;
        sib_ijk->nfnb[f] = 1;
        //printelm(sib_ijk);
      }
    } /* end f-loop */
  }
}


/* write number of eplocs in var amr_elm_nbinfo into
   elm->dat->info->nnbinfo */
void amr_elm_nbinfo_set_nnbinfo(tElm *elm, int positive)
{
  tDat *dat = elm->dat;
  int f;

  if(!dat) return; /* do nothing if no dat */

  for(f=0; f<6; f++)
  {
    int i_nbinfo = amr->elm_nbinfo0 + f;
    tArray *nbinfo = VarA(elm, i_nbinfo);
    int nnbinfo_new;
    int nnbinfo_old = elm->dat->info->nnbinfo[f];

    /* if there is no nbinfo there are no nbs */
    /* is nbinfo enabled? */
    if(!nbinfo)
      nnbinfo_new = 0; /* if there is no nbinfo there are no nbs */
    else
      nnbinfo_new = array_Neplocs(nbinfo); /* num. of nbs we have */

    /* write back into nnbinfo[f] */
    /* if it was non-negative we keep it so */
    if(nnbinfo_old >= 0 || positive)
      elm->dat->info->nnbinfo[f] = nnbinfo_new;
    else /* we keep it negative */
        elm->dat->info->nnbinfo[f] = -nnbinfo_new - 1;
  }
}

/* call amr_elm_nbinfo_set_nnbinfo for all elms in mesh */
int amr_elm_nbinfo_set_nnbinfo_mesh(tMesh *mesh, int positive)
{
  formyelms(mesh)
  {
    tElm *elm = MyElm;
    amr_elm_nbinfo_set_nnbinfo(elm, positive);
  }
  return 0;
}

/* redim var amr_elm_nbinfo according to elm->dat->info->nnbinfo */
void amr_elm_nbinfo_redim_according_to_nnbinfo(tElm *elm)
{
  int f;
  for(f=0; f<6; f++)
  {
    int i_nbinfo = amr->elm_nbinfo0 + f;
    tArray *nbinfo = VarA(elm, i_nbinfo);
    int nnb;

    /* if there is no nbinfo there should be no nbs */
    if(nbinfo==NULL)
    {
      if(elm->dat->info->nnbinfo[f] > 0) //check if there should be nbs
        errorexit("nbinfo=NULL contradicts nnbinfo>0");
      else        //nnbinfo[f]=0 is consistent, nnbinfo[f]<0 means invalid
        continue; //for nnbinfo[f]<=0 we do nothing
    }

    nnb = elm->dat->info->nnbinfo[f];
    /* since a need to update is signaled by nnbinfo -> -nnbinfo-1 */
    if(nnb<0) nnb = -nnb-1;
    redim_array_Neplocs(nbinfo, nnb);
    //PRF;printf(": eid%lu f%d nnb=%d ", Elm_eid(elm), f, nnb);
    //printf("  ");print_amr_elm_nbinfo(elm, f);printf("\n");
  }
}

/* update eid part of eplocs in nbinfo vars */
void amr_elm_nbinfo_update_eid_locally_using_fnb(tElm *elm)
{
  int f;
  for(f=0; f<6; f++)
  {
    int i_nbinfo = amr->elm_nbinfo0 + f;
    tArray *nbinfo = VarA(elm, i_nbinfo);
    int nfnb = elm->nfnb[f];
    tElm **fnb = elm->fnb[f];
    int nnb, Neplocs, ni;
    //PRFs(":\n");
    //printf("nbinfo=%p fnb=%p", nbinfo, fnb);

    /* if there is no nbinfo or if it is outdated do nothing */
    if(!nbinfo) continue;
    nnb = elm->dat->info->nnbinfo[f];
    //NOTE and CHECK: We had this:
    //if(nnb <= 0) continue;
    if(nnb<0) nnb = -nnb-1;
    Neplocs = array_Neplocs(nbinfo);
    if(nnb!=Neplocs)
      errorexiti("amr_elm_nbinfo and nnbinfo disagree on face%d", f);

    //PRFs(": ");printelm(elm);
    //printf(" f%d Neplocs=%d nfnb=%d\n", f, Neplocs, nfnb);

    if( (Neplocs!=nfnb) || (!fnb && Neplocs>0) )
      errorexiti("nbinfo and fnb disagree on face%d", f);

    /* set eid in each amr_elm_nbinfo entry to the actual eid of the nb */
    for(ni=0; ni<Neplocs; ni++)
    {
      tElm *nb = fnb[ni];
      tEploc *eploc = &(nbinfo->eploc[ni]);
      if(nb) eploc->eid = Elm_eid(nb);
    }
  }
}
/* call amr_elm_nbinfo_update_eid_locally_using_fnb on entire new
   mesh->myelm */
void amr_elm_nbinfo_update_eid_locally_using_fnb_mesh(tMesh *mesh)
{
  formyelms(mesh)
  {
    tElm *elm = MyElm;
    amr_elm_nbinfo_update_eid_locally_using_fnb(elm);
  }
}


/* Look in elm-array arr (in [arr+off,arr+num-1]) to find the nb of
   elm on face elmface.
   In: elm,elmface, narr,arr  =>  Out: fnb_head, Returns list length
   Note: elems in fnb_head have to be freed later! */
int amr_make_fnb_list(tElm *elm, int elmface, long narr, tElm **arr,
                      struct list_head *fnb_head)
{
  //int patface[6];
  const tEploc *eploc = elm->eploc;
  tEloc eloc[1];

  eloc_from_eploc(eloc, eploc);
  //if(elmname_is(elm, "1_7124"))// && elmface==1)
  //{ PRFs(": ");printeloc(eloc);printf(" elmface=%d\n", elmface); }

  /* simple case where elmface is inside the patch */
  if(!connections_loc_on_patchface_f(eloc->l, eloc->loc, elmface))
  {
    tEloc nbeloc[1];
    int nbface, l0;
    nbeloc->p = eloc->p;
    l0 = connections_get_nbloc_InsidePat(eloc->l, eloc->loc, elmface,
                                         nbeloc->loc, &nbface);
    nbeloc->l = eloc->l;
    //if(elocname_is(nbeloc, "1_7124"))// && elmface==1)
    //{
    //  printf(" -> l0=%d nbeloc=", l0);printeloc(nbeloc);printf(" nbface=%d\n", nbface);
    //  printf("nbeloc->l=%d\n", nbeloc->l);
    //  printf("find nbeloc,nbface in: ");printelmarray(narr, arr);
    //}
    //if(elmname_is(elm, "1_7124"))// && elmface==1)
    //{
    //  char s[100];
    //  printf(" -> l0=%d nbeloc=%s", l0, elocname(nbeloc,s,100));printf(" nbface=%d\n", nbface);
    //  printf("nbeloc->l=%d\n", nbeloc->l);
    //  printf("find nbeloc,nbface in: ");printelmarray(narr, arr);
    //}
    /* Look for nbs at level l0.
       l0 = max(nb_l-4, 0) would be ok if we allow for only a
       difference of 4 in refinement of nbs... */
    //l0 = 0; //For now we set l0 = 0.
    /* but the best may be to use the output we got above ??? */
    amr_make_elms_on_eloc_face_list(narr, arr, 0, narr, nbeloc, nbface,
                                    l0, fnb_head);
  }
  else /* complicated case where elmface is on patch surface */
  {
    amr_make_patchface_fnb_list(elm,elmface, narr,arr, fnb_head);
  }

  //if(elmname_is(elm, "1_7124"))// && elmface==1)
  //{
  //  PRF;printf(": myrank=%d elmface=%d patface[elmface]=%d ",
  //             nMPI_rank(), elmface, connections_loc_on_patchface_f(eloc->l, eloc->loc, elmface));
  //  printeploc_s(elm->eploc, " \n");
  //  printelmglist(fnb_head);
  //}
  //  //printf("on=%d\n", connections_loc_on_patchface_f(0,"530", 2));
  //  //printf("on=%d\n", connections_loc_on_patchface_f(1,"530", 1));
  //  //printf("on=%d\n", connections_loc_on_patchface_f(2,"530", 1));
  //  //printf("on=%d\n", connections_loc_on_patchface_f(3,"530", 1));
  //  //printf("on=%d\n", connections_loc_on_patchface_f(4,"5305", 1));
  //  //printf("on=%d\n", connections_loc_on_patchface_f(8,"53351537", 1));
  //  //int pf[6];
  //  //printf("on=%d\n", connections_loc_on_patchface(0,"530", pf));
  //  //printf("on=%d\n", connections_loc_on_patchface(1,"530", pf));
  //  //printf("on=%d\n", connections_loc_on_patchface(2,"530", pf));
  //  //printf("on=%d\n", connections_loc_on_patchface(3,"530", pf));
  //  //printf("on=%d\n", connections_loc_on_patchface(4,"5305", pf));
  //  //printf("on=%d\n", connections_loc_on_patchface(8,"53351537", pf));
  //  errorexit("stop");
  //}
  return list_count_nodes(fnb_head);
}


/* Update amr_elm_nbinfo vars on all faces where:
   elm->dat->info->nnbinfo[f] < 0.
   This is done for all elms on all ranks */
/* NOTE: amr_update_elm_nbinfo_if_nnbinfo_negative is very general as it
   assumes no prior knowledge of which MPI rank I am actually in contact with.
   BUT amr_update_elm_nbinfo_if_nnbinfo_negative is very slow!!!
   This happens because every rank sends its ef0 to ALL other ranks:
     nMPI_Bcast(ef0, nef0[f], nMPIvars->TELM0, rk);
   So we should keep this func, but also make an improved version, where:
   +I should only send to the ranks that I am in contact with, i.e. the ones
    I had nbs with on this face before. We could destill a contacts list for
    each face-dir out of mesh->nbelms and store it in mesh.
   +Or MAYBE we should save a contacts list for each face in each elm. Then
    we can do local updates (similar to the surface exchanges). */
int amr_update_elm_nbinfo_if_nnbinfo_negative(tMesh *mesh)
{
  struct list_head *pos;
  struct list_head ef0_head[6]; // one list for each face
  ulong nmyef0[6] = {0};        // total number of elms
  int f, rk;
  int rank=nMPI_rank();
  int size=nMPI_size();

  /* init */
  for(f=0; f<6; f++) { INIT_LIST_HEAD(&ef0_head[f]); nmyef0[f]=0; }

  /* find all my elmfaces that have nnbinfo<0 (all ranks do this) */
  list_for_each(pos, &mesh->myelm_head)
  {
    tElm *elm = list_entry(pos, tElm, list);
    tDat *dat = elm->dat;

    /* go over elm-faces */
    if(dat)
    {
      tNodeInfo *info = elm->dat->info;
      for(f=0; f<6; f++)
      {
        /* if nnbinfo<0 nb info is not there yet */
        if(info->nnbinfo[f] < 0)
        {
          /* erase all nb info in var amr_elm_nbinfo[f] */
          disablevarcomp_innode(elm, amr->elm_nbinfo0+f);
          /* add elm to list of elms that need nb-info */
          glist_entry_add_tail(elm, &ef0_head[f]);
          nmyef0[f] += 1;
        }
      }
    }
  }
  /* now we have 6 lists ef0_head[f] that contain elms where the fnb info
     needs to be updated. */

  /* print what we have so far */
  /*
  for(f=0; f<6; f++)
  {
    PRF;printf(": nmyef0[f]=%lu  &ef0_head[%d]:\n", nmyef0[f], f);
    list_for_each(pos, &ef0_head[f])
    {
      tElm *elm = glist_entry(pos);
      printeploc_s(elm->eploc, " ");
    }
    printf("\n");
  }
  */

  //printf("XXXXXX count=%lu\n", list_count_nodes(&ef0_head[0]));
  //printf("XXXXXX rank%d nmyef0[0]=%lu\n", rank, nmyef0[0]);

  /* send my lists to the other ranks */
  for(rk=0; rk<size; rk++)
  {
    ulong nef0[6];
    /* ef0_nbs is a large array, where we will store all nbs of all the
       nef0[f] elms rank rk needs nb info about, for all faces f. Layout is:
        ef0_nbs = |nef0[0]|nnb0|nb_eploc[0...nnb0-1]|
                          |nnb1|nb_eploc[0...nnb1-1]| <--all entries have
                          ...                            sizeof(tEploc) bytes
                  |nef0[5]|nnb0|nb_eploc[0...nnb0-1]|
                          |nnb1|nb_eploc[0...nnb1-1]|
                          ... */
    tArray *ef0_nbs = alloc_array1d((18*sizeof(tEploc))/sizeof(double));
    //tArray *ef0_nbs = alloc_array1d((1*sizeof(tEploc))/sizeof(double));
    ulong ef0_nbs_idx = 0; /* index of next entry to add */
    ulong nmyEplocs;       /* number of tEploc sized entries in ef0_nbs */

    //printf("YYYYYYY rank%d nmyef0[0]=%lu\n", rank, nmyef0[0]);

    /* rank rk copies his nmyef0 into nef0 */
    if(rank == rk)
      for(f=0; f<6; f++) nef0[f] = nmyef0[f];

    /* rank rk sends his nef0[f] to all others, to tell how many elms he
       wants to find face nbs of */
    //printf("111111 rank%d rk%d nef0[0]=%lu\n", rank, rk, nef0[0]);
    nMPI_Bcast(&nef0[0],6, nMPI_UNSIGNED_LONG, rk);
    //printf("222222 rank%d rk%d nef0[0]=%lu\n", rank, rk, nef0[0]);

    /* init ef0_nbs index counter */
    ef0_nbs_idx = 0;

    /* make list for each face and send them... */
    for(f=0; f<6; f++)
    {
      /* put the number nef0[f] into ef0_nbs array */
      memcpy_to_array_redim(ef0_nbs, sizeof(tEploc), ef0_nbs_idx,
                            &(nef0[f]), sizeof(nef0[f]));
      ef0_nbs_idx++;

      //printf("ADD nef0[f] int, next ef0_nbs_idx=%lu\n", ef0_nbs_idx);
      //printarray_eploc(ef0_nbs, 1);
      //exit(16);

      /* there is something to do only if nef0[f]>0 */
      if(nef0[f])
      {
        // get &ef0_head[f] from rk into elmhead or eploc array ef0
        //tElm0 *ef0 = checked_calloc(nef0[f], sizeof(ef0[0]));
        //hey maybe we should just send tEploc not tElm0!!!???
        tEploc *ef0 = checked_calloc(nef0[f], sizeof(ef0[0]));
        struct list_head *pos0;
        ulong i;

        /* rank rk now fills the ef0 array */
        if(rank == rk)
        {
          i=0;
          list_for_each(pos0, &ef0_head[f])
          {
            tElm *elm = glist_entry(pos0);
            //memcpy(&ef0[i], elm, sizeof(ef0[0]));
            memcpy(&ef0[i], elm->eploc, sizeof(ef0[0]));
            i++;
          }
          if(nef0[f]!=i) errorexit("nef0[f]!=i");
        }

        /* broadcast all elmheaders in ef0 from rank rk to all */
        //nMPI_Bcast(ef0, nef0[f], nMPIvars->TELM0, rk);
         /* broadcast all eplocs in ef0 from rank rk to all */
        nMPI_Bcast(ef0, nef0[f], nMPIvars->TEPLOC, rk);
                               //^^^^^^^^^^^^^^^-is this right???

        /* all ranks do work on ef0 array and find all nbs of all in ef0 */
        for(i=0; i<nef0[f]; i++)
        {
          //tElm *elmi = alloc_elm_of_elmheader(mesh, &ef0[i]);
          tElm *elmi = alloc_elm_of_eploc(mesh, &ef0[i]);
          struct list_head *pos1, *sav;
          struct list_head fnb_head;
          ulong j, nnb;

          INIT_LIST_HEAD(&fnb_head);

          /* put the nbs of elmi into fnb_head list */
          nnb = amr_make_fnb_list(elmi, f, mesh->nmyelm, mesh->myelm,
                                  &fnb_head);

          /* put the number nnb into ef0_nbs array */
          memcpy_to_array_redim(ef0_nbs, sizeof(tEploc), ef0_nbs_idx,
                                &(nnb), sizeof(nnb));
          ef0_nbs_idx++;

          /* get nb eploc into ef0_nbs array */
          j=0;
          list_for_each_safe(pos1, sav, &fnb_head)
          {
            tGlist *elem = list_entry(pos1, tGlist, list);
            tElm *nb = elem->entry;
            /* put nb->eploc into ef0_nbs array */
            memcpy_to_array_redim(ef0_nbs, sizeof(tEploc), ef0_nbs_idx,
                                  nb->eploc, sizeof(tEploc));
            ef0_nbs_idx++;
            j++;

            /* once nb->eploc is in ef0_nbs, del elem with nb */
            glist_elem_del(elem);
          }
          if(nnb!=j) errorexit("nnb!=j");

          //if(elmname_is(elmi, "2_7") && f==2)
          //{
          //  PRF;printf(": rank%d rk=%d: ", rank, rk);
          //  printeploc(elmi->eploc);
          //  printf(": 2_7 f%d\n", f);
          //  for(int ni=0; ni<nnb; ni++)
          //    printeploc_s(ef0_nbs->eploc+ef0_nbs_idx-j+ni, " ");
          //  printf("\n");
          //}

          /* now &fnb_head is freed, so just free the elmi */
          free_elm(elmi);
        }
        free(ef0);
      }
    } /* end loop over face f */

    /* number of tEploc sized entries in ef0_nbs */
    nmyEplocs = ef0_nbs_idx;
    //printf("nmyEplocs=%lu\n", nmyEplocs);
    //printarray_eploc(ef0_nbs, 1);

    if(rank != rk) /* send to rank rk */
    {
      /* first send number of tEploc sized entries in ef0_nbs */
      nMPI_Send(&nmyEplocs,1, nMPI_UNSIGNED_LONG, rk, 1000);

      /* now send contents of ef0_nbs */
      nMPI_Send(ef0_nbs->d,nmyEplocs, nMPIvars->TEPLOC, rk, 2000);

    }
    else /* rank=rk: i.e. I am rank rk and will revc from all others */
    {
      ulong *N_eplocs = checked_calloc(size, sizeof(N_eplocs[0]));
      tEploc **eplocs;
      int r;

      /* revc number of tEploc sized entries from each rank r */
      for(r=0; r<size; r++)
      {
        if(r != rk)
          nMPI_Recv(&N_eplocs[r],1, nMPI_UNSIGNED_LONG, r, 1000);
        else
          N_eplocs[r] = nmyEplocs;
      }

      /* make recv buffers for data that is recvd */
      eplocs = rows_calloc(size, N_eplocs, sizeof(eplocs[0][0]));
      /* but rank rk does not need a recv buffer, because it has all
         in ef0_nbs->d already, so we use that here */
      free(eplocs[rk]);
      /* transfer ef0_nbs->d from ef0_nbs to eplocs[rk] */
      eplocs[rk] = ef0_nbs->eploc; //eplocs[rk] is in my tArray ef0_nbs
      ef0_nbs->d_nofree=1;
      free_array(ef0_nbs);
      ef0_nbs = alloc_array1d(1); /* dummy that will be freed below */
      ef0_nbs_idx = 0;

      /* revc contents of ef0_nbs from each rank r */
      for(r=0; r<size; r++)
      {
        if(r != rk)
          nMPI_Recv(eplocs[r], N_eplocs[r], nMPIvars->TEPLOC, r, 2000);
        /* Note: eplocs[rk] already has what was in ef0_nbs->d before */
      }


      /**********************************************/
      /* read eplocs[r] to build nb info on rank rk */
      /**********************************************/
      {
        for(r=0; r<size; r++)
        {
          ulong epi;
          /* each eplocs[r] is a tEploc array, where we will store all nbs of
             all the nef0[f] elms rank rk needs nb info about, for all faces f.
             Layout is:
             eplocs[r] = |nelms[0]|nnb0|nb_eploc[0...nnb0-1]|
                                  |nnb1|nb_eploc[0...nnb1-1]|
                                  ...
                         |nelms[5]|nnb0|nb_eploc[0...nnb0-1]|
                                  |nnb1|nb_eploc[0...nnb1-1]|
                                  ... */
          epi = 0;
          for(f=0; f<6; f++)
          {
            union { tEploc e; ulong ul; } e2ul;
            struct list_head *pos1;
            ulong nelms, ei;

            /* pos of 1st elm in ef0_head[f] list */
            pos1 = ef0_head[f].next;

            /* get number of elms nelms out of eplocs[r] */
            e2ul.e = eplocs[r][epi++];
            nelms  = e2ul.ul;

            //printf("f=%d ef0_head[%d] count=%lu\n", f,f, list_count_nodes(&ef0_head[f]));
            //printf("f=%d: epi-1=%lu nelms=%lu ", f, epi-1, nelms);
            //printeploc_s(&(eplocs[r][epi-1]),"\n");

            for(ei=0; ei<nelms; ei++)
            {
              tGlist *elem;
              tElm *elm;
              ulong nnb;

              /* get elm from list ef0_head[f] */
              elem = list_entry(pos1, tGlist, list);
              elm  = elem->entry;
              pos1 = pos1->next; /* go forward now, because we del below */

              /* get number of nbs out of eplocs[r] */
              e2ul.e = eplocs[r][epi++];
              nnb    = e2ul.ul;

              //printeploc_s(elm->eploc, " ");
              //printf("f%d r=%d ", f, r);
              //printf("ei=%lu  nnb=%lu epi=%lu ", ei, nnb, epi);
              //if(nnb) printeploc_s(&(eplocs[r][epi]), " ...");
              //printf("\n");

              /* add all nbs in eplocs[r] to var amr_elm_nbinfo */
              amr_elm_nbinfo_add_nbeploc(elm, f, nnb, &(eplocs[r][epi]));
              epi += nnb;

              /* cannot del elem from ef0_head[f] list here, because it is
                 needed for every r */
            }
          } /* end for f */
        }
        /* now clear the 6 ef0_head[f] lists */
        for(f=0; f<6; f++) glist_free_elems(&(ef0_head[f]));
      } /* end func that builds nb-info from elocs[r] */

      /* the eplocs has been all read now, so free it */
      rows_free(eplocs, size);
      free(N_eplocs);
    }

    /* we could reuse the large tArray, but for now we just free it */
    free_array(ef0_nbs);
  } /* end loop over rk */


  /* finally set nnbinfo according to the new nb-info we have now,
     but we keep them negative for now */
  amr_elm_nbinfo_set_nnbinfo_mesh(mesh, 0);

  return 0;
}

/* clear amr_elm_nbinfo[f] if nnbinfo[f] < 0 */
void amr_disable_elm_nbinfo_if_nnbinfo_negative(tMesh *mesh)
{
  formyelms(mesh)
  {
    tElm *elm = MyElm;
    tDat *dat = elm->dat;
    if(dat)
    {
      tNodeInfo *info = elm->dat->info;
      int f;
      for(f=0; f<6; f++)  /* go over elm-faces */
        if(info->nnbinfo[f] < 0) /* if nnbinfo<0 nb info is not there yet */
        {
          /* erase all nb info in var amr_elm_nbinfo[f] */
          disablevarcomp_innode(elm, amr->elm_nbinfo0+f);
          info->nnbinfo[f] = -1;
        }
    }
  }
}

/* find all my elmfaces that have nnbinfo<0 */
void amr_update_elm_nbinfo_locally_if_nnbinfo_negative(tMesh *mesh)
{
  formyelms(mesh)
  {
    tElm *elm = MyElm;
    tDat *dat = elm->dat;
    if(dat)
    {
      tNodeInfo *info = elm->dat->info;
      int f;
      for(f=0; f<6; f++)  /* go over elm-faces */
        if(info->nnbinfo[f] < 0)
        {
          struct list_head *pos1, *sav;
          struct list_head fnb_head;
          /* put the nbs of elm into fnb_head list */
          INIT_LIST_HEAD(&fnb_head);
          amr_make_fnb_list(elm, f, mesh->nmyelm, mesh->myelm, &fnb_head);
          /* add all nbs to amr_elm_nbinfo */
          list_for_each_safe(pos1, sav, &fnb_head)
          {
            tGlist *elem = list_entry(pos1, tGlist, list);
            tElm *nb = elem->entry;
            amr_elm_nbinfo_add_nbeploc(elm, f, 1, nb->eploc);
            /* once nb->eploc is in amr_elm_nbinfo, del elem with nb */
            glist_elem_del(elem);
          }
        }
    }
  }
}


/* Update amr_elm_nbinfo vars on all faces where:
   elm->dat->info->nnbinfo[f] < 0.
   This is done for all elms on all ranks, but it uses the info in nbranks
   and ef to communicate only with nb ranks. */
int amr_update_elm_nbinfo_if_nnbinfo_negative_ef(tMesh *mesh,
                                                 khash_t(u32) *nbranks,
                                                 khash_t(u32_gptr) *ef)
{
  int f, mywork;
  //int rank=nMPI_rank();
  //int size=nMPI_size();
  tCom *com0, *com1, *com2, *scom;
  int rq0, rq, rq2, srq;
  khiter_t ki;

  /* send/recv com to/from other ranks */
  com0 = alloc_com(sizeof(ulong), 1);

  /* send number of elms that I want info about from the other ranks and
     recv number of elms that other ranks need info about */
  forkhiter(nbranks, ki)
  {
    unsigned rk = kh_key(nbranks, ki);
    khiter_t ki2 = kh_get(u32_gptr, ef, rk); /* get kh iter in ef */
    ulong *ns = checked_calloc(1, sizeof(ns[0]));
    ulong *nr = checked_calloc(1, sizeof(nr[0]));

    /* count how many are in ef for rank rk */
    ns[0] = 0;
    if(ki2!=kh_end(ef)) /* if rk is in ef, we count */
    {
      struct list_head *fhead = kh_val(ef, ki2);
      for(f=0; f<6; f++)
        ns[0] += 1 + list_count_nodes(&(fhead[f]));
        /* we add 1 to also send the number of eplocs on each face */
    }
    else
    {
      for(f=0; f<6; f++) ns[0] += 1;
    }

    /* tell how many eplocs (ns[0]) I want to send to rank rk, and also
       find out how many it wants to send to me (nr[0]) */
    rq = append_buffers_to_com(com0, ns,1, nr,1);
    nMPI_Isend_Irecv_com(com0, rq, nMPI_UNSIGNED_LONG, rk, 10,10, WORLD,WORLD);
  }

  /* wait for sends and recvs in com0 */
  nMPI_Waitall_com_send(com0);
  nMPI_Waitall_com_recv(com0);

  /* send/recv com to/from other ranks */
  com1 = alloc_com(sizeof(ulong), 1);

  /* send elms that I want info about from the other ranks and
     recv elms that other ranks need info about */
  rq0=0;
  forkhiter(nbranks, ki)
  {
    unsigned rk = kh_key(nbranks, ki);
    khiter_t ki2 = kh_get(u32_gptr, ef, rk); /* get kh iter in ef */
    ulong i;
    ulong *ns = get_com_send_buf(com0, rq0); //num.of eplocs to send to rk
    ulong *nr = get_com_recv_buf(com0, rq0); //num.of eplocs to recv from rk
    tEploc *sef = checked_calloc(ns[0], sizeof(sef[0]));
    tEploc *ref = checked_calloc(nr[0], sizeof(ref[0]));

PRFs(":1r ");printf("ns[0]=%lu nr[0]=%lu\n", ns[0], nr[0]);

    if(ki2!=kh_end(ef)) /* if rk is in ef */
    {
      struct list_head *fhead = kh_val(ef, ki2);

      /* now fill the sef array, i is index into sef */
      for(i=0, f=0; f<6; f++)
      {
        ulong num;
        struct list_head *pos0;

        /* put in number on face f */
        num = list_count_nodes(&(fhead[f]));
        memcpy(&sef[i], &num, sizeof(num));
        i++;

        /* fill sef on face f */
        list_for_each(pos0, &(fhead[f]))
        {
          tElm *elm = glist_entry(pos0);
          memcpy(&sef[i], elm->eploc, sizeof(sef[0]));
          i++;
        }
      }
      if(i!=ns[0]) errorexit("i!=ns");
    }

    rq = append_buffers_to_com(com1, sef,ns[0], ref,nr[0]);
    nMPI_Isend_Irecv_com(com1, rq, nMPIvars->TEPLOC, rk, 20,20, WORLD,WORLD);

PRFs(":1b ");printf("ns[0]=%lu nr[0]=%lu\n", ns[0], nr[0]);
PRFs(":2a sef=");
for(int ii=0; ii<ns[0]; ii++) printeploc_s(&(sef[ii]), " ");
//PRFs(":2b ");
//for(int ii=0; ii<nr[0]; ii++) printeploc_s(&(ref[ii]), " ");
fflush(stdout);
//abort();
    rq0++;
  }

  /* MPI is still transferring stuff, but by now we have finished the
     Sends and Recvs to exchange all elm eplocs about which we need info */

  /* more coms to exchange nb results with other ranks */
  com2 = alloc_com(sizeof(ulong), 0);
  scom = alloc_com(sizeof(tEploc), 1);

  /* wait for each ref and then find nbs and put them into ef0_nbs array,
     and also send the number of entries (nmyEplocs) in ef0_nbs */
  rq=0;
  forkhiter(nbranks, ki)
  {
    unsigned rk = kh_key(nbranks, ki);
    //khiter_t ki2 = kh_get(u32_gptr, ef, rk); /* get kh iter in ef */
    ulong i;
    //tEploc *sef = get_com_send_buf(com, rq);
    tEploc *ref = get_com_recv_buf(com1, rq);

    /* ef0_nbs is a large array, where we will store all nbs of all the
       nef_f elms rank rk needs nb info about, for all faces f. Layout is:
        ef0_nbs = |nef_0|nnb0|nb_eploc[0...nnb0-1]|
                        |nnb1|nb_eploc[0...nnb1-1]| <--all entries have
                        ...                            sizeof(tEploc) bytes
                  |nef_5|nnb0|nb_eploc[0...nnb0-1]|
                        |nnb1|nb_eploc[0...nnb1-1]|
                        ... */
    tArray *ef0_nbs = alloc_array1d((18*sizeof(tEploc))/sizeof(double));
    ulong ef0_nbs_idx = 0; /* index of next entry to add */
    /* number of tEploc sized entries in ef0_nbs */
    ulong *nsE = checked_calloc(1, sizeof(nsE[0]));
    ulong *nrE = checked_calloc(1, sizeof(nrE[0]));

    /* process eplocs that others want to know about */
    nMPI_Wait_com_recv(com1, rq); /* wait for request number rq */

PRFs(":3 ref=");
for(int ii=0; ii<16; ii++) printeploc_s(&(ref[ii]), " ");
fflush(stdout);

    /* do work on ref array for request rq and find all nbs of all in ref */
    for(i=0, f=0; f<6; f++)
    {
      union { tEploc eploc; ulong num; } eploc2ulong;
      ulong nef, efi;

      /* read nef */
      eploc2ulong.eploc = ref[i];
      nef = eploc2ulong.num;
      i++;

      /* put the number nef into ef0_nbs array */
      memcpy_to_array_redim(ef0_nbs, sizeof(tEploc), ef0_nbs_idx,
                            &nef, sizeof(nef));
      ef0_nbs_idx++;

      for(efi=0; efi<nef; efi++)
      {
        struct list_head *pos1, *sav;
        struct list_head fnb_head;
        ulong nnb, j;
        tElm *elmi = alloc_elm_of_eploc(mesh, &ref[i]);
        i++;

        /* put the nbs of elmi into fnb_head list */
        INIT_LIST_HEAD(&fnb_head);
        nnb = amr_make_fnb_list(elmi, f, mesh->nmyelm, mesh->myelm,
                                &fnb_head);

        /* put the number nnb into ef0_nbs array */
        memcpy_to_array_redim(ef0_nbs, sizeof(tEploc), ef0_nbs_idx,
                              &(nnb), sizeof(nnb));
        ef0_nbs_idx++;

        /* put the nb eplocs into ef0_nbs array */
        j=0;
        list_for_each_safe(pos1, sav, &fnb_head)
        {
          tGlist *elem = list_entry(pos1, tGlist, list);
          tElm *nb = elem->entry;
          /* put nb->eploc into ef0_nbs array */
          memcpy_to_array_redim(ef0_nbs, sizeof(tEploc), ef0_nbs_idx,
                                nb->eploc, sizeof(tEploc));
          ef0_nbs_idx++;
          j++;

          /* once nb->eploc is in ef0_nbs, del elem with nb */
          glist_elem_del(elem);
        }
        if(nnb!=j) errorexit("nnb!=j");

        /* now &fnb_head is freed, so just free the elmi */
        free_elm(elmi);
      }
    } /* end loop over ref */

    /* save number of tEploc sized entries in ef0_nbs */
    nsE[0] = ef0_nbs_idx;

    /* wait for sef array */
    nMPI_Wait_com_send(com1, rq); /* wait for request number rq */

    /* send/recv my nsE/nrE */
    rq2 = append_buffers_to_com(com2, nsE,1, nrE,1);
    nMPI_Isend_Irecv_com(com2, rq2, nMPI_UNSIGNED_LONG, rk, 30,30, WORLD,WORLD);

    /* send ef0_nbs */
    srq = append_buffers_to_com(scom, ef0_nbs->eploc,nsE[0], NULL,0);
    nMPI_Isend_com(scom, srq, nMPIvars->TEPLOC, rk, 40, WORLD);

PRFs(":4 ");printf("nsE[0]=%lu nrE[0]=%lu\n ef0_nbs=", nsE[0], nrE[0]);
printarray_eploc(ef0_nbs, 0);
fflush(stdout);
//abort();

    rq++;
  }

  /* we are now done with all in com0 and com1 */
  free_com(com1);
  free_com(com0);

  /* clear nnbinfo for new entries */
  amr_disable_elm_nbinfo_if_nnbinfo_negative(mesh);

  /* work on my own elms, while MPI is busy:
     find all my elmfaces that have nnbinfo<0 */
  //if(kh_size(nbranks)==0)
  amr_update_elm_nbinfo_locally_if_nnbinfo_negative(mesh);

  /* recv results about my nbs */
  rq=0;
  mywork=0;
  forkhiter(nbranks, ki)
  {
    unsigned rk = kh_key(nbranks, ki);
    khiter_t ki2 = kh_get(u32_gptr, ef, rk); /* get kh iter in ef */
    /*
    //if( (rk>rank) && (mywork==0) ) //this is in the correct order
    if(mywork==0) // this does all mine first
    {
      mywork = 1;
      // work on my own elms, while MPI is busy:
      // find all my elmfaces that have nnbinfo<0
      amr_update_elm_nbinfo_locally_if_nnbinfo_negative(mesh);
    }
    */

    nMPI_Wait_com_recv(com2, rq); /* wait for request number rq */

    {
      ulong *nr = get_com_recv_buf(com2, rq); //num.of eplocs to recv from rk
      tEploc *ef0_nbs = checked_calloc(nr[0], sizeof(ef0_nbs[0]));
      /* ef0_nbs is a large array, where we have stored all nbs of all the
         nef_f elms rank rk needs nb info about, for all faces f. Layout is:
          ef0_nbs = |nef_0|nnb0|nb_eploc[0...nnb0-1]|
                          |nnb1|nb_eploc[0...nnb1-1]| <--all entries have
                          ...                            sizeof(tEploc) bytes
                    |nef_5|nnb0|nb_eploc[0...nnb0-1]|
                          |nnb1|nb_eploc[0...nnb1-1]|
                          ... */
      //ulong ef0_nbs_idx = 0; /* index of next entry to add */
      //ulong nmyEplocs;       /* number of tEploc sized entries in ef0_nbs */

     if(ki2!=kh_end(ef)) /* if rk is in ef */
       if(nr==0) errorexit("what??? how can nr be zero?");

PRFs(":5 ");printf("nr[0]=%lu\n", nr[0]);
fflush(stdout);
//abort();

      /* recv ef0_nbs from rk */
      nMPI_Recv(ef0_nbs,nr[0], nMPIvars->TEPLOC, rk, 40);
      //FIXME: can we pair Isend with Recv???

PRFs(":6 ");
for(int ii=0; ii<nr[0]; ii++) printeploc_s(&(ef0_nbs[ii]), " ");
fflush(stdout);
//abort();

      /* add all in ef0_nbs to my elms as nbs */
      /*******************************************/
      /* read ef0_nbs to build info about my nbs */
      /*******************************************/
      if(ki2!=kh_end(ef)) /* if rk is in ef */
      {
        struct list_head *fhead = kh_val(ef, ki2);
        ulong epi;
        /* each ef0_nbs is a tEploc array, where we will store all nbs of
           all the nef0[f] elms rank rk needs nb info about, for all faces f.
           Layout is:
           ef0_nbs = |nelms[0]|nnb0|nb_eploc[0...nnb0-1]|
                              |nnb1|nb_eploc[0...nnb1-1]|
                              ...
                     |nelms[5]|nnb0|nb_eploc[0...nnb0-1]|
                              |nnb1|nb_eploc[0...nnb1-1]|
                              ... */
        epi = 0;
        for(f=0; f<6; f++)
        {
          union { tEploc e; ulong ul; } e2ul;
          struct list_head *pos1;
          ulong nelms, ei;

          /* pos of 1st elm in ef list */
          pos1 = fhead[f].next;

          /* get number of elms nelms out of ef0_nbs */
          e2ul.e = ef0_nbs[epi++];
          nelms  = e2ul.ul;

          for(ei=0; ei<nelms; ei++)
          {
            tGlist *elem;
            tElm *elm;
            ulong nnb;

            /* get elm at pos1 */
            elem = list_entry(pos1, tGlist, list);
            elm  = elem->entry;
            pos1 = pos1->next; /* go forward now, so we could del below */

            /* get number of nbs out of ef0_nbs */
            e2ul.e = ef0_nbs[epi++];
            nnb    = e2ul.ul;

            //printeploc_s(elm->eploc, " ");
            //printf("f%d r=%d ", f, r);
            //printf("ei=%lu  nnb=%lu epi=%lu ", ei, nnb, epi);
            //if(nnb) printeploc_s(&(ef0_nbs[epi]), " ...");
            //printf("\n");

            /* add all nbs in ef0_nbs to var amr_elm_nbinfo */
            amr_elm_nbinfo_add_nbeploc(elm, f, nnb, &(ef0_nbs[epi]));
            epi += nnb;
            PRFs(":7 elm=");printelm0(elm,"\n");print_amr_elm_nbinfo(elm, f);
          }
        } /* end for f */
      } /* end func that builds nb-info from ef0_nbs */
      free(ef0_nbs);
    }
    rq++;
  }

  /* wait for sends in com2 */
  nMPI_Waitall_com_send(com2);

  /* we are now done with all in com2 */
  free_com(com2);

  /* wait for sends in scom */
  nMPI_Waitall_com_send(scom);
  /* we are now done with all in scom */
  free_com(scom);
fflush(stdout);
//abort();

  /* finally set nnbinfo according to the new nb-info we have now,
     but we keep them negative for now */
  amr_elm_nbinfo_set_nnbinfo_mesh(mesh, 0);

  return 0;
}


/* Add elm to nb->fnb, where nb=elm->fnb[f][ni] */
void amr_add_elm_to_nbelm_fnb(tElm *elm, int f, int ni)
{
  tElm *nb = elm->fnb[f][ni];
  int nb_nfnb;
  tElm **nb_fnb;
  /* figure out face on nb */
  int nb_f = amr_get_nbface(elm, f, nb);
  if(nb_f<0) errorexit("nb_f not found");

  nb_nfnb = nb->nfnb[nb_f];
  nb_fnb = nb->fnb[nb_f];

  /* make room for one more */
  nb_fnb = realloc(nb_fnb, (nb_nfnb+1) * sizeof(nb_fnb[0]));
  nb->fnb[nb_f]  = nb_fnb;
  nb->nfnb[nb_f] = nb_nfnb+1;

  /* add elm to nb->fnb */
  nb_fnb[nb_nfnb] = elm;
}

/* Add elm to nb->fnb, but only if elm is not there yet.
   Here nb=elm->fnb[f][ni] */
void amr_unionadd_elm_to_nbelm_fnb(tElm *elm, int f, int ni)
{
  tElm *nb = elm->fnb[f][ni];
  int nb_nfnb;
  tElm **nb_fnb;
  /* figure out face on nb */
  int nb_f = amr_get_nbface(elm, f, nb);
  if(nb_f<0) errorexit("nb_f not found");

  nb_nfnb = nb->nfnb[nb_f];
  nb_fnb = nb->fnb[nb_f];

  /* see if elm is already in nb_fnb */
  tElm **f_elm = amr_elmarray_linsearch(nb_nfnb, nb_fnb, elm);
  /* if elm is not there yet add it to nb_fnb = nb->fnb[nb_f] */
  if(!f_elm)
  {
    /* make room for one more */
    nb_fnb = realloc(nb_fnb, (nb_nfnb+1) * sizeof(nb_fnb[0]));
    nb->fnb[nb_f]  = nb_fnb;
    nb->nfnb[nb_f] = nb_nfnb+1;

    /* add elm to nb->fnb */
    nb_fnb[nb_nfnb] = elm;
  }
  else
  {
    if(f_elm[0] != elm)
      errorexit("found elm with same loc that is not elm");
  }
}


/* erase all that is in elm->fnb */
void amr_erase_all_elm_fnb(tMesh *mesh)
{
  formyelms(mesh)
  {
    tElm *elm = MyElm;
    int f;
    for(f=0; f<6; f++)
    {
      free(elm->fnb[f]);
      elm->fnb[f]  = NULL;
      elm->nfnb[f] = 0;
    }
  }
}

/* Update elm->fnb[f] from the amr_elm_nbinfo0+f var if elm->nfnb[f] < 0,
   and also add nb elms to mesh->nbelm if this rank does not have them */
int amr_elm_nbinfo_to_elm_fnb(tMesh *mesh)
{
  int rank = nMPI_rank();

  formyelms(mesh)
  {
    tElm *elm = MyElm;
    int f;

    /* go over all 6 faces of elm */
    for(f=0; f<6; f++)
    {
      int i_nbinfo = amr->elm_nbinfo0 + f;
      tArray *nbinfo = VarA(elm, i_nbinfo);
      int nnbinfo_f;
      int i, neplocs;

      /* if there is no nbinfo there are no nbs */
      if(!nbinfo)
      {
        free(elm->fnb[f]);
        elm->nfnb[f] = 0;
        elm->fnb[f]  = NULL;
        continue;
      }

      /* we only set the elm->fnb[f] if an update is needed, i.e.
         If nnbinfo[f]<0 or elm->nfnb[f]!=nnbinfo[f] */
      nnbinfo_f = elm->dat->info->nnbinfo[f];
      if( (nnbinfo_f >= 0) && (elm->nfnb[f] == nnbinfo_f) )
        continue; /* do nothing */

      /* If we get here, we add the nbs in nbinfo to fnb */
      neplocs = array_Neplocs(nbinfo); //num. of nbs we have

      /* first free and then allocate room for neighbors */
      free(elm->fnb[f]);
      elm->nfnb[f] = neplocs;
      elm->fnb[f]  = NULL;
      if(neplocs)
        elm->fnb[f] = checked_calloc(neplocs, sizeof(elm->fnb[f][0]));

      /* now set fnb on face f */
      for(i=0; i<neplocs; i++)
      {
        tEploc *eploc = &(nbinfo->eploc[i]);
        int datrank;
        ulong nbidx;

        amr_elmindex_and_elmrank_of_eid(mesh, eploc->eid, &nbidx, &datrank);

        //printf("QQQQQQ ");printeploc(elm->eploc);printf(" f%d\t",f);printeploc(eploc);
        //printf(" nbidx=%lu datrank=%d\n", nbidx, datrank);

        if(datrank == rank) /* get elm of eploc from mesh->myelm */
        {
          tElm *nb = mesh->myelm[nbidx];
          elm->fnb[f][i] = nb;
        }
        else /* get elm of eploc from mesh->myelm of rank datrank */
        {
          /* something like move_node_to_rank would not work here
             because it needs to be called also by the sending rank */
          tElm **f_elm;
          /* make a new empty elm that is missing some info, like elm->n */
          tElm *nb;
          tElm0 nb0[1] = {0}; /* initialze all to zero */

          /* copy eploc and datrank into nb0, and set bbox */
          nb0->eploc[0] = eploc[0];
          nb0->datrank  = datrank;
          amr_set_elm0_bbox(mesh, nb0);

          /* make new nb-elm */
          nb = alloc_elm_of_elmheader(mesh, nb0);

          /* Is nb in mesh->nbelm already? */
          f_elm = amr_elmarray_bsearch(mesh->nnbelm, mesh->nbelm, nb);
          /* if yes we use the nb from mesh->nbelm */
          if(f_elm)
          {
            free_elm(nb);
            nb = f_elm[0];
          }
          else /* otherwise we add nb to mesh->nbelm */
          {
            amr_elmarray_add_sort(&(mesh->nnbelm), &(mesh->nbelm), nb);
          }
          /* NOTE: nb->n and nb->np need to be set later!!! */

          /* finally also point at this nb */
          elm->fnb[f][i] = nb;

          /* add elm also to nb->fnb[nb_f][nb_i] */
          amr_unionadd_elm_to_nbelm_fnb(elm, f, i);
          /* This is done locally and may miss nbs of the elms in nbelm,
             that are on yet other ranks! */
        }
      }
    } /* end loop over f */
  }

  return 0;
}


/* get the full elmheader for all elms in mesh->nbelm from the other rank */
int amr_get_nbelm_elmheaders(tMesh *mesh)
{
  int size = nMPI_size();
  int rank = nMPI_rank();
  /* numbers of elms we send to or recv from rank r: */
  ulong *ns_elm0 = checked_calloc(size, sizeof(ns_elm0[0]));
  ulong *nr_elm0 = checked_calloc(size, sizeof(nr_elm0[0]));
  // s_elm0[r][i] is elm_i that is sent to rank r
  // r_elm0[r][i] is elm_i that is revcd from rank r
  tElm0 **s_elm0 = checked_calloc(size, sizeof(s_elm0[0]));
  tElm0 **r_elm0 = checked_calloc(size, sizeof(r_elm0[0]));
  tCom *scom, *rcom;
  ulong ei;
  int rk;

  /* mesh->nbelm has all nb-elms about which we have to exchange info */
  for(ei=0; ei<mesh->nnbelm; ei++)
  {
    tElm *nb = mesh->nbelm[ei];
    int nbrank = nb->datrank;
    union { tElm *elm; tElm0 *elm0; } e2e0;
    int f, ni;

    /* add nb to list we want to recv from rank nbrank */
    e2e0.elm = nb;
    amr_elm0array_add_sort(&(nr_elm0[nbrank]),&(r_elm0[nbrank]),
                           e2e0.elm0);
    /* the neighbors of nb are my elms that I have to send to rank nbrank */
    for(f=0; f<6; f++)
      for(ni=0; ni<nb->nfnb[f]; ni++)
      {
        tElm *elm = nb->fnb[f][ni];

        /* add elm to list we want to recv */
        e2e0.elm = elm;
        amr_elm0array_unionadd_sort(&(ns_elm0[nbrank]),&(s_elm0[nbrank]),
                                    e2e0.elm0);
      }
  }

  /*
  PRF;printf(": mesh->nnbelm=%lu\n", mesh->nnbelm);
  for(rk=0; rk<size; rk++)
  {
    printf("rk%d: nr=%lu ns=%lu\n", rk, nr_elm0[rk], ns_elm0[rk]);
    printf("r_elm0 = ");
    for(ei=0; ei<nr_elm0[rk]; ei++)
      printeploc_s(r_elm0[rk][ei].eploc, " ");
    printf("\n");

    printf("s_elm0 = ");
    for(ei=0; ei<ns_elm0[rk]; ei++)
      printeploc_s(s_elm0[rk][ei].eploc, " ");
    printf("\n");
  }
  */

  /* send and recv from rank rk */
  scom = alloc_com(sizeof(s_elm0[0][0]), 0);
  rcom = alloc_com(sizeof(r_elm0[0][0]), 0);
  for(rk=0; rk<size; rk++)
    if(rk != rank)
    {
      int rq;

      /* send s_elm0[rk] buffer  */
      rq = append_buffers_to_com(scom, s_elm0[rk],ns_elm0[rk], NULL,0);
      nMPI_Isend_com(scom, rq, nMPIvars->TELM0, rk, 31, WORLD);

      /* recv in r_elm0[rk] */
      rq = append_buffers_to_com(rcom, NULL,0, r_elm0[rk],nr_elm0[rk]);
      nMPI_Irecv_com(rcom, rq, nMPIvars->TELM0, rk, 31, WORLD);
    }

  /* wait for recvs in rcom */
  nMPI_Waitall_com_recv(rcom);
  free_com(rcom);

  /* write the now revcd r_elm0 contents into mesh->nbelm */
  for(rk=0; rk<size; rk++)
    for(ei=0; ei<nr_elm0[rk]; ei++)
    {
      tElm0 *nb0;
      tElm *nb, *nb_header, **f_nb;
      union { tElm *elm; tElm0 *elm0; } e2e0;
      e2e0.elm0 = nb0 = &(r_elm0[rk][ei]);
      nb_header = e2e0.elm; //is a pointer to tElm that points to a tElm0

      /* find nb corresponding to nb_header in mesh->nbelm */
      f_nb = amr_elmarray_bsearch(mesh->nnbelm, mesh->nbelm, nb_header);
      nb = *f_nb;

      /* write nb_header into nb */
      memcpy(nb, nb0, sizeof(nb0[0]));
    }

  //printf("mesh->nbelm:\n");
  //printnbelms(mesh);

  /* wait for sends in scom, then free all send related stuff */
  nMPI_Waitall_com_send(scom);
  free_com(scom);

  /* free all arrays */
  rows_free(r_elm0, size);
  rows_free(s_elm0, size);
  free(nr_elm0);
  free(ns_elm0);
  return 0;
}


/****************************************************************************/
/* functions to get elm-headers from another ranks */
/****************************************************************************/

/* Get the full elmheader for all elms with eids in eidarr.
   In: neids,eidarr  Out: elm0-array
   This needs to be called on all MPI processes! */
int amr_get_elm0_for_eids(tMesh *mesh, ulong neids, ulong *eidarr,
                          tElm0 *elm0)
{
  int size = nMPI_size();
  int rank = nMPI_rank();

  /* numbers of eids that we send to rank r to get elm0 */
  ulong *ns_deseid = checked_calloc(size, sizeof(ns_deseid[0]));
  // s_deseid[r][i] is eid_i that is sent to rank r
  // s_des_ei[r][i] is index of eid in eidarr
  ulong **s_deseid = checked_calloc(size, sizeof(s_deseid[0]));
  ulong **s_des_ei = checked_calloc(size, sizeof(s_des_ei[0]));

  /* numbers of eids that we recv from rank r so we can send it elm0 */
  ulong *nr_deseid = checked_calloc(size, sizeof(nr_deseid[0]));
  // r_deseid[r][i] is eid_i that is sent to rank r
  ulong **r_deseid;

  /* numbers of elms we recv/send from rank r: */
  ulong *nr_elm0 = checked_calloc(size, sizeof(nr_elm0[0]));
  ulong *ns_elm0 = checked_calloc(size, sizeof(ns_elm0[0]));
  // r_elm0[r][i] is elm_i that is revcd from rank r
  // s_elm0[r][i] is elm_i that is sent to rank r
  tElm0 **r_elm0;
  tElm0 **s_elm0;
  tCom *scom, *rcom;
  ulong ei, k;
  int rk;

  /* sort eidarr into s_deseid */
  for(ei=0; ei<neids; ei++)
  {
    ulong eid = eidarr[ei];
    int datrank = amr_rank_of_eid(mesh, eid);
    int num = ns_deseid[datrank];

    if(datrank!=rank)
    {
      s_deseid[datrank] = checked_realloc(s_deseid[datrank],
                                          (num+1) * sizeof(s_deseid[0][0]));
      s_des_ei[datrank] = checked_realloc(s_des_ei[datrank],
                                          (num+1) * sizeof(s_des_ei[0][0]));
      s_deseid[datrank][num] = eid;
      s_des_ei[datrank][num] = ei;
      ns_deseid[datrank] = num+1;
    }
    else
    {
      union { tElm *elm; tElm0 *elm0; } e2e0;
      ulong elmindex;
      amr_elmindex_and_elmrank_of_eid(mesh, eid,
                                      &elmindex, &datrank);
      /* just fill in elm0 */
      e2e0.elm = mesh->myelm[elmindex];
      elm0[ei] = *(e2e0.elm0);
    }
  }

  /* we now know the size of r_elm0 */
  r_elm0 = rows_calloc(size, ns_deseid, sizeof(r_elm0[0][0]));
  //printf("r_elm0 from ns_deseid:\n");
  //rows_print_sizes(size, ns_deseid, sizeof(r_elm0[0][0]));

  ////nr_deseid[0] = nr_deseid[2] = 66;
  ////printf("fake nr_deseid:\n");
  ////rows_print_sizes(size, nr_deseid, sizeof(r_deseid[0][0]));

  /* send and recv coms */
  scom = alloc_com(sizeof(s_elm0[0][0]), 0);
  rcom = alloc_com(sizeof(r_elm0[0][0]), 0);

  /* find out how much needs to be exchanged */
  for(rk=0; rk<size; rk++)
    if(rk != rank)
    {
      int rq;
      /* tell how many I want from each rank rk */
      {
        rq = append_buffers_to_com(scom, &(ns_deseid[rk]),1, NULL,0);
        nMPI_Isend_com(scom, rq, nMPI_UNSIGNED_LONG, rk, 10, WORLD);
      }
      /* get how many the other ranks want */
      {
        rq = append_buffers_to_com(rcom, NULL,0, &(nr_deseid[rk]),1);
        nMPI_Irecv_com(rcom, rq, nMPI_UNSIGNED_LONG, rk, 10, WORLD);
      }
    }

  /* wait for nr_deseid[rk] */
  nMPI_Waitall_com_recv(rcom);
  realloc_com_reqs(rcom, 0);

  /* we now know the size of r_deseid[rk] and s_elm0 */
  r_deseid = rows_calloc(size, nr_deseid, sizeof(r_deseid[0][0]));
  s_elm0 =   rows_calloc(size, nr_deseid, sizeof(s_elm0[0][0]));

  //printf("r_deseid from nr_deseid:\n");
  //rows_print_sizes(size, nr_deseid, sizeof(r_deseid[0][0]));
  //printf("s_elm0 from nr_deseid:\n");
  //rows_print_sizes(size, nr_deseid, sizeof(s_elm0[0][0]));

  /* now send/recv s_deseid and r_deseid */
  for(rk=0; rk<size; rk++)
    if(rk != rank)
    {
      int rq;

      /* send s_deseid buffer, i.e. the eids for which I want elm0 from rk */
      if(ns_deseid[rk])
      {
        rq = append_buffers_to_com(scom, s_deseid[rk],ns_deseid[rk], NULL,0);
        nMPI_Isend_com(scom, rq, nMPI_UNSIGNED_LONG, rk, 20, WORLD);
      }
      /* recv in r_deseid, i.e. the eids for which rk wants to know elm0 */
      if(nr_deseid[rk])
      {
        rq = append_buffers_to_com(rcom, NULL,0, r_deseid[rk],nr_deseid[rk]);
        nMPI_Irecv_com(rcom, rq, nMPI_UNSIGNED_LONG, rk, 20, WORLD);
      }
    }

  /* wait for r_deseid[rk] */
  nMPI_Waitall_com_recv(rcom);
  realloc_com_reqs(rcom, 0);

  /* use r_deseid[rk] to fill in s_elm0[rk] arrays */
  for(rk=0; rk<size; rk++)
    if(rk != rank)
      for(k=0; k<nr_deseid[rk]; k++)
      {
        union { tElm *elm; tElm0 *elm0; } e2e0;
        ulong elmindex;
        int datrank;

        amr_elmindex_and_elmrank_of_eid(mesh, r_deseid[rk][k],
                                        &elmindex, &datrank);
        /* fill in s_elm0 */
        e2e0.elm = mesh->myelm[elmindex];
        s_elm0[rk][k] = *(e2e0.elm0);
      }

  /* wait for sends in scom, then free all send related stuff */
  nMPI_Waitall_com_send(scom);
  realloc_com_reqs(scom, 0);


  /* finally exchange elm0 */
  for(rk=0; rk<size; rk++)
    if(rk != rank)
    {
      int rq;

      /* send elm0 that the others want */
      if(nr_deseid[rk])
      {
        rq = append_buffers_to_com(scom, s_elm0[rk],nr_deseid[rk], NULL,0);
        nMPI_Isend_com(scom, rq, nMPIvars->TELM0, rk, 30, WORLD);
      }
      /* recv elm0 we want from the others */
      if(ns_deseid[rk])
      {
        rq = append_buffers_to_com(rcom, NULL,0, r_elm0[rk],ns_deseid[rk]);
        nMPI_Irecv_com(rcom, rq, nMPIvars->TELM0, rk, 30, WORLD);
      }
    }

  /* wait for r_elm0[rk] in rcom */
  nMPI_Waitall_com_recv(rcom);
  realloc_com_reqs(rcom, 0);


  /* fill in elm0 array from r_elm0[rk] */
  for(rk=0; rk<size; rk++)
    if(rk != rank)
      for(k=0; k<ns_deseid[rk]; k++)
      {
        ei = s_des_ei[rk][k];
        elm0[ei] = r_elm0[rk][k];
      }

  /* wait for final sends */
  nMPI_Waitall_com_send(scom);
  free_com(rcom);
  free_com(scom);

  /* free all arrays */
  rows_free(s_elm0, size);
  rows_free(r_deseid, size);
  rows_free(r_elm0, size);

  rows_free(s_des_ei, size);
  rows_free(s_deseid, size);

  free(ns_elm0);
  free(nr_elm0);

  free(nr_deseid);
  free(ns_deseid);

  return 0;
}

/* allocate and fill an array with all elm0 on rank rk,
   Returns: pointer to array, Out: nelm0s <- size of array,
   returned array needs to be freed by caller */
tElm0 *amr_alloc_get_elm0array_of_rank(tMesh *mesh, int rk, ulong *nelm0s)
{
  ulong nelm0 = amr_nelms_on_rank(mesh, rk);
  tElm0 *elm0ar = checked_calloc(nelm0, sizeof(elm0ar[0]));
  ulong *eidar  = checked_calloc(nelm0, sizeof(eidar[0]));
  ulong eid0 = amr_1st_eid_on_rank(mesh, rk);
  ulong i;

  /* fill eid array */
  for(i=0; i<nelm0; i++) eidar[i] = eid0 + i;

  /* now put all elm0 for the eid array eidar into elm0ar */
  amr_get_elm0_for_eids(mesh, nelm0, eidar, elm0ar);

  free(eidar);
  *nelm0s = nelm0;
  return elm0ar;
}


/****************************************************************************/
/* functions to partially set nb-info after refine or unrefine */
/****************************************************************************/

/* Invalidate nbinfo for all nbs of elm on elmface.
   Return: 0 if all nbs have dat (are on my rank), 1 if one nb has no dat */
int amr_invalidate_nbinfo_of_nbs(tElm *elm, int elmface,
                                 int Keep_nbs_fnb_to_elm)
{
  int ni;
  int nbs_on_other_rank = 0;

  /* go over nbs of elm,elmface and invalidate their nbinfo */
  for(ni=0; ni<elm->nfnb[elmface]; ni++)
  {
    tElm *nb = elm->fnb[elmface][ni];
    int nb_f, nb_ni;

    if(!nb) continue; /* do nothing if there is no nb */

    /* face of nb */
    nb_f = amr_get_nbface(elm,elmface, nb);
    if(nb_f<0) errorexit("nb_f not found");

    if(nb->dat)
    {
      int nnb = nb->dat->info->nnbinfo[nb_f];
      /* invalidate nbinfo */
      if(nnb>=0) nb->dat->info->nnbinfo[nb_f] = -nnb-1;
    }
    else
    {
      nbs_on_other_rank = 1;
    }

    //CHECK:
    /* Set pointers in nb that point back at elm to NULL.
       This is needed if we remove elm and we then want to remove nb soon
       after and then call:
       amr_invalidate_nbinfo_of_nbs(nb,...); */
    if(!Keep_nbs_fnb_to_elm)
      for(nb_ni=0; nb_ni<nb->nfnb[nb_f]; nb_ni++)
      {
        tElm *nbnb =  nb->fnb[nb_f][nb_ni];
        if(nbnb==elm) nb->fnb[nb_f][nb_ni] = NULL;
      }
  }
  return nbs_on_other_rank;
}

/* Invalidate nbinfo for all nbs of elm on elmface.
   Return: 0 if all nbs have dat (are on my rank), # of nbs witout dat */
int amr_invalidate_nbinfo_of_all_nbs(tElm *elm, int Keep_nbs_fnb)
{
  int f;
  int nbs_on_other_rank = 0;
  for(f=0; f<6; f++)
    nbs_on_other_rank += amr_invalidate_nbinfo_of_nbs(elm, f, Keep_nbs_fnb);
  return nbs_on_other_rank;
}

/* Go over mesh->nbelm list and invalidate nbinfo for all my elms that
   are nbs of any elm in mesh->nbelm. */
void amr_invalidate_nbinfo_of_mesh_nbelm_nbs(tMesh *mesh, int Keep_nbs_fnb)
{
  ulong ei;
  for(ei=0; ei < mesh->nnbelm; ei++)
  {
    tElm *elm = mesh->nbelm[ei];
    amr_invalidate_nbinfo_of_all_nbs(elm, Keep_nbs_fnb);
  }
}

/* Remove mesh->nbelm and make sure all nbinfo about it is deleted */
void amr_remove_mesh_nbelm(tMesh *mesh, int Keep_nbs_fnb)
{
  ulong ei;

  /* first make sure nobody has info about elms in nbelm */
  amr_invalidate_nbinfo_of_mesh_nbelm_nbs(mesh, Keep_nbs_fnb);

  for(ei=0; ei < mesh->nnbelm; ei++)
  {
    tElm *elm = mesh->nbelm[ei];
    free_elm(elm);
  }
  free(mesh->nbelm);
  mesh->nbelm  = NULL;
  mesh->nnbelm = 0;
}


/****************************************************************************/
/* functions to record info about neighboring ranks */
/****************************************************************************/

/* Record all the nbranks I am in contact with, i.e. add datranks of all nbs
   to the hash set called nbranks. */
/* needs both:   khash_t(u32) *nbranks = kh_init(u32);
                 kh_destroy(u32, nbranks);              */
int amr_khset_add_nb_ranks(tMesh *mesh, khash_t(u32) *nbranks)
{
  int nadded = 0;
  ulong ei;
  for(ei=0; ei<mesh->nnbelm; ei++)
  {
    tElm *elm = mesh->nbelm[ei];
    int is_missing;
    kh_put(u32, nbranks, elm->datrank, &is_missing);
    if(is_missing) nadded++;
  }
  return nadded;
}

/* Record elm,face for the key rank in the hash table ef. The value of key
   rank is 6 lists (one for each face) to which we append the elm. */
void amr_khmap_add_elm_face_for_rank(khash_t(u32_gptr) *ef, int rank,
                                     tElm *elm, int face)
{
  int is_missing;
  khiter_t ki;
  struct list_head *fhead;

  ki = kh_put(u32_gptr, ef, rank, &is_missing);
  if(is_missing)
  {
    int f;
    /* alloc space for 6 list heads, one for each face */
    kh_val(ef, ki) = calloc(6, sizeof(struct list_head));
    fhead = kh_val(ef, ki);
    for(f=0; f<6; f++) INIT_LIST_HEAD(&(fhead[f]));
  }

  /* add elm to list for this face */
  fhead = kh_val(ef, ki);
  //PRF;printf("1: rank%d ki=%u f%d %p\n", rank, ki, face, &(fhead[face]));
  //printelmglist(&(fhead[face]));
  glist_entry_add_tail(elm, &(fhead[face]));
  //PRF;printf("2: rank%d ki=%u f%d %p\n", rank, ki, face, &(fhead[face]));
  //printelmglist(&(fhead[face]));
}

/* free mem allocated by amr_khmap_add_elm_face_for_rank allocs
   for the 6 lists */
void amr_khmap_free_all_lists(khash_t(u32_gptr) *ef)
{
  khiter_t ki;
  forkhiter(ef, ki)
  {
    struct list_head *fhead = kh_val(ef, ki);
    int f;
    /* clear the 6 lists in val */
    for(f=0; f<6; f++) glist_free_elems(&(fhead[f]));
    free(fhead);
    kh_val(ef, ki) = NULL;
  }
}

/* Find all ranks that elm touches on face, and save elm,face once for each
   touching rank. (Note: fnbranks is only there to track if elm was already
   added once before.) */
void amr_khmap_add_elm_forface(khash_t(u32_gptr) *ef, tElm *elm, int face)
{
  khash_t(u32) *fnbranks = kh_init(u32); /* empty nb ranks set for face */
  int ni;
  //printelm(elm);
  for(ni=0; ni<elm->nfnb[face]; ni++)
  {
    tElm *nb = elm->fnb[face][ni];
    int nb_rk = nb->datrank;
    int is_missing;

    kh_put(u32, fnbranks, nb_rk, &is_missing); /* record nb rank */
    /* if this is the 1st time we find this rank on this face, add elm */
    if(is_missing)
      amr_khmap_add_elm_face_for_rank(ef, nb_rk, elm, face);
  }
  kh_destroy(u32, fnbranks);
}

/* Find all ranks that elm touches, and save elm,face once for each
   touching rank. (Note: fnbranks is only there to track if elm was already
   added once before.) */
void amr_khmap_add_negelm_forallfaces(khash_t(u32_gptr) *ef, tElm *elm)
{
  int f;
  for(f=0; f<6; f++)
    if(elm->dat->info->nnbinfo[f] < 0)
      amr_khmap_add_elm_forface(ef, elm, f);
}

/* Find all ranks that parent elm touches, and save child0-7,face once for
   each touching rank. (Note: fnbranks is only there to track if children
   were already added once before.) */
void amr_khmap_add_negchildren_forallparentfaces(khash_t(u32_gptr) *ef,
                                                 tElm *child0, tElm *parent)
{
  khash_t(u32) *fnbranks = kh_init(u32);
  int p_rk = parent->datrank;

  PRF;printelm(parent);
  //printelm(child0);

  int f, ni;
  for(f=0; f<6; f++)
  {
    struct list_head *pos_ijk;
    int ijk;
    int nnbinfo_neg;

    /* check if any child has nnbinfo<0 */
    nnbinfo_neg = 0;
    pos_ijk = &child0->list; /* pos of 1st child */
    for(ijk=0; ijk<8; ijk++)
    {
      tElm *child = list_entry(pos_ijk, tElm, list);
      if(child->dat->info->nnbinfo[f] < 0) { nnbinfo_neg = 1; break; }
      pos_ijk = pos_ijk->next;  /* pos of next child */
    }

    /* if yes, try to add them */
    if(nnbinfo_neg)
    {
      kh_clear(u32, fnbranks); /* empty nb ranks set for face f */

      for(ni=0; ni<parent->nfnb[f]; ni++)
      {
        tElm *nb = parent->fnb[f][ni];
        int nb_rk = nb->datrank;

        if(nb_rk != p_rk)
        {
          int is_missing;
          kh_put(u32, fnbranks, nb_rk, &is_missing); /* record nb rank */
          /* if this is the 1st time we find this rank on this face,
             add the 4 children on this face */
          if(is_missing)
          {
            pos_ijk = &child0->list; /* pos of 1st child */
            for(ijk=0; ijk<8; ijk++)
            {
              if(connections_ijk_is_at_parentface(ijk, f))
              {
                tElm *child = list_entry(pos_ijk, tElm, list);
                amr_khmap_add_elm_face_for_rank(ef, nb_rk, child, f);
              }
              pos_ijk = pos_ijk->next;  /* pos of next child */
            }
          }
        }
      }
    } /* end if(nnbinfo_neg) */
  }
  kh_destroy(u32, fnbranks);
}

/* Find all ranks that children touch, and save parent,face once for
   each touching rank. (Note: fnbranks is only there to track if parent
   was already added once before.) */
void amr_khmap_add_negparent_forallchildrenfaces(khash_t(u32_gptr) *ef,
                                                 tElm *parent,
                                                 struct list_head *ch_head)
{
  khash_t(u32) *fnbranks = kh_init(u32);
  int p_rk = parent->datrank;
  int f, ni;
  for(f=0; f<6; f++)
  {
    if(parent->dat->info->nnbinfo[f] < 0)
    {
      struct list_head *pos_ijk;

      kh_clear(u32, fnbranks); /* empty nb ranks set for face f */
      list_for_each(pos_ijk, ch_head)
      {
        tElm *ch_ijk = list_entry(pos_ijk, tElm, list);

        for(ni=0; ni<ch_ijk->nfnb[f]; ni++)
        {
          tElm *nb = ch_ijk->fnb[f][ni];
          int nb_rk = nb->datrank;

          if(nb_rk != p_rk)
          {
            int is_missing;
            kh_put(u32, fnbranks, nb_rk, &is_missing); /* record nb rank */
            /* if this is the 1st time we find this rank on this face,
               add parent */
            if(is_missing)
              amr_khmap_add_elm_face_for_rank(ef, nb_rk, parent, f);
          }
        }
      }
    }
  }
  kh_destroy(u32, fnbranks);
}

/* Go over mesh->nbelm list, invalidate nbinfo for all my elms that
   are nbs of any elm in mesh->nbelm, and record them in ef */
void amr_invalidate_nbinfo_of_mesh_nbelm_nbs_ef(tMesh *mesh,
                                                khash_t(u32_gptr) *ef)
{
  khash_t(u64) *rk_elm_f = kh_init(u64);
  ulong key_n[] = {nMPI_size(), mesh->nnbelm, 6};
  ulong ei;

  /* make nnbinf0<0 but keep all fnb pointers */
  amr_invalidate_nbinfo_of_mesh_nbelm_nbs(mesh, 1);

  /* go over nbelms */
  for(ei=0; ei < mesh->nnbelm; ei++)
  {
    tElm *elm = mesh->nbelm[ei];
    int elm_rk = elm->datrank;
    int f, ni;

    for(f=0; f<6; f++)
      for(ni=0; ni<elm->nfnb[f]; ni++)
      {
        tElm *nb = elm->fnb[f][ni];
        int nb_f, nb_ni;

        if(!nb) continue; /* do nothing if there is no nb */

        /* face of nb */
        nb_f = amr_get_nbface(elm,f, nb);
        if(nb_f<0) errorexit("nb_f not found");

        /* go over nbs of nb on its face nb_f, and look for elm */
        for(nb_ni=0; nb_ni<nb->nfnb[nb_f]; nb_ni++)
        {
          tElm *nbnb = nb->fnb[nb_f][nb_ni];
          if(nbnb==elm) /* if it points back at elm we need to record it */
          {
            int is_missing;
            ulong nb_lid = calc_local_elm_id(nb);
            ulong key = Ind_n(elm_rk, nb_lid, nb_f,  key_n);
            kh_put(u64, rk_elm_f, key, &is_missing); /* record */
            /* if this is the 1st time we find this rank on this face, add elm */
            if(is_missing)
              amr_khmap_add_elm_face_for_rank(ef, elm_rk, nb, nb_f);
          }
        }
      } /* end ni-loop */
  }
  kh_destroy(u64, rk_elm_f);
}

/* Remove mesh->nbelm, make sure all nbinfo about it is deleted, and
   record what is missing in ef */
void amr_remove_mesh_nbelm_ef(tMesh *mesh, int Keep_nbs_fnb,
                              khash_t(u32_gptr) *ef)
{
  /* record in ef, but keep all pointers */
  amr_invalidate_nbinfo_of_mesh_nbelm_nbs_ef(mesh, ef);

  /* make sure nobody has info about elms in nbelm */
  amr_remove_mesh_nbelm(mesh, Keep_nbs_fnb);
}



/****************************************************************************/

/* We probably do NOT NEED all this stuff BELOW: */

/****************************************************************************/

//FIXME: maybe remove
/* struct that has has both eloc and a face */
typedef struct tELOCFACE {
  tEloc eloc[1];
  int face;
} tElocFace;

void printelocface(const tElocFace *ef)
{
  printeloc(ef->eloc);
  printf(" f=%d\n", ef->face);
}

/****************************************************************************/
/* functions to exchange info between rank_i and rank_{i \pm 1} */
/****************************************************************************/

/* for MPI exchange between neighboring ranks */
typedef struct tELMFL {
  tElm elm_fl[2]; /* first and last elm on rank */
  int nelms;      /* number of elms on rank */
} tElmfl;
typedef struct tNBR {
  tElmfl fl_m1[1]; /* first and last elm on rank-1 */
  tElmfl fl_p1[1]; /* first and last elm on rank+1 */
} tNbr;



/* fill in tElmfl myfl[1] with my first and last elm */
void get_my_Elmfl(tMesh *mesh, tElmfl myfl[1])
{
  long nelms = mesh->nmyelm;

  myfl->nelms = nelms;
  if(nelms <= 0) return;
  myfl->elm_fl[0] = *(mesh->myelm[0]);        /* shallow copies */
  myfl->elm_fl[1] = *(mesh->myelm[nelms-1]);
}

/* exchange first and last elms with rank+1 and rank-1 */
void get_nbr_rank_info(tMesh *mesh)
{
  int size = nMPI_size();
  int rank = nMPI_rank();
  tCom *com;
  int rq;
  tElmfl myfl[1];
  //this is correct:
  //tElmfl *fl_m1 = mesh->nbr->fl_m1;
  //tElmfl *fl_p1 = mesh->nbr->fl_p1;
  //but for now:
  tElmfl *fl_m1 = NULL;
  tElmfl *fl_p1 = NULL;

  /* get my first and last elm from mesh->myelm */
  get_my_Elmfl(mesh, myfl);

  /* for MPI data transfers */
  //FIXME: why sizeof(double)???, should 2nd arg be 0 or 1????
  //com = alloc_com(sizeof(double), 0);
  com = alloc_com(sizeof(char), 0);

  //alloc_com is stupid!!! Its 1st arg should always be sizeof(void *)

  /* send myfl to rank-1 and also receive fl_m1 from rank-1 */
  if(rank>0)
  {
    rq = append_buffers_to_com(com, myfl,sizeof(myfl[0]),
                                    fl_m1,sizeof(fl_m1[0]));
    nMPI_Isend_Irecv_com(com, rq, nMPI_CHAR, rank-1, -1,+1, WORLD, WORLD);
  }

  /* send myfl to rank+1 and also receiv fl_p1 from rank+1 */
  if(rank < size-1)
  {
    rq = append_buffers_to_com(com, myfl,sizeof(myfl[0]),
                                    fl_p1,sizeof(fl_p1[0]));
    nMPI_Isend_Irecv_com(com, rq, nMPI_CHAR, rank+1, +1,-1, WORLD, WORLD);
  }

  /* wait until all sent and received */
  nMPI_Waitall_com(com);
  free_com(com);
}



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
