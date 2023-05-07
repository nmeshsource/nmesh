/* connections.c */
/* Wolfgang Tichy, 3/2023 */

#include "nmesh.h"
#include "amr.h"


extern tnMPIvars nMPIvars[1];


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
  const tEploc *pelc = elm_arr[0]->eploc;
  tEloc elc[1];
  int cmp;

  /* we could optimize this by caching un unpacked loc in each elm: */
  eloc_from_eploc(elc, pelc);

  cmp = loccmp(lc, elc);
  PRFs(": ");printeloc_s(lc, " ");printeloc_s(elc, " ");
  printf("--> cmp=%d\n", cmp);
  //printelm(elm_arr[0]);
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


/* Out: return value: number of faces l,loc is on
        patface[f] = 1 if l,loc is on patch face f */
int connections_loc_on_patchface(int l, const char loc[NLOCS],
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
int connections_get_nbloc_InsidePat(int l, const char loc[NLOCS], int face,
                                    char nbloc[NLOCS], int *nb_f)
{
  int patface[6];
  int ijk, nb_ijk, lret;

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
    //char pnbloc[NLOCS]; /* location of parent nb */
    PRF;printf(" at parentface\n");
    //printf(" pnbloc=%s\n", pnbloc);

    /* l-1,loc is parent, write parent nb loc into nbloc */
    lret = connections_get_nbloc_InsidePat(l-1,loc, face, nbloc, nb_f);
    //pnbloc[l-1] = 0; /* add string-end marker */
    //printf(" pnbloc=%s\n", pnbloc);
    //strncpy(nbloc, pnbloc, NLOCS);
    printf("  nbloc=%s\n", nbloc);
    //nbloc[l-1] = nb_ijk + '0';
    //printf(" nbloc=%s\n", nbloc);
    return lret;
  }
  else
  {
    strncpy(nbloc, loc, NLOCS);
    nbloc[l-1] = nb_ijk + '0';
    if(l<NLOCS) nbloc[l] = 0; /* add string-end marker */
    return l;
  }
  return l;
}



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

  for(i=0; i<NPBYTES-1; i++) eploc->ploc[i] = 255;

  for(i=0; i<NLOCS-1; i++) eloc->loc[i] = 7+'0';
  eloc->loc[NLOCS-1]= 5+'0';
  eloc->p=9;
  eloc->l=19;
  eloc->eid=98;
  strcpy(eloc->loc, "1234567564321012345");
  eloc->loc[44] = '1';
  printeloc_s(eloc, "\n");

  eloc_to_eploc(eloc, eploc);
  eloc_from_eploc(eloc2, eploc);

  eloc_to_eploc(eloc2, eploc);
  eloc_from_eploc(eloc3, eploc);

  printeloc_s(eloc2, "\n");
  printeloc_s(eloc3, "\n");
}

/* put data (of n bytes) into a tArray at position of eploc i */
void memcpy_to_array_at_bytestridepos(tArray *ar, ulong i, size_t bytestride,
                                      const void *src, size_t n)
{
  ulong sz = bytestride;
  ulong sd = sizeof(ar->d[0]); // sizeof double
  ulong nd = (sz+sd-1)/sd;     // num. of doubles in bytestride bytes
  ulong ff = nd*sd;            // num. of bytes equiv to nd doubles
  ulong nd_in_n = (n+sd-1)/sd; // num. of doubles in n
  ulong len = nd*i + nd_in_n;

  redim_array(ar, len,1,1);

  memcpy(ar->d + ff*i, src, n);
}

/* get data (of n bytes) from a tArray at position of eploc i */
void memcpy_from_array_at_bytestridepos(tArray *ar, ulong i,
                                        size_t bytestride,
                                        void *dest, size_t n)
{
  ulong sz = bytestride;
  ulong sd = sizeof(ar->d[0]); // sizeof double
  ulong nd = (sz+sd-1)/sd;     // num. of doubles in bytestride bytes
  ulong ff = nd*sd;            // num. of bytes equiv to nd doubles
  memcpy(ar->d + ff*i, dest, n);
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
  char *e = (char *) elm;
  tElm0 *elm0 = (tElm0 *) e;
  amr_set_elm0_bbox(mesh, elm0);
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

/* return the rank that an elm with eploc is on */
int amr_rank_of_elm_eploc(tMesh* mesh, tEploc *eploc)
{
  int size = nMPI_size();
  ulong *eidlim = mesh->eidlim;
  ulong eid = eploc->eid;
  ulong *li;
  size_t off, num;

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

  return off+1;
}

/* init elm0 data as far as possible from eploc */
void amr_init_elm0_from_eploc(tMesh* mesh, tEploc *eploc, tElm0 *elm0)
{
  /* set eploc in elm0 */
  elm0->eploc[0] = eploc[0];

  /* set bbox and datrank, all else in elm0 remains unchanged */
  amr_set_elm0_bbox(mesh, elm0);
  elm0->datrank = amr_rank_of_elm_eploc(mesh, eploc);
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

/* check if elm and nb has common points on faces f and nb_f
   since res in find_elmfacepoints_in_nbface is low, try it both ways */
int elm_common_facepoints(tElm *elm, int f, tElm *nb, int nb_f)
{
  int f1, f2;

  f1 = find_elmfacepoints_in_nbface(elm,f, nb,nb_f);
  if(f1) return 1;

  f2 = find_elmfacepoints_in_nbface(nb,nb_f, elm,f);

  return f2 || f1;
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



/* Look in elm-array arr (in [arr+off,arr+num-1]) to find the elm
   with loc s_eloc and face s_f.
   *s_eloc is a loc where we start searching
   *But we start searching 1st for s_eloc's ancestor on level l0.
   *Return list with elms on face s_f */
int amr_elms_on_eloc_face(long narr, tElm **arr,
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
  printeloc(s_eloc);
  printf(" s_f=%d\n  ", s_f);
  printelmarray(narr, arr);

  for(l = l0; l <= lmax; l++)
  {
Yo(l);
    /* search for ancestor of s_eloc of level l */
    f_eloc[0] = s_eloc[0];
    f_eloc->l = l;
    f_elm = binarysearch(f_eloc, arr, &off, &num, sizeof(*arr), lecmp, NULL);

//printf("off=%zu num=%zu  f_elm pos=%zu\n",
//off, num, (size_t) ((const tElm **)f_elm - arr));
//printf("got ");printelm(*f_elm);

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
  for(ijk = 0; ijk<8; ijk++) /* loop over all children */
  {
    if(connections_ijk_is_at_parentface(ijk, s_f)) /* only 4 are relevant */
    {
      int ret;
      /* child ijk */
      cheloc->loc[l-1] = ijk + '0';
      if(l<NLOCS) cheloc->loc[l] = 0;
      ret = amr_elms_on_eloc_face(narr, arr, off, num, cheloc, s_f, l,
                                  f_elms_head);
      if(ret >=0) lret = l; /* record success */
    }
  }

  /* finally signal failure or success with at least one nb child */
  printf("final lret=%d\n", lret);
  return lret;
}




/* find a neighbors on patchface of elm,elmface in elmarray narr,arr
   In:  elm,elmface, narr,arr
   Out: will append to list nbelocface_head */
int amr_set_patchface_fnb_list(tElm *elm, int elmface,
                               long narr, tElm **arr,
                               struct list_head *fnb_head)
{
  int patface[6]; //, nfaces;
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
  PRFs(": ");printeloc(eloc);printf(" elmface=%d\n", elmface);

  /* sanity check */
  connections_loc_on_patchface(eloc->l, eloc->loc, patface);
  if(!patface[elmface])
    errorexit("call this for patch faces only!!!");

printbfaces_on_f(pat, elmface);

  /* loop over all bfaces on face and find nb */
  forbfacesonface(pat, elmface, bface)
  {
    tBface *obface = bface->obface;

    printf("bface\n");
    printbface(bface);
    printf("obface\n");
    printbface(obface);

    /* do nothing if no other patch face */
    if(!obface) continue;

    /* eloc and face of root elm in other patch */
    nbeloc->l      = 0;
    nbeloc->loc[0] = 0;
    nbeloc->p      = obface->pat->p;
    nb_f = obface->f;

    /* set pos to last entry to start of list_for_each_safe_continue below */
    pos = fnb_head->prev;

    printeloc(nbeloc);printf(" nb_f=%d\n", nb_f);
    /* add all elms in arr on face nb_f to list fnb_head */
    l2=amr_elms_on_eloc_face(narr, arr, 0, narr, nbeloc, nb_f, 0, fnb_head);

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

  PRFs(": fnb list\n");
  list_for_each(pos, fnb_head)
  {
    tElm *fnb = glist_entry(pos);
    tEploc *fnbeploc = fnb->eploc;
    tEloc fnbeloc[1];
    eloc_from_eploc(fnbeloc, fnbeploc);
    printeloc_s(fnbeloc, " ");
  }
  printf("\n");

  return l;
}




/* Look in elm-array arr (in [arr+off,arr+num-1]) to find the nb of
   elm on face elmface.
   In: elm,elmface, narr,arr  =>  Out: fnb_head, Returns list length */
int amr_set_fnb_list(tElm *elm, int elmface, long narr, tElm **arr,
                     struct list_head *fnb_head)
{
  int patface[6];
  const tEploc *eploc = elm->eploc;
  tEloc eloc[1];

  eloc_from_eploc(eloc, eploc);

  PRFs(": ");printeloc(eloc);printf(" elmface=%d\n", elmface);

  connections_loc_on_patchface(eloc->l, eloc->loc, patface);
  /* simple case where elmface is inside the patch */
  if(!patface[elmface])
  {
    tEloc nbeloc[1];
    int nbface, l0;
    nbeloc->p = eloc->p;
for(l0=0; l0<10; l0++) nbeloc->loc[l0] = 'X';
printf("nbeloc->l=%d\n", nbeloc->l);
    l0 = connections_get_nbloc_InsidePat(eloc->l, eloc->loc, elmface,
                                         nbeloc->loc, &nbface);
    nbeloc->l = l0;
    printf(" -> l0=%d nbeloc=", l0);printeloc(nbeloc);printf(" nbface=%d\n", nbface);
    printf("nbeloc->l=%d\n", nbeloc->l);
    //FIXME: we should pass in a better l0 (not just 0), i.e. the one coming
    //       from connections_get_nbloc_InsidePat
    //l0=0;
    amr_elms_on_eloc_face(narr, arr, 0, narr, nbeloc, nbface, l0, fnb_head);
  }
  else /* complicated case where elmface is on patch surface */
  {
    amr_set_patchface_fnb_list(elm,elmface, narr,arr, fnb_head);
  }

  return list_count_nodes(fnb_head);
}

/* set all fnb for all elms on all ranks */
int amr_set_all_fnbs(tMesh *mesh)
{
  struct list_head *pos;
  struct list_head ef0_head[6]; // one list for each face
  ulong nmyef0[6];              // total number of elms
  int f, rk;
  int rank=nMPI_rank();
  int size=nMPI_size();

  /* init */
  for(f=0; f<6; f++) { INIT_LIST_HEAD(&ef0_head[f]); nmyef0[f]=0; }

  /* find all my elmfaces that have fnb=NULL (all ranks do this) */
  list_for_each(pos, &mesh->myelm_head)
  {
    tElm *elm = list_entry(pos, tElm, list);

    /* go over elm-faces */
    for(f=0; f<6; f++)
    {
      /* if fnb=NULL nb info is not there yet */
      if(!elm->fnb[f])
      {
        glist_entry_add_tail(elm, &ef0_head[f]);
        nmyef0[f] += 1;
      }
    }
  }
  /* now we have 6 lists ef0_head[f] that contain elms where the fnb info
     needs to be updated. */

  /* send my lists to the other ranks */
  for(rk=0; rk<size; rk++)
  {
    int nef0[6];

    /* rank rk copies his nmyef0 into nef0 */
    if(rank == rk)
      for(f=0; f<6; f++) nef0[f] = nmyef0[f];

    /* rank rk sends his nef0[f] to all others, to tell how many elms he
       wants to find face nbs of */
    nMPI_Bcast(&nef0[0],6, nMPI_UNSIGNED_LONG, rk);

    /* make list for each face and send them... */
    for(f=0; f<6; f++)
    {
      /* large array where we will store all nbs of all the nef0[f] elms */
      tArray *ef0_nbs = alloc_array1d(sizeof(tEploc));
      ulong ef0_nbs_idx = 0;

      /* put a zero at the start  */
      memcpy_to_array_at_bytestridepos(ef0_nbs, ef0_nbs_idx, sizeof(tEploc),
                                       &(ef0_nbs_idx), sizeof(ulong));

      /* there is something to do only if nef0[f]>0 */
      if(nef0[f])
      {
        // get &ef0_head[f] from rk into elmhead array ef0
        tElm0 *ef0 = checked_calloc(nef0[f], sizeof(ef0[0]));
        //hey maybe we should just send tEploc not tElm0!!!???
        struct list_head *pos0;
        ulong i;

        /* rank rk now fills the ef0 array */
        if(rank == rk)
        {
          i=0;
          list_for_each(pos0, &ef0_head[f])
          {
            tElm *elm = glist_entry(pos0);
            memcpy(&ef0[i], elm, sizeof(ef0[0]));
            i++;
          }
          if(nef0[f]!=i) errorexit("nef0[f]!=i");
        }

        /* broadcast all elmheaders in ef0 from rank rk to all */
        nMPI_Bcast(ef0, nef0[f], nMPIvars->TELM0, rk);
                               //^^^^^^^^^^^^^^^-is this right???

        /* init ef0_nbs and put 1st value into array */
        ef0_nbs_idx = 0;
        memcpy_to_array_at_bytestridepos(ef0_nbs, ef0_nbs_idx, sizeof(tEploc),
                                         &(nef0[f]), sizeof(nef0[f]));
        ef0_nbs_idx++;

        /* all ranks do work on ef0 array and find all nbs of all in ef0 */
        for(i=0; i<nef0[f]; i++)
        {
          tElm *elmi = alloc_elm_of_elmheader(mesh, &ef0[i]);
          struct list_head *pos1;
          struct list_head fnb_head;
          int j, nnb;
          tEploc *nb_eploc;

          INIT_LIST_HEAD(&fnb_head);

          /* put nbs of elmi into fnb_head list */
          nnb = amr_set_fnb_list(elmi, f, mesh->nmyelm, mesh->myelm,
                                 &fnb_head);

          /* memory for eploc of each nb in nb_head list */
          nb_eploc = checked_calloc(nnb, sizeof(nb_eploc[0]));

          /* get nb eploc into nb_eploc array */
          j=0;
          list_for_each(pos1, &fnb_head)
          {
            tElm *nb = glist_entry(pos);
            nb_eploc[j] = nb->eploc[0];
            j++;
          }
          if(nnb!=j) errorexit("nnb!=j");

          /* fill the ef0_nbs array, layout is:
          ef0_nbs = |nef0[f]|nnb0|nb_eploc[0...nnb0-1]|
                            |nnb1|nb_eploc[0...nnb1-1]|... */
          memcpy_to_array_at_bytestridepos(ef0_nbs, ef0_nbs_idx,
                                           sizeof(tEploc),
                                           &(nnb), sizeof(nnb));
          ef0_nbs_idx++;
          memcpy_to_array_at_bytestridepos(ef0_nbs, ef0_nbs_idx,
                                           sizeof(tEploc),
                                           &(nb_eploc[0]),
                                           sizeof(nb_eploc[0])*nnb);
          ef0_nbs_idx += nnb;



          /* now from the nbs of this elmi  and also elmi */
          free(nb_eploc);
          free_elm(elmi);
        }

        //...  ???

        //we now need to send the ef0_nbs arrays of each rank to rank rk
        //do NOT use: nMPI_Bcast(ef0_nbs->d, len???, nMPI_DOUBLE, rk);
        // all need to send, and only rk receives

        // ???

        free(ef0);
      } /* end if(nef0[f]) */

      free_array(ef0_nbs);
    } /* end loop over face f */

  } /* end loop over rk */











  char *buf;
  int sz1 = sizeof(buf[0]);
  nMPI_Win win;

  /* make a RMA window through which we communicate all */
  nMPI_Win_allocate(1000*sz1, sz1, nMPI_INFO_NULL, WORLD, &buf, &win);

  /* send elm0 to all MPI ranks */

  nMPI_Win_free(&win);
  return 0;
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
