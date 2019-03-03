/* bfaces.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "amr.h"





/*************************************************************************/
/* funcs to add and remove bfaces */
/*************************************************************************/


/* get mem for 1 bface */
tBface *alloc_bface(tPat *pat, int f)
{
  tBface *bface = calloc(1, sizeof(bface[0]));
  if(!bface) errorexit("not enough memory");
  bface->pat = pat;
  bface->f   = f;
  return bface;
}

/* add a bface to a pat, f denotes the pat face (0 to 5),
   return the new bface */
tBface *add_empty_bface(tPat *pat, int f)
{
  tBface *bface = alloc_bface(pat, f);
  tBface *bface0 = pat->bface0;
  tBface *bf;

  /* add bface to end of list in pat */
  if(bface0)
  {
    for(bf=bface0; bf->next; bf=bf->next) ;
    bf->next = bface;
    bface->prev = bf;
  }
  else
    bface0 = bface;

  /* set some bface info */
  bface->oXi  = -1; /* var indices of other coords not known yet */
  bface->oYi  = -1; /* var indices of other coords not known yet */
  bface->oZi  = -1; /* var indices of other coords not known yet */
  return bface;
}

/* remove a bface with index fi from a pat, return number of bfaces removed */
void remove_bface(tBface *bface)
{
  tPat *pat = bface->pat;
  tBface *obface = bface->obface;
  tBface *bn, *bp;

  if(!bface) return;

  /* remove bface from list */
  bp = bface->prev;
  bn = bface->next;
  if(bn) bn->prev = bp;
  if(bp) bp->next = bn;
  else   pat->bface0 = bn;

  /* now free bface */
  free(bface);

  /* remove pointer to this bface in other bface */
  obface->obface = NULL;
  remove_bface(obface);
}

/* free all bfaces in pat */
void remove_all_bfaces(tPat *pat)
{
  tBface *bface0 = pat->bface0;
  tBface *bf, *bft;

  for(bf=bface0; bf;)
  {
    bft = bf;
    bf  = bf->next;
    remove_bface(bft);
  }
}
