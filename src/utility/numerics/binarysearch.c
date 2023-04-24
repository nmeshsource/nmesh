#include <stdio.h>

/* Perform a binary search, in way very similar to bsearch_r.

   The function below is a bit sneaky. After a comparison fails, we
   divide the work in half by moving either left or right. If lim
   is odd, moving left simply involves halving lim: e.g., when lim
   is 5 we look at item 2, so we change lim to 2 so that we will
   look at items 0 & 1.  If lim is even, the same applies.  If lim
   is odd, moving right again involes halving lim, this time moving
   base0offset up one item past pos: e.g., when lim is 5 we change base
   to item 3 and make lim 2 so that we will look at items 3 and 4.
   If lim is even, however, we have to shrink it by one before
   halving: e.g., when lim is 4, we still looked at item 2, so we
   have to make lim 3, then halve, obtaining 1, so that we will only
   look at item 3. */

/* Look for key in base0[i] for i in [*base0offset, *base0offset + *num-1].
   +Both *base0offset and *num get modified during the search, so that after
    return they contain the last interval that we searched in. This can be
    useful if there are many base0[i] that compare equal to the key, because
    the last interval should contain them all, even though it may be too wide.
   +Returns pointer to the base0[i] that has the key if the key was found,
    otherwise it returns NULL. */
void *binarysearch(const void *key, const void *base0,
                   size_t *base0offset, size_t *num, size_t size,
                   int (*compar)(const void *, const void *, void *),
                   void *arg)
{
    int cmp;
    const void *p;
    size_t lim, pos;
    size_t off = *base0offset;

    for(lim = *num; lim != 0; lim >>= 1)
    {
        /* save last reasonble search pars */
        *base0offset = off;
        *num = lim;
        //printf("*base0offset=%zu *num=%zu off=%zu lim=%zu\n",
        //       *base0offset, *num, off, lim);

        pos = off + (lim >> 1);
        p = (const char *) base0 + pos* size;
        cmp = (*compar)(key, p, arg);
        if(cmp == 0)
        {
            return (void *)p;
        }
        if(cmp > 0)  /* key > p: move right, by increasing off */
        {
            off = pos + 1;
            lim--;
        }
        /* else: move left, by halving lim (done in for loop) */
    }
    return NULL;
}

/* check if key also matches left and right of result of binarysearch
   Return vals:
   1 left of result also matches
   2 right of result also matches
   3 left and right of result also match
   0 neither left nor right match */
int binarysearchmore(const void *key, const void *base0,
                     size_t nmemb, size_t size, const void *result,
                     int (*compar)(const void *, const void *, void *),
                     void *arg)
{
    int cmp, more;
    const char *res = (const char *) result;
    size_t pos;

    if(!result) return 0;

    pos = (res - (const char *) base0)/size;
    //printf("result=%p base0=%p pos=%ld\n", result, base0, pos);
    more = 0;

    /* check on left */
    if(pos>0)
    {
        cmp = (*compar)(key, res-size, arg);
        if(cmp == 0) more |= 1;
    }

    /* check on right */
    if(pos<nmemb-1)
    {
        cmp = (*compar)(key, res+size, arg);
        if(cmp == 0) more |= 2;
    }

    return more;
}


/* Look for key in base0[i] for i in [*base0offset, *base0offset + *num-1].
   +If key is found exactly:
     sets *base0offset to the i where key is, sets *num=1
     Returns pointer to the base0[i] that has the key
   +If key is between 2 values:
     sets *base0offset to the i left of key , sets *num=2
     Returns pointer to the base0[i]
   +If key is not bracketed at all:
     Returns NULL */
void *bisectionsearch(const void *key, const void *base0,
                      size_t *base0offset, size_t *num, size_t size,
                      int (*compar)(const void *, const void *, void *),
                      void *arg)
{
    int cmp, cmp_a, cmp_b;
    const void *p;
    size_t pos, pos_a, pos_b;
    size_t off = *base0offset;
    size_t nf;

    pos_a = off;
    pos_b = off + (*num)-1;

    nf = 1;

    pos = pos_a;
    p = (const char *) base0 + pos * size;
    cmp_a = (*compar)(key, p, arg);
    if(cmp_a == 0) goto FoundKey;

    pos = pos_b;
    p = (const char *) base0 + pos * size;
    cmp_b = (*compar)(key, p, arg);
    if(cmp_b == 0) goto FoundKey;

    /* no bracket */
    if(cmp_a*cmp_b > 0)
        return NULL;

    /* tighten bracket */
    while((pos = (pos_a + pos_b)/2) > pos_a)
    {
        p = (const char *) base0 + pos * size;
        cmp = (*compar)(key, p, arg);
        if(cmp == 0) goto FoundKey;
        if(cmp*cmp_b > 0) { pos_b = pos;  cmp_b = cmp; }
        else              { pos_a = pos;  cmp_a = cmp; }
    }

    /* we still have a bracket */
    nf = 2;
    pos = pos_a;
    p = (const char *) base0 + pos * size;

FoundKey:
    *base0offset = pos;
    *num = nf;
    return (void *)p;
}
