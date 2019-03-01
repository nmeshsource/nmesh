/* main/main/list_templates.h */
/* Wolfgang Tichy, May 2017 */

#include "list_templates.h"
/* the list template is defined in list_templates.h */

/* To compile list functions with lists of type int, we need to use:
#define TYP int
#include "list_templates.c"
#undef TYP

To use them in another file we need to add
#define TYP int
#include "list_templates.h"
#undef TYP
to this file. */

/* If we need list of pointers do e.g. this for a list of pointers to tVarList:
typedef tVarList *pVL;        // list_templates.h only works with numbers
#define TYP pVL               // the pointer pVL is a number
#include "list_templates.c"
#undef TYP
*/

/************************************************************************/
/* utility functions for these lists */


/* allocate an empty list */
LIST(TYP) *FN(alloc,LIST(TYP))(void)
{
  LIST(TYP) *u;

  u = calloc(1, sizeof(LIST(TYP)));
  return u;
}

/* free a list */
void FN(free,LIST(TYP))(LIST(TYP) *u)
{
  if(u)
  {
    if(u->e) free(u->e);
    free(u);
    u = NULL;
  }
}

/* clear a list */
void FN(clear,LIST(TYP))(LIST(TYP) *u)
{
  if(u)
  {
    if(u->e) free(u->e);
    u->e = NULL;
    u->n = 0;
  }
}

/* print list u */
void FN(pr,LIST(TYP))(LIST(TYP) *u)
{
  int i;
  long ui;
  printf("%s: n=%d  e =", __func__, u->n);
  for(i=0; i<u->n; i++)
  {
    ui = (long) u->e[i];
    printf(" %ld" , ui);
  }
  printf("\n");
}

/* add an entry to a list */
void FN(push,LIST(TYP))(LIST(TYP) *v, TYP vi)
{
  v->n += 1;
  v->e = realloc(v->e, sizeof(TYP) * v->n);
  v->e[v->n-1] = vi;
}

/* append a list to a list */
void FN(pushlist,LIST(TYP))(LIST(TYP) *v, LIST(TYP) *u)
{
  int i;

  if(!v || !u) return;
  v->n += u->n;
  v->e = realloc(v->e, sizeof(TYP) * v->n);
  for (i = 0; i < u->n; i++)
    v->e[v->n - u->n + i] = u->e[i];
}

/* add to list (if not already in it) */
void FN(unionpush,LIST(TYP))(LIST(TYP) *v, TYP vi)
{
  int i;
  int addvi=1;
  /* add vi only if it is not already in blist */
  for(i=0; i<v->n; i++) if(v->e[i]==vi) { addvi=0; break; }
  if(addvi) FN(push,LIST(TYP))(v, vi);
}

/* v = union(v, u): add all of u to list v (if not already in v) */
void FN(unionpushlist,LIST(TYP))(LIST(TYP) *v, LIST(TYP) *u)
{
  int i;
  for(i=0; i<u->n; i++) FN(unionpush,LIST(TYP))(v, u->e[i]);
}

/* drop an entry from a list */
void FN(dropindex,LIST(TYP))(LIST(TYP) *v, int ind)
{
  int i;
  if(ind<0 || ind >= v->n) return;
  for(i = ind; i < v->n-1; i++)  v->e[i] = v->e[i+1];
  v->n -= 1;
}

/* drop an entry from a list */
void FN(drop,LIST(TYP))(LIST(TYP) *v, TYP vi)
{
  int i;

  for(i = 0; i < v->n; i++)
    if (v->e[i] == vi)
    {
      v->n -= 1;
      for (; i < v->n; i++)
        v->e[i] = v->e[i+1];
      break;
    }
}

/* drop last n entries from a list */
void FN(droplastn,LIST(TYP))(LIST(TYP) *v, int n)
{
  if(n <= 0)
    return;
  if(n >= v->n)
    v->n = 0;
  else
    v->n -= n;
}

/* drop all in u from v */
void FN(droplist,LIST(TYP))(LIST(TYP) *v, LIST(TYP) *u)
{
  int i;
  for(i=0; i<u->n; i++) FN(drop,LIST(TYP))(v, u->e[i]);
}

/* duplicate a list */
LIST(TYP) *FN(duplicate,LIST(TYP))(LIST(TYP) *v)
{
  int i;
  LIST(TYP) *u = FN(alloc,LIST(TYP))();

  for(i = 0; i < v->n; i++)
    FN(push,LIST(TYP))(u, v->e[i]);

  return u;
}

/* return 1 if vi is in list v */
int FN(in,LIST(TYP))(LIST(TYP) *v, TYP vi)
{
  int i;
  int in=0;
  for(i=0; i<v->n; i++) if(v->e[i]==vi) { in=1; break; }
  return in;
}

/* return index of first element vi in list v, returns -1 if not in list */
int FN(index,LIST(TYP))(LIST(TYP) *v, TYP vi)
{
  int i;
  int in=-1; /* is not in list */
  for(i=0; i<v->n; i++) if(v->e[i]==vi) { in=i; break; }
  return in;
}

/* return index of first element in list v that has prop returning 1,
   we start checking with element i0, returns -1 if prop returns 0 for all in v */
/* the function prop could be as simple as:
   int prop(void *p, int vi)
   {
     int *pi = (int *) p;
     return (vi == *pi);
   }
*/
int FN(index_prop,LIST(TYP))(LIST(TYP) *v, int i0,
                             int (*prop)(), //(void *obj, TYP vi),
                             void *obj)
{
  int i;
  int in=-1; /* is not in list */
  if(i0<0) i0=0;
  for(i=i0; i<v->n; i++) if(prop(obj, v->e[i])) { in=i; break; }
  return in;
}

/* copy contents of src into dest, so that dest = src
   We need to pass in a func copy that know how to copy list elements.
   e.g. for ints it could be just:
   void copy(int d, int s){ d = s; }
*/
void FN(copy,LIST(TYP))(LIST(TYP) *dest, LIST(TYP) *src,
                        void (*copy)(), //(void *obj, TYP d, TYP s),
                        void *obj)
{
  int i;
  for(i=0; i<dest->n; i++)
    copy(obj, dest->e[i], src->e[i]); /* func that know how to copy */
}

/* add contents: r = ca*a + cb*b */
void FN(add,LIST(TYP))(LIST(TYP) *r, double ca, LIST(TYP) *a,
                       double cb, LIST(TYP) *b,
                       void (*add)(), //(void *obj, TYP r, double ca, TYP a, double cb, TYP b),
                       void *obj)
{
  int i;
  for(i=0; i<r->n; i++)
    add(obj, r->e[i], ca,a->e[i], cb,b->e[i]); /* func that adds */
}

/* add to contents: r += ca*a */
void FN(addto,LIST(TYP))(LIST(TYP) *r, double ca, LIST(TYP) *a,
                         void (*addto)(), //(void *obj, TYP r, double ca, TYP a),
                         void *obj)
{
  int i;
  for(i=0; i<r->n; i++)
    addto(obj, r->e[i], ca,a->e[i]); /* func that adds to r */
}

/* free contents of r */
void FN(freeclear,LIST(TYP))(LIST(TYP) *r, void (*Free)(), void *obj)
{
  if(Free && r)
  {
    int i;
    for(i=0; i<r->n; i++)
    {
      if(obj) Free(obj, r->e[i]); /* func that frees r */
      else    Free(r->e[i]);
    }
  }
  FN(clear,LIST(TYP))(r);
}

/* free contents of r and then r itself */
void FN(freeall,LIST(TYP))(LIST(TYP) *r, void (*Free)(), void *obj)
{
  FN(freeclear,LIST(TYP))(r, Free, obj);
  FN(free,LIST(TYP))(r);
}

/***************************************************************************/
/* undef all from list_templates.h */
#undef DUMMY
#undef PASTE
#undef PASTE_
#undef FN
#undef LIST
#undef typedefLIST
#undef TYP
