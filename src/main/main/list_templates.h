/* main/main/list_templates.h */
/* Wolfgang Tichy, May 2017 */
/* this header file contains macros & prototypes for generic lists */

/*******************************************************************/
/* essential helper macros. Only using them, the stuff below works */
#define DUMMY(a)      a
#define PASTE(a, b)   a ## b
#define PASTE_(a, b)  a ## _ ## b

/* make a function name */
#define FN(LTYPE,name) PASTE_(LTYPE,name)

/*******************************************************************/

/*********************/
/* for generic lists */
/*********************/
/* the new type name can be obtained like this. E.g. LIST(float) gives foatList */
#define LIST(TYP) PASTE(TYP,List)

/* define a new list with entries of a certain type TYP, 
   e.g. typedefLIST(double);
   would be expanded to:
   typedef struct doubleLIST {int n; double *e;} doubleList */
#define typedefLIST(TYP) typedef struct PASTE(TYP,LIST) { \
  int n;      /* number of entries */ \
  TYP *e;     /* array of list entries */ \
} PASTE(TYP,List)

/* generate a typedef to define the list type we need */
typedefLIST(TYP);

/************************************************************************/
/* utility functions for these lists */

/* from list_templates.c */
LIST(TYP) *FN(LIST(TYP),alloc)(void);
void FN(LIST(TYP),free)(LIST(TYP) *u);
void FN(LIST(TYP),clear)(LIST(TYP) *u);
void FN(LIST(TYP),pr)(LIST(TYP) *u);
void FN(LIST(TYP),push)(LIST(TYP) *v, TYP vi);
void FN(LIST(TYP),pushlist)(LIST(TYP) *v, LIST(TYP) *u);
void FN(LIST(TYP),unionpush)(LIST(TYP) *v, TYP vi);
void FN(LIST(TYP),unionpushlist)(LIST(TYP) *u, LIST(TYP) *v);
void FN(LIST(TYP),dropindex)(LIST(TYP) *v, int ind);
void FN(LIST(TYP),drop)(LIST(TYP) *v, TYP vi);
void FN(LIST(TYP),droplastn)(LIST(TYP) *v, int n);
void FN(LIST(TYP),droplist)(LIST(TYP) *v, LIST(TYP) *u);
void FN(LIST(TYP),setatindex)(LIST(TYP) *v, int ind, TYP vi);
LIST(TYP) *FN(LIST(TYP),duplicate)(LIST(TYP) *v);
int FN(LIST(TYP),in)(LIST(TYP) *v, TYP vi);
int FN(LIST(TYP),index)(LIST(TYP) *v, TYP vi);
int FN(LIST(TYP),index_prop)(LIST(TYP) *v, int i0,
                             int (*prop)(), //(const void *obj, TYP vi),
                             const void *obj);
void FN(LIST(TYP),copy)(LIST(TYP) *dest, LIST(TYP) *src,
                        void (*copy)(), //(const void *obj, TYP d, TYP s),
                        const void *obj);
void FN(LIST(TYP),add)(LIST(TYP) *r, double ca, LIST(TYP) *a,
                       double cb, LIST(TYP) *b,
                       void (*add)(), //(const void *obj, TYP r, double ca, TYP a, double cb, TYP b),
                       const void *obj);
void FN(LIST(TYP),addto)(LIST(TYP) *r, double ca, LIST(TYP) *a,
                         void (*addto)(), //(const void *obj, TYP r, double ca, TYP a),
                         const void *obj);
void FN(LIST(TYP),freeclear)(LIST(TYP) *r, void (*Free)(), const void *obj);
void FN(LIST(TYP),freeall)(LIST(TYP) *r, void (*Free)(), const void *obj);

/* e.g.: if TYP = int 
   LIST(TYP) FN(LIST(TYP),alloc)(void)
   becomes
   intList intList_alloc(void) */

/************************************************************************/
/* macros to do things with lists like looping */
#include "list_loops.h"
