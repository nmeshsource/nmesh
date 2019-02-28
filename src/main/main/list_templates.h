/* main/main/list_templates.h */
/* Wolfgang Tichy, May 2017 */
/* this header file contains macros & prototypes for generic lists */

/*******************************************************************/
/* essential helper macros. Only using them, the stuff below works */
#define DUMMY(a)      a
#define PASTE(a, b)   a ## b
#define PASTE_(a, b)  a ## _ ## b

/* make a function name */
#define FN(name,LTYPE) PASTE_(name,LTYPE)

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
LIST(TYP) *FN(alloc,LIST(TYP))(void);
void FN(free,LIST(TYP))(LIST(TYP) *u);
void FN(clear,LIST(TYP))(LIST(TYP) *u);
void FN(pr,LIST(TYP))(LIST(TYP) *u);
void FN(push,LIST(TYP))(LIST(TYP) *v, TYP vi);
void FN(pushlist,LIST(TYP))(LIST(TYP) *v, LIST(TYP) *u);
void FN(unionpush,LIST(TYP))(LIST(TYP) *v, TYP vi);
void FN(unionpushlist,LIST(TYP))(LIST(TYP) *u, LIST(TYP) *v);
void FN(dropindex,LIST(TYP))(LIST(TYP) *v, int ind);
void FN(drop,LIST(TYP))(LIST(TYP) *v, TYP vi);
void FN(droplastn,LIST(TYP))(LIST(TYP) *v, int n);
void FN(droplist,LIST(TYP))(LIST(TYP) *v, LIST(TYP) *u);
LIST(TYP) *FN(duplicate,LIST(TYP))(LIST(TYP) *v);
int FN(in,LIST(TYP))(LIST(TYP) *v, TYP vi);
int FN(index,LIST(TYP))(LIST(TYP) *v, TYP vi);
int FN(index_prop,LIST(TYP))(LIST(TYP) *v, int i0,
                             int (*prop)(void *obj, TYP vi), void *obj);
void FN(copy,LIST(TYP))(LIST(TYP) *dest, LIST(TYP) *src,
                        void (*copy)(void *obj, TYP d, TYP s), void *obj);
void FN(add,LIST(TYP))(LIST(TYP) *r, double ca, LIST(TYP) *a,
                       double cb, LIST(TYP) *b,
                       void (*add)(void *obj, TYP r, double ca, TYP a,
                                   double cb, TYP b), void *obj);
void FN(addto,LIST(TYP))(LIST(TYP) *r, double ca, LIST(TYP) *a,
                         void (*addto)(void *obj, TYP r, double ca, TYP a),
                         void *obj);

/* e.g.: if TYP = int 
   LIST(TYP) FN(alloc,LIST(TYP))(void)
   becomes
   intList alloc_intList(void) */

/************************************************************************/
/* macros to do things with lists like looping */
#include "list_loops.h"
