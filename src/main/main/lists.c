/* lists.c */
/* Wolfgang Tichy, May 2017 */

#include "sources.h"

#include <stdlib.h>
#include <stdio.h>


/***************************************************************************/
/* compile lists with entries of type int */
/***************************************************************************/
/* compile list functions with lists of type int
   e.g. intList *alloc_intList(void);  */
#define TYP int
#include "list_templates.c"
#undef TYP
/* to use them in a file we need to add
#define TYP int
#include "list_templates.h"
#undef TYP
to this file. */

/* if we need the same lists but with double entries, do this:
#define TYP double
#include "list_templates.c"
#undef TYP
*/

/***************************************************************************/
/* compile lists with entries of type (tVarList *) */
/***************************************************************************/
#include "../amr/nmesh_amr.h" /* get def of tVarList and tArray */
typedef tVarList *pVL;        /* list_templates.h only works with numbers */
#define TYP pVL               /* the pointer pVL is a number */
#include "list_templates.c"
#undef TYP

/***************************************************************************/
/* compile lists with entries of type func. pointer */
/***************************************************************************/
typedef int (*EvoFuncPtr)();  /* a func. pointer */
#define TYP EvoFuncPtr        /* the pointer EvoFuncPtr is a number */
#include "list_templates.c"
#undef TYP

/***************************************************************************/
/* compile lists with entries of type (tArray *) */
/***************************************************************************/
typedef tArray *pArr;        /* list_templates.h only works with numbers */
#define TYP pArr             /* the pointer pArr is a number */
#include "list_templates.c"
#undef TYP

/***************************************************************************/
/* compile lists with entries of type (const char *), i.e. strings */
/***************************************************************************/
typedef const char *constString; /* list_templates.h only works with numbers */
#define TYP constString          /* the pointer String is a number */
#include "list_templates.c"
#undef TYP
