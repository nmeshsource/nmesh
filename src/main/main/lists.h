/* lists.h */
/* Wolfgang Tichy, 3/2019 */

/**************************************************************************/
/* header to use the functions compiled in lists.c */
/**************************************************************************/

/* use lists with int entries */
#define TYP int
#include "list_templates.h"
#undef TYP

/* use lists with (tVarList *) entries */
typedef tVarList *pVL;        /* list_templates.h only works with numbers */
#define TYP pVL               /* the pointer pVL is a number */
#include "list_templates.h"
#undef TYP

/* use lists with entries of type func. pointer */
typedef int (*FuncPointer)();  /* a func. pointer */
#define TYP FuncPointer        /* the pointer FuncPointer is a number */
#include "list_templates.h"
#undef TYP

/* use lists with (tArray *) entries */
typedef tArray *pArr;        /* list_templates.h only works with numbers */
#define TYP pArr             /* the pointer pArr is a number */
#include "list_templates.h"
#undef TYP

/* use lists with (const char *) entries */
typedef const char *constString; /* list_templates.h only works with numbers */
#define TYP constString          /* the pointer String is a number */
#include "list_templates.h"
#undef TYP
