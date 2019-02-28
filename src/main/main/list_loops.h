/* list_loops.h */
/* Wolfgang Tichy, 2/2019 */

/****************************************************************************/
/* Loops are performed by macros so that the user has to know very little
   about the implementation details. */
/****************************************************************************/

/* loop over List */
#define forList(List,k) for(k=0; k<List->n; k++)

/* get entry number i in List */
#define ListEntry(List, i) List->e[i]
