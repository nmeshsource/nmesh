/* glist.c */
/* Wolfgang Tichy, 4/2023 */


#include "nmesh.h"


/**********************************************************************/
/* support for generic lists using linux lists */
/**********************************************************************/

/* alloc one element in a generic list, can be freed with just free() */
tGlist *glist_elem_alloc(void)
{
  tGlist *elem = malloc(sizeof(elem[0]));
  if(!elem) errorexit("no memory for elem");
  return elem;
}

/* make a new tGlist element, add it to list begin, and add entry into it */
void glist_entry_add(void *entry, struct list_head *head)
{
  tGlist *elem = glist_elem_alloc();
  elem->entry = entry;
  list_add(&elem->list, head);
}

/* make a new tGlist element, add it to list end, and add entry into it */
void glist_entry_add_tail(void *entry, struct list_head *head)
{
  tGlist *elem = glist_elem_alloc();
  elem->entry = entry;
  list_add_tail(&elem->list, head);
}

/* remove one tGlist element from list, and also free the element,
   but do not free its entry */
void glist_elem_del(tGlist *elem)
{
  list_del(&elem->list);
  free(elem);
}
