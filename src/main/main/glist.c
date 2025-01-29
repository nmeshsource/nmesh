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

/* remove and free the tGlist element from list which contains entry,
   but do not free entry */
void glist_entry_del(void *entry)
{
  tGlist *elem = container_of(&entry, tGlist, entry);
  errorexit("Too dangerous: It works ONLY if we use the pointer entry that "
            "is actually in tGlist. It fails if we use another one"
            " (e.g. elm)!!!");
  list_del(&elem->list);
  free(elem);
}

/* free all entries in glist andremove one tGlist element from list, and also free the element,
   but do not free its entry */
//was: void glist_free_elems_and_entries(struct list_head *head, void (*Free)())
void glist_free_elems_and_entries(struct list_head *head, void (*Free)(void *e))
{
  struct list_head *pos, *sav;
  list_for_each_prev_safe(pos, sav, head)
  {
    tGlist *elem = list_entry(pos, tGlist, list);
    void *entry = elem->entry;
    if(Free) Free(entry);
    glist_elem_del(elem);
  }
}

/* free all tGlist elements as to clear list, but do not free entries */
void glist_free_elems(struct list_head *head)
{
  glist_free_elems_and_entries(head, NULL);
}
