/* skeleton.c */
/* Wolfgang Tichy, 2019 */


#include "nmesh.h"
#include "main.h"

/* global var that tells us how often we have restarted nmesh */
extern int nmesh_restarts;


void AddFun(tMesh *mesh, int step, int (*f)(tMesh *), char *name)
{
  tTodo *t;
  tTodo *skel = mesh->skel;

  if (1) printf("  function  %s\n", name);

  if (!skel[step]) skel[step] = (tTodo *) calloc(1, sizeof(tTodo));

  for (t = skel[step]; t->next; t = t->next);
  t->next = (tTodo *) calloc(1, sizeof(tTodo));
  t->f = f;
  t->name = (char *) calloc(strlen(name)+1, sizeof(char));
  strcpy(t->name, name);
}


void RunFun(tMesh *mesh, int step) 
{
  tTodo *t;
  tTodo *skel = mesh->skel;

  if (!skel[step]) return;

  for (t = skel[step]; t->next; t = t->next) {
    (*(t->f))(mesh);
  }
}
