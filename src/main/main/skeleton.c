/* skeleton.c */
/* Wolfgang Tichy, 2019 */


#include "nmesh.h"
#include "main.h"

/* global var that tells us how often we have restarted nmesh */
extern int nmesh_restarts;


void AddMeshFun(tMesh *mesh, int step, int (*f)(tMesh *), char *name)
{
  tTodo *t;
  tTodo **skel = mesh->skel;

  if(1) printf("  func_%03d  %s\n", step, name);

  if(!skel[step]) skel[step] = calloc(1, sizeof(tTodo));

  for(t = skel[step]; t->next; t = t->next);
  t->next = calloc(1, sizeof(tTodo)); /* allocate next one for later use */
  t->f = f;
  t->name = strdup(name);
}

void remove_all_MeshFuns(tMesh *mesh)
{
  tTodo **skel = mesh->skel;
  int step;

  for(step=0; step<NFUNCBINS; step++)
  {
    tTodo *t, *p=NULL;

    for(t = skel[step]; t; )
    {
      //printf("step%d %p %s\n", step, t, t->name);
      p = t;         /* previous entry */
      t = t->next;   /* get next entry */
      free(p->name); /* free previous name */
      free(p);       /* free previous entry */
    }
  }
}

void RunMeshFun(tMesh *mesh, int step) 
{
  tTodo *t;
  tTodo **skel = mesh->skel;

  if(!skel[step]) return;

  for(t = skel[step]; t->next; t = t->next)
    (*(t->f))(mesh);
}
