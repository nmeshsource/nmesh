/* skeleton.c */
/* Wolfgang Tichy 1/2019, Bernd Bruegmann 12/99 */


#include "nmesh.h"
#include "main.h"

/* global var that tells us how often we have restarted nmesh */
extern int nmesh_restarts;

/* add a function to a function bin (or time bin) */
void AddMeshFun(tMesh *mesh, int step, int (*f)(tMesh *), const char *name)
{
  tTodo *t;
  tTodo **skel = mesh->skel;

  if(1) printf("  func_T%02d  %s\n", step, name);

  if(!skel[step]) skel[step] = calloc(1, sizeof(tTodo));

  for(t = skel[step]; t->next; t = t->next);
  t->next = calloc(1, sizeof(tTodo)); /* allocate next one for later use */
  t->f = f;
  t->name = strdup(name);
}

/* remove all functions from all func bins */
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

/* run functions in func bin step */
void RunMeshFun(tMesh *mesh, int step) 
{
  tTodo *t;
  tTodo **skel = mesh->skel;

  if(!skel[step]) return;

  for(t = skel[step]; t->next; t = t->next)
    (*(t->f))(mesh);
}

/* print functions in func bin step */
void PrintMeshFun(tMesh *mesh, int step)
{
  tTodo *t;
  tTodo **skel = mesh->skel;

  if(!skel[step]) return;

  printf("func_T%02d:\n", step);
  for(t = skel[step]; t->next; t = t->next)
  {
    printf("          %s\n", t->name);
    if(t->f == evolve_myln) evolve_print_evosys(mesh);
  }
}

/* print functions in all function bins */
int PrintMeshFuncs(tMesh *mesh)
{
  int step;

  prdivider(0);
  printf("Function bin skeleton: functions are called in this order:\n");
  for(step=0; step<NFUNCBINS; step++) PrintMeshFun(mesh, step);
  prdivider(0);
  return 0;
}
