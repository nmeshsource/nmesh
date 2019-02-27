/* evolve.c */
/* Wolfgang Tichy, 2/2019 */

#include "nmesh.h"
#include "evolve.h"

#define PR 1


/* Runge-Kutta 4 */
void evolve_RK4(tNode *node)
{
  tMesh *mesh = node->pat->mesh;
  double  t = mesh->time;
  double dt = mesh->dt;
  pVLList *u_p = evosys->u_p;
  pVLList *r   = evosys->uk[0];
  pVLList *w   = evosys->uk[1];

  copy_pVLList(u_p, u, vlcopy);              /* u_p = u */

  evosys->setrhs(r, u, t);                   /* r  = RHS(u, t)       */
  addto_pVLList(u, dt/6.0, r, vladdto);      /* u += r dt/6          */

  add_pVLList(w, 1., u_p, dt/2., r, vladd);  /* w  = u_p + r dt/2    */
  evosys->setrhs(r, w, t+0.5*dt);            /* r  = RHS(w, t+dt/2)  */
  addto_pVLList(u, dt/3., r, vladdto);       /* u += r dt/3          */

  add_pVLList(w, 1., u_p, dt/2., r, vladd);  /* w  = u_p + r dt/2    */
  evosys->setrhs(r, w, t+0.5*dt);            /* r  = RHS(w, t+dt/2)  */
  addto_pVLList(u, dt/3., r, vladdto);       /* u += r dt/3          */
 
  add_pVLList(w, 1., u_p, dt, r, vladd);     /* w  = u_p + r dt      */
  evosys->setrhs(r, w, t+dt);                /* r  = RHS(w, t+dt)    */
  addto_pVLList(u, dt/6., r, vladdto);       /* u += r dt/6          */
}
