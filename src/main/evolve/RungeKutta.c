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
  pVLList *r   = evosys->rhs;
  pVLList *w   = evosys->w;

  copy_pVLList(u_p, u, vlcopy_node, node);             // u_p = u

  node->time = t;
  evosys->setrhs(node, r, u);                          // r  = RHS(u, t)
  addto_pVLList(u, dt/6.0, r, vladdto_node, node);     // u += r dt/6

  add_pVLList(w, 1., u_p, dt/2., r, vladd_node, node); // w  = u_p + r dt/2
  node->time = t+0.5*dt;
  evosys->setrhs(node, r, w);                          // r  = RHS(w, t+dt/2)
  addto_pVLList(u, dt/3., r, vladdto_node, node);      // u += r dt/3

  add_pVLList(w, 1., u_p, dt/2., r, vladd_node, node); // w  = u_p + r dt/2
  node->time = t+0.5*dt;
  evosys->setrhs(node, r, w);                          // r  = RHS(w, t+dt/2)
  addto_pVLList(u, dt/3., r, vladdto_node, node);      // u += r dt/3

  add_pVLList(w, 1., u_p, dt, r, vladd_node, node);    // w  = u_p + r dt
  node->time = t+dt;
  evosys->setrhs(node, r, w);                          // r  = RHS(w, t+dt)
  addto_pVLList(u, dt/6., r, vladdto_node, node);      // u += r dt/6
}
