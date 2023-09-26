/* evolve.c */
/* Wolfgang Tichy, 2/2019 */

#include "nmesh.h"
#include "evolve.h"

#define PR 0


/*************************************************************************/
/* functions to evolve all the mesh with a uniform time step
*/
/*************************************************************************/

/* Runge-Kutta 4 */
void evolve_RK4_mesh(tMesh *mesh)
{
  tEvoSys *evosys = mesh->evosys;
  double  t = mesh->time;
  double dt = mesh->dt;
  pVLList *u   = evosys->u;
  pVLList *u_p = evosys->u_p;
  pVLList *r   = evosys->rhs;
  pVLList *w   = evosys->w;

  copy_pVLList(u_p, u, vlcopy,0);             // u_p = u
  mesh->time = t;
  evolve_setrhs_mesh(mesh, r, u);             // r  = RHS(u, t)
  addto_pVLList(u, dt/6.0, r, vladdto,0);     // u += r dt/6

  add_pVLList(w, 1., u_p, dt/2., r, vladd,0); // w  = u_p + r dt/2
  mesh->time = t+0.5*dt;
  evolve_limiter_mesh(mesh, w, 0);
  evolve_setrhs_mesh(mesh, r, w);             // r  = RHS(w, t+dt/2)
  addto_pVLList(u, dt/3., r, vladdto,0);      // u += r dt/3

  add_pVLList(w, 1., u_p, dt/2., r, vladd,0); // w  = u_p + r dt/2
  mesh->time = t+0.5*dt;
  evolve_limiter_mesh(mesh, w, 0);
  evolve_setrhs_mesh(mesh, r, w);             // r  = RHS(w, t+dt/2)
  addto_pVLList(u, dt/3., r, vladdto,0);      // u += r dt/3

  add_pVLList(w, 1., u_p, dt, r, vladd,0);    // w  = u_p + r dt
  mesh->time = t+dt;
  evolve_limiter_mesh(mesh, w, 0);
  evolve_setrhs_mesh(mesh, r, w);             // r  = RHS(w, t+dt)
  addto_pVLList(u, dt/6., r, vladdto,0);      // u += r dt/6
  mesh->time = t+dt;                          // we are now at t+dt
  /* The new u is not limited yet!
     A final evolve_limiter_mesh(mesh, u, 0) is called in evolve_myln */
}

/* Euler step */
void evolve_Euler_mesh(tMesh *mesh)
{
  tEvoSys *evosys = mesh->evosys;
  double  t = mesh->time;
  double dt = mesh->dt;
  pVLList *u   = evosys->u;
  pVLList *u_p = evosys->u_p;
  pVLList *r   = evosys->rhs;
  //pVLList *w   = evosys->w;

  copy_pVLList(u_p, u, vlcopy,0);         // u_p = u
  mesh->time = t;
  evolve_setrhs_mesh(mesh, r, u);         // r  = RHS(u, t)
  addto_pVLList(u, dt, r, vladdto,0);     // u += r dt
  mesh->time = t+dt;                      // we are now at t+dt
  /* The new u is not limited yet!
     A final evolve_limiter_mesh(mesh, u, 0) is called in evolve_myln */
}

/* third order strong stability preserving Runge-Kutta scheme from
   arXiv:1804.02003v2 .
   Since arXiv:1804.02003v2 does not mention it, mesh->time = t + c_i*dt
   is set using the 1st column of the Butcher tableau for SSPRK3
             c_i   a_{ij}       i,j \in [1,2,3,...]
              /       \
             0   | 0     0     0
             1   | 1     0     0
             1/2 | 1/4   1/4   0
             ----------------------
     b_i ->      | 1/6   1/6   2/3
   in https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods
*/
void evolve_sspRK3_mesh(tMesh *mesh)
{
  tEvoSys *evosys = mesh->evosys;
  double  t = mesh->time;
  double dt = mesh->dt;
  pVLList *u   = evosys->u;
  pVLList *u_p = evosys->u_p;
  pVLList *r   = evosys->rhs;
  pVLList *w   = evosys->w;

  //tNode *node = node_from_nodename(mesh, "0_366");
  //tPoint pt[] =  {{.node=node, .ijk=17}};

  if(PR) PRFs(": 0\n");
  //printvarlist_atpoint(pt, ListEntry(u,0), "");
  //printvarlist_atpoint(pt, ListEntry(u,1), "");
  copy_pVLList(u_p, u, vlcopy,0);              // u_p = u
  mesh->time = t;
  evolve_setrhs_mesh(mesh, r, u);              // r  = RHS(u, t)
  addto_pVLList(u, dt/6., r, vladdto,0);       // u += r dt/6

  add_pVLList(w, 1., u_p, dt, r, vladd,0);     // w  = u_p + r dt
  mesh->time = t+dt;                           // c_2=1 from Butcher tab.
  //printvarlist_atpoint(pt, ListEntry(w,0), "");
  evolve_limiter_mesh(mesh, w, 0);
  if(PR) PRFs(": 1\n");
  //printvarlist_atpoint(pt, ListEntry(w,0), "");
  evolve_setrhs_mesh(mesh, r, w);              // r  = RHS(w, t+dt)
  addto_pVLList(u, dt/6., r, vladdto,0);       // u += r dt/6

  addto_pVLList(w, dt, r, vladdto,0);          // w += r dt
  add_pVLList(w, 0.75, u_p, 0.25, w, vladd,0); // w = 0.75*u_p + 0.25*w
  mesh->time = t+0.5*dt;                       // c_3=1/2 from Butcher tab.
  //printvarlist_atpoint(pt, ListEntry(w,0), "");
  evolve_limiter_mesh(mesh, w, 0);
  if(PR) PRFs(": 2\n");
  //printvarlist_atpoint(pt, ListEntry(w,0), "");
  evolve_setrhs_mesh(mesh, r, w);              // r  = RHS(w, t+dt/2)
  addto_pVLList(u, dt*2./3., r, vladdto,0);    // u += r dt*2/3
  mesh->time = t+dt;                           // we are now at t+dt
  if(PR) PRFs(": 3\n");
  //printvarlist_atpoint(pt, ListEntry(u,0), "");
  /* The new u is not limited yet!
     A final evolve_limiter_mesh(mesh, u, 0) is called in evolve_myln */
}



/*************************************************************************/
/* NOTE: functions below do not work yet !!! */
/*************************************************************************/

/*************************************************************************/
/* functions to evolve on just one node
   will work only once request_all_vl_surfaces and such in
   main/amr/surface.c start working
*/
/*************************************************************************/

/* Runge-Kutta 4 */
void evolve_RK4(tNode *node)
{
  tMesh *mesh = node->pat->mesh;
  tEvoSys *evosys = mesh->evosys;
  double  t = mesh->time;
  double dt = mesh->dt;
  pVLList *u   = evosys->u;
  pVLList *u_p = evosys->u_p;
  pVLList *r   = evosys->rhs;
  pVLList *w   = evosys->w;

  copy_pVLList(u_p, u, vlcopy_node, node);             // u_p = u
  node->time = t;
  evolve_setrhs(node, r, u, 1);                        // r  = RHS(u, t)
  addto_pVLList(u, dt/6.0, r, vladdto_node, node);     // u += r dt/6

  add_pVLList(w, 1., u_p, dt/2., r, vladd_node, node); // w  = u_p + r dt/2
  node->time = t+0.5*dt;
  evolve_setrhs(node, r, w, 1);                        // r  = RHS(w, t+dt/2)
  addto_pVLList(u, dt/3., r, vladdto_node, node);      // u += r dt/3

  add_pVLList(w, 1., u_p, dt/2., r, vladd_node, node); // w  = u_p + r dt/2
  node->time = t+0.5*dt;
  evolve_setrhs(node, r, w, 1);                        // r  = RHS(w, t+dt/2)
  addto_pVLList(u, dt/3., r, vladdto_node, node);      // u += r dt/3

  add_pVLList(w, 1., u_p, dt, r, vladd_node, node);    // w  = u_p + r dt
  node->time = t+dt;
  evolve_setrhs(node, r, w, 1);                        // r  = RHS(w, t+dt)
  addto_pVLList(u, dt/6., r, vladdto_node, node);      // u += r dt/6
}

/* Euler step */
void evolve_Euler(tNode *node)
{
  tMesh *mesh = node->pat->mesh;
  tEvoSys *evosys = mesh->evosys;
  double  t = mesh->time;
  double dt = mesh->dt;
  pVLList *u   = evosys->u;
  pVLList *u_p = evosys->u_p;
  pVLList *r   = evosys->rhs;
  //pVLList *w   = evosys->w;

  copy_pVLList(u_p, u, vlcopy_node, node);         // u_p = u
  node->time = t;
  evolve_setrhs(node, r, u, 1);                    // r  = RHS(u, t)
  addto_pVLList(u, dt, r, vladdto_node, node);     // u += r dt
}
