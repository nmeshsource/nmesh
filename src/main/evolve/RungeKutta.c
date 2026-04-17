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
  int trouble_score;

  pVLList_copy(u_p, u, vlcopy,0);             // u_p = u
  mesh->time = t;
  evolve_setrhs_mesh(mesh, r, u);             // r  = RHS(u, t)
  pVLList_addto(u, dt/6.0, r, vladdto,0);     // u += r dt/6
  trouble_score = evolve_set_trouble_score_mesh(mesh); // scores u
  if(trouble_score>0) evolve_trouble_redo_u_step_mesh(mesh, dt/6.0);

  pVLList_add(w, 1., u_p, dt/2., r, vladd,0); // w  = u_p + r dt/2
  mesh->time = t+0.5*dt;
  evolve_limiter_mesh(mesh, w, 0);





  evolve_setrhs_mesh(mesh, r, w);             // r  = RHS(w, t+dt/2)
  pVLList_addto(u, dt/3., r, vladdto,0);      // u += r dt/3

  pVLList_add(w, 1., u_p, dt/2., r, vladd,0); // w  = u_p + r dt/2
  mesh->time = t+0.5*dt;
  evolve_limiter_mesh(mesh, w, 0);
  evolve_setrhs_mesh(mesh, r, w);             // r  = RHS(w, t+dt/2)
  pVLList_addto(u, dt/3., r, vladdto,0);      // u += r dt/3

  pVLList_add(w, 1., u_p, dt, r, vladd,0);    // w  = u_p + r dt
  mesh->time = t+dt;
  evolve_limiter_mesh(mesh, w, 0);
  evolve_setrhs_mesh(mesh, r, w);             // r  = RHS(w, t+dt)
  pVLList_addto(u, dt/6., r, vladdto,0);      // u += r dt/6
  mesh->time = t+dt;                          // we are now at t+dt

  /* switch from fv to dg if we had NOTROUBLES*3 RK substeps without trouble */
  evolve_switch_nontroubled_nodes_mesh(mesh, NOTROUBLES*3);
  // apply limiter
  // FIXME: call limiter


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

  pVLList_copy(u_p, u, vlcopy,0);         // u_p = u
  mesh->time = t;
  evolve_setrhs_mesh(mesh, r, u);         // r  = RHS(u, t)
  pVLList_addto(u, dt, r, vladdto,0);     // u += r dt
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
  pVLList_copy(u_p, u, vlcopy,0);              // u_p = u
  mesh->time = t;
  evolve_setrhs_mesh(mesh, r, u);              // r  = RHS(u, t)
  pVLList_addto(u, dt/6., r, vladdto,0);       // u += r dt/6

  pVLList_add(w, 1., u_p, dt, r, vladd,0);     // w  = u_p + r dt
  mesh->time = t+dt;                           // c_2=1 from Butcher tab.
  //printvarlist_atpoint(pt, ListEntry(w,0), "");
  evolve_limiter_mesh(mesh, w, 0);
  if(PR) PRFs(": 1\n");
  //printvarlist_atpoint(pt, ListEntry(w,0), "");
  evolve_setrhs_mesh(mesh, r, w);              // r  = RHS(w, t+dt)
  pVLList_addto(u, dt/6., r, vladdto,0);       // u += r dt/6

  pVLList_addto(w, dt, r, vladdto,0);          // w += r dt
  pVLList_add(w, 0.75, u_p, 0.25, w, vladd,0); // w = 0.75*u_p + 0.25*w
  mesh->time = t+0.5*dt;                       // c_3=1/2 from Butcher tab.
  //printvarlist_atpoint(pt, ListEntry(w,0), "");
  evolve_limiter_mesh(mesh, w, 0);
  if(PR) PRFs(": 2\n");
  //printvarlist_atpoint(pt, ListEntry(w,0), "");
  evolve_setrhs_mesh(mesh, r, w);              // r  = RHS(w, t+dt/2)
  pVLList_addto(u, dt*2./3., r, vladdto,0);    // u += r dt*2/3
  mesh->time = t+dt;                           // we are now at t+dt
  if(PR) PRFs(": 3\n");
  //printvarlist_atpoint(pt, ListEntry(u,0), "");
  /* The new u is not limited yet!
     A final evolve_limiter_mesh(mesh, u, 0) is called in evolve_myln */
}


/*************************************************************************/
/* functions to redo an RK substep */
/*************************************************************************/

void evolve_trouble_redo_u_step_mesh(tMesh *mesh, double rfac)
{
  tEvoSys *evosys = mesh->evosys;
  pVLList *u   = evosys->u;
  pVLList *u_p = evosys->u_p;
  pVLList *r   = evosys->rhs;

  /* alloc list to store old dg elms */
  tElm **elm_old = checked_calloc(mesh->nmyelm, sizeof(elm_old[0]));

  formyelms(mesh)
  {
    tElm *elm = MyElm;
    int n[3], pt_typ[3];
    int trb = elm->dat->info->trbl_score;
    if(trb > 0)
    {
      int li;
      tRef *ref = elm->dat->info->trbl_ref;

      // take substep back
      pVLList_addto(u, -rfac, r, vladdto, elm);  // e.g. u -= r dt/6

      // pick n, pt_typ
      hp_refine_set_n_pt_typ(elm, ref, n, pt_typ);

      // make u_p a DATAVAR so that it will be interp'd on p-refine
      forList(u_p, li) VLSetType(ListEntry(u_p,li), DATAVAR);

      // p-refine locally
      elm_old[MyID] = update_node_n_pt_typ_return_node_old(elm, n, pt_typ);
      // ^-need to keep old elm in case something is pointing to its data
      //   So we store old elm in list elm_old and free it later.

      // make u_p an AUXVAR again:
      forList(u_p, li) VLSetType(ListEntry(u_p,li), AUXVAR);

      // remove all ajsurf of elm_old
      free_all_ajsurf_only(elm_old[MyID]);

      // copy all surface data pointers from elm_old to elm
      surface_copy_all_pointers(elm_old[MyID], elm);

      // interp nb surfs to adj for our new elm
      set_all_ajsurf(elm);

      // run RHS funcs again
      /*
      ???
      todo: 1. rewrite evolve_setrhs_mesh(mesh), s.t. it calls a funcs
               for each elm in each of the formylnodes(mesh) loops
            2. make evolve_setrhs(elm) that calls these same elm based funcs
      */

      // set u again
      pVLList_addto(u, rfac, r, vladdto, elm);  // e.g. u += r dt/6
    }
  }
  /* Finally we call p-refine so that nb-info such as fnb gets updated. But
     interp should not happen again since update_node_n_pt_typ (called form
     hp_refine_elms_if_rflag) will notice that n,pt_typ are set already. */
  evolve_switch_troubled_nodes_mesh(mesh);

  /* now free all the elms in list elm_old */
  formyelms(mesh)
    if(elm_old[MyID]) update_node_n_pt_typ_free_node_old(elm_old[MyID]);
  /* free list itself */
  free(elm_old);
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

  pVLList_copy(u_p, u, vlcopy, node);                  // u_p = u
  node->time = t;
  evolve_setrhs(node, r, u, 1);                        // r  = RHS(u, t)
  pVLList_addto(u, dt/6.0, r, vladdto, node);          // u += r dt/6

  pVLList_add(w, 1., u_p, dt/2., r, vladd, node);      // w  = u_p + r dt/2
  node->time = t+0.5*dt;
  evolve_setrhs(node, r, w, 1);                        // r  = RHS(w, t+dt/2)
  pVLList_addto(u, dt/3., r, vladdto, node);           // u += r dt/3

  pVLList_add(w, 1., u_p, dt/2., r, vladd, node);      // w  = u_p + r dt/2
  node->time = t+0.5*dt;
  evolve_setrhs(node, r, w, 1);                        // r  = RHS(w, t+dt/2)
  pVLList_addto(u, dt/3., r, vladdto, node);           // u += r dt/3

  pVLList_add(w, 1., u_p, dt, r, vladd, node);         // w  = u_p + r dt
  node->time = t+dt;
  evolve_setrhs(node, r, w, 1);                        // r  = RHS(w, t+dt)
  pVLList_addto(u, dt/6., r, vladdto, node);           // u += r dt/6
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

  pVLList_copy(u_p, u, vlcopy, node);              // u_p = u
  node->time = t;
  evolve_setrhs(node, r, u, 1);                    // r  = RHS(u, t)
  pVLList_addto(u, dt, r, vladdto, node);          // u += r dt
}
