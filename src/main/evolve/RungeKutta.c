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
  int redo_substep = Getv(Par("evolve_redo_troubled"), "substep");
  int nrs = 0; //counts how many substep we redone

  pVLList_copy(u_p, u, vlcopy,0);             // u_p = u
  mesh->time = t;
  evolve_setrhs_mesh(mesh, r, u);             // r  = RHS(u, t)
  pVLList_addto(u, dt/6., r, vladdto,0);      // u += r dt/6
  if(redo_substep)
  {
    pVLList_copy(w, u_p, vlcopy,0);           // w = u_p
    if(evolve_set_trouble_score_mesh(mesh)>0) //score u after step
    { evolve_trouble_redo_u_step_mesh(mesh, dt/6.);  nrs++; }
  }

  pVLList_add(w, 1., u_p, dt/2., r, vladd,0); // w  = u_p + r dt/2
  mesh->time = t+0.5*dt;
  evolve_limiter_mesh(mesh, w, 0);
  evolve_setrhs_mesh(mesh, r, w);             // r  = RHS(w, t+dt/2)
  pVLList_addto(u, dt/3., r, vladdto,0);      // u += r dt/3
  if(redo_substep)
  {
    if(evolve_set_trouble_score_mesh(mesh)>0) //score u after step
    { evolve_trouble_redo_u_step_mesh(mesh, dt/3.);  nrs++; }
  }

  pVLList_add(w, 1., u_p, dt/2., r, vladd,0); // w  = u_p + r dt/2
  mesh->time = t+0.5*dt;
  evolve_limiter_mesh(mesh, w, 0);
  evolve_setrhs_mesh(mesh, r, w);             // r  = RHS(w, t+dt/2)
  pVLList_addto(u, dt/3., r, vladdto,0);      // u += r dt/3
  if(redo_substep)
  {
    if(evolve_set_trouble_score_mesh(mesh)>0) //score u after step
    { evolve_trouble_redo_u_step_mesh(mesh, dt/3.);  nrs++; }
  }

  pVLList_add(w, 1., u_p, dt, r, vladd,0);    // w  = u_p + r dt
  mesh->time = t+dt;
  evolve_limiter_mesh(mesh, w, 0);
  evolve_setrhs_mesh(mesh, r, w);             // r  = RHS(w, t+dt)
  pVLList_addto(u, dt/6., r, vladdto,0);      // u += r dt/6
  if(redo_substep)
  {
    if(evolve_set_trouble_score_mesh(mesh)>0) //score u after step
    { evolve_trouble_redo_u_step_mesh(mesh, dt/6.);  nrs++; }
  }

  mesh->time = t+dt;                          // we are now at t+dt

  if(nrs>0) {PRF;printf(": switched troubled elms & redid %d substeps\n",nrs);}

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
  int redo_substep = Getv(Par("evolve_redo_troubled"), "substep");
  int nrs = 0; //counts how many substep we redone

  pVLList_copy(u_p, u, vlcopy,0);         // u_p = u
  mesh->time = t;
  evolve_setrhs_mesh(mesh, r, u);         // r  = RHS(u, t)
  pVLList_addto(u, dt, r, vladdto,0);     // u += r dt
  if(redo_substep)
  {
    pVLList_copy(w, u_p, vlcopy,0);       // w = u_p
    if(evolve_set_trouble_score_mesh(mesh)>0) //score u after step
    { evolve_trouble_redo_u_step_mesh(mesh, dt);  nrs++; }
  }

  mesh->time = t+dt;                      // we are now at t+dt

  if(nrs>0) {PRF;printf(": switched troubled elms & redid %d substeps\n",nrs);}

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
  int redo_substep = Getv(Par("evolve_redo_troubled"), "substep");
  int nrs = 0; //counts how many substeps were redone

  //tNode *node = node_from_nodename(mesh, "0_366");
  //tPoint pt[] =  {{.node=node, .ijk=17}};
  //if(PR) PRFs(": 0\n");
  //printvarlist_atpoint(pt, ListEntry(u,0), "");
  //printvarlist_atpoint(pt, ListEntry(u,1), "");
  pVLList_copy(u_p, u, vlcopy,0);              // u_p = u
  mesh->time = t;
  evolve_setrhs_mesh(mesh, r, u);              // r  = RHS(u, t)
  pVLList_addto(u, dt/6., r, vladdto,0);       // u += r dt/6
  if(redo_substep)
  {
    pVLList_copy(w, u_p, vlcopy,0);            // w = u_p
    if(evolve_set_trouble_score_mesh(mesh)>0)  //score u after step
    { evolve_trouble_redo_u_step_mesh(mesh, dt/6.);  nrs++; }
  }

  pVLList_add(w, 1., u_p, dt, r, vladd,0);     // w  = u_p + r dt
  mesh->time = t+dt;                           // c_2=1 from Butcher tab.
  evolve_limiter_mesh(mesh, w, 0);
  evolve_setrhs_mesh(mesh, r, w);              // r  = RHS(w, t+dt)
  pVLList_addto(u, dt/6., r, vladdto,0);       // u += r dt/6
  if(redo_substep)
  {
    if(evolve_set_trouble_score_mesh(mesh)>0) //score u after step
    { evolve_trouble_redo_u_step_mesh(mesh, dt/6.);  nrs++; }
  }


  pVLList_addto(w, dt, r, vladdto,0);          // w += r dt
  pVLList_add(w, 0.75, u_p, 0.25, w, vladd,0); // w = 0.75*u_p + 0.25*w
  mesh->time = t+0.5*dt;                       // c_3=1/2 from Butcher tab.
  evolve_limiter_mesh(mesh, w, 0);
  evolve_setrhs_mesh(mesh, r, w);              // r  = RHS(w, t+dt/2)
  pVLList_addto(u, dt*2./3., r, vladdto,0);    // u += r dt*2/3
  if(redo_substep)
  {
    if(evolve_set_trouble_score_mesh(mesh)>0) //score u after step
    { evolve_trouble_redo_u_step_mesh(mesh, dt*2./3.);  nrs++; }
  }

  mesh->time = t+dt;                           // we are now at t+dt

  if(nrs>0) {PRF;printf(": switched troubled elms & redid %d substeps\n",nrs);}

  /* The new u is not limited yet!
     A final evolve_limiter_mesh(mesh, u, 0) is called in evolve_myln */
}

/* return number of substeps an RK scheme takes */
int evolve_nsubsteps(tMesh *mesh)
{
  /* check evo method */
  int evolve_method = Par("evolve_method");
  if(Getv(evolve_method,      "RK4"))    return 4;
  else if(Getv(evolve_method, "Euler"))  return 1;
  else if(Getv(evolve_method, "sspRK3")) return 3;
  else
    errorexits("unknown value:   evolve_method = %s", Gets(evolve_method));
}

/* Return limit after which we consider switching back to dg. If we check
   for trouble only after full evo step we want to say all is great if the
   trouble score is less than -NOTROUBLES, but if we check after each
   substep it's great if it's less then -NOTROUBLES*evolve_nsubsteps(mesh). */
int evolve_notroubles(tMesh *mesh)
{
  int redo_substep = Getv(Par("evolve_redo_troubled"), "substep");
  if(redo_substep) return NOTROUBLES*evolve_nsubsteps(mesh);
  else             return NOTROUBLES;
}


/*************************************************************************/
/* functions to redo an RK substep */
/*************************************************************************/

/* Redo substep on all elms with elm->dat->info->trbl_score > 0.
   We 1st take pre substep of amount rfac back. Here rfac is e.g. dt/6 */
void evolve_trouble_redo_u_step_mesh(tMesh *mesh, double rfac)
{
  tEvoSys *evosys = mesh->evosys;
  pVLList *u   = evosys->u;
  pVLList *u_p = evosys->u_p;
  pVLList *r   = evosys->rhs;
  pVLList *w   = evosys->w;
  int notroubles = evolve_notroubles(mesh);

  /* alloc list to store old dg elms */
  tElm **elm_new = checked_calloc(mesh->nmyelm, sizeof(elm_new[0]));

  formyelms(mesh)
  {
    tElm *elm = MyElm;
    int trb = elm->dat->info->trbl_score;
    if(trb > 0)
    {
      tElm *elm_sav;
      int li;
      int n[3], pt_typ[3];

      // take substep back
      pVLList_addto(u, -rfac, r, vladdto, elm);  // u -= r rfac

      // pick n, pt_typ
      hp_refine_set_n_pt_typ(elm, elm->dat->info->trbl_ref, n, pt_typ);

      // make u_p, w DATAVARs so that they will be interp'd on p-refine
      forList(u_p, li)
      {
        VLSetType(ListEntry(u_p,li), DATAVAR);
        VLSetType(ListEntry(w,li), DATAVAR);
      }

      // p-refine locally
      elm_sav = update_node_n_pt_typ_return_node_old(elm, n, pt_typ);
      // ^-need to keep old elm in case something is pointing to its data
      //   So we store old elm in elm_sav.

      // Since elm is now refined we reset its ref method to PARENT_n, which
      // is a no-op. Then evolve_switch_troubled_nodes_mesh below will not
      // refine again.
      elm->dat->info->trbl_ref->method = PARENT_n;

      // make u_p, w AUXVAR again:
      forList(u_p, li)
      {
        VLSetType(ListEntry(u_p,li), AUXVAR);
        VLSetType(ListEntry(w,li), AUXVAR);
      }

      // remove all ajsurf of old elm
      free_all_ajsurf_only(elm_sav);

      // init surfaces in new elm, could use MPIexchange_init(elm)
      init_all_surfaces(elm);

      // set mysurf on new elm, could use MPIexchange_set_localdata(elm)
      set_all_mysurf(elm);

      // Let surfaces in new elm point to nbsurf of old elm.
      surface_copy_nbsurf_pointers(elm_sav, elm);

      // interp nb surfs to adj for our new elm
      set_all_ajsurf(elm);

      /* NOTE: elm->dat now has surfaces, but no indic */

      // init indicators for w in new elm
      init_myindc_for_evosys_u_or_w(elm, w);

      // let all indicators in new elm point to indicators of old elm.
      indic_copy_nbindc_pointers(elm_sav, elm);

      // run RHS funcs again
      evolve_limiter(elm, w, 0, notroubles, 0);//like evolve_limiter_mesh(mesh, w, 0);
      evolve_setrhs(elm, r, w, 0);

      // set u again
      pVLList_addto(u, rfac, r, vladdto, elm);  // u += r rfac

      // make elm_sav the new elm and save it in elm_new list
      elm_swap_shallow(elm, elm_sav);
      elm_new[MyID] = elm_sav;
      // Now elm is the old elm again. We need this for the next iter of
      // this loop, otherwise set_all_ajsurf gets confused.
      // Below, we swap the new elm back in!
    }
  }

  /* now insert all the elms in list elm_new, and free the old elms */
  formyelms(mesh)
  {
    tElm *elm = MyElm;
    tElm *elm2 = elm_new[MyID];
    if(elm2)
    {
      elm_swap_shallow(elm, elm2); //new elm is in elm, old elm is in elm2
      update_node_n_pt_typ_free_node_old(elm2); // free old elm2
    }
  }
  /* free list itself */
  free(elm_new);

  /* Finally we call p-refine so that nb-info such as fnb gets updated. But
     interp should not happen again since update_node_n_pt_typ (called form
     hp_refine_elms_if_rflag) will notice that n,pt_typ are set already. */
  evolve_switch_troubled_nodes_mesh(mesh);
  // ^-This calls evolve_init_communication_structs and prefine_elms_if_rflag
  //   prefine_elms_if_rflag updates nb-info such as fnb.
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
