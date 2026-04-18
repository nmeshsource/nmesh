/* evolve.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for coupled evolution systems */



/* evolve.c */
void evolve_setrhs_mesh(tMesh *mesh, pVLList *rhs, pVLList *u);
void evolve_setrhs_PRESURF(tElm *node, pVLList *rhs, pVLList *u);
void evolve_setrhs_VOLRHS(tElm *node, pVLList *rhs, pVLList *u);
void evolve_setrhs_SURFRHS(tElm *node, pVLList *rhs, pVLList *u);
void evolve_setrhs(tElm *elm, pVLList *rhs, pVLList *u, int MPI_exchange);
void evolve_setsrc_again_nontroubled_nodes_mesh(tMesh *mesh,
                                                pVLList *rhs, pVLList *u,
                                                int notroubles);
int evolve_filter_evosys_mesh(tMesh *mesh);
int EVOLVE_timer_start(tMesh *mesh);
int EVOLVE_timer_stop(tMesh *mesh);
int evolve_output_timers(tMesh *mesh);


/* evosys.c */
int evolve_free_evosys(tMesh *mesh);
int evolve_init_evosys(tMesh *mesh);
void evolve_request_surfaces(tNode *node, pVLList *u);


/* evolve_test.c */
int evolve_test_init(tMesh *mesh);
int evolve_test_analyze(tMesh *mesh);


/* RungeKutta.c */
void evolve_RK4_mesh(tMesh *mesh);
void evolve_Euler_mesh(tMesh *mesh);
void evolve_sspRK3_mesh(tMesh *mesh);
int evolve_notroubles(tMesh *mesh);
void evolve_trouble_redo_u_step_mesh(tMesh *mesh, double rfac);
void evolve_RK4(tNode *node);
void evolve_Euler(tNode *node);


/* trouble.c */
#define NOTROUBLES 11 /* if accum. trouble<=-NOTROUBLES, we switch to dg */
int evolve_set_trouble_score_mesh(tMesh *mesh);
int evolve_read_trouble_score_mesh(tMesh *mesh);
void evolve_prepare_do_over_mesh(tMesh *mesh);
void evolve_switch_troubled_nodes_mesh(tMesh *mesh);
void evolve_switch_nontroubled_nodes_mesh(tMesh *mesh, int notroubles);
void evolve_collect_u_p_data_mesh(tMesh *mesh, pVLList *u_p);
int trouble_reset_evo_troubled_mesh(tMesh *mesh);
