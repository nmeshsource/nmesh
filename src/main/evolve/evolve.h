/* evolve.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for coupled evolution systems */



/* evolve.c */
int evolve_myln(tMesh *mesh);
void evolve_setrhs_mesh(tMesh *mesh, pVLList *rhs, pVLList *u);
void evolve_limiter_mesh(tMesh *mesh, pVLList *u);
void evolve_setrhs(tNode *node, pVLList *rhs, pVLList *u, int request_surfs);
int evolve_filter_evosys_mesh(tMesh *mesh);


/* evosys.c */
int evolve_free_evosys(tMesh *mesh);
void evolve_print_evosys(tMesh *mesh);
int evolve_init_evosys(tMesh *mesh);
void evolve_request_surfaces(tNode *node, pVLList *u);


/* evolve_test.c */
int evolve_test_init(tMesh *mesh);
int evolve_test_analyze(tMesh *mesh);


/* RungeKutta.c */
void evolve_RK4(tNode *node);
void evolve_Euler(tNode *node);
void evolve_RK4_mesh(tMesh *mesh);
void evolve_Euler_mesh(tMesh *mesh);
void evolve_sspRK3_mesh(tMesh *mesh);
