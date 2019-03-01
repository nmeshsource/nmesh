/* evolve.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for coupled evolution systems */




/* evolve.c */
void evolve_register_subsys_u_rhs_src(tMesh *mesh, tVarList *u,
                                      FuncPointer rhs, FuncPointer src);
int evolve_finalize(tMesh *mesh);
void evolve_setrhs(tNode *node, pVLList *rhs, pVLList *u);
int evolve_myln(tMesh *mesh);
void evolve(tNode *node);

int evolve_test_init(tMesh *mesh);
int evolve_test_analyze(tMesh *mesh);


/* RungeKutta.c */
void evolve_RK4(tNode *node);
