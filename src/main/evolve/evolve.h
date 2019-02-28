/* evolve.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for coupled evolution systems */




/* evolve.c */
int evolve_myln(tMesh *mesh);
void evolve(tNode *node);
int evolve_finalize(tMesh *mesh);

int evolve_test(tMesh *mesh);


/* RungeKutta.c */
void evolve_RK4(tNode *node);
