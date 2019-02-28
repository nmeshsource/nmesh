/* evolve.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for coupled evolution systems */




/* evolve.c */
void evolve(tNode *node);

int evolve_test(tMesh *mesh);


/* RungeKutta.c */
void evolve_RK4(tNode *node);
