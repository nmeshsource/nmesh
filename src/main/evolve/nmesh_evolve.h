/* nmesh_evolve.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for global functions */



/* evolve.c */
void evolve_register_subsys_u_rhs_src(tMesh *mesh, tVarList *u,
                                      void (*rhs)(), void (*src)());
void evolve_print_evosys(tMesh *mesh);
