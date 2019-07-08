/* nmesh_evolve.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for global functions */



/* evolve.c */
void evolve_register_subsys(tMesh *mesh, tVarList *u,
                FuncPointer prelim, FuncPointer limdata, FuncPointer limiter,
                FuncPointer presurf, FuncPointer setsrc, FuncPointer volrhs,
                FuncPointer surfrhs);
void evolve_register_subsys_u_rhs_lim(tMesh *mesh, tVarList *u,
                                      FuncPointer volrhs, FuncPointer surfrhs,
                                      FuncPointer limdata,
                                      FuncPointer limiter);
void evolve_print_evosys(tMesh *mesh);
int evolve_init_communication_structs(tMesh *mesh);
int evolve_free_communication_structs(tMesh *mesh);
