/* nmesh_evolve.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for global functions */



/* evolve.c */
void evolve_register_subsys(tMesh *mesh, tVarList *u,
                FuncPointer prelim, FuncPointer limdata, FuncPointer limiter,
                FuncPointer presurf, FuncPointer setsrc, FuncPointer setrhs);
void evolve_register_subsys_u_rhs_src_lim(tMesh *mesh, tVarList *u,
                                          FuncPointer rhs, FuncPointer src,
                                          FuncPointer limdata,
                                          FuncPointer limiter);
void evolve_print_evosys(tMesh *mesh);
