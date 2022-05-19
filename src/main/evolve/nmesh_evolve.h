/* nmesh_evolve.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for global functions */


/* macro to set a func for a varlist in an evolution bin */
#define evolve_SetFun(bin, f, vl) evolve_SetEvoFun(bin, f, vl, MSTR(f))


/* evolve.c */
int evolve_myln(tMesh *mesh);

/* evosys.c */
void evolve_register_vl(tVarList *vl);
void evolve_SetEvoFun(int bin, FuncPointer f, tVarList *vl, const char *name);
void evolve_register_subsys(tMesh *mesh, tVarList *u,
                FuncPointer prelim, FuncPointer limdata, FuncPointer limiter,
                FuncPointer presurf, FuncPointer setsrc, FuncPointer volrhs,
                FuncPointer surfrhs);
void evolve_register_subsys_u_rhs_lim(tMesh *mesh, tVarList *u,
                                      FuncPointer volrhs, FuncPointer surfrhs,
                                      FuncPointer limdata,
                                      FuncPointer limiter);
tVarList *evolve_get_rhs_vl(tVarList *vl);
void evolve_print_evosys(tMesh *mesh);
int var_added_by_evolve_init_evosys(tMesh *mesh, int vi);
int evolve_init_communication_structs(tMesh *mesh);
int evolve_free_communication_structs(tMesh *mesh);
void init_all_myln_myindc_in_evosys(tMesh *mesh);
void free_all_myln_myindc_in_evosys(tMesh *mesh);

/* trouble.c */
int evolve_RDMP_trouble(tNode *node, tVarList *vlu, tVarList *vlu_p);
