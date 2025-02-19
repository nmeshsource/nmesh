/* nmesh_evolve.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for global functions */


/* macro to set a func for a varlist in an evolution bin */
#define evolve_SetFun(bin, f, vlu) evolve_SetEvoFun(bin, f, vlu, MSTR(f))


/* evolve.c */
int evolve_myln(tMesh *mesh);
void evolve_limiter_mesh(tMesh *mesh, pVLList *u, int opt);

/* evosys.c */
void evolve_register_vl(tVarList *vl);
void evolve_SetEvoFun(int bin, EvoFuncPtr f, tVarList *vlu, const char *name);
void evolve_register_subsys(tMesh *mesh, tVarList *u,
                EvoFuncPtr prelim, EvoFuncPtr limdata, EvoFuncPtr limiter,
                EvoFuncPtr presurf, EvoFuncPtr setsrc, EvoFuncPtr volrhs,
                EvoFuncPtr surfrhs);
void evolve_register_subsys_u_rhs_lim(tMesh *mesh, tVarList *u,
                                      EvoFuncPtr volrhs, EvoFuncPtr surfrhs,
                                      EvoFuncPtr limdata,
                                      EvoFuncPtr limiter);
void evolve_SetVLx(tVarList *vlx, tVarList *vlu);
tVarList *evolve_get_rhs_vl(tVarList *vl);
void evolve_print_evosys(tMesh *mesh);
int var_added_by_evolve_init_evosys(tMesh *mesh, int vi);
int evolve_init_communication_structs(tMesh *mesh);
int evolve_free_communication_structs(tMesh *mesh);
void init_all_myln_myindc_in_evosys(tMesh *mesh);
void free_all_myln_myindc_in_evosys(tMesh *mesh);

/* trouble.c */
int evolve_RDMP_trouble(tNode *node, tVarList *vlu, tVarList *vlu_p,
                        double deltafac, double delta0, double epsilon);
double evolve_Persson_indicator_ncoeffs(tNode *node, int iu, double u_scale,
                                        int ncoeffs[3]);
int evolve_Persson_trouble_ncoeffs(tNode *node, int iu, double u_scale,
                                   int ncoeffs[3],
                                   double alpha, double alpha_fv);
int evolve_Persson_trouble_ncoeffs_dg(tNode *node, int iu, double u_scale,
                                      int ncoeffs[3],
                                      double alpha, double alpha_fv);
int evolve_Persson_trouble(tNode *node, int iu, double u_scale,
                           double alpha, double alpha_fv);
int trouble_score(tNode *node, int troubled);
void trouble_print_if_name(tNode *node, const char *nname,
                           int trbl, const char *text);
int evolve_evosteps_troubled(tNode *node);
