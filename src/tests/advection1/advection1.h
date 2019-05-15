/* advection1.h */
/* (c) Wolfgang Tichy 1/2019 */
/* header file for advection1 local functions */



/* advection1.c */
int advection1_init_global_pars(tMesh *mesh);
void advection1_fluxes_pt(tDGinfo *d);
int advection1_init(tMesh *mesh);
int advection1_analyze(tMesh *mesh);
