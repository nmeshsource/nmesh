/* advectionGhoFDy.h */
/* (c) Wolfgang Tichy 12/2019 */
/* header file for advectionGhoFDy local functions */


/* structure that holds global vars and pars */
typedef struct {
  int sin_profile;
  int square_profile;
} tadvectionGhoFDy;


/* advectionGhoFDy.c */
int advectionGhoFDy_init_global_pars(tMesh *mesh);
void advectionGhoFDy_set_profile_pt(double xyz[3], double t, int nv, double *u);
int advectionGhoFDy_init(tMesh *mesh);
int advectionGhoFDy_analyze(tMesh *mesh);
