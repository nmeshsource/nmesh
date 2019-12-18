/* advectionFDy.h */
/* (c) Wolfgang Tichy 12/2019 */
/* header file for advectionFDy local functions */


/* structure that holds global vars and pars */
typedef struct {
  int sin_profile;
  int square_profile;
} tadvectionFDy;


/* advectionFDy.c */
int advectionFDy_init_global_pars(tMesh *mesh);
void advectionFDy_set_profile_pt(double xyz[3], double t, int nv, double *u);
int advectionFDy_init(tMesh *mesh);
int advectionFDy_analyze(tMesh *mesh);
