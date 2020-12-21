/* scalarwave1.h */
/* (c) Wolfgang Tichy 9/2019 */
/* header file for scalarwave1 local functions */


/* structure that holds global scalarwave1 vars and pars */
typedef struct {
  void (*numflux)(tDGinfo *d); // func pointer for numerical flux
  /* frequently used pars */
  double k[3];
  int sin_profile;
  int square_profile;
} tscalarwave1;


/* scalarwave1.c */
int scalarwave1_init_global_pars(tMesh *mesh);
void scalarwave1_fluxes_pt(tDGinfo *d);
void scalarwave1_set_profile_pt(double xyz[3], double t, int nv, double *u);
int scalarwave1_init(tMesh *mesh);
int scalarwave1_analyze(tMesh *mesh);
void scalarwave1_divf_FV(tNode *node, tVarList *vlu);
int scalarwave1_set_use_fv_flag(tMesh *mesh);
void scalarwave1_rec_fluxes_lam(tFVinfo *fv, tDGinfo *d);
