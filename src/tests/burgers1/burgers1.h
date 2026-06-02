/* burgers1.h */
/* (c) Wolfgang Tichy 1/2019 */
/* header file for burgers1 local functions */


/* structure that holds global scalarwave1 vars and pars */
typedef struct {
  /* saved indices of some vars */
  int idivf; /* Ind("burgers1_divf") */
  void (*numflux)(tDGinfo *d); // func pointer for numerical flux
  /* frequently used pars */
  double direction[3];
  int profile; /* Par("burgers1_profile") */
} tburgers1;


/* burgers1.c */
int burgers1_init_global_pars(tMesh *mesh);
int burgers1_init(tMesh *mesh);
void burgers1_set_profile_pt(tMesh *mesh,
                              double xyz[3], double t, double *u);
int burgers1_analyze(tMesh *mesh);
void burgers1_rec_u_f_lam(tFVinfo *fv, tDGinfo *d);
