/* burgers1.h */
/* (c) Wolfgang Tichy 1/2019 */
/* header file for burgers1 local functions */


/* structure that holds global scalarwave1 vars and pars */
typedef struct {
  /* saved indices of some vars */
  //int if_pix ; /* Ind("scalarwave1_f_pix") */
  void (*numflux)(tDGinfo *d); // func pointer for numerical flux
  /* frequently used pars */
  double direction[3];
} tburgers1;


/* burgers1.c */
int burgers1_init_global_pars(tMesh *mesh);
int burgers1_init(tMesh *mesh);
int burgers1_analyze(tMesh *mesh);
