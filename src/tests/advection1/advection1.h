/* advection1.h */
/* (c) Wolfgang Tichy 1/2019 */
/* header file for advection1 local functions */


/* structure that holds global vars and pars */
typedef struct {
  void (*numflux)(tDGinfo *d);  /* func pointer for numerical flux */
  double direction[3];
  int outerBC_influxes;
  int sin_profile;
  int square_profile;
} tadvection1;


/* advection1.c */
int advection1_init_global_pars(tMesh *mesh);
void advection1_fluxes_pt(tDGinfo *d);
int advection1_init(tMesh *mesh);
int advection1_analyze(tMesh *mesh);
int advection1_refine(tMesh *mesh);
