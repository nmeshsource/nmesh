/* advection1.h */
/* (c) Wolfgang Tichy 1/2019 */
/* header file for advection1 local functions */


/* structure that holds global vars and pars */
typedef struct {
  /* saved indices of some vars */
  int ifx;   /* Ind("advection1_fx"); */
  int idivf; /* Ind("advection1_divf"); */
  void (*numflux)(tDGinfo *d);  /* func pointer for numerical flux */
  double direction[3];
  int outerBC_influxes; /* Par("advection1_outerBC_influxes") */
  int sin_profile;
  int square_profile;
  int fd_dissfac;       /* Par("advection1_fd_dissfac") */
  int fd_dissorder;     /* Par("advection1_fd_dissorder") */
} tadvection1;


/* advection1.c */
int advection1_init_global_pars(tMesh *mesh);
void advection1_fluxes_pt(tDGinfo *d);
int advection1_init(tMesh *mesh);
int advection1_analyze(tMesh *mesh);
int advection1_refine(tMesh *mesh);
