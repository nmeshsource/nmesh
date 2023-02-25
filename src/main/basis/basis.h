/* basis.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for basis local functions */


/* structure that holds global coordinates vars and pars */
typedef struct {
  int expfilter_JacobianPower; /* Par("basis_expfilter_JacobianPower") */
  int filter_fv;               /* Par("basis_filter_fv") */
} tbasis;

/* basis.c */
int basis_init_globals(tMesh *mesh);

/* gridpoints.c */
int init_gridpoints(tMesh *mesh);
int free_gridpoints(tMesh *mesh);

/* Fourier.c */
void Fourier_coeffs(int N, const double u[], double c[]);
