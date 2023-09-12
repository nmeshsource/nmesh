/* basis.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for basis local functions */


/* basis.c */
int basis_init_globals(tMesh *mesh);

/* gridpoints.c */
int init_gridpoints(tMesh *mesh);
int free_gridpoints(tMesh *mesh);

/* Fourier.c */
void Fourier_coeffs(int N, const double u[], double c[]);
void set_TrafoArray(tArray *At,
                    void (*get_coeffs)(int N, const double *u, double *c));
