/* center.h */
/* Wolfgang Tichy, Jan 2024 */

/* center.c */
int center_update(tMesh *mesh);
int center_init_globals(tMesh *mesh);
double average_grid_spacing(tMesh *mesh, double x[3]);
int center_track_extremum(tMesh *mesh, double h, int var, int findMax,
                          const double xold[3], double minmove,
                          double xnew[3]);
