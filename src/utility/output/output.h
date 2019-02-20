/* output.h */
/* Wolfgang Tichy 4/2005 */

//#include <sys/stat.h>
//#include <sys/types.h>


/* output2d.c */
void output2d_meshvar(tMesh *mesh, char *name, int It, double T);

/* gnuplot2d.c */
void write_plane_ascii(tNode *node, FILE *fp, int normal, int plane[], int iv,
                       int Iter, double Time);
