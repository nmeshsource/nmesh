/* output.h */
/* Wolfgang Tichy 4/2005 */


/* object to hold pars we need to pass around */
typedef struct tOUTPARS {
  char *name;
  int p;
  char *nodeloc;
  int text;
  int binary;
  int arrange_as_1d;
  int flt;
  int dbl;
} tOutpars;


/* output.c */
int TimeForMeshOutput_di_dt(tMesh *mesh, int di, double dt) ;

/* output2d.c */
void output2d_meshvar(tMesh *mesh, char *name, int It, double T, int wr_xyz);

/* output1d.c */
void output1d_meshvar(tMesh *mesh, char *name, int It, double T, int wr_xyz);

/* gnuplot2d.c */
void write_plane_ascii(tNode *node, FILE *fp, int normal, int plane[],
                       tArray *va, int Iter, double Time);
void write_line_ascii(tNode *node, FILE *fp, int dir, int axis[],
                      tArray *va, int Iter, double Time);

/* output3d.c */
FILE *fopen_vtk(char *varname, char *outdir, char *suffix,
                int p, char *nstr, int series);
void write3d_vtk(tNode *node, FILE *fp, tArray *va, int Iter, double Time,
                 int series, tOutpars *par);
void output3d_meshvar(tMesh *mesh, char *name, int It, double T, int wr_xyz);

/* XDMF_out.c */
void write_plane_xdmf(tNode *node, int normal, int plane[], tArray *va,
                      char *vname, char *outdir, char *suffix, double Time,
                      int write_xyz);
size_t write_buffer_idx(const double *buf, intList *idx, int dbl, FILE *fp);
size_t write_3buffers_idx(const double *b1, const double *b2, const double *b3,
                          intList *idx, int dbl, FILE *fp);
size_t fwrite_buffer_idx(const void *ptr, size_t size, 
                         intList *idx, FILE *fp);
size_t fwrite_3buffers_idx(const void *p1, const void *p2, const void *p3,
                           size_t size, intList *idx, FILE *fp);
