/* output.h */
/* Wolfgang Tichy 4/2005 */


/* object to hold pars we need to pass around */
typedef struct tOUTPARS {
  char *name;
  int p;
  char *nodeloc;
  int text;
  int arrange_as_1d;
  int flt;
  int dbl;
} tOutpars;


/* output.c */
int TimeForMeshOutput_di_dt(tMesh *mesh, int di, double dt) ;

/* output2d.c */
void output2d_vl(tVarList *vl, int It, double T);

/* output1d.c */
void output1d_vl(tVarList *vl, int It, double T);

/* output3d.c */
void output3d_vl(tVarList *vl, int It, double T);

/* gnuplot2d.c */
void write_plane_ascii(tNode *node, FILE *fp, int normal, int plane[],
                       tArray *va, int Iter, double Time);
void write_line_ascii(tNode *node, FILE *fp, int dir, int axis[],
                      tArray *va, int Iter, double Time);
void gnuplot_output2d_meshvar(tMesh *mesh, char *name, int It, double T);
void gnuplot_output1d_meshvar(tMesh *mesh, char *name, int It, double T);

/* VTK_out.c */
FILE *fopen_vtk(char *varname, char *outdir, char *suffix,
                int p, char *nstr, int series);
void write3d_vtk(tNode *node, FILE *fp, tArray *va, int Iter, double Time,
                 int series, tOutpars *par);
void vtk_output3d_meshvar(tMesh *mesh, char *name, int It, double T);
void vtk_output2d_meshvar(tMesh *mesh, char *name, int It, double T);

/* XDMF_out.c */
void output2d_xdmf(tVarList *vl, int It, double T);
void write_plane_xdmf(tVarList *vl, int norm, char *outdir, double Time);
void output3d_xdmf(tVarList *vl, int It, double Time);
size_t write_buffer(const double *buf, int buflen, int dbl, FILE *fp);
size_t write_3buffers(const double *b1, const double *b2, const double *b3,
                      int buflen, int dbl, FILE *fp);
size_t write_buffer_idx(const double *buf, intList *idx, int dbl, FILE *fp);
size_t write_3buffers_idx(const double *b1, const double *b2, const double *b3,
                          intList *idx, int dbl, FILE *fp);
size_t fwrite_buffer_idx(const void *ptr, size_t size, 
                         intList *idx, FILE *fp);
size_t fwrite_3buffers_idx(const void *p1, const void *p2, const void *p3,
                           size_t size, intList *idx, FILE *fp);

/* output0d.c */
void output0d_vl(tVarList *vl, int It, double T);
void output0d_mesh_vl(tVarList *vl, tPat *pat, int It, double T);
void output0d_filename(tMesh *mesh, char *filename, int len,
                       char *name, char *type, tPat *pat);
void output0d_value(char *filename, double time, double val);
