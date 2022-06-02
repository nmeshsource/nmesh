/* nmesh_output.h */
/* Wolfgang Tichy, February 2019 */



/* global output pars */
typedef struct tOUTPUT {
  int Noutpt;               /* 0doutput at this number of specific points */
#define Noutptmax 2         /* max Noutpt */
  double xpt[Noutptmax][3]; /* x,y,z coords of specific 0doutput points */
} tOutput;

    
/* output.c */
int mesh_output(tMesh *mesh);
int write_mesh(tMesh *mesh, int Iteration, double Time);
int TimeForMeshOutput_vl(tMesh *mesh, tVarList *vl);

/* VTK_out.c: for quick output for debugging */
void write_array(tNode *node, tArray *va, char *name, int as_1d,
                 int fake_it, double fake_t);
void write_var(tNode *node, char *name, int as_1d,
               int fake_it, double fake_t);
void write_vl(tNode *node, tVarList *vl, int as_1d,
              int fake_it, double fake_t);
void write_var_nodenamelist(tMesh *mesh, char *nodenamelist, char *varname,
                            int as_1d, int fake_it, double fake_t);
void write_vl_nodenamelist(char *nodenamelist, tVarList *vl,
                           int as_1d, int fake_it, double fake_t);

/* gnuplot.c */
void outputPatchPlanes_meshvar(tMesh *mesh, char *name, int It, double T);
