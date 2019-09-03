/* nmesh_output.h */
/* Wolfgang Tichy, February 2019 */

    
/* output.c */
int mesh_output(tMesh *mesh);
int write_mesh(tMesh *mesh, int Iteration, double Time);
int TimeForMeshOutput_vl(tMesh *mesh, tVarList *vl);

/* quick output for debugging */
void write_array(tNode *node, tArray *va, char *name, int as_1d,
                 int fake_it, double fake_t);
void write_var(tNode *node, char *name, int as_1d,
               int fake_it, double fake_t);

/* quick varlist output in vtk format */
void write_vl(tNode *node, tVarList *vl, int as_1d,
              int fake_it, double fake_t);

/* gnuplot.c */
void outputPatchPlanes_meshvar(tMesh *mesh, char *name, int It, double T);
