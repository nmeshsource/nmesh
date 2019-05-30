/* nmesh_output.h */
/* Wolfgang Tichy, February 2019 */

    
/* output.c */
int mesh_output(tMesh *mesh);
int write_mesh(tMesh *mesh, int Iteration, double Time);
int TimeForMeshOutput_vl(tMesh *mesh, tVarList *vl);

/* quick_out.c for debugging */
/*
void quick_Vars_output(tBox *box, char *names, double fake_t, int fake_i);
void quick_Array_output(tBox *box, double *Ar, char *name,
                        double fake_t, int fake_i);
void quick_VarList_output(tBox *box, tVarList *vl, double fake_t, int fake_i);
*/

/* gnuplot.c */
void outputPatchPlanes_meshvar(tMesh *mesh, char *name, int It, double T);
