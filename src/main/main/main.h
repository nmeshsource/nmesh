/* main.h */
/* Wolfgang Tichy, 1/2019 */

/* main.c */
int read_command_line(tMesh *mesh, int argc, char **argv);
int parse_command_line_options(tMesh *mesh);
int make_output_directory(tMesh *mesh);
int redirect_stdout_and_stderr(tMesh *mesh, const char *mode);
int check_compiledphysics(tMesh *mesh);
int inidata_mesh(tMesh *mesh);
int evolve_mesh(tMesh *mesh);
int finalize_mesh(tMesh *mesh);
void finalize_all(tMesh *mesh);

/* parameters.c */
void parse_parameter_file(tMesh *mesh, char *parfile);
int iterate_parameters(tMesh *mesh, int next);
int nmesh_update_parameters(tMesh *mesh);
int CheckForBannedPars(tMesh *mesh);

/* sysmon.c */
int sysmon(tMesh *mesh);

/* utilities.c */
int enable_FPEEXCEPTIONS(void);
int file_sync_mode(FILE *fp, int mode);
int fclose_sync_mode(FILE *fp, int mode);
int fclose_buf_sync_mode(FILE *fp, char **buf, int mode);
