/* main.h */
/* Wolfgang Tichy, 1/2019 */

int iterate_parameters(tMesh *mesh, int next);

int read_command_line(tMesh *mesh, int argc, char **argv);
int parse_command_line_options(tMesh *mesh);
int make_output_directory(tMesh *mesh);
int inidata_mesh(tMesh *mesh);
int evolve_mesh(tMesh *mesh);
int finalize_mesh(tMesh *mesh);
