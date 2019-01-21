/* nmesh_main.h */
/* Wolfgang Tichy, 1/2019 */

/* constants */
#ifdef PI
#undef PI
#endif
#define PI  3.14159265358979323846264338327950
#define PIh 1.57079632679489661923132169163975
#define PIq 0.785398163397448309615660845819876

/* Indices */
#define Index(i,j,k)   ((i)+n[0]*((j)+n[1]*(k)))
#define Ind_n(i,j,k,n) ((i)+(n[0])*((j)+(n[1])*(k)))
/* ijk = i + n0*j + n0*n1*k, thus:
   ijk/(n0*n1) = k 
   (ijk - n0*n1*k)/n0 = j
   (ijk - n0*n1*k - n0*j ) = i   */
#define kOfInd_n(ijk,n)        ((ijk)/((n[0])*(n[1])))
#define jOfInd_n_k(ijk,n,k)    (((ijk) - (n[0])*(n[1])*(k))/(n[0]))
#define iOfInd_n_jk(ijk,n,j,k) ((ijk) - (n[0])*(n[1])*(k) - (n[0])*(j))
#define forallijk(i,j,k) \
  for (k = 0; k < n[2]; k++) \
  for (j = 0; j < n[1]; j++) \
  for (i = 0; i < n[0]; i++)    


/* NOTE: In C99 these two have the same effect:
   #pragma omp parallel for
   _Pragma ( "omp parallel for" )
*/
/* To parallelize with OpenMP we need _Pragma ( "omp parallel for" ) 
   in many places. But for different applications we want to switch 
   them on or off depending on where they are.
   SGRID_LEVEL2_Pragma used for omp loops over a plane in a pat (2d)
   SGRID_LEVEL3_Pragma used for omp loops over all points in a pat (3d)
   SGRID_LEVEL4_Pragma used for omp loops over all pates
   SGRID_LEVEL6_Pragma used for 6d omp loops (e.g. loop over pat while interpolating onto each point)
   more can be defined easily.
   */
/* define SGRID_LEVEL2_Pragma macros that allow us to include
   certain pragmas only if certain things like LEVEL2_Pragmas are defined */
#ifdef LEVEL2_Pragmas
#define SGRID_LEVEL2_Pragma(x)  _Pragma ( #x )
#else
#define SGRID_LEVEL2_Pragma(x)
#endif

#ifdef LEVEL3_Pragmas
#define SGRID_LEVEL3_Pragma(x)  _Pragma ( #x )
#else
#define SGRID_LEVEL3_Pragma(x)
#endif

#ifdef LEVEL4_Pragmas
#define SGRID_LEVEL4_Pragma(x)  _Pragma ( #x )
#else
#define SGRID_LEVEL4_Pragma(x)
#endif

#ifdef LEVEL6_Pragmas
#define SGRID_LEVEL6_Pragma(x)  _Pragma ( #x )
#else
#define SGRID_LEVEL6_Pragma(x)
#endif

#ifdef TOPLEVEL_Pragmas
#define SGRID_TOPLEVEL_Pragma(x)  _Pragma ( #x )
#else
#define SGRID_TOPLEVEL_Pragma(x)
#endif

#if defined(LEVEL6_Pragmas) || defined(TOPLEVEL_Pragmas)
#define SGRID_LEVEL6orTOP_Pragma(x)  _Pragma ( #x )
#else
#define SGRID_LEVEL6orTOP_Pragma(x)
#endif


/* snap effect for mesh coordinates */
#define dequaleps 1e-10
#define dless(a,b) ((a)<(b)-dequaleps)
#define dequal(a,b) (!(dless(a,b)||dless(b,a)))
#define dgreater(a,b) ((a)>(b)+dequaleps)

/* approx. <= and >= with dequaleps tolerance */
#define dlesseq(a,b) ( (a)<(b)+dequaleps )
#define dgreatereq(a,b) ( (a)>(b)-dequaleps )

#define signum(v) ((v) > 0.0 ? (1.0) : ((v) < 0.0 ? (-1.0) : (0.0)))


/* skeleton.c */
void AddFun(tMesh *mesh, int step, int (*f)(tMesh *), char *name);
void RunFun(tMesh *mesh, int step);

/* parameters.c */
/* parameter data base structure */
typedef struct {
  char *name;
  char *value;
  char *description;
  double numericalvalue; /* some pars are pure numbers, we cache them here */
  int booleanvalue; /* some pars are true/false, we cache them here as 1/0 */
} tParameter;
void makeparameter(char *name, char *value, char *description);
void AddPar(char *name, char *value, char *description);
void AddOrModifyPar(char *name, char *value, char *description);
void Sets(char *name, char *value);
void Seti(char *name, int i);
void Setd(char *name, double d);
char *Gets(char *name);
char *GetsLax(char *name);
int Geti(char *name);
double Getd(char *name);
int Getv(char *name, char *value);
int GetvLax(char *name, char *value);
char *NextEntry(char *list);
void Appends(char *name, char *value);
char *GetsInd(int i);
char *GetnameInd(int i);
tParameter *GetPointerTo_pbd(void);
int GetnParameters(void);
double GetCachedNumValByParIndex(int i);
int GetCachedBoolValByParIndex(int i);
int GetParIndex(char *name);
int Set_pdb_iStart_AtPar(char *name);
void print_pdb_i1_i2(tParameter *pdb, int i1, int i2, int pr_ind, int pr_cache);
void print_parameter_database(void);
void create_copy_of_pdb1_in_pdb2(tParameter *pdb1, int npdb1, int npdb1max,
                                 tParameter **pdb2);
tParameter *make_empty_pdb(int npdb1max);
void copy_pdb(tParameter *pdb1, int npdb1, tParameter *pdb2);
void free_pdb(tParameter *pdb1, int npdb1);
void free_global_parameter_database_contents(void);


/* tensors.c */
#define NINDEXLIST 100
void tensorindexlist(char *tensorindices, int *nilist, char **ilist, int *sym);

/* variables.c */
/* variable data base structures */
typedef struct {
  double iotime[4];
  int ioiter[4];
  int ioflag[4];
  int ioatall;
} tIO;
typedef struct tVAR {
  char *name;
  char *tensorindices;
  char *description;
  int index;
  int ncomponents;
  int component;
  tIO *io;
  double farlimit;
  double falloff;
  double propspeed;
  int sym[3];
  int constant;
} tVar;
int IndLax(char *name);
int Ind(char *name);
int Set_vdb_iStart_AtPar(char *name);
void AddVar(char *name, char *indices, char *description);
void AddConstantVar(char *name, char *tensorindices, char *description);
void AddVarToGrid(tMesh *mesh, char *name, char *tensorindices,
                  char *description);
tVarList *AddDuplicate(tVarList *vl, char *postfix);
tVarList *AddDuplicateEnable(tVarList *vl, char *postfix);

char *VarName(int i);
int VarNComponents(int i);
int VarComponent(int i);
int IndComponent0(int i);
char *VarNameComponent0(char *name);
char *VarTensorIndices(int i);
void VarNameSetBoundaryInfo(char *name, 
			    double farlimit, double falloff, double propspeed);
double VarFallOff(int i);
double VarFarLimit(int i);
double VarPropSpeed(int i);
int VarSymmetry(int i, int dir);
void VarNameSetConstantFlag(char *name);
int VarConstantFlag(int i);

void prvarlist(tVarList *v);
void prvarlist_inpat(tPat *pat, tVarList *v);
tVarList *vlalloc(tMesh *mesh);
void vlenable(tVarList *v);
void vlenablemesh(tMesh *mesh, tVarList *v);
void vldisable(tVarList *v);
void vlfree(tVarList *u);
void vlpushone(tVarList *v, int vi);
void vlpush(tVarList *v, int vi);
void vlpushvl(tVarList *v, tVarList *u);
void vldropone(tVarList *v, int vi);
void vldrop(tVarList *v, int vi);
void vldropn(tVarList *v, int n);
tVarList *vlduplicate(tVarList *v);
void vlsetconstant(tVarList *u, const double c);
void vlcopy(tVarList *v, tVarList *u);
void vlcopymesh(tMesh *mesh, tVarList *v, tVarList *u);
void varcopy(tMesh *mesh, int iv, int iu);
void vlswap(tVarList *v, tVarList *u);
void varswap(tMesh *mesh, int iv, int iu);
void vlaverage(tVarList *r, tVarList *a, tVarList *b);
void vlsubtract(tVarList *r, tVarList *a, tVarList *b);
void vladd(tVarList *r, double ca, tVarList *a, double cb, tVarList*b);
void varadd(tMesh *mesh, int ir, double ca, int ia, double cb, int ib);
void vladdto(tVarList *r, const double ca, tVarList *a);

tVarList *VLPtrEnable1(tMesh *mesh, char *varname);
void VLDisableFree(tVarList *vl);

/* utilities.c */
void  errorexit(char *file, int line, const char *func, char *s);
void errorexits(char *file, int line, const char *func, char *s, char *t);
void errorexiti(char *file, int line, const char *func, char *s, int i);
#define errorexit(s)     errorexit(__FILE__, __LINE__, __func__, (s))
#define errorexits(s,t) errorexits(__FILE__, __LINE__, __func__, (s), (t))
#define errorexiti(s,i) errorexiti(__FILE__, __LINE__, __func__, (s), (i))
#define PRF     printf("%s", __func__)
#define PRFs(s) printf("%s%s", __func__, s)

void Yo(double x);
void prdivider(int n);
void initTimeIn_s(void);
double getTimeIn_s(void);
void prTimeIn_s(char *comment);
double min2(double x, double y);
double min3(double x, double y, double z);
double max2(double x, double y);
double max3(double x, double y, double z);
double min_in_1d_array(double *f, int n, int *imin);
double max_in_1d_array(double *f, int n, int *imax);
double min2_in_1d_array(double *f0, int n0, double *f1, int n1, 
                        int *ai, int *imin);
double max2_in_1d_array(double *f0, int n0, double *f1, int n1, 
                        int *ai, int *imax);
double min3_in_1d_array(double *f0, int n0, double *f1, int n1, double *f2, int n2,
                        int *ai, int *imin);
double max3_in_1d_array(double *f0, int n0, double *f1, int n1, double *f2, int n2,
                        int *ai, int *imax);
int copy_file_into_dir(char *fname, char *dir);
int system2(char *s1, char *s2);
int system3(char *s1, char *s2, char *s3);
int system_emu(const char *command);
int construct_argv(char *str, char ***argv);
double *dmalloc(int n);
int *imalloc(int n);
char *cmalloc(int n);
void *pmalloc(int n);

/* endianIO.c */
size_t fwrite_double_little(const double *buf, size_t nmemb, FILE *fp);
size_t fread_double_little(double *buf, size_t nmemb, FILE *fp);

/* nmesh_MPI.c */
int nmesh_MPI_Init(int *pargc, char ***pargv);
int nmesh_MPI_Finalize(void);
int nmesh_MPI_rank(void);
int nmesh_MPI_size(void);  
int nmesh_MPI_barrier(void);
