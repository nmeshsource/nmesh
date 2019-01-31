/* nmesh_main.h */
/* Wolfgang Tichy, 1/2019 */

/* use MSTR(x) to transform x into a string */
#define XMSTR(x) STRING(x)
#define MSTR(x ) #x

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
void AddMeshFun(tMesh *mesh, int step, int (*f)(tMesh *), char *name);
void RunMeshFun(tMesh *mesh, int step);
/* conveniece macros for functions */
#define AddFun(step, f)  AddMeshFun(mesh, step, f, MSTR(f))
#define RunFun(step)     RunMeshFun(mesh, step)


/* parameters.c */
/* parameter data base structure */
typedef struct tPAR {
  char *name;
  char *value;
  char *description;
  double numericalvalue; /* some pars are pure numbers, we cache them here */
  int booleanvalue; /* some pars are true/false, we cache them here as 1/0 */
} tParameter;
void AddMeshPar(tMesh *mesh, char *name, char *value, char *description);
void AddOrModifyMeshPar(tMesh *mesh, char *name, char *value, char *description);
void makeparameter(tMesh *mesh, char *name, char *value, char *description);
int findparameterindex(tMesh *mesh, char *name, int fatal);
void MeshParSets(tMesh *mesh, int pi, char *value);
void MeshParSeti(tMesh *mesh, int pi, int i);
void MeshParSetd(tMesh *mesh, int pi, double d);
void MeshParAppends(tMesh *mesh, int pi, char *value);
char *MeshParGets(tMesh *mesh, int i);
char *MeshParGetsLax(tMesh *mesh, int i);
int MeshParGeti(tMesh *mesh, int i);
double MeshParGetd(tMesh *mesh, int i);
int MeshParGetb(tMesh *mesh, int i);
int MeshParGetv_fatal(tMesh *mesh, int i, char *value, int fatal);
/* conveniece macros to query and set  pars */
#define AddPar(name, val, desc) AddMeshPar(mesh, (name), (val), (desc))
#define AddOrModifyPar(name, val, desc) AddOrModifyMeshPar(mesh, (name), \
                                                           (val), (desc))
#define GetMeshParIndex(mesh, name, fatal) findparameterindex((mesh), (name), (fatal))
#define Par(name)     findparameterindex(mesh, (name), 0)
#define ParF(name)    findparameterindex(mesh, (name), 1)
#define Geti(ip)      MeshParGeti(mesh, (ip))
#define Getd(ip)      MeshParGetd(mesh, (ip))
#define Getb(ip)      MeshParGetb(mesh, (ip))
#define Gets(ip)      MeshParGets(mesh, (ip))
#define GetsLax(ip)   MeshParGetsLax(mesh, (ip))
#define Getv(ip, val) MeshParGetv_fatal(mesh, (ip), (val), 1)
#define GetvLax(ip, val) MeshParGetv_fatal(mesh, (ip), (val), 0)
// use like this: if(Getv(Par("parname"), "value")) x = Setd(Par("parname2"))
#define Seti(ip, i)   MeshParSeti(mesh, (ip), (i))
#define Setd(ip, x)   MeshParSetd(mesh, (ip), (x))
#define Sets(ip, s)   MeshParSets(mesh, (ip), (s))
#define Appends(ip, s)  MeshParAppends(mesh, (ip), (s))


/* tensors.c */
#define NINDEXLIST 100
void tensorindexlist(char *tensorindices, int *nilist, char **ilist, int *sym);

/* variables.c */
/* variable data base structure */
typedef struct tVAR {
  char *name;
  char *tensorindices;
  char *description;
  int index;
  int ncomponents;
  int component;
  double farlimit;
  double falloff;
  int sym[3];
  int constant;
} tVar;
/* functions to create and access variables */
void AddMeshVar(tMesh *mesh, char *name, char *tensorindices, char *description);
void AddConstantMeshVar(tMesh *mesh, char *name, 
                        char *tensorindices, char *description);
int MeshVarIndLax(tMesh *mesh, char *name) ;
int MeshVarInd(tMesh *mesh, char *name);
int Set_vdb_iStart_AtVar(tMesh *mesh, char *name);

char *MeshVarName(tMesh *mesh, int i);
int MeshVarNComponents(tMesh *mesh, int i);
int MeshVarComponent(tMesh *mesh, int i);
int MeshVarIndComponent0(tMesh *mesh, int i);
char *MeshVarNameComponent0(tMesh *mesh, char *name);
char *MeshVarTensorIndices(tMesh *mesh, int i);
void MeshVarNameSetBoundaryInfo(tMesh *mesh, char *name,
			        double farlimit, double falloff);
void MeshVarNameSetConstantFlag(tMesh *mesh, char *name);
double MeshVarFallOff(tMesh *mesh, int i);
double MeshVarFarLimit(tMesh *mesh, int i);
int MeshVarSymmetry(tMesh *mesh, int i, int dir);
int MeshVarConstantFlag(tMesh *mesh, int i);
/* conveniece macros for vars */
#define VarName(i) MeshVarName(mesh, (i))
#define Ind(name)  MeshVarInd(mesh, (name))
#define AddVar(name, tensorindices, description) \
  AddMeshVar(mesh, (name), (tensorindices), (description))

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
