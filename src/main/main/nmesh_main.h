/* nmesh_main.h */
/* Wolfgang Tichy, 1/2019 */

/* use MSTR(x) to transform x into a string, MSTR_OFVAL(x) first evaluates x
   and then makes the result into a string */
#define MSTR(x) #x
#define MSTR_OFVAL(x) MSTR(x)

/* constants */
#ifdef PI
#undef PI
#endif
#define PI  3.14159265358979323846264338327950
#define PIh 1.57079632679489661923132169163975
#define PIq 0.785398163397448309615660845819876

/* min x we pass into log(x) */
#define LOGARGFLOOR 1e-50

/* define NORET as _Noreturn, unless we use -DNO_C11 */
#ifdef NO_C11
#define NORET
#else
#define NORET _Noreturn
#endif

/* snap effect for mesh coordinates */
#define dequaleps 1e-10
#define dless(a,b) ((a)<(b)-dequaleps)
#define dequal(a,b) (!(dless(a,b)||dless(b,a)))
#define dgreater(a,b) ((a)>(b)+dequaleps)
#define dless_tol(a,b, tol) ((a)<(b)-(tol))
#define dequal_tol(a,b, tol) (!(dless_tol(a,b,tol)||dless_tol(b,a,tol)))

/* approx. <= and >= with dequaleps tolerance */
#define dlesseq(a,b) ( (a)<(b)+dequaleps )
#define dgreatereq(a,b) ( (a)>(b)-dequaleps )

#define signum(v) ((v) > 0.0 ? (1.0) : ((v) < 0.0 ? (-1.0) : (0.0)))


/* skeleton.c */
void AddMeshFun(tMesh *mesh, int step, int (*f)(tMesh *), const char *name);
void remove_all_MeshFuns(tMesh *mesh);
void RunMeshFun(tMesh *mesh, int step);
void PrintMeshFuncs(tMesh *mesh);
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
  int valuelen; /* cache strlen(value) */
} tPar;
int nmesh_load_parameters(tMesh *mesh, char *fname, int fatal, int pr);
void AddMeshPar(tMesh *mesh, const char *name, const char *value,
                const char *description);
void AddOrModifyMeshPar(tMesh *mesh, const char *name, const char *value,
                        const char *description);
void makeparameter(tMesh *mesh, const char *name, const char *value,
                   const char *description);
void free_mesh_pdb_contents(tMesh *mesh);
int findparameterindex(tMesh *mesh, const char *name, int fatal);
void MeshParSets(tMesh *mesh, int pi, const char *value);
void MeshParSeti(tMesh *mesh, int pi, int i);
void MeshParSetd(tMesh *mesh, int pi, double d);
void MeshParAppends(tMesh *mesh, int pi, const char *value);
char *MeshParGets(tMesh *mesh, int i);
char *MeshParGetsLax(tMesh *mesh, int i);
int MeshParGeti(tMesh *mesh, int i);
double MeshParGetd(tMesh *mesh, int i);
int MeshParGetb(tMesh *mesh, int i);
int MeshParGetLen(tMesh *mesh, int i);
int MeshParGetv_fatal(tMesh *mesh, int i, const char *value, int fatal);
/* conveniece macros to query and set  pars */
#define AddPar(name, val, desc) AddMeshPar(mesh, (name), (val), (desc))
#define AddOrModifyPar(name, val, desc) AddOrModifyMeshPar(mesh, (name), \
                                                           (val), (desc))
#define GetMeshParIndex(mesh, name, fatal) findparameterindex((mesh), (name), (fatal))
#define Par(name)     findparameterindex(mesh, (name), 1)
#define ParLax(name)  findparameterindex(mesh, (name), 0)
#define Geti(ip)      MeshParGeti(mesh, (ip))
#define Getd(ip)      MeshParGetd(mesh, (ip))
#define Getb(ip)      MeshParGetb(mesh, (ip))
#define Gets(ip)      MeshParGets(mesh, (ip))
#define GetLen(ip)    MeshParGetLen(mesh, (ip))
#define GetsLax(ip)   MeshParGetsLax(mesh, (ip))
#define Getv(ip, val) MeshParGetv_fatal(mesh, (ip), (val), 1)
#define GetvLax(ip, val) MeshParGetv_fatal(mesh, (ip), (val), 0)
// use like this: if(Getv(Par("parname"), "value")) x = Setd(Par("parname2"))
#define Seti(ip, i)   MeshParSeti(mesh, (ip), (i))
#define Setd(ip, x)   MeshParSetd(mesh, (ip), (x))
#define Sets(ip, s)   MeshParSets(mesh, (ip), (s))
#define Appends(ip, s)  MeshParAppends(mesh, (ip), (s))
/* NOTE: If a parvalue contains several words (e.g. 3), we can do this:
num = sscanf(Gets(Par("parname")), "%s %s %s", str1, str2, str3);
*/

/* banparameters.c */
void BanPar(const char *name, const char *ban_reason);
int ExitIfParBanned(const char *name);
int FreeBannedParList(tMesh *mesh);

/* tensors.c */
#define NINDEXLIST 100
void tensorindexlist(const char *tensind, int *nilist, char **ilist, int *sym);

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
  int sym[3];       /* var symmetries */
  int surfacezones; /* surfacezone number we have for this var on each side */
  int n_special[3]; /* if >0, use this dim in dirs 0,1,2 */
  int type;         /* 0: evo. var., 1: aux. var. no need to copy or interp. */
  int Nextra;       /* extra space allocated in var array */
} tVar;

/* Extra possibilities for n_special[i] in tVar, besides any positive number.
   These are checked in enablevarcomp_innode */
enum
{
  NOT_USED=0,   /* everything else in here needs to be negative! */
  NODE_n=-1,    /* use node->n[d] in dir d */
  NODE_nM1=-2,  /* use node->n[d]-1 in dir d */
  NODE_nP1=-3   /* use node->n[d]+1 in dir d */
};

/* functions to create and access variables */
void AddMeshVar(tMesh *mesh, const char *name, const char *tensorindices,
                const char *description);
void free_mesh_vdb_contents(tMesh *mesh);
void AddEvoMeshVar(tMesh *mesh, const char *name,
                   const char *tensorindices, const char *description);
void AddAuxMeshVar(tMesh *mesh, const char *name,
                   const char *tensorindices, const char *description);
void AddMeshVarDim(tMesh *mesh, const char *name,
                   const char *tensorindices, const char *description,
                   int n_special0, int n_special1, int n_special2);
void AddAuxMeshVarDim(tMesh *mesh, const char *name,
                      const char *tensorindices, const char *description,
                      int n_special0, int n_special1, int n_special2);
int MeshVarIndLax(tMesh *mesh, const char *name);
int MeshVarInd(tMesh *mesh, const char *name);
int Set_vdb_iStart_AtVar(tMesh *mesh, const char *name);
char *MeshVarName(tMesh *mesh, int i);
int MeshVarNComponents(tMesh *mesh, int i);
int MeshVarComponent(tMesh *mesh, int i);
int MeshVarIndComponent0(tMesh *mesh, int i);
char *MeshVarNameComponent0(tMesh *mesh, const char *name);
char *MeshVarTensorIndices(tMesh *mesh, int i);
void MeshVarNameSetBoundaryInfo(tMesh *mesh, const char *name,
			        double farlimit, double falloff);
void MeshVarSetType(tMesh *mesh, int i, int type);
void MeshVarSetSurfInfo(tMesh *mesh, int i, int surfacezones);
void MeshVarSetSpecial(tMesh *mesh, int i,  int ns0, int ns1, int ns2);
void MeshVarSetNextra(tMesh *mesh, int i, int Nextra);
double MeshVarFallOff(tMesh *mesh, int i);
double MeshVarFarLimit(tMesh *mesh, int i);
int MeshVarSymmetry(tMesh *mesh, int i, int dir);
int MeshVarSurfacezones(tMesh *mesh, int i);
int MeshVarType(tMesh *mesh, int i);
int *MeshVar_n_special(tMesh *mesh, int i);
int MeshVar_Nextra(tMesh *mesh, int i);
/* convenience macros for vars */
#define VarName(i) MeshVarName(mesh, (i))
#define Ind(name)  MeshVarInd(mesh, (name))
#define AddEvoVar(name, tensorindices, description) \
  AddEvoMeshVar(mesh, (name), (tensorindices), (description))
#define AddAuxVar(name, tensorindices, description) \
  AddAuxMeshVar(mesh, (name), (tensorindices), (description))

#define AddVarDim(name, tensorindices, description, ns0,ns1,ns2) \
  AddMeshVarDim(mesh, (name), (tensorindices), (description), \
                (ns0),(ns1),(ns2))
#define AddAuxVarDim(name, tensorindices, description, ns0,ns1,ns2) \
  AddAuxMeshVarDim(mesh, (name), (tensorindices), (description), \
                  (ns0),(ns1),(ns2))
/* variable lists in variables.c*/
void prvarlist(tVarList *v);
tVarList *vlalloc(tMesh *mesh);
void vlfree(tVarList *u);
void vlpushone(tVarList *v, int vi);
void vlpush(tVarList *v, int vi);
void vlpushvl(tVarList *v, tVarList *u);
void vldropn(tVarList *v, int n);
int vlpushone_index(tVarList *v, int vi);
int vlpush_index(tVarList *v, int vi);
int vlpushvl_index(tVarList *v, tVarList *u);
int vlindex(tVarList *v, int vi);
void vlsort(tVarList *v);
int vlindex_if_sorted(tVarList *v, int vi);
tVarList *AddDuplicate(tVarList *vl, const char *postfix,
                       int type, int surfacezones);
tVarList *AddDuplicateEnable(tVarList *vl, const char *postfix,
                             int type, int surfacezones);
void vlsetconstant_node(tNode *node, tVarList *u, const double c);
void vlsetconstant(tVarList *u, const double c);
void vlcopy_node(tNode *node, tVarList *v, tVarList *u);
void vlcopy(tVarList *v, tVarList *u);
void vladd_node(tNode *node,
                tVarList *r, double ca, tVarList *a, double cb, tVarList *b);
void vladd(tVarList *r, double ca, tVarList *a, double cb, tVarList *b);
void vladdto_node(tNode *node, tVarList *r, const double ca, tVarList *a);
void vladdto(tVarList *r, const double ca, tVarList *a);
void vladdto_onfaces_node(tNode *node, tVarList *r,
                          const double ca, tVarList *a);
void vladdto_onfaces(tVarList *r, const double ca, tVarList *a);
intList *vl2intList(tVarList *v);

/* utilities.c */
NORET void  errorexit(const char *file, int line, const char *func,
                      const char *s);
NORET void errorexits(const char *file, int line, const char *func,
                      const char *s, const char *t);
NORET void errorexiti(const char *file, int line, const char *func,
                      const char *s, int i);
#define errorexit(s)     errorexit(__FILE__, __LINE__, __func__, (s))
#define errorexits(s,t) errorexits(__FILE__, __LINE__, __func__, (s), (t))
#define errorexiti(s,i) errorexiti(__FILE__, __LINE__, __func__, (s), (i))
#define PRF     printf("%s", __func__)
#define PRFs(s) printf("%s%s", __func__, s)

#define VarMem(mesh, cat,li, vi) Vard_((mesh)->myln->ln[cat][li]->node, vi)
#define prVarMem(mesh, cat,li, vi,ijk) \
  printf("cat%d node%d var%d at %p+%d = %g\n", cat,(int) (li), vi, \
         (void *) VarMem(mesh, cat,li, vi), ijk, \
         VarMem(mesh, cat,li, vi)[ijk])

void Yo(double x);
void prdivider(int n);
int getRealTime(struct timespec *tp);
void initTimeIn_s(void);
void getTimeDiff(struct timespec dtp[1],
                 struct timespec tp1[1], struct timespec tp0[1]);
double getTimeDiffIn_s(struct timespec tp1[1], struct timespec tp0[1]);
double getTimeIn_s(void);
void prTimeIn_s(const char *comment);
void wait_for_debugger_if_NMESH_MPI_DEBUG(void);
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
int finit(double x);
int remove_chars_from_str(char *str, const char *del);
void trim_whitespace(char *str);
int get_par_from_str(const char *str, char *name, const char *delim,
                     char *value, int n);
int str_to_intList(const char *str, const char *delim, intList *il);
long str_from_buf(const char *buffer, long nbuffer, long offset,
                  char delim, char *str, long nstr, long *strlen);
long nbytes_infile(FILE *fp);
FILE *fopen_buf(const char *pathname, const char *mode,
                char **buf, size_t bufsiz);
int fclose_buf(FILE *fp, char **buf);
int copy_file_into_dir(char *fname, char *dir);
int system1(const char *s1);
int system2(const char *s1, const char *s2);
int system3(const char *s1, const char *s2, const char *s3);
void print_system_info(void);
int system_emu(const char *command);
int construct_argv(char *str, char ***argv);
void sort_int_array(int n, int *ar);
int search_sorted_int_array(int n, int *ar, int key);
double *dmalloc(int n);
int *imalloc(int n);
char *cmalloc(int n);
void *pmalloc(int n);
void *dtensor(size_t size);
void *pcalloc(size_t n);
void *checked_calloc(size_t nmemb, size_t size);
void *checked_realloc(void *ptr, size_t size);
void *rows_calloc(size_t nx, unsigned long ny[nx], size_t size);
void rows_free(void *g, size_t nx);
void rows_print_sizes(size_t nx, unsigned long ny[nx], size_t size);

/* endianIO.c */
size_t fread_little(void *ptr, size_t size, size_t nmemb, FILE *fp);
size_t fwrite_little(const void *ptr, size_t size, size_t nmemb, FILE *fp);
size_t fwrite_big(const void *ptr, size_t size, size_t nmemb, FILE *fp);
size_t fread_big(void *ptr, size_t size, size_t nmemb, FILE *fp);
int return_BYTE_ORDER_LITTLE(void);
int print_endian_info(tMesh *mesh);

/* glist.c */
void glist_entry_add(void *entry, struct list_head *head);
void glist_entry_add_tail(void *entry, struct list_head *head);
void glist_elem_del(tGlist *elem);
void glist_free_elems_and_entries(struct list_head *head, void (*Free)());
void glist_free_elems(struct list_head *head);

/* timer.c */
int write_all_timers(tMesh *mesh);
int free_all_timers(tMesh *mesh);
struct tTIMER *timer_start(const char *name);
struct tTIMER *timer_stop(const char *name);
#define TIMER_START timer_start(__func__)
#define TIMER_STOP  timer_stop(__func__)

/* nan_checker.c */
int array_finite(tArray *a, char *name, int ijk[3]);
int var_finite(tNode *node, int vi);
int vl_finite(tNode *node, tVarList *vl);
int vl_finite_mesh(tVarList *vl);
int nan_checker(tMesh *mesh);

/* norms.c */
double MeshVolumeIntegral(tMesh *mesh, tPat *pat, int vind,
                          double power, int mode);
double MeshMax(tMesh *mesh, tPat *pat, int vind);
double MeshMin(tMesh *mesh, tPat *pat, int vind);
double MeshExtremumLoc(tMesh *mesh, tPat *pat, int vind, int findMax,
                       int *Mp, char Mnodeloc[104], int *Mijk,
                       double *MX, double *Mx);

/* main.c */
FILE *fopen_bufsize(tMesh *mesh, const char *pathname, const char *mode,
                    char **buf);
int fs_sync(tMesh *mesh);
void finalize_all_and_exit(tMesh *mesh, int ec);
