/* utilities.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "main.h"

#include <time.h>
#include <ctype.h>       /* for isspace */

/* for POSIX.1-2001 mkdir, opendir, fork, wait functions */
#include <unistd.h>      /* for fork */
#include <sys/stat.h>
#include <sys/types.h>   /* for pid_t */
#include <sys/wait.h>    /* for wait */
#include <dirent.h>      /* for opendir */
#include <fenv.h>        /* for feenableexcept */


/* global vars for timing */
struct timespec tp_at_nmesh_start[1]; /* initTimeIn_s is called in main.c */
extern tMesh *main_mesh;              /* mesh created in main */


/* quick debug */
void Yo(double x) {fprintf(stdout, "Yo:%g\n", x);fflush(stdout);}


/* print divider lines */
void prdivider(int n)
{
  int i, c;
  switch(n)
  {
  case 1:
    printf("==============================================================================\n");
    break;
  case 2:
    printf("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n");
    break;
  case 3:
    printf("******************************************************************************\n");
    break;
  default:
    if(n>=' ' && n<='~') c = n;
    else                 c = '-';
    for(i=0; i<78; i++) printf("%c", c);
    printf("\n");
  }
  fflush(stdout);
}


/* wrapper around clock_gettime */
int getRealTime(struct timespec *tp)
{
  return clock_gettime(CLOCK_REALTIME, tp);
}

/* initialize time_in_s_at_nmesh_start */
void initTimeIn_s(void)
{
  struct timespec tp[1];

  getRealTime(tp);
  tp_at_nmesh_start->tv_sec  = tp->tv_sec;
  tp_at_nmesh_start->tv_nsec = tp->tv_nsec;
}

/* get time difference: dtp = tp1 - tp0 */
void getTimeDiff(struct timespec dtp[1],
                 struct timespec tp1[1], struct timespec tp0[1])
{
  dtp->tv_sec  = (tp1->tv_sec  - tp0->tv_sec);
  dtp->tv_nsec = (tp1->tv_nsec - tp0->tv_nsec);
}

/* add time dtp to tp: tp = tp + dtp */
void addtoTime(struct timespec tp[1], struct timespec dtp[1])
{
  tp->tv_sec  += dtp->tv_sec;
  tp->tv_nsec += dtp->tv_nsec;
}

/* get time difference in seconds */
double getTimeDiffIn_s(struct timespec tp1[1], struct timespec tp0[1])
{
  double t_in_s;

  t_in_s  = (tp1->tv_sec - tp0->tv_sec);
  t_in_s += 1e-9 * (tp1->tv_nsec - tp0->tv_nsec);

  return t_in_s;
}

/* get current time in seconds */
double getTimeIn_s(void)
{
  struct timespec tp[1];
  getRealTime(tp);
  return getTimeDiffIn_s(tp, tp_at_nmesh_start);
}

/* print current time */
void prTimeIn_s(const char *comment)
{
  double t_in_s = getTimeIn_s();
  printf("%s%gs\n", comment, t_in_s);
  fflush(stdout);
}

/* get current CPU clock time in seconds */
double getClockTimeIn_s(void)
{
  double t_in_s = clock();
  t_in_s = t_in_s/CLOCKS_PER_SEC;
  return t_in_s;
}

/* print current CPU clock time */
void prClockTimeIn_s(char *comment)
{
  double t_in_s = clock();
  t_in_s = t_in_s/CLOCKS_PER_SEC;
  printf("%s%gs\n", comment, t_in_s);
  fflush(stdout);
}

/* for MPI debugging */
/* wait until var i is set to 1 in debugger */
void wait_for_debugger_if_NMESH_MPI_DEBUG(void)
{
#ifdef USEMPI
  if(getenv("NMESH_MPI_DEBUG") != NULL)
  {
    long pid = getpid();
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank==0)
    {
      volatile int i=0;

      sleep(4);
      printf("nmesh rank%d has pid %ld and is waiting for debugger.\n",
             rank, pid);
      printf("Use the following commands in gdb to stop waiting:\n");
      printf(" attach %ld\n", pid);
      printf(" set var i=1\n");
      fflush(stdout);

      /* infinite loop */
      while(i==0) { /* set var i=1 to leave this infinite loop */ }
    }
    else
    {
      printf("nmesh rank%d has pid %ld.\n", rank, pid);
      fflush(stdout);
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }
#endif
}


/* minimum and maximum funcs, works for integers if they are not too big */
double min2(double x, double y)
{
  return (x < y) ? x : y;
}
double min3(double x, double y, double z)
{
  return min2(min2(x, y), z);
}
double max2(double x, double y)
{
  return (x > y) ? x : y;
}
double max3(double x, double y, double z)
{
  return max2(max2(x, y), z);
}


/* find min in a 1d array f, return min and set imin to index of min */
double min_in_1d_array(double *f, int n, int *imin)
{
  int i;   
  double fmin=f[0];

  *imin=0;
  for(i=1; i<n; i++)
    if(f[i]<fmin) 
    {
      fmin = f[i];
      *imin= i;
    }
  return fmin;
}
/* find max in a 1d array f, return max and set imax to index of max */
double max_in_1d_array(double *f, int n, int *imax)
{
  int i;   
  double fmax=f[0];

  *imax=0;
  for(i=1; i<n; i++)
    if(f[i]>fmax) 
    {
      fmax = f[i];
      *imax= i;
    }
  return fmax;
}
/* find min in two 1d arrays f0, f1, 
   return min, set ai to 0 or 1 depending on which array the min is in,
   set imin to index of min in array ai */
double min2_in_1d_array(double *f0, int n0, double *f1, int n1, 
                        int *ai, int *imin)
{
  int imin0, imin1;   
  double fmin0, fmin1, fmin;
  
  fmin0 = min_in_1d_array(f0, n0, &imin0);
  fmin1 = min_in_1d_array(f1, n1, &imin1);

  if(fmin1<fmin0) { *ai=1;  fmin=fmin1;  *imin=imin1; }
  else            { *ai=0;  fmin=fmin0;  *imin=imin0; }  
  return fmin;
}
/* find max in two 1d arrays f0, f1, 
   return max, set ai to 0 or 1 depending on which array the max is in,
   set imax to index of max in array ai */
double max2_in_1d_array(double *f0, int n0, double *f1, int n1, 
                        int *ai, int *imax)
{
  int imax0, imax1;   
  double fmax0, fmax1, fmax;
  
  fmax0 = max_in_1d_array(f0, n0, &imax0);
  fmax1 = max_in_1d_array(f1, n1, &imax1);

  if(fmax1>fmax0) { *ai=1;  fmax=fmax1;  *imax=imax1; }
  else            { *ai=0;  fmax=fmax0;  *imax=imax0; }  
  return fmax;
}
/* find min in three 1d arrays f0, f1, f2, 
   return min, set ai to 0,1,2 depending on which array the min is in,
   set imin to index of min in array ai */
double min3_in_1d_array(double *f0, int n0, double *f1, int n1, double *f2, int n2,
                        int *ai, int *imin)
{
  double fmin;

  min2_in_1d_array(f1,n1, f2,n2, ai, imin);
  if(*ai==1) { fmin = min2_in_1d_array(f0,n0, f2,n2, ai, imin); *ai = (*ai)*2; }
  else       { fmin = min2_in_1d_array(f0,n0, f1,n1, ai, imin); }
  return fmin;
}
/* find max in three 1d arrays f0, f1, f2, 
   return max, set ai to 0,1,2 depending on which array the max is in,
   set imax to index of max in array ai */
double max3_in_1d_array(double *f0, int n0, double *f1, int n1, double *f2, int n2,
                        int *ai, int *imax)
{
  double fmax;

  max2_in_1d_array(f1,n1, f2,n2, ai, imax);
  if(*ai==1) { fmax = max2_in_1d_array(f0,n0, f2,n2, ai, imax); *ai = (*ai)*2; }
  else       { fmax = max2_in_1d_array(f0,n0, f1,n1, ai, imax); }
  return fmax;
}

/* Read bit pattern in x to determine if it is finite, because
   isfinite(x) does not always work with -Ofast or -ffast-math */
int finit(double x)
{
  union {
    uint64_t bits; /* must be same number of bytes as double */
    double d;
  } ud;
  ud.d = x;

  //printf(" %lx ", ud.bits);

  /* According to IEEE 754:
     NaN starts with 7FF0, 7FF8 or 7FFF. -NaN starts with FFF0, FFF8 or FFFF.
     Inf starts with 7FF0. -Inf starts with FFF0. */
  if((ud.bits | 0x800FFFFFFFFFFFFF) == 0xFFFFFFFFFFFFFFFF)
    return 0;
  else
    return 1;
}

/* if MyConfig contains, e.g.
   DFLAGS += -D_GNU_SOURCE -DFPEEXCEPTIONS="FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW"
   we enable the corresponding floating point exceptions */
int enable_FPEEXCEPTIONS(void)
{
  PRF;
#ifdef FPEEXCEPTIONS
  printf(": feenableexcept(%s);\n", MSTR_OFVAL(FPEEXCEPTIONS));
  return feenableexcept(FPEEXCEPTIONS);
#else
  printf(": none.\n");
  return 0;
#endif
}


/* remove all chars that occur in del from string str,
   return number of removed chars */
int remove_chars_from_str(char *str, const char *del)
{
  unsigned writ = 0, read = 0;
  unsigned nstr = strlen(str);
  unsigned ndel = strlen(del);
  unsigned n;
  int keep;

  for(read=writ=0; read<nstr; read++)
  {
    /* do we keep char? */
    keep = 1;
    for(n=0; n<ndel; n++)
      if(str[read] == del[n]) { keep=0; break; }

    /* copy char we keep */
    if(keep)
      str[writ++] = str[read];
  }
  /* add nul at end of string */
  str[writ]=0;
  return nstr - writ;
}

/* remove all trailing and leading whitespaces */
void trim_whitespace(char *str)
{
  int len, f,i;

  if(!str) return;

  len = strlen(str);
  //printf("str=|%s| len=%d\n", str, len);

  /* remove all trailing spaces */
  for(f=len-1; f>=0; f--) if(!isspace(str[f])) break;
  len = f+1;
  str[len] = 0;
  //printf("str=|%s| len=%d\n", str, len);

  /* find first no-space char */
  for(f=0; f<len; f++) if(!isspace(str[f])) break;
  len = len-f;

  /* shift all chars left by f */
  if(f) for(i=0; i<=len; i++) str[i] = str[i+f];
  //printf("str=|%s| len=%d\n", str, len);
}

/* parse a string to find a parname and its value,
   returns 1 if str contains the delimiter delim (e.g. "=") otherwise 0. */
int get_par_from_str(const char *str, char *name, const char *delim,
                     char *value, int n)
{
  char *str2, *saveptr, *nam, *val;
  int ret;

  /* duplicate, because strtok_r writes into str2 */
  str2 = strdup(str);

  /* find parname and its value */
  nam = strtok_r(str2, delim, &saveptr);
  val = strtok_r(NULL, delim, &saveptr);

  if( (!(*nam)) || (strlen(str2)==strlen(str)) )
  {
    ret = 0;
  }
  else
  {
    /* trim spaces */
    trim_whitespace(nam);
    trim_whitespace(val);

    strncpy(name,  nam, n);
    if(val) strncpy(value, val, n);
    else    value[0] = 0;
    ret = 1;
  }

  free(str2);
  return ret;
}

/* convert a string to an intList */
int str_to_intList(const char *str, const char *delim, intList *il)
{
  char *str2, *sav, *val;
  int ret = 0;

  /* duplicate, because strtok_r writes into str2 */
  str2 = strdup(str);
  for(val=strtok_r(str2, delim, &sav); val!=NULL;
      val=strtok_r(NULL, delim, &sav))
  {
    int i = atoi(val);
    push_intList(il, i);
    ret++;
  }
  free(str2);
  return ret;
}

/* start at buffer + offset and read into array mem until first delim,
   return offset to mem after last read, write length of mem read in len */
long mem_from_buf(const char *buffer, long nbuffer, long offset,
                  char delim, char *mem, long memsize, long *len)
{
  const char *buf;
  long nbuf, i, im, newoffset;

  if(offset>=nbuffer) return -1; /* signal end of buffer */
  if(offset<0)        return -1; /* if offset<0 stay at end of buffer */

  buf = buffer + offset;
  nbuf = nbuffer - offset;

  if(nbuf < memsize) im = nbuf;
  else               im = memsize;

  for(i=0; i<im; i++)
  {
    mem[i] = buf[i];
    if(buf[i] == delim) { i++;  break; }
  }
  newoffset = offset + i;
  *len = i;

  return newoffset;
}

/* start at buffer + offset and read into string str until first delim,
   return offset to mem after last read, write length of mem read in len:
   To be used like e.g. this:
     while((off = str_from_buf(buffer,nbuffer, off, '\n', buf,999, &len))>=0)
   or:
     off = str_from_buf(buffer,nbuffer, off, '\n', buf,999, &len);  */
long str_from_buf(const char *buffer, long nbuffer, long offset,
                  char delim, char *str, long nstr, long *strlen)
{
  long newoffset;
  newoffset = mem_from_buf(buffer,nbuffer,offset, delim, str, nstr-1, strlen);
  str[*strlen] = 0; /* add string terminator */
  return newoffset;
}


/* return the number of bytes in a file */
long nbytes_infile(FILE *fp)
{
  long nbytes;

  /* find number of bytes in file */
  fseek(fp, 0L, SEEK_END);
  nbytes = ftell(fp);
  fseek(fp, 0L, SEEK_SET);
  return nbytes;
}


/* Open a file where we use a buffer of size bufsiz (set with setvbuf).
   The buffer buf is allocated here and needs to be freed later, e.g.
   by calling fclose_buf. If bufsiz=0 no buffer is allocated and the default
   buffers of the C library are used. */
FILE *fopen_buf(const char *pathname, const char *mode,
                char **buf, size_t bufsiz)
{
  FILE *fp = fopen(pathname, mode);
  if(fp)
  {
    if(bufsiz)
    {
      *buf = cmalloc(bufsiz);
      setvbuf(fp, *buf, _IOFBF, bufsiz);
    }
    else
    {
      *buf = NULL;
    }
  }
  return fp;
}

/* counterpart to fopen_buf, closes file and frees buf */
int fclose_buf(FILE *fp, char **buf)
{
  if(fp)
  {
    int ret = fclose(fp);
    free(*buf);
    *buf = NULL;
    return ret;
  }
  return EOF;
}


/* make copy of a file: cp fname newname */
int copy_file(char *fname, char *newname)
{
  FILE *in, *out;
  void *buffer;
  size_t BUFSIZE=16777216; /* 16MiB */
  size_t bcount;
  size_t bufsiz = BUFSIZE/16;
  char *inbuf;
  char *outbuf;

  printf("copy_file(\"%s\", \"%s\");\n", fname, newname);

  /* open source file */
  in = fopen_buf(fname, "rb", &inbuf,bufsiz);
  if(!in) errorexits("failed opening %s", fname);

  /* open destination file */
  out = fopen_buf(newname, "wb", &outbuf,bufsiz);
  if(!out) errorexits("failed opening %s", newname);

  /* copy char by char */
  /* while( (ch=fgetc(in)) != EOF)
       fputc(ch, out);             */
  buffer = malloc(BUFSIZE * sizeof(char));
  if(!buffer)
    errorexiti("copy_file: out of memory for buffer (%d chars)", BUFSIZE);
  do
  {
    bcount = fread(buffer, sizeof(char), BUFSIZE, in);
    fwrite(buffer, sizeof(char), bcount, out);
  } while(bcount==BUFSIZE);
  free(buffer);

  fclose_buf(out, &outbuf);
  fclose_buf(in, &inbuf);
  return 0;
}

/* make a copy of a file in some dir */
int copy_file_into_dir(char *fname, char *dir)
{
  char newname[10000];
  int i;

  /* find / in fname, to determine filename without dirname */
  for(i=strlen(fname)-1; i>=0; i--) if(fname[i]=='/') break;
  snprintf(newname, 9999, "%s/%s", dir, fname+i+1);

  return copy_file(fname, newname);
}

/* use opendir to scan through dir and remove the entire dir */
int remove_dir(const char *dirname)
{
  DIR *d;
  struct dirent *dir;
  char file[8192];

  /* open the dir dirname */
  d = opendir(dirname);
  if(d)
  {
    while((dir = readdir(d)) != NULL)
    {
      DIR *sd;

      /* exclude . and .. directories */
      if( strcmp(dir->d_name, ".") == 0 || strcmp(dir->d_name, "..") == 0 )
        continue;

      snprintf(file, 8192, "./%s/%s", dirname, dir->d_name);

      sd = opendir(file);
      if(sd!=NULL)
      {
        remove_dir(file);
        closedir(sd);
      }
      else
      {
        if(remove(file) != 0)
        {
          //errorexits("remove(%s) failed!", file);
          closedir(d);
          return -2;
        }
      }
    } /* end of while loop */

    /* now close the dir dirname */
    closedir(d);

    /* delete the now empty dirname directory */
    if(remove(dirname) != 0)
    {
      //errorexits("remove(%s) failed!", dirname);
      return -1;
    }
  }
  return 0;
}

/* call system with 1, 2 or 3 arguments, also call c-functions directly
   in some special case to avoid some unix shell commands */
int system1(const char *s1)
{
  return system2(s1, "");
}
int system2(const char *s1, const char *s2)
{
  return system3(s1, s2, "");
}
int system3(const char *s1, const char *s2, const char *s3)
{
  char command[10000];
  int status = 0;

  /* check for special cases where we can use c-functions */
  if( strcmp(s1,"mv")==0 || strcmp(s1,"mv -f")==0 ) /* use rename */
  {
    sprintf(command, "rename(\"%s\", \"%s\");", s2, s3);
    printf("ANSI C call: %s\n", command);
    status = rename(s2, s3);
  }
  else if( strcmp(s1,"rm -rf")==0 ) /* use remove */
  {
    if(strlen(s2)>0)
    {
      printf("remove_dir(\"%s\");\n", s2);
      status = remove_dir(s2);
    }
    if(strlen(s3)>0)
    {
      printf("remove_dir(\"%s\");\n", s3);
      status = remove_dir(s3);
    }
  }
  else if( strcmp(s1,"mkdir")== 0 ) /* use POSIX.1-2001 mkdir function */
  {
    sprintf(command, "mkdir(\"%s\", S_IRWXU | S_IRWXG);", s2);
    printf("POSIX.1-2001 call: %s\n", command);
    status = mkdir(s2, S_IRWXU | S_IRWXG);
  }
  else /* use system */
  { 
    sprintf(command, "%s %s %s", s1, s2, s3);
    printf("System call: %s\n", command);
    fflush(stdout);
    status = system(command);
    fflush(stdout);
  }
  
  if(status!=0) printf(" -> WARNING: Return value = %d\n", status);
  return status;
}

/* print some system info */
void print_system_info(void)
{
  char str[1024];
  long pid = getpid();

  prdivider(0);
  printf("print_system_info: calling some shell commands\n");
  system1("hostname");
  system1("uname -a");
  system1("uptime");
  system1("lscpu");
  snprintf(str,1023, "cat /proc/%ld/status", pid);
  system1(str);
  system1("free -h");
}

/* construct an argv array from a string and return number of args */
/* NOTE: str is modified and used as mem for argv! */
int construct_argv(char *str, char ***argv)
{
  char *str1, *token, *saveptr;
  int count;

  *argv = NULL;
  for(count=0, str1=str; ; count++, str1=NULL)
  {
    *argv = (char **) realloc(*argv, sizeof(char *)*(count+1));
    token = strtok_r(str1, " ", &saveptr);
    //printf("token=%p:%s\n", token,token);
    (*argv)[count] = token;
    if(token == NULL) break;
  }
  //printf("saveptr=%p:%s\n", saveptr,saveptr);
  return count;
}

/* run a command, without a shell */
int system_emu(const char *command)
{
  char *com = strdup(command); /* duplicate since construct_argv modifies its args */
  int ret, status;
  pid_t cpid;
  printf("system_emu: running command:\n%s\n", command);

  /* Spawn a child to run the program. */
  cpid = fork();
  if(cpid<0) /* fork failed */
  {
    printf("*** WARNING: fork failed! ***\n");
    status = ret = -911;
  }
  else if(cpid==0) /* child process */
  {
    char **argv;
    construct_argv(com, &argv);
    ret = execv(argv[0], argv);
    printf("*** WARNING: command not found, (execv returned %d) ***\n", ret);
    exit(127); /* exit child, only if execv fails */
  }
  else /* cpid!=0; parent process */
  {
    int wret = waitpid(cpid, &ret, 0); /* wait for child to exit */
    if(wret<0)
    {
      printf("*** WARNING: waitpid failed! ***\n");
      status = ret = -42;
    }  
    else
      status = ret;
    //printf("wret=%d  ret=%d  status=%d\n", wret, ret, status);
  }
  if(status!=0) printf(" -> WARNING: Return value = %d\n", status);
  free(com);
  return status;
}

/* Lock a file from current file position to the end. The lock will be
   released when the file is closed.
   fd is a file descriptor open for writing. */
int lock_curr_til_EOF(FILE *out)
{
  int fd = fileno(out); /* get file dscriptor */
  if(fd==-1) return fd; /* return -1 on error */
  return lockf(fd, F_LOCK, 0);
}
/* Unlock a file from current file position to the end.
   fd is a file descriptor open for writing. */
int unlock_curr_til_EOF(FILE *out)
{
  int fd = fileno(out); /* get file dscriptor */
  if(fd==-1) return fd; /* return -1 on error */
  return lockf(fd, F_ULOCK, 0);
}


/* function to compare two ints for qsort */
int qsort_compar_int(const void *x1, const void *x2)
{
  const int *p1 = x1;
  const int *p2 = x2;
  const int i1 = *p1;
  const int i2 = *p2;

  /* without overflow we could use: return i1 - i2; */
  if(i1 < i2) return -1;
  if(i1 > i2) return +1;
  return 0;
}
/* sort array a of n integers */
void sort_int_array(int n, int *ar)
{
  qsort(ar, n, sizeof(ar[0]), qsort_compar_int);
}

/* search for key in a sorted int array, returns position in array */
int search_sorted_int_array(int n, int *ar, int key)
{
  int *p = bsearch(&key, ar, n, sizeof(ar[0]), qsort_compar_int);
  if(p) return (int) (p-ar);
  else  return -1;
}


/* malloc memory (and check out of memory) for some simple cases */
double *dmalloc(int n)
{
  double *p = malloc(sizeof(double) * n);
  
  if(!p) errorexiti("out of memory (%d double)", n);
  return p;
}

int *imalloc(int n)
{
  int *p = malloc(sizeof(int) * n);
  
  if(!p) errorexiti("out of memory (%d int)", n);
  return p;
}

char *cmalloc(int n)
{
  char *p = malloc(sizeof(char) * n);
  
  if(!p) errorexiti("out of memory (%d char)", n);
  return p;
}

void *pmalloc(int n)
{
  void *p = malloc(sizeof(void *) * n);
  
  if(!p) errorexiti("out of memory (%d void *)", n);
  return p;
}

/* malloc memory (and check out of memory) for VLA matrix or tensor.
   use it like this:
   double (*M)[ny] = dtensor(nx*ny);        // gives M[nx][ny]
   double (*T)[ny][nz] = dtensor(nx*ny*nz); // gives T[nx][ny][nz]
   //...
   free(T); free(M); */
void *dtensor(size_t size)
{
  void *p = malloc(sizeof(double) * size);

  if(!p) errorexiti("out of memory (%d double)", size);
  return p;
}

/********************************************/
/* Functions that have to do with errorexit */
/********************************************/

/* function that selects how we exit inside errorexit */
NORET void finalexit(int ec)
{
  tMesh *mesh = main_mesh;
  fflush(stderr);
  fflush(stdout);
  sync();
  nMPI_Abort(ec);
  if(mesh && GetvLax(Par("errorexit"), "abort"))  abort();
  else                                            exit(ec);
}

/* errorexit functions */
/* note that nmesh_main.h defines a macro so that the user does not have
   to specify __FILE__ , __LINE__ and __func__ to describe where the
   error occured */
#undef errorexit
#undef errorexits
#undef errorexiti

NORET void errorexit(const char *file, int line, const char *func,
                     const char *s)
{
  fprintf(stdout, "%s:%d: error in %s\n", file, line, func);
  fprintf(stdout, "Error: %s\n", s);
  finalexit(1);
}

NORET void errorexits(const char *file, int line, const char *func,
                      const char *s, const char *t)
{
  fprintf(stdout, "%s:%d: error in %s\n", file, line, func);
  fprintf(stdout, "Error: ");
  fprintf(stdout, s, t);
  fprintf(stdout, "\n");
  finalexit(1);
}

NORET void errorexiti(const char *file, int line, const char *func,
                      const char *s, int i)
{
  fprintf(stdout, "%s:%d: error in %s\n", file, line, func);
  fprintf(stdout, "Error: ");
  fprintf(stdout, s, i);
  fprintf(stdout, "\n");
  finalexit(1);
}
/************************************************************************/
/* NOTE: DO NOT WRITE ANYTHING BELOW THIS LINE!!!
   Reason: The method of using macros and undef/define for the errorexit
   functions works only if they are at the end. */
/************************************************************************/
