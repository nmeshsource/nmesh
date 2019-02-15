/* main.c */
/* Wolfgang Tichy, Jan. 2019 */

#include "nmesh.h"
#include "main.h"
#include <time.h> // for time and ctime functions

/* global var that contain the mesh made in main */
tMesh *main_mesh;


/* initialize libraries 
   the automatically generated file calls the initializers for each module */
void initialize_libraries(struct tMESH *mesh)
{
  prdivider(1);
  printf("Initializing libraries\n");

#include "nmesh_automatic_initialize.c"
}

/**************************************************************************/
/* main */
int main(int argc, char **argv) 
{
  tMesh *mesh;
  
  nMPI_Init(&argc, &argv);

  /* make first mesh in which we store pars, vars and funs */
  mesh = make_empty_mesh(0);
  main_mesh = mesh; /* save this mesh in a global var as well */

  /* start the main parts of nmesh: */
  initTimeIn_s();
  read_command_line(mesh, argc, argv);
  parse_parameter_file(mesh, Gets(Par("parameterfile")));
  parse_command_line_options(mesh);
  make_output_directory(mesh);
  initialize_libraries(mesh);

  iterate_parameters(mesh, 0); /* start of new iteration */
  while(iterate_parameters(mesh, 1))
  {
    RunFun(POST_PARAMETERS); //hook for funs right after iterate_parameters
    RunFun(INITMESH); // here we schedule funcs to programatically set up the mesh
    inidata_mesh(mesh);
    evolve_mesh(mesh);
    makeparameter(mesh, "outdir_previous_iteration", "",
                  "outdir of previous iteration");
    Sets(Par("outdir_previous_iteration"), Gets(Par("outdir")));
  }

  finalize_mesh(mesh);
  nMPI_Finalize();
  return 0;
}


/* read command line */
int read_command_line(tMesh *mesh, int argc, char **argv)
{
  int i; 

  if(0) 
    for(i = 0; i < argc; i++) printf("argv[%d] = %s\n", i, argv[i]);

  if(Rank0)
  {
    prdivider(1);
    printf("Welcome to nmesh, compiled on %s at %s\n", __DATE__, __TIME__);
    prdivider(1);
  }

  if(argc < 2)
  {
    if(Rank0)
    {
      printf("Usage:  nmesh name.par\n");
      printf("or:     nmesh name.par options and extra arguments\n");
      printf("\n");
      printf("options: --keep_previous          do not touch name_previous\n");
      printf("         --modify-par:\"P=v\"       set par P to value v\n");
      printf(" all options must start with --\n");
    }
    nMPI_Finalize();
    exit(0);
  }

  /* this par is needed here already, so that finalexit can check it */
  if(Rank0) printf("Making first parameters\n");
  makeparameter(mesh, "errorexit", "exit", "how we exit");
  /* this is about how we output */
  makeparameter(mesh, "logfile_creation", "append", "how we create logfile");

  /* got two or more arguments? */
  if (argc >= 2)
  {
    int nopts, nargs;
    char argi[1000];
    char descr[1000];
    char options[1000];
    char *parfile = (char *) calloc(strlen(argv[1])+40, sizeof(char));
    char *outdir  = (char *) calloc(strlen(argv[1])+40, sizeof(char));
    int parnamelen;

    /* determine name of parameter file and output directory */
    strcpy(parfile, argv[1]);
    parnamelen = strlen(parfile);
    if(!strstr(parfile, ".par") || parnamelen<5) strcat(parfile, ".par");
    parnamelen = strlen(parfile);
    if(!strstr(parfile+parnamelen-4, ".par")) strcat(parfile, ".par");
    strcpy(outdir, parfile);
    outdir[strlen(outdir)-4]=0; /* remove .par */

    /* first parameter initializes parameter data base */
    if(Rank0) printf("Adding command line parameters\n");
    makeparameter(mesh, "outdir", outdir, "output directory");
    makeparameter(mesh, "parameterfile", parfile, 
                  "name of parameter file given on command line");

    /* add other args */
    nargs=nopts=0;
    options[0]=0;
    for(i=2; i<argc; i++)
    {
      if(argv[i][0]=='-')
      {
        if(strlen(argv[i])==1) errorexit("- is not a valid option");
        if(argv[i][1]!='-') errorexit("all options must start with --");
        /* save all options in options */
        if(nopts>0) strncat(options, " ", 999);
        strncat(options, argv[i], 999);
        nopts++;
      }
      else
      {
        snprintf(argi, 999, "nmesh_arg%d", nargs+2);
        snprintf(descr, 999, "nmesh command line argument%d", nargs+2);
        makeparameter(mesh, argi, argv[i], descr);
        nargs++;
      }
    }
    /* add nmesh command line options */
    if(nopts>0)
      makeparameter(mesh, "nmesh_options", options,
                    "nmesh command line options");
    free(parfile);
    free(outdir);
  }

  /* more initialization */
  
  return 0;
}


int make_output_directory(tMesh *mesh)
{
  char *outdir  = Gets(Par("outdir"));
  char *outdirp = (char *) calloc(strlen(outdir)+40, sizeof(char));
  char so[1000];
  time_t mytime = time(NULL);       // get time
  char *time_str = ctime(&mytime);

  /* set outdirp to outdir_previous */
  strcpy(outdirp, outdir);
  strcat(outdirp, "_previous");

  /* check if a shell is available to execute commands later */
  /* NOTE: system2 and system3 are smart enough to do "mkdir", "rm -rf"
           and "mv" even without a shell by using POSIX calls. */
  if(system(NULL)==0 && Rank0)
  {
    printf("WARNING: system(NULL)=0 => cannot execute shell commands!\n");
    printf("         Consider using system_emu.\n");
  }

  /* check if we remove outdir_previous */
  if(!GetvLax(Par("nmesh_options"), "--keep_previous"))
  {
    /*
    char *prev = checkpoint_filename("_previous", "");
    char *curr = checkpoint_filename("", "");
    FILE *fpprev = fopen(prev, "r");
    FILE *fpcurr = fopen(curr, "r");

    //printf("prev: %s, curr: %s\n", prev, curr);
    if(0) printf("Safety first: checking data before overwriting "
                 "previous directory.\n");
    if(fpprev && !fpcurr)
    {
        printf("Warning: %s exists, while %s does not.\n", prev, curr);
        printf("Warning: This could result in the loss of the checkpoint!\n");
        errorexit("There may be important data in the previous directory!");
    }
    if(fpprev) fclose(fpprev);
    if(fpcurr) fclose(fpcurr);
    free(prev);
    free(curr);
    */
    /* remove outdir_previous and move outdir to outdir_previous */
    if(Rank0)
    {
      system2("rm -rf", outdirp);
      system3("mv", outdir, outdirp);
    }
  }

  /* make output directory, save parfile */
  if(Rank0)
  {
    system2("mkdir", outdir);
    copy_file_into_dir(Gets(Par("parameterfile")), outdir);
  }
  /* all wait here until mkdir is done. */
  nMPI_barrier();

  /* redirect stdout and stderr. Do it for all MPI ranks>0
     all output is collected in outdir/stdout.001 etc  */
  if(nMPI_rank()>0)
  {
    char f[100];
    snprintf(f,99, "%%s/stdout.%%0%dd", (int) log10(nMPI_size())+1);
    snprintf(so,999, f, outdir, nMPI_rank());
    if(nMPI_rank()==1)
    {
      prdivider(3);
      printf("*** NOTE ***  Output from other procs redirected to e.g.:\n %s\n",
             so);
      prdivider(3);
    }
    freopen(so, "w", stdout);
    freopen(so, "w", stderr);
  }
  else if(!Getv(Par("logfile_creation"),"no"))
  {
    char *opt;
    snprintf(so,999, "%s.log", outdir);
    prdivider(3);
    printf("*** NOTE ***  Output redirected to:\n %s\n", so);
    prdivider(3);
    if(Getv(Par("logfile_creation"),"append"))
      opt = "a";
    else
      opt = "w";
    freopen(so, opt, stdout);
    freopen(so, opt, stderr);
  }
  /* say what we have after redirection: */
  prdivider(1);
  time_str[strlen(time_str)-1] = '\0';
  printf("The current time is %s\n", time_str);
  printf("nmesh was compiled on %s at %s\n", __DATE__, __TIME__);
  printf("  outdir = %s \n", outdir);

  free(outdirp);
  return 0;
}


/* go through options and act accordingly */
int parse_command_line_options(tMesh *mesh)
{
  char *optionstr;
  char *str1;
  char *str2;
  char *par;
  char *val;

  //PRF;printf(":\nPar(\"nmesh_options\")=%i\n", Par("nmesh_options"));
  /* get length of nmesh_options string */
  if(GetsLax(Par("nmesh_options"))==0) return 0;
  if(Rank0) printf("Parsing command line options\n");

  /* parse for all --modify-par: */
  optionstr = strdup((GetsLax(Par("nmesh_options"))));
  str1 = optionstr;
  while( (str1=strstr(str1, "--modify-par:"))!=NULL )
  {
    str1+=13;
    str2=strstr(str1, "=");
    if(str2==NULL) break;
    str2[0]=0; /* replace = with 0 */
    par=str1;
    val=str2+1;
    str2=strstr(val, " --");
    if(str2!=NULL) str2[0]=0; /* now val is 0 terminated for sure */
    if(0) printf("par=%s|val=%s|\n", par,val);
    
    /* set par to new value */
    AddOrModifyPar(par, val, "set with --modify-par option");

    /* move forward in str1 */
    if(str2!=NULL) str2[0]=' '; /* restore space before -- */
    str1=val;
    if(0) printf("str1=%s|\n", str1);
  }
  free(optionstr);
  return 0;
} 


/* get initial data for mesh */
int inidata_mesh(tMesh *mesh)
{
  if(1)
  {
    prdivider(0);
    printf("Initializing mesh\n");
  }

  /* compute initial data */
  RunFun(INITIALDATA);

  /* initial data is just another new time slice */
  RunFun(POST_EVOLVE);

  /* initial data complete */
  prdivider(0);
  printf("Done with initialization\n");
  printf(" iteration %d, time=%g\n", mesh->iteration, mesh->time);

  /* analyze initial data */
  RunFun(ANALYZE);

  /* output for permanent variables */
  RunFun(OUTPUT);

  /* checkpoint, in case we init from checkpoint files */
  //checkpoint(mesh);

  return 0;
}


/* evolve mesh */
int evolve_mesh(tMesh *mesh)
{
  int iterationmax = Geti(Par("iterations"));
  double timemax   = Getd(Par("finaltime"));

  prdivider(0);

  if(timemax > 0)
    iterationmax = timemax/mesh->dt + 0.5;

  if(iterationmax > 0) 
    printf("Evolving mesh for %d iterations to time %g\n", 
	   iterationmax, iterationmax * mesh->dt);

  if(iterationmax <= 0) return 0;

  /* outermost evolution loop */
  while(mesh->iteration < iterationmax)
  { 
    /* pre evolve */
    RunFun(PRE_EVOLVE); 

    /* make one evolution step */
    RunFun(EVOLVE); 

    /* post evolve */
    RunFun(POST_EVOLVE); 

    /* the evolution step is complete now */
    mesh->iteration++;
    mesh->time = mesh->iteration * mesh->dt;

    /* print some info */
    printf(" iteration %d, time=%g\n", mesh->iteration, mesh->time);
    fflush(stdout); 

    /* call analyze functions */
    RunFun(ANALYZE);

    /* call output functions, say for variable output */
    RunFun(OUTPUT);

    /* post output */
    RunFun(POST_OUTPUT);

    /* checkpoint */
    //checkpoint(mesh);
  
    /* update since this may change during evolution, say when checkpointing */
    timemax = Getd(Par("finaltime"));
    iterationmax
      = (timemax > 0) ? timemax/mesh->dt + 0.5 : Geti(Par("iterations"));
  }
  return 0;
}



/* finalize mesh */
int finalize_mesh(tMesh *mesh)
{
  free_mesh(mesh);
  return 0;
}
