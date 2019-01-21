/* main.c */
/* Wolfgang Tichy, Jan. 2019 */

#include "nmesh.h"
#include "main.h"


/* initialize libraries 
   the automatically generated file calls the initializers for each module */
void initialize_libraries(struct tMESH *mesh)
{
  prdivider(0);
  printf("Initializing libraries\n");

#include "nmesh_automatic_initialize.c"
}

/**************************************************************************/
/* main */
int main(int argc, char **argv) 
{
  tMesh *m;
  
  nmesh_MPI_Init(&argc, &argv);

  /* make first mesh in which we store pars, vars and funs */
  m = make_empty_mesh(1);

  /* start the main parts of nmesh: */
  initTimeIn_s();
  read_command_line(m, argc, argv);
  parse_parameter_file(m, Gets("parameterfile"));
  parse_command_line_options(m);
  make_output_directory(m);
  initialize_libraries(m);

  iterate_parameters(m, 0); /* start of new iteration */
  while(iterate_parameters(m, 1))
  {
    RunFun(m, POST_PARAMETERS); //hook for funs right after iterate_parameters
    RunFun(m, INITMESH); // here we schedule funcs to programatically set up the mesh
    inidata_mesh(m);
    evolve_mesh(m);
    finalize_mesh(m);
    RunFun(m, POST_FINALIZE_MESH); //hook after finalize_mesh, e.g. for special cleanup
    makeparameter("outdir_previous_iteration", "", "outdir of previous iteration");
    Sets("outdir_previous_iteration", Gets("outdir"));
  }

  nmesh_MPI_Finalize();
  return 0;
}


/* read command line */
int read_command_line(int argc, char **argv)
{
  int i; 

  if (0) 
    for (i = 0; i < argc; i++)
      printf("argv[%d] = %s\n", i, argv[i]);

  prdivider(0);
  printf("Welcome to nmesh, compiled on %s at %s\n", __DATE__, __TIME__);
  prdivider(0);

  if (argc < 2)
  {
    printf("Usage:  nmesh name.par\n");
    printf("or:     nmesh name.par options and extra arguments\n");
    printf("or:     nmesh --argsfile args.txt\n");
    printf("\n");
    printf("options: --keep_previous           do not touch name_previous\n");
    printf("         --modify-par:\"P=v\"        set par P to value v\n");
    printf(" all options must start with --\n"); 
    exit(0);
  }

  /* got two parameters */
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
    printf("Adding command line parameters\n");
    AddPar("outdir", outdir, "output directory");
    AddPar("parameterfile", parfile, 
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
        AddPar(argi, argv[i], descr);
        nargs++;
      }
    }
    /* add nmesh command line options */
    if(nopts>0) AddPar("nmesh_options", options, "nmesh command line options");

    free(parfile);
    free(outdir);
  }

  /* more initialization */
  
  return 0;
}

int make_output_directory(void)
{
  char *outdir  = Gets("outdir");
  char *outdirp = (char *) calloc(strlen(outdir)+40, sizeof(char));
  char so[1000];

  /* set outdirp to outdir_previous */
  strcpy(outdirp, outdir);
  strcat(outdirp, "_previous");

  /* check if a shell is available to execute commands later */
  /* NOTE: system2 and system3 are smart enough to do "mkdir", "rm -rf"
           and "mv" even without a shell by using POSIX calls. */
  if(system(NULL)==0)
  {
    printf("WARNING: system(NULL)=0 => cannot execute shell commands!\n");
    printf("         Consider using system_emu.\n");
  }

  /* check if we remove outdir_previous */
  /*
  if(!GetvLax("nmesh_options", "--keep_previous"))
  {
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
    system2("rm -rf", outdirp);
    system3("mv", outdir, outdirp);
  }

  /* make output directory, save parfile */
  system2("mkdir", outdir);
  /* system3("cp", Gets("parameterfile"), outdir); */
  copy_file_into_dir(Gets("parameterfile"), outdir);

  /* redirect stdout and stderr for MPI jobs 
     all output is collected in outdir/stdout.001 etc  */
  if (nmesh_MPI_rank()>0)
  {
    f[100];
    snprintf(f,99, "%%s/stdout.%%0%dd", (int) log10(nmesh_MPI_size())+1);
    snprintf(so,999, f, outdir, nmesh_MPI_rank());  
    freopen(so, "w", stdout);   
    freopen(so, "w", stderr);
  }
  else /* always redirect output of proc0 to outdir.log */
  {
    snprintf(s,999, "%s.log", outdir);
    freopen(so, "a", stdout);   
    freopen(so, "a", stderr);
  }
  /* say what we have after redirection: */
  prdivider(0);
  PRF;printf(": stdout redirected to:\n %s\n", so);
  printf("nmesh was compiled on %s at %s\n", __DATE__, __TIME__);
  prdivider(0);

  free(outdirp);
  return 0;
}


/* go through options and act accordingly */
int parse_command_line_options(void)
{
  char *optionstr;
  char *str1;
  char *str2;
  char *par;
  char *val;
  
  /* get length of nmesh_options string */
  if(GetsLax("nmesh_options")==0) return 0;
  printf("Parsing command line options\n");

  /* parse for all --modify-par: */
  optionstr = strdup((GetsLax("nmesh_options")));
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


/* initialize mesh */
int inidata_mesh(tMesh *m)
{
  if (1) {
    prdivider(0);
    printf("Initializing mesh\n");
  }

  /* compute initial data */
  RunFun(m, INITIALDATA);

  /* initial data is just another new time slice */
  RunFun(m, POST_EVOLVE);

  /* initial data complete */
  prdivider(0);
  printf("Done with initialization\n");
  printf(" iteration %d, time=%g\n", m->iteration, m->time);

  /* analyze initial data */
  RunFun(m, ANALYZE);

  /* output for permanent variables */
  RunFun(m, OUTPUT);

  /* checkpoint, just in case we need it here already */
  //checkpoint(m);

  return 0;
}




/* evolve mesh */
int evolve_mesh(tMesh *mesh)
{
  int iterationmax = Geti("iterations");
  double timemax = Getd("finaltime");

  prdivider(0);

  if (timemax > 0)
    iterationmax = timemax/mesh->dt + 0.5;

  if (iterationmax > 0) 
    printf("Evolving mesh for %d top mesh iterations to time %.3f\n", 
	   iterationmax, iterationmax * mesh->dt);

  if (iterationmax <= 0) return 0;

  /* outermost evolution loop */
  while (mesh->iteration < iterationmax)
  { 
    /* pre evolve */
    RunFun(mesh, PRE_EVOLVE); 

    /* evolve */
    RunFun(mesh, EVOLVE); 

    /* post evolve */
    RunFun(mesh, POST_EVOLVE); 

    /* evolution step complete */
    mesh->iteration++;
    mesh->time = mesh->iteration * mesh->dt;

    /* print info */
    printf(" iteration %d, time=%g\n", mesh->iteration, mesh->time);
    fflush(stdout); 

    /* analyze */
    RunFun(mesh, ANALYZE);

    /* output for permanent variables */
    RunFun(mesh, OUTPUT);

    /* post output */
    RunFun(mesh, POST_OUTPUT);

    /* checkpoint */
    //checkpoint(mesh);
  
    /* update since this may change during evolution, say when checkpointing */
    timemax= Getd("finaltime");
    iterationmax= (timemax > 0) ? timemax/mesh->dt + 0.5 : Geti("iterations");
  }
  return 0;
}



/* finalize mesh */
int finalize_mesh(tMesh *m)
{
  prdivider(0);
  free_mesh(m);
  return 0;
}
