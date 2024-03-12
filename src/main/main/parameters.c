/* parameters.c */
/* Wolfgang Tichy, 1/2019 & Bernd Bruegmann, 12/99 */

/* The parameter file has this format:

parameter_name = parameter_value1 parameter_value2 ...

   Note that if there are lots of values it can also be written like this
   in the parameter file:

parameter_name = parameter_value1 parameter_value2 parameter_value3
  parameter_value4 parameter_value5

   Comments:
   Anything after the character # is ignored until the end of the line.
   I.e. comments can look like these:

# my fav parameter:
parameter_name = parameter_value1  # this is just a comment

   For processing the parameters and their values, all comments are removed
   and then the file content is converted to the form
   #par=val ... val#par=val ... val ...
   Here # now indicates the start of a par and not a comment!
*/


#include "nmesh.h"
#include <ctype.h>  // needed for isspace


/* maximum number of pars allowed */
int npdbmax = 2000;

/* functions */
void makeparameter(tMesh *mesh, const char *name, const char *value,
                   const char *description);
int findparameterindex(tMesh *mesh, const char *name, int fatal);
tPar *findparameter(tMesh *mesh, const char *name, int fatal);
void setparameter(tMesh *mesh, int i, const char *value);
void printparameter(tPar *p);
void printparameters(tMesh *mesh);
void translatevalue(char **value);
int set_numericalvalue_byIndex(tPar *pdb1, int ind, int npdb1max);
int set_booleanvalue_byIndex(tPar *pdb1, int ind, int npdb1max);
int set_valuelen_byIndex(tPar *pdb1, int ind, int npdb1max);
int make_output_directory(tMesh *mesh);



/* open and parse a given parameter file */
int nmesh_load_parameters(tMesh *mesh, char *parfile, int fatal, int pr)
{
  FILE *fp;
  int c, i, j;
  int nbuffer, file_exists;
  char *buffer;
  char *par, *val;
  int lpar, lval;

  if(pr) printf("Reading parameter file \"%s\"\n", parfile);

  /* read file into memory, and also add one space at end and beginning */
  file_exists = 1;
  buffer = NULL;
  nbuffer = 0;
  if(Rank0_or_NoMPI)
  {
    int IObufsz = 1048576;
    char *IObuf; /* larger buffer for read */

    fp = fopen_buf(parfile, "r", &IObuf,IObufsz);
    if(!fp)
    {
      if(pr) printf("  parameter file \"%s\" does not exist!\n", parfile);
      if(fatal)
        errorexits("Could not open parameter file \"%s\"", parfile);
      else
        file_exists = 0; /* signals error reading file */
    }

    if(file_exists)
    {
      for(i = nbuffer = 0;; i++)
      {
        if(i >= nbuffer-2)
        {
          if(nbuffer > 1000000)
            errorexit("Parameter files bigger than 1MB are not allowed!");
          buffer = (char *) realloc(buffer, sizeof(char)*(nbuffer += 1000));
          if(!buffer)
            errorexit("Out of memory for buffer while reading parameter file!");
        }
        if(i == 0) buffer[i++] = ' ';
        if((c = fgetc(fp)) == EOF) break;
        buffer[i] = c;
      }
      fclose_buf(fp, &IObuf);
      buffer[i++] = ' ';
      buffer[i] = '\0';
      nbuffer = strlen(buffer);
    }
    else
    {
      nbuffer = 0;
    }
  }

  /* broadcast whether parfile exists */
  MCK( nMPI_Bcast(&file_exists,1, nMPI_INT, 0) );
  if(!file_exists) return 0;

  /* broadcast buffer to all MPI ranks */
  MCK( nMPI_Bcast(&nbuffer,1, nMPI_INT, 0) );
  if(!Rank0_or_NoMPI) buffer = calloc(nbuffer+2, sizeof(char));
  MCK( nMPI_Bcast(buffer,nbuffer+1, nMPI_CHAR, 0) );
  if(0) { printf("%s", buffer); Yo(1); }

  /* replace comments by spaces */
  for(i = 0; i < nbuffer; i++)
  {
    if(buffer[i] == '#')
      while (i < nbuffer && buffer[i] != '\n')
        buffer[i++] = ' ';
  }

  /* collapse all white space into single space */
  for(i = j = 1; i < nbuffer; i++) {
    if(!isspace(buffer[i])) buffer[j++] = buffer[i];
    else if(!isspace(buffer[j-1])) buffer[j++] = ' ';
  }
  buffer[j] = '\0';
  nbuffer = strlen(buffer);
  if(0) printf("|%s|\n", buffer);

  /* put a # in front of each par name */
  for(i=0; i<nbuffer-2; i++)
  {
    if(buffer[i]==' ') j=i; /* save pos of newest space found in j */
    if(buffer[i+1]=='=' || buffer[i+2]=='=') /* there's a '=' after pos i */
    {
      buffer[j]='#';  /* put a '#' in front of name, i.e. at pos j */
      i++;
    }
  }
  if(0) printf("|%s|\n", buffer);

  /* now remove spaces around = */
  for(i = j = 1; i < nbuffer; i++) {
    if(buffer[i] != ' ' || (buffer[i-1] != '=' && buffer[i+1] != '='))
      buffer[j++] = buffer[i];
  }
  buffer[j] = '\0';
  if(buffer[j-1] == ' ') buffer[j-1] = '\0';
  nbuffer = strlen(buffer);
  if(0) printf("|%s|\n", buffer);

  /* now the buffer is
     |#par=val ... val#par=val ... val| */

  /* split parameter names and values by replacing '#' and '=' by zero */
  for(i=0; i<nbuffer; i++)
  {
    if(buffer[i] == '#' || buffer[i] == '=')
      buffer[i] = '\0';
  }

  /* loop over all parameter/value pairs */
  for(i = 1; i < nbuffer; i += lpar + lval + 2)
  {
    int pari;
    par = buffer+i;
    lpar = strlen(par);
    val = par + lpar + 1;
    lval = strlen(val);
    if(0) printf("%s = |%s|\n", par, val);

    pari = findparameterindex(mesh, par, 0);
    if(pari<0)
    {
      /* errorexits("%s in parameterfile is not yet in mesh->pdb", par); */
      makeparameter(mesh, par, val, "parameter found only in parfile");
    }
    else
    {
      setparameter(mesh, pari, val);
    }
  }

  /* print parameters */
  if(0)
  {
    printf("parameters after reading the parameterfile:\n");
    printparameters(mesh);
  }
  free(buffer);
  return 1;
}

/* parse a parameter file */
void parse_parameter_file(tMesh *mesh, char *parfile)
{
  nmesh_load_parameters(mesh, parfile, 1, Rank0);
}

/* update parameters from file outdir/nmesh_update_parameters.par */
int nmesh_update_parameters(tMesh *mesh)
{
  static double last_update_time = 0.;
  double hours = Getd(Par("update_parameters_hours"));
  double time  = getTimeIn_s()/3600.;
  double time_since_update;
  int do_update = 0;

  /* test if it is time */
  time_since_update = time - last_update_time;
  if(Rank0)
  {
    /* test based on walltime */
    if(hours >= 0. && hours <= time_since_update)
      do_update = 1; /* yes, we want to update pars */
  }

  /* broadcast do_update from rank0 to all others */
  MCK( nMPI_Bcast(&do_update, 1, nMPI_INT, 0) );

  /* now do it if needed */
  if(do_update)
  {
    char update_file[] = "nmesh_update_parameters.par";
    char *outdir = Gets(Par("outdir"));
    int pl = strlen(outdir) + 80;
    char *pars = cmalloc(pl);
    int parsread;

    /* read pars in file nmesh_update_parameters.par */
    snprintf(pars,pl, "%s/%s", outdir, update_file);
    parsread = nmesh_load_parameters(mesh, pars, 0, 0);
    if(parsread)
    {
      PRF;printf(":\n read \"%s\"\n", pars);
      if(Rank0)
      {
        char *pars2 = cmalloc(pl);
        /* rename file, so that we do not update over and over again */
        snprintf(pars2,pl, "%s/%s%s", outdir, update_file,"_done");
        rename(pars, pars2);
        free(pars2);
      }
    }

    /* update times */
    last_update_time = time;
    free(pars);
  }
  return 0;
}

/* save all current pars in file fname */
void nmesh_save_parameters(tMesh *mesh, const char *fname)
{
  FILE *fp;
  int i;

  fp = fopen(fname, "w");
  if(!fp) errorexits("could not open %s", fname);

  for(i=0; i<mesh->npdb; i++)
  {
    char *name = MeshParGetName(mesh, i);
    char *val  = Gets(i);

    fprintf(fp, "%s = %s\n", name, val);
  }

  fclose(fp);
}


/***************************************************************************/
/* parameter data base */

/* make new parameter in parameter data base, merge if already there */
void makeparameter(tMesh *mesh, const char *name, const char *value,
                   const char *description)
{
  tPar *p;

  if(0) {PRF;printf(" mesh=%p\n", (void *) mesh);}
  if(0) {PRF;printf(" %s = %s,  %s\n", name, value, description);}

  /* check if name is banned (i.e. blacklisted) */
  ExitIfParBanned(name);

  mesh->pdb = realloc(mesh->pdb, npdbmax*sizeof(tPar));
  if(!mesh->pdb) errorexit("out of memory for mesh->pdb");

  p = findparameter(mesh, name, 0);
  if(!p)
  {
    int i = mesh->npdb;
    p = &(mesh->pdb[i]);
    mesh->npdb++;
    p->name  = strdup(name);
    p->value = NULL;
    setparameter(mesh, i, value);
  }
  else
  {
    free(p->description);
  }
  p->description = (char *) calloc(strlen(description)+1, sizeof(char));
  strcpy(p->description, description);

  if(mesh->npdb >= npdbmax)
    errorexit("no space for more parameters because npdbmax is hard coded");

  if(0) printparameters(mesh);
}

/* free strings in vdb */
void free_mesh_pdb_contents(tMesh *mesh)
{
  tPar *par;
  int i;

  for(i = 0; i < mesh->npdb; i++)
  {
    par = &(mesh->pdb[i]);
    free(par->name);
    free(par->value);
    free(par->description);
  }
}



/* set parameter */
void setparameter(tMesh *mesh, int i, const char *value)
{
  tPar *p;

  if(i<0 || i>=mesh->npdb) errorexit("this parameter index does not exist");

  p = &(mesh->pdb[i]);
  free(p->value);
  p->value = strdup(value);
  translatevalue(&p->value);
  set_numericalvalue_byIndex(p, 0, 1);
  set_booleanvalue_byIndex(p, 0, 1);
  set_valuelen_byIndex(p, 0, 1);
}


/* set the mesh pdb_iStart to the index of par name */
void Set_pdb_iStart(tMesh *mesh, int i)
{
  if(i<0 || i>=mesh->npdb) errorexit("this index does not exist");
  mesh->pdb_iStart = i;
}


/* find parameter index */
int findparameterindex(tMesh *mesh, const char *name, int fatal)
{
  tPar *pdb = mesh->pdb;
  int npdb = mesh->npdb;
  int iS   = mesh->pdb_iStart;
  int i;

  if(!name) errorexit("no parameter name");

  if( (iS < 0) || (iS >= npdb) ) iS = 0; /* make sure first i is in range */

  for(i = iS; i < npdb; i++)
    if(!strcmp(pdb[i].name, name))
      return i;

  for(i = 0; i < iS; i++)
    if(!strcmp(pdb[i].name, name))
      return i;

  if(fatal) errorexits("Could not find parameter \"%s\"", name);

  /* -1 means par was not found */
  return -1;
}
tPar *parameterfromindex(tMesh *mesh, int i)
{
  if(i<0 || i>=mesh->npdb)
    errorexit("parameter with this index does not exist");

  return &(mesh->pdb[i]);
}
/* find parameter */
tPar *findparameter(tMesh *mesh, const char *name, int fatal)
{
  int i = findparameterindex(mesh, name, fatal);

  if(i<0) return 0;

  return &(mesh->pdb[i]);
}


/* translate parameter value in some simple cases */
void translatevalue(char **value)
{
  double x = 0;

  if(strcmp(*value, "pi")    == 0) x =    PI;
  if(strcmp(*value, "-pi")   == 0) x =   -PI;
  if(strcmp(*value, "pi/2")  == 0) x =    PI/2;
  if(strcmp(*value, "-pi/2") == 0) x =   -PI/2;
  if(strcmp(*value, "2*pi")  == 0) x =  2*PI;
  if(strcmp(*value, "2pi")   == 0) x =  2*PI;
  if(strcmp(*value, "-2*pi") == 0) x = -2*PI;
  if(strcmp(*value, "-2pi")  == 0) x = -2*PI;

  if(strcmp(*value,  "DBL_MAX")  == 0) x =  DBL_MAX;
  if(strcmp(*value, "-DBL_MAX")  == 0) x = -DBL_MAX;

  if(x)
  {
    char newvalue[100];
    snprintf(newvalue,99, "%.18e", x);
    free(*value);
    *value = strdup(newvalue);
  }
}

/* Write the numerical (double) value of the par with index ind into the
   par cache.
   Note: we keep the numerical (double) values of each par in a cache */
int set_numericalvalue_byIndex(tPar *pdb1, int ind, int npdb1max)
{
  if(pdb1!=NULL && ind>=0 && ind<npdb1max)
    pdb1[ind].numericalvalue = atof(pdb1[ind].value);
  else
    errorexit("parameter index out of range");

  return 1;
}

/* Write the boolean (int) value of the par with index ind into the
   par cache.
   Note: we keep the boolean values of each par in a cache */
int set_booleanvalue_byIndex(tPar *pdb1, int ind, int npdb1max)
{
  if(pdb1!=NULL && ind>=0 && ind<npdb1max)
  {
    char *par = pdb1[ind].value;
    int boolval=atoi(par); /* try to read integer value */

    if(strstr(par, "yes")) boolval=1;
    if(strstr(par, "Yes")) boolval=1;
    if(strstr(par, "YES")) boolval=1;
    if(strstr(par, "on")) boolval=1;
    if(strstr(par, "On")) boolval=1;
    if(strstr(par, "ON")) boolval=1;
    if(strstr(par, "true")) boolval=1;
    if(strstr(par, "True")) boolval=1;
    if(strstr(par, "TRUE")) boolval=1;

    pdb1[ind].booleanvalue = boolval;
  }
  else
    errorexit("parameter index out of range");

  return 1;
}

/* Write the length of value of the par with index ind into the
   par cache. */
int set_valuelen_byIndex(tPar *pdb1, int ind, int npdb1max)
{
  if(pdb1!=NULL && ind>=0 && ind<npdb1max)
    pdb1[ind].valuelen = strlen(pdb1[ind].value);
  else
    errorexit("parameter index out of range");

  return 1;
}

void printparameter(tPar *p)
{
  printf("%16s = %-16s,  %s\n                 = %g  ->  %d\n",
  p->name, p->value, p->description,  p->numericalvalue, p->booleanvalue);
}

/* print parameters */
void printparameters(tMesh *mesh)
{
  tPar *pdb = mesh->pdb;
  int npdb = mesh->npdb;
  int i;

  for(i = 0; i < npdb; i++)
    printf("pdb[%2d]:  %12s = %-16s,  %s\n",
           i, pdb[i].name, pdb[i].value, pdb[i].description);
}

/* check all parameters */
int CheckForBannedPars(tMesh *mesh)
{
  tPar *pdb = mesh->pdb;
  int npdb = mesh->npdb;
  int i;

  for(i = 0; i < npdb; i++) ExitIfParBanned(pdb[i].name);

  return 0;
}


/***************************************************************************/
/* functions for external calls, we have macros for many of them, e.g.
   AddPar(name, value, description);
   instead of AddMeshPar(mesh, name, value, description); */

/* creation functions */
void AddMeshPar(tMesh *mesh, const char *name, const char *value,
                const char *description)
{
  int i;

  makeparameter(mesh, name, value, description);
  i = Par(name);
  printf("  par_%04d  %-25s  =  %s\n", i, name, Gets(i));
}

void AddOrModifyMeshPar(tMesh *mesh, const char *name, const char *value,
                        const char *description)
{
  int pari = findparameterindex(mesh, name, 0);

  if(pari<0)
    makeparameter(mesh, name, value, description);
  else
    setparameter(mesh, pari, value);
  if(Rank0) printf("  parameter %-25s  =  %s\n", name, Gets(Par(name)));
}

/* functions for setting pars */
void MeshParSets(tMesh *mesh, int pi, const char *value)
{
  setparameter(mesh, pi, value);
}

void MeshParSeti(tMesh *mesh, int pi, int i)
{
  char value[100];
  snprintf(value,99, "%d", i);
  setparameter(mesh, pi, value);
}

void MeshParSetd(tMesh *mesh, int pi, double d)
{
  char value[100];
  snprintf(value,99, "%.20e", d);
  setparameter(mesh, pi, value);
}

void MeshParAppends(tMesh *mesh, int pi, const char *value)
{
  if(MeshParGetv_fatal(mesh, pi, value, 1)) return;
  {
    char *oldvalue = MeshParGets(mesh, pi);
    char *newvalue = cmalloc(strlen(oldvalue) + strlen(value) + 2);
    if(oldvalue[0]) sprintf(newvalue, "%s %s", oldvalue, value);
    else            sprintf(newvalue, "%s", value);
    setparameter(mesh, pi, newvalue);
    free(newvalue);
  }
}


/* functions for getting par values */
char *MeshParGets(tMesh *mesh, int i)
{
  if(i<0 || i>=mesh->npdb)
    errorexit("parameter with this index does not exist");
  return mesh->pdb[i].value;
}

char *MeshParGetsLax(tMesh *mesh, int i)
{
  if(i<0 || i>=mesh->npdb) return 0;
  return mesh->pdb[i].value;
}

int MeshParGeti(tMesh *mesh, int i)
{
  if(i<0 || i>=mesh->npdb)
    errorexit("parameter with this index does not exist");
  return mesh->pdb[i].numericalvalue;
}

double MeshParGetd(tMesh *mesh, int i)
{
  if(i<0 || i>=mesh->npdb)
    errorexit("parameter with this index does not exist");
  return mesh->pdb[i].numericalvalue;
}

/* return 0 or 1 */
int MeshParGetb(tMesh *mesh, int i)
{
  if(i<0 || i>=mesh->npdb)
    errorexit("parameter with this index does not exist");
  return mesh->pdb[i].booleanvalue;
}

/* return length of par value string */
int MeshParGetLen(tMesh *mesh, int i)
{
  if(i<0 || i>=mesh->npdb)
    errorexit("parameter with this index does not exist");
  return mesh->pdb[i].valuelen;
}

/* "get value?" returns 1 if value is in the list of values and 0 else
   (not equivalent to value being a substring of string parameter) */
int MeshParGetv_fatal(tMesh *mesh, int i, const char *value, int fatal)
{
  tPar *p;
  char *s=NULL;
  char *parval;
  int lv, ls, lp, startok, endok;

  if(i<0 || i>=mesh->npdb)
  {
    if(fatal) errorexit("parameter with this index does not exist");
    else      return 0;
  }

  p = &(mesh->pdb[i]);

  if(!p) return 0;
  parval = p->value;
  while( (s = strstr(parval, value)) )
  {
    lv = strlen(value);
    //printf("ls=%d  value=%s| s=%s| p->value=%s|\n", ls, value, s, p->value);
    if( s[lv]==' ' || s[lv]==0 ) break;
    parval = s+1;
  }
  if(!s) return 0;
  ls = strlen(s);
  lp = strlen(p->value);
  startok = (s == p->value || *(s-1) == ' ');  /* how robust is this? */
  endok   = (s+ls == p->value+lp || *(s+ls) == ' ');
  return startok && endok ? 1 : 0;
}

/* get name of par with index i */
char *MeshParGetName(tMesh *mesh, int i)
{
  if(i<0 || i>=mesh->npdb)
    errorexit("parameter with this index does not exist");
  return mesh->pdb[i].name;
}



/**************************************************************************/
/* we can iterate over pars */
int iterate_parameters(tMesh* mesh, int next)
{
  static int iter = 0; // should make this a value in mesh or a global var
  tPar *p;
  char *list, *name, *newvalue=NULL, *value, *saveptr;
  char iterpar[100] = "iterate_parameter1";
  char newoutdir[10000];
  int i, j, l;

  /* reset iter to zero if next=0 */
  if(next==0) { iter=0;  return 0; }

  /* the default is that we don't want to iterate */
  if(!Getv(Par("iterate_parameters"), "yes"))
  {
    /* return 1 for first call, but 0 for second call, which exits nmesh */
    if(iter < 0)
      return 0;
    iter = -1;
    return 1;
  }

  /* If we get here, we want to iterate */
  if(iter<0) iter=0;
  printf("\n");
  prdivider(0);
  printf("Iterating parameters:\n");

  sprintf(newoutdir, "%s", Gets(Par("parameterfile")));
  *strstr(newoutdir, ".par") = '\0';

  p = findparameter(mesh, iterpar, 0);
  if(!p)
    errorexit("nothing to iterate, specify at least iterate_parameter1\n");

  j = 1;
  while(p)
  {
    /* put val of iterate_parameter? in list,
       duplicate Gets(Par(iterpar)) because strtok_r will modify list */
    list = strdup( Gets(Par(iterpar)) );
    name = strtok_r(list, " ", &saveptr);

    if(!findparameter(mesh, name, 0))
      errorexit("iterate_parameterN has to start with name of "
                "existing parameter");

    value = strtok_r(0, " ", &saveptr);
    if(!value)
      errorexit("iterate_parameterN needs at least one value");

    /* loop until newvalue is value number iter */
    for(i=0, newvalue=0; value; i++)
    {
      if(i == iter) newvalue = value;
      value = strtok_r(0, " ", &saveptr);
    }

    /* if we get a new value we set the par and make a new output dir */
    if(newvalue)
    {
      Sets(Par(name), newvalue);
      printf("%38s = %s\n", name, Gets(Par(name)));

      l = strlen(newoutdir);
      sprintf(newoutdir+l, "_%s", newvalue);
      Sets(Par("outdir"), newoutdir);
    }

    /* get the next iterate_parameterN */
    j++;
    l = strlen("iterate_parameter");
    sprintf(iterpar+l, "%d", j);
    p = findparameter(mesh, iterpar, 0);

    free(list);
  }

  //printparameters(mesh);
  if(newvalue)
  {
    printf("Starting parameter iteration %d:\n", iter);
    if(1) printf("outdir = %s\n", Gets(Par("outdir")));
    make_output_directory(mesh);
    iter++;
    return 1;
  }
  printf("  finished iterating.\n");
  prdivider(0);
  return 0;
}

/* print pars with index i1 to i2 in pdb */
void print_pdb_i1_i2(tPar *pdb, int i1, int i2, int pr_ind, int pr_cache)
{
  int i;

  printf("print_pdb_i1_i2: pdb=%p i1=%d i2=%d\n", (void *) pdb, i1,i2);
  for(i=i1; i<=i2; i++)
  {
    if(pr_ind) printf("#%d# ", i);
    if(pr_cache) printf("|%g|%d| ", pdb[i].numericalvalue, pdb[i].booleanvalue);
    printf("%s:\n", pdb[i].description);
    printf("%s = %s\n", pdb[i].name, pdb[i].value);
  }
}
/* print entire parameter database */
void print_parameter_database(tMesh *mesh)
{
  print_pdb_i1_i2(mesh->pdb, 0, mesh->npdb-1, 0,0);
}

/**********************************************/
/* functions to copy the entire par data base */
/**********************************************/

/* Create a copy of parameter database pdb1 in pdb2.
   This allocates all memory needed for pdb2.
   The caller has to free pdb2 later on its own, e.g. with free_pdb */
/* if we have: tPar *pdb2;
   call like:  create_copy_of_pdb1_in_pdb2(pdb,npdb,npdbmax, &pdb2); */
void create_copy_of_pdb1_in_pdb2(tPar *pdb1, int npdb1, int npdb1max,
                                 tPar **pdb2)
{
  int i;
  tPar *p2;

  /* allocate array for p2 */
  p2 = (tPar *) calloc(npdb1max, sizeof(tPar));
  if(!p2) errorexit("create_copy_of_pdb1_in_pdb2: out of memory");

  /* go over pars in pdb1 and use strdup to create copies in p2 */
  for(i=0; i<npdb1; i++)
  {
    p2[i].name  = strdup(pdb1[i].name);
    p2[i].value = strdup(pdb1[i].value);
    p2[i].description = strdup(pdb1[i].description);
    p2[i].numericalvalue = pdb1[i].numericalvalue;
    p2[i].booleanvalue = pdb1[i].booleanvalue;
    if(!p2[i].name || !p2[i].value || !p2[i].description)
      errorexit("create_copy_of_pdb1_in_pdb2: out of memory");
  }
  /* set pdb2 */
  *pdb2=p2;
}

/* Make an empty par database of max length npdb1max.
   Usage: pdb2 = make_empty_pdb(npdbmax); */
tPar *make_empty_pdb(int npdb1max)
{
  int i;
  tPar *pdb1;
  /* allocate array for pdb1 */
  pdb1 = (tPar *) calloc(npdb1max, sizeof(tPar));
  if(!pdb1) errorexit("make_empty_pdb: out of memory");

  /* set all entries to NULL */
  for(i=0; i<npdb1max; i++)
  {
    pdb1[i].name  = NULL;
    pdb1[i].value = NULL;
    pdb1[i].description = NULL;
    pdb1[i].numericalvalue = 0.0;
    pdb1[i].booleanvalue = 0;
  }

  return pdb1;
}

/* copy pdb1 into pdb2 */
void copy_pdb(tPar *pdb1, int npdb1, tPar *pdb2)
{
  int i;

  /* Go over pars in pdb1. Use realloc and strcpy to create copies in pdb2 */
  for(i=0; i<npdb1; i++)
  {
    pdb2[i].name =
     (char *) realloc(pdb2[i].name,
                      sizeof(char)*(strlen(pdb1[i].name)+1));
    pdb2[i].value =
     (char *) realloc(pdb2[i].value,
                      sizeof(char)*(strlen(pdb1[i].value)+1));
    pdb2[i].description =
     (char *) realloc(pdb2[i].description,
                      sizeof(char)*(strlen(pdb1[i].description)+1));
    if(!pdb2[i].name || !pdb2[i].value || !pdb2[i].description)
      errorexit("copy_pdb: out of memory");
    strcpy(pdb2[i].name, pdb1[i].name);
    strcpy(pdb2[i].value, pdb1[i].value);
    strcpy(pdb2[i].description, pdb1[i].description);
    pdb2[i].numericalvalue = pdb1[i].numericalvalue;
    pdb2[i].booleanvalue = pdb1[i].booleanvalue;
  }
}

/* free the parameter database content in pdb1 */
void free_pdb_contents(tPar *pdb1, int npdb1)
{
  int i;

  /* go over pars in pdb1 and free them */
  for(i=0; i<npdb1; i++)
  {
    free(pdb1[i].name);
    free(pdb1[i].value);
    free(pdb1[i].description);
  }
}

/* free the parameter database pdb1 */
void free_pdb(tPar *pdb1, int npdb1)
{
  free_pdb_contents(pdb1, npdb1);
  free(pdb1);
}
