/*
Below is an example on how we could use X-macros to replace
Par("GHG_parname") and Ind("GHG_varcompname)" by
GHG->parname and GHG->varname.
In this way we could avoid string searches altogether!!!
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/***************************************************************************/
/* put into other files: */
void AddPar(char *parname, char *val, char *comm)
{
  //this is just an example for this function
  printf("  %s = %s\t\t/* %s */\n", parname, val, comm);
}

void AddVar(char *varname, char *inds, char *comm)
{
  //this is just an example for this function
  printf("  add %s\"%s\"\t\t/* %s */\n", varname, inds, comm);
}

int Par(char *parname)
{
  //this is just an example for this function
  return parname[4];
}

int Ind(char *varname)
{
  //this is just an example for this function
  return strlen(varname);
}
void varname_Comp0(const char *name, const char *indices, char *varname0)
{
  #define NINDEXLIST 100
  //int nilist;
  char *ilist[NINDEXLIST];
  //int sym[3*NINDEXLIST];
  /* use
       tensorindexlist(indices, &nilist, ilist, sym);
     to find ilist[0] */
  ilist[0] = malloc(666);
  strcpy(ilist[0], indices);
  sprintf(varname0, "%s%s", name, ilist[0]);
  free(ilist[0]);
}

/* useful macros */
#define ADDPAR(module, name, value, comment) \
  AddPar(#module"_"#name, value, comment); \
  module->name = Par(#module"_"#name);

#define ADDVAR(module, name, indices, comment) \
  AddVar(#module"_"#name, indices, comment); \
  { char varname[100]; \
    varname_Comp0(#module"_"#name, indices, varname); \
    module->name = Ind(varname); }

#define INITPARINDEX(module, name, value, comment) \
  module->name = Par(#module"_"#name);
/***************************************************************************/

/***************************************************************************/
/* put into GHG header file GHG_vars_and_pars.h */
#define GHG_VARS \
X (g, "ab",  "4-metric") \
X (d, "iab", "deriv of 4-metric") \

#define GHG_PARS \
X (numflux,   "LLF", "num flux we use [LLF,upwind]") \
X (gamma1,    "-1",  "constraint damping par that we use in the evo system " \
                     "[some number]") \
X (dissfac,   "0",   "diss fac we use") \

/* struct with global vars and pars */
typedef struct {
  #define X(name,val,comment) int name; /* index of par or var called name */
  /* all vars */
  GHG_VARS
  /* all pars */
  GHG_PARS
  #undef X
  /* add custom fields here */
  //...
} tGHG;
/***************************************************************************/
//#include GHG header file GHG_vars_and_pars.h in GHG.h or nmesh_GHG.h
/***************************************************************************/

/***************************************************************************/
//...

/* global vars */
tGHG GHG[1];

/* int nmesh_GHG(tMesh *mesh) */
int main()
{
  /* add vars from header file to GHG */
  #define X(name,indices,comment) ADDVAR(GHG, name,indices,comment)
  GHG_VARS
  #undef X
  /* add pars from header file to GHG */
  #define X(name,value,comment) ADDPAR(GHG, name,value,comment)
  GHG_PARS
  #undef X


  /* some print */
  printf("print var indices with X-macro:\n");
  #define X(name,indices,comment) \
    printf("GHG->"#name" = %d\n", GHG->name);
  GHG_VARS
  #undef X

  /* some more print */
  printf("print par indices with X-macro:\n");
  #define X(name,value,comment) \
    printf("GHG->"#name" = %d\n", GHG->name);
  GHG_PARS
  #undef X


  /* test use */
  printf("some globals are:\n");
  printf("GHG->numflux = %d\n", GHG->numflux);
  printf("GHG->g = %d\n", GHG->g);
}
