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
/* use main part of nmesh for these: */
void AddPar(char *parname, char *val, char *comm)
{
  //this is just an example for this function
  printf("  %s = %s\t\t/* %s */\n", parname, val, comm);
}

void AddEvoVar(char *varname, char *inds, char *comm)
{
  //this is just an example for this function
  printf("  AddEvoVar %s\"%s\"\t\t/* %s */\n", varname, inds, comm);
}
void AddAuxVar(char *varname, char *inds, char *comm)
{
  //this is just an example for this function
  printf("  AddAuxVar %s\"%s\"\t\t/* %s */\n", varname, inds, comm);
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
  //this is just an example for this function
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

/* useful macros to be added to main part of nmesh */
#define ADDPAR(module, type, name, value, comment) \
  AddPar(#module"_"#name, value, comment); \
  module->name = Par(#module"_"#name);

#define ADDVAR(module, type, name, indices, comment) \
  Add##type(#module"_"#name, indices, comment); \
  { char varname[100]; \
    varname_Comp0(#module"_"#name, indices, varname); \
    module->name = Ind(varname); }

#define INITPARINDEX(module, type, name, value, comment) \
  module->name = Par(#module"_"#name);
/***************************************************************************/

/***************************************************************************/
/* this goes into GHG_vars_and_pars.h */
#define GHG_VARS \
X(EvoVar, g, "ab",  "4-metric") \
X(AuxVar, d, "iab", "deriv of 4-metric") \

#define GHG_PARS \
X(Par, numflux,   "LLF", "num flux we use [LLF,upwind]") \
X(Par, gamma1,    "-1",  "constraint damping par that we use in the evo " \
"system [some number]") \
X(Par, dissfac,   "0",   "diss fac we use") \

/* struct with global vars and pars */
typedef struct {
  #define X(type,name,val,comment) int name; /* index of par or var called name */
  /* all vars */
  GHG_VARS
  /* all pars */
  GHG_PARS
  #undef X
  /* add custom fields here */
  //...
} tGHG;
/***************************************************************************/
// #include GHG_vars_and_pars.h in GHG.h or nmesh_GHG.h
/***************************************************************************/

/***************************************************************************/
/* This goes into nmesh_GHG.c */
//...

/* global vars */
tGHG GHG[1];

/* int nmesh_GHG(tMesh *mesh) */
int main()
{
  printf("Adding GHG\n");

  /* add vars from header file to GHG */
  #define X(type,name,indices,comment) ADDVAR(GHG, type,name,indices,comment)
  GHG_VARS
  #undef X
  /* add pars from header file to GHG */
  #define X(type,name,value,comment) ADDPAR(GHG, type,name,value,comment)
  GHG_PARS
  #undef X


  /* some print */
  printf("print var indices with X-macro:\n");
  #define X(type,name,indices,comment) \
    printf("  GHG->"#name" = %d\n", GHG->name);
  GHG_VARS
  #undef X

  /* some more print */
  printf("print par indices with X-macro:\n");
  #define X(type,name,value,comment) \
    printf("  GHG->"#name" = %d\n", GHG->name);
  GHG_PARS
  #undef X


  /* test use */
  printf("some globals are:\n");
  printf("  GHG->numflux = %d\n", GHG->numflux);
  printf("  GHG->g = %d\n", GHG->g);
}
