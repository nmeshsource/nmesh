/* variables.h */
/* Wolfgang Tichy, 5/2019 */

/* for variable types used in variables.c */

enum
{
  EVOVAR=0,  //loadbal transfer via MPI, interp on refine, exchange surfaces
  AUXVAR=1,  //no loadbal transfer, no interp on refine, no surfaces
  DATAVAR=2, //loadbal transfer via MPI and interp on refine, but no surfaces
  LBTVAR=4   //loadbal transfer via MPI, but no interp. and no surfaces
};
//Note: we may use 0,1 instead of EVOVAR,AUXVAR somwhere
