/* variables.h */
/* Wolfgang Tichy, 5/2019 */

/* for variable types used in variables.c */
enum
{
  EVOVAR=0,  //legacy warning: we may use 0 instead of EVOVAR somwhere
  AUXVAR=1,  //legacy warning: we may use 1 instead of AUXVAR somwhere
  DATAVAR=2,
  EXCHVAR=4
};
