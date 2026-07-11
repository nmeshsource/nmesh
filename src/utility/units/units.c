/* units.c */
/* Wolfgang Tichy 6/2022 */

#include "nmesh.h"
#include "units.h"


/* set global unit struct */
int units_set_SItoGMc1_GMc1toSI(tMesh *mesh)
{
  PRFs(": set SItoGMc1 struct for c = G = G*Msun = 1.\n");
  units_set_SItoGMc1();
  units_test_tSItoGMc1();
  return 0;
}

/* print some unit info */
int print_unit_conversion_factors(tMesh *mesh)
{
  PRFs(":\n");
  printf("In nmesh we use units with c = G = G*Msun = 1.\n");
  printf("In SI units these constants are:\n");
  print_natureconsts();
  printf("Conversion factors to obtain SI units from nmesh units:\n");
  print_conversion_factors();
  return 0;
}
