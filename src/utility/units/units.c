/* units.c */
/* Wolfgang Tichy 6/2022 */

#include "nmesh.h"
#include "units.h"


/* set global unit struct */
int unit_set_for_GMc1_mesh(tMesh *mesh)
{
  PRFs(": set unit struct for c = G = G*Msun = 1.\n");
  unit_set_for_GMc1();
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
