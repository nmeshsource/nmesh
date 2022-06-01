/* convert_units.c */
/* Wolfgang Tichy 6/2022 */

#include <stdio.h>

/***************************************************************************/
/* struct with constans of nature */
/***************************************************************************/

/* constans of nature*/
struct tNATURECONSTS
{
  const double c;
  const double G;
  const double GMsun;
};

/* values of consts of nature we use */
struct tNATURECONSTS natureconsts =
{
  .c     = 299792458,        // [m/s],     c is defined to have this value
  .G     = 6.67430e-11,      // [m^3/(kg s^2)], from https://physics.nist.gov/cgi-bin/cuu/Value?bg
  .GMsun = 1.32712440018e20, // [m^3/s^2], standard gravitational parameter, measured using Kepler's law
};


/***************************************************************************/
/* convert from units with G=c=Msun=1 to SI */
/***************************************************************************/

/* Length in SI units */
double Length_GMc1_to_SI(double x)
{
  double c     = natureconsts.c;
  double GMsun = natureconsts.GMsun;
  return x * GMsun/(c*c);
}

/* Time in SI units */
double Time_GMc1_to_SI(double t)
{
  double c     = natureconsts.c;
  double GMsun = natureconsts.GMsun;
  return t * GMsun/(c*c*c);
}

/* frequency in SI units */
double Frequency_GMc1_to_SI(double f)
{
  return 1. / Time_GMc1_to_SI(1./f);
}

/* speed in SI units */
double Speed_GMc1_to_SI(double v)
{
  return  v * natureconsts.c;
}

/* acceleration in SI units */
double Acceleration_GMc1_to_SI(double a)
{
  double L = Length_GMc1_to_SI(1.);
  double T = Time_GMc1_to_SI(1.);
  return a * L/(T*T);
}

/* Mass in SI units */
double Mass_GMc1_to_SI(double m)
{
  double G     = natureconsts.G;
  double GMsun = natureconsts.GMsun;
  return m * GMsun/G;
}

/* mass density in SI */
double MassDensity_GMc1_to_SI(double rho)
{
  double L = Length_GMc1_to_SI(1.);
  double L3 = L*L*L;
  double Msun = Mass_GMc1_to_SI(1.);
  return rho * Msun/L3;
}

/* energy density in SI */
double EnergyDensity_GMc1_to_SI(double rho)
{
  double c = natureconsts.c;
  return MassDensity_GMc1_to_SI(rho) * c*c;
}

/* pressure density in SI */
double Pressure_GMc1_to_SI(double P)
{
  double L = Length_GMc1_to_SI(1.);
  double a = Acceleration_GMc1_to_SI(1.);
  double Msun = Mass_GMc1_to_SI(1.);
  return P * Msun*a/(L*L);
}


/***************************************************************************/
/* print some info */
/***************************************************************************/
void print_natureconsts(void)
{
  printf("  c = %.15g m/s\n", natureconsts.c);
  printf("  G = %.15g m^3/(kg s^2)\n", natureconsts.G);
  printf("  GMsun = %.15g m^3/s^2\n", natureconsts.GMsun);
}

void print_conversion_factors(void)
{
  double L = Length_GMc1_to_SI(1.);
  double T = Time_GMc1_to_SI(1.);
  double c = Speed_GMc1_to_SI(1.);
  double f = Frequency_GMc1_to_SI(1.);
  double Msun = Mass_GMc1_to_SI(1.);
  double rhoM = MassDensity_GMc1_to_SI(1.);
  double rhoE = EnergyDensity_GMc1_to_SI(1.);
  double Pc   = Pressure_GMc1_to_SI(1.);

  printf("  L = %.15g m\n", L);
  printf("  T = %.15g s\n", T);
  printf("  c = %.15g m/s\n", c);
  printf("  f = %.15g Hz\n", f);
  printf("  Msun = %.15g kg\n", Msun);
  printf("  rhoM = %.15g kg/m^3\n", rhoM);
  printf("  rhoE = %.15g J/m^3\n", rhoE);
  printf("  Pc = %.15g N/m^2\n", Pc);
}
