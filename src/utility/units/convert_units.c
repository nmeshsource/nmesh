/* convert_units.c */
/* Wolfgang Tichy 6/2022 */

#include <stdio.h>
#include "convert_units.h"



/***************************************************************************/
/* init struct with constants of nature */
/***************************************************************************/

/* SI values of consts of nature we use */
tNatureconsts natureconsts =
{
  .c     = 299792458,     // [m/s],          c is defined to have this value
  .G     = 6.67430e-11,   // [m^3/(kg s^2)], from https://physics.nist.gov/cgi-bin/cuu/Value?bg
  .GMsun = 1.32712440041279419e20, // [m^3/s^2], standard gravitational parameter
  /* Note: GMsun can be very precisely measured using Kepler's law.
     Also, any Msun listed anywhere is computed from:  Msun = GMsun/G.
     Unfortunately G itself is not known as precisely.
     Some older sources give:
   GMsun = 1.32712440018e20 \pm 8e9 m^3/s^2
     But https://ssd.jpl.nasa.gov/astro_par.html says:
   GMsun = 1.32712440041279419e20 m^3/s^2
     This page also cites:
      Park, R.S., et al., 2021, "The JPL Planetary and Lunar Ephemerides
      DE440 and DE441", Astronomical Journal, 161:105.
      https://ssd.jpl.nasa.gov/doc/Park.2021.AJ.DE440.pdf
     In Tab. 2 it has:
   GMsun = 1.32712440041279419e20 m^3/s^2 (estimated from DE440)
     There is also:
      E. V. Pitjeva, "Determination of the Value of the Heliocentric
      Gravitational Constant (GMsun) from Modern Observations of Planets and
      Spacecraft", Journal of Physical and Chemical Reference Data 44,
      031210 (2015); https://doi.org/10.1063/1.4921980
     In its abstract it says:
   GMsun = 132 712 440 042 \pm 10 (km^3/s^2)
     Within the error this is the same as the 1.32712440041279419e20 m^3/s^2
     we use here.
     On https://ssd.jpl.nasa.gov/planets/phys_par.html it lists only planets.
  */
  /* From: CODATA values of the fundamental physical constants 2022
     ( see https://physics.nist.gov/cuu/Constants ): */
  .h     = 6.62607015e-34, // [J s],         h is defined to have this value
  .k_B   = 1.380649e-23,   // [J/K],       k_B is defined to have this value
  .e     = 1.602176634e-19,// [C],           e is defined to have this value
  .m_e   = 9.1093837139e-31,  // [kg],     measured, last 2 digits +-28
  .m_u   = 1.66053906892e-27, // [kg],     measured, last 2 digits +-52
  .m_n   = 1.67492750056e-27, // [kg],     measured, last 2 digits +-85
  .m_p   = 1.67262192595e-27, // [kg],     measured, last 2 digits +-52
};


/***************************************************************************/
/* convert from SI to nmesh units */
/***************************************************************************/

/* global struct that contains the values of various SI units */
tSItoGMc1 SItoGMc1;

/* set the values of SI units for G=Msun=c=1 */
void units_set_SItoGMc1(void)
{
  double c     = natureconsts.c;
  double G     = natureconsts.G;
  double GMsun = natureconsts.GMsun;
  double Msun  = GMsun/G;
  double e     = natureconsts.e;

  SItoGMc1.m  = (c*c)/GMsun;
  SItoGMc1.m2 = SItoGMc1.m * SItoGMc1.m;
  SItoGMc1.m3 = SItoGMc1.m * SItoGMc1.m2;
  SItoGMc1.s  = (c*c*c)/GMsun;
  SItoGMc1.s2 = SItoGMc1.s * SItoGMc1.s;
  SItoGMc1.kg = 1./Msun;
  SItoGMc1.J  = 1./(Msun*c*c);
  SItoGMc1.W  = SItoGMc1.J / SItoGMc1.s;
  SItoGMc1.N  = SItoGMc1.J / SItoGMc1.m;
  SItoGMc1.Pa = SItoGMc1.N / SItoGMc1.m2; // = SItoGMc1.J / SItoGMc1.m3;
  SItoGMc1.MeV= SItoGMc1.J * e * 1e6;
}


/***************************************************************************/
/* convert from units with G=c=Msun=1 to SI */
/***************************************************************************/

/* global struct that contains the values of various SI units */
tGMc1toSI GMc1toSI;

/* set the values of G=Msun=c=1 units in SI units */
void units_set_GMc1toSI(void)
{
  double c     = natureconsts.c;
  double c2    = c*c;
  double G     = natureconsts.G;
  double GMsun = natureconsts.GMsun;
  double Msun  = GMsun/G;
  double e     = natureconsts.e;

  GMc1toSI.Length = GMsun / c2;
  GMc1toSI.Area   = GMc1toSI.Length * GMc1toSI.Length;
  GMc1toSI.Volume = GMc1toSI.Area   * GMc1toSI.Length;
  GMc1toSI.Time   = GMc1toSI.Length / c;
  GMc1toSI.Mass   = Msun;
  GMc1toSI.Energy = GMc1toSI.Mass *c2;
  GMc1toSI.Power  = GMc1toSI.Energy / GMc1toSI.Time;
  GMc1toSI.Force  = GMc1toSI.Energy / GMc1toSI.Length;
  GMc1toSI.Press  = GMc1toSI.Force / GMc1toSI.Area;
  GMc1toSI.Edens  = GMc1toSI.Energy / GMc1toSI.Volume;
  GMc1toSI.Mdens  = GMc1toSI.Mass / GMc1toSI.Volume;
  GMc1toSI.Ndens  = 1. / GMc1toSI.Volume;
  GMc1toSI.Energy_MeV = GMc1toSI.Energy / (e * 1e6);
  GMc1toSI.Edens_MeV  = GMc1toSI.Edens  / (e * 1e6);
}


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

/* mass density in SI units */
double MassDensity_GMc1_to_SI(double rho)
{
  double L = Length_GMc1_to_SI(1.);
  double L3 = L*L*L;
  double Msun = Mass_GMc1_to_SI(1.);
  return rho * Msun/L3;
}

/* energy density in SI units */
double EnergyDensity_GMc1_to_SI(double rho)
{
  double c = natureconsts.c;
  return MassDensity_GMc1_to_SI(rho) * c*c;
}

/* pressure in SI units */
double Pressure_GMc1_to_SI(double P)
{
  double L = Length_GMc1_to_SI(1.);
  double a = Acceleration_GMc1_to_SI(1.);
  double Msun = Mass_GMc1_to_SI(1.);
  return P * Msun*a/(L*L);
}


/* test some entries of the SItoGMc1 and GMc1toSI structs */
void units_test_SItoGMc1_GMc1toSI(void)
{
  printf("1m in SI is %.19g\n", Length_GMc1_to_SI(SItoGMc1.m));
  printf("1s in SI is %.19g\n", Time_GMc1_to_SI(SItoGMc1.s));
  printf("1kg in SI is %.19g\n", Mass_GMc1_to_SI(SItoGMc1.kg));
  printf("1kg/m^3 in SI is %.19g\n", MassDensity_GMc1_to_SI(SItoGMc1.kg/SItoGMc1.m3));
  printf("1J/m^3 in SI is %.19g\n", EnergyDensity_GMc1_to_SI(SItoGMc1.J/SItoGMc1.m3));
  printf("1Pa in SI is %.19g\n", Pressure_GMc1_to_SI(SItoGMc1.Pa));

  printf("GMc1toSI.Length = %.19g\n", GMc1toSI.Length);
  printf("GMc1toSI.Area   = %.19g\n", GMc1toSI.Area);
  printf("GMc1toSI.Volume = %.19g\n", GMc1toSI.Volume);
  printf("GMc1toSI.Time   = %.19g\n", GMc1toSI.Time);
  printf("GMc1toSI.Mass   = %.19g\n", GMc1toSI.Mass);
  printf("GMc1toSI.Energy = %.19g\n", GMc1toSI.Energy);
  printf("GMc1toSI.Power  = %.19g\n", GMc1toSI.Power);
  printf("GMc1toSI.Force  = %.19g\n", GMc1toSI.Force);
  printf("GMc1toSI.Press  = %.19g\n", GMc1toSI.Press);
  printf("GMc1toSI.Edens  = %.19g\n", GMc1toSI.Edens);
  printf("GMc1toSI.Mdens  = %.19g\n", GMc1toSI.Mdens);
  printf("GMc1toSI.Ndens  = %.19g\n", GMc1toSI.Ndens);
  printf("GMc1toSI.Energy_MeV = %.19g\n", GMc1toSI.Energy_MeV);
  printf("GMc1toSI.Edens_MeV  = %.19g\n", GMc1toSI.Edens_MeV);

  printf("h in G=Msun=c=1 units is %.19g\n",
         natureconsts.h * SItoGMc1.J * SItoGMc1.s);
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
