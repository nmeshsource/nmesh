/* convert_units.c */
/* Wolfgang Tichy 6/2022 */

#include <stdio.h>

/***************************************************************************/
/* struct with constants of nature */
/***************************************************************************/

/* constants of nature */
struct tNATURECONSTS
{
  const double c;      /* speed of light */
  const double G;      /* gravitational constant */
  const double GMsun;  /* G * solar Mass */
  const double h;      /* Planck's constant */
  const double k_B;    /* Boltzmann constant */
  const double e;      /* electron charge */
  const double m_u;    /* unified atomic mass unit, (unbound C atom mass)/12 */
};

/* SI values of consts of nature we use */
struct tNATURECONSTS natureconsts =
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
  /* From: CODATA values of the fundamental physical constants 2022: */
  .h     = 6.62607015e-34, // [J s],         h is defined to have this value
  .k_B   = 1.380649e-23,   // [J/K],       k_B is defined to have this value
  .e     = 1.602176634e-19,// [C],           e is defined to have this value
  .m_u   = 1.6605390689252e-27, //[kg],  measured, last 2 digits uncertain
};


/***************************************************************************/
/* struct with conversion factors from SI to nmesh units */
/***************************************************************************/

/* struct with values of SI units */
struct tUNIT
{
  double m;    /* 1 m in nmesh units */
  double m2;   /* 1 m^2 in nmesh units */
  double m3;   /* 1 m^3 in nmesh units */
  double s;    /* 1 s in nmesh units */
  double s2;   /* 1 s^2 in nmesh units */
  double kg;   /* 1 kg in nmesh units */
  double J;    /* 1 J in nmesh units */
  double W;    /* 1 W in nmesh units */
  double N;    /* 1 N in nmesh units */
  double Pa;   /* 1 Pa = 1 N/m^2 = 1 J/m^3 in nmesh units */
  double eV;   /* 1 eV in nmesh units */
};


/***************************************************************************/
/* convert from SI to nmesh units */
/***************************************************************************/

/* global struct that contains the values of various SI units */
struct tUNIT unit;

/* set the values of SI units for G=Msun=c=1 */
void unit_set_for_GMc1(void)
{
  double c     = natureconsts.c;
  double G     = natureconsts.G;
  double GMsun = natureconsts.GMsun;
  double Msun  = GMsun/G;
  double e     = natureconsts.e;

  unit.m  = (c*c)/GMsun;
  unit.m2 = unit.m * unit.m;
  unit.m3 = unit.m * unit.m2;
  unit.s  = (c*c*c)/GMsun;
  unit.s2 = unit.s * unit.s;
  unit.kg = 1./Msun;
  unit.J  = 1./(Msun*c*c);
  unit.W  = unit.J / unit.s;
  unit.N  = unit.J / unit.m;
  unit.Pa = unit.N / unit.m2; // = unit.J / unit.m3;
  unit.eV = unit.J * e;
}


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
