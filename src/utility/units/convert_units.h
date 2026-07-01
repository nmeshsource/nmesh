/* convert_units.h */
/* Wolfgang Tichy, 7/2026 */

/***************************************************************************/
/* struct with constants of nature */
/***************************************************************************/

/* constants of nature in SI units */
typedef struct
{
  const double c;      /* speed of light */
  const double G;      /* gravitational constant */
  const double GMsun;  /* G * solar Mass */
  const double h;      /* Planck's constant */
  const double k_B;    /* Boltzmann constant */
  const double e;      /* electron charge */
  const double m_u;    /* unified atomic mass unit, (unbound C atom mass)/12 */
} tNatureconsts;


/***************************************************************************/
/* struct with conversion factors from SI to nmesh units */
/***************************************************************************/

/* struct with values of SI units */
typedef struct
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
} tUnits;


/***************************************************************************/
/* functions in convert_units.c */
/***************************************************************************/
void units_set_for_GMc1(void);
void units_test_values(void);
double Length_GMc1_to_SI(double x);
double Time_GMc1_to_SI(double t);
double Frequency_GMc1_to_SI(double f);
double Speed_GMc1_to_SI(double v);
double Acceleration_GMc1_to_SI(double a);
double Mass_GMc1_to_SI(double m);
double MassDensity_GMc1_to_SI(double rho);
double EnergyDensity_GMc1_to_SI(double rho);
double Pressure_GMc1_to_SI(double P);
void print_natureconsts(void);
void print_conversion_factors(void);
