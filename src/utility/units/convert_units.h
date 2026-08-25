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
  const double m_e;    /* electron mass */
  const double m_u;    /* unified atomic mass unit, (unbound C atom mass)/12 */
  const double m_n;    /* neutron mass */
  const double m_p;    /* proton mass */
} tNatureconsts;


/***************************************************************************/
/* struct with conversion factors from SI to G=Msun=c=1 units */
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
  double MeV;  /* 1 MeV in nmesh units */
} tSItoGMc1;


/***************************************************************************/
/* struct with conversion factors from G=Msun=c=1 units to SI */
/***************************************************************************/

/* struct with values of SI units */
typedef struct
{
  double Length;  /* Length=1 in m */
  double Area;    /* Area=1 in m^2 */
  double Volume;  /* Volume=1 in m^3 */
  double Time;    /* Time=1 in s */
  double Mass ;   /* Mass=1 in kg */
  double Energy;  /* Energy=1 in J */
  double Power;   /* Power=1 in W=J/s */
  double Force;   /* Force=1 in N */
  double Press;   /* Press=1 in Pa=N/m^2=J/m^3 */
  double Edens;   /* Edens=1 in J/m^3 */
  double Mdens;   /* Mdens=1 in kg/m^3 */
  double Ndens;   /* Ndens=1 in 1/m^3 */
  double Energy_MeV;  /* Energy=1 in MeV */
  double Edens_MeV;   /* Edens=1 in MeV/m^3 */
} tGMc1toSI;


/***************************************************************************/
/* functions in convert_units.c */
/***************************************************************************/
void units_set_SItoGMc1(void);
void units_test_SItoGMc1_GMc1toSI(void);
void units_set_GMc1toSI(void);
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
