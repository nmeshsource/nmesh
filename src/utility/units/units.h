/* units.h */
/* Wolfgang Tichy, June 2022 */

/* convert_units.c */
void unit_set_for_GMc1(void);
void unit_test_values(void);
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

/* units.c */
int unit_set_for_GMc1_mesh(tMesh *mesh);
int print_unit_conversion_factors(tMesh *mesh);
