#include "dynamical_system.h"

// only for the time being
#define G 1.0
#define m 1.0
#define gamma 0.01

int field_rigid(double t, const double y[], double f[], 
				void *params)
{
	double e = *(double *)params;
	double u = kepler_equation(e,t);
    double r = 1 - e * cos(u);
    double f_e = 2.0 * atan(sqrt((1.0 + e)/(1.0 - e)) 
			* tan(0.5 * u));
	double aux = (-3.0/2.0) * gamma * G * m / (r * r * r);

	f[0] = y[1];
	f[1] = aux * sin(2.0 * (f_e - y[0]));

	return GSL_SUCCESS;
}

int jacobian_rigid(double t, const double y[], double *dfdy, 
					double dfdt[], void *params)
{
	double e = *(double *)params;
	double u = kepler_equation(e,t);
    double r = 1 - e * cos(u);
    double f_e = 2.0 * atan(sqrt((1.0 + e)/(1.0 - e)) 
			* tan(0.5 * u));
	double aux = (-3.0/2.0) * gamma * G * m / (r * r * r);

	gsl_matrix_view dfdy_mat 
		= gsl_matrix_view_array(dfdy, 4, 4);
	gsl_matrix *mat = &dfdy_mat.matrix;
	gsl_matrix_set(mat, 0, 0, 0.0);
	gsl_matrix_set(mat, 0, 1, 1.0);
	gsl_matrix_set(mat, 1, 0, -2.0 * aux 
					*cos(2.0 * (f_e - y[0])) );
	gsl_matrix_set(mat, 1, 1, 0.0);

	dfdt[0] = 0.0;
	dfdt[1] = 0.0;

	return GSL_SUCCESS;
}