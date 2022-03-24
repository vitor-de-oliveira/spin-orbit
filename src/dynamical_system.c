#include "dynamical_system.h"

int field_two_body(double t, const double y[], double f[], 
				void *params)
{
	double G = 1.0;
	double m_secondary = 1.215e-2; //Moon
	double m_primary = 1.0 - m_secondary;
	double total_mass = m_primary + m_secondary;

	(void)t;
	double e = *(double *)params;
    double r = sqrt((y[0] * y[0]) + (y[1] * y[1]));
	double r_cube = r * r * r;

	f[0] = y[2];
	f[1] = y[3];
	f[2] = -1.0 * G * total_mass * y[0] / r_cube;
	f[3] = -1.0 * G * total_mass * y[1] / r_cube;

	return GSL_SUCCESS;
}

int jacobian_two_body(double t, const double y[], double *dfdy, 
					double dfdt[], void *params)
{
	double G = 1.0;
	double m_secondary = 1.215e-2; //Moon
	double m_primary = 1.0 - m_secondary;
	double total_mass = m_primary + m_secondary;

	(void)t;
	double e = *(double *)params;
    double r = sqrt((y[0] * y[0]) + (y[1] * y[1]));
	double r_fifth = r * r * r * r * r;

	double alpha_1 = G * total_mass * 
			(2.0 * y[0] * y[0] - y[1] * y[1]) / r_fifth;

	double alpha_2 = 3.0 * G * total_mass * y[0] * y[1] 
					/ r_fifth;

	double alpha_3 = alpha_2;

	double alpha_4 = G * total_mass *  
			(2.0 * y[1] * y[1] - y[0] * y[0]) / r_fifth;

	gsl_matrix_view dfdy_mat 
		= gsl_matrix_view_array(dfdy, 4, 4);
	gsl_matrix *mat = &dfdy_mat.matrix;
	gsl_matrix_set(mat, 0, 0, 0.0);
	gsl_matrix_set(mat, 0, 1, 0.0);
	gsl_matrix_set(mat, 0, 2, 1.0);
	gsl_matrix_set(mat, 0, 3, 0.0);
	gsl_matrix_set(mat, 1, 0, 0.0);
	gsl_matrix_set(mat, 1, 1, 0.0);
	gsl_matrix_set(mat, 1, 2, 0.0);
	gsl_matrix_set(mat, 1, 3, 1.0);
	gsl_matrix_set(mat, 2, 0, alpha_1);
	gsl_matrix_set(mat, 2, 1, alpha_2);
	gsl_matrix_set(mat, 2, 2, 0.0);
	gsl_matrix_set(mat, 2, 3, 0.0);
	gsl_matrix_set(mat, 3, 0, alpha_3);
	gsl_matrix_set(mat, 3, 1, alpha_4);
	gsl_matrix_set(mat, 3, 2, 1.0);
	gsl_matrix_set(mat, 3, 3, 0.0);

	dfdt[0] = 0.0;
	dfdt[1] = 0.0;
	dfdt[2] = 0.0;
	dfdt[3] = 0.0;

	return GSL_SUCCESS;
}

int field_rigid(double t, const double y[], double f[], 
				void *params)
{
	double G = 1.0;
	double m = 1.0;
	double gamma = 0.01;

	double m_secondary = 1.215e-2; //Moon
	double m_primary = 1.0 - m_secondary;
	double total_mass = m_primary + m_secondary;

	(void)t;
	double e = *(double *)params;
    double r = sqrt((y[2] * y[2]) + (y[3] * y[3]));
	double r_cube = r * r * r;
    double f_e = atan2(y[3], y[2]) + M_PI;
	double aux = (-3.0/2.0) * gamma * G * m;

	// y[0] = theta
	// y[1] = theta_dot
	// y[2] = x
	// y[3] = y
	// y[4] = x_dot
	// y[5] = y_dot

	f[0] = y[1];
	f[1] = aux * sin(2.0 * (f_e - y[0])) / r_cube;

	f[2] = y[4];
	f[3] = y[5];
	f[4] = -1.0 * G * total_mass * y[2] / r_cube;
	f[5] = -1.0 * G * total_mass * y[3] / r_cube;

	return GSL_SUCCESS;
}

int jacobian_rigid(double t, const double y[], double *dfdy, 
					double dfdt[], void *params)
{
	double G = 1.0;
	double m = 1.0;
	double gamma = 0.01;

	double m_secondary = 1.215e-2; //Moon
	double m_primary = 1.0 - m_secondary;
	double total_mass = m_primary + m_secondary;

	double e = *(double *)params;
    double r = sqrt((y[2] * y[2]) + (y[3] * y[3]));
	double r_cube = r * r * r;
    double f_e = atan2(y[3], y[2]) + M_PI;
	double aux = (-3.0/2.0) * gamma * G * m;

	double r_fifth = r * r * r * r * r;
	double alpha_1 = G * total_mass * 
			(2.0 * y[2] * y[2] - y[3] * y[3]) / r_fifth;
	double alpha_2 = 3.0 * G * total_mass * y[2] * y[3] 
					/ r_fifth;
	double alpha_3 = alpha_2;
	double alpha_4 = G * total_mass *  
			(2.0 * y[3] * y[3] - y[2] * y[2]) / r_fifth;

	gsl_matrix_view dfdy_mat 
		= gsl_matrix_view_array(dfdy, 2, 2);
	gsl_matrix *mat = &dfdy_mat.matrix;
	gsl_matrix_set(mat, 0, 0, 0.0);
	gsl_matrix_set(mat, 0, 1, 1.0);
	gsl_matrix_set(mat, 1, 0, -2.0 * aux 
					*cos(2.0 * (f_e - y[0])) );
	gsl_matrix_set(mat, 1, 1, 0.0);

	gsl_matrix_set(mat, 2, 2, 0.0);
	gsl_matrix_set(mat, 2, 3, 0.0);
	gsl_matrix_set(mat, 2, 4, 1.0);
	gsl_matrix_set(mat, 2, 5, 0.0);
	gsl_matrix_set(mat, 3, 2, 0.0);
	gsl_matrix_set(mat, 3, 3, 0.0);
	gsl_matrix_set(mat, 3, 4, 0.0);
	gsl_matrix_set(mat, 3, 5, 1.0);
	gsl_matrix_set(mat, 4, 2, alpha_1);
	gsl_matrix_set(mat, 4, 3, alpha_2);
	gsl_matrix_set(mat, 4, 4, 0.0);
	gsl_matrix_set(mat, 4, 5, 0.0);
	gsl_matrix_set(mat, 5, 2, alpha_3);
	gsl_matrix_set(mat, 5, 3, alpha_4);
	gsl_matrix_set(mat, 5, 4, 1.0);
	gsl_matrix_set(mat, 5, 5, 0.0);

	// THINK ABOUT THIS! I AM NOT SURE
	gsl_matrix_set(mat, 0, 2, 0.0);
	gsl_matrix_set(mat, 0, 3, 0.0);
	gsl_matrix_set(mat, 0, 4, 0.0);
	gsl_matrix_set(mat, 0, 5, 0.0);
	gsl_matrix_set(mat, 1, 2, 0.0);
	gsl_matrix_set(mat, 1, 3, 0.0);
	gsl_matrix_set(mat, 1, 4, 0.0);
	gsl_matrix_set(mat, 1, 5, 0.0);
	gsl_matrix_set(mat, 2, 0, 0.0);
	gsl_matrix_set(mat, 2, 1, 0.0);
	gsl_matrix_set(mat, 3, 0, 0.0);
	gsl_matrix_set(mat, 3, 1, 0.0);
	gsl_matrix_set(mat, 4, 0, 0.0);
	gsl_matrix_set(mat, 4, 1, 0.0);
	gsl_matrix_set(mat, 5, 0, 0.0);
	gsl_matrix_set(mat, 5, 1, 0.0);

	dfdt[0] = 0.0;
	dfdt[1] = 0.0;
	dfdt[2] = 0.0;
	dfdt[3] = 0.0;
	dfdt[4] = 0.0;
	dfdt[5] = 0.0;

	return GSL_SUCCESS;
}

int field_rigid_kepler(double t, const double y[], 
							double f[], void *params)
{
	double G = 1.0;
	double m = 1.0;
	double gamma = 0.01;

	double e = *(double *)params;
	double u = kepler_equation(e,t);
    double r = 1.0 - e * cos(u);
    double f_e = 2.0 * atan(sqrt((1.0 + e)/(1.0 - e)) 
			* tan(0.5 * u));
	double aux = (-3.0/2.0) * gamma * G * m / (r * r * r);

	f[0] = y[1];
	f[1] = aux * sin(2.0 * (f_e - y[0]));

	return GSL_SUCCESS;
}

int jacobian_rigid_kepler(double t, const double y[], 
				double *dfdy, double dfdt[], void *params)
{
	double G = 1.0;
	double m = 1.0;
	double gamma = 0.01;

	double e = *(double *)params;
	double u = kepler_equation(e,t);
    double r = 1.0 - e * cos(u);
    double f_e = 2.0 * atan(sqrt((1.0 + e)/(1.0 - e)) 
			* tan(0.5 * u));
	double aux = (-3.0/2.0) * gamma * G * m / (r * r * r);

	gsl_matrix_view dfdy_mat 
		= gsl_matrix_view_array(dfdy, 2, 2);
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