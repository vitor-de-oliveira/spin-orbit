#include "dynamical_system.h"

int field_rigid(double t, const double y[], double f[], void *params)
{
	(void)(t);
	double mu = *(double *)params;
	double mu_1 = 1. - mu;
	double mu_2 = mu;
	double r_1 = sqrt((y[0] + mu_2) * (y[0] + mu_2) + y[1] * y[1]);
	double r_2 = sqrt((y[0] - mu_1) * (y[0] - mu_1) + y[1] * y[1]);
	double omega_x = y[0] - mu_1 * (y[0] + mu_2) / (r_1 * r_1 * r_1) - mu_2 * (y[0] - mu_1) / (r_2 * r_2 * r_2);
	double omega_y = y[1] - mu_1 * y[1] / (r_1 * r_1 * r_1) - mu_2 * y[1] / (r_2 * r_2 * r_2);
	f[0] = y[2];
	f[1] = y[3];
	f[2] = omega_x + 2.0 * y[3];
	f[3] = omega_y - 2.0 * y[2];
	return GSL_SUCCESS;
}

int jacobian_rigid(double t, const double y[], double *dfdy, double dfdt[], void *params)
{
	(void)(t);
	double mu = *(double *)params;
	double mu_1 = 1. - mu;
	double mu_2 = mu;
	double r_1 = sqrt((y[0] + mu_2) * (y[0] + mu_2) + y[1] * y[1]);
	double r_2 = sqrt((y[0] - mu_1) * (y[0] - mu_1) + y[1] * y[1]);
	double K1 = (1. - mu) / (r_1 * r_1 * r_1);
	double K2 = mu / (r_2 * r_2 * r_2);
	double A = (y[0] + mu) * (y[0] + mu) / (r_1 * r_1);
	double B = (y[0] + mu - 1.) * (y[0] + mu - 1.) / (r_2 * r_2);
	double C = y[1] * y[1] / (r_1 * r_1);
	double D = y[1] * y[1] / (r_2 * r_2);
	double E = (1. - mu) * (y[0] + mu) / (r_1 * r_1 * r_1 * r_1 * r_1);
	double F = mu * (y[0] + mu - 1.) / (r_2 * r_2 * r_2 * r_2 * r_2);
	double omega_xx = 1. - K1 * (1. - 3. * A) - K2 * (1. - 3. * B);
	double omega_yy = 1. - K1 * (1. - 3. * C) - K2 * (1. - 3. * D);
	double omega_xy = 3. * y[1] * (E + F);
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 4, 4);
	gsl_matrix *m = &dfdy_mat.matrix;
	gsl_matrix_set(m, 0, 0, 0.0);
	gsl_matrix_set(m, 0, 1, 0.0);
	gsl_matrix_set(m, 0, 2, 1.0);
	gsl_matrix_set(m, 0, 3, 0.0);
	gsl_matrix_set(m, 1, 0, 0.0);
	gsl_matrix_set(m, 1, 1, 0.0);
	gsl_matrix_set(m, 1, 2, 0.0);
	gsl_matrix_set(m, 1, 3, 1.0);
	gsl_matrix_set(m, 2, 0, omega_xx);
	gsl_matrix_set(m, 2, 1, omega_xy);
	gsl_matrix_set(m, 2, 2, 0.0);
	gsl_matrix_set(m, 2, 3, 2.0);
	gsl_matrix_set(m, 3, 0, omega_xy);
	gsl_matrix_set(m, 3, 1, omega_yy);
	gsl_matrix_set(m, 3, 2, -2.0);
	gsl_matrix_set(m, 3, 3, 0.0);
	dfdt[0] = 0.0;
	dfdt[1] = 0.0;
	dfdt[2] = 0.0;
	dfdt[3] = 0.0;
	return GSL_SUCCESS;
}