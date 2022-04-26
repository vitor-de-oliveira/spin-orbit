#include "dynamical_system.h"

int field_two_body(double t, const double y[], double f[], 
				void *params)
{
	(void)t;

	double *par = (double *)params;

	double e 			= par[1];
	double m_primary 	= par[2]; 
	double m_secondary 	= par[3];
	double G 			= par[4];

	double total_mass = m_primary + m_secondary;

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
	(void)t;

	double *par = (double *)params;

	double e 			= par[1];
	double m_primary 	= par[2]; 
	double m_secondary 	= par[3];
	double G 			= par[4];

	double total_mass = m_primary + m_secondary;

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
	(void)t;

	double *par = (double *)params;

	double gamma 		= par[0];
	double e 			= par[1];
	double m_primary 	= par[2]; 
	double m_secondary 	= par[3];
	double G 			= par[4];

	double total_mass = m_primary + m_secondary;

    double r = sqrt((y[2] * y[2]) + (y[3] * y[3]));
    double f_e = atan2(y[3], y[2]);

	double aux = (-3.0/2.0) * gamma * G * m_primary;
	double r_cube = r * r * r;

	// y[0] = theta
	// y[1] = theta_dot
	// y[2] = x
	// y[3] = y
	// y[4] = x_dot
	// y[5] = y_dot

	f[0] = y[1];
	f[1] = aux * sin(2.0 * (y[0] - f_e)) / r_cube;
	f[2] = y[4];
	f[3] = y[5];
	f[4] = -1.0 * G * total_mass * y[2] / r_cube;
	f[5] = -1.0 * G * total_mass * y[3] / r_cube;

	return GSL_SUCCESS;
}

int jacobian_rigid(double t, const double y[], double *dfdy, 
					double dfdt[], void *params)
{
	(void)t;

	double *par = (double *)params;

	double gamma 		= par[0];
	double e 			= par[1];
	double m_primary 	= par[2]; 
	double m_secondary 	= par[3];
	double G 			= par[4];

	double total_mass = m_primary + m_secondary;

    double r = sqrt((y[2] * y[2]) + (y[3] * y[3]));
    double f_e = atan2(y[3], y[2]);

	double aux = (-3.0/2.0) * gamma * G * m_primary;
	double r_cube = r * r * r;

	double r_fifth = r * r * r * r * r;
	double alpha_1 = G * total_mass * 
			(2.0 * y[2] * y[2] - y[3] * y[3]) / r_fifth;
	double alpha_2 = 3.0 * G * total_mass * y[2] * y[3] 
					/ r_fifth;
	double alpha_3 = alpha_2;
	double alpha_4 = G * total_mass *  
			(2.0 * y[3] * y[3] - y[2] * y[2]) / r_fifth;

	double mixed_1 = (aux / r_fifth) * 
			(2.0*y[3]*cos(2.0*(y[0]-f_e)) - 
			 3.0*y[2]*sin(2.0*(y[0]-f_e)));

	double mixed_2 = (aux / r_fifth) * 
			(-3.0*y[3]*sin(2.0*(y[0]-f_e)) 
			 -2.0*y[2]*cos(2.0*(y[0]-f_e)));

	gsl_matrix_view dfdy_mat 
		= gsl_matrix_view_array(dfdy, 2, 2);
	gsl_matrix *mat = &dfdy_mat.matrix;

	gsl_matrix_set(mat, 0, 0, 0.0);
	gsl_matrix_set(mat, 0, 1, 1.0);
	gsl_matrix_set(mat, 1, 0, 2.0 * aux 
					*cos(2.0 * (y[0] - f_e)) / r_cube );
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

	gsl_matrix_set(mat, 0, 2, 0.0);
	gsl_matrix_set(mat, 0, 3, 0.0);
	gsl_matrix_set(mat, 0, 4, 0.0);
	gsl_matrix_set(mat, 0, 5, 0.0);
	gsl_matrix_set(mat, 1, 2, mixed_1);
	gsl_matrix_set(mat, 1, 3, mixed_2);
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

	double *par = (double *)params;

	double gamma 		= par[0];
	double e 			= par[1];
	double m_primary 	= par[2]; 
	double m_secondary 	= par[3];
	double G 			= par[4];

	double u = kepler_equation(e,t);
    double r = 1.0 - e * cos(u);
    double f_e = 2.0 * atan(sqrt((1.0 + e)/(1.0 - e)) 
			* tan(0.5 * u));

	double aux = (-3.0/2.0) * gamma * G * m_primary;
	double r_cube = r * r * r;
	
	f[0] = y[1];
	f[1] = aux * sin(2.0 * (y[0] - f_e)) / r_cube;

	return GSL_SUCCESS;
}

int jacobian_rigid_kepler(double t, const double y[], 
				double *dfdy, double dfdt[], void *params)
{
	double *par = (double *)params;

	double gamma 		= par[0];
	double e 			= par[1];
	double m_primary 	= par[2]; 
	double m_secondary 	= par[3];
	double G 			= par[4];

	double u = kepler_equation(e,t);
    double r = 1.0 - e * cos(u);
    double f_e = 2.0 * atan(sqrt((1.0 + e)/(1.0 - e)) 
			* tan(0.5 * u));

	double aux = (-3.0/2.0) * gamma * G * m_primary;
	double r_cube = r * r * r;

	gsl_matrix_view dfdy_mat 
		= gsl_matrix_view_array(dfdy, 2, 2);
	gsl_matrix *mat = &dfdy_mat.matrix;
	gsl_matrix_set(mat, 0, 0, 0.0);
	gsl_matrix_set(mat, 0, 1, 1.0);
	gsl_matrix_set(mat, 1, 0, 2.0 * aux 
					*cos(2.0 * (y[0] - f_e)) / r_cube );
	gsl_matrix_set(mat, 1, 1, 0.0);

	dfdt[0] = 0.0;
	dfdt[1] = 0.0;

	return GSL_SUCCESS;
}

double angular_momentum_two_body(double y[4])
{
	double h;

	h = y[0] * y[3] - y[1] * y[2];

	return h;
}

double vis_viva_two_body(double y[4], double T, double a)
{
	double n, mu, r, v, C;

	n = 2.0 * M_PI / T;

	mu = n * n * a * a * a;

	r = sqrt(y[0] * y[0] + y[1] * y[1]);

	v = sqrt(y[2] * y[2] + y[3] * y[3]);

	C = 0.5 * v * v - mu / r;

	return C;
}

int init_orbital(double orb[4], double e)
{
	double a = 1.0;
	
	double x = a * (1.0 - e * e) / (1.0 + e);
	double y = 0.0;

	double f_e = atan2(y, x);
	
	double x_dot = 0.0;
	double y_dot = (e + cos(f_e))/sqrt(1.0-e*e);
	// double y_dot = (e + 1.0 / 
	// 	sqrt(y[3]*y[3]/(y[2]*y[2])+1.0)) /
	// 	sqrt(1.0-e*e);

	orb[0] = x;
	orb[1] = y;
	orb[2] = x_dot;
	orb[3] = y_dot;	

	return 0;
}

dynsys init_rigid(void *params)
{
	dynsys rigid;
    rigid.name = "rigid";
    rigid.dim = 6;
	rigid.field = &field_rigid;
	rigid.jac = &jacobian_rigid;
	rigid.params = params;
    return rigid;
}

dynsys init_rigid_kepler(void *params)
{
	dynsys rigid_kepler;
    rigid_kepler.name = "rigid_kepler";
    rigid_kepler.dim = 2;
	rigid_kepler.field = &field_rigid_kepler;
	rigid_kepler.jac = &jacobian_rigid_kepler;
	rigid_kepler.params = params;
    return rigid_kepler;
}

dynsys init_two_body(void *params)
{
	dynsys two_body;
    two_body.name = "two_body";
    two_body.dim = 4;
	two_body.field = &field_two_body;
	two_body.jac = &jacobian_two_body;
	two_body.params = params;
    return two_body;
}

int evolve_cycle(double *y,	double cycle_period, 
				 double *t, dynsys system)
{
	// declare variables
	int status;
	double h, h_max, h_min;
	double error_abs, error_rel;
	char *method, *control;

	// initialize control variables
	h = 1e-3 * sign(cycle_period);
	error_abs = 1e-14;
	error_rel = 0.0;
	h_max = 1e-1;
	h_min = 1e-11;
	method = "rk8pd";
	control = "adaptive";

	// set system, integrator and driver
	gsl_odeiv2_system sys;
	set_system(&sys, system);

	const gsl_odeiv2_step_type *T;
	set_integrator(&T, method);

	gsl_odeiv2_driver *d;
	set_driver(&d, &sys, T, h, h_max, h_min,
			   error_abs, error_rel, control);

	// cycle evolution
	status = gsl_odeiv2_driver_apply (d, t, 
		*t + cycle_period, y);

	// check if integration was successfull
	if (status != GSL_SUCCESS)
	{
		printf("Warning: %s\n", gsl_strerror(status));
		exit(2);
	}

	// free memory
	gsl_odeiv2_driver_free(d);

	return 0;
}

int evolve_orbit(double *ic, double cycle_period, 
				int number_of_cycles, 
				double ***orbit, int *orbit_size,
				dynsys system)
{
	// declare variables
	double y[system.dim];
	double box = 1e5;
	double t = 0.0;

	// allocate memory and initializes exit data
	if (orbit != NULL)
	{
		// takes into consideration initial condition
		alloc_2d_double(orbit, number_of_cycles + 1, 
			system.dim);
		copy((*orbit)[0], ic, system.dim);
	}
	else
	{
		printf("Error: NULL orbit.");
		exit(2);
	}
	
	// takes into consideration initial condition
	int counter = 1;

	// orbit evolution
	copy(y, ic, system.dim);
	for (int i = 0; i < number_of_cycles; i++)
	{
		evolve_cycle(y, cycle_period, &t, system);
	
		// check if orbit diverges
		for (int j = 0; j < system.dim; j++)
		{
			if (fabs(y[j]) > box)
			{
				printf("Warning: box limit reached\n");
				goto out;
			}
		}

		counter++;

		// write orbit element
		copy((*orbit)[i + 1], y, system.dim);
	}

	out:;

	*orbit_size = counter;

	return 0;
}

int set_driver(gsl_odeiv2_driver **d, 
			   gsl_odeiv2_system *sys, 
			   const gsl_odeiv2_step_type *T, 
			   double h, double h_max, double h_min,
			   double error_abs, double error_rel, 
			   char *control)
{
	*d = gsl_odeiv2_driver_alloc_y_new(sys, 
		T, h, error_abs, error_rel);
		
	if (strcmp(control, "adaptive") == 0)
	{
		gsl_odeiv2_driver_set_hmax(*d, h_max);
		gsl_odeiv2_driver_set_hmin(*d, h_min);
	}
	else if (strcmp(control, "fixed") != 0)
	{
		printf("Warning: invalid GSL error control\n");
		exit(2);
	}
	return 0;
}

int set_integrator(const gsl_odeiv2_step_type **T, 
				   char *integrator)
{
    if (strcmp(integrator, "rk2") == 0)
	{
		*T = gsl_odeiv2_step_rk2;
	}
	else if (strcmp(integrator, "rk4") == 0)
	{
		*T = gsl_odeiv2_step_rk4;
	}
	else if (strcmp(integrator, "rk45") == 0)
	{
		*T = gsl_odeiv2_step_rkf45;
	}
	else if (strcmp(integrator, "rkck") == 0)
	{
		*T = gsl_odeiv2_step_rkck;
	}
	else if (strcmp(integrator, "rk8pd") == 0)
	{
		*T = gsl_odeiv2_step_rk8pd;
	}
	else if (strcmp(integrator, "rk4imp") == 0)
	{
		*T = gsl_odeiv2_step_rk4imp;
	}
	else if (strcmp(integrator, "bsimp") == 0)
	{
		*T = gsl_odeiv2_step_bsimp;
	}
	else if (strcmp(integrator, "msadams") == 0)
	{
		*T = gsl_odeiv2_step_msadams;
	}
	else if (strcmp(integrator, "msbdf") == 0)
	{
		*T = gsl_odeiv2_step_msbdf;	
	}
	else
	{
		printf("Warning: invalid GSL integrator\n");
		exit(2);
	}
	return 0;
}

int set_system(gsl_odeiv2_system *sys, 
			   dynsys system)
{
	sys->function = system.field;
	sys->jacobian = system.jac;
	sys->dimension = system.dim;
	sys->params = system.params;
	return 0;
}

double root_function_kepler(double u, void *params)
{
	struct root_params *p = (struct root_params *)params;
	double e = p->e;
	double t = p->t;
	return u - e * sin(u) - t;
}

double root_derivative_kepler(double u, void *params)
{
	struct root_params *p = (struct root_params *)params;
	double e = p->e;
	return 1.0 - e * cos(u);
}

void root_fdf_kepler(double u, void *params, double *y, 
                    double *dy)
{
	struct root_params *p = (struct root_params *)params;
	double e = p->e;
	double t = p->t;
	*y = u - e * sin(u) - t;
	*dy = 1.0 - e * cos(u);
}

double kepler_equation(double e, double t)
{
	int status, iter = 0, max_iter = 100;
	double u0, u;
	struct root_params params_root = {e, t};

    // initial guess (works up to e = 0.99)
    u = t + e * sin(t);

    // stonger initial guess
    // u = t + e * sin(t) + 0.5 * e * e * sin(2.0 * t)
        // + (e * e * e / 8.0) * (3.0 * sin (3.0 * t) 
        // - sin(t));

    const gsl_root_fdfsolver_type *T 
        = gsl_root_fdfsolver_steffenson;
    gsl_root_fdfsolver *s = gsl_root_fdfsolver_alloc(T);
    gsl_function_fdf FDF;
    FDF.f = &root_function_kepler;
    FDF.df = &root_derivative_kepler;
    FDF.fdf = &root_fdf_kepler;
    FDF.params = &params_root;
    gsl_root_fdfsolver_set(s, &FDF, u);

    do
    {
        iter++;
        gsl_root_fdfsolver_iterate(s);
        u0 = u;
        u = gsl_root_fdfsolver_root(s);
        status = gsl_root_test_delta(u, u0, 1e-15, 0);
        if (iter == max_iter)
        {
            printf("Warning: reached maximum iterate\n");
            return 1;
        }
    } while (status == GSL_CONTINUE);

    gsl_root_fdfsolver_free(s);

    return u;
}
