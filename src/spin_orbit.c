#include "spin_orbit.h"

int field_two_body(double t, const double y[], double f[], 
				void *params)
{
	(void)t;

	double *par = (double *)params;

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
	double a 			= par[5];

	double u = kepler_equation(e,t);
    double r = a * (1.0 - e * cos(u));
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
	double a 			= par[5];

	double u = kepler_equation(e,t);
    double r = a * (1.0 - e * cos(u));
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

int field_linear(double t, const double y[], 
				double f[], void *params)
{
	(void)t;

	double *par = (double *)params;

	double gamma 		= par[0];
	double e 			= par[1];
	double m_primary 	= par[2]; 
	double m_secondary 	= par[3];
	double G 			= par[4];
	double a			= par[5];
	double K			= par[6];
	double T			= par[7];

	double total_mass = m_primary + m_secondary;

    double r = sqrt((y[2] * y[2]) + (y[3] * y[3]));
    double f_e = atan2(y[3], y[2]);

	double aux = (-3.0/2.0) * gamma * G * m_primary;
	double r_cube = r * r * r;

	double a_over_r = a / r;
	double aux_2 = pow(a_over_r, 6.0);
	
	double n = 2.0 * M_PI / T;
	double f_dot = a_over_r * a_over_r * n * sqrt(1.0 - (e * e));

	double L = aux_2;
	double N = aux_2 * f_dot;

	// y[0] = theta
	// y[1] = theta_dot
	// y[2] = x
	// y[3] = y
	// y[4] = x_dot
	// y[5] = y_dot

	f[0] = y[1];
	f[1] = aux * sin(2.0 * (y[0] - f_e)) / r_cube 
			- K * (L * y[1] - N);
	f[2] = y[4];
	f[3] = y[5];
	f[4] = -1.0 * G * total_mass * y[2] / r_cube;
	f[5] = -1.0 * G * total_mass * y[3] / r_cube;

	return GSL_SUCCESS;
}

int field_linear_average(double t, const double y[], 
						double f[], void *params)
{
	(void)t;

	double *par = (double *)params;

	double gamma 		= par[0];
	double e 			= par[1];
	double m_primary 	= par[2]; 
	double m_secondary 	= par[3];
	double G 			= par[4];

	double K			= par[6];

	double total_mass = m_primary + m_secondary;

    double r = sqrt((y[2] * y[2]) + (y[3] * y[3]));
    double f_e = atan2(y[3], y[2]);

	double aux = (-3.0/2.0) * gamma * G * m_primary;
	double r_cube = r * r * r;

	double e_2 = e * e;
	double e_4 = e * e * e * e;
	double e_6 = e * e * e * e * e * e;

	double L_avg = (1.+3.*e_2+(3./8.)*e_4) / 
					pow(1.-e_2,(9./2.));

	double N_avg = (1.+(15./2.)*e_2+(45./8.)*e_4+
					(5./16.)*e_6) / pow(1.-e_2,6.);

	// y[0] = theta
	// y[1] = theta_dot
	// y[2] = x
	// y[3] = y
	// y[4] = x_dot
	// y[5] = y_dot

	f[0] = y[1];
	f[1] = aux * sin(2.0 * (y[0] - f_e)) / r_cube 
			- K * (L_avg * y[1] - N_avg);
	f[2] = y[4];
	f[3] = y[5];
	f[4] = -1.0 * G * total_mass * y[2] / r_cube;
	f[5] = -1.0 * G * total_mass * y[3] / r_cube;

	return GSL_SUCCESS;
}

double root_function_kepler(double u, void *params)
{
	struct root_params_kepler *p = 
            (struct root_params_kepler *)params;
	double e = p->e;
	double t = p->t;
	return u - e * sin(u) - t;
}

double root_derivative_kepler(double u, void *params)
{
	struct root_params_kepler *p = 
            (struct root_params_kepler *)params;
	double e = p->e;
	return 1.0 - e * cos(u);
}

void root_fdf_kepler(double u, void *params, double *y, 
                    double *dy)
{
	struct root_params_kepler *p = 
            (struct root_params_kepler *)params;
	double e = p->e;
	double t = p->t;
	*y = u - e * sin(u) - t;
	*dy = 1.0 - e * cos(u);
}

double kepler_equation(double e, double t)
{
	int status, iter = 0, max_iter = 100;
	double u0, u;
	struct root_params_kepler params_root = {e, t};

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

double angular_momentum_two_body(double y[4])
{
	double h;

	h = y[0] * y[3] - y[1] * y[2];

	return h;
}

double vis_viva_two_body(double y[4],
                         dynsys system)
{
	double *par = (double *)system.params;
	double a = par[5];
	double T = par[7];
	double n, mu, r, v, C;

	n = 2.0 * M_PI / T;

	mu = n * n * a * a * a;

	r = sqrt(y[0] * y[0] + y[1] * y[1]);

	v = sqrt(y[2] * y[2] + y[3] * y[3]);

	C = 0.5 * v * v - mu / r;

	return C;
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

dynsys init_linear(void *params)
{
	dynsys linear;
    linear.name = "linear";
    linear.dim = 6;
	linear.field = &field_linear;
	linear.jac = NULL;
	linear.params = params;
    return linear;
}

dynsys init_linear_average(void *params)
{
	dynsys linear_average;
    linear_average.name = "linear_average";
    linear_average.dim = 6;
	linear_average.field = &field_linear_average;
	linear_average.jac = NULL;
	linear_average.params = params;
    return linear_average;
}

int init_orbital(double orb[4],
                 dynsys system)
{
	double	*par = (double *)system.params;
	double	gamma = par[0];
	double	e = par[1];
	double	a = par[5];
	double	T = par[7];

	double	n = 2.0 * M_PI / T;
	
	// position and velocity at periapsis from Murray
	double	x = a * (1.0 - e);
	double	y = 0.0;	
	double	x_dot = 0.0;
	double	y_dot = n * a * sqrt((1.0 + e)/(1.0 - e));

	orb[0] = x;
	orb[1] = y;
	orb[2] = x_dot;
	orb[3] = y_dot;

	return 0;
}

int complete_orbital_part   (double y[],
                             dynsys system)
{
	double orb_ini[4];
	init_orbital(orb_ini, system);
	if (system.dim == 6)
	{
		for (int k = 0; k < 4; k++)
		{
			y[k+2] = orb_ini[k];
		}
	}
	return 0;
}

double kepler_period(double m1,
					 double m2,
					 double G,
					 double a)
{
	return sqrt(4.0 * M_PI * M_PI * a * a * a / (G * (m1 + m2)));
}

int orbit_two_body	(double *ic,
					 dynsys system,
					 anlsis analysis,
					 rngkta rk)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/orbit", &st) == -1) {
		mkdir("output/orbit", 0700);
	}

	if (system.name != "two_body")
	{
		printf("Warning: cannot use this function.\n");
		exit(2);
	}

	printf("Calculating orbit\n");

	double *par = (double *)system.params;
	double e = par[1];

	// prepare and open exit files 
	FILE	*out_orb, *out_orb_ic, 
			*out_orb_err;
	char	filename[150];

	sprintf(filename, "output/orbit/orbit_two_body_e_%1.3f.dat", e);
	out_orb = fopen(filename, "w");

	sprintf(filename, "output/orbit/orbit_ic_two_body_e_%1.3f.dat", e);
	out_orb_ic = fopen(filename, "w");

	sprintf(filename, "output/orbit/orbit_orbital_error_angular_momentum_two_body_e_%1.3f.dat", e);
	out_orb_err = fopen(filename, "w");
	
	// declare variables
	int orbit_size;
	double **orbit;
	double orb[4], orb_ini[4];
	double last_angle_dif;

	for (int i = 0; i < 4; i++)
	{
		orb_ini[i] = ic[i];
	}

	// evolve system
	evolve_orbit(ic, &orbit, &orbit_size, system, analysis, rk);

	// write orbit and constant error to file
	fprintf(out_orb_ic, "%1.15e %1.15e\n", 
			orbit[0][0], orbit[0][1]);

	for (int i = 0; i < orbit_size; i++)
	{
		fprintf(out_orb, "%1.15e %1.15e\n", 
				orbit[i][0], orbit[i][1]);
		
		for (int j = 0; j < 4; j++)
		{
			orb[j] = orbit[i][j];
		}
		fprintf(out_orb_err, "%d %1.15e\n", 
				i, fabs(angular_momentum_two_body(orb)-
				angular_momentum_two_body(orb_ini)));
	}

	// free memory
	dealloc_2d_double(&orbit, analysis.number_of_cycles);

	// close files
	fclose(out_orb);
	fclose(out_orb_ic);
	fclose(out_orb_err);

	printf("Data written in output/orbit/\n");

	return 0;
}

int orbit_map(double *ic, dynsys system,
					anlsis analysis, rngkta rk)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/orbit", &st) == -1) {
		mkdir("output/orbit", 0700);
	}

	printf("Calculating orbit\n");

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];
	double T = par[7];

	// prepare and open exit files 
	FILE	*out_orb, *out_orb_ic, 
			*out_orb_err;
	char	filename[150];

	sprintf(filename, "output/orbit/orbit_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f.dat", 
		gamma, e, system.name, K);
	out_orb = fopen(filename, "w");

	sprintf(filename, "output/orbit/orbit_ic_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f.dat", 
		gamma, e, system.name, K);
	out_orb_ic = fopen(filename, "w");

	sprintf(filename, "output/orbit/orbit_orbital_error_angular_momentum_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f.dat", 
		gamma, e, system.name, K);
	out_orb_err = fopen(filename, "w");
	
	// declare variables
	int orbit_size;
	int gcd, numerator, denominator;
	double **orbit;
	double orb[4], orb_ini[4];
	double last_angle_dif;

	for (int i = 0; i < 4; i++)
	{
		orb_ini[i] = ic[i+2];
	}

	// evolve system
	evolve_orbit(ic, &orbit, &orbit_size, system, analysis, rk);

	// write orbit and constant error to file
	fprintf(out_orb_ic, "%1.15e %1.15e\n", 
			angle_mod(orbit[0][0]), orbit[0][1]);

	for (int i = 0; i < orbit_size; i++)
	{
		// fprintf(out_orb, "%1.15e %1.15e\n", 
		// 		angle_mod_pos(orbit[i][0]), orbit[i][1]);
		fprintf(out_orb, "%1.15e %1.15e\n", 
				angle_mod(orbit[i][0]), orbit[i][1]);
		
		// if (strcmp(system.name, "rigid") == 0)
		if (system.dim == 6)
		{
			for (int j = 0; j < 4; j++)
			{
				orb[j] = orbit[i][j+2];
			}
			fprintf(out_orb_err, "%d %1.15e\n", 
					i, fabs(angular_momentum_two_body(orb)-
					angular_momentum_two_body(orb_ini)));
		}
	}

	// if (analysis.convergence_transient_wn + analysis.convergence_window_wn < orbit_size)
	// {
	// 	FILE 	*out_wn = fopen("output/tests/orbit_progress_wn.dat", "w");
	// 	double 	winding_number_window[analysis.convergence_window_wn];
	// 	for (int i = analysis.convergence_transient_wn + 1; i < orbit_size; i++)
	// 	{
	// 		int 	i_local = i - (analysis.convergence_transient_wn + 1);
	// 		double 	delta_theta = orbit[i][0] - orbit[analysis.convergence_transient_wn][0];
	// 		double 	wn = (delta_theta / (double) (i_local + 1));
	// 		fprintf(out_wn, "%d %f\n", i, wn);

	// 		if(i_local < analysis.convergence_window_wn)
	// 		{
	// 			winding_number_window[i_local] = wn;
	// 		}
	// 		else
	// 		{
	// 			double max_winding_number = wn;
	// 			double min_winding_number = wn;
	// 			for (int j = 0; j < analysis.convergence_window_wn - 1; j++)
	// 			{
	// 				winding_number_window[j] = winding_number_window[j+1];
	// 				if(winding_number_window[j] > max_winding_number)
	// 				{
	// 					max_winding_number = winding_number_window[j];
	// 				}
	// 				if(winding_number_window[j] < min_winding_number)
	// 				{
	// 					min_winding_number = winding_number_window[j];
	// 				}
	// 			}
	// 			winding_number_window[analysis.convergence_window_wn - 1] = wn;
	// 			if (fabs(max_winding_number-min_winding_number) < analysis.convergence_precision_wn)
	// 			{
	// 				printf("Winding number converged at iterate number %d\n", i);
	// 				printf("wn = %1.10e\n", wn);
	// 				double precise_wn, pwn_numerator = 0.0, pwn_denominator = 0.0;
	// 				for (int l = analysis.convergence_transient_wn + 2; l < i; l++)
	// 				{
	// 					double twp = ((double) (l - analysis.convergence_transient_wn - 1)) / ((double) (i - analysis.convergence_transient_wn - 1));
	// 					double factor = 1.0/exp(1.0/(twp*(1.0-twp)));
	// 					pwn_numerator += (orbit[l][0] - orbit[l-1][0]) * factor;
	// 					pwn_denominator += factor;
	// 				}
	// 				precise_wn = pwn_numerator / pwn_denominator;
	// 				printf("precise_wn = %1.10e\n", precise_wn);
	// 				break;
	// 			}
	// 		}
	// 	}
	// 	fclose(out_wn);
	// }

	// free memory
	dealloc_2d_double(&orbit, analysis.number_of_cycles);

	// close files
	fclose(out_orb);
	fclose(out_orb_ic);
	fclose(out_orb_err);

	printf("Data written in output/orbit/\n");

	return 0;
}

int phase_space(dynsys system, anlsis analysis, rngkta rk)
{
	// little warning
	if (strcmp(system.name, "two_body") == 0)
	{
		printf("Warning: cant draw phase space\n");
		printf("for two-body system\n");
		exit(2);
	}

	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/phase_space", &st) == -1) {
		mkdir("output/phase_space", 0700);
	}

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];

	// prepare and open exit files 
	FILE 	*out_psp, *out_ic,
			*out_orb_ang_mom_err, *out_vis_viva_err;
	char	filename[100];

	sprintf(filename, 
		"output/phase_space/phase_space_gamma_%1.6f_e_%1.3f.dat", gamma, e);
	out_psp = fopen(filename, "w");

	sprintf(filename, 
		"output/phase_space/phase_space_initial_conditions_gamma_%1.6f_e_%1.3f.dat", 
		gamma, e);
	out_ic = fopen(filename, "w");

	sprintf(filename, 
		"output/phase_space/phase_space_orbital_angular_momentum_error_gamma_%1.6f_e_%1.3f.dat",
		gamma, e);
	out_orb_ang_mom_err = fopen(filename, "w");

	sprintf(filename, 
		"output/phase_space/phase_space_vis_viva_error_gamma_%1.6f_e_%1.3f.dat",
		gamma, e);
	out_vis_viva_err = fopen(filename, "w");

	// declare variables
	double y[system.dim];
	double coordinate, velocity;
	int orbit_fw_size, orbit_bw_size;
	double **orbit_fw, **orbit_bw;
	double orb[4], orb_ini[4];
	init_orbital(orb_ini, system);
	anlsis analysis_fw, analysis_bw;

	analysis_fw = analysis;
	analysis_bw = analysis;
	analysis_bw.cycle_period *= -1.0;

	printf("Calculating phase space for e = %1.3f and gamma = %1.6f\n", e, gamma);

	// loop over coordinate values
	for (int i = 0; i < analysis.nc; i++)
	{
		// print progress on coordinate
		printf("Calculating set %d of %d\n", i + 1, analysis.nc);

		#pragma omp parallel private(y, coordinate, velocity, \
				orbit_fw_size, orbit_bw_size, orbit_fw, orbit_bw, orb) //num_threads(THREADS_NUM)
		{

		if (analysis.nc == 1)
		{
			coordinate = analysis.coordinate_min;
		}
		else
		{
			coordinate = analysis.coordinate_min + 
					(double)(i) * (analysis.coordinate_max - 
					analysis.coordinate_min) / 
					(double)(analysis.nc - 1);
		}

		#pragma omp for schedule(dynamic)
			// loop over velocity values
			for (int j = 0; j < analysis.nv; j++)
			{
				// print progress on velocity
				printf("Calculating subset %d of %d\n", 
							j + 1, analysis.nv);

				if (analysis.nv == 1)
				{
					velocity = analysis.velocity_min;
				}
				else
				{
					velocity = analysis.velocity_min + 
							(double)(j) * (analysis.velocity_max - 
							analysis.velocity_min) / 
							(double)(analysis.nv - 1);
				}

				y[0] = coordinate;
				y[1] = velocity;

				if (system.dim == 6)
				{
					for (int k = 0; k < 4; k++)
					{
						y[k+2] = orb_ini[k];
					}
				}

				// calculate forward integration
				evolve_orbit(y, &orbit_fw, 
					&orbit_fw_size, system, analysis_fw, rk);

				// calculate backward integration
				evolve_orbit(y, &orbit_bw, 
					&orbit_bw_size, system, analysis_bw, rk);

				#pragma omp critical
				{
					// write initial condition to file
					fprintf(out_ic, "%1.15e %1.15e\n", coordinate, velocity);

					// write orbit and error to file
					for (int k = 0; k < orbit_fw_size; k++)
					{
						// fprintf(out_psp, "%1.15e %1.15e\n", 
						// 	angle_mod_pos(orbit_fw[k][0]), 
						// 	orbit_fw[k][1]);

						fprintf(out_psp, "%1.15e %1.15e\n", 
							angle_mod(orbit_fw[k][0]), 
							orbit_fw[k][1]);

						if (system.dim == 6)
						{
							for (int l = 0; l < 4; l++)
							{
								orb[l] = orbit_fw[k][l+2];
							}
							fprintf(out_orb_ang_mom_err, "%d %1.15e\n", 
									k, fabs(angular_momentum_two_body(orb)-
									angular_momentum_two_body(orb_ini)));
							fprintf(out_vis_viva_err, "%d %1.15e\n", 
									k, fabs(vis_viva_two_body(orb, system)-
									vis_viva_two_body(orb_ini, system)));
						}
					}

					for (int k = 0; k < orbit_bw_size; k++)
					{
						// fprintf(out_psp, "%1.15e %1.15e\n", 
						// 	angle_mod_pos(orbit_bw[k][0]), 
						// 	orbit_bw[k][1]);

						fprintf(out_psp, "%1.15e %1.15e\n", 
							angle_mod(orbit_bw[k][0]), 
							orbit_bw[k][1]);

						if (system.dim == 6)
						{
							for (int l = 0; l < 4; l++)
							{
								orb[l] = orbit_bw[k][l+2];
							}
							fprintf(out_orb_ang_mom_err, "%d %1.15e\n", 
									k, fabs(angular_momentum_two_body(orb)-
									angular_momentum_two_body(orb_ini)));
							fprintf(out_vis_viva_err, "%d %1.15e\n", 
									k, fabs(vis_viva_two_body(orb, system)-
									vis_viva_two_body(orb_ini, system)));
						}
					}
				} // end pragma omp critical

				// free memory
				dealloc_2d_double(&orbit_fw, 
						analysis.number_of_cycles);
				dealloc_2d_double(&orbit_bw, 
						analysis.number_of_cycles);
			
			}
		} // end pragma
		// new line on terminal
		printf("\n");
	}

	// close exit files
	fclose(out_psp);
	fclose(out_ic);
	fclose(out_orb_ang_mom_err);
	fclose(out_vis_viva_err);

	printf("Data written in output folder\n");

	return 0;
}

int time_series(dynsys system,
				anlsis analysis,
				rngkta rk)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/time_series", &st) == -1) {
		mkdir("output/time_series", 0700);
	}

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	// prepare and open exit files 
	FILE *out;
	char	filename[100];

	// indexed file
	sprintf(filename, "output/time_series/time_series_e_%1.3f_K_%1.5f.dat", e, K);
	out = fopen(filename, "w");

	printf("Writting time series with e = %1.3f and K = %1.5f\n", e, K);

	// declare variables
	int orbit_size;
	double res_mean;
	double **orbit;
	double ic[system.dim], orb[4];

	ic[0] = 0.0, ic[1] = 1000.;

	if (system.dim == 6)
	{
		init_orbital(orb, system);
		for (int j = 0; j < 4; j++) ic[j+2] = orb[j];
	}

	// evolve system
	evolve_orbit(ic, &orbit, &orbit_size, system, analysis, rk);

	for (int j = 0; j < orbit_size; j++)
	{
		fprintf(out, "%d %1.15e\n", j, orbit[j][1]);
	}

	fprintf(out, "\n");

	// free memory
	dealloc_2d_double(&orbit, analysis.number_of_cycles);

	// close files
	fclose(out);

	printf("Data written in output/time_series/ folder\n");

	return 0;
}

int multiple_time_series(dynsys system,
						 anlsis analysis,
						 rngkta rk)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/time_series", &st) == -1) {
		mkdir("output/time_series", 0700);
	}

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	// prepare and open exit files 
	FILE *out;
	char	filename[100];

	// indexed file
	sprintf(filename, "output/time_series/multiple_time_series_e_%1.3f_K_%1.5f.dat", e, K);
	out = fopen(filename, "w");

	// declare variables
	int orbit_size;
	double **orbit;
	double ic[system.dim], orb[4];

	#pragma omp parallel private(orbit_size, orbit, ic, orb) //num_threads(THREADS_NUM)
	{

	#pragma omp for schedule(dynamic)
		for (int i = 0; i < 20; i++)
		{
			printf("Calculating set %d of %d\n", i, 20);

			if (i < 10)
			{
				ic[0] = 0.0; ic[1] = 10. + 10. * (double)(i);
			}
			else
			{
				ic[0] = 0.0; ic[1] = 100. + 100. * (double)(i-10);
			}

			if (system.dim == 6)
			{
				init_orbital(orb, system);
				for (int j = 0; j < 4; j++) ic[j+2] = orb[j];
			}

			// evolve system
			evolve_orbit(ic, &orbit, &orbit_size, system, analysis, rk);

			#pragma omp critical
			{
				for (int j = 0; j < orbit_size; j++)
				{
					fprintf(out, "%d %1.15e\n", 
							j, orbit[j][1]);
				}
				fprintf(out, "\n");
			}

			// free memory
			dealloc_2d_double(&orbit, analysis.number_of_cycles);
		
		}

	} // end pragma
	
	// close files
	fclose(out);

	printf("Data written in output/time_series/ folder\n");

	return 0;
}

int multiple_time_series_delta_theta_dot(dynsys system,
										 anlsis analysis,
										 rngkta rk)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/time_series", &st) == -1) {
		mkdir("output/time_series", 0700);
	}

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	// prepare and open exit files 
	FILE *out;
	char	filename[300];

	// indexed file
	sprintf(filename, "output/time_series/multiple_time_series_delta_theta_dot_e_%1.3f_K_%1.5f_nic_%d_delta_%1.2f_system_%s.dat", 
				e, K, analysis.number_of_time_series, analysis.time_series_delta, system.name);
	out = fopen(filename, "w");

	// declare variables
	int orbit_size;
	double **orbit;
	double ic[system.dim], orb[4];

	#pragma omp parallel private(orbit_size, orbit, ic, orb) //num_threads(THREADS_NUM)
	{

	#pragma omp for schedule(dynamic)
		for (int i = -analysis.number_of_time_series/2; i <= analysis.number_of_time_series/2; i++) 
		{
			printf("Calculating time series number %d\n", i + 1 + analysis.number_of_time_series/2);

			ic[0] = 0.0;
			ic[1] = 1000.0 + analysis.time_series_delta * (double)(i);

			if (system.dim == 6)
			{
				init_orbital(orb, system);
				for (int j = 0; j < 4; j++) ic[j+2] = orb[j];
			}

			// evolve system
			evolve_orbit(ic, &orbit, &orbit_size, system, analysis, rk);

			#pragma omp critical
			{
				for (int j = 0; j < orbit_size; j++)
				{
					fprintf(out, "%d %1.15e\n", 
							j, orbit[j][1]);
				}
				fprintf(out, "\n\n");
			}

			// free memory
			dealloc_2d_double(&orbit, analysis.number_of_cycles);
		}

	} // end pragma

	// close files
	fclose(out);

	printf("Data written in output/time_series/ folder\n");

	return 0;
}

int multiple_time_series_delta_theta(dynsys system,
									 anlsis analysis,
									 rngkta rk)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/time_series", &st) == -1) {
		mkdir("output/time_series", 0700);
	}

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	// prepare and open exit files 
	FILE *out;
	char	filename[100];

	// indexed file
	sprintf(filename, "output/time_series/multiple_time_series_delta_theta_e_%1.3f_K_%1.5f_delta_%1.2f.dat", 
				e, K, analysis.time_series_delta);
	out = fopen(filename, "w");

	// declare variables
	int orbit_size;
	double **orbit;
	double ic[system.dim], orb[4];

	#pragma omp parallel private(orbit_size, orbit, ic, orb) //num_threads(THREADS_NUM)
	{

	#pragma omp for schedule(dynamic)
		// for (int i = 0; i < 20; i++)
		// {
		// 	ic[0] = 0.0 + 0.01 * (double)(i); ic[1] = 1000.0;
		for (int i = 0; i <= 20; i++) 
		{
			printf("Calculating time series number %d\n", i + 11);

			ic[0] = 0.0 + analysis.time_series_delta * (double)(i);
			ic[1] = 1000.0;

			if (system.dim == 6)
			{
				init_orbital(orb, system);
				for (int j = 0; j < 4; j++) ic[j+2] = orb[j];
			}

			// evolve system
			evolve_orbit(ic, &orbit, &orbit_size, system, analysis, rk);

			#pragma omp critical
			{
				for (int j = 0; j < orbit_size; j++)
				{
					fprintf(out, "%d %1.15e\n", 
							j, orbit[j][1]);
				}
				fprintf(out, "\n\n");
			}

			// free memory
			dealloc_2d_double(&orbit, analysis.number_of_cycles);
		
		}

	} // end pragma

	// close files
	fclose(out);

	printf("Data written in output/time_series/ folder\n");

	return 0;
}

int evolve_n_cycles_po  (double y0[2],
						 int n,
                         dynsys system,
                         anlsis analysis,
						 rngkta rk)
{
    double t = 0;
	double y[system.dim], orbital[4];
	double *par = (double *)system.params;
	double e = par[1];

	copy (y, y0, 2);

	if (system.dim == 6)
	{
		init_orbital(orbital, system);
		for (int i = 0; i < 4; i++)
		{
			y[i+2] = orbital[i];
		}
	}

    // loop over cycles
	for (int i = 0; i < n; i++)
	{
		evolve_cycle(y, &t, system, analysis, rk);
	
		// check if orbit diverges
		for (int j = 0; j < system.dim; j++)
		{
			if (fabs(y[j]) > analysis.evolve_box_size)
			{
                printf("Warning: evolve periodic orbit failed\n");
				printf("box limit reached\n");
				printf("y[%d] = %1.10e\n", j, y[j]);
				exit(2);
			}
		}
	}

	copy (y0, y, 2);

	return 0;
}

int periodic_orbit	(perorb *po,
                     dynsys system,
                     anlsis analysis,
					 rngkta rk)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/periodic_orbit", &st) == -1)
    {
		mkdir("output/periodic_orbit", 0700);
	}

	FILE    *out_orb, *out_orb_res;
	char    filename[300];
	double 	t = 0;
    double 	y[system.dim], orbital[4];
	double	po_ic_after_one_period[system.dim];
	double 	*par = (double *)system.params;
	double 	gamma = par[0];
	double 	e = par[1];
	double 	K = par[6];
	double	number_of_spins;
	double	one_period_angular_diff;

	(*po).dist_on_phase_space = &dist_from_ref;
	(*po).evolve_n_cycles = &evolve_n_cycles_po;

	if (calculate_periodic_orbit_ic(po, system, analysis, rk) == -1)
	{
		return -1;
	}

    copy(y, (*po).initial_condition, 2);
	if (system.dim == 6)
	{
		init_orbital(orbital, system);
		for (int i = 0; i < 4; i++)
		{
			y[i+2] = orbital[i];
		}
	}

    for (int i = 0; i < (*po).period; i++)
    {
        copy((*po).orbit[i], y, 2);
		evolve_cycle(y, &t, system, analysis, rk);
    }

    // finds the index for the smallest theta value
    // inside the periodic orbit to use as a reference
    // on the file name
    int index_theta_min = 0;
    for (int i = 1; i < (*po).period; i++)
    {
        if(angle_mod((*po).orbit[i][0]) < angle_mod((*po).orbit[i-1][0]))
        {
            index_theta_min = i;
        }
    }

	(*po).initial_condition[0] = angle_mod((*po).orbit[index_theta_min][0]);
	(*po).initial_condition[1] = (*po).orbit[index_theta_min][1];

	// indexed file
	// sprintf(filename, "output/periodic_orbit/periodic_orbit_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_period_%d_ic_%1.3f_%1.3f.dat", 
    //         gamma, e, system.name, K, (*po).period, angle_mod((*po).initial_condition[0]), (*po).initial_condition[1]);
	// out_orb = fopen(filename, "w");

	copy(y, (*po).initial_condition, 2);
	if (system.dim == 6)
	{
		init_orbital(orbital, system);
		for (int i = 0; i < 4; i++)
		{
			y[i+2] = orbital[i];
		}
	}

    for (int i = 0; i < (*po).period; i++)
    {
        copy((*po).orbit[i], y, system.dim);

		// fprintf(out_orb, "%1.10e %1.10e\n", 
        //     angle_mod((*po).orbit[i][0]), (*po).orbit[i][1]);

		evolve_cycle(y, &t, system, analysis, rk);
    }

	copy(po_ic_after_one_period, y, system.dim);

	// sprintf(filename, "output/periodic_orbit/periodic_orbit_resonance_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_period_%d_ic_%1.3f_%1.3f.dat", 
    //         gamma, e, system.name, K, (*po).period, angle_mod((*po).initial_condition[0]), (*po).initial_condition[1]);
	// out_orb_res = fopen(filename, "w");

	// calculating the resonance
	one_period_angular_diff = 
		po_ic_after_one_period[0] - (*po).initial_condition[0];

	number_of_spins = one_period_angular_diff / (2.*M_PI);

	(*po).winding_number_numerator = (int) round(number_of_spins);
	(*po).winding_number_denominator = (*po).period;

	// fprintf(out_orb_res, "Orbit:\n\n");

	// for (int i = 0; i < (*po).period; i++)
    // {
    //     fprintf(out_orb_res, "%1.5e %1.5e\n", 
    //         angle_mod((*po).orbit[i][0]), (*po).orbit[i][1]);
    // }

	// fprintf(out_orb_res, "\n\n");

	// fprintf(out_orb_res, "Orbit period:\n\n%d\n\n\n", 
	// 	(*po).period);

	// fprintf(out_orb_res, "Angular difference after 1 period:\n\n%1.10e\n\n\n", 
	// 	one_period_angular_diff);

	// fprintf(out_orb_res, "Number of spins:\n\n%1.10e\n\n\n", 
	// 	number_of_spins);

	// fprintf(out_orb_res, "Resonance:\n\n%d / %d\n\n\n", 
	// 	(*po).winding_number_numerator, 
	// 	(*po).winding_number_denominator);

    // fclose(out_orb);
	// fclose(out_orb_res);

	return 0;
}

int look_for_resonance	(int number_of_candidates,
						 double candidates[][2],
						 int spin_period,
                         int orbit_period,
                         dynsys system, 
						 anlsis analysis,
						 rngkta rk)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/periodic_orbit", &st) == -1) {
		mkdir("output/periodic_orbit", 0700);
	}

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	// prepare and open exit files
	FILE	*out_dist, *out_cand;
	char	filename[300];

	// sprintf(filename, "output/periodic_orbit/resonance_distances_gamma_%1.6f_e_%1.3f_spin_orbit_%d_%d.dat", 
	// 	gamma, e, spin_period, orbit_period);
	// out_dist = fopen(filename, "w");
	// sprintf(filename, "output/periodic_orbit/resonance_candidates_gamma_%1.6f_e_%1.3f_spin_orbit_%d_%d.dat", 
	// 	gamma, e, spin_period, orbit_period);
	// out_cand = fopen(filename, "w");

	// declare variables
	double y[system.dim], y0[system.dim];
	double rot_ini[2], orb_ini[4];
	init_orbital(orb_ini, system);
	int grid[2];
	double ic[2];
	double t;

	int safety_factor = 5;	// looks for a multiple of the p.o instead of the p.o itself

	int **candidates_ij;
	alloc_2d_int(&candidates_ij, number_of_candidates, 2);

	double **orbit_distance, **spin_distance;

	alloc_2d_double(&orbit_distance, 
		analysis.grid_resolution, analysis.grid_resolution);
	alloc_2d_double(&spin_distance, 
		analysis.grid_resolution, analysis.grid_resolution);
	for (int i = 0; i < analysis.grid_resolution; i++)
	{
		for (int j = 0; j < analysis.grid_resolution; j++)
		{
			orbit_distance[i][j] = NAN;
			spin_distance[i][j] = NAN;
		}
	}

	// loop over coordinate values
	for (int i = 0; i < analysis.grid_resolution; i++)
	{
		// print progress on coordinate
		printf("Calculating set %d of %d\n", 
					i + 1, analysis.grid_resolution);

		#pragma omp parallel private(y, y0, grid, rot_ini, t) shared(orbit_distance, spin_distance) //num_threads(THREADS_NUM)
		{

		#pragma omp for schedule(dynamic)
			// loop over velocity values
			for (int j = 0; j < analysis.grid_resolution; j++)
			{
				// print progress on velocity
				// printf("Calculating subset %d of %d\n", 
				// 			j + 1, analysis.grid_resolution);

				grid[0] = i;
				grid[1] = j;

				grid_to_double(grid, rot_ini, analysis);

				copy(y, rot_ini, 2);
				if (system.dim == 6)
				{
					for (int k = 0; k < 4; k++)
					{
						y[k+2] = orb_ini[k];
					}
				}

				t = 0.0;
				copy(y0, y, system.dim);
				for (int i = 0; i < safety_factor * orbit_period; i++)
				{
					evolve_cycle(y, &t, system, analysis, rk);
				}

				orbit_distance[i][j] = dist_from_ref(y,y0);

				spin_distance[i][j] = 
					fabs((y[0] - y0[0]) / (2.*M_PI) - (double) (safety_factor * spin_period));
			}
		} // end pragma

		// new line on terminal
		// printf("\n");
	}

	// for (int i = 0; i < analysis.grid_resolution; i++)
	// {
	// 	for (int j = 0; j < analysis.grid_resolution; j++)
	// 	{
	// 		if (orbit_distance[i][j] != orbit_distance[i][j] ||
	// 			 spin_distance[i][j] != spin_distance[i][j])
	// 		{
	// 			printf("Looking for the resonance gave NAN value.\n");
	// 		}

	// 		grid[0] = i;
	// 		grid[1] = j;
	// 		grid_to_double(grid, ic, analysis);

	// 		fprintf(out_dist, "%1.10f %1.10f %1.10e %1.10e\n", 
	// 			ic[0], ic[1], orbit_distance[i][j],	spin_distance[i][j]);
	// 	}
	// 	fprintf(out_dist, "\n");
	// }

	int counter = 0;
	for (int i = 0; i < analysis.grid_resolution; i++)
	{
		for (int j = 0; j < analysis.grid_resolution; j++)
		{
			if (counter < number_of_candidates)
			{
				candidates_ij[counter][0] = i;
				candidates_ij[counter][1] = j;
				counter++;
			}
			else
			{
				for (int k = 0; k < number_of_candidates - 1; k++)
				{
					for (int l = k + 1; l < number_of_candidates; l++)
					{
						if ((spin_distance[candidates_ij[k][0]][candidates_ij[k][1]] < 
							 spin_distance[candidates_ij[l][0]][candidates_ij[l][1]]) &&
							(orbit_distance[candidates_ij[k][0]][candidates_ij[k][1]] < 
							 orbit_distance[candidates_ij[l][0]][candidates_ij[l][1]]))
						{
							int hold[2];
							copy_int(hold, candidates_ij[k], 2);
							copy_int(candidates_ij[k], candidates_ij[l], 2);
							copy_int(candidates_ij[l], hold, 2);
						}
					}
				}
				for (int k = 0; k < number_of_candidates; k++)
				{
					if ((spin_distance[i][j] < spin_distance[candidates_ij[k][0]][candidates_ij[k][1]]) &&
						(orbit_distance[i][j] < orbit_distance[candidates_ij[k][0]][candidates_ij[k][1]]))
					{
						candidates_ij[k][0] = i;
						candidates_ij[k][1] = j;
						break;
					}
				}
			}
		}
	}

	for (int i = 0; i < number_of_candidates; i++)
	{
		grid[0] = candidates_ij[i][0];
		grid[1] = candidates_ij[i][1];
		grid_to_double(grid, ic, analysis);
		// fprintf(out_cand, "%1.5f %1.5f\n", ic[0], ic[1]);
		copy(candidates[i], ic, 2);
	}
				
	dealloc_2d_int(&candidates_ij, number_of_candidates);

	dealloc_2d_double(&orbit_distance, 
					analysis.grid_resolution);
	dealloc_2d_double(&spin_distance, 
					analysis.grid_resolution);

	// close exit files
	// fclose(out_dist);
	// fclose(out_cand);

	// printf("Data written in output/periodic_orbit/\n");

	return 0;
}

int find_all_periodic_attractors(int *number_of_pos,
							 	 perorb **multiple_pos,
							 	 dynsys system,
                         	 	 anlsis analysis,
								 rngkta rk)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/periodic_orbit", &st) == -1) {
		mkdir("output/periodic_orbit", 0700);
	}
	
	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	// prepare and open exit files
	FILE	*out;
	char	filename[300];

	sprintf(filename, "output/periodic_orbit/all_periodic_orbits_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f.dat", 
		gamma, e, system.name, K);

	int spin_period, orbit_period;
	int number_of_candidates;
	double tol_same_po, tol_between_pos;

	// distance value for which we say a periodic orbit of period n
	// has actually period m where n is a multiple of m
	tol_same_po = 1e-6;

	// distance value for which we say two pos are actually the same one
	tol_between_pos = 1e-3;		

	*number_of_pos = 0;

	for (orbit_period = analysis.orbit_period_min; orbit_period <= analysis.orbit_period_max; orbit_period++)
	{
		for (spin_period = analysis.spin_period_min; spin_period <= analysis.spin_period_max; spin_period++)
		{
			printf("SPIN = %d ORBIT = %d\n", spin_period, orbit_period);
			number_of_candidates = 10 * orbit_period; // 2 * orbit_period
			
			double candidates[number_of_candidates][2];
			look_for_resonance (number_of_candidates, candidates, 
				spin_period, orbit_period, system, analysis, rk);

			perorb multiple_candidates[number_of_candidates];
			int successful_candidates_indices[number_of_candidates];
			for (int i = 0; i < number_of_candidates; i++)
			{
				successful_candidates_indices[i] = -1;
			}

			for (int i = 0; i < number_of_candidates; i++)
			{
				multiple_candidates[i].period = orbit_period;
				copy (multiple_candidates[i].seed, candidates[i], 2);
				alloc_2d_double(&multiple_candidates[i].orbit, multiple_candidates[i].period, system.dim);
				if (periodic_orbit(&multiple_candidates[i], system, analysis, rk) == 0)
				{
					// check if spin exceeds maximum defined spin
					if (multiple_candidates[i].winding_number_numerator > analysis.spin_period_max)
					{
						goto skip;
					}
					// check if p.o. already registered
					for (int k = 0; k < *number_of_pos; k++)
					{
						for (int l = 0; l < (*multiple_pos)[k].period; l++)
						{
							for (int m = 0; m < orbit_period; m++)
							{
								if (dist_from_ref((*multiple_pos)[k].orbit[l], multiple_candidates[i].orbit[m]) 
									< tol_between_pos)
								{
									goto skip;
								}
							}
						}
					}
					// check if p.o is stable
					jacobian_eigenvalues_magnitude(&multiple_candidates[i], system, analysis, rk);
					if ((multiple_candidates[i].eigenvalues_absolute_value[0] > 1.0) ||
						(multiple_candidates[i].eigenvalues_absolute_value[1] > 1.0))
					{
						goto skip;
					}
					// check if orbit period is actually lower
					int number_of_similar_points = 1;
					for(int k = 0; k < multiple_candidates[i].period - 1; k++)
					{
						for(int l = k + 1; l < multiple_candidates[i].period; l++)
						{
							if (dist_from_ref(multiple_candidates[i].orbit[k], multiple_candidates[i].orbit[l]) 
									< tol_same_po)
							{
								number_of_similar_points++;
							}
						}
						if (number_of_similar_points > 1)
						{
							break;
						}
					}
					successful_candidates_indices[i] = 1;
					*number_of_pos = *number_of_pos + 1;
					if (*number_of_pos == 1)
					{
						*multiple_pos = (perorb*) malloc(*number_of_pos * sizeof(perorb));
					}
					else
					{
						*multiple_pos = realloc(*multiple_pos, *number_of_pos * sizeof(perorb));
					}

					(*multiple_pos)[*number_of_pos-1].period = multiple_candidates[i].period/number_of_similar_points;
					(*multiple_pos)[*number_of_pos-1].seed[0] = multiple_candidates[i].seed[0];
					(*multiple_pos)[*number_of_pos-1].seed[1] = multiple_candidates[i].seed[1];
					(*multiple_pos)[*number_of_pos-1].initial_condition[0] = multiple_candidates[i].initial_condition[0];
					(*multiple_pos)[*number_of_pos-1].initial_condition[1] = multiple_candidates[i].initial_condition[1];
					(*multiple_pos)[*number_of_pos-1].eigenvalues_absolute_value[0] = pow(multiple_candidates[i].eigenvalues_absolute_value[0],1.0/number_of_similar_points);
					(*multiple_pos)[*number_of_pos-1].eigenvalues_absolute_value[1] = pow(multiple_candidates[i].eigenvalues_absolute_value[1],1.0/number_of_similar_points);
					(*multiple_pos)[*number_of_pos-1].winding_number_numerator = multiple_candidates[i].winding_number_numerator/number_of_similar_points;
					(*multiple_pos)[*number_of_pos-1].winding_number_denominator = multiple_candidates[i].winding_number_denominator/number_of_similar_points;

					alloc_2d_double(&(*multiple_pos)[*number_of_pos-1].orbit, (*multiple_pos)[*number_of_pos-1].period, system.dim);
					for (int j = 0; j < (*multiple_pos)[*number_of_pos-1].period; j++)
					{
						copy((*multiple_pos)[*number_of_pos-1].orbit[j], multiple_candidates[i].orbit[j], system.dim);
					}
				}
				skip:;
				dealloc_2d_double(&multiple_candidates[i].orbit, multiple_candidates[i].period);
			}
		}
	}

	if (*number_of_pos == 0)
	{
		printf("Could not find any po.\n");
	}
	else
	{
		out = fopen(filename, "w");
		for (int j = 0; j < *number_of_pos; j++)
		{
			for (int k = 0; k < (*multiple_pos)[j].period; k++)
			{
				fprintf(out, "%1.15e %1.15e %d %d\n", 
					angle_mod((*multiple_pos)[j].orbit[k][0]),
					(*multiple_pos)[j].orbit[k][1],
					(*multiple_pos)[j].winding_number_numerator,
					(*multiple_pos)[j].winding_number_denominator);
			}
			fprintf(out, "\n");
		}
		fclose(out);
	}

	return 0;
}

int fill_attractor_array(int *number_of_pos,
						 perorb **multiple_pos,
                         dynsys system,
                         anlsis analysis,
						 rngkta rk)
{
	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	FILE	*in;
	char	filename[300];

	sprintf(filename, "output/periodic_orbit/all_periodic_orbits_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f.dat", 
		gamma, e, system.name, K);

	// sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_ref_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat", 
	// 	gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	
	in = fopen(filename, "r");
	if (in == NULL)
	{
		find_all_periodic_attractors(number_of_pos, multiple_pos, system, analysis, rk);
	}
	else
	{
		int 	wind_x, wind_y;
		double	po_x, po_y;
	
		printf("Getting pos from printed file\n");

		*number_of_pos = 0;
		while(fscanf(in, "%lf %lf %d %d", &po_x, &po_y, &wind_x, &wind_y) != EOF)
		{
			*number_of_pos = *number_of_pos + 1;
			
			if (*number_of_pos == 1)
			{
				*multiple_pos = (perorb*) malloc(*number_of_pos * sizeof(perorb));
			}
			else
			{
				*multiple_pos = realloc(*multiple_pos, *number_of_pos * sizeof(perorb));
			}
			
			(*multiple_pos)[*number_of_pos-1].period = wind_y;
			(*multiple_pos)[*number_of_pos-1].seed[0] = po_x;
			(*multiple_pos)[*number_of_pos-1].seed[1] = po_y;
			(*multiple_pos)[*number_of_pos-1].initial_condition[0] = po_x;
			(*multiple_pos)[*number_of_pos-1].initial_condition[1] = po_y;
			(*multiple_pos)[*number_of_pos-1].winding_number_numerator = wind_x;
			(*multiple_pos)[*number_of_pos-1].winding_number_denominator = wind_y;
			
			alloc_2d_double(&(*multiple_pos)[*number_of_pos-1].orbit, (*multiple_pos)[*number_of_pos-1].period, 2);
			(*multiple_pos)[*number_of_pos-1].orbit[0][0] = po_x;
			(*multiple_pos)[*number_of_pos-1].orbit[0][1] = po_y;
			for (int l = 1; l < (*multiple_pos)[*number_of_pos-1].period; l++)
			{
				fscanf(in, "%lf %lf %d %d", &po_x, &po_y, &wind_x, &wind_y);
				(*multiple_pos)[*number_of_pos-1].orbit[l][0] = po_x;
				(*multiple_pos)[*number_of_pos-1].orbit[l][1] = po_y;
			}
		}
		fclose(in);
		printf("Done\n");
	}
	return 0;
}

int fill_basin_matrix	(double	***basin_matrix,
						 dynsys system,
                         anlsis analysis)
{
	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	FILE	*in;
	char	filename[300];

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	
	in = fopen(filename, "r");
	if (in == NULL)
	{
		printf("Error, basin file not found.\n");
		exit(4);
	}
	else
	{
		alloc_2d_double(basin_matrix, 
			analysis.grid_resolution, analysis.grid_resolution);

		int 	fool_spin, fool_orbit, i, j, fool_control;
		double	fool_basin_x, fool_basin_y, cantor, fool_time;
	
		printf("Getting basin matrix from printed file\n");

		while(fscanf(in, "%lf %lf %lf %lf %d %d %d %d %d\n", 
			&fool_basin_x, &fool_basin_y, &cantor, &fool_time, &fool_spin, 
			&fool_orbit, &i, &j, &fool_control) != EOF)
		{
			(*basin_matrix)[i][j] = cantor;
		}
		fclose(in);
	}
	
	return 0;
}

int fill_control_matrix	(int	***control_matrix,
						 dynsys system,
                         anlsis analysis)
{
	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	FILE	*in;
	char	filename[300];

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);

	in = fopen(filename, "r");
	if (in == NULL)
	{
		printf("Error, basin file not found.\n");
		exit(4);
	}
	else
	{
		alloc_2d_int(control_matrix, 
			analysis.grid_resolution, analysis.grid_resolution);

		int 	fool_spin, fool_orbit, i, j, control;
		double	fool_basin_x, fool_basin_y, fool_cantor, fool_time;
	
		printf("Getting control matrix from printed file\n");

		while(fscanf(in, "%lf %lf %lf %lf %d %d %d %d %d\n", 
			&fool_basin_x, &fool_basin_y, &fool_cantor, &fool_time, &fool_spin, 
			&fool_orbit, &i, &j, &control) != EOF)
		{
			(*control_matrix)[i][j] = control;
		}
		fclose(in);
	}

	return 0;
}

int fill_control_monte_carlo(int	**control,
						 	 dynsys system,
                         	 anlsis analysis)
{
	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	FILE	*in;
	char	filename[300];

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_control_monte_carlo_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc);

	in = fopen(filename, "r");
	if (in == NULL)
	{
		printf("Error, control monte carlo file not found.\n");
		exit(4);
	}
	else
	{
		alloc_1d_int(control, analysis.number_of_rand_orbits_mc);

		int	i, control_i;
	
		printf("Getting control monte carlo from printed file\n");

		while(fscanf(in, "%d %d\n", &i, &control_i) != EOF)
		{
			(*control)[i] = control_i;
		}
		fclose(in);
	}

	return 0;
}

int fill_control_monte_carlo_with_break(int *control_size,
										int	**control,
						 	 			dynsys system,
                         	 			anlsis analysis)
{
	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	FILE	*in;
	char	filename[300];

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_control_monte_carlo_with_break_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc);

	in = fopen(filename, "r");
	if (in == NULL)
	{
		printf("Error, control monte carlo file not found.\n");
		exit(4);
	}
	else
	{
		alloc_1d_int(control, analysis.number_of_rand_orbits_mc);

		int	i, control_i;
	
		printf("Getting control monte carlo from printed file\n");

		while(fscanf(in, "%d %d\n", &i, &control_i) != EOF)
		{
			(*control)[i] = control_i;
		}
		*control_size = i + 1;
		fclose(in);
	}

	return 0;
}

double dist_from_ref(double x[2],
                    double ref[2])
{
	double dist_from_ref_theta = 
			angular_dist(x[0], ref[0]);

    double dist_from_ref_theta_dot = 
            fabs(x[1]-ref[1]);

    double dist_ref = 
        sqrt(dist_from_ref_theta*dist_from_ref_theta
            +dist_from_ref_theta_dot*dist_from_ref_theta_dot);

    return dist_ref;
}

int evolve_basin(double *ic,
				 bool *converged,
                 int *convergence_time,
                 perorb po,
                 dynsys system,
				 anlsis analysis,
				 rngkta rk)
{
	// declare variables
	bool is_close_to;
	int orbit_counter, close_time_counter;
	double y[system.dim], rot[2];
	double t = 0.0;

	// takes into consideration initial condition
	orbit_counter = 1;

	// starts proximity counter
	close_time_counter = 0;
	
	*converged = false;

	copy(y, ic, system.dim);
	for (int i = 0; i < analysis.number_of_cycles; i++)
	{
		copy(rot, y, 2);

		is_close_to = false;
		for (int j = 0; j < po.period; j++)
		{
			if(dist_from_ref(rot, po.orbit[j]) < analysis.evolve_basin_eps)
			{
				is_close_to = true;
			}
		}
		
		if (is_close_to == true)
		{
			close_time_counter++;
		}
		else
		{
			close_time_counter = 0;
		}

		if (close_time_counter > analysis.evolve_basin_time_tol)
		{
			*converged = true;
			goto out;
		}

		evolve_cycle(y, &t, system, analysis, rk);
	
		// check if orbit diverges
		for (int j = 0; j < system.dim; j++)
		{
			if (fabs(y[j]) > analysis.evolve_box_size)
			{
				printf("Warning: box limit reached\n");
				printf("y[%d] = %1.10e\n", j, y[j]);
				goto out;
			}
		}

		orbit_counter++;

	}

	out:;

	if (*converged == true)
	{
		*convergence_time = orbit_counter;
	}
	else
	{
		*convergence_time = -1;
	}

	return 0;
}

int basin_of_attraction (perorb po,
                         dynsys system,
                         anlsis analysis,
						 rngkta rk)
{
	// little warning
	if (strcmp(system.name, "two_body") == 0)
	{
		printf("Warning: cant draw basins\n");
		printf("for two-body system\n");
		exit(2);
	}

	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/basin_of_attraction", &st) == -1) {
		mkdir("output/basin_of_attraction", 0700);
	}

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	// prepare and open exit files
	FILE	*out_boa, *out_ref, *out_size;
	char	filename[300];

	sprintf(filename, "output/basin_of_attraction/basin_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat", 
		gamma, e, system.name, K, po.initial_condition[0], po.initial_condition[1], po.period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	out_boa = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/basin_ref_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat", 
		gamma, e, system.name, K, po.initial_condition[0], po.initial_condition[1], po.period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	out_ref = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/basin_size_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat", 
		gamma, e, system.name, K, po.initial_condition[0], po.initial_condition[1], po.period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	out_size = fopen(filename, "w");

	// declare variables
	double y[system.dim];
	double coordinate, velocity;
	double rot_ini[2];
	double orb[4], orb_ini[4];
	init_orbital(orb_ini, system);
	int basin_size = 0;
	double basin_size_fraction;
	int grid[2];
	double basin[2];
	bool converged;
	int convergence_time;

	for (int i = 0; i < po.period; i++)
	{
		fprintf(out_ref, "%1.15e %1.15e\n", po.orbit[i][0], po.orbit[i][1]);
	}

	int **basin_matrix, **control_matrix;
	double **time_matrix;

	alloc_2d_int(&basin_matrix, 
		analysis.grid_resolution, analysis.grid_resolution);
	alloc_2d_int(&control_matrix, 
		analysis.grid_resolution, analysis.grid_resolution);
	alloc_2d_double(&time_matrix, 
		analysis.grid_resolution, analysis.grid_resolution);
	for (int i = 0; i < analysis.grid_resolution; i++)
	{
		for (int j = 0; j < analysis.grid_resolution; j++)
		{
			basin_matrix[i][j] = 0;
			control_matrix[i][j] = 0;
			time_matrix[i][j] = NAN;
		}
	}

	// loop over coordinate values
	for (int i = 0; i < analysis.grid_resolution; i++)
	{
		// print progress on coordinate
		printf("Calculating set %d of %d\n", 
					i + 1, analysis.grid_resolution);

		#pragma omp parallel private(y, coordinate, velocity, basin, grid, \
				converged, convergence_time, orb, rot_ini) shared(basin_matrix, \
				control_matrix, time_matrix) //num_threads(THREADS_NUM)
		{

		#pragma omp for schedule(dynamic)
			// loop over velocity values
			for (int j = 0; j < analysis.grid_resolution; j++)
			{
				// print progress on velocity
				// printf("Calculating subset %d of %d\n", 
				// 			j + 1, analysis.grid_resolution);

				if (control_matrix[i][j] == 0)
				{
					grid[0] = i;
					grid[1] = j;

					grid_to_double(grid, rot_ini, analysis);

					copy(y, rot_ini, 2);

					if (system.dim == 6)
					{
						for (int k = 0; k < 4; k++)
						{
							y[k+2] = orb_ini[k];
						}
					}

					// calculate forward integration
					evolve_basin(y, &converged, &convergence_time,
						po, system, analysis, rk);

					if(converged == true)
					{
						basin_matrix[i][j] = 1;
						control_matrix[i][j] = 1;
						time_matrix[i][j] = (double)(convergence_time);
						basin_size++;
					}
					else
					{
						control_matrix[i][j] = -1;
					}

				}
			
			}
		} // end pragma

		// new line on terminal
		// printf("\n");
	}

	for (int i = 0; i < analysis.grid_resolution; i++)
	{
		for (int j = 0; j < analysis.grid_resolution; j++)
		{
			grid[0] = i;
			grid[1] = j;
			grid_to_double(grid, basin, analysis);

			fprintf(out_boa, "%1.5f %1.5f %d %1.5e\n", 
				basin[0], basin[1], 
				basin_matrix[grid[0]][grid[1]],
				time_matrix[grid[0]][grid[1]]);
		}
		fprintf(out_boa, "\n");
	}

	basin_size_fraction = ((double)basin_size) / 
		(((double)analysis.grid_resolution) * ((double)analysis.grid_resolution));

	printf("Basin size = %d\n", basin_size);
	printf("Basin size fraction = %1.4f\n", basin_size_fraction);
	fprintf(out_size, "%1.4f", basin_size_fraction);

	// close exit files
	fclose(out_boa);
	fclose(out_ref);
	fclose(out_size);

	dealloc_2d_int(&basin_matrix, 
					analysis.grid_resolution);

	dealloc_2d_int(&control_matrix, 
					analysis.grid_resolution);

	dealloc_2d_double(&time_matrix, 
					analysis.grid_resolution);

	printf("Data written in output/basin_of_attraction/\n");

	return 0;
}

int evolve_multiple_basin_determined(double *ic,
									 int number_of_po,
									 int *converged_po_id,
									 int *convergence_time,
									 perorb po[],
									 dynsys system,
									 anlsis analysis,
									 rngkta rk)
{
	// declare variables
	bool 	is_close_to_po;
	int 	orbit_counter, close_time_counter;
	int 	internal_converged_po_id;
	double 	y[system.dim], rot[2];
	double 	t = 0.0;

	bool	check_winding_number = true;
	double 	y_ref[system.dim];
	double 	winding_number_window[analysis.convergence_window_wn];

	double *angle_progress;
	alloc_1d_double(&angle_progress, 1);

	// takes into consideration initial condition
	orbit_counter = 1;

	// starts proximity counters
	close_time_counter = 0;
	
	*converged_po_id = -1;

	copy(y, ic, system.dim);
	for (int i = 0; i < analysis.number_of_cycles; i++)
	{
		copy(rot, y, 2);

		// check if orbit is close to an indexed po
		is_close_to_po = false;
		for (int j = 0; j < number_of_po; j++)
		{
			for (int k = 0; k < po[j].period; k++)
			{
				if(dist_from_ref(rot, po[j].orbit[k]) < analysis.evolve_basin_eps)
				{
					is_close_to_po = true;
					internal_converged_po_id = j;
					goto stop_checking_po;
				}
			}
		}
		stop_checking_po:;
		if (is_close_to_po == true)
		{
			close_time_counter++;
		}
		else
		{
			close_time_counter = 0;
		}
		if (close_time_counter > analysis.evolve_basin_time_tol)
		{
			*converged_po_id = internal_converged_po_id;
			goto out;
		}

		// check if orbit is close to a limit cycle by checking
		// if its rotation number converged
		if (check_winding_number == true)
		{
			int i_winding = i - (analysis.convergence_transient_wn + 1);
			if (i == analysis.convergence_transient_wn) // i_winding = -1
			{
				copy(y_ref, y, system.dim);
				angle_progress[0] = y[0];
			}
			else if (i > analysis.convergence_transient_wn)	// i_winding >= 0
			{
				double 	delta_theta = y[0] - y_ref[0];
				double 	wn = delta_theta / (double) (i_winding + 1);

				angle_progress = (double*) realloc(angle_progress, (i_winding + 2) * sizeof(double));
				angle_progress[i_winding + 1] = y[0];

				if(i_winding < analysis.convergence_window_wn)
				{
					winding_number_window[i_winding] = wn;
				}
				else if (i_winding >= analysis.convergence_window_wn)
				{
					double max_winding_number = wn;
					double min_winding_number = wn;
					for (int j = 0; j < analysis.convergence_window_wn - 1; j++)
					{
						winding_number_window[j] = winding_number_window[j+1];
						if(winding_number_window[j] > max_winding_number)
						{
							max_winding_number = winding_number_window[j];
						}
						if(winding_number_window[j] < min_winding_number)
						{
							min_winding_number = winding_number_window[j];
						}
					}
					winding_number_window[analysis.convergence_window_wn - 1] = wn;
					if (fabs(max_winding_number-min_winding_number) < analysis.convergence_precision_wn)
					{
						double pwn_numerator = 0.0, pwn_denominator = 0.0;
						for (int l = 1; l <= i_winding; l++)
						{
							double twp = (((double) l) / ((double) (i_winding + 1)));
							double factor = 1.0 / exp(1.0 / (twp * (1.0 - twp)));
							pwn_numerator += (angle_progress[l] - angle_progress[l-1]) * factor;
							pwn_denominator += factor;
						}
						double precise_wn = pwn_numerator / pwn_denominator;

						wn = precise_wn;

						double 	wn_mod = angle_mod_pos(wn);
						if(fabs(wn_mod) < 5e-2) wn_mod = 2.0*M_PI;
						double 	dist_from_int = fabs((2.0*M_PI/wn_mod) - round(2.0*M_PI/wn_mod));
						if(dist_from_int > 5e-2)
						{
							*converged_po_id = -2;
							goto out; // irrational winding number found
						}
					}
				}
			}
		}

		evolve_cycle(y, &t, system, analysis, rk);

		// check if orbit diverges
		for (int j = 0; j < system.dim; j++)
		{
			if (fabs(y[j]) > analysis.evolve_box_size)
			{
				printf("Warning: box limit reached\n");
				printf("y[%d] = %1.10e\n", j, y[j]);
				goto out;
			}
		}

		orbit_counter++;

	}

	out:;

	if (*converged_po_id == -1)
	{
		*convergence_time = -1;
	}
	else if (*converged_po_id == -2)
	{
		*convergence_time = orbit_counter - analysis.convergence_window_wn;
	}
	else
	{
		*convergence_time = orbit_counter - analysis.evolve_basin_time_tol;
	}

	dealloc_1d_double(&angle_progress);

	return 0;
}

int multiple_basin_of_attraction_determined (int number_of_po,
											 perorb po[],
                         					 dynsys system,
                         					 anlsis analysis,
											 rngkta rk)
{
	// little warning
	if (strcmp(system.name, "two_body") == 0)
	{
		printf("Warning: cant draw basins\n");
		printf("for two-body system\n");
		exit(2);
	}

	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/basin_of_attraction", &st) == -1) {
		mkdir("output/basin_of_attraction", 0700);
	}

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	// prepare and open exit files
	FILE	*out_boa, *out_ref;
	char	filename[300];

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	out_boa = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_ref_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	out_ref = fopen(filename, "w");

	// declare variables
	double y[system.dim];
	double coordinate, velocity;
	double rot_ini[2];
	double orb[4];

	int grid[2];
	double basin[2];
	
	int converged_po_id;
	int convergence_time;

	for (int i = 0; i < number_of_po; i++)
	{
		for (int j = 0; j < po[i].period; j++)
		{
			fprintf(out_ref, "%1.15e %1.15e %d %d\n",
				angle_mod(po[i].orbit[j][0]), 
				po[i].orbit[j][1],
				po[i].winding_number_numerator,
				po[i].winding_number_denominator);
		}
		fprintf(out_ref, "\n");
	}
	fclose(out_ref);

	int **control_matrix;
	double **basin_matrix, **time_matrix;

	alloc_2d_int(&control_matrix, 
		analysis.grid_resolution, analysis.grid_resolution);
	alloc_2d_double(&basin_matrix, 
		analysis.grid_resolution, analysis.grid_resolution);
	alloc_2d_double(&time_matrix, 
		analysis.grid_resolution, analysis.grid_resolution);
	for (int i = 0; i < analysis.grid_resolution; i++)
	{
		for (int j = 0; j < analysis.grid_resolution; j++)
		{
			control_matrix[i][j] = 0;
			basin_matrix[i][j] = NAN;
			time_matrix[i][j] = NAN;
		}
	}

	// loop over coordinate values
	for (int i = 0; i < analysis.grid_resolution; i++)
	{
		// print progress on coordinate
		printf("Calculating set %d of %d\n", 
					i + 1, analysis.grid_resolution);

		#pragma omp parallel private(y, coordinate, velocity, basin, grid, \
				converged_po_id, convergence_time, orb, rot_ini) shared(basin_matrix, \
				control_matrix, time_matrix) //num_threads(THREADS_NUM)
		{

		#pragma omp for schedule(dynamic)
			// loop over velocity values
			for (int j = 0; j < analysis.grid_resolution; j++)
			{
				if (control_matrix[i][j] == 0)
				{
					grid[0] = i;
					grid[1] = j;

					// grid_to_double(grid, rot_ini, analysis);
					grid_to_double_v2(grid, rot_ini, analysis);

					copy(y, rot_ini, 2);
					complete_orbital_part(y, system);

					// calculate forward integration
					evolve_multiple_basin_determined(y, number_of_po, 
						&converged_po_id, &convergence_time,
						po, system, analysis, rk);

					if (converged_po_id == -1)
					{
						control_matrix[i][j] = -1;
						printf("Orbit did not converge\n");
						// exit(42);
					}
					else if (converged_po_id == -2)
					{
						control_matrix[i][j] = -2;
						basin_matrix[i][j] = 0.0;
						time_matrix[i][j] = (double)(convergence_time);
					}
					else
					{
						basin_matrix[i][j] = 
							(double) cantor_pairing_function(po[converged_po_id].winding_number_numerator,po[converged_po_id].winding_number_denominator);
						control_matrix[i][j] = converged_po_id + 1;
						time_matrix[i][j] = (double)(convergence_time);
					}

				}
			
			}
		} // end pragma

		// new line on terminal
		// printf("\n");
	}

	int spin_period, orbit_period;
	for (int i = 0; i < analysis.grid_resolution; i++)
	{
		for (int j = 0; j < analysis.grid_resolution; j++)
		{
			grid[0] = i;
			grid[1] = j;
			// grid_to_double(grid, basin, analysis);
			grid_to_double_v2(grid, basin, analysis);

			if (control_matrix[i][j] != -1)
			{
				spin_period = po[control_matrix[i][j]-1].winding_number_numerator;
				orbit_period = po[control_matrix[i][j]-1].winding_number_denominator;
			}
			else
			{
				spin_period = 0;
				orbit_period = 0;
			}

			fprintf(out_boa, "%1.5f %1.5f %f %1.5e %d %d %d %d %d\n", 
				basin[0], basin[1], 
				basin_matrix[grid[0]][grid[1]],
				time_matrix[grid[0]][grid[1]],
				spin_period,
				orbit_period,
				grid[0],
				grid[1],
				control_matrix[i][j]);
		}
		fprintf(out_boa, "\n");
	}
	fclose(out_boa);

	// free memory
	dealloc_2d_int(&control_matrix, 
					analysis.grid_resolution);
	dealloc_2d_double(&basin_matrix, 
					analysis.grid_resolution);
	dealloc_2d_double(&time_matrix, 
					analysis.grid_resolution);

	printf("Data written in output/basin_of_attraction/\n");

	return 0;
}

int multiple_basin_of_attraction_determined_monte_carlo	(int number_of_po,
											 			 perorb po[],
                         					 			 dynsys system,
                         					 			 anlsis analysis,
														 rngkta rk)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/basin_of_attraction", &st) == -1) {
		mkdir("output/basin_of_attraction", 0700);
	}

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	// prepare and open exit files
	FILE	*out_control, *out_ref, *out_converged_time;
	char	filename[300];

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_control_monte_carlo_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc);
	out_control = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_ref_monte_carlo_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc);
	out_ref = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_converged_time_monte_carlo_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc);
	out_converged_time = fopen(filename, "w");

	// declare variables
	double y[system.dim];
	int converged_po_id, convergence_time;
	int *basin_size;
	alloc_1d_int(&basin_size, number_of_po);
	for (int i = 0; i < number_of_po; i++) basin_size[i] = 0;
	double *entropy_progress;
	alloc_1d_double(&entropy_progress, analysis.number_of_rand_orbits_mc);
	for (int i = 0; i < analysis.number_of_rand_orbits_mc; i++) entropy_progress[i] = NAN;
	int orbits_counter = 0;
	time_t t;	srand((unsigned) time(&t));
	bool test_convergence = false;
	int converged_orbit_number;
	int *control;
	alloc_1d_int(&control, analysis.number_of_rand_orbits_mc);

	for (int i = 0; i < number_of_po; i++)
	{
		for (int j = 0; j < po[i].period; j++)
		{
			fprintf(out_ref, "%1.15e %1.15e %d %d\n",
				angle_mod(po[i].orbit[j][0]), 
				po[i].orbit[j][1],
				po[i].winding_number_numerator,
				po[i].winding_number_denominator);
		}
		fprintf(out_ref, "\n");
	}
	fclose(out_ref);

	#pragma omp parallel private(y,	converged_po_id, convergence_time) shared(basin_size, test_convergence) //num_threads(THREADS_NUM)
	{
	#pragma omp for schedule(dynamic)
		for (int i = 0; i < analysis.number_of_rand_orbits_mc; i++)
		{
			y[0] = rand_number_in_interval(analysis.grid_coordinate_min, analysis.grid_coordinate_max);
			y[1] = rand_number_in_interval(analysis.grid_velocity_min, analysis.grid_velocity_max);
			complete_orbital_part(y, system);

			evolve_multiple_basin_determined(y, number_of_po, 
				&converged_po_id, &convergence_time,
				po, system, analysis, rk);

			#pragma omp critical
			{
				// print progress
				print_prog((double)++orbits_counter/(double)analysis.number_of_rand_orbits_mc);

				if (converged_po_id != -1)
				{
					basin_size[converged_po_id]++;
					control[orbits_counter-1] = converged_po_id + 1;
				}
				else
				{
					control[orbits_counter-1] = -1;
				}

				entropy_progress[orbits_counter-1] = 
					basin_entropy(orbits_counter, number_of_po, basin_size, po, analysis); 
				
				// test if method converged
				if (test_convergence)
				{
					if (orbits_counter > analysis.convergence_window_mc)
					{
						double max_entropy_progress = entropy_progress[orbits_counter-1];
						double min_entropy_progress = entropy_progress[orbits_counter-1];
						for (int m = 2; m < analysis.convergence_window_mc + 1; m++)
						{
							if(entropy_progress[orbits_counter-m] > max_entropy_progress)
							{
								max_entropy_progress = entropy_progress[orbits_counter-m];
							}
							else if(entropy_progress[orbits_counter-m] < min_entropy_progress)
							{
								min_entropy_progress = entropy_progress[orbits_counter-m];
							}
						}
						if (fabs(max_entropy_progress-min_entropy_progress) < analysis.convergence_precision_mc)
						{
							printf("Method converged\n");
							test_convergence = false;
							converged_orbit_number = orbits_counter-1;
						}
					}
				}
			}
		}
	} // end pragma
	printf("\n");

	for (int i = 0; i < analysis.number_of_rand_orbits_mc; i++)
	{
		fprintf(out_control, "%d %d\n", i, control[i]);
	}
	fclose(out_control);

	fprintf(out_converged_time, "%1.3f %d %d\n", 
		e, converged_orbit_number, converged_orbit_number-analysis.convergence_window_mc);
	fclose(out_converged_time);

	// free memory
	dealloc_1d_int(&basin_size);
	dealloc_1d_int(&control);
	dealloc_1d_double(&entropy_progress);

	printf("Data written in output/basin_of_attraction/\n");

	return 0;
}

int multiple_basin_of_attraction_determined_monte_carlo_with_break	(int number_of_po,
											 			 			 perorb po[],
                         					 			 			 dynsys system,
                         					 			 			 anlsis analysis,
																	 rngkta rk)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/basin_of_attraction", &st) == -1) {
		mkdir("output/basin_of_attraction", 0700);
	}

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	// prepare and open exit files
	FILE	*out_control, *out_ref, *out_converged_number, *out_times;
	char	filename[300];

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_control_monte_carlo_with_break_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc);
	out_control = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_ref_monte_carlo_with_break_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc);
	out_ref = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_method_converged_number_monte_carlo_with_break_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc);
	out_converged_number = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_times_monte_carlo_with_break_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc);
	out_times = fopen(filename, "w");

	// declare variables
	double y[system.dim];
	int converged_po_id, convergence_time;
	int *basin_size;
	alloc_1d_int(&basin_size, number_of_po);
	for (int i = 0; i < number_of_po; i++) basin_size[i] = 0;
	double *entropy_progress;
	alloc_1d_double(&entropy_progress, analysis.number_of_rand_orbits_mc);
	for (int i = 0; i < analysis.number_of_rand_orbits_mc; i++) entropy_progress[i] = NAN;
	int orbits_counter = 0;
	time_t t;	srand((unsigned) time(&t));
	int converged_orbit_number;
	int *control;
	double *times;
	alloc_1d_int(&control, analysis.number_of_rand_orbits_mc);
	alloc_1d_double(&times, analysis.number_of_rand_orbits_mc);

	for (int i = 0; i < number_of_po; i++)
	{
		for (int j = 0; j < po[i].period; j++)
		{
			fprintf(out_ref, "%1.15e %1.15e %d %d\n",
				angle_mod(po[i].orbit[j][0]), 
				po[i].orbit[j][1],
				po[i].winding_number_numerator,
				po[i].winding_number_denominator);
		}
		fprintf(out_ref, "\n");
	}
	fclose(out_ref);

	const int N = analysis.number_of_rand_orbits_mc;
	bool go = true;
	unsigned int give = 0;

	#pragma omp parallel private(y,	converged_po_id, convergence_time) shared(basin_size, go) //num_threads(THREADS_NUM)
	{
		unsigned int loop_variable, stop;

		#pragma omp critical
		{
			loop_variable = give;
			give += N/omp_get_num_threads();
			stop = give;

			if(omp_get_thread_num() == omp_get_num_threads()-1)
				stop = N;
		} 

		while(loop_variable < stop && go)
		{
			y[0] = rand_number_in_interval(analysis.grid_coordinate_min, analysis.grid_coordinate_max);
			y[1] = rand_number_in_interval(analysis.grid_velocity_min, analysis.grid_velocity_max);
			complete_orbital_part(y, system);

			evolve_multiple_basin_determined(y, number_of_po, 
				&converged_po_id, &convergence_time,
				po, system, analysis, rk);

			#pragma omp critical
			{
				// print progress
				print_prog((double)++orbits_counter/(double)analysis.number_of_rand_orbits_mc);

				if (converged_po_id != -1)
				{
					basin_size[converged_po_id]++;
					control[orbits_counter-1] = converged_po_id + 1;
					times[orbits_counter-1] = (double)(convergence_time);
				}
				else
				{
					control[orbits_counter-1] = -1;
					times[orbits_counter-1] = NAN;
				}

				entropy_progress[orbits_counter-1] = 
					basin_entropy(orbits_counter, number_of_po, basin_size, po, analysis); 
				
				// test if method converged
				if (orbits_counter > analysis.convergence_window_mc)
				{
					double max_entropy_progress = entropy_progress[orbits_counter-1];
					double min_entropy_progress = entropy_progress[orbits_counter-1];
					for (int m = 2; m < analysis.convergence_window_mc + 1; m++)
					{
						if(entropy_progress[orbits_counter-m] > max_entropy_progress)
						{
							max_entropy_progress = entropy_progress[orbits_counter-m];
						}
						else if(entropy_progress[orbits_counter-m] < min_entropy_progress)
						{
							min_entropy_progress = entropy_progress[orbits_counter-m];
						}
					}
					if (fabs(max_entropy_progress-min_entropy_progress) < analysis.convergence_precision_mc)
					{
						printf("Method converged\n");
						go = false;
						converged_orbit_number = orbits_counter-1;
					}
				}
			}
			loop_variable++;
		}
	} // end pragma
	printf("\n");

	for (int i = 0; i < converged_orbit_number + 1; i++)
	{
		fprintf(out_control, "%d %d\n", i, control[i]);
		fprintf(out_times, "%d %1.5e\n", i, times[i]);
	}
	fclose(out_control);
	fclose(out_times);

	fprintf(out_converged_number, "%1.3f %d %d\n", 
		e, converged_orbit_number, converged_orbit_number-analysis.convergence_window_mc);
	fclose(out_converged_number);

	// free memory
	dealloc_1d_int(&basin_size);
	dealloc_1d_int(&control);
	dealloc_1d_double(&times);
	dealloc_1d_double(&entropy_progress);

	printf("Data written in output/basin_of_attraction/\n");

	return 0;

}

int evolve_multiple_basin_undetermined_winding	(double *ic,
									 			 bool *converged,
									 			 int *convergence_time,
									 			 atrtor *A,
									 			 dynsys system,
									 			 anlsis analysis,
												 rngkta rk)
{
	double *par = (double *)system.params;
	double T = par[7];
	
	// declare variables
	int 	orbit_counter;
	double 	y[system.dim], y_ref[system.dim];
	double	t = 0.0;

	double 	winding_number_window[analysis.convergence_window_wn];

	double *angle_progress;
	alloc_1d_double(&angle_progress, 1);

	// takes into consideration initial condition
	orbit_counter = 1;

	*converged = false;

	copy(y, ic, system.dim);
	for (int i = 0; i < analysis.number_of_cycles; i++)
	{
		// check if the rotation number converged
		int i_winding = i - (analysis.convergence_transient_wn + 1);
		if (i == analysis.convergence_transient_wn) // i_winding = -1
		{
			copy(y_ref, y, system.dim);
			angle_progress[0] = y[0];
		}
		else if (i > analysis.convergence_transient_wn)	// i_winding >= 0
		{
			double 	delta_theta = y[0] - y_ref[0];
			double 	wn = delta_theta / (double) (i_winding + 1);

			angle_progress = (double*) realloc(angle_progress, (i_winding + 2) * sizeof(double));
			angle_progress[i_winding + 1] = y[0];

			if(i_winding < analysis.convergence_window_wn)
			{
				winding_number_window[i_winding] = wn;
			}
			else if (i_winding >= analysis.convergence_window_wn)
			{
				double max_winding_number = wn;
				double min_winding_number = wn;
				for (int j = 0; j < analysis.convergence_window_wn - 1; j++)
				{
					winding_number_window[j] = winding_number_window[j+1];
					if(winding_number_window[j] > max_winding_number)
					{
						max_winding_number = winding_number_window[j];
					}
					if(winding_number_window[j] < min_winding_number)
					{
						min_winding_number = winding_number_window[j];
					}
				}
				winding_number_window[analysis.convergence_window_wn - 1] = wn;
				if (fabs(max_winding_number-min_winding_number) < analysis.convergence_precision_wn)
				{
					*converged = true;

					double pwn_numerator = 0.0, pwn_denominator = 0.0;
					for (int l = 1; l <= i_winding; l++)
					{
						double twp = (((double) l) / ((double) (i_winding + 1)));
						double factor = 1.0 / exp(1.0 / (twp * (1.0 - twp)));
						pwn_numerator += (angle_progress[l] - angle_progress[l-1]) * factor;
						pwn_denominator += factor;
					}
					double precise_wn = pwn_numerator / pwn_denominator;
					// printf("precise_wn = %1.10e\n", precise_wn);

					wn = precise_wn;

					A->winding_number = wn;
					A->basin_size = 1;
					A->theta = y[0];
					A->theta_dot = y[1];

					/* all this next part is wrong. I could not find a why to diferentiate 
					 * between resonances with the same winding number */

					double 	wn_mod = angle_mod_pos(wn);
					if(fabs(wn_mod) < 5e-2) wn_mod = 2.0*M_PI;

					double 	dist_from_int = fabs((2.0*M_PI/wn_mod) - round(2.0*M_PI/wn_mod));

					if(dist_from_int < 1e-3)
					{
						int		wn_period = (int) round(2.0*M_PI/wn_mod);
						double 	one_period_angular_diff = angle_progress[i_winding + 1] - angle_progress[i_winding + 1 - wn_period];
						int 	number_of_spins = (int) round(one_period_angular_diff / T);

						A->res_spin = number_of_spins;
						A->res_orbit = wn_period;
					}
					else
					{
						A->res_spin = 0;
						A->res_orbit = 0;
					}
					goto out;

					// A->res_spin = 1;	// dummy value
					// A->res_orbit = 1;	// dummy value

					// bool resonance_found = false;
					// double 	wn_mod = angle_mod_pos(wn);
					// if(fabs(wn_mod) < 5e-2) wn_mod = 2.0*M_PI;
					// double 	dist_from_int = fabs((2.0*M_PI/wn_mod) - round(2.0*M_PI/wn_mod));
					// if(dist_from_int < 5e-2)
					// {
					// 	double	y_new[system.dim];
					// 	double	t_new = 0.0;
					// 	copy(y_new, y, system.dim);
					// 	for (int number_of_orbits = analysis.orbit_period_min; number_of_orbits <= analysis.orbit_period_max; number_of_orbits++)
					// 	{
					// 		evolve_cycle(y_new, &t_new, system, analysis);
					// 		double angular_dif = angular_dist(y_new[0], y[0]);
					// 		// printf("%f %f %f %d\n", angular_dif, y[0], y_new[0], number_of_orbits);
					// 		if (angular_dif < 1e-1)
					// 		{
					// 			double 	one_period_angular_diff = y_new[0] - y[0];
					// 			int 	number_of_spins = (int) round(one_period_angular_diff / T);
					// 			A->res_spin = number_of_spins;
					// 			A->res_orbit = number_of_orbits;
					// 			resonance_found = true;
					// 			break;
					// 		}
					// 	}
					// }
					// else
					// {
					// 	A->res_spin = 0;
					// 	A->res_orbit = 0;
					// 	resonance_found = true;
					// }

					// if (resonance_found == true)
					// {
					// 	goto out;
					// }
					// else
					// {
					// 	printf("Ressonance not found inside given range\n");
					// 	// printf("x = %1.10e y = %1.10e\n", ic[0], ic[1]);
					// 	// exit(184);
					// 	A->res_spin = -2;
					// 	A->res_orbit = -2;
					// 	goto out;
					// }
				}
			}

		}

		evolve_cycle(y, &t, system, analysis, rk);

		// check if orbit diverges
		for (int j = 0; j < system.dim; j++)
		{
			if (fabs(y[j]) > analysis.evolve_box_size)
			{
				printf("Warning: box limit reached\n");
				printf("y[%d] = %1.10e\n", j, y[j]);
				goto out;
			}
		}

		orbit_counter++;

	}

	out:;

	if (*converged == true)
	{
		*convergence_time = orbit_counter - analysis.convergence_window_wn;
	}
	else
	{
		printf("Orbit that did not converge found\n");
		// printf("x = %1.10e y = %1.10e\n", ic[0], ic[1]);
		// exit(183);
		A->basin_size = 1;
		A->winding_number = NAN;
		A->res_spin = -1;
		A->res_orbit = -1;
		A->theta = NAN;
		A->theta_dot = NAN;
		*convergence_time = -1;
	}

	dealloc_1d_double(&angle_progress);

	return 0;
}

int multiple_basin_of_attraction_undetermined_monte_carlo_with_break(dynsys system,
                         					 			 			 anlsis analysis,
																	 rngkta rk)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/basin_of_attraction", &st) == -1) {
		mkdir("output/basin_of_attraction", 0700);
	}

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	// prepare and open exit files
	FILE	*out_control, *out_ref, *out_converged_number, *out_times;
	FILE	*out_size_full, *out_size, *out_entropy_prog, *out_entropy;
	char	filename[300];

	sprintf(filename, "output/basin_of_attraction/multiple_basin_undetermined_control_monte_carlo_with_break_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_rand_%d_window_mc_%d_precision_mc_%1.3f_transient_wn_%d_window_wn_%d_precision_wn_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc, analysis.convergence_transient_wn, analysis.convergence_window_wn, analysis.convergence_precision_wn);
	out_control = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/multiple_basin_undetermined_ref_monte_carlo_with_break_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_rand_%d_window_mc_%d_precision_mc_%1.3f_transient_wn_%d_window_wn_%d_precision_wn_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc, analysis.convergence_transient_wn, analysis.convergence_window_wn, analysis.convergence_precision_wn);
	out_ref = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/multiple_basin_undetermined_method_converged_number_monte_carlo_with_break_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_rand_%d_window_mc_%d_precision_mc_%1.3f_transient_wn_%d_window_wn_%d_precision_wn_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc, analysis.convergence_transient_wn, analysis.convergence_window_wn, analysis.convergence_precision_wn);
	out_converged_number = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/multiple_basin_undetermined_times_monte_carlo_with_break_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_rand_%d_window_mc_%d_precision_mc_%1.3f_transient_wn_%d_window_wn_%d_precision_wn_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc, analysis.convergence_transient_wn, analysis.convergence_window_wn, analysis.convergence_precision_wn);
	out_times = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/multiple_basin_undetermined_size_full_monte_carlo_with_break_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_rand_%d_window_mc_%d_precision_mc_%1.3f_transient_wn_%d_window_wn_%d_precision_wn_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc, analysis.convergence_transient_wn, analysis.convergence_window_wn, analysis.convergence_precision_wn);
	out_size_full = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/multiple_basin_undetermined_size_monte_carlo_with_break_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_rand_%d_window_mc_%d_precision_mc_%1.3f_transient_wn_%d_window_wn_%d_precision_wn_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc, analysis.convergence_transient_wn, analysis.convergence_window_wn, analysis.convergence_precision_wn);
	out_size = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/multiple_basin_undetermined_entropy_progress_monte_carlo_with_break_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_rand_%d_window_mc_%d_precision_mc_%1.3f_transient_wn_%d_window_wn_%d_precision_wn_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc, analysis.convergence_transient_wn, analysis.convergence_window_wn, analysis.convergence_precision_wn);
	out_entropy_prog = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/multiple_basin_undetermined_entropy_monte_carlo_with_break_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_rand_%d_window_mc_%d_precision_mc_%1.3f_transient_wn_%d_window_wn_%d_precision_wn_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc, analysis.convergence_transient_wn, analysis.convergence_window_wn, analysis.convergence_precision_wn);
	out_entropy = fopen(filename, "w");


	// declare variables
	double y[system.dim];
	int converged_po_id, convergence_time;
	double *entropy_progress;
	alloc_1d_double(&entropy_progress, analysis.number_of_rand_orbits_mc);
	for (int i = 0; i < analysis.number_of_rand_orbits_mc; i++) entropy_progress[i] = NAN;
	int orbits_counter = 0;
	time_t t;
	srand((unsigned) time(&t));
	int converged_orbit_number = analysis.number_of_rand_orbits_mc;
	int *control;
	double *times;
	alloc_1d_int(&control, analysis.number_of_rand_orbits_mc);
	alloc_1d_double(&times, analysis.number_of_rand_orbits_mc);

	atrtor 		A, *A_all;
	int			number_of_attractors = 0;

	bool		converged;

	int			index_for_not_converged = -1;

	volatile bool flag = false;

	#pragma omp parallel private(y,	convergence_time, converged, A) shared(flag, A_all, number_of_attractors, index_for_not_converged, orbits_counter) //num_threads(THREADS_NUM)
	{
	#pragma omp for schedule(dynamic)
		for (int i = 0; i < analysis.number_of_rand_orbits_mc; i++)
		{
			if(flag) continue;

			y[0] = rand_number_in_interval(analysis.grid_coordinate_min, analysis.grid_coordinate_max);
			y[1] = rand_number_in_interval(analysis.grid_velocity_min, analysis.grid_velocity_max);
			complete_orbital_part(y, system);

			double y_keep[system.dim];
			copy(y_keep, y, system.dim);

			evolve_multiple_basin_undetermined_winding(y,
				&converged, &convergence_time, &A, 
				system, analysis, rk);

			#pragma omp critical
			{
				if (flag == false)
				{
					// print progress
					print_prog((double)++orbits_counter/(double)analysis.number_of_rand_orbits_mc);
					// print_prog_float((double)++orbits_counter/(double)analysis.number_of_rand_orbits_mc);

					if (converged == true)
					{
						if (number_of_attractors == 0)
						{
							number_of_attractors++;
							A_all = (atrtor*) malloc(number_of_attractors * sizeof(atrtor));
							A_all[number_of_attractors-1] = A;
							control[orbits_counter-1] = number_of_attractors;
						}
						else
						{
							bool attractor_already_logged = false;
							for (int j = 0; j < number_of_attractors; j++)
							{
								if ((A.res_spin == A_all[j].res_spin) &&
									(A.res_orbit == A_all[j].res_orbit))
								{
									A_all[j].basin_size++;
									attractor_already_logged = true;
									control[orbits_counter-1] = j + 1;
									break;
								}
							}
							if(attractor_already_logged == false)
							{
								number_of_attractors++;
								A_all = (atrtor*) realloc(A_all, number_of_attractors * sizeof(atrtor));
								A_all[number_of_attractors-1] = A;
								control[orbits_counter-1] = number_of_attractors;
							}
						}
						times[orbits_counter-1] = (double)(convergence_time);
					}
					else
					{
						if (index_for_not_converged != -1)
						{
							A_all[index_for_not_converged].basin_size++;
						}
						else
						{
							if (number_of_attractors == 0)
							{
								number_of_attractors++;
								A_all = (atrtor*) malloc(number_of_attractors * sizeof(atrtor));
							}
							else
							{
								number_of_attractors++;
								A_all = (atrtor*) realloc(A_all, number_of_attractors * sizeof(atrtor));
							}
							A_all[number_of_attractors-1] = A;
							index_for_not_converged = number_of_attractors-1;
						}
						control[orbits_counter-1] = -1;
						times[orbits_counter-1] = NAN;
					}

					double entropy = 0.0;
					for (int j = 0; j < number_of_attractors; j++)
					{
						if(A_all[j].basin_size > 0)
						{
							double size_frac = 
								(double)A_all[j].basin_size / (double)orbits_counter;
							entropy += size_frac * log (1.0 / size_frac);
						}
					}
					if(number_of_attractors > 1) entropy /= log(number_of_attractors);

					entropy_progress[orbits_counter-1] = entropy; 
					
					// test if method converged
					if (orbits_counter > analysis.convergence_window_mc)
					{
						double max_entropy_progress = entropy_progress[orbits_counter-1];
						double min_entropy_progress = entropy_progress[orbits_counter-1];
						for (int m = 2; m < analysis.convergence_window_mc + 1; m++)
						{
							if(entropy_progress[orbits_counter-m] > max_entropy_progress)
							{
								max_entropy_progress = entropy_progress[orbits_counter-m];
							}
							else if(entropy_progress[orbits_counter-m] < min_entropy_progress)
							{
								min_entropy_progress = entropy_progress[orbits_counter-m];
							}
						}
						if (fabs(max_entropy_progress-min_entropy_progress) < analysis.convergence_precision_mc
							&& flag == false)
						{
							flag = true;
							printf(" Method converged");
							converged_orbit_number = orbits_counter-1;
						}
					}
				}
			}
		}
	} // end pragma
	printf("\n");

	for (int i = 0; i < converged_orbit_number + 1; i++)
	{
		fprintf(out_control, "%d %d\n", i, control[i]);
		fprintf(out_times, "%d %1.5e\n", i, times[i]);
	}
	fclose(out_control);
	fclose(out_times);

	fprintf(out_converged_number, "%1.3f %d %d\n", 
		e, converged_orbit_number, converged_orbit_number-analysis.convergence_window_mc);
	fclose(out_converged_number);

	for (int i = 0; i < number_of_attractors; i++)
	{
		if (A_all[i].res_orbit == 0)
		{
			fprintf(out_ref, "%1.15e %1.15e %d %d %1.5e %1.5e\n",
				angle_mod(A_all[i].theta), 
				A_all[i].theta_dot,
				A_all[i].res_spin,
				A_all[i].res_orbit,
				A_all[i].winding_number, 
				angle_mod(A_all[i].winding_number));
		}
		else
		{
			double y_local[system.dim], t_local = 0.0;
			y_local[0] = A_all[i].theta;
			y_local[1] = A_all[i].theta_dot;
			complete_orbital_part(y_local, system);

			for (int j = 0; j < A_all[i].res_orbit; j++)
			{
				fprintf(out_ref, "%1.15e %1.15e %d %d %1.5e %1.5e\n",
					angle_mod(y_local[0]), 
					y_local[1],
					A_all[i].res_spin,
					A_all[i].res_orbit,
					A_all[i].winding_number, 
					angle_mod(A_all[i].winding_number));

				evolve_cycle(y_local, &t_local, system, analysis, rk);
			}
		}
		fprintf(out_ref, "\n");
	}
	fclose(out_ref);

	for (int i = 0; i < number_of_attractors; i++)
	{
		fprintf(out_size_full, "w = %f s = %d o = %d size = %f\n", 
			A_all[i].winding_number,
			A_all[i].res_spin,
			A_all[i].res_orbit,
			A_all[i].basin_size/((double)orbits_counter));
	}
	fclose(out_size_full);

	bool empty_res;
	for (int orbit_period = analysis.orbit_period_min; orbit_period <= analysis.orbit_period_max; orbit_period++)
	{
		for (int spin_period = analysis.spin_period_min; spin_period <= analysis.spin_period_max; spin_period++)
		{
			empty_res = true;
			for (int k = 0; k < number_of_attractors; k++)
			{
				if (A_all[k].res_spin == spin_period &&
					A_all[k].res_orbit == orbit_period)
				{
					fprintf(out_size, "%1.3f %d %d %1.10f %d\n", 
						e, spin_period, orbit_period, A_all[k].basin_size/((double)orbits_counter),
						cantor_pairing_function(spin_period, orbit_period));
					empty_res = false;
				}
			}
			if (empty_res == true)
			{
				fprintf(out_size, "%1.3f %d %d %1.10f %d\n", 
					e, spin_period, orbit_period, 0.0,
					cantor_pairing_function(spin_period, orbit_period));
			}
		}
	}
	empty_res = true;
	for (int k = 0; k < number_of_attractors; k++)
	{
		if (A_all[k].res_spin == 0 &&
			A_all[k].res_orbit == 0)
		{
			fprintf(out_size, "%1.3f %d %d %1.10f %d\n", 
				e, A_all[k].res_spin, A_all[k].res_orbit, A_all[k].basin_size/((double)orbits_counter),
				cantor_pairing_function(A_all[k].res_spin, A_all[k].res_orbit));
			empty_res = false;
		}
	}
	if (empty_res == true)
	{
		fprintf(out_size, "%1.3f %d %d %1.10f %d\n", 
			e, 0, 0, 0.0,
			cantor_pairing_function(0, 0));
	}
	fclose(out_size);

	for (int i = 0; i < converged_orbit_number; i++)
	{
		fprintf(out_entropy_prog, "%1.10e\n", entropy_progress[i]);
	}
	fclose(out_entropy_prog);

	fprintf(out_entropy, "%1.3f %1.10e\n", e, entropy_progress[converged_orbit_number-1]);
	fclose(out_entropy);	

	// free memory
	free(A_all);
	dealloc_1d_int(&control);
	dealloc_1d_double(&times);
	dealloc_1d_double(&entropy_progress);

	printf("Data written in output/basin_of_attraction/\n");

	return 0;

}

int comparison_entropy_grid_vs_monte_carlo	(int number_of_po,
											 perorb po[],
                         					 dynsys system,
                         					 anlsis analysis)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/basin_of_attraction", &st) == -1) {
		mkdir("output/basin_of_attraction", 0700);
	}

	double	*par = (double *)system.params;
	double	gamma = par[0];
	double	e = par[1];
	double	K = par[6];

	int		**control_matrix;
	int 	basin_size[number_of_po];

	fill_control_matrix(&control_matrix, system, analysis);

	FILE	*out_grid, *out_mc, *in_entropy_grid, *in_entropy_mc;
	char	filename[300];

	sprintf(filename, "output/basin_of_attraction/comparison_entropy_grid_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	out_grid = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/comparison_entropy_monte_carlo_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc);
	out_mc = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_entropy_progress_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc);
	in_entropy_grid = fopen(filename, "r");
	if (in_entropy_grid == NULL)
	{
		printf("File not found\n");
		exit(10);
	}
	int counter_grid;
	double *entropy_progress_grid, entropy_grid;
	alloc_1d_double(&entropy_progress_grid, analysis.number_of_rand_orbits_mc);
	for (int i = 0; i < analysis.number_of_rand_orbits_mc; i++)
	{
		entropy_progress_grid[i] = NAN;
	}
	while(fscanf(in_entropy_grid, "%d %lf", &counter_grid, &entropy_grid) != EOF)
	{
		entropy_progress_grid[counter_grid] = entropy_grid;
	}
	fclose(in_entropy_grid);

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_entropy_progress_monte_carlo_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.dat", 
	gamma, e, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc);
	in_entropy_mc = fopen(filename, "r");
	if (in_entropy_mc == NULL)
	{
		printf("File not found\n");
		exit(10);
	}
	int counter_mc;
	double *entropy_progress_mc, entropy_mc;
	alloc_1d_double(&entropy_progress_mc, analysis.number_of_rand_orbits_mc);
	while(fscanf(in_entropy_mc, "%d %lf", &counter_mc, &entropy_mc) != EOF)
	{
		entropy_progress_mc[counter_mc] = entropy_mc;
	}
	fclose(in_entropy_mc);

	int last_value_grid;
	for (int i = analysis.number_of_rand_orbits_mc - 1; i > 0; i--)
	{
		if (entropy_progress_grid[i] == entropy_progress_grid[i])
		{
			last_value_grid = i;
			break;
		}
	}

	for (int i = 0; i < analysis.number_of_rand_orbits_mc; i++)
	{
		if (entropy_progress_grid[i] == entropy_progress_grid[i])
		{
			fprintf(out_grid, "%d %1.10f %1.10f\n", i+1, entropy_progress_grid[i], 
				fabs(entropy_progress_grid[i]-entropy_progress_grid[last_value_grid]));
			fprintf(out_mc, "%d %1.10f %1.10f\n", i+1, entropy_progress_mc[i], 
				fabs(entropy_progress_mc[i]-entropy_progress_mc[last_value_grid]));
		}
	}
	fclose(out_grid);
	fclose(out_mc);

	dealloc_1d_double(&entropy_progress_grid);
	dealloc_1d_double(&entropy_progress_mc);
	dealloc_2d_int(&control_matrix, analysis.grid_resolution);

	return 0;
}

int basin_entropy_vs_box_size	(int number_of_po,
								 perorb po[],
                         		 dynsys system,
                         		 anlsis analysis)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/basin_of_attraction", &st) == -1) {
		mkdir("output/basin_of_attraction", 0700);
	}

	double	*par = (double *)system.params;
	double	gamma = par[0];
	double	e = par[1];
	double	K = par[6];

	double 	**basin_matrix;
	int		**control_matrix;
	fill_basin_matrix(&basin_matrix, system, analysis);
	fill_control_matrix(&control_matrix, system, analysis);

	FILE	*out_entropy;
	char	filename[300];

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_entropy_vs_box_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	out_entropy = fopen(filename, "w");

	for (int step = 1; step < (analysis.grid_resolution / analysis.sqrt_orbits_on_box); step++)
	{
		if ((analysis.grid_resolution % (analysis.sqrt_orbits_on_box * step)) == 0)
		{
			double gibbs_entropy = 0.0;
			for (int n = 0; n <= (analysis.grid_resolution / analysis.sqrt_orbits_on_box) - step; n = n + step)
			{
				for (int m = 0; m <= (analysis.grid_resolution / analysis.sqrt_orbits_on_box) - step; m = m + step)
				{
					int non_zero_prob[number_of_po];
					double prob[number_of_po];
					for (int p = 0; p < number_of_po; p++)
					{
						non_zero_prob[p] = -1;
						prob[p] = 0.0;
					}
					for (int i = n * analysis.sqrt_orbits_on_box; i < (n + step) * analysis.sqrt_orbits_on_box; i = i + step)
					{
						for (int j = m * analysis.sqrt_orbits_on_box; j < (m + step) * analysis.sqrt_orbits_on_box; j = j + step)
						{
							if (control_matrix[i][j] != -1)
							{
								prob[control_matrix[i][j] - 1] += 1.0;
								non_zero_prob[control_matrix[i][j] - 1] = 1;
							}
						}
					}
					double gibbs_entropy_on_box = 0.0;
					int number_of_at = 0;
					for (int p = 0; p < number_of_po; p++)
					{
						if (non_zero_prob[p] != -1)
						{
							prob[p] /= ((double)(analysis.sqrt_orbits_on_box * analysis.sqrt_orbits_on_box));
							gibbs_entropy_on_box += prob[p] * log (1.0 / prob[p]);
							number_of_at++;
						}
					}
					gibbs_entropy_on_box /= (double)number_of_at; // don't know if this helps or not
					gibbs_entropy += gibbs_entropy_on_box;
				}
			}
			double sqrt_number_of_boxes = (double) (analysis.grid_resolution / (analysis.sqrt_orbits_on_box * step));
			double gibbs_entropy_per_box = gibbs_entropy / (sqrt_number_of_boxes * sqrt_number_of_boxes);
			double gibbs_entropy_per_box_2 = gibbs_entropy / (double)(analysis.grid_resolution * analysis.grid_resolution);
			fprintf(out_entropy, "%1.5e %1.5e %1.5e %1.5e %1.5e\n",  
				log(1.0 / sqrt_number_of_boxes),
				log(1.0 / ((double)(analysis.sqrt_orbits_on_box * step))), 
				log(gibbs_entropy),
				log(gibbs_entropy_per_box),
				log(gibbs_entropy_per_box_2));
		}
	}
	fclose(out_entropy);

	dealloc_2d_double(&basin_matrix, analysis.grid_resolution);
	dealloc_2d_int(&control_matrix, analysis.grid_resolution);

	return 0;
}

int basin_size_from_data(int number_of_po,
						 perorb po[],
                         dynsys system,
                         anlsis analysis)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/basin_of_attraction", &st) == -1) {
		mkdir("output/basin_of_attraction", 0700);
	}

	double	*par = (double *)system.params;
	double	gamma = par[0];
	double	e = par[1];
	double	K = par[6];

	int		**control_matrix;
	fill_control_matrix(&control_matrix, system, analysis);

	FILE	*out_size;
	char	filename[300];

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_size_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	out_size = fopen(filename, "w");

	int 	basin_size[number_of_po];

	for (int i = 0; i < number_of_po; i++) basin_size[i] = 0;

	for (int i = 0; i < analysis.grid_resolution; i++)
	{
		for (int j = 0; j < analysis.grid_resolution; j++)
		{
			if (control_matrix[i][j] > 0)
			{
				basin_size[control_matrix[i][j] - 1]++;
			}
		}
	}

	double basin_size_fraction_sum = 0.0;
	for (int orbit_period = analysis.orbit_period_min; orbit_period <= analysis.orbit_period_max; orbit_period++)
	{
		for (int spin_period = analysis.spin_period_min; spin_period <= analysis.spin_period_max; spin_period++)
		{
			int basin_size_sum = 0;
			double basin_size_fraction = 0.0;
			for (int k = 0; k < number_of_po; k++)
			{
				if (po[k].winding_number_numerator == spin_period &&
					po[k].winding_number_denominator == orbit_period)
				{
					basin_size_sum += basin_size[k];
				}
			}
			if(basin_size_sum > 0)
			{
				basin_size_fraction = (double)basin_size_sum / 
					((double)(analysis.grid_resolution * analysis.grid_resolution));
			}

			fprintf(out_size, "%1.3f %d %d %1.10f %d\n", 
				e, spin_period, orbit_period, basin_size_fraction,
				cantor_pairing_function(spin_period, orbit_period));	
			
			basin_size_fraction_sum += basin_size_fraction;
		}
	}
	fprintf(out_size, "%1.3f %d %d %1.10f nan\n", e, 0, 0, 1.0 - basin_size_fraction_sum);
	fclose(out_size);

	dealloc_2d_int(&control_matrix, analysis.grid_resolution);

	return 0;
}

int basin_size_from_data_monte_carlo(int number_of_po,
						 			 perorb po[],
                         			 dynsys system,
                         			 anlsis analysis)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/basin_of_attraction", &st) == -1) {
		mkdir("output/basin_of_attraction", 0700);
	}

	double	*par = (double *)system.params;
	double	gamma = par[0];
	double	e = par[1];
	double	K = par[6];

	int		*control;
	fill_control_monte_carlo(&control, system, analysis);

	FILE	*out_size;
	char	filename[300];

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_size_monte_carlo_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc);
	out_size = fopen(filename, "w");

	int 	basin_size[number_of_po];

	for (int i = 0; i < number_of_po; i++) basin_size[i] = 0;

	for (int i = 0; i < analysis.number_of_rand_orbits_mc; i++)
	{
		if (control[i] > 0)
		{
			basin_size[control[i] - 1]++;
		}
	}

	double basin_size_fraction_sum = 0.0;
	for (int orbit_period = analysis.orbit_period_min; orbit_period <= analysis.orbit_period_max; orbit_period++)
	{
		for (int spin_period = analysis.spin_period_min; spin_period <= analysis.spin_period_max; spin_period++)
		{
			int basin_size_sum = 0;
			double basin_size_fraction = 0.0;
			for (int k = 0; k < number_of_po; k++)
			{
				if (po[k].winding_number_numerator == spin_period &&
					po[k].winding_number_denominator == orbit_period)
				{
					basin_size_sum += basin_size[k];
				}
			}
			if(basin_size_sum > 0)
			{
				basin_size_fraction = (double)basin_size_sum / 
					((double)(analysis.grid_resolution * analysis.grid_resolution));
			}

			fprintf(out_size, "%1.3f %d %d %1.10f %d\n", 
				e, spin_period, orbit_period, basin_size_fraction,
				cantor_pairing_function(spin_period, orbit_period));	
			
			basin_size_fraction_sum += basin_size_fraction;
		}
	}
	fprintf(out_size, "%1.3f %d %d %1.10f nan\n", e, 0, 0, 1.0 - basin_size_fraction_sum);
	fclose(out_size);

	dealloc_1d_int(&control);

	return 0;
}

int basin_size_from_data_monte_carlo_with_break(int number_of_po,
						 			 			perorb po[],
                         			 			dynsys system,
                         			 			anlsis analysis)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/basin_of_attraction", &st) == -1) {
		mkdir("output/basin_of_attraction", 0700);
	}

	double	*par = (double *)system.params;
	double	gamma = par[0];
	double	e = par[1];
	double	K = par[6];

	int		control_size, *control;
	fill_control_monte_carlo_with_break(&control_size, &control, system, analysis);

	FILE	*out_size;
	char	filename[300];

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_size_monte_carlo_with_break_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc);
	out_size = fopen(filename, "w");

	int 	basin_size[number_of_po];

	for (int i = 0; i < number_of_po; i++) basin_size[i] = 0;

	for (int i = 0; i < control_size; i++)
	{
		if (control[i] > 0)
		{
			basin_size[control[i] - 1]++;
		}
	}

	double basin_size_fraction_sum = 0.0;
	for (int orbit_period = analysis.orbit_period_min; orbit_period <= analysis.orbit_period_max; orbit_period++)
	{
		for (int spin_period = analysis.spin_period_min; spin_period <= analysis.spin_period_max; spin_period++)
		{
			int basin_size_sum = 0;
			double basin_size_fraction = 0.0;
			for (int k = 0; k < number_of_po; k++)
			{
				if (po[k].winding_number_numerator == spin_period &&
					po[k].winding_number_denominator == orbit_period)
				{
					basin_size_sum += basin_size[k];
				}
			}
			if(basin_size_sum > 0)
			{
				basin_size_fraction = (double)basin_size_sum / 
					((double)(control_size));
			}

			fprintf(out_size, "%1.3f %d %d %1.10f %d\n", 
				e, spin_period, orbit_period, basin_size_fraction,
				cantor_pairing_function(spin_period, orbit_period));	
			
			basin_size_fraction_sum += basin_size_fraction;
		}
	}
	fprintf(out_size, "%1.3f %d %d %1.10f nan\n", e, 0, 0, 1.0 - basin_size_fraction_sum);
	fclose(out_size);

	dealloc_1d_int(&control);

	return 0;
}

int basin_entropy_from_data (dynsys system,
                         	 anlsis analysis)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/basin_of_attraction", &st) == -1) {
		mkdir("output/basin_of_attraction", 0700);
	}

	double	*par = (double *)system.params;
	double	gamma = par[0];
	double	e = par[1];
	double	K = par[6];

	FILE	*out_entropy, *in_basin_size_grid;
	char	filename[300];

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_entropy_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	out_entropy = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_size_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat", 
			gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	in_basin_size_grid = fopen(filename, "r");
	if (in_basin_size_grid == NULL)
	{
		printf("File not found\n");
		exit(10);
	}
	
	int 	fool_spin, fool_orbit, number_attractors;
	double	ec, size_frac, entropy, fool_cantor;

	printf("Getting size of basins from printed file\n");

	number_attractors = 0;
	while(fscanf(in_basin_size_grid, "%lf %d %d %lf %lf\n", &ec, &fool_spin, &fool_orbit, &size_frac, &fool_cantor) != EOF)
	{
		if (size_frac > 0.0)
		{
			entropy += size_frac * log (1.0 / size_frac);
			number_attractors++;
		}
	}
	fclose(in_basin_size_grid);

	if (number_attractors > 1) entropy /= log(number_attractors);

	fprintf(out_entropy, "%1.3f %1.10f\n", 
		e, entropy);
	fclose(out_entropy);

	return 0;
}

int basin_entropy_from_data_monte_carlo (dynsys system,
                         	 			 anlsis analysis)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/basin_of_attraction", &st) == -1) {
		mkdir("output/basin_of_attraction", 0700);
	}

	double	*par = (double *)system.params;
	double	gamma = par[0];
	double	e = par[1];
	double	K = par[6];

	FILE	*out_entropy, *in_basin_size_mc;
	char	filename[300];

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_entropy_monte_carlo_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc);
	out_entropy = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_size_monte_carlo_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc);
	in_basin_size_mc = fopen(filename, "r");
	if (in_basin_size_mc == NULL)
	{
		printf("File not found\n");
		exit(10);
	}
	
	int 	fool_spin, fool_orbit, number_attractors;
	double	ec, size_frac, entropy, fool_cantor;

	printf("Getting size of basins monte carlo from printed file\n");

	number_attractors = 0;
	while(fscanf(in_basin_size_mc, "%lf %d %d %lf %lf\n", &ec, &fool_spin, &fool_orbit, &size_frac, &fool_cantor) != EOF)
	{
		if (size_frac > 0.0)
		{
			entropy += size_frac * log (1.0 / size_frac);
			number_attractors++;
		}
	}
	fclose(in_basin_size_mc);

	if (number_attractors > 1) entropy /= log(number_attractors);

	fprintf(out_entropy, "%1.3f %1.10f\n", 
		e, entropy);
	fclose(out_entropy);

	return 0;
}

int basin_entropy_from_data_monte_carlo_with_break	(dynsys system,
                         	 			 			 anlsis analysis)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/basin_of_attraction", &st) == -1) {
		mkdir("output/basin_of_attraction", 0700);
	}

	double	*par = (double *)system.params;
	double	gamma = par[0];
	double	e = par[1];
	double	K = par[6];

	FILE	*out_entropy, *in_basin_size_mc;
	char	filename[300];

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_entropy_monte_carlo_with_break_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc);
	out_entropy = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_size_monte_carlo_with_break_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc);
	in_basin_size_mc = fopen(filename, "r");
	if (in_basin_size_mc == NULL)
	{
		printf("File not found\n");
		exit(10);
	}
	
	int 	fool_spin, fool_orbit, number_attractors;
	double	ec, size_frac, entropy, fool_cantor;

	printf("Getting size of basins monte carlo from printed file\n");

	number_attractors = 0;
	while(fscanf(in_basin_size_mc, "%lf %d %d %lf %lf\n", &ec, &fool_spin, &fool_orbit, &size_frac, &fool_cantor) != EOF)
	{
		if (size_frac > 0.0)
		{
			entropy += size_frac * log (1.0 / size_frac);
			number_attractors++;
		}
	}
	fclose(in_basin_size_mc);

	if (number_attractors > 1) entropy /= log(number_attractors);

	fprintf(out_entropy, "%1.3f %1.10f\n", 
		e, entropy);
	fclose(out_entropy);

	return 0;
}

int basin_entropy_progress_from_data(int number_of_po,
						 			 perorb po[],
                         			 dynsys system,
                         			 anlsis analysis)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/basin_of_attraction", &st) == -1) {
		mkdir("output/basin_of_attraction", 0700);
	}

	double	*par = (double *)system.params;
	double	gamma = par[0];
	double	e = par[1];
	double	K = par[6];

	int		**control_matrix;
	int 	basin_size[number_of_po];

	fill_control_matrix(&control_matrix, system, analysis);

	FILE	*out_progress;
	char	filename[300];

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_entropy_progress_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc);
	out_progress = fopen(filename, "w");

	for (int n = analysis.grid_resolution; n > 0; n--)
	{
		if (analysis.grid_resolution % n == 0)
		{
			for (int p = 0; p < number_of_po; p++)
			{
				basin_size[p] = 0;
			}
			int local_counter = 0;
			for (int i = 0; i < analysis.grid_resolution; i = i + n)
			{
				for (int j = 0; j < analysis.grid_resolution; j = j + n)
				{
					if (control_matrix[i][j] > 0)
					{
						basin_size[control_matrix[i][j] - 1]++;
					}
					local_counter++;
				}
			}
			double entropy = basin_entropy(local_counter, number_of_po, basin_size, po, analysis);
			fprintf(out_progress, "%d %1.10e\n", local_counter - 1, entropy);
		}
	}
	fclose(out_progress);

	dealloc_2d_int(&control_matrix, analysis.grid_resolution);

	return 0;
}

int basin_entropy_progress_from_data_monte_carlo(int number_of_po,
						 			 			 perorb po[],
                         			 			 dynsys system,
                         			 			 anlsis analysis)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/basin_of_attraction", &st) == -1) {
		mkdir("output/basin_of_attraction", 0700);
	}

	double	*par = (double *)system.params;
	double	gamma = par[0];
	double	e = par[1];
	double	K = par[6];

	double	*entropy_progress;
	alloc_1d_double(&entropy_progress, analysis.number_of_rand_orbits_mc);

	int		*control;
	fill_control_monte_carlo(&control, system, analysis);

	FILE	*out_progress;
	char	filename[300];

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_entropy_progress_monte_carlo_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc);
	out_progress = fopen(filename, "w");

	int 	basin_size[number_of_po];

	for (int i = 0; i < number_of_po; i++) basin_size[i] = 0;

	for (int i = 0; i < analysis.number_of_rand_orbits_mc; i++)
	{
		if (control[i] > 0)
		{
			basin_size[control[i] - 1]++;
		}

		entropy_progress[i] = basin_entropy(i+1, number_of_po, basin_size, po, analysis); 
	}

	for (int i = 0; i < analysis.number_of_rand_orbits_mc; i++)
	{
		fprintf(out_progress, "%d %1.10e\n", i, entropy_progress[i]);
	}
	fclose(out_progress);

	dealloc_1d_int(&control);
	dealloc_1d_double(&entropy_progress);

	return 0;
}

int basin_entropy_progress_from_data_monte_carlo_with_break(int number_of_po,
						 			 			 			perorb po[],
                         			 			 			dynsys system,
                         			 			 			anlsis analysis)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/basin_of_attraction", &st) == -1) {
		mkdir("output/basin_of_attraction", 0700);
	}

	double	*par = (double *)system.params;
	double	gamma = par[0];
	double	e = par[1];
	double	K = par[6];

	int		control_size, *control;
	fill_control_monte_carlo_with_break(&control_size, &control, system, analysis);

	double	*entropy_progress;
	alloc_1d_double(&entropy_progress, control_size);

	FILE	*out_progress;
	char	filename[300];

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_entropy_progress_monte_carlo_with_break_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc);
	out_progress = fopen(filename, "w");

	int 	basin_size[number_of_po];

	for (int i = 0; i < number_of_po; i++) basin_size[i] = 0;

	for (int i = 0; i < control_size; i++)
	{
		if (control[i] > 0)
		{
			basin_size[control[i] - 1]++;
		}

		entropy_progress[i] = basin_entropy(i+1, number_of_po, basin_size, po, analysis); 
	}

	for (int i = 0; i < control_size; i++)
	{
		fprintf(out_progress, "%d %1.10e\n", i, entropy_progress[i]);
	}
	fclose(out_progress);

	dealloc_1d_int(&control);
	dealloc_1d_double(&entropy_progress);

	return 0;
}

int evolve_multiple_basin_undetermined  (double *ic,
									     bool *converged,
                                         int *attractor_period,
									     int *convergence_time,
									     dynsys system,
									     anlsis analysis,
										 rngkta rk)
{
	// declare variables
	bool close_encounter;
	int close_time_counter;
	int time_to_look_for_po = analysis.number_of_cycles / 10;
	double y[system.dim], y_copy[system.dim];
	double t = 0.0;
	int max_po_period;
	int internal_attractor_period;
	int internal_convergence_time;
	double close_encounter_tolerance;

	// max po period 
 	max_po_period = 2;

	// starts proximity counter
	close_time_counter = 0;
	
	close_encounter = false;
	*converged = false;
	close_encounter_tolerance = 1e-8;

	double window[max_po_period][system.dim];

	copy(y, ic, system.dim);
	for (int i = 0; i < analysis.number_of_cycles; i++)
	{
		// fill the window
		if (i < max_po_period)
		{
			copy(window[i], y, system.dim);
		}
		else
		{
			if (close_encounter == false)
			{
				if (i > time_to_look_for_po)
				{
					for (int j = max_po_period - 1; j >= 0; j--)
					{
						if (dist_from_ref(window[j], y) 
							< close_encounter_tolerance)
						{
							internal_attractor_period = max_po_period - j;
							close_encounter = true;
							internal_convergence_time = i + 1;
							copy(y_copy, y, system.dim);
						}
					}
				}
			}
			else
			{
				if (i % internal_attractor_period == 0)
				{
					if (dist_from_ref(y, y_copy) < close_encounter_tolerance)
					{
						close_time_counter += internal_attractor_period;
					}
					else
					{
						close_time_counter = 0;
						close_encounter == false;
					}
				}
				if (close_time_counter > analysis.evolve_basin_time_tol)
				{
					*converged = true;
					goto out;
				}
			}
			
			// shift the window forward
			for (int j = 0; j < max_po_period - 1; j++)
			{
				copy(window[j], window[j + 1], system.dim);
			}
			copy(window[max_po_period - 1], y, system.dim);
		}

		evolve_cycle(y, &t, system, analysis, rk);
	
		// check if orbit diverges
		for (int j = 0; j < system.dim; j++)
		{
			if (fabs(y[j]) > analysis.evolve_box_size)
			{
				printf("Warning: box limit reached\n");
				printf("y[%d] = %1.10e\n", j, y[j]);
				goto out;
			}
		}
	}

	out:;
	
	if (*converged == true)
	{
		*convergence_time = internal_convergence_time;
		*attractor_period = internal_attractor_period;
	}
	else
	{
		*convergence_time = -1;
		*attractor_period = -1;
	}

	return 0;
}

int multiple_basin_of_attraction_undetermined	(dynsys system,
                         					 	 anlsis analysis,
												 rngkta rk)
{
	printf("Warning: multiple_basin_of_attraction_undetermined is very outdated\n");
	printf("Compare it to multiple_basin_of_attraction_determined before usage\n");
	exit(11);

	// little warning
	if (strcmp(system.name, "two_body") == 0)
	{
		printf("Warning: cant draw basins\n");
		printf("for two-body system\n");
		exit(2);
	}

	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/basin_of_attraction", &st) == -1) {
		mkdir("output/basin_of_attraction", 0700);
	}

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	// prepare and open exit files
	FILE	*out_boa, *out_ref, *out_size;
	char	filename[300];

	sprintf(filename, "output/basin_of_attraction/multiple_basin_undetermined_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d.dat", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles);
	out_boa = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/multiple_basin_undetermined_ref_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d.dat", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles);
	out_ref = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/multiple_basin_undetermined_size_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d.dat", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles);
	out_size = fopen(filename, "w");

	// declare variables
	double y[system.dim];
	double rot_ini[2];
	double orb[4], orb_ini[4];
	init_orbital(orb_ini, system);
	int grid[2];
	double basin[2];
	bool converged;
	bool po_already_found;
	int converged_po_id;
	int convergence_time;
	int attractor_period;
	int number_of_po;
	double eps;
	int *basin_size;
	perorb *multiple_po;
	int **basin_matrix, **control_matrix;
	double **time_matrix;

	alloc_2d_int(&basin_matrix, 
		analysis.grid_resolution, analysis.grid_resolution);
	alloc_2d_int(&control_matrix, 
		analysis.grid_resolution, analysis.grid_resolution);
	alloc_2d_double(&time_matrix, 
		analysis.grid_resolution, analysis.grid_resolution);
	for (int i = 0; i < analysis.grid_resolution; i++)
	{
		for (int j = 0; j < analysis.grid_resolution; j++)
		{
			basin_matrix[i][j] = 0;
			control_matrix[i][j] = 0;
			time_matrix[i][j] = NAN;
		}
	}

	// tolerance for comparing p.o.s
	eps = 1e-5;

	number_of_po = 0;
	// loop over coordinate values
	for (int i = 0; i < analysis.grid_resolution; i++)
	{
		// print progress on coordinate
		printf("Calculating set %d of %d\n", 
					i + 1, analysis.grid_resolution);

		#pragma omp parallel private(y, grid, converged, po_already_found, attractor_period, \
				converged_po_id, convergence_time, orb, rot_ini) shared(basin_matrix, \
				control_matrix, time_matrix, number_of_po, basin_size, multiple_po) //num_threads(THREADS_NUM)
		{

		#pragma omp for schedule(dynamic)
			// loop over velocity values
			for (int j = 0; j < analysis.grid_resolution; j++)
			{
				// print progress on velocity
				// printf("Calculating subset %d of %d\n", 
				// 			j + 1, analysis.grid_resolution);

				if (control_matrix[i][j] == 0)
				{
					grid[0] = i;
					grid[1] = j;

					grid_to_double(grid, rot_ini, analysis);

					copy(y, rot_ini, 2);

					if (system.dim == 6)
					{
						for (int k = 0; k < 4; k++)
						{
							y[k+2] = orb_ini[k];
						}
					}

					evolve_multiple_basin_undetermined(y, &converged,
						&attractor_period, &convergence_time,
						system, analysis, rk);

					// #pragma omp critical
					// {
						if (converged == true)
						{
							// po_already_found = false;
							// for (int k = 0; k < number_of_po; k++)
							// {
							// 	if (multiple_po[k].period == attractor_period)
							// 	{
							// 		for (int l = 0; l < multiple_po[k].period; l++)
							// 		{
							// 			if (dist_from_ref(y, multiple_po[k].orbit[l]) < eps)
							// 			{
							// 				po_already_found = true;
							// 				converged_po_id = k;
							// 			}
							// 		}
							// 	}
							// }
							// if(po_already_found == false)
							{
								perorb po_try;
								copy(po_try.seed, y, 2);
								po_try.period = attractor_period;
								alloc_2d_double(&po_try.orbit, po_try.period, system.dim);
								if (periodic_orbit(&po_try, system, analysis, rk) == 0)
								{
									for (int k = 0; k < number_of_po; k++)
									{
										if (multiple_po[k].period == attractor_period)
										{
											for (int l = 0; l < multiple_po[k].period; l++)
											{
												for (int m = 0; m < po_try.period; m++)
												{
													if (dist_from_ref(po_try.orbit[m], multiple_po[k].orbit[l]) < eps)
													{
														// po_already_found = true;
														converged_po_id = k;
														// printf("po was already found\n");
														goto po_found;
													}
												}
											}
										}
									}
									number_of_po++;
									converged_po_id = number_of_po - 1;
									if (number_of_po == 1)
									{
										multiple_po = (perorb*) malloc(number_of_po * sizeof(perorb));
										basin_size = (int*) malloc(number_of_po * sizeof(int));
									}
									else
									{
										multiple_po = realloc(multiple_po, number_of_po * sizeof(perorb));
										basin_size = realloc(basin_size, number_of_po * sizeof(int));
									}
									copy(multiple_po[converged_po_id].seed, po_try.seed, 2);
									copy(multiple_po[converged_po_id].initial_condition, po_try.seed, 2);
									multiple_po[converged_po_id].period = po_try.period;
									multiple_po[converged_po_id].winding_number_numerator = po_try.winding_number_numerator;
									multiple_po[converged_po_id].winding_number_denominator = po_try.winding_number_denominator;
									alloc_2d_double(&multiple_po[converged_po_id].orbit, multiple_po[converged_po_id].period, system.dim);
									for(int k = 0; k < po_try.period; k++)
									{
										copy(multiple_po[converged_po_id].orbit[k], po_try.orbit[k], system.dim);
									}
									dealloc_2d_double(&po_try.orbit, po_try.period);
								}
								else
								{
									dealloc_2d_double(&po_try.orbit, po_try.period);
									goto po_not_found;
								}
								
							}
							po_found:;
							// printf("w = %d / %d\n", multiple_po[converged_po_id].winding_number_numerator,multiple_po[converged_po_id].winding_number_denominator);
							basin_matrix[i][j] = 
								cantor_pairing_function(multiple_po[converged_po_id].winding_number_numerator,multiple_po[converged_po_id].winding_number_denominator);
							control_matrix[i][j] = converged_po_id + 1;
							time_matrix[i][j] = (double)(convergence_time);
							basin_size[converged_po_id]++;
						}
						else
						{
							po_not_found:;
							control_matrix[i][j] = -1;
						}
					// } // end pragma critical
				}
			}
		} // end pragma

		// new line on terminal
		// printf("\n");
		printf("number_of_po = %d\n", number_of_po);
	}

	if (number_of_po > 0)
	{

		// print p.o.s to file
		for (int i = 0; i < number_of_po; i++)
		{
			for (int j = 0; j < multiple_po[i].period; j++)
			{
				fprintf(out_ref, "%1.15e %1.15e %d %d\n",
					angle_mod(multiple_po[i].orbit[j][0]), multiple_po[i].orbit[j][1],
					multiple_po[i].winding_number_numerator,
					multiple_po[i].winding_number_denominator);
			}
			fprintf(out_ref, "\n");
		}

		int res_spin, res_orbit;

		for (int i = 0; i < analysis.grid_resolution; i++)
		{
			for (int j = 0; j < analysis.grid_resolution; j++)
			{
				grid[0] = i;
				grid[1] = j;
				grid_to_double(grid, basin, analysis);

				if (control_matrix[i][j] != -1)
				{
					res_spin = multiple_po[control_matrix[i][j]-1].winding_number_numerator;
					res_orbit = multiple_po[control_matrix[i][j]-1].winding_number_denominator;
				}
				else
				{
					res_spin = 0;
					res_orbit = 0;
				}

				fprintf(out_boa, "%1.5f %1.5f %d %1.5e %d %d\n", 
					basin[0], basin[1], 
					basin_matrix[grid[0]][grid[1]],
					time_matrix[grid[0]][grid[1]],
					res_spin,
					res_orbit);

			}
			fprintf(out_boa, "\n");
		}

		int winding_number_numerator_max = 0;
		int winding_number_denominator_max = 0;

		for (int i = 0; i < number_of_po; i++)
		{
			if (multiple_po[i].winding_number_numerator > winding_number_numerator_max)
			{
				winding_number_numerator_max = multiple_po[i].winding_number_numerator;
			}
			if (multiple_po[i].winding_number_denominator > winding_number_denominator_max)
			{
				winding_number_denominator_max = multiple_po[i].winding_number_denominator;
			}
		}

		double basin_size_fraction;
		double basin_size_fraction_sum;

		basin_size_fraction_sum = 0.0;
		for (int i = 1; i <= winding_number_numerator_max; i++)
		{
			for (int j = 1; j <= winding_number_denominator_max; j++)
			{
				basin_size_fraction = 0.0;
				for (int k = 0; k < number_of_po; k++)
				{
					if (multiple_po[k].winding_number_numerator == i &&
						multiple_po[k].winding_number_denominator == j)
					{
						basin_size_fraction += ((double)basin_size[k]) / 
							(((double)analysis.grid_resolution) * ((double)analysis.grid_resolution));
					}
				}

				fprintf(out_size, "%d %d %1.4f\n", 
					i, j, basin_size_fraction);	
				
				basin_size_fraction_sum += basin_size_fraction;
			}
		}

		fprintf(out_size, "n m %1.4f\n", basin_size_fraction_sum);

		dealloc_1d_int(&basin_size);
		for (int i = 0; i < number_of_po; i++)
		{
			dealloc_2d_double(&multiple_po[i].orbit, multiple_po[i].period);
		}
		free(multiple_po);
	}

	// close exit files
	fclose(out_boa);
	fclose(out_ref);
	fclose(out_size);

	dealloc_2d_int(&basin_matrix, 
					analysis.grid_resolution);

	dealloc_2d_int(&control_matrix, 
					analysis.grid_resolution);

	dealloc_2d_double(&time_matrix, 
					analysis.grid_resolution);

	printf("Data written in output/basin_of_attraction/\n");

	return 0;
}

double basin_entropy(int number_of_orbits,
					 int number_of_po,
                     int *basin_size,
                     perorb po[],
					 anlsis analysis)
{
	int number_of_attractors = 0;
	int total_basin_size_sum = 0;
	double entropy = 0.0;
	for (int k = analysis.orbit_period_min; k <= analysis.orbit_period_max; k++)
	{
		for (int l = analysis.spin_period_min; l <= analysis.spin_period_max; l++)
		{
			int basin_size_sum = 0;
			for (int j = 0; j < number_of_po; j++)
			{
				if (po[j].winding_number_numerator == l &&
					po[j].winding_number_denominator == k)
				{
					basin_size_sum += basin_size[j];
				}
			}
			if(basin_size_sum > 0)
			{
				double size_frac = 
					(double)basin_size_sum / (double)number_of_orbits;
				entropy += size_frac * log (1.0 / size_frac);
				number_of_attractors++;
			}
			total_basin_size_sum += basin_size_sum;
		}
	}
	if (total_basin_size_sum != number_of_orbits)
	{
		double size_frac = 
			(double)(number_of_orbits - total_basin_size_sum) / (double)number_of_orbits;
		entropy += size_frac * log (1.0 / size_frac);
		number_of_attractors++;
	}
	if(number_of_attractors > 1)
	{
		entropy /= log(number_of_attractors);
	}
	return entropy;
}

int trace_ellipse()
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/tests", &st) == -1) {
		mkdir("output/tests", 0700);
	}

	FILE *out = fopen("output/tests/ellipse.dat","w");

	double a = 2.0;
	double e = 0.5;
	double b = a * sqrt(1.0 - e*e);
	double delta = 1e-3;
	double k = 0.0;
	double h = -1.0 * e * a;
	double x_max = h + a;
	double x_min = h - a;

	for (double x = x_max; x > x_min; x -= delta)
	{
		double y = k + sqrt(b*b * (1.0 - (x-h)*(x-h)/(a*a)));
		fprintf(out, "%1.5e %1.5e\n", x, y);
	}
	for (double x = x_min; x < x_max; x += delta)
	{
		double y = k - sqrt(b*b * (1.0 - (x-h)*(x-h)/(a*a)));
		fprintf(out, "%1.5e %1.5e\n", x, y);
	}
	double x = x_max;
	double y = k + sqrt(b*b * (1.0 - (x-h)*(x-h)/(a*a)));
	fprintf(out, "%1.5e %1.5e\n", x, y);

	fclose(out);

	return 0;
}

int draw_orbit_map(dynsys system)
{
	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	printf("Drawing orbit of system %s with gamma = %1.3f, e = %1.3f and K = %1.5f\n", 
		system.name, gamma, e, K);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/orbit\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/orbit/fig_orbit_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f.png\"\n", 
		gamma, e, system.name, K);
	fprintf(gnuplotPipe, "set xlabel \"{/Symbol q}\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "set xrange[-3.1415:3.1415]\n");
	fprintf(gnuplotPipe, "set yrange [0.0:3.0]\n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, 
		"set title \"gamma = %1.3f    e = %1.3f    K = %1.5f\"\n", 
		gamma, e, K);
	fprintf(gnuplotPipe, "plot 'orbit_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f.dat' w p pt 7 ps 1.5 palette notitle, 'orbit_ic_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f.dat' w p pt 7 ps 1.5 notitle",
		gamma, e, system.name, K, gamma, e, system.name, K);
	fclose(gnuplotPipe);

	printf("Done!\n");
	printf("Data written in output/orbit/\n");

	return 0;
}

int draw_phase_space(dynsys system,
					 anlsis analysis)
{
	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];

	printf("Drawing phase space with gamma = %1.6f and e = %1.3f\n", 
		gamma, e);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/phase_space\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/phase_space/fig_phase_space_gamma_%1.6f_e_%1.3f.png\"\n", gamma, e);
	fprintf(gnuplotPipe, "set xlabel \"{/Symbol q}\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "set xrange[-3.1415:3.1415]\n");
	fprintf(gnuplotPipe, "set yrange [%f:%f]\n", analysis.velocity_min, analysis.velocity_max);
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, 
		"set key title \"~{/Symbol g}{0.5-} = %1.6f e = %1.3f\" box opaque top right width 2\n", 
		gamma, e);
	fprintf(gnuplotPipe, "plot 'phase_space_gamma_%1.6f_e_%1.3f.dat' w d notitle",
		gamma, e);
	fclose(gnuplotPipe);

	printf("Done!\n");

	return 0;
}

int draw_phase_space_clean(dynsys system)
{

	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/clean_figures", &st) == -1) {
		mkdir("output/clean_figures", 0700);
	}

	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];

	printf("Drawing phase space with gamma = %1.3f and e = %1.3f\n", 
		gamma, e);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set loadpath \"output/phase_space\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/clean_figures/fig_phase_space_gamma_%1.6f_e_%1.3f.png\"\n", gamma, e);
	fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,50\"\n");
	fprintf(gnuplotPipe, "set size square \n");
	fprintf(gnuplotPipe, "set lmargin at screen 0.05\n");
	fprintf(gnuplotPipe, "set bmargin at screen 0.05\n");
	fprintf(gnuplotPipe, "set rmargin at screen 0.95\n");
	fprintf(gnuplotPipe, "set tmargin at screen 0.95\n");
	fprintf(gnuplotPipe, "set border lw 2 \n");
	fprintf(gnuplotPipe, "set xrange[-3.1415:3.1415]\n");
	fprintf(gnuplotPipe, "set yrange [0.0:3.0]\n");
	fprintf(gnuplotPipe, "set xtics format \" \"\n");
	fprintf(gnuplotPipe, "set ytics format \" \"\n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, "unset title\n");
	// fprintf(gnuplotPipe, "plot 'phase_space_gamma_%1.6f_e_%1.3f.dat' w d notitle", gamma, e);
	fprintf(gnuplotPipe, "plot 'phase_space_gamma_%1.6f_e_%1.3f.dat' w p pt 7 ps 0.2 notitle", gamma, e);
	fclose(gnuplotPipe);

	printf("Done!\n");

	return 0;
}

int draw_phase_space_latex(dynsys system)
{

	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/tests", &st) == -1) {
		mkdir("output/tests", 0700);
	}

	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];

	printf("Drawing phase space with gamma = %1.3f and e = %1.3f\n", 
		gamma, e);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set loadpath \"output/phase_space\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/tests/fig_phase_space_gamma_%1.6f_e_%1.3f.png\"\n", gamma, e);
	fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,50\"\n");
	fprintf(gnuplotPipe, "set size square \n");
	fprintf(gnuplotPipe, "set border lw 2 \n");
	fprintf(gnuplotPipe, "set xrange[-3.1415:3.1415]\n");
	fprintf(gnuplotPipe, "set yrange [0.0:3.0]\n");
	fprintf(gnuplotPipe, "set key font \"fonts/cmr10.ttf,35\" \n");
	fprintf(gnuplotPipe, "set xlabel \"{/Symbol q}\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, "unset title\n");
	// fprintf(gnuplotPipe, "plot 'phase_space_gamma_%1.6f_e_%1.3f.dat' w d lc rgb \"black\" notitle",
	// 	gamma, e);
	// fprintf(gnuplotPipe, "plot 'phase_space_gamma_%1.6f_e_%1.3f.dat' w d notitle",
	// 	gamma, e);
	fprintf(gnuplotPipe, "plot 'phase_space_gamma_%1.6f_e_%1.3f.dat' w p pt 7 ps 0.2 notitle", gamma, e);
	fclose(gnuplotPipe);

	printf("Done!\n");

	return 0;
}

int draw_orbit_on_phase_space(dynsys system)
{
	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	printf("Drawing orbit on phase space of system %s with gamma = %1.3f, e = %1.3f and K = %1.5f\n", 
		system.name, gamma, e, K);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/orbit/fig_orbit_on_phase_space_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f.png\"\n", 
		gamma, e, system.name, K);
	fprintf(gnuplotPipe, "set xlabel \"{/Symbol q}\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "set xrange[-3.1415:3.1415]\n");
	fprintf(gnuplotPipe, "set yrange [0.0:3.0]\n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, 
		"set title \"gamma = %1.3f    e = %1.3f    K = %1.5f\"\n", 
		gamma, e, K);
	fprintf(gnuplotPipe, "plot 'phase_space/phase_space_gamma_%1.6f_e_%1.3f.dat' w d lc rgb \"gray40\" notitle ,'orbit/orbit_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f.dat' w p pt 7 ps 1.5 palette notitle, 'orbit/orbit_ic_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f.dat' w p pt 7 ps 1.5 notitle",
		gamma, e, gamma, e, system.name, K, gamma, e, system.name, K);
	fclose(gnuplotPipe);

	printf("Done!\n");
	printf("Data written in output/orbit/\n");

	return 0;
}

int draw_orbit_on_phase_space_latex(dynsys system)
{
	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	printf("Drawing orbit on phase space of system %s with gamma = %1.3f, e = %1.3f and K = %1.5f\n", 
		system.name, gamma, e, K);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 2200,2000 font \"fonts/cmr10.ttf,50\"\n");
		fprintf(gnuplotPipe, "set key font \"fonts/cmr10.ttf,35\" \n");
	fprintf(gnuplotPipe, "set loadpath \"output\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/orbit/fig_orbit_on_phase_space_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_latex.png\"\n", 
		gamma, e, system.name, K);
	fprintf(gnuplotPipe, "set size square \n");
	fprintf(gnuplotPipe, "set border lw 2 \n");
	fprintf(gnuplotPipe, "set xlabel \"{/Symbol q}\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "set xrange[-3.1415:3.1415]\n");
	fprintf(gnuplotPipe, "set yrange [0.0:3.0]\n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, "unset title\n");
	fprintf(gnuplotPipe, "plot 'phase_space/phase_space_gamma_%1.6f_e_%1.3f.dat' w d notitle ,'orbit/orbit_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f.dat' w p pt 7 ps 3 palette notitle", 
		gamma, e, gamma, e, system.name, K);
	fclose(gnuplotPipe);

	printf("Done!\n");
	printf("Data written in output/orbit/\n");

	return 0;
}

int draw_time_series(dynsys system)
{
	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double e = par[1];
	double K = par[6];

	printf("Drawing time series with e = %1.3f and K = %1.5f\n", e, K);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/time_series\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/time_series/fig_time_series_e_%1.3f_K_%1.5f.png\"\n", e, K);
	fprintf(gnuplotPipe, "set xlabel \"n\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	// fprintf(gnuplotPipe, "set xrange [60:200] \n");
	// fprintf(gnuplotPipe, "set yrange [0.5:2] \n");
	fprintf(gnuplotPipe, "set log y\n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, 
		"set key title \"e = %1.3f K = %1.5f\" box opaque top right width 2\n", e, K);
	fprintf(gnuplotPipe, "plot 'time_series_e_%1.3f_K_%1.5f.dat' u 1:2 w l lw 2 notitle", e, K);
	fclose(gnuplotPipe);

	printf("Done!\n");

	return 0;
}

int draw_time_series_latex(dynsys system)
{
	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double e = par[1];
	double K = par[6];

	printf("Drawing time series with e = %1.3f and K = %1.5f\n", e, K);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 1000,1000 font \"fonts/cmr10.ttf,25\"\n");
	fprintf(gnuplotPipe, "set key font \"fonts/cmr10.ttf,20\" \n");
	fprintf(gnuplotPipe, "set loadpath \"output/time_series\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/time_series/fig_time_series_e_%1.3f_K_%1.5f_latex.png\"\n", e, K);
	fprintf(gnuplotPipe, "set size square \n");
	fprintf(gnuplotPipe, "set border lw 2 \n");
	fprintf(gnuplotPipe, "set xlabel \"n\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "unset title\n");
	// fprintf(gnuplotPipe, "set xrange [60:200] \n");
	// fprintf(gnuplotPipe, "set yrange [0.5:2] \n");
	fprintf(gnuplotPipe, "set log y\n");
	fprintf(gnuplotPipe, "unset key\n");
	// fprintf(gnuplotPipe, 
	// 	"set key title \"e = %1.3f K = %1.5f\" box opaque top right width 2\n", e, K);
	fprintf(gnuplotPipe, "plot 'time_series_e_%1.3f_K_%1.5f.dat' u 1:2 w l lw 2 notitle", e, K);
	fclose(gnuplotPipe);

	printf("Done!\n");

	return 0;
}

int draw_time_series_union_e(dynsys system)
{
	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double K = par[6];

	int color;
	double e;

	printf("Drawing union of time series\n");

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/time_series\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/time_series/fig_time_series_union_K_%1.5f.png\"\n", K);
	fprintf(gnuplotPipe, "set xlabel \"n\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "set title \"K = %1.5f\"\n", K);
	fprintf(gnuplotPipe, "set log y\n");
	fprintf(gnuplotPipe, "set linetype cycle 20\n");
	color = 1;
	e = 0.00;
	fprintf(gnuplotPipe, "plot 'time_series_e_%1.3f_K_%1.5f.dat' u 1:2 w l lw 2 lc %d title \"e = %1.3f\"", e, K, color, e);
	for (e = 0.02; e < 0.205; e += 0.02)
	{
		color++;
		fprintf(gnuplotPipe, ", 'time_series_e_%1.3f_K_%1.5f.dat' u 1:2 w l lw 2 lc %d title \"e = %1.3f\"", e, K, color, e);
	}
	fclose(gnuplotPipe);

	// gnuplotPipe = popen("gnuplot -persistent", "w");
	// fprintf(gnuplotPipe, "reset\n");
	// fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	// fprintf(gnuplotPipe, "set loadpath \"output/time_series\"\n");
	// fprintf(gnuplotPipe, 
	// 	"set output \"output/time_series/fig_time_series_union_K_%1.5f_zoom.png\"\n", K);
	// fprintf(gnuplotPipe, "set xlabel \"n\"\n");
	// fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	// fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	// fprintf(gnuplotPipe, "set title \"K = %1.5f\"\n", K);
	// fprintf(gnuplotPipe, "set log y\n");
	// fprintf(gnuplotPipe, "set xrange[0:120]\n");
	// fprintf(gnuplotPipe, "set yrange[0.8:1000]\n");
	// fprintf(gnuplotPipe, "set linetype cycle 20\n");
	// color = 1;
	// e = 0.00;
	// fprintf(gnuplotPipe, "plot 'time_series_e_%1.3f_K_%1.5f.dat' u 1:2 w l lw 2 lc %d title \"e = %1.3f\"", e, K, color, e);
	// for (e = 0.02; e < 0.205; e += 0.02)
	// {
	// 	color++;
	// 	fprintf(gnuplotPipe, ", 'time_series_e_%1.3f_K_%1.5f.dat' u 1:2 w l lw 2 lc %d title \"e = %1.3f\"", e, K, color, e);
	// }
	// fclose(gnuplotPipe);

	printf("Done!\n");

	return 0;
}

int draw_time_series_union_e_latex(dynsys system)
{
	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double K = par[6];

	int color;
	double e;

	printf("Drawing union of time series\n");

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 1000,1000 font \"fonts/cmr10.ttf,25\"\n");
	fprintf(gnuplotPipe, "set key font \"fonts/cmr10.ttf,20\" \n");
	fprintf(gnuplotPipe, "set loadpath \"output/time_series\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/time_series/fig_time_series_union_K_%1.5f_latex.png\"\n", K);
	fprintf(gnuplotPipe, "set size square \n");
	fprintf(gnuplotPipe, "set border lw 2 \n");
	fprintf(gnuplotPipe, "set xlabel \"n\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "unset title\n");
	fprintf(gnuplotPipe, "set key box lw 2 opaque\n");
	fprintf(gnuplotPipe, "set xrange[0:200]\n");
	// fprintf(gnuplotPipe, "set title \"K = %1.5f\"\n", K);
	fprintf(gnuplotPipe, "set log y\n");
	fprintf(gnuplotPipe, "set linetype cycle 20\n");
	color = 1;
	e = 0.00;
	fprintf(gnuplotPipe, "plot 'time_series_e_%1.3f_K_%1.5f.dat' u 1:2 w l lw 2 lc %d title \"e = %1.3f\"", e, K, color, e);
	for (e = 0.02; e < 0.205; e += 0.02)
	{
		color++;
		fprintf(gnuplotPipe, ", 'time_series_e_%1.3f_K_%1.5f.dat' u 1:2 w l lw 2 lc %d title \"e = %1.3f\"", e, K, color, e);
	}
	fclose(gnuplotPipe);

	printf("Done!\n");

	return 0;
}

int draw_time_series_union_e_eps(dynsys system)
{
	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double K = par[6];

	int color;
	double e;

	printf("Drawing union of time series\n");

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal postscript eps size 9.20,8.00 enhanced color font \"fonts/cmr10.ttf,50\" linewidth 2\n");
	fprintf(gnuplotPipe, "set key font \"fonts/cmr10.ttf,40\" \n");
	fprintf(gnuplotPipe, "set loadpath \"output/time_series\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/time_series/fig_time_series_union_K_%1.5f_eps.eps\"\n", K);
	fprintf(gnuplotPipe, "set size square \n");
	fprintf(gnuplotPipe, "set xlabel \"n\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "unset title\n");
	fprintf(gnuplotPipe, "set key box opaque\n");
	fprintf(gnuplotPipe, "set xrange[0:200]\n");
	// fprintf(gnuplotPipe, "set title \"K = %1.5f\"\n", K);
	fprintf(gnuplotPipe, "set log y\n");
	fprintf(gnuplotPipe, "set linetype cycle 20\n");
	color = 1;
	e = 0.00;
	fprintf(gnuplotPipe, "plot 'time_series_e_%1.3f_K_%1.5f.dat' u 1:2 w l lw 2 lc %d title \"e = %1.3f\"", e, K, color, e);
	for (e = 0.02; e < 0.205; e += 0.02)
	{
		color++;
		fprintf(gnuplotPipe, ", 'time_series_e_%1.3f_K_%1.5f.dat' u 1:2 w l lw 2 lc %d title \"e = %1.3f\"", e, K, color, e);
	}
	fclose(gnuplotPipe);

	printf("Done!\n");

	return 0;
}

int draw_time_series_union_K(dynsys system)
{
	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double e = par[1];

	int color;
	double K;

	printf("Drawing union of time series\n");

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/time_series\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/time_series/fig_time_series_union_e_%1.3f.png\"\n", e);
	fprintf(gnuplotPipe, "set xlabel \"n\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "set title \"e = %1.3f\"\n", e);
	fprintf(gnuplotPipe, "set log y\n");
	fprintf(gnuplotPipe, "set linetype cycle 20\n");
	color = 1;
	K = 0.01;
	fprintf(gnuplotPipe, "plot 'time_series_e_%1.3f_K_%1.5f.dat' u 1:2 w l lw 2 lc %d title \"K = %1.3f\"", e, K, color, K);
	for (K = 0.009; K > 0.00095; K -= 0.001)
	{
		color++;
		fprintf(gnuplotPipe, ", 'time_series_e_%1.3f_K_%1.5f.dat' u 1:2 w l lw 2 lc %d title \"K = %1.3f\"", e, K, color, K);
	}
	fclose(gnuplotPipe);

	printf("Done!\n");

	return 0;
}

int draw_time_series_union_K_latex(dynsys system)
{
	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double e = par[1];

	int color;
	double K;

	printf("Drawing union of time series\n");

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 1000,1000 font \"fonts/cmr10.ttf,25\"\n");
	fprintf(gnuplotPipe, "set key font \"fonts/cmr10.ttf,20\" \n");
	fprintf(gnuplotPipe, "set loadpath \"output/time_series\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/time_series/fig_time_series_union_e_%1.3f_latex.png\"\n", e);
	fprintf(gnuplotPipe, "set size square \n");
	fprintf(gnuplotPipe, "set border lw 2 \n");
	fprintf(gnuplotPipe, "set xlabel \"n\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "unset title\n");
	fprintf(gnuplotPipe, "set key box lw 2 opaque\n");
	// fprintf(gnuplotPipe, "set title \"e = %1.3f\"\n", e);
	fprintf(gnuplotPipe, "set log y\n");
	fprintf(gnuplotPipe, "set linetype cycle 20\n");
	color = 1;
	K = 0.01;
	fprintf(gnuplotPipe, "plot 'time_series_e_%1.3f_K_%1.5f.dat' u 1:2 w l lw 2 lc %d title \"K = %1.3f\"", e, K, color, K);
	for (K = 0.009; K > 0.00095; K -= 0.001)
	{
		color++;
		fprintf(gnuplotPipe, ", 'time_series_e_%1.3f_K_%1.5f.dat' u 1:2 w l lw 2 lc %d title \"K = %1.3f\"", e, K, color, K);
	}
	fclose(gnuplotPipe);

	printf("Done!\n");

	return 0;
}

int draw_multiple_time_series(dynsys system)
{
	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double e = par[1];
	double K = par[6];

	printf("Drawing multiple time series with e = %1.3f and K = %1.5f\n", e, K);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/time_series\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/time_series/fig_multiple_time_series_e_%1.3f_K_%1.5f.png\"\n", e, K);
	fprintf(gnuplotPipe, "set xlabel \"n\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "set log y\n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, 
		"set key title \"e = %1.3f K = %1.5f\" box opaque top right width 2\n", e, K);
	fprintf(gnuplotPipe, "plot 'multiple_time_series_e_%1.3f_K_%1.5f.dat' u 1:2 w l lw 2 notitle", e, K);
	fclose(gnuplotPipe);

	printf("Done!\n");

	return 0;
}

int draw_multiple_time_series_delta_theta_dot	(dynsys system,
												 anlsis analysis)
{
	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double e = par[1];
	double K = par[6];

	printf("Drawing multiple time series for e = %1.3f K = %1.5f and system = %s\n", e, K, system.name);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/time_series\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/time_series/fig_multiple_time_series_delta_theta_dot_e_%1.3f_K_%1.5f_nic_%d_delta_%1.2f_system_%s.png\"\n", e, K, analysis.number_of_time_series, analysis.time_series_delta, system.name);
	fprintf(gnuplotPipe, "set xlabel \"n\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "set yrange [0.01:1000] \n");
	// fprintf(gnuplotPipe, "set xrange [0.0:50] \n");
	fprintf(gnuplotPipe, "set log y\n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, 
		"set key title \"e = %1.3f  K = %1.5f  nic = %d  delta = %1.2f  system = %s\" box opaque top right\n", e, K, analysis.number_of_time_series, analysis.time_series_delta, system.name);
	fprintf(gnuplotPipe, "plot 'multiple_time_series_delta_theta_dot_e_%1.3f_K_%1.5f_nic_%d_delta_%1.2f_system_%s.dat' u 1:2:-2 w l lc var lw 2 notitle", e, K, analysis.number_of_time_series, analysis.time_series_delta, system.name);
	fclose(gnuplotPipe);

	printf("Done!\n");

	return 0;
}

int draw_multiple_time_series_delta_theta_dot_latex	(dynsys system,
												 	 anlsis analysis)
{
	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double e = par[1];
	double K = par[6];

	printf("Drawing multiple time series with e = %1.3f and K = %1.5f\n", e, K);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 1000,1000 font \"fonts/cmr10.ttf,25\"\n");
	fprintf(gnuplotPipe, "set key font \"fonts/cmr10.ttf,20\" \n");
	fprintf(gnuplotPipe, "set loadpath \"output/time_series\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/time_series/fig_multiple_time_series_delta_theta_dot_e_%1.3f_K_%1.5f_delta_%1.2f_latex.png\"\n", e, K, analysis.time_series_delta);
	fprintf(gnuplotPipe, "set size square \n");
	fprintf(gnuplotPipe, "set border lw 2 \n");
	fprintf(gnuplotPipe, "set xlabel \"n\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "unset title\n");
	fprintf(gnuplotPipe, "set yrange [0.01:1000] \n");
	// fprintf(gnuplotPipe, "set xrange [0.0:50] \n");
	fprintf(gnuplotPipe, "set log y\n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, "plot 'multiple_time_series_delta_theta_dot_e_%1.3f_K_%1.5f_delta_%1.2f.dat' u 1:2:-2 w l lc var lw 2 notitle", e, K, analysis.time_series_delta);
	fclose(gnuplotPipe);

	printf("Done!\n");

	return 0;
}


int draw_multiple_time_series_delta_theta   (dynsys system,
                                             anlsis analysis)
{
	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double e = par[1];
	double K = par[6];

	printf("Drawing multiple time series with e = %1.3f and K = %1.5f\n", e, K);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/time_series\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/time_series/fig_multiple_time_series_delta_theta_e_%1.3f_K_%1.5f_delta_%1.2f.png\"\n", e, K, analysis.time_series_delta);
	fprintf(gnuplotPipe, "set xlabel \"n\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "set log y\n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, 
		"set key title \"e = %1.3f K = %1.5f delta = %1.2f\" box opaque top right width 2\n", e, K, analysis.time_series_delta);
	fprintf(gnuplotPipe, "plot 'multiple_time_series_delta_theta_e_%1.3f_K_%1.5f_delta_%1.2f.dat' u 1:2:-2 w l lc var lw 2 notitle", e, K, analysis.time_series_delta);
	fclose(gnuplotPipe);

	printf("Done!\n");

	return 0;
}

int draw_periodic_orbit_on_phase_space  (perorb po,
                                         dynsys system)
{
	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	printf("Drawing periodic orbit of period %d on phase space of system %s with gamma = %1.3f, e = %1.3f and K = %1.5f\n", 
		po.period, system.name, gamma, e, K);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/periodic_orbit/fig_periodic_orbit_on_phase_space_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_period_%d_ic_%1.3f_%1.3f.png\"\n", 
		gamma, e, system.name, K, po.period, po.initial_condition[0], po.initial_condition[1]);
	fprintf(gnuplotPipe, "set xlabel \"{/Symbol q}\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "set xrange[-3.1416:3.1416]\n");
	fprintf(gnuplotPipe, "set yrange [0.0:3.0]\n");
	fprintf(gnuplotPipe, "unset key\n");
	if(system.name == "rigid")
	{
		fprintf(gnuplotPipe, 
			"set title \"Periodic orbit for system %s and gamma = %1.3f e = %1.3f. Resonance = %d / %d\"\n", 
			system.name, gamma, e, po.winding_number_numerator, po.winding_number_denominator);
	}
	else
	{
		fprintf(gnuplotPipe, 
			"set title \"Periodic orbit for system %s and gamma = %1.3f e = %1.3f K = %1.5f. Resonance = %d / %d\"\n", 
			system.name, gamma, e, K, po.winding_number_numerator, po.winding_number_denominator);
	}
	fprintf(gnuplotPipe, "plot 'phase_space/phase_space_gamma_%1.6f_e_%1.3f.dat' w d lc rgb \"gray40\" notitle ,'periodic_orbit/periodic_orbit_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_period_%d_ic_%1.3f_%1.3f.dat' w p pt 7 ps 1.5 lc rgb \"black\" notitle",
		gamma, e, gamma, e, system.name, K, po.period, po.initial_condition[0], po.initial_condition[1]);
	fclose(gnuplotPipe);

	printf("Done!\n");
	printf("Data written in output/periodic_orbit/\n");

	return 0;
}

int draw_periodic_orbit_on_phase_space_clean(perorb po,
                                         	 dynsys system)
{

	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/clean_figures", &st) == -1) {
		mkdir("output/clean_figures", 0700);
	}

	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	printf("Drawing periodic orbit of period %d on phase space of system %s with gamma = %1.3f, e = %1.3f and K = %1.5f\n", 
		po.period, system.name, gamma, e, K);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set loadpath \"output\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/clean_figures/fig_periodic_orbit_on_phase_space_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_period_%d_ic_%1.3f_%1.3f.png\"\n", 
		gamma, e, system.name, K, po.period, po.initial_condition[0], po.initial_condition[1]);
	fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,50\"\n");
	fprintf(gnuplotPipe, "set size square \n");
	fprintf(gnuplotPipe, "set lmargin at screen 0.05\n");
	fprintf(gnuplotPipe, "set bmargin at screen 0.05\n");
	fprintf(gnuplotPipe, "set rmargin at screen 0.95\n");
	fprintf(gnuplotPipe, "set tmargin at screen 0.95\n");
	fprintf(gnuplotPipe, "set border lw 2 \n");
	fprintf(gnuplotPipe, "set xrange[-3.1416:3.1416]\n");
	fprintf(gnuplotPipe, "set yrange [0.0:3.0]\n");
	fprintf(gnuplotPipe, "set xtics format \" \"\n");
	fprintf(gnuplotPipe, "set ytics format \" \"\n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, "unset title\n");
	// fprintf(gnuplotPipe, "plot 'phase_space/phase_space_gamma_%1.6f_e_%1.3f.dat' w d notitle ,'periodic_orbit/periodic_orbit_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_period_%d_ic_%1.3f_%1.3f.dat' w p pt 7 ps 4 lc rgb \"black\" notitle",
	// 	gamma, e, gamma, e, system.name, K, po.period, po.initial_condition[0], po.initial_condition[1]);
	fprintf(gnuplotPipe, "plot 'phase_space/phase_space_gamma_%1.6f_e_%1.3f.dat' w p pt 7 ps 0.2 notitle ,'periodic_orbit/periodic_orbit_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_period_%d_ic_%1.3f_%1.3f.dat' w p pt 7 ps 5 lc rgb \"black\" notitle",
		gamma, e, gamma, e, system.name, K, po.period, po.initial_condition[0], po.initial_condition[1]);
	fclose(gnuplotPipe);

	printf("Done!\n");
	printf("Data written in output/periodic_orbit/\n");

	return 0;
}

int draw_basin_of_attraction(perorb po,
                             dynsys system,
                             anlsis analysis)
{
	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	printf("Drawing basin of attraction of system %s with gamma = %1.3f, e = %1.3f and K = %1.5f and reference theta = %1.3f theta_dot = %1.3f with period = %d\n", 
		system.name, gamma, e, K, po.initial_condition[0], po.initial_condition[1], po.period);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/basin_of_attraction/fig_basin_of_attraction_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.png\"\n", 
		gamma, e, system.name, K, po.initial_condition[0], po.initial_condition[1], po.period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "set xlabel \"{/Symbol q}\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "unset colorbox\n");
	fprintf(gnuplotPipe, "set xrange[-3.1415:3.1415]\n");
	fprintf(gnuplotPipe, "set yrange [0.0:3.0]\n");
	fprintf(gnuplotPipe, "unset key\n");
	// fprintf(gnuplotPipe, 
	// 	"set title \"Basin of attraction  for  gamma = %1.3f  e = %1.3f  K = %1.0e  res = %d  n = %1.0e  eps = %1.0e\"\n", 
	// 	gamma, e, K, analysis.grid_resolution, (double)analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "FILE = \"basin_size_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat\"\n", 
		gamma, e, system.name, K, po.initial_condition[0], po.initial_condition[1], po.period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "stats FILE u (myTitle=strcol(1),0) nooutput\n");
	fprintf(gnuplotPipe, "set title \"gamma = %1.3f  e = %1.3f  K = %1.0e  res = %d  n = %1.0e  eps = %1.0e  size = \".myTitle\n", 
		gamma, e, K, analysis.grid_resolution, (double)analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "plot 'basin_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat' u 1:2:3 w image notitle, 'basin_ref_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat' w p pt 7 ps 1.5 lc rgb \"green\" notitle",
		gamma, e, system.name, K, po.initial_condition[0], po.initial_condition[1], po.period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps, gamma, e, system.name, K, po.initial_condition[0], po.initial_condition[1], po.period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fclose(gnuplotPipe);

	printf("Done!\n");
	printf("Data written in output/basin_of_attraction/\n");

	printf("Drawing convergence times of system %s with gamma = %1.3f, e = %1.3f and K = %1.5f and reference theta = %1.3f theta_dot = %1.3f with period = %d\n", 
		system.name, gamma, e, K, po.initial_condition[0], po.initial_condition[1], po.period);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/basin_of_attraction/fig_convergence_times_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.png\"\n", 
		gamma, e, system.name, K, po.initial_condition[0], po.initial_condition[1], po.period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "set xlabel \"{/Symbol q}\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "set xrange[-3.1415:3.1415]\n");
	fprintf(gnuplotPipe, "set yrange [0.0:3.0]\n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, 
		"set title \"        Convergence times  for  gamma = %1.3f  e = %1.3f  K = %1.0e  res = %d  n = %1.0e  eps = %1.0e\"\n", 
		gamma, e, K, analysis.grid_resolution, (double)analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "plot 'basin_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat' u 1:2:4 w image notitle, 'basin_ref_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat' w p pt 7 ps 1.5 lc rgb \"green\" notitle",
		gamma, e, system.name, K, po.initial_condition[0], po.initial_condition[1], po.period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps, gamma, e, system.name, K, po.initial_condition[0], po.initial_condition[1], po.period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fclose(gnuplotPipe);

	printf("Done!\n");
	printf("Data written in output/basin_of_attraction/\n");

	return 0;
}

int draw_basin_of_attraction_clean	(int ref_period, double ref[][2],
                            		 dynsys system, anlsis analysis)
{

	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/clean_figures", &st) == -1) {
		mkdir("output/clean_figures", 0700);
	}

	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	printf("Drawing basin of attraction of system %s with gamma = %1.3f, e = %1.3f and K = %1.5f and reference theta = %1.3f theta_dot = %1.3f with period = %d\n", 
		system.name, gamma, e, K, ref[0][0], ref[0][1], ref_period);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/clean_figures/fig_basin_of_attraction_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.png\"\n", 
		gamma, e, system.name, K, ref[0][0], ref[0][1], ref_period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,50\"\n");
	fprintf(gnuplotPipe, "set size square \n");
	fprintf(gnuplotPipe, "set lmargin at screen 0.05\n");
	fprintf(gnuplotPipe, "set bmargin at screen 0.05\n");
	fprintf(gnuplotPipe, "set rmargin at screen 0.95\n");
	fprintf(gnuplotPipe, "set tmargin at screen 0.95\n");
	fprintf(gnuplotPipe, "set border lw 2 \n");
	fprintf(gnuplotPipe, "set xrange[-3.1415:3.1415]\n");
	fprintf(gnuplotPipe, "set yrange [0.0:3.0]\n");
	fprintf(gnuplotPipe, "set xtics format \" \"\n");
	fprintf(gnuplotPipe, "set ytics format \" \"\n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, "unset title\n");
	fprintf(gnuplotPipe, "unset colorbox\n");
	fprintf(gnuplotPipe, "plot 'basin_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat' u 1:2:3 w image notitle, 'basin_ref_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat' w p pt 7 ps 1.5 lc rgb \"green\" notitle",
		gamma, e, system.name, K, ref[0][0], ref[0][1], ref_period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps, gamma, e, system.name, K, ref[0][0], ref[0][1], ref_period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fclose(gnuplotPipe);

	printf("Done!\n");
	printf("Data written in output/clean_figures/\n");

	printf("Drawing convergence times of system %s with gamma = %1.3f, e = %1.3f and K = %1.5f and reference theta = %1.3f theta_dot = %1.3f with period = %d\n", 
		system.name, gamma, e, K, ref[0][0], ref[0][1], ref_period);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/clean_figures/fig_convergence_times_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.png\"\n", 
		gamma, e, system.name, K, ref[0][0], ref[0][1], ref_period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "set terminal pngcairo size 2400,2000 font \"fonts/cmr10.ttf,50\"\n");
	fprintf(gnuplotPipe, "set size square \n");
	fprintf(gnuplotPipe, "set lmargin at screen 0.05\n");
	fprintf(gnuplotPipe, "set bmargin at screen 0.05\n");
	fprintf(gnuplotPipe, "set tmargin at screen 0.95\n");
	fprintf(gnuplotPipe, "set border lw 2 \n");
	fprintf(gnuplotPipe, "set xrange[-3.1415:3.1415]\n");
	fprintf(gnuplotPipe, "set yrange [0.0:3.0]\n");
	fprintf(gnuplotPipe, "set xtics format \" \"\n");
	fprintf(gnuplotPipe, "set ytics format \" \"\n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, "unset title\n");
	fprintf(gnuplotPipe, "plot 'basin_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat' u 1:2:4 w image notitle, 'basin_ref_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat' w p pt 7 ps 1.5 lc rgb \"green\" notitle",
		gamma, e, system.name, K, ref[0][0], ref[0][1], ref_period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps, gamma, e, system.name, K, ref[0][0], ref[0][1], ref_period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fclose(gnuplotPipe);

	printf("Done!\n");
	printf("Data written in output/clean_figures/\n");

	return 0;
}

int draw_multiple_basin_of_attraction_determined(dynsys system,
                                        		 anlsis analysis)
{
	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	int cb_range_min = cantor_pairing_function(analysis.spin_period_min, analysis.orbit_period_min);
	int cb_range_max = cantor_pairing_function(analysis.spin_period_max, analysis.orbit_period_max);

	printf("Drawing multtiple basin of attraction of system %s with gamma = %1.6f, e = %1.3f and K = %1.5f\n", 
		system.name, gamma, e, K);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/basin_of_attraction/fig_multiple_basin_of_attraction_determined_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.png\"\n", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "set xlabel \"{/Symbol q}\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "unset colorbox\n");
	fprintf(gnuplotPipe, "set xrange[-3.15:3.15]\n");
	fprintf(gnuplotPipe, "set yrange [%f:%f]\n", analysis.grid_velocity_min, analysis.grid_velocity_max);
	// fprintf(gnuplotPipe, "set cbrange [%d:%d]\n", cb_range_min, cb_range_max);
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, 
		"set title \"Multiple basin of attraction (D) for {/Symbol g} = %1.6f e = %1.3f K = %1.0e res = %d n = %1.0e {/Symbol e} = %1.0e\"\n", 
		gamma, e, K, analysis.grid_resolution, (double)analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "plot 'multiple_basin_determined_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat' u 1:2:3 w image notitle",
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, ", 'multiple_basin_determined_ref_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat' u 1:2 w p pt 7 ps 2 lc rgb \"green\" title \"spin-orbit resonances\"",
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fclose(gnuplotPipe);

	printf("Done!\n");
	printf("Data written in output/basin_of_attraction/\n");

	// printf("Drawing convergence times for system %s with gamma = %1.3f, e = %1.3f and K = %1.5f\n", 
	// 	system.name, gamma, e, K);

	// gnuplotPipe = popen("gnuplot -persistent", "w");
	// fprintf(gnuplotPipe, "reset\n");
	// fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	// fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	// fprintf(gnuplotPipe, 
	// 	"set output \"output/basin_of_attraction/fig_multiple_convergence_time_determined_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.png\"\n", 
	// 	gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	// fprintf(gnuplotPipe, "set xlabel \"{/Symbol q}\"\n");
	// fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	// fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	// fprintf(gnuplotPipe, "unset colorbox\n");
	// fprintf(gnuplotPipe, "set xrange[-3.15:3.15]\n");
	// fprintf(gnuplotPipe, "set yrange [0.0:3.0]\n");
	// fprintf(gnuplotPipe, "unset key\n");
	// fprintf(gnuplotPipe, 
	// 	"set title \"Convergence times (D) for {/Symbol g} = %1.3f e = %1.3f K = %1.0e res = %d n = %1.0e {/Symbol e} = %1.0e\"\n", 
	// 	gamma, e, K, analysis.grid_resolution, (double)analysis.number_of_cycles, analysis.evolve_basin_eps);
	// fprintf(gnuplotPipe, "plot 'multiple_basin_determined_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat' u 1:2:4 w image notitle",
	// 	gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	// fprintf(gnuplotPipe, ", 'multiple_basin_determined_ref_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat' u 1:2 w p pt 7 ps 2 lc rgb \"green\" title \"spin-orbit resonances\"",
	// 	gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	// fclose(gnuplotPipe);

	// printf("Done!\n");
	// printf("Data written in output/basin_of_attraction/\n");

	return 0;
}

int draw_multiple_basin_of_attraction_determined_clean	(dynsys system,
                                        		 		 anlsis analysis)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/clean_figures", &st) == -1) {
		mkdir("output/clean_figures", 0700);
	}
	
	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	int cb_range_min = cantor_pairing_function(analysis.spin_period_min, analysis.orbit_period_min);
	int cb_range_max = cantor_pairing_function(analysis.spin_period_max, analysis.orbit_period_max);

	printf("Drawing multtiple basin of attraction of system %s with gamma = %1.3f, e = %1.3f and K = %1.5f\n", 
		system.name, gamma, e, K);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
		fprintf(gnuplotPipe, "set terminal pngcairo size 2000,2000 font \"fonts/cmr10.ttf,50\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/clean_figures/fig_multiple_basin_of_attraction_determined_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.png\"\n", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "set lmargin at screen 0.05\n");
	fprintf(gnuplotPipe, "set bmargin at screen 0.05\n");
	fprintf(gnuplotPipe, "set rmargin at screen 0.95\n");
	fprintf(gnuplotPipe, "set tmargin at screen 0.95\n");
	fprintf(gnuplotPipe, "set border lw 2 \n");
	fprintf(gnuplotPipe, "set xrange[-3.1415:3.1415]\n");
	fprintf(gnuplotPipe, "set yrange [0.0:3.0]\n");
	fprintf(gnuplotPipe, "set xtics format \" \"\n");
	fprintf(gnuplotPipe, "set ytics format \" \"\n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, "unset title\n");
	fprintf(gnuplotPipe, "unset colorbox\n");
	fprintf(gnuplotPipe, "set cbrange [%d:%d]\n", cb_range_min, cb_range_max);
	fprintf(gnuplotPipe, "plot 'multiple_basin_determined_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat' u 1:2:3 w image notitle",
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, ", 'multiple_basin_determined_ref_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat' u 1:2 w p pt 7 ps 2 lc rgb \"green\" title \"spin-orbit resonances\"",
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fclose(gnuplotPipe);

	printf("Done!\n");
	printf("Data written in output/basin_of_attraction/\n");

	return 0;
}

int plot_size_multiple_basin_of_attraction_determined_range_e	(int number_of_e,
																 double e_initial,
																 double e_final,
																 dynsys system,
																 anlsis analysis)
{
	FILE 	*combineFiles;
	FILE 	*gnuplotPipe;
	int		size_filename = 200;
	char 	local_filename[size_filename];
	char	filename[size_filename * number_of_e + 50];
	int		spin_period;
	int		orbit_period;
	int		cb_range_min;
	int		cb_range_max;
	double 	*par = (double *)system.params;
	double 	gamma = par[0];
	double 	K = par[6];
	double 	ec, e_step;

	if (number_of_e < 2)
	{
		printf("Invalid range\n");
		exit(2);
	}

	e_step = (e_final - e_initial) / (double)(number_of_e); 

	sprintf(filename, "paste -d \"\n\"");

	ec = e_initial;
	for(int i = 0; i <= number_of_e; i++)
	{
		sprintf(local_filename, " output/basin_of_attraction/multiple_basin_determined_size_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat", 
		gamma, ec, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
		strcat(filename, local_filename);
		ec += e_step;
	}

	strcat(filename, " > output/basin_of_attraction/basins_size_combined.dat");

	combineFiles = popen(filename, "w");

	fclose(combineFiles);

	cb_range_min = cantor_pairing_function(analysis.spin_period_min, analysis.orbit_period_min);
	cb_range_max = cantor_pairing_function(analysis.spin_period_max, analysis.orbit_period_max);

	gnuplotPipe = popen("gnuplot -persistent", "w");

	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/basin_of_attraction/fig_size_multiple_basin_of_attraction_determined_range_e_gamma_%1.6f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.png\"\n", 
		gamma, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "set xlabel \"Orbital eccentricity e\"\n");
	fprintf(gnuplotPipe, "set ylabel \"Basin size\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	// fprintf(gnuplotPipe, 
	// 	"set title \"       Size of multiple basin of attraction (D) for {/Symbol g} = %1.3f K = %1.0e res = %d n = %1.0e {/Symbol e} = %1.0e\"\n", 
	// 	gamma, K, analysis.grid_resolution, (double)analysis.number_of_cycles, analysis.evolve_basin_eps);

	fprintf(gnuplotPipe, "set key opaque box outside top right width 1.1\n");

	fprintf(gnuplotPipe, "unset colorbox\n");
	fprintf(gnuplotPipe, "set cbrange [%d:%d]\n", cb_range_min, cb_range_max);

	fprintf(gnuplotPipe, "plot ");

	for (orbit_period = analysis.orbit_period_min; orbit_period <= analysis.orbit_period_max; orbit_period++)
	{
		for (spin_period = analysis.spin_period_min; spin_period <= analysis.spin_period_max; spin_period++)
		{
			fprintf(gnuplotPipe, "'basins_size_combined.dat' u 1:($2==%d&&$3==%d?$4:1/0) w lp pt %d ps 2 title \"%d/%d\", ", 
				spin_period, orbit_period, orbit_period + 4, spin_period, orbit_period);
		}
	}

	// fprintf(gnuplotPipe, "'basins_size_combined.dat' u 1:(strcol(2) eq \"s\"?$4:1/0) w lp pt 7 ps 2 lc rgb \"black\" title \"sum\"");

	fprintf(gnuplotPipe, "'basins_size_combined.dat' u 1:($2==0&&$3==0?$4:1/0) w lp pt 7 ps 2 lc rgb \"black\" title \"HO,QP\"");

	fclose(gnuplotPipe);

	printf("Done!\n");

	return 0;
}

int plot_size_multiple_basin_of_attraction_determined_range_e_latex	(int number_of_e,
																 	 double e_initial,
																 	 double e_final,
																 	 dynsys system,
																 	 anlsis analysis)
{
	FILE 	*combineFiles;
	FILE 	*gnuplotPipe;
	int		size_filename = 200;
	char 	local_filename[size_filename];
	char	filename[size_filename * number_of_e + 50];
	int		spin_period;
	int		orbit_period;
	int		cb_range_min;
	int		cb_range_max;
	double 	*par = (double *)system.params;
	double 	gamma = par[0];
	double 	K = par[6];
	double 	ec, e_step;

	if (number_of_e < 2)
	{
		printf("Invalid range\n");
		exit(2);
	}

	e_step = (e_final - e_initial) / (double)(number_of_e); 

	sprintf(filename, "paste -d \"\n\"");

	ec = e_initial;
	for(int i = 0; i <= number_of_e; i++)
	{
		sprintf(local_filename, " output/basin_of_attraction/multiple_basin_determined_size_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat", 
		gamma, ec, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
		strcat(filename, local_filename);
		ec += e_step;
	}

	strcat(filename, " > output/basin_of_attraction/basins_size_combined.dat");

	combineFiles = popen(filename, "w");

	fclose(combineFiles);

	cb_range_min = cantor_pairing_function(analysis.spin_period_min, analysis.orbit_period_min);
	cb_range_max = cantor_pairing_function(analysis.spin_period_max, analysis.orbit_period_max);

	gnuplotPipe = popen("gnuplot -persistent", "w");

	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 1000,1000 font \"fonts/cmr10.ttf,25\"\n");
	fprintf(gnuplotPipe, "set key font \"fonts/cmr10.ttf,18\" \n");
	fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/basin_of_attraction/fig_size_multiple_basin_of_attraction_determined_range_e_gamma_%1.6f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f_latex.png\"\n", 
		gamma, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "set size square \n");
	fprintf(gnuplotPipe, "set border lw 2 \n");
	fprintf(gnuplotPipe, "set xlabel \"Orbital eccentricity e\"\n");
	fprintf(gnuplotPipe, "set ylabel \"Basin size\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "set key opaque box lw 2 outside top right width 1.1\n");
	fprintf(gnuplotPipe, "set bmargin at screen 0.15\n");
	fprintf(gnuplotPipe, "unset title\n");

	fprintf(gnuplotPipe, "unset colorbox\n");
	fprintf(gnuplotPipe, "set cbrange [%d:%d]\n", cb_range_min, cb_range_max);

	fprintf(gnuplotPipe, "plot ");

	for (orbit_period = analysis.orbit_period_min; orbit_period <= analysis.orbit_period_max; orbit_period++)
	{
		for (spin_period = analysis.spin_period_min; spin_period <= analysis.spin_period_max; spin_period++)
		{
			fprintf(gnuplotPipe, "'basins_size_combined.dat' u 1:($2==%d&&$3==%d?$4:1/0) w lp pt %d ps 2 title \"%d/%d\", ", 
				spin_period, orbit_period, orbit_period + 4, spin_period, orbit_period);
		}
	}

	fprintf(gnuplotPipe, "'basins_size_combined.dat' u 1:(strcol(2) eq \"s\"?$4:1/0) w lp pt 7 ps 2 lc rgb \"black\" title \"sum\"");

	fclose(gnuplotPipe);

	printf("Done!\n");

	return 0;
}

int draw_multiple_basin_of_attraction_undetermined	(dynsys system,
                                        		 	 anlsis analysis)
{
	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	printf("Drawing basin of attraction of system %s with gamma = %1.3f, e = %1.3f and K = %1.5f\n", 
		system.name, gamma, e, K);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/basin_of_attraction/fig_multiple_basin_of_attraction_undetermined_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d.png\"\n", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles);
	fprintf(gnuplotPipe, "set xlabel \"{/Symbol q}\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "unset colorbox\n");
	fprintf(gnuplotPipe, "set xrange[-3.1415:3.1415]\n");
	fprintf(gnuplotPipe, "set yrange [0.0:3.0]\n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, 
		"set title \"Multiple basin of attraction (U) for {/Symbol g} = %1.3f e = %1.3f K = %1.0e res = %d n = %1.0e\"\n", 
		gamma, e, K, analysis.grid_resolution, (double)analysis.number_of_cycles);
	fprintf(gnuplotPipe, "plot 'multiple_basin_undetermined_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d.dat' u 1:2:3 w image notitle",
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles);
	fprintf(gnuplotPipe, ", 'multiple_basin_undetermined_ref_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d.dat' u 1:2 w p pt 7 ps 2 lc rgb \"green\" title \"spin-orbit resonances\"",
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles);
	fclose(gnuplotPipe);

	printf("Done!\n");
	printf("Data written in output/basin_of_attraction/\n");

	printf("Drawing convergence times for system %s with gamma = %1.3f, e = %1.3f and K = %1.5f\n", 
		system.name, gamma, e, K);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/basin_of_attraction/fig_multiple_convergence_time_undetermined_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d.png\"\n", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles);
	fprintf(gnuplotPipe, "set xlabel \"{/Symbol q}\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "unset colorbox\n");
	fprintf(gnuplotPipe, "set xrange[-3.1415:3.1415]\n");
	fprintf(gnuplotPipe, "set yrange [0.0:3.0]\n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, 
		"set title \"Convergence times (U) for {/Symbol g} = %1.3f e = %1.3f K = %1.0e res = %d n = %1.0e\"\n", 
		gamma, e, K, analysis.grid_resolution, (double)analysis.number_of_cycles);
	fprintf(gnuplotPipe, "plot 'multiple_basin_undetermined_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d.dat' u 1:2:4 w image notitle",
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles);
	fprintf(gnuplotPipe, ", 'multiple_basin_undetermined_ref_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d.dat' u 1:2 w p pt 7 ps 2 lc rgb \"green\" title \"spin-orbit resonances\"",
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles);
	fclose(gnuplotPipe);

	printf("Done!\n");
	printf("Data written in output/basin_of_attraction/\n");

	return 0;
}

int plot_basin_entropy_vs_box_size	(dynsys system,
                         			 anlsis analysis)
{
	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	int cb_range_min = cantor_pairing_function(analysis.spin_period_min, analysis.orbit_period_min);
	int cb_range_max = cantor_pairing_function(analysis.spin_period_max, analysis.orbit_period_max);

	printf("Drawing basin entropy for system %s with gamma = %1.3f, e = %1.3f and K = %1.5f\n", 
		system.name, gamma, e, K);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/basin_of_attraction/fig_basin_entropy_vs_box_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.png\"\n", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "set xlabel \"log {/Symbol e}\"\n");
	fprintf(gnuplotPipe, "set ylabel \"log S/N\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	// fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, "set key opaque box top left\n");
	fprintf(gnuplotPipe, "set fit quiet\n");
	fprintf(gnuplotPipe, "set fit logfile '/dev/null'\n");
	fprintf(gnuplotPipe, "f(x) = a * x + b\n");
	fprintf(gnuplotPipe, "fit f(x) 'multiple_basin_determined_entropy_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat' u 1:2 via a,b\n",
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "set print \"output/basin_of_attraction/basin_entropy_slope_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat\"\n",
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "print sprintf(\"%1.3f %%1.3f\", a)\n", e);
	fprintf(gnuplotPipe, "set print\n");
	fprintf(gnuplotPipe, 
		"set title \"Basin entropy for {/Symbol g} = %1.3f e = %1.3f K = %1.0e res = %d n = %1.0e {/Symbol e} = %1.0e\"\n", 
		gamma, e, K, analysis.grid_resolution, (double)analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "plot 'multiple_basin_determined_entropy_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat' u 1:2 w p pt 3 ps 2 notitle, f(x) lw 2 title sprintf(\"Slope = %%1.3f\", a)",
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fclose(gnuplotPipe);

	printf("Done!\n");
	printf("Data written in output/basin_of_attraction/\n");

	return 0;
}

int plot_slope_basin_entropy_range_e(int number_of_e,
									 double e_initial,
									 double e_final,
									 dynsys system,
									 anlsis analysis)
{
	FILE 	*combineFiles;
	FILE 	*gnuplotPipe;
	int		size_filename = 200;
	char 	local_filename[size_filename];
	char	filename[size_filename * number_of_e + 50];
	int		spin_period;
	int		orbit_period;
	double 	*par = (double *)system.params;
	double 	gamma = par[0];
	double 	K = par[6];
	double 	ec, e_step;

	printf("Drawing basin entropy slope\n");

	if (number_of_e < 2)
	{
		printf("Invalid range\n");
		exit(2);
	}

	e_step = (e_final - e_initial) / (double)(number_of_e); 

	sprintf(filename, "paste -s ");

	ec = e_initial;
	for(int i = 0; i <= number_of_e; i++)
	{
		sprintf(local_filename, " output/basin_of_attraction/basin_entropy_slope_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat", 
		gamma, ec, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
		strcat(filename, local_filename);
		ec += e_step;
	}

	strcat(filename, " > output/basin_of_attraction/basin_entropy_slope_combined.dat");

	combineFiles = popen(filename, "w");

	fclose(combineFiles);

	gnuplotPipe = popen("gnuplot -persistent", "w");

	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/basin_of_attraction/fig_basin_entropy_slope_range_e_gamma_%1.6f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.png\"\n", 
		gamma, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "set xlabel \"Orbital eccentricity e\"\n");
	fprintf(gnuplotPipe, "set ylabel \"Basin entropy slope\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, 
		"set title \"       Basin entropy slope for {/Symbol g} = %1.3f K = %1.0e res = %d n = %1.0e {/Symbol e} = %1.0e\"\n", 
		gamma, K, analysis.grid_resolution, (double)analysis.number_of_cycles, analysis.evolve_basin_eps);

	fprintf(gnuplotPipe, "unset key\n");

	fprintf(gnuplotPipe, "plot 'basin_entropy_slope_combined.dat' w lp pt 7 ps 1.5");

	fclose(gnuplotPipe);

	printf("Done!\n");

	return 0;
}

int plot_basin_entropy_range_e	(int number_of_e,
								 double e_initial,
								 double e_final,
								 dynsys system,
								 anlsis analysis)
{
	FILE 	*combineFiles;
	FILE 	*gnuplotPipe;
	int		size_filename = 200;
	char 	local_filename[size_filename];
	char	filename[size_filename * number_of_e + 50];
	int		spin_period;
	int		orbit_period;
	double 	*par = (double *)system.params;
	double 	gamma = par[0];
	double 	K = par[6];
	double 	ec, e_step;

	printf("Drawing basin entropy range e\n");

	if (number_of_e < 2)
	{
		printf("Invalid range\n");
		exit(2);
	}

	e_step = (e_final - e_initial) / (double)(number_of_e); 

	sprintf(filename, "paste -s");

	ec = e_initial;
	for(int i = 0; i <= number_of_e; i++)
	{
		sprintf(local_filename, " output/basin_of_attraction/multiple_basin_determined_entropy_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat", 
		gamma, ec, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
		strcat(filename, local_filename);
		ec += e_step;
	}

	strcat(filename, " > output/basin_of_attraction/entropy_size_combined.dat");

	combineFiles = popen(filename, "w");

	fclose(combineFiles);

	gnuplotPipe = popen("gnuplot -persistent", "w");

	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/basin_of_attraction/fig_basin_entropy_size_range_e_gamma_%1.6f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.png\"\n", 
		gamma, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "set xlabel \"Orbital eccentricity e\"\n");
	fprintf(gnuplotPipe, "set ylabel \"Basin entropy\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, 
		"set title \"       Basin entropy for {/Symbol g} = %1.3f K = %1.0e res = %d n = %1.0e {/Symbol e} = %1.0e\"\n", 
		gamma, K, analysis.grid_resolution, (double)analysis.number_of_cycles, analysis.evolve_basin_eps);

	fprintf(gnuplotPipe, "unset key\n");

	fprintf(gnuplotPipe, "plot 'entropy_size_combined.dat' w lp pt 7 ps 1.5,");
	fprintf(gnuplotPipe, " 'entropy_size_combined.dat' u 1:3 w lp pt 7 ps 1.5");

	fclose(gnuplotPipe);

	printf("Done!\n");

	return 0;
}

int plot_size_multiple_basin_of_attraction_determined_plus_basin_entropy_range_e(int number_of_e,
																 	 			 double e_initial,
																 	 			 double e_final,
																 	 			 dynsys system,
																 	 			 anlsis analysis)
{
	FILE 	*combineFiles;
	FILE 	*gnuplotPipe;
	int		size_filename = 200;
	char 	local_filename[size_filename];
	char	filename[size_filename * number_of_e + 50];
	int		spin_period;
	int		orbit_period;
	int		cb_range_min;
	int		cb_range_max;
	double 	*par = (double *)system.params;
	double 	gamma = par[0];
	double 	K = par[6];
	double 	ec, e_step;

	if (number_of_e < 2)
	{
		printf("Invalid range\n");
		exit(2);
	}

	e_step = (e_final - e_initial) / (double)(number_of_e); 

	sprintf(filename, "paste -d \"\n\"");

	ec = e_initial;
	for(int i = 0; i <= number_of_e; i++)
	{
		sprintf(local_filename, "output/basin_of_attraction/multiple_basin_determined_size_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat", 
		gamma, ec, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
		FILE *file_verify = fopen(local_filename, "r");
		if (file_verify != NULL)
		{
			strcat(filename, " ");
			strcat(filename, local_filename);
			fclose(file_verify);
		}
		ec += e_step;
	}

	strcat(filename, " > output/basin_of_attraction/basins_size_combined.dat");

	combineFiles = popen(filename, "w");

	fclose(combineFiles);

	sprintf(filename, "paste -s");

	ec = e_initial;
	for(int i = 0; i <= number_of_e; i++)
	{
		sprintf(local_filename, " output/basin_of_attraction/multiple_basin_determined_entropy_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat", 
		gamma, ec, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
		strcat(filename, local_filename);
		ec += e_step;
	}

	strcat(filename, " > output/basin_of_attraction/entropy_size_combined.dat");

	combineFiles = popen(filename, "w");

	fclose(combineFiles);

	cb_range_min = cantor_pairing_function(analysis.spin_period_min, analysis.orbit_period_min);
	cb_range_max = cantor_pairing_function(analysis.spin_period_max, analysis.orbit_period_max);

	gnuplotPipe = popen("gnuplot -persistent", "w");

	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	// fprintf(gnuplotPipe, 
	// 	"set output \"output/basin_of_attraction/fig_size_multiple_basin_of_attraction_determined_with_basin_entropy_range_e_gamma_%1.6f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.png\"\n", 
	// 	gamma, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, 
		"set output \"output/basin_of_attraction/fig_size_multiple_basin_of_attraction_determined_with_basin_entropy_synchronous_range_e_gamma_%1.6f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.png\"\n", 
		gamma, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "set xlabel \"Orbital eccentricity e\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	// fprintf(gnuplotPipe, 
	// 	"set title \"       Size of multiple basin of attraction (D) for {/Symbol g} = %1.3f K = %1.0e res = %d n = %1.0e {/Symbol e} = %1.0e\"\n", 
	// 	gamma, K, analysis.grid_resolution, (double)analysis.number_of_cycles, analysis.evolve_basin_eps);

	fprintf(gnuplotPipe, "set key opaque box outside top right width 1.1\n");

	fprintf(gnuplotPipe, "unset colorbox\n");
	fprintf(gnuplotPipe, "set cbrange [%d:%d]\n", cb_range_min, cb_range_max);

	fprintf(gnuplotPipe, "plot ");

	double M[analysis.spin_period_max + 1][analysis.orbit_period_max + 1];
	for (int j = 0; j <= analysis.orbit_period_max; j++)
	{
		for (int i = 0; i <= analysis.spin_period_max; i++)
		{
			M[i][j] = 0.0;
		}
	}

	FILE *in_combined = fopen("output/basin_of_attraction/basins_size_combined.dat", "r");
	int sp, or, fool_cantor;
	double fool_e, basin_size;
	while(fscanf(in_combined, "%lf %d %d %lf %d\n", &fool_e, &sp, &or, &basin_size, &fool_cantor) != EOF)
	{
		M[sp][or] += basin_size;
	}
	fclose(in_combined);

	// orbit_period = 1;
	// spin_period = 1;
	for (orbit_period = analysis.orbit_period_min; orbit_period <= analysis.orbit_period_max; orbit_period++)
	{
		for (spin_period = analysis.spin_period_min; spin_period <= analysis.spin_period_max; spin_period++)
		{
			// printf("%d %d %f\n", spin_period, orbit_period, M[spin_period][orbit_period]);
			if (M[spin_period][orbit_period] > 0.0)
			{
				fprintf(gnuplotPipe, "'basins_size_combined.dat' u 1:($2==%d&&$3==%d?$4:1/0) w lp pt %d ps 2 title \"%d/%d\", ", 
					spin_period, orbit_period, orbit_period + 4, spin_period, orbit_period);
			}
		}
	}

	fprintf(gnuplotPipe, "'basins_size_combined.dat' u 1:($2==0&&$3==0?$4:1/0) w lp pt 7 ps 2 lc rgb \"red\" title \"HO,QP\", ");

	fprintf(gnuplotPipe, "'entropy_size_combined.dat' u 1:2 w l lw 2 lc rgb \"black\" title \"NBE\", ");

	fclose(gnuplotPipe);

	printf("Done!\n");

	return 0;
}

int plot_size_multiple_basin_of_attraction_determined_plus_basin_entropy_monte_carlo_with_break_range_e(int number_of_e,
																										double e_step_ext,
																 	 			 						double e_initial,
																 	 			 						double e_final,
																 	 			 						dynsys system,
																 	 			 						anlsis analysis)
{
	FILE 	*combineFiles;
	FILE 	*gnuplotPipe;
	int		size_filename = 300;
	char 	local_filename[size_filename];
	char	filename[50000];
	int		spin_period;
	int		orbit_period;
	int		cb_range_min;
	int		cb_range_max;
	double 	*par = (double *)system.params;
	double 	gamma = par[0];
	double 	K = par[6];
	double 	ec, e_step;

	if (number_of_e > 0)
	{
		e_step = (e_final - e_initial) / (double)(number_of_e);
	}
	else
	{
		e_step = e_step_ext;
	}
	ec = e_initial;
	sprintf(filename, "paste -d \"\n\"");
	while(ec < e_final + e_step/2.0)
	{
		sprintf(local_filename, " output/basin_of_attraction/multiple_basin_determined_size_monte_carlo_with_break_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.dat", 
		gamma, ec, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc);
		strcat(filename, local_filename);
		ec += e_step;
	}
	strcat(filename, " > output/basin_of_attraction/basins_size_combined_monte_carlo_with_break.dat");

	combineFiles = popen(filename, "w");

	fclose(combineFiles);

	sprintf(filename, "paste -s");

	ec = e_initial;
	while(ec < e_final + e_step/2.0)
	{
		sprintf(local_filename, " output/basin_of_attraction/multiple_basin_determined_entropy_monte_carlo_with_break_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.dat", 
		gamma, ec, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc);
		strcat(filename, local_filename);
		ec += e_step;
	}

	strcat(filename, " > output/basin_of_attraction/entropy_size_combined_monte_carlo_with_break.dat");

	combineFiles = popen(filename, "w");

	fclose(combineFiles);

	cb_range_min = cantor_pairing_function(analysis.spin_period_min, analysis.orbit_period_min);
	cb_range_max = cantor_pairing_function(analysis.spin_period_max, analysis.orbit_period_max);

	gnuplotPipe = popen("gnuplot -persistent", "w");

	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/basin_of_attraction/fig_size_multiple_basin_of_attraction_determined_with_basin_entropy__monte_carlo_with_break_range_e_gamma_%1.6f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.png\"\n", 
		gamma, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	// fprintf(gnuplotPipe, 
	// 	"set output \"output/basin_of_attraction/fig_size_multiple_basin_of_attraction_determined_with_basin_entropy_monte_carlo_with_break_synchronous_range_e_gamma_%1.6f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.png\"\n", 
	// 	gamma, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "set xlabel \"Orbital eccentricity e\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	// fprintf(gnuplotPipe, 
	// 	"set title \"       Size of multiple basin of attraction (D) for {/Symbol g} = %1.3f K = %1.0e res = %d n = %1.0e {/Symbol e} = %1.0e\"\n", 
	// 	gamma, K, analysis.grid_resolution, (double)analysis.number_of_cycles, analysis.evolve_basin_eps);

	fprintf(gnuplotPipe, "set key opaque box outside top right width 1.1\n");

	fprintf(gnuplotPipe, "unset colorbox\n");
	fprintf(gnuplotPipe, "set cbrange [%d:%d]\n", cb_range_min, cb_range_max);

	fprintf(gnuplotPipe, "plot ");

	// orbit_period = 1;
	// spin_period = 1;
	for (orbit_period = analysis.orbit_period_min; orbit_period <= analysis.orbit_period_max; orbit_period++)
	{
		for (spin_period = analysis.spin_period_min; spin_period <= analysis.spin_period_max; spin_period++)
		{
			fprintf(gnuplotPipe, "'basins_size_combined_monte_carlo_with_break.dat' u 1:($2==%d&&$3==%d?$4:1/0) w lp pt %d ps 2 title \"%d/%d\", ", 
				spin_period, orbit_period, orbit_period + 4, spin_period, orbit_period);
		}
	}

	fprintf(gnuplotPipe, "'basins_size_combined_monte_carlo_with_break.dat' u 1:($2==0&&$3==0?$4:1/0) w lp pt 7 ps 2 lc rgb \"black\" title \"HO,QP\", ");

	fprintf(gnuplotPipe, "'entropy_size_combined_monte_carlo_with_break.dat' u 1:2 w l lw 2 title \"NBE\", ");

	fclose(gnuplotPipe);

	printf("Done!\n");

	return 0;
}

int plot_size_multiple_basin_of_attraction_undetermined_plus_basin_entropy_monte_carlo_with_break_range_e	(int number_of_e,
																											 double e_step_ext,
																 	 			 							 double e_initial,
																 	 			 							 double e_final,
																 	 			 							 dynsys system,
																 	 			 							 anlsis analysis)
{
	FILE 	*combineFiles;
	FILE 	*gnuplotPipe;
	int		size_filename = 300;
	char 	local_filename[size_filename];
	char	filename[50000];
	int		spin_period;
	int		orbit_period;
	int		cb_range_min;
	int		cb_range_max;
	double 	*par = (double *)system.params;
	double 	gamma = par[0];
	double 	K = par[6];
	double 	ec, e_step;

	if (number_of_e > 0)
	{
		e_step = (e_final - e_initial) / (double)(number_of_e);
	}
	else
	{
		e_step = e_step_ext;
	}
	ec = e_initial;
	sprintf(filename, "paste -d \"\n\"");
	while(ec < e_final + e_step/2.0)
	{
		sprintf(local_filename, " output/basin_of_attraction/multiple_basin_undetermined_size_monte_carlo_with_break_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_rand_%d_window_mc_%d_precision_mc_%1.3f_transient_wn_%d_window_wn_%d_precision_wn_%1.3f.dat", 
		gamma, ec, system.name, K, analysis.number_of_cycles, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc, analysis.convergence_transient_wn, analysis.convergence_window_mc, analysis.convergence_precision_mc);
		strcat(filename, local_filename);
		ec += e_step;
	}
	strcat(filename, " > output/basin_of_attraction/basins_size_undetermined_combined_monte_carlo_with_break.dat");

	combineFiles = popen(filename, "w");

	fclose(combineFiles);

	sprintf(filename, "paste -s");

	ec = e_initial;
	while(ec < e_final + e_step/2.0)
	{
		sprintf(local_filename, " output/basin_of_attraction/multiple_basin_undetermined_entropy_monte_carlo_with_break_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_rand_%d_window_mc_%d_precision_mc_%1.3f_transient_wn_%d_window_wn_%d_precision_wn_%1.3f.dat", 
		gamma, ec, system.name, K, analysis.number_of_cycles, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc, analysis.convergence_transient_wn, analysis.convergence_window_mc, analysis.convergence_precision_mc);
		strcat(filename, local_filename);
		ec += e_step;
	}

	strcat(filename, " > output/basin_of_attraction/entropy_size_undetermined_combined_monte_carlo_with_break.dat");

	combineFiles = popen(filename, "w");

	fclose(combineFiles);

	cb_range_min = cantor_pairing_function(analysis.spin_period_min, analysis.orbit_period_min);
	cb_range_max = cantor_pairing_function(analysis.spin_period_max, analysis.orbit_period_max);

	gnuplotPipe = popen("gnuplot -persistent", "w");

	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/basin_of_attraction/fig_basin_size_with_basin_entropy_range_e_multiple_basin_undetermined_monte_carlo_with_break_gamma_%1.6f_system_%s_K_%1.5f_n_%d_rand_%d_window_mc_%d_precision_mc_%1.3f_transient_wn_%d_window_wn_%d_precision_wn_%1.3f.dat", 
		gamma, system.name, K, analysis.number_of_cycles, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc, analysis.convergence_transient_wn, analysis.convergence_window_mc, analysis.convergence_precision_mc);
	// fprintf(gnuplotPipe, 
	// 	"set output \"output/basin_of_attraction/fig_basin_size_with_basin_entropy_range_e_multiple_basin_undetermined_monte_carlo_with_break_synchronous_gamma_%1.6f_system_%s_K_%1.5f_n_%d_rand_%d_window_%d_transient_%d_precision_%1.3f.dat", 
	// 	gamma, system.name, K, analysis.number_of_cycles, analysis.number_of_rand_orbits_mc, analysis.convergence_window, analysis.convergence_transient, analysis.convergence_precision);
	fprintf(gnuplotPipe, "set xlabel \"Orbital eccentricity e\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	// fprintf(gnuplotPipe, 
	// 	"set title \"       Size of multiple basin of attraction (D) for {/Symbol g} = %1.3f K = %1.0e res = %d n = %1.0e {/Symbol e} = %1.0e\"\n", 
	// 	gamma, K, analysis.grid_resolution, (double)analysis.number_of_cycles, analysis.evolve_basin_eps);

	fprintf(gnuplotPipe, "set key opaque box outside top right width 1.1\n");

	fprintf(gnuplotPipe, "unset colorbox\n");
	fprintf(gnuplotPipe, "set cbrange [%d:%d]\n", cb_range_min, cb_range_max);

	fprintf(gnuplotPipe, "plot ");

	// orbit_period = 1;
	// spin_period = 1;
	for (orbit_period = analysis.orbit_period_min; orbit_period <= analysis.orbit_period_max; orbit_period++)
	{
		for (spin_period = analysis.spin_period_min; spin_period <= analysis.spin_period_max; spin_period++)
		{
			fprintf(gnuplotPipe, "'basins_size_undetermined_combined_monte_carlo_with_break.dat' u 1:($2==%d&&$3==%d?$4:1/0) w lp pt %d ps 2 title \"%d/%d\", ", 
				spin_period, orbit_period, orbit_period + 4, spin_period, orbit_period);
		}
	}

	fprintf(gnuplotPipe, "'basins_size_undetermined_combined_monte_carlo_with_break.dat' u 1:($2==0&&$3==0?$4:1/0) w lp pt 7 ps 2 lc rgb \"black\" title \"QP\", ");

	fprintf(gnuplotPipe, "'entropy_size_undetermined_combined_monte_carlo_with_break.dat' u 1:2 w l lw 2 title \"NBE\", ");

	fclose(gnuplotPipe);

	printf("Done!\n");

	return 0;
}

int plot_entropy_comparison_monte_carlo_range_e	(int number_of_e,
												 double e_initial,
												 double e_final,
												 dynsys system,
												 anlsis analysis)
{
	FILE 	*combineFiles;
	FILE 	*gnuplotPipe;
	int		size_filename = 200;
	char 	local_filename[size_filename];
	char	filename[size_filename * number_of_e + 50];
	int		spin_period;
	int		orbit_period;
	int		cb_range_min;
	int		cb_range_max;
	double 	*par = (double *)system.params;
	double 	gamma = par[0];
	double 	K = par[6];
	double 	ec, e_step;

	if (number_of_e < 2)
	{
		printf("Invalid range\n");
		exit(2);
	}

	e_step = (e_final - e_initial) / (double)(number_of_e); 

	sprintf(filename, "paste -s");

	ec = e_initial;
	for(int i = 0; i <= number_of_e; i++)
	{
		sprintf(local_filename, " output/basin_of_attraction/multiple_basin_determined_entropy_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat", 
		gamma, ec, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
		strcat(filename, local_filename);
		ec += e_step;
	}

	strcat(filename, " > output/basin_of_attraction/entropy_size_combined.dat");

	combineFiles = popen(filename, "w");

	fclose(combineFiles);

	sprintf(filename, "paste -s");

	ec = e_initial;
	for(int i = 0; i <= number_of_e; i++)
	{
		sprintf(local_filename, " output/basin_of_attraction/multiple_basin_determined_converged_time_monte_carlo_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.dat", 
		gamma, ec, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc);
		strcat(filename, local_filename);
		ec += e_step;
	}

	strcat(filename, " > output/basin_of_attraction/entropy_converged_time_monte_carlo_combined.dat");

	combineFiles = popen(filename, "w");

	fclose(combineFiles);

	// sprintf(filename, "paste -s");

	// ec = e_initial;
	// for(int i = 0; i <= number_of_e; i++)
	// {
	// 	sprintf(local_filename, " output/basin_of_attraction/multiple_basin_determined_entropy_monte_carlo_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.dat", 
	// 	gamma, ec, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window, analysis.convergence_precision);
	// 	strcat(filename, local_filename);
	// 	ec += e_step;
	// }

	// strcat(filename, " > output/basin_of_attraction/entropy_size_monte_carlo_combined.dat");

	// combineFiles = popen(filename, "w");

	// fclose(combineFiles);

	// sprintf(filename, "paste -s");

	// ec = e_initial;
	// for(int i = 0; i <= number_of_e; i++)
	// {
	// 	sprintf(local_filename, " output/basin_of_attraction/multiple_basin_determined_entropy_converged_monte_carlo_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.dat", 
	// 	gamma, ec, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window, analysis.convergence_precision);
	// 	strcat(filename, local_filename);
	// 	ec += e_step;
	// }

	// strcat(filename, " > output/basin_of_attraction/entropy_size_converged_monte_carlo_combined.dat");

	combineFiles = popen(filename, "w");

	fclose(combineFiles);

	sprintf(filename, "paste -d \"\n\"");

	ec = e_initial;
	for(int i = 0; i <= number_of_e; i++)
	{
		sprintf(local_filename, " output/basin_of_attraction/multiple_basin_determined_size_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat", 
		gamma, ec, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
		strcat(filename, local_filename);
		ec += e_step;
	}

	strcat(filename, " > output/basin_of_attraction/basins_size_combined.dat");

	combineFiles = popen(filename, "w");

	fclose(combineFiles);

	// gnuplotPipe = popen("gnuplot -persistent", "w");

	// fprintf(gnuplotPipe, "reset\n");
	// fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	// fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	// fprintf(gnuplotPipe, 
	// 	"set output \"output/basin_of_attraction/fig_entropy_comparison_monte_carlo_range_e_gamma_%1.6f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.png\"\n", 
	// 	gamma, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window, analysis.convergence_precision);
	// fprintf(gnuplotPipe, "set xlabel \"Orbital eccentricity e\"\n");
	// fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	// // fprintf(gnuplotPipe, 
	// // 	"set title \"       Size of multiple basin of attraction (D) for {/Symbol g} = %1.3f K = %1.0e res = %d n = %1.0e {/Symbol e} = %1.0e\"\n", 
	// // 	gamma, K, analysis.grid_resolution, (double)analysis.number_of_cycles, analysis.evolve_basin_eps);

	// fprintf(gnuplotPipe, "set key opaque box top right width 1.1\n");

	// fprintf(gnuplotPipe, "plot ");

	// fprintf(gnuplotPipe, "'entropy_size_combined.dat' u 1:($2!=0?$2/$3:$2) w l lw 2 title \"NBE - Grid\", ");
	// fprintf(gnuplotPipe, "'entropy_size_monte_carlo_combined.dat' u 1:2 w l lw 2 title \"NBE - MC\", ");

	// fclose(gnuplotPipe);

	// gnuplotPipe = popen("gnuplot -persistent", "w");

	// fprintf(gnuplotPipe, "reset\n");
	// fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	// fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	// fprintf(gnuplotPipe, 
	// 	"set output \"output/basin_of_attraction/fig_entropy_comparison_converged_monte_carlo_range_e_gamma_%1.6f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.png\"\n", 
	// 	gamma, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window, analysis.convergence_precision);
	// fprintf(gnuplotPipe, "set xlabel \"Orbital eccentricity e\"\n");
	// fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	// // fprintf(gnuplotPipe, 
	// // 	"set title \"       Size of multiple basin of attraction (D) for {/Symbol g} = %1.3f K = %1.0e res = %d n = %1.0e {/Symbol e} = %1.0e\"\n", 
	// // 	gamma, K, analysis.grid_resolution, (double)analysis.number_of_cycles, analysis.evolve_basin_eps);

	// fprintf(gnuplotPipe, "set key opaque box top right width 1.1\n");

	// fprintf(gnuplotPipe, "plot ");

	// fprintf(gnuplotPipe, "'entropy_size_combined.dat' u 1:($2!=0?$2/$3:$2) w l lw 2 title \"NBE - Grid\", ");
	// fprintf(gnuplotPipe, "'entropy_size_converged_monte_carlo_combined.dat' u 1:2 w l lw 2 title \"NBE - MC\", ");

	// fclose(gnuplotPipe);

	// gnuplotPipe = popen("gnuplot -persistent", "w");

	// fprintf(gnuplotPipe, "reset\n");
	// fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	// fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	// fprintf(gnuplotPipe, 
	// 	"set output \"output/basin_of_attraction/fig_entropy_convergence_converged_monte_carlo_range_e_gamma_%1.6f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.png\"\n", 
	// 	gamma, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window, analysis.convergence_precision);
	// fprintf(gnuplotPipe, "set xlabel \"Orbital eccentricity e\"\n");
	// fprintf(gnuplotPipe, "set ylabel \"Number of realizations\"\n");
	// fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	// // fprintf(gnuplotPipe, 
	// // 	"set title \"       Size of multiple basin of attraction (D) for {/Symbol g} = %1.3f K = %1.0e res = %d n = %1.0e {/Symbol e} = %1.0e\"\n", 
	// // 	gamma, K, analysis.grid_resolution, (double)analysis.number_of_cycles, analysis.evolve_basin_eps);

	// fprintf(gnuplotPipe, "unset key\n");

	// fprintf(gnuplotPipe, "plot ");

	// fprintf(gnuplotPipe, "'entropy_size_converged_monte_carlo_combined.dat' u 1:3 w l lw 2 notitle");

	// fclose(gnuplotPipe);

	// gnuplotPipe = popen("gnuplot -persistent", "w");

	// fprintf(gnuplotPipe, "reset\n");
	// fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	// fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	// fprintf(gnuplotPipe, 
	// 	"set output \"output/basin_of_attraction/fig_size_multiple_basin_of_attraction_determined_with_entropy_convergence_monte_carlo_range_e_gamma_%1.6f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.png\"\n", 
	// 	gamma, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	// fprintf(gnuplotPipe, "set xlabel \"Orbital eccentricity e\"\n");
	// fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	// // fprintf(gnuplotPipe, 
	// // 	"set title \"       Size of multiple basin of attraction (D) for {/Symbol g} = %1.3f K = %1.0e res = %d n = %1.0e {/Symbol e} = %1.0e\"\n", 
	// // 	gamma, K, analysis.grid_resolution, (double)analysis.number_of_cycles, analysis.evolve_basin_eps);

	// fprintf(gnuplotPipe, "set key opaque box outside top right width 1.1\n");

	// fprintf(gnuplotPipe, "unset colorbox\n");
	// fprintf(gnuplotPipe, "set cbrange [%d:%d]\n", cb_range_min, cb_range_max);

	// fprintf(gnuplotPipe, "plot ");

	// orbit_period = 1;
	// spin_period = 1;
	// for (orbit_period = analysis.orbit_period_min; orbit_period <= analysis.orbit_period_max; orbit_period++)
	// {
	// 	for (spin_period = analysis.spin_period_min; spin_period <= analysis.spin_period_max; spin_period++)
	// 	{
	// 		fprintf(gnuplotPipe, "'basins_size_combined.dat' u 1:($2==%d&&$3==%d?$4:1/0) w lp pt %d ps 2 title \"%d/%d\", ", 
	// 			spin_period, orbit_period, orbit_period + 4, spin_period, orbit_period);
	// 	}
	// }

	// // fprintf(gnuplotPipe, "'basins_size_combined.dat' u 1:(strcol(2) eq \"s\"?$4:1/0) w lp pt 7 ps 2 lc rgb \"black\" title \"sum\", ");

	// fprintf(gnuplotPipe, "'entropy_size_converged_monte_carlo_combined.dat' u 1:($4/8000) w l lw 2 title \"\\# of orbits\"");

	// fclose(gnuplotPipe);

	// gnuplotPipe = popen("gnuplot -persistent", "w");

	// fprintf(gnuplotPipe, "reset\n");
	// fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	// fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	// fprintf(gnuplotPipe, 
	// 	"set output \"output/basin_of_attraction/fig_entropy_comparison_converged_with_convergence_monte_carlo_range_e_gamma_%1.6f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.png\"\n", 
	// 	gamma, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window, analysis.convergence_precision);
	// fprintf(gnuplotPipe, "set xlabel \"Orbital eccentricity e\"\n");
	// fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	// // fprintf(gnuplotPipe, 
	// // 	"set title \"       Size of multiple basin of attraction (D) for {/Symbol g} = %1.3f K = %1.0e res = %d n = %1.0e {/Symbol e} = %1.0e\"\n", 
	// // 	gamma, K, analysis.grid_resolution, (double)analysis.number_of_cycles, analysis.evolve_basin_eps);

	// fprintf(gnuplotPipe, "set key opaque box top right width 1.1\n");

	// fprintf(gnuplotPipe, "plot ");

	// fprintf(gnuplotPipe, "'entropy_size_converged_monte_carlo_combined.dat' u 1:2 w l lw 2 title \"NBE - MC\", ");
	// fprintf(gnuplotPipe, "'entropy_size_converged_monte_carlo_combined.dat' u 1:($4/8000) w l lw 2 title \"\\# of orbits\"");

	// fclose(gnuplotPipe);

	gnuplotPipe = popen("gnuplot -persistent", "w");

	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/basin_of_attraction/fig_entropy_comparison_converged_with_convergence_monte_carlo_range_e_gamma_%1.6f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.png\"\n", 
		gamma, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc);
	fprintf(gnuplotPipe, "set xlabel \"Orbital eccentricity e\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	// fprintf(gnuplotPipe, 
	// 	"set title \"       Size of multiple basin of attraction (D) for {/Symbol g} = %1.3f K = %1.0e res = %d n = %1.0e {/Symbol e} = %1.0e\"\n", 
	// 	gamma, K, analysis.grid_resolution, (double)analysis.number_of_cycles, analysis.evolve_basin_eps);

	fprintf(gnuplotPipe, "set key opaque box top right width 1.1\n");

	fprintf(gnuplotPipe, "plot ");

	fprintf(gnuplotPipe, "'entropy_size_combined.dat' u 1:2 w l lw 2 title \"NBE\", ");
	fprintf(gnuplotPipe, "'entropy_converged_time_monte_carlo_combined.dat' u 1:($3/10000) w l lw 2 title \"\\# of orbits (x10000)\"");

	fclose(gnuplotPipe);

	printf("Done!\n");

	return 0;
}

int plot_comparison_entropy_grid_vs_monte_carlo	(dynsys system,
												 anlsis analysis)
{
	FILE 	*gnuplotPipe;
	double 	*par = (double *)system.params;
	double 	gamma = par[0];
	double	e = par[1];
	double 	K = par[6];

	gnuplotPipe = popen("gnuplot -persistent", "w");

	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/basin_of_attraction/fig_comparison_entropy_grid_vs_monte_carlo_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f_rand_%d.png\"\n", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc);
	fprintf(gnuplotPipe, "set xlabel \"Number of orbits\"\n");
	fprintf(gnuplotPipe, "set ylabel \"Basin entropy\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	// fprintf(gnuplotPipe, 
	// 	"set title \"Basin entropy comparison for {/Symbol g} = %1.3f K = %1.0e res = %d n = %1.0e {/Symbol e} = %1.0e\"\n", 
	// 	gamma, K, analysis.grid_resolution, (double)analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "set title \"e = %1.3f\"\n", e);

	fprintf(gnuplotPipe, "set key box opaque top right\n");

	fprintf(gnuplotPipe, "plot 'comparison_entropy_grid_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat' u 1:2 w lp pt 7 ps 1.5 lw 1.5 title \"Grid\",", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, " 'comparison_entropy_monte_carlo_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d.dat'u 1:2 w lp pt 7 ps 1.5 lw 1.5 title \"Monte Carlo\",", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc);

	fclose(gnuplotPipe);

	gnuplotPipe = popen("gnuplot -persistent", "w");

	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/basin_of_attraction/fig_comparison_entropy_error_grid_vs_monte_carlo_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f_rand_%d.png\"\n", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc);
	fprintf(gnuplotPipe, "set xlabel \"Number of orbits\"\n");
	fprintf(gnuplotPipe, "set ylabel \"Basin entropy error\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	// fprintf(gnuplotPipe, 
	// 	"set title \"Basin entropy comparison for {/Symbol g} = %1.3f K = %1.0e res = %d n = %1.0e {/Symbol e} = %1.0e\"\n", 
	// 	gamma, K, analysis.grid_resolution, (double)analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "set title \"e = %1.3f\"\n", e);

	fprintf(gnuplotPipe, "set key box opaque top right\n");

	// fprintf(gnuplotPipe, "set log y\n");
	fprintf(gnuplotPipe, "set yrange [0.0:0.01]\n");

	fprintf(gnuplotPipe, "plot 'comparison_entropy_grid_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat' u 1:3 w lp pt 7 ps 1.5 lw 1.5 title \"Grid\",", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, " 'comparison_entropy_monte_carlo_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d.dat'u 1:3 w lp pt 7 ps 1.5 lw 1.5 title \"Monte Carlo\",", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc);

	fclose(gnuplotPipe);

	printf("Done!\n");

	return 0;
}

int plot_histogram_python	(dynsys system,
							 anlsis analysis)
{
	FILE 	*pythonPipe;
	int		size_filename = 400;
	char 	filename[size_filename];
	char 	filename_input[size_filename];
	char 	filename_output[size_filename];
	char 	filename_parameter[size_filename];
	double 	*par = (double *)system.params;
	double 	gamma = par[0];
	double	e = par[1];
	double 	K = par[6];

	printf("Plotting histogram for e = %1.3f and gamma = %1.3f\n", e, gamma);

	sprintf(filename, "python3 python_tools/plot_histogram.py");
	sprintf(filename_input, " --input output/basin_of_attraction/multiple_basin_determined_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.dat", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	strcat(filename, filename_input);
	sprintf(filename_output, " --output output/basin_of_attraction/fig_histogram_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f.png", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	strcat(filename, filename_output);
	sprintf(filename_parameter, " --parameter e=%1.3f", e);
	strcat(filename, filename_parameter);
	pythonPipe = popen(filename, "w");
	fclose(pythonPipe);

	printf("Done!\n");

	return 0;
}

int plot_histogram_python_monte_carlo_with_break(dynsys system,
							 					 anlsis analysis)
{
	FILE 	*pythonPipe;
	int		size_filename = 400;
	char 	filename[size_filename];
	char 	filename_input[size_filename];
	char 	filename_output[size_filename];
	char 	filename_parameter[size_filename];
	char 	filename_method[size_filename];
	double 	*par = (double *)system.params;
	double 	gamma = par[0];
	double	e = par[1];
	double 	K = par[6];

	printf("Plotting histogram for e = %1.3f and gamma = %1.3f\n", e, gamma);

	sprintf(filename, "python3 python_tools/plot_histogram.py");
	sprintf(filename_input, " --input output/basin_of_attraction/multiple_basin_determined_times_monte_carlo_with_break_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.dat", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc);
	strcat(filename, filename_input);
	sprintf(filename_output, " --output output/basin_of_attraction/fig_histogram_monte_carlo_with_break_gamma_%1.6f_e_%1.3f_system_%s_K_%1.5f_n_%d_basin_eps_%1.3f_rand_%d_window_%d_precision_%1.3f.png", 
		gamma, e, system.name, K, analysis.number_of_cycles, analysis.evolve_basin_eps, analysis.number_of_rand_orbits_mc, analysis.convergence_window_mc, analysis.convergence_precision_mc);
	strcat(filename, filename_output);
	sprintf(filename_parameter, " --parameter e=%1.3f", e);
	strcat(filename, filename_parameter);
	sprintf(filename_method, " --method MC");
	strcat(filename, filename_method);
	pythonPipe = popen(filename, "w");
	fclose(pythonPipe);

	printf("Done!\n");

	return 0;
}
