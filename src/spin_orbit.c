#include "spin_orbit.h"

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
	printf("linear field not implemented yet\n");
	exit(2);

	return 0;
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

double vis_viva_two_body(double y[4])
{
	double n, mu, r, v, C;

	double a = 1.0;

	double T = 2.0 * M_PI;

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

int orbit_map(double *ic, dynsys system,
					anlsis analysis)
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

	// prepare and open exit files 
	FILE	*out_orb, *out_orb_ic, 
			*out_orb_res, *out_orb_err;
	char	filename[150];

	sprintf(filename, "output/orbit/orbit_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f.dat", 
		gamma, e, system.name, K);
	out_orb = fopen(filename, "w");

	sprintf(filename, "output/orbit/orbit_ic_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f.dat", 
		gamma, e, system.name, K);
	out_orb_ic = fopen(filename, "w");

	sprintf(filename, "output/orbit/orbit_resonance_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f.dat", 
		gamma, e, system.name, K);
	out_orb_res = fopen(filename, "w");

	sprintf(filename, "output/orbit/orbit_orbital_error_angular_momentum_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f.dat", 
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
	evolve_orbit(ic, &orbit, &orbit_size, system, analysis);

	// write orbit and constant error to file
	fprintf(out_orb_ic, "%1.15e %1.15e\n", 
			angle_mod_pos(orbit[0][0]), orbit[0][1]);

	for (int i = 0; i < orbit_size; i++)
	{
		// fprintf(out_orb, "%1.15e %1.15e\n", 
		// 		angle_mod_pos(orbit[i][0]), orbit[i][1]);
		fprintf(out_orb, "%1.15e %1.15e %d\n", 
				angle_mod(orbit[i][0]), orbit[i][1], i);
		
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

	// calculating the resonance
	last_angle_dif = 
		orbit[orbit_size-1][0] - orbit[orbit_size-2][0];

	fprintf(out_orb_res, "Angular difference:\n\n%1.10e\n\n", 
		last_angle_dif);

	fprintf(out_orb_res, "Number of spins:\n\n%1.10e\n\n", 
		last_angle_dif / (2.*M_PI));

	fprintf(out_orb_res, "Fractional approximations:\n\n"); 
	denominator = 10;
	for (int i = 0; i < 5; i++)
	{
		numerator = 
		(int) round((last_angle_dif / (2.*M_PI)) * (double) denominator);

		gcd = greatest_common_divisor (numerator, denominator);
	
		if (gcd != 1)
		{
			fprintf(out_orb_res, "%d / %d\n", 
				numerator / gcd, denominator / gcd);
		}
		else
		{
			fprintf(out_orb_res, "%d / %d (GCD not found)\n", 
				numerator / gcd, denominator / gcd);
		}
		denominator *= 10;
	}

	// free memory
	dealloc_2d_double(&orbit, analysis.number_of_cycles);

	// close files
	fclose(out_orb);
	fclose(out_orb_ic);
	fclose(out_orb_res);
	fclose(out_orb_err);

	printf("Data written in output/orbit/\n");

	return 0;
}

int phase_space(dynsys system, anlsis analysis)
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
		"output/phase_space/phase_space_gamma_%1.3f_e_%1.3f.dat", gamma, e);
	out_psp = fopen(filename, "w");

	sprintf(filename, 
		"output/phase_space/phase_space_initial_conditions_gamma_%1.3f_e_%1.3f.dat", 
		gamma, e);
	out_ic = fopen(filename, "w");

	sprintf(filename, 
		"output/phase_space/phase_space_orbital_angular_momentum_error_gamma_%1.3f_e_%1.3f.dat",
		gamma, e);
	out_orb_ang_mom_err = fopen(filename, "w");

	sprintf(filename, 
		"output/phase_space/phase_space_vis_viva_error_gamma_%1.3f_e_%1.3f.dat",
		gamma, e);
	out_vis_viva_err = fopen(filename, "w");

	// declare variables
	double y[system.dim], y0[system.dim];
	double coordinate, velocity;
	int orbit_fw_size, orbit_bw_size;
	double **orbit_fw, **orbit_bw;
	double orb[4], orb_ini[4];
	init_orbital(orb_ini, e);

	// loop over coordinate values
	for (int i = 0; i < analysis.nc; i++)
	{
		// print progress on coordinate
		printf("Calculating set %d of %d\n", i + 1, analysis.nc);

		omp_set_dynamic(0);     // Explicitly disable dynamic teams
		omp_set_num_threads(12); // Use 4 threads for all consecutive parallel regions

		#pragma omp parallel private(y, y0, coordinate, velocity, \
				orbit_fw_size, orbit_bw_size, orbit_fw, orbit_bw, orb)
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

		#pragma omp for
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

				// keep IC for backward integration
				copy(y0, y, system.dim);

				// calculate forward integration
				evolve_orbit(y, &orbit_fw, 
					&orbit_fw_size, system, analysis);

				// calculate backward integration
				analysis.cycle_period *= -1.0;
				evolve_orbit(y0, &orbit_bw, 
					&orbit_bw_size, system, analysis);
				analysis.cycle_period *= -1.0;

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
									k, fabs(vis_viva_two_body(orb)-
									vis_viva_two_body(orb_ini)));
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
									k, fabs(vis_viva_two_body(orb)-
									vis_viva_two_body(orb_ini)));
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
				anlsis analysis)
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
		init_orbital(orb, e);
		for (int j = 0; j < 4; j++) ic[j+2] = orb[j];
	}

	// evolve system
	evolve_orbit(ic, &orbit, &orbit_size, system, analysis);

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
						anlsis analysis)
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

	#pragma omp parallel private(orbit_size, orbit, ic, orb)
	{

	#pragma omp for
		for (int i = 0; i < 20; i++)
		{
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
				init_orbital(orb, e);
				for (int j = 0; j < 4; j++) ic[j+2] = orb[j];
			}

			// evolve system
			evolve_orbit(ic, &orbit, &orbit_size, system, analysis);

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
										anlsis analysis)
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
	sprintf(filename, "output/time_series/multiple_time_series_delta_theta_dot_e_%1.3f_K_%1.5f_delta_%1.2f.dat", 
				e, K, analysis.time_series_delta);
	out = fopen(filename, "w");

	// declare variables
	int orbit_size;
	double **orbit;
	double ic[system.dim], orb[4];

	#pragma omp parallel private(orbit_size, orbit, ic, orb)
	{

	#pragma omp for
		// for (int i = 0; i < 20; i++)
		// {
		// 	ic[0] = 0.0; ic[1] = 999.0 + 0.1 * (double)(i);
		for (int i = -10; i <= 10; i++) 
		{
			printf("Calculating time series number %d\n", i + 11);

			ic[0] = 0.0;
			ic[1] = 1000.0 + analysis.time_series_delta * (double)(i);

			if (system.dim == 6)
			{
				init_orbital(orb, e);
				for (int j = 0; j < 4; j++) ic[j+2] = orb[j];
			}

			// evolve system
			evolve_orbit(ic, &orbit, &orbit_size, system, analysis);

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
									anlsis analysis)
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

	#pragma omp parallel private(orbit_size, orbit, ic, orb)
	{

	#pragma omp for
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
				init_orbital(orb, e);
				for (int j = 0; j < 4; j++) ic[j+2] = orb[j];
			}

			// evolve system
			evolve_orbit(ic, &orbit, &orbit_size, system, analysis);

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

int evolve_basin(double *ic, double ref[][2], int ref_period, bool *converged,
                 double ***orbit, int *orbit_size,
                 dynsys system, anlsis analysis)
{
	// declare variables
	double y[system.dim], rot[2];
	double t = 0.0;

	// allocate memory and initializes exit data
	if (orbit != NULL)
	{
		// takes into consideration initial condition
		alloc_2d_double(orbit, 
			analysis.number_of_cycles + 1, 
			system.dim);
		copy((*orbit)[0], ic, system.dim);
	}
	else
	{
		printf("Error: NULL orbit.\n");
		exit(2);
	}
	
	// takes into consideration initial condition
	int counter = 1;

	// orbit evolution
	copy(y, ic, system.dim);
	for (int i = 0; i < analysis.number_of_cycles; i++)
	{
		copy(rot, y, 2);

		for (int j = 0; j < ref_period; j++)
		{
			if(dist_from_ref(rot, ref[j]) < analysis.evolve_basin_eps)
			{
				*converged = true;
				goto out;
			}
		}

		evolve_cycle(y, &t, system, analysis);
	
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

		counter++;

		// write orbit element
		copy((*orbit)[i + 1], y, system.dim);

	}

	out:;

	*orbit_size = counter;

	return 0;
}

int basin_of_attraction(double ref[][2], int ref_period, dynsys system, 
						anlsis analysis)
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
	char	filename[200];

	sprintf(filename, "output/basin_of_attraction/basin_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat", 
		gamma, e, system.name, K, ref[0][0], ref[0][1], ref_period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	out_boa = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/basin_ref_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat", 
		gamma, e, system.name, K, ref[0][0], ref[0][1], ref_period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	out_ref = fopen(filename, "w");

	// declare variables
	double y[system.dim];
	double coordinate, velocity;
	int orbit_fw_size;
	double **orbit_fw;
	double rot_ini[2];
	double orb[4], orb_ini[4];
	init_orbital(orb_ini, e);
	int basin_counter = 0;
	int grid[2];
	double basin[2];
	bool converged;
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

	for (int i = 0; i < ref_period; i++)
	{
		fprintf(out_ref, "%1.15e %1.15e\n", 
						ref[i][0], ref[i][1]);
	}	

	// loop over coordinate values
	for (int i = 0; i < analysis.grid_resolution; i++)
	{
		// print progress on coordinate
		printf("Calculating set %d of %d\n", 
					i + 1, analysis.grid_resolution);

		omp_set_dynamic(0);     // Explicitly disable dynamic teams
		omp_set_num_threads(12); // Use 4 threads for all consecutive parallel regions

		#pragma omp parallel private(y, coordinate, velocity, basin, grid, \
				orbit_fw_size, orbit_fw, converged, orb, rot_ini) shared(basin_matrix, \
				control_matrix, time_matrix)
		{

		#pragma omp for
			// loop over velocity values
			for (int j = 0; j < analysis.grid_resolution; j++)
			{
				// print progress on velocity
				printf("Calculating subset %d of %d\n", 
							j + 1, analysis.grid_resolution);

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

					converged = false;

					// calculate forward integration
					evolve_basin(y, ref, ref_period, &converged,
						&orbit_fw, &orbit_fw_size,
						system, analysis);

					if(converged == true)
					{
						basin_matrix[i][j] = 1;
						control_matrix[i][j] = 1;
						time_matrix[i][j] = (double)(orbit_fw_size);
						basin_counter++;
					}
					else
					{
						control_matrix[i][j] = -1;
					}

					// free memory
					dealloc_2d_double(&orbit_fw, 
							analysis.number_of_cycles);
				}
			
			}
		} // end pragma

		// new line on terminal
		printf("\n");
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

	printf("Basin counter = %d\n", basin_counter);

	// close exit files
	fclose(out_boa);
	fclose(out_ref);

	dealloc_2d_int(&basin_matrix, 
					analysis.grid_resolution);

	dealloc_2d_int(&control_matrix, 
					analysis.grid_resolution);

	dealloc_2d_double(&time_matrix, 
					analysis.grid_resolution);

	printf("Data written in output/basin_of_attraction/\n");

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
		"set output \"output/orbit/fig_orbit_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f.png\"\n", 
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
	fprintf(gnuplotPipe, "plot 'orbit_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f.dat' w p pt 7 ps 1.5 palette notitle, 'orbit_ic_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f.dat' w p pt 7 ps 1.5 notitle",
		gamma, e, system.name, K, gamma, e, system.name, K);
	fclose(gnuplotPipe);

	printf("Done!\n");
	printf("Data written in output/orbit/\n");

	return 0;
}

int draw_phase_space(dynsys system)
{
	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];

	printf("Drawing phase space with gamma = %1.3f and e = %1.3f\n", 
		gamma, e);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/phase_space\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/phase_space/fig_phase_space_gamma_%1.3f_e_%1.3f.png\"\n", gamma, e);
	fprintf(gnuplotPipe, "set xlabel \"{/Symbol q}\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "set xrange[-3.1415:3.1415]\n");
	fprintf(gnuplotPipe, "set yrange [0.0:3.0]\n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, 
		"set key title \"gamma = %1.3f e = %1.3f\" box opaque top right width 2\n", 
		gamma, e);
	fprintf(gnuplotPipe, "plot 'phase_space_gamma_%1.3f_e_%1.3f.dat' w d notitle",
		gamma, e);
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
		"set output \"output/orbit/fig_orbit_on_phase_space_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f.png\"\n", 
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
	fprintf(gnuplotPipe, "plot 'phase_space/phase_space_gamma_%1.3f_e_%1.3f.dat' w d lc rgb \"gray40\" notitle ,'orbit/orbit_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f.dat' w p pt 7 ps 1.5 palette notitle, 'orbit/orbit_ic_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f.dat' w p pt 7 ps 1.5 notitle",
		gamma, e, gamma, e, system.name, K, gamma, e, system.name, K);
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
	fprintf(gnuplotPipe, "set xrange [60:200] \n");
	fprintf(gnuplotPipe, "set yrange [0.5:2] \n");
	fprintf(gnuplotPipe, "set log y\n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, 
		"set key title \"e = %1.3f K = %1.5f\" box opaque top right width 2\n", e, K);
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

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/time_series\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/time_series/fig_time_series_union_K_%1.5f_zoom.png\"\n", K);
	fprintf(gnuplotPipe, "set xlabel \"n\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "set title \"K = %1.5f\"\n", K);
	fprintf(gnuplotPipe, "set log y\n");
	fprintf(gnuplotPipe, "set xrange[0:120]\n");
	fprintf(gnuplotPipe, "set yrange[0.8:1000]\n");
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

	printf("Drawing multiple time series with e = %1.3f and K = %1.5f\n", e, K);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/time_series\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/time_series/fig_multiple_time_series_delta_theta_dot_e_%1.3f_K_%1.5f_delta_%1.2f.png\"\n", e, K, analysis.time_series_delta);
	fprintf(gnuplotPipe, "set xlabel \"n\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "set yrange [0.01:1000] \n");
	// fprintf(gnuplotPipe, "set xrange [0.0:50] \n");
	fprintf(gnuplotPipe, "set log y\n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, 
		"set key title \"e = %1.3f K = %1.5f delta = %1.2f\" box opaque top right width 2\n", e, K, analysis.time_series_delta);
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

int draw_basin_of_attraction(double ref[][2], int ref_period,
                            dynsys system, anlsis analysis)
{
	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	printf("Drawing basin of attraction of system %s with gamma = %1.3f, e = %1.3f and K = %1.5f and reference theta = %1.3f theta_dot = %1.3f with period = %d\n", 
		system.name, gamma, e, K, ref[0][0], ref[0][1], ref_period);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/basin_of_attraction/fig_basin_of_attraction_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.png\"\n", 
		gamma, e, system.name, K, ref[0][0], ref[0][1], ref_period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "set xlabel \"{/Symbol q}\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "set xrange[-3.1415:3.1415]\n");
	fprintf(gnuplotPipe, "set yrange [0.0:3.0]\n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, 
		"set title \"        Basin of attraction  for  gamma = %1.3f  e = %1.3f  K = %1.0e  res = %d  n = %1.0e  eps = %1.0e\"\n", 
		gamma, e, K, analysis.grid_resolution, (double)analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "plot 'basin_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat' u 1:2:3 w image notitle, 'basin_ref_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat' w p pt 7 ps 1.5 lc rgb \"green\" notitle",
		gamma, e, system.name, K, ref[0][0], ref[0][1], ref_period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps, gamma, e, system.name, K, ref[0][0], ref[0][1], ref_period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fclose(gnuplotPipe);

	printf("Done!\n");
	printf("Data written in output/basin_of_attraction/\n");

	printf("Drawing convergence times of system %s with gamma = %1.3f, e = %1.3f and K = %1.5f and reference theta = %1.3f theta_dot = %1.3f with period = %d\n", 
		system.name, gamma, e, K, ref[0][0], ref[0][1], ref_period);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/basin_of_attraction/fig_convergence_times_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.png\"\n", 
		gamma, e, system.name, K, ref[0][0], ref[0][1], ref_period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "set xlabel \"{/Symbol q}\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "set xrange[-3.1415:3.1415]\n");
	fprintf(gnuplotPipe, "set yrange [0.0:3.0]\n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, 
		"set title \"        Convergence times  for  gamma = %1.3f  e = %1.3f  K = %1.0e  res = %d  n = %1.0e  eps = %1.0e\"\n", 
		gamma, e, K, analysis.grid_resolution, (double)analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "plot 'basin_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat' u 1:2:4 w image notitle, 'basin_ref_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat' w p pt 7 ps 1.5 lc rgb \"green\" notitle",
		gamma, e, system.name, K, ref[0][0], ref[0][1], ref_period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps, gamma, e, system.name, K, ref[0][0], ref[0][1], ref_period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fclose(gnuplotPipe);

	printf("Done!\n");
	printf("Data written in output/basin_of_attraction/\n");

	return 0;
}

