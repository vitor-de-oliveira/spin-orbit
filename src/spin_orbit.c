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
	(void)t;

	double *par = (double *)params;

	double gamma 		= par[0];
	double e 			= par[1];
	double m_primary 	= par[2]; 
	double m_secondary 	= par[3];
	double G 			= par[4];
	double a			= par[5];
	double K			= par[6];

	double total_mass = m_primary + m_secondary;

    double r = sqrt((y[2] * y[2]) + (y[3] * y[3]));
    double f_e = atan2(y[3], y[2]);

	double aux = (-3.0/2.0) * gamma * G * m_primary;
	double r_cube = r * r * r;

	double a_over_r = a / r;
	double aux_2 = pow(a_over_r, 6.0);
	
	double n = 1.0;
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
			*out_orb_err;
	char	filename[150];

	sprintf(filename, "output/orbit/orbit_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f.dat", 
		gamma, e, system.name, K);
	out_orb = fopen(filename, "w");

	sprintf(filename, "output/orbit/orbit_ic_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f.dat", 
		gamma, e, system.name, K);
	out_orb_ic = fopen(filename, "w");

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

	// free memory
	dealloc_2d_double(&orbit, analysis.number_of_cycles);

	// close files
	fclose(out_orb);
	fclose(out_orb_ic);
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
	double y[system.dim];
	double coordinate, velocity;
	int orbit_fw_size, orbit_bw_size;
	double **orbit_fw, **orbit_bw;
	double orb[4], orb_ini[4];
	init_orbital(orb_ini, e);
	anlsis analysis_fw, analysis_bw;

	analysis_fw = copy_anlsis(analysis);
	analysis_bw = copy_anlsis(analysis);
	analysis_bw.cycle_period *= -1.0;

	// loop over coordinate values
	for (int i = 0; i < analysis.nc; i++)
	{
		// print progress on coordinate
		printf("Calculating set %d of %d\n", i + 1, analysis.nc);

		omp_set_dynamic(0);     // Explicitly disable dynamic teams
		omp_set_num_threads(12); // Use 4 threads for all consecutive parallel regions

		#pragma omp parallel private(y, coordinate, velocity, \
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

				// calculate forward integration
				evolve_orbit(y, &orbit_fw, 
					&orbit_fw_size, system, analysis_fw);

				// calculate backward integration
				evolve_orbit(y, &orbit_bw, 
					&orbit_bw_size, system, analysis_bw);

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
	char	filename[200];

	// indexed file
	sprintf(filename, "output/time_series/multiple_time_series_delta_theta_dot_e_%1.3f_K_%1.5f_nic_%d_delta_%1.2f_system_%s.dat", 
				e, K, analysis.number_of_time_series, analysis.time_series_delta, system.name);
	out = fopen(filename, "w");

	// declare variables
	int orbit_size;
	double **orbit;
	double ic[system.dim], orb[4];

	#pragma omp parallel private(orbit_size, orbit, ic, orb)
	{

	#pragma omp for
		for (int i = -analysis.number_of_time_series/2; i <= analysis.number_of_time_series/2; i++) 
		{
			printf("Calculating time series number %d\n", i + 1 + analysis.number_of_time_series/2);

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

int evolve_basin(double *ic,
				 bool *converged,
                 int *convergence_time,
                 perorb po,
                 dynsys system,
				 anlsis analysis)
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
	FILE	*out_boa, *out_ref, *out_size;
	char	filename[200];

	sprintf(filename, "output/basin_of_attraction/basin_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat", 
		gamma, e, system.name, K, po.initial_condition[0], po.initial_condition[1], po.period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	out_boa = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/basin_ref_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat", 
		gamma, e, system.name, K, po.initial_condition[0], po.initial_condition[1], po.period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	out_ref = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/basin_size_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat", 
		gamma, e, system.name, K, po.initial_condition[0], po.initial_condition[1], po.period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	out_size = fopen(filename, "w");

	// declare variables
	double y[system.dim];
	double coordinate, velocity;
	double rot_ini[2];
	double orb[4], orb_ini[4];
	init_orbital(orb_ini, e);
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

		omp_set_dynamic(0);     // Explicitly disable dynamic teams
		omp_set_num_threads(12); // Use 4 threads for all consecutive parallel regions

		#pragma omp parallel private(y, coordinate, velocity, basin, grid, \
				converged, convergence_time, orb, rot_ini) shared(basin_matrix, \
				control_matrix, time_matrix)
		{

		#pragma omp for
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
						po, system, analysis);

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
									 anlsis analysis)
{
	// declare variables
	bool is_close_to;
	int orbit_counter, close_time_counter;
	int internal_converged_po_id;
	double y[system.dim], rot[2];
	double t = 0.0;

	// takes into consideration initial condition
	orbit_counter = 1;

	// starts proximity counter
	close_time_counter = 0;
	
	*converged_po_id = -1;

	copy(y, ic, system.dim);
	for (int i = 0; i < analysis.number_of_cycles; i++)
	{
		copy(rot, y, 2);

		is_close_to = false;
		for (int j = 0; j < number_of_po; j++)
		{
			for (int k = 0; k < po[j].period; k++)
			{
				if(dist_from_ref(rot, po[j].orbit[k]) < analysis.evolve_basin_eps)
				{
					is_close_to = true;
					internal_converged_po_id = j;	// tem um problema aqui
					// as vezes o calculo da orbita periodica encontra a mesma
					// orbita para 1/1 e 2/2. Nessa linha, é atribuído o ponto
					// da bacia de 1/1 para 2/2
				}
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
			*converged_po_id = internal_converged_po_id;
			goto out;
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

		orbit_counter++;

	}

	out:;

	if (*converged_po_id != -1)
	{
		*convergence_time = orbit_counter;
	}
	else
	{
		*convergence_time = -1;
	}

	return 0;
}

int multiple_basin_of_attraction_determined (int number_of_po,
											 perorb po[],
                         					 dynsys system,
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
	FILE	*out_boa, *out_ref, *out_size;
	char	filename[200];

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f_number_of_po_%d.dat", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps, number_of_po);
	out_boa = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_ref_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f_number_of_po_%d.dat", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps, number_of_po);
	out_ref = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/multiple_basin_determined_size_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f_number_of_po_%d.dat", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps, number_of_po);
	out_size = fopen(filename, "w");

	// declare variables
	double y[system.dim];
	double coordinate, velocity;
	double rot_ini[2];
	double orb[4], orb_ini[4];
	init_orbital(orb_ini, e);

	int grid[2];
	double basin[2];
	
	int converged_po_id;
	int convergence_time;

	for (int i = 0; i < number_of_po; i++)
	{
		for (int j = 0; j < po[i].period; j++)
		{
			fprintf(out_ref, "%1.15e %1.15e %d %d\n",
				po[i].orbit[j][0], po[i].orbit[j][1],
				po[i].winding_number_numerator,
				po[i].winding_number_denominator);
		}
		fprintf(out_ref, "\n");
	}

	int *basin_size;
	int **basin_matrix, **control_matrix;
	double **time_matrix;

	alloc_1d_int(&basin_size, number_of_po);
	for (int i = 0; i < number_of_po; i++)
	{
		basin_size[i] = 0;
	}

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

		omp_set_dynamic(0);     // Explicitly disable dynamic teams
		omp_set_num_threads(12); // Use 4 threads for all consecutive parallel regions

		#pragma omp parallel private(y, coordinate, velocity, basin, grid, \
				converged_po_id, convergence_time, orb, rot_ini) shared(basin_matrix, \
				control_matrix, time_matrix)
		{

		#pragma omp for
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
					evolve_multiple_basin_determined(y, number_of_po, 
						&converged_po_id, &convergence_time,
						po, system, analysis);

					if(converged_po_id != -1)
					{
						basin_matrix[i][j] = 
							cantor_pairing_function(po[converged_po_id].winding_number_numerator,po[converged_po_id].winding_number_denominator);
						control_matrix[i][j] = converged_po_id + 1;
						time_matrix[i][j] = (double)(convergence_time);
						basin_size[converged_po_id]++;
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
				res_spin = po[control_matrix[i][j]-1].winding_number_numerator;
				res_orbit = po[control_matrix[i][j]-1].winding_number_denominator;
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
		if (po[i].winding_number_numerator > winding_number_numerator_max)
		{
			winding_number_numerator_max = po[i].winding_number_numerator;
		}
		if (po[i].winding_number_denominator > winding_number_denominator_max)
		{
			winding_number_denominator_max = po[i].winding_number_denominator;
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
				if (po[k].winding_number_numerator == i &&
					po[k].winding_number_denominator == j)
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

	// close exit files
	fclose(out_boa);
	fclose(out_ref);
	fclose(out_size);

	dealloc_1d_int(&basin_size);

	dealloc_2d_int(&basin_matrix, 
					analysis.grid_resolution);

	dealloc_2d_int(&control_matrix, 
					analysis.grid_resolution);

	dealloc_2d_double(&time_matrix, 
					analysis.grid_resolution);

	printf("Data written in output/basin_of_attraction/\n");

	return 0;
}
int evolve_multiple_basin_undetermined  (double *ic,
									     bool *converged,
                                         int *attractor_period,
									     int *convergence_time,
									     dynsys system,
									     anlsis analysis)
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
	FILE	*out_boa, *out_ref, *out_size;
	char	filename[200];

	sprintf(filename, "output/basin_of_attraction/multiple_basin_undetermined_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d.dat", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles);
	out_boa = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/multiple_basin_undetermined_ref_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d.dat", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles);
	out_ref = fopen(filename, "w");

	sprintf(filename, "output/basin_of_attraction/multiple_basin_undetermined_size_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d.dat", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles);
	out_size = fopen(filename, "w");

	// declare variables
	double y[system.dim];
	double rot_ini[2];
	double orb[4], orb_ini[4];
	init_orbital(orb_ini, e);
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

		omp_set_dynamic(0);     // Explicitly disable dynamic teams
		omp_set_num_threads(12); // Use 4 threads for all consecutive parallel regions

		#pragma omp parallel private(y, grid, converged, po_already_found, attractor_period, \
				converged_po_id, convergence_time, orb, rot_ini) shared(basin_matrix, \
				control_matrix, time_matrix, number_of_po, basin_size, multiple_po)
		{

		#pragma omp for
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
						system, analysis);

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
								if (periodic_orbit(&po_try, system, analysis) == 0)
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

int look_for_resonance	(int number_of_candidates,
						 double candidates[][2],
						 int spin_period,
                         int orbit_period,
                         dynsys system, 
						 anlsis analysis)
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
	char	filename[200];

	sprintf(filename, "output/periodic_orbit/resonance_distances_gamma_%1.3f_e_%1.3f_spin_orbit_%d_%d.dat", 
		gamma, e, spin_period, orbit_period);
	out_dist = fopen(filename, "w");
	sprintf(filename, "output/periodic_orbit/resonance_candidates_gamma_%1.3f_e_%1.3f_spin_orbit_%d_%d.dat", 
		gamma, e, spin_period, orbit_period);
	out_cand = fopen(filename, "w");

	// declare variables
	double y[system.dim], y0[system.dim];
	double rot_ini[2], orb_ini[4];
	init_orbital(orb_ini, e);
	int grid[2];
	double ic[2];
	double t;

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

		omp_set_dynamic(0);     // Explicitly disable dynamic teams
		omp_set_num_threads(12); // Use 4 threads for all consecutive parallel regions

		#pragma omp parallel private(y, y0, grid, rot_ini, t) shared(orbit_distance, spin_distance)
		{

		#pragma omp for
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
				for (int i = 0; i < orbit_period; i++)
				{
					evolve_cycle(y, &t, system, analysis);
				}

				orbit_distance[i][j] = dist_from_ref(y,y0);

				spin_distance[i][j] = 
					fabs((y[0] - y0[0]) / (2.*M_PI) - (double) spin_period);
			}
		} // end pragma

		// new line on terminal
		// printf("\n");
	}

	for (int i = 0; i < analysis.grid_resolution; i++)
	{
		for (int j = 0; j < analysis.grid_resolution; j++)
		{
			if (orbit_distance[i][j] != orbit_distance[i][j] ||
				 spin_distance[i][j] != spin_distance[i][j])
			{
				printf("Looking for the resonance gave NAN value.\n");
			}

			grid[0] = i;
			grid[1] = j;
			grid_to_double(grid, ic, analysis);

			fprintf(out_dist, "%1.10f %1.10f %1.10e %1.10e\n", 
				ic[0], ic[1], orbit_distance[i][j],	spin_distance[i][j]);
		}
		fprintf(out_dist, "\n");
	}

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
		fprintf(out_cand, "%1.5f %1.5f\n", ic[0], ic[1]);
		copy(candidates[i], ic, 2);
	}
				
	dealloc_2d_int(&candidates_ij, number_of_candidates);

	dealloc_2d_double(&orbit_distance, 
					analysis.grid_resolution);
	dealloc_2d_double(&spin_distance, 
					analysis.grid_resolution);

	// close exit files
	fclose(out_dist);
	fclose(out_cand);

	printf("Data written in output/periodic_orbit/\n");

	return 0;
}

int evolve_n_cycles_po  (double y0[2],
						 int n,
                         dynsys system,
                         anlsis analysis)
{
    double t = 0;
	double y[system.dim], orbital[4];
	double *par = (double *)system.params;
	double e = par[1];

	copy (y, y0, 2);

	if (system.dim == 6)
	{
		init_orbital(orbital, e);
		for (int i = 0; i < 4; i++)
		{
			y[i+2] = orbital[i];
		}
	}

    // loop over cycles
	for (int i = 0; i < n; i++)
	{
		evolve_cycle(y, &t, system, analysis);
	
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
                     anlsis analysis)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/periodic_orbit", &st) == -1)
    {
		mkdir("output/periodic_orbit", 0700);
	}

	FILE    *out_orb, *out_orb_res;
	char    filename[200];
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

	if (calculate_periodic_orbit_ic(po, system, analysis) == -1)
	{
		return -1;
	}

    copy(y, (*po).initial_condition, 2);
	if (system.dim == 6)
	{
		init_orbital(orbital, e);
		for (int i = 0; i < 4; i++)
		{
			y[i+2] = orbital[i];
		}
	}

    for (int i = 0; i < (*po).period; i++)
    {
        copy((*po).orbit[i], y, 2);
		evolve_cycle(y, &t, system, analysis);
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
	sprintf(filename, "output/periodic_orbit/periodic_orbit_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_period_%d_ic_%1.3f_%1.3f.dat", 
            gamma, e, system.name, K, (*po).period, angle_mod((*po).initial_condition[0]), (*po).initial_condition[1]);
	out_orb = fopen(filename, "w");

	copy(y, (*po).initial_condition, 2);
	if (system.dim == 6)
	{
		init_orbital(orbital, e);
		for (int i = 0; i < 4; i++)
		{
			y[i+2] = orbital[i];
		}
	}

    for (int i = 0; i < (*po).period; i++)
    {
        copy((*po).orbit[i], y, system.dim);

		fprintf(out_orb, "%1.10e %1.10e\n", 
            angle_mod((*po).orbit[i][0]), (*po).orbit[i][1]);

		evolve_cycle(y, &t, system, analysis);
    }

	copy(po_ic_after_one_period, y, system.dim);

	sprintf(filename, "output/periodic_orbit/periodic_orbit_resonance_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_period_%d_ic_%1.3f_%1.3f.dat", 
            gamma, e, system.name, K, (*po).period, angle_mod((*po).initial_condition[0]), (*po).initial_condition[1]);
	out_orb_res = fopen(filename, "w");

	// calculating the resonance
	one_period_angular_diff = 
		po_ic_after_one_period[0] - (*po).initial_condition[0];

	number_of_spins = one_period_angular_diff / (2.*M_PI);

	(*po).winding_number_numerator = (int) round(number_of_spins);
	(*po).winding_number_denominator = (*po).period;

	fprintf(out_orb_res, "Orbit:\n\n");

	for (int i = 0; i < (*po).period; i++)
    {
        fprintf(out_orb_res, "%1.5e %1.5e\n", 
            angle_mod((*po).orbit[i][0]), (*po).orbit[i][1]);
    }

	fprintf(out_orb_res, "\n\n");

	fprintf(out_orb_res, "Orbit period:\n\n%d\n\n\n", 
		(*po).period);

	fprintf(out_orb_res, "Angular difference after 1 period:\n\n%1.10e\n\n\n", 
		one_period_angular_diff);

	fprintf(out_orb_res, "Number of spins:\n\n%1.10e\n\n\n", 
		number_of_spins);

	fprintf(out_orb_res, "Resonance:\n\n%d / %d\n\n\n", 
		(*po).winding_number_numerator, 
		(*po).winding_number_denominator);

    fclose(out_orb);
	fclose(out_orb_res);

	return 0;
}

int linear_average_benchmark()
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output/test", &st) == -1) {
		mkdir("output/test", 0700);
	}

    double gamma;			// equatorial flattening
    double e;				// eccentricity
	double m_primary;		// mass of primary
	double m_secondary;		// mass of secondary
	double G;				// gravitational constant
	double a;				// semimajor axis
	double K;				// dissipation parameter

	double *params[7] = {&gamma,
						 &e,
						 &m_primary,
						 &m_secondary,
						 &G,
						 &a,
						 &K};

	anlsis analysis;

	dynsys system = init_two_body(*params);

	// gamma = ((.89 * .89) / 3.); //~0.264
	m_secondary = 0.;
	m_primary = 1.0 - m_secondary;
	G = 1.0;
	a = 1.0;
	K = 1e-2;

	analysis.number_of_cycles = 1e3;
	analysis.cycle_period = 2.0 * M_PI * 1e-3;
	analysis.evolve_box_size = 1e8;
	analysis.evolve_basin_eps = 1e-1;

	printf("Calculating orbit\n");

	// prepare and open exit files 
	FILE	*out, *out_avg, *out_num_avg, *out_avg_dist;
	char	filename[150];

	e = 0.2;

	sprintf(filename, 
		"output/test/benchmark_linear_e_%1.3f.dat", e);
	out = fopen(filename, "w");
	sprintf(filename, 
		"output/test/benchmark_linear_average_e_%1.3f.dat", e);
	out_avg = fopen(filename, "w");
	sprintf(filename, 
		"output/test/benchmark_linear_numerical_average_e_%1.3f.dat", e);
	out_num_avg = fopen(filename, "w");
	sprintf(filename, 
		"output/test/benchmark_linear_average_distance.dat");
	out_avg_dist = fopen(filename, "w");
	
	// declare variables
	int orbit_size;
	double **orbit;
	double ic[4];
	double e_2, e_4, e_6, L_avg, N_avg;
	double L_num_avg, N_num_avg;
	double r, a_over_r, aux_2, n, f_dot, L, N;
	double distance_L, distance_N;

	e_2 = e * e;
	e_4 = e * e * e * e;
	e_6 = e * e * e * e * e * e;
	L_avg = (1.+3.*e_2+(3./8.)*e_4) / 
			pow(1.-e_2,(9./2.));
	N_avg = (1.+(15./2.)*e_2+(45./8.)*e_4+
			(5./16.)*e_6) / pow(1.-e_2,6.);

	init_orbital(ic, e);

	// evolve system
	evolve_orbit(ic, &orbit, &orbit_size, system, analysis);

	L_num_avg = 0.0;
	N_num_avg = 0.0;
	for (int i = 0; i < orbit_size; i++)
	{
		r = sqrt((orbit[i][2] * orbit[i][2]) + (orbit[i][3] * orbit[i][3]));
		a_over_r = a / r;
		aux_2 = pow(a_over_r, 6.0);
		n = 1.0;
		f_dot = a_over_r * a_over_r * n * sqrt(1.0 - (e * e));
		L = aux_2;
		N = aux_2 * f_dot;
		L_num_avg += L;
		N_num_avg += N;
	}
	L_num_avg /= (double) orbit_size;
	N_num_avg /= (double) orbit_size;

	for (int i = 0; i < orbit_size; i++)
	{
		r = sqrt((orbit[i][2] * orbit[i][2]) + (orbit[i][3] * orbit[i][3]));
		a_over_r = a / r;
		aux_2 = pow(a_over_r, 6.0);
		n = 1.0;
		f_dot = a_over_r * a_over_r * n * sqrt(1.0 - (e * e));
		L = aux_2;
		N = aux_2 * f_dot;
		fprintf(out, "%1.10e %1.10e\n", L, N);
		fprintf(out_avg, "%1.10e %1.10e\n", L_avg, N_avg);
		fprintf(out_num_avg, "%1.10e %1.10e\n", L_num_avg, N_num_avg);
	}

	// free memory
	dealloc_2d_double(&orbit, analysis.number_of_cycles);

	for (e = 0.0; e < 0.201; e += 0.001)
	{
		printf("Calculating distance for e = %1.3e\n", e);
	
		init_orbital(ic, e);
		evolve_orbit(ic, &orbit, &orbit_size, system, analysis);

		L_num_avg = 0.0;
		N_num_avg = 0.0;
		for (int i = 0; i < orbit_size; i++)
		{
			r = sqrt((orbit[i][2] * orbit[i][2]) + (orbit[i][3] * orbit[i][3]));
			a_over_r = a / r;
			aux_2 = pow(a_over_r, 6.0);
			n = 1.0;
			f_dot = a_over_r * a_over_r * n * sqrt(1.0 - (e * e));
			L = aux_2;
			N = aux_2 * f_dot;
			L_num_avg += L;
			N_num_avg += N;
		}
		L_num_avg /= (double) orbit_size;
		N_num_avg /= (double) orbit_size;
		e_2 = e * e;
		e_4 = e * e * e * e;
		e_6 = e * e * e * e * e * e;
		L_avg = (1.+3.*e_2+(3./8.)*e_4) / 
				pow(1.-e_2,(9./2.));
		N_avg = (1.+(15./2.)*e_2+(45./8.)*e_4+
				(5./16.)*e_6) / pow(1.-e_2,6.);
		distance_L = fabs(L_num_avg - L_avg);
		distance_N = fabs(N_num_avg - N_avg);
		fprintf(out_avg_dist, "%1.3e %1.10e %1.10e\n", 
				e, distance_L, distance_N);

		dealloc_2d_double(&orbit, analysis.number_of_cycles);
		fprintf(out_avg_dist, "\n");
	}

	// close files
	fclose(out);
	fclose(out_avg);
	fclose(out_num_avg);
	fclose(out_avg_dist);

	printf("Data written in output/test/\n");

	return 0;

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
		"set key title \"~{/Symbol g}{0.5-} = %1.3f e = %1.3f\" box opaque top right width 2\n", 
		gamma, e);
	fprintf(gnuplotPipe, "plot 'phase_space_gamma_%1.3f_e_%1.3f.dat' w d notitle",
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
		"set output \"output/clean_figures/fig_phase_space_gamma_%1.3f_e_%1.3f.png\"\n", gamma, e);
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
	// fprintf(gnuplotPipe, "plot 'phase_space_gamma_%1.3f_e_%1.3f.dat' w d notitle", gamma, e);
	fprintf(gnuplotPipe, "plot 'phase_space_gamma_%1.3f_e_%1.3f.dat' w p pt 7 ps 0.2 notitle", gamma, e);
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
		"set output \"output/tests/fig_phase_space_gamma_%1.3f_e_%1.3f.png\"\n", gamma, e);
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
	// fprintf(gnuplotPipe, "plot 'phase_space_gamma_%1.3f_e_%1.3f.dat' w d lc rgb \"black\" notitle",
	// 	gamma, e);
	// fprintf(gnuplotPipe, "plot 'phase_space_gamma_%1.3f_e_%1.3f.dat' w d notitle",
	// 	gamma, e);
	fprintf(gnuplotPipe, "plot 'phase_space_gamma_%1.3f_e_%1.3f.dat' w p pt 7 ps 0.2 notitle", gamma, e);
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
		"set output \"output/orbit/fig_orbit_on_phase_space_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_latex.png\"\n", 
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
	fprintf(gnuplotPipe, "plot 'phase_space/phase_space_gamma_%1.3f_e_%1.3f.dat' w d notitle ,'orbit/orbit_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f.dat' w p pt 7 ps 3 palette notitle", 
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
		"set output \"output/basin_of_attraction/fig_basin_of_attraction_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.png\"\n", 
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
	fprintf(gnuplotPipe, "FILE = \"basin_size_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat\"\n", 
		gamma, e, system.name, K, po.initial_condition[0], po.initial_condition[1], po.period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "stats FILE u (myTitle=strcol(1),0) nooutput\n");
	fprintf(gnuplotPipe, "set title \"gamma = %1.3f  e = %1.3f  K = %1.0e  res = %d  n = %1.0e  eps = %1.0e  size = \".myTitle\n", 
		gamma, e, K, analysis.grid_resolution, (double)analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "plot 'basin_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat' u 1:2:3 w image notitle, 'basin_ref_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat' w p pt 7 ps 1.5 lc rgb \"green\" notitle",
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
		"set output \"output/basin_of_attraction/fig_convergence_times_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.png\"\n", 
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
	fprintf(gnuplotPipe, "plot 'basin_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat' u 1:2:4 w image notitle, 'basin_ref_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat' w p pt 7 ps 1.5 lc rgb \"green\" notitle",
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
		"set output \"output/clean_figures/fig_basin_of_attraction_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.png\"\n", 
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
	fprintf(gnuplotPipe, "plot 'basin_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat' u 1:2:3 w image notitle, 'basin_ref_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat' w p pt 7 ps 1.5 lc rgb \"green\" notitle",
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
		"set output \"output/clean_figures/fig_convergence_times_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.png\"\n", 
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
	fprintf(gnuplotPipe, "plot 'basin_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat' u 1:2:4 w image notitle, 'basin_ref_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_ref_%1.3f_%1.3f_period_%d_res_%d_n_%d_basin_eps_%1.3f.dat' w p pt 7 ps 1.5 lc rgb \"green\" notitle",
		gamma, e, system.name, K, ref[0][0], ref[0][1], ref_period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps, gamma, e, system.name, K, ref[0][0], ref[0][1], ref_period, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps);
	fclose(gnuplotPipe);

	printf("Done!\n");
	printf("Data written in output/clean_figures/\n");

	return 0;
}

int draw_multiple_basin_of_attraction_determined(int number_of_po,
                                        		 dynsys system,
                                        		 anlsis analysis)
{
	FILE *gnuplotPipe;

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];
	double K = par[6];

	printf("Drawing basin of attraction of system %s with gamma = %1.3f, e = %1.3f and K = %1.5f for %d resonances\n", 
		system.name, gamma, e, K, number_of_po);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/basin_of_attraction/fig_multiple_basin_of_attraction_determined_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f_number_of_po_%d.png\"\n", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps, number_of_po);
	fprintf(gnuplotPipe, "set xlabel \"{/Symbol q}\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "unset colorbox\n");
	fprintf(gnuplotPipe, "set xrange[-3.1415:3.1415]\n");
	fprintf(gnuplotPipe, "set yrange [0.0:3.0]\n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, 
		"set title \"Multiple basin of attraction (D) for {/Symbol g} = %1.3f e = %1.3f K = %1.0e res = %d n = %1.0e {/Symbol e} = %1.0e\"\n", 
		gamma, e, K, analysis.grid_resolution, (double)analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "plot 'multiple_basin_determined_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f_number_of_po_%d.dat' u 1:2:3 w image notitle",
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps, number_of_po);
	fprintf(gnuplotPipe, ", 'multiple_basin_determined_ref_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f_number_of_po_%d.dat' u 1:2 w p pt 7 ps 2 lc rgb \"green\" title \"spin-orbit resonances\"",
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps, number_of_po);
	fclose(gnuplotPipe);

	printf("Done!\n");
	printf("Data written in output/basin_of_attraction/\n");

	printf("Drawing convergence times for system %s with gamma = %1.3f, e = %1.3f and K = %1.5f for %d resonances\n", 
		system.name, gamma, e, K, number_of_po);

	gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "reset\n");
	fprintf(gnuplotPipe, "set terminal pngcairo size 920,800 font \"Helvetica,15\"\n");
	fprintf(gnuplotPipe, "set loadpath \"output/basin_of_attraction\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/basin_of_attraction/fig_multiple_convergence_time_determined_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f_number_of_po_%d.png\"\n", 
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps, number_of_po);
	fprintf(gnuplotPipe, "set xlabel \"{/Symbol q}\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set ylabel offset 0.8 \n");
	fprintf(gnuplotPipe, "unset colorbox\n");
	fprintf(gnuplotPipe, "set xrange[-3.1415:3.1415]\n");
	fprintf(gnuplotPipe, "set yrange [0.0:3.0]\n");
	fprintf(gnuplotPipe, "unset key\n");
	fprintf(gnuplotPipe, 
		"set title \"Convergence times (D) for {/Symbol g} = %1.3f e = %1.3f K = %1.0e res = %d n = %1.0e {/Symbol e} = %1.0e\"\n", 
		gamma, e, K, analysis.grid_resolution, (double)analysis.number_of_cycles, analysis.evolve_basin_eps);
	fprintf(gnuplotPipe, "plot 'multiple_basin_determined_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f_number_of_po_%d.dat' u 1:2:4 w image notitle",
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps, number_of_po);
	fprintf(gnuplotPipe, ", 'multiple_basin_determined_ref_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d_basin_eps_%1.3f_number_of_po_%d.dat' u 1:2 w p pt 7 ps 2 lc rgb \"green\" title \"spin-orbit resonances\"",
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles, analysis.evolve_basin_eps, number_of_po);
	fclose(gnuplotPipe);

	printf("Done!\n");
	printf("Data written in output/basin_of_attraction/\n");

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
		"set output \"output/basin_of_attraction/fig_multiple_basin_of_attraction_undetermined_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d.png\"\n", 
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
	fprintf(gnuplotPipe, "plot 'multiple_basin_undetermined_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d.dat' u 1:2:3 w image notitle",
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles);
	fprintf(gnuplotPipe, ", 'multiple_basin_undetermined_ref_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d.dat' u 1:2 w p pt 7 ps 2 lc rgb \"green\" title \"spin-orbit resonances\"",
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
		"set output \"output/basin_of_attraction/fig_multiple_convergence_time_undetermined_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d.png\"\n", 
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
	fprintf(gnuplotPipe, "plot 'multiple_basin_undetermined_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d.dat' u 1:2:4 w image notitle",
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles);
	fprintf(gnuplotPipe, ", 'multiple_basin_undetermined_ref_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_res_%d_n_%d.dat' u 1:2 w p pt 7 ps 2 lc rgb \"green\" title \"spin-orbit resonances\"",
		gamma, e, system.name, K, analysis.grid_resolution, analysis.number_of_cycles);
	fclose(gnuplotPipe);

	printf("Done!\n");
	printf("Data written in output/basin_of_attraction/\n");

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
		"set output \"output/periodic_orbit/fig_periodic_orbit_on_phase_space_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_period_%d_ic_%1.3f_%1.3f.png\"\n", 
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
	fprintf(gnuplotPipe, "plot 'phase_space/phase_space_gamma_%1.3f_e_%1.3f.dat' w d lc rgb \"gray40\" notitle ,'periodic_orbit/periodic_orbit_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_period_%d_ic_%1.3f_%1.3f.dat' w p pt 7 ps 1.5 lc rgb \"black\" notitle",
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
		"set output \"output/clean_figures/fig_periodic_orbit_on_phase_space_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_period_%d_ic_%1.3f_%1.3f.png\"\n", 
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
	// fprintf(gnuplotPipe, "plot 'phase_space/phase_space_gamma_%1.3f_e_%1.3f.dat' w d notitle ,'periodic_orbit/periodic_orbit_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_period_%d_ic_%1.3f_%1.3f.dat' w p pt 7 ps 4 lc rgb \"black\" notitle",
	// 	gamma, e, gamma, e, system.name, K, po.period, po.initial_condition[0], po.initial_condition[1]);
	fprintf(gnuplotPipe, "plot 'phase_space/phase_space_gamma_%1.3f_e_%1.3f.dat' w p pt 7 ps 0.2 notitle ,'periodic_orbit/periodic_orbit_gamma_%1.3f_e_%1.3f_system_%s_K_%1.5f_period_%d_ic_%1.3f_%1.3f.dat' w p pt 7 ps 5 lc rgb \"black\" notitle",
		gamma, e, gamma, e, system.name, K, po.period, po.initial_condition[0], po.initial_condition[1]);
	fclose(gnuplotPipe);

	printf("Done!\n");
	printf("Data written in output/periodic_orbit/\n");

	return 0;
}