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

int trace_orbit_map(double *ic, dynsys system,
					anlsis analysis)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output", &st) == -1) {
		mkdir("output", 0700);
	}

	// declare variables
	FILE *out_orb, *out_orb_ic, *out_orb_res,
		 *out_orb_res_mean,
		 *out_orb_ang_mom_err, *out_vis_viva_err;
	int orbit_size;
	double res_mean;
	double **orbit;
	double orb[4], orb_ini[4];

	for (int i = 0; i < 4; i++)
	{
		orb_ini[i] = ic[i+2];
	}

	// open exit files
	out_orb = fopen("output/orbit.dat", "w");
	out_orb_ic = fopen("output/orbit_ic.dat", "w");
	out_orb_res = fopen("output/orbit_resonance.dat", "w");
	out_orb_res_mean = 
		fopen("output/orbit_resonance_mean.dat", "w");
	out_orb_ang_mom_err = 
		fopen("output/orbit_orbital_angular_momentum_error.dat", "w");
	out_vis_viva_err = fopen("output/orbit_vis_viva_error.dat", "w");

	// evolve system
	evolve_orbit(ic, analysis.cycle_period, 
				analysis.number_of_cycles, 
				&orbit, &orbit_size, system);

	// write orbit and constant error to file
	fprintf(out_orb_ic, "%1.15e %1.15e\n", 
			angle_mod_pos(orbit[0][0]), orbit[0][1]);

	res_mean = 0.0;
	for (int i = 0; i < orbit_size; i++)
	{
		fprintf(out_orb, "%1.15e %1.15e\n", 
				angle_mod_pos(orbit[i][0]), orbit[i][1]);
		
		if (i > 0) 	
		{
			double angle_dif 
					= orbit[i][0] - orbit[i-1][0];
			fprintf(out_orb_res, "%d %1.5e\n", 
				i, angle_dif / analysis.cycle_period);
			res_mean += angle_dif / analysis.cycle_period;
		}

		// if (strcmp(system.name, "rigid") == 0)
		if (system.dim == 6)
		{
			for (int j = 0; j < 4; j++)
			{
				orb[j] = orbit[i][j+2];
			}
			fprintf(out_orb_ang_mom_err, "%d %1.15e\n", 
					i, fabs(angular_momentum_two_body(orb)-
					angular_momentum_two_body(orb_ini)));
			fprintf(out_vis_viva_err, "%d %1.15e\n", 
					i, fabs(vis_viva_two_body(orb)-
					vis_viva_two_body(orb_ini)));
		}
	}

	printf("w = %1.10e\n", 
		angular_dist(orbit[orbit_size-1][0], orbit[0][0])
		/ (double)analysis.number_of_cycles);

	printf("resonance mean = %1.3f\n", 
			res_mean / (double)(orbit_size - 1));

	// free memory
	dealloc_2d_double(&orbit, analysis.number_of_cycles);

	// close files
	fclose(out_orb);
	fclose(out_orb_ic);
	fclose(out_orb_res);
	fclose(out_orb_res_mean);
	fclose(out_orb_ang_mom_err);
	fclose(out_vis_viva_err);

	printf("Data written in output folder\n");

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
	if (stat("output", &st) == -1) {
		mkdir("output", 0700);
	}

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];

	// prepare and open exit files 
	FILE 	*out_psp, *out_ic,
			*out_orb_ang_mom_err, *out_vis_viva_err;
	char	filename[100];

	sprintf(filename, 
		"output/phase_space_gamma_%1.3f_e_%1.3f.dat", gamma, e);
	out_psp = fopen(filename, "w");

	sprintf(filename, 
		"output/phase_space_initial_conditions_gamma_%1.3f_e_%1.3f.dat", 
		gamma, e);
	out_ic = fopen(filename, "w");

	sprintf(filename, 
		"output/phase_space_orbital_angular_momentum_error_gamma_%1.3f_e_%1.3f.dat",
		gamma, e);
	out_orb_ang_mom_err = fopen(filename, "w");

	sprintf(filename, 
		"output/phase_space_vis_viva_error_gamma_%1.3f_e_%1.3f.dat",
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
		printf("Calculating set %d of %ld\n", i + 1, analysis.nc);

		#pragma omp parallel private(y, y0, coordinate, velocity, \
				orbit_fw_size, orbit_bw_size, orbit_fw, orbit_bw)
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
				printf("Calculating subset %d of %ld\n", 
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

				if (strcmp(system.name, "rigid") == 0)
				{
					for (int k = 0; k < 4; k++)
					{
						y[k+2] = orb_ini[k];
					}
				}

				// keep IC for backward integration
				copy(y0, y, system.dim);

				// calculate forward integration
				evolve_orbit(y, analysis.cycle_period, 
					analysis.number_of_cycles, &orbit_fw, 
					&orbit_fw_size, system);

				// calculate backward integration
				evolve_orbit(y0, -1.0 * analysis.cycle_period, 
					analysis.number_of_cycles, &orbit_bw, 
					&orbit_bw_size, system);

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

						if (strcmp(system.name, "rigid") == 0)
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

						if (strcmp(system.name, "rigid") == 0)
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
	fprintf(gnuplotPipe, "set loadpath \"output\"\n");
	fprintf(gnuplotPipe, 
		"set output \"output/fig_phase_space_gamma_%1.3f_e_%1.3f.png\"\n", gamma, e);
	fprintf(gnuplotPipe, "set xlabel \"{/Symbol q}\"\n");
	fprintf(gnuplotPipe, "set ylabel \"~{/Symbol q}{1.1.}\"\n");
	fprintf(gnuplotPipe, "set xrange [0.0:6.29]\n");
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

int mean_resonance(dynsys system, anlsis analysis)
{
	// little warning
	if (strcmp(system.name, "two_body") == 0)
	{
		printf("Warning: cant draw resonance\n");
		printf("for two-body system\n");
		exit(2);
	}

	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output", &st) == -1) {
		mkdir("output", 0700);
	}

	double *par = (double *)system.params;
	double gamma = par[0];
	double e = par[1];

	// prepare and open exit files 
	FILE 	*out_orb_mean_res;
	char	filename[100];

	sprintf(filename, 
		"output/mean_resonance_gamma_%1.3f_e_%1.3f.dat", 
		gamma, e);
	out_orb_mean_res = fopen(filename, "w");

	// declare variables
	double y[system.dim];
	double coordinate, velocity;
	double mean_res;
	int orbit_fw_size;
	double **orbit_fw;
	double orb[4], orb_ini[4];
	init_orbital(orb_ini, e);
	int grid[2];
	double resonance[2];
	double **mean_res_matrix;
	alloc_2d_double(&mean_res_matrix, 
		analysis.grid_resolution, analysis.grid_resolution);

	// loop over coordinate values
	for (int i = 0; i < analysis.grid_resolution; i++)
	{
		// print progress on coordinate
		printf("Calculating set %d of %d\n", 
			i + 1, analysis.grid_resolution);

		omp_set_dynamic(0);     // Explicitly disable dynamic teams
		omp_set_num_threads(12); // Use 4 threads for all consecutive parallel regions

		#pragma omp parallel private(y, coordinate, velocity, \
				mean_res, orbit_fw_size, orbit_fw)
		{

		coordinate = analysis.grid_coordinate_min + 
				(double)(i) * (analysis.grid_coordinate_max - 
				analysis.grid_coordinate_min) / 
				(double)(analysis.grid_resolution - 1);

		#pragma omp for
			// loop over velocity values
			for (int j = 0; j < analysis.grid_resolution; j++)
			{
				// print progress on velocity
				printf("Calculating subset %d of %d\n", 
							j + 1, analysis.grid_resolution);

				velocity = analysis.grid_velocity_min + 
						(double)(j) * (analysis.grid_velocity_max - 
						analysis.grid_velocity_min) / 
						(double)(analysis.grid_resolution - 1);
				
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
				evolve_orbit(y, analysis.cycle_period, 
					analysis.number_of_cycles, &orbit_fw, 
					&orbit_fw_size, system);

				mean_res = 0.0;
				for (int k = 1; k < orbit_fw_size; k++)
				{
					double angle_dif 
							= orbit_fw[k][0] - orbit_fw[k-1][0];
					mean_res += angle_dif / analysis.cycle_period;
				}

				mean_res_matrix[i][j] = 
					mean_res / (double)(orbit_fw_size - 1);

				// free memory
				dealloc_2d_double(&orbit_fw, 
						analysis.number_of_cycles);
			
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
			grid_to_double(grid, resonance, analysis);

			fprintf(out_orb_mean_res, "%1.5f %1.5f %1.10e\n", 
				resonance[0], resonance[1],
				mean_res_matrix[grid[0]][grid[1]]);
		}
		fprintf(out_orb_mean_res, "\n");
	}

	// close exit files
	fclose(out_orb_mean_res);

	dealloc_2d_double(&mean_res_matrix, 
		analysis.grid_resolution);

	printf("Data written in output folder\n");

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

int evolve_basin(double *ic, double *ref, bool *converged,
                 double ***orbit, int *orbit_size,
                 dynsys system, anlsis analysis)
{
	// declare variables
	double y[system.dim], rot[2];
	double box = 1e6;
	double t = 0.0;
	// tolerance for which we say the orbit converged
	double eps = 1e-1;

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
		printf("Error: NULL orbit.");
		exit(2);
	}
	
	// takes into consideration initial condition
	int counter = 1;

	// orbit evolution
	copy(y, ic, system.dim);
	for (int i = 0; i < analysis.number_of_cycles; i++)
	{
		copy(rot, y, 2);
	
		if(dist_from_ref(rot, ref) < eps)
		{
			*converged = true;
			goto out;
		}

		evolve_cycle(y, analysis.cycle_period, &t, system);
	
		// check if orbit diverges
		for (int j = 0; j < system.dim; j++)
		{
			if (fabs(y[j]) > box)
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

int evolve_basin_union (double *ic, double ref[][2],
						int num_of_basins, bool *converged,
                 		double ***orbit, int *orbit_size,
                 		dynsys system, anlsis analysis)
{
	// declare variables
	double y[system.dim];
	double box = 1e6;
	double t = 0.0;
	double dist_from_ref_theta;
	double dist_from_ref_theta_dot;
	double dist_from_ref;
	// tolerance for which we say the orbit converged
	double eps = 1e-1;

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
		printf("Error: NULL orbit.");
		exit(2);
	}
	
	// takes into consideration initial condition
	int counter = 1;

	// orbit evolution
	copy(y, ic, system.dim);
	for (int i = 0; i < analysis.number_of_cycles; i++)
	{
		for (int j = 0; j < num_of_basins; j++)
		{
			dist_from_ref_theta = 
					angular_dist(y[0], ref[j][0]);
			
			dist_from_ref_theta_dot = 
					fabs(y[1]-ref[j][1]);

			dist_from_ref = 
				sqrt(dist_from_ref_theta*dist_from_ref_theta
					+dist_from_ref_theta_dot*dist_from_ref_theta_dot);
		
			if(dist_from_ref < eps)
			{
				*converged = true;
				goto out;
			}
		}

		evolve_cycle(y, analysis.cycle_period, &t, system);
	
		// check if orbit diverges
		for (int j = 0; j < system.dim; j++)
		{
			if (fabs(y[j]) > box)
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

int basin_of_attraction(double *ref, dynsys system, 
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
	if (stat("output", &st) == -1) {
		mkdir("output", 0700);
	}

	// declare variables
	FILE 	*out_boa, *out_ref;
	double y[system.dim];
	double coordinate, velocity;
	int orbit_fw_size;
	double **orbit_fw;
	double *par = (double *)system.params;
	double e = par[1];
	double rot_ini[2];
	double orb[4], orb_ini[4];
	init_orbital(orb_ini, e);
	int basin_counter = 0;
	int grid[2];
	double basin[2];
	int **basin_matrix, **control_matrix;
	double **time_matrix, **rotation_matrix;
	alloc_2d_int(&basin_matrix, 
		analysis.grid_resolution, analysis.grid_resolution);
	alloc_2d_int(&control_matrix, 
		analysis.grid_resolution, analysis.grid_resolution);
	alloc_2d_double(&time_matrix, 
		analysis.grid_resolution, analysis.grid_resolution);
	alloc_2d_double(&rotation_matrix, 
		analysis.grid_resolution, analysis.grid_resolution);
	for (int i = 0; i < analysis.grid_resolution; i++)
	{
		for (int j = 0; j < analysis.grid_resolution; j++)
		{
			basin_matrix[i][j] = 0;
			control_matrix[i][j] = 0;
			time_matrix[i][j] = NAN;
			rotation_matrix[i][j] = NAN;
		}
	}
	bool converged;

	// open exit files
	out_boa = fopen("output/basin.dat", "w");
	out_ref = fopen("output/basin_ref.dat", "w");

	fprintf(out_ref, "%1.15e %1.15e\n", ref[0], ref[1]);		

	// loop over coordinate values
	for (int i = 0; i < analysis.grid_resolution; i++)
	{
		// print progress on coordinate
		printf("Calculating set %d of %d\n", 
					i + 1, analysis.grid_resolution);

		omp_set_dynamic(0);     // Explicitly disable dynamic teams
		omp_set_num_threads(12); // Use 4 threads for all consecutive parallel regions

		#pragma omp parallel private(y, coordinate, velocity, basin, grid, \
				orbit_fw_size, orbit_fw, converged) shared(basin_matrix, \
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
					evolve_basin(y, ref, &converged,
						&orbit_fw, &orbit_fw_size,
						system, analysis);

					if(converged == true)
					{
						basin_matrix[i][j] = 1;
						control_matrix[i][j] = 1;
						time_matrix[i][j] = (double)(orbit_fw_size);
						rotation_matrix[i][j] = 
							orbit_fw[orbit_fw_size-1][0] / 2.*M_PI;
						if (analysis.grid_coordinate_min < 0.0)
						{
							rotation_matrix[i][j] += 1.0;
						}
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

			fprintf(out_boa, "%1.5f %1.5f %d %1.5e %1.5e\n", 
				basin[0], basin[1], 
				basin_matrix[grid[0]][grid[1]],
				time_matrix[grid[0]][grid[1]],
				rotation_matrix[grid[0]][grid[1]]);
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

	dealloc_2d_double(&rotation_matrix, 
					analysis.grid_resolution);

	printf("Data written in output folder\n");

	return 0;
}

int union_basin_of_attraction(double ref[][2], 
	int num_of_basins, dynsys system, anlsis analysis)
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
	if (stat("output", &st) == -1) {
		mkdir("output", 0700);
	}

	// declare variables
	FILE 	*out_boa, *out_ref;
	double y[system.dim];
	double coordinate, velocity;
	int orbit_fw_size;
	double **orbit_fw;
	double *par = (double *)system.params;
	double e = par[1];
	double rot_ini[2];
	double orb[4], orb_ini[4];
	init_orbital(orb_ini, e);
	int basin_counter = 0;
	int grid[2];
	double basin[2];
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
	bool converged;

	// open exit files
	out_boa = fopen("output/basin_union.dat", "w");
	out_ref = fopen("output/basin_ref_union.dat", "w");

	for (int i = 0; i < num_of_basins; i++)
	{
		fprintf(out_ref, "%1.15e %1.15e\n\n", 
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
				orbit_fw_size, orbit_fw, converged) shared(basin_matrix, \
				control_matrix)
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
					evolve_basin_union(y, ref, 
						num_of_basins, &converged,
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

			fprintf(out_boa, "%1.5f %1.5f %d %1.10e\n", 
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

	printf("Data written in output folder\n");

	return 0;
}
