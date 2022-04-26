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

	// declare and open exit files
	FILE *orb = fopen("output/orbit.dat", "w");

	// declare variables
	int orbit_size;
	double **orbit;

	// evolve system
	evolve_orbit(ic, analysis.cycle_period, 
				analysis.number_of_cycles, 
				&orbit, &orbit_size, system);

	// write orbit and constant error to file
	for (int i = 0; i < orbit_size; i++)
	{
		// fprintf(orb, "%1.15e %1.15e\n", 
		// 		angle_mod(orbit[i][0]), orbit[i][1]);
		fprintf(orb, "%1.15e %1.15e\n", 
			fmod(orbit[i][0], 2.0*M_PI), orbit[i][1]);
	}

	printf("w = %1.10e\n", 
		angular_dist(orbit[orbit_size-1][0], orbit[0][0])
		/ (double)analysis.number_of_cycles);

	// free memory
	dealloc_2d_double(&orbit, analysis.number_of_cycles);

	// close files
	fclose(orb);

	printf("Data written in output folder\n");

	return 0;
}

int draw_phase_space(dynsys system, anlsis analysis)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output", &st) == -1) {
		mkdir("output", 0700);
	}

	// declare and open exit files
	FILE *psp = fopen("output/phase_space.dat", "w");
	FILE *inc 
		= fopen("output/phase_space_initial_conditions.dat", 
				"w");

	if (strcmp(system.name, "two_body") == 0)
	{
		printf("Warning: cant draw phase space\n");
		printf("for two-body system\n");
		exit(2);
	}

	// declare variables
	double y[system.dim], y0[system.dim];
	double coordinate, velocity;
	int orbit_fw_size, orbit_bw_size;
	double **orbit_fw, **orbit_bw;
	double *par = (double *)system.params;
	double e = par[1];
	double orb[4];
	init_orbital(orb, e);

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
							y[k+2] = orb[k];
						}
					}

					// keep IC for backward integration
					copy(y0, y, system.dim);

					// calculate forward integration
					evolve_orbit(y, analysis.cycle_period, 
						analysis.number_of_cycles, &orbit_fw, 
						&orbit_fw_size, system);

					// calculate backward integration
					evolve_orbit(y, -1.0 * analysis.cycle_period, 
						analysis.number_of_cycles, &orbit_bw, 
						&orbit_bw_size, system);

					#pragma omp critical
					{
						// write initial condition to file
						fprintf(inc, "%1.15e %1.15e\n", coordinate, 
								velocity);
						// write orbit and constant error to file
						for (int k = 0; k < orbit_fw_size; k++)
						{
							fprintf(psp, "%1.15e %1.15e\n", 
								angle_mod(orbit_fw[k][0]), 
								orbit_fw[k][1]);
						}
						for (int k = 0; k < orbit_bw_size; k++)
						{
							fprintf(psp, "%1.15e %1.15e\n", 
								angle_mod(orbit_bw[k][0]), 
								orbit_bw[k][1]);
						}
						// create new line on exit file
						fprintf(psp, "\n");

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
	fclose(psp);
	fclose(inc);

	printf("Data written in output folder\n");

	return 0;
}
