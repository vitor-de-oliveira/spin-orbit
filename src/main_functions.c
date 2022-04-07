#include "main_functions.h"

int trace_orbit_map(double *ic, void *params, 
		double cycle_period, int number_of_cycles, 
		char *system)
{
	// declare and open exit files
	FILE *orb = fopen("output/orbit.dat", "w");

	// declare variables
	int orbit_size;
	double **orbit;

	// evolve system
	evolve_orbit(params, ic, cycle_period, number_of_cycles, 
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
		/ (double)number_of_cycles);

	// free memory
	dealloc_2d_double(&orbit, number_of_cycles);

	// close files
	fclose(orb);

	printf("Data written in folder output\n");

	return 0;
}

int draw_phase_space(void *params, double cycle_period, 
		int number_of_cycles, double coordinate_min, 
		double coordinate_max, double velocity_min, 
		double velocity_max, int nc, int nv, char *system)
{
	// declare and open exit files
	FILE *psp = fopen("output/phase_space.dat", "w");
	FILE *inc 
		= fopen("output/phase_space_initial_conditions.dat", 
				"w");

	int system_dimension;

	if (strcmp(system, "rigid") == 0)
	{
		system_dimension = 6;
	}
	else if (strcmp(system, "rigid_kepler") == 0)
	{
		system_dimension = 2;
	}
	else if (strcmp(system, "two_body") == 0)
	{
		system_dimension = 4;
	}
	else
	{
		printf("Warning: undefined system\n");
		exit(2);
	}

	// declare variables
	double y[system_dimension], y0[system_dimension];
	double coordinate, velocity;
	int orbit_fw_size, orbit_bw_size;
	double **orbit_fw, **orbit_bw;
	double *par = (double *)params;
	double e = par[1];

		// loop over coordinate values
		for (int i = 0; i < nc; i++)
		{
			// print progress on coordinate
			printf("Calculating set %d of %d\n", i + 1, nc);

			if (nc == 1)
			{
				coordinate = coordinate_min;
			}
			else
			{
				coordinate = coordinate_min + 
						(double)(i) * (coordinate_max - 
						coordinate_min) / (double)(nc - 1);
			}

			// #pragma omp parallel default(none) private(params, 
			// 		cycle_period, number_of_cycles, 
			// 		coordinate_min, coordinate_max, 
			// 		velocity_min, velocity_max, nc, nv, system,
			// 		system_dimension, y, y0, 
			// 		coordinate, velocity, orbit_size, 
			// 		orbit_fw, orbit_bw, e) shared(psp, inc)
			// {
			// #pragma omp for
				// loop over velocity values
				for (int j = 0; j < nv; j++)
				{
					// print progress on velocity
					// print_prog((double)(j + 1) / (double)nv);
					printf("Calculating subset %d of %d\n", 
								j + 1, nv);

					if (nv == 1)
					{
						velocity = velocity_min;
					}
					else
					{
						velocity = velocity_min + 
								(double)(j) * (velocity_max - 
								velocity_min) / (double)(nv - 1);
					}

					y[0] = coordinate;
					y[1] = velocity;

					if (strcmp(system, "rigid") == 0)
					{
						y[2] = 0.9451; // Moon
						y[3] = 0.0;
						y[4] = 0.0;
						y[5] = sqrt((1.0 + e) / 
						0.9451);
					}

					// keep IC for backward integration
					copy(y0, y, system_dimension);

					// calculate forward integration
					evolve_orbit(params, y, cycle_period, 
						number_of_cycles, &orbit_fw, 
						&orbit_fw_size, system);

					// calculate backward integration
					evolve_orbit(params, y, -1.0 * cycle_period, 
						number_of_cycles, &orbit_bw, 
						&orbit_bw_size, system);

					// #pragma omp critical
					// {
						// write initial condition to file
						fprintf(inc, "%1.15e %1.15e\n", coordinate, 
								velocity);
						// write orbit and constant error to file
						for (int i = 0; i < orbit_fw_size; i++)
						{
							fprintf(psp, "%1.15e %1.15e\n", 
								angle_mod(orbit_fw[i][0]), 
								orbit_fw[i][1]);
						}
						for (int i = 0; i < orbit_bw_size; i++)
						{
							fprintf(psp, "%1.15e %1.15e\n", 
								angle_mod(orbit_bw[i][0]), 
								orbit_bw[i][1]);
						}
					// } // end pragma omp critical

					// write orbit and constant error to file


					// free memory
					dealloc_2d_double(&orbit_fw, number_of_cycles);
					dealloc_2d_double(&orbit_bw, number_of_cycles);

					// create new line on exit file
					fprintf(psp, "\n");
				
				}
			// } // end pragma

			// new line on terminal
			printf("\n");
		}

	// close exit files
	fclose(psp);
	fclose(inc);

	printf("Data written in folder output\n");

	return 0;
}