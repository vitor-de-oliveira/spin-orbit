#include "pendulum.h"

int field_pendulum  (double t, const double y[],
                    double f[], void *params)
{
	(void)t;

	f[0] = y[1];
	f[1] = -1.0 * sin(y[0]);

	return GSL_SUCCESS;
}

dynsys init_pendulum(void *params)
{
	dynsys pendulum;
    pendulum.name = "pendulum";
    pendulum.dim = 2;
	pendulum.field = &field_pendulum;
	pendulum.jac = NULL;
	pendulum.params = params;
    return pendulum;
}

int draw_phase_space_grid(dynsys system, anlsis analysis)
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

	// declare variables
	FILE 	*out_psp;
	double y[system.dim], y0[system.dim];
	double coordinate, velocity;
	int orbit_fw_size, orbit_bw_size;
	double **orbit_fw, **orbit_bw;
	double *par = (double *)system.params;
	double e = par[1];
	double orb[4], orb_ini[4];
	// init_orbital(orb_ini, e);

	int grid[2];
	double ps[2];
	int **ps_matrix;
	alloc_2d_int(&ps_matrix, 
		analysis.grid_resolution, analysis.grid_resolution);
	for (int i = 0; i < analysis.grid_resolution; i++)
	{
		for (int j = 0; j < analysis.grid_resolution; j++)
		{
			ps_matrix[i][j] = 0;
		}
	}

	// open exit files
	out_psp = fopen("output/test_phase_space_grid.dat", "w");

	// loop over coordinate values
	for (int i = 0; i < analysis.nc; i++)
	{
		// print progress on coordinate
		printf("Calculating set %d of %ld\n", i + 1, analysis.nc);

		// #pragma omp parallel private(y, y0, coordinate, velocity, \
		// 		orbit_fw_size, orbit_bw_size, orbit_fw, orbit_bw)
		// {

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

		// #pragma omp for
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

				// if (strcmp(system.name, "rigid") == 0)
				// {
				// 	for (int k = 0; k < 4; k++)
				// 	{
				// 		y[k+2] = orb_ini[k];
				// 	}
				// }

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

				// write orbit and error to file
				for (int k = 0; k < orbit_fw_size; k++)
				{
					ps[0] = angle_mod(orbit_fw[k][0]);
					ps[1] = orbit_fw[k][1];
					double_to_grid(grid, ps, analysis);
					
					if ((grid[0] >= 0) && 
						 (grid[0] < analysis.grid_resolution) && 
						(grid[1] >= 0) && 
						 (grid[1] < analysis.grid_resolution))
					{
						// printf("Here 1\n");
						// printf("grid[0] = %d\n", grid[0]);
						// printf("grid[1] = %d\n", grid[1]);
						ps_matrix[grid[0]][grid[1]] = 1;
						// printf("Here 2\n");
					}
				}
				for (int k = 0; k < orbit_bw_size; k++)
				{
					ps[0] = angle_mod(orbit_bw[k][0]);
					ps[1] = orbit_bw[k][1];
					double_to_grid(grid, ps, analysis);
					
					if ((grid[0] >= 0 && 
						 grid[0] < analysis.grid_resolution) && 
						(grid[1] >= 0 && 
						 grid[1] < analysis.grid_resolution))
					{
						ps_matrix[grid[0]][grid[1]] = 1;
					}
				}

				// free memory
				dealloc_2d_double(&orbit_fw, 
						analysis.number_of_cycles);
				dealloc_2d_double(&orbit_bw, 
						analysis.number_of_cycles);
			
			}
		// } // end pragma

		// new line on terminal
		printf("\n");
	}

	for (int i = 0; i < analysis.grid_resolution; i++)
	{
		for (int j = 0; j < analysis.grid_resolution; j++)
		{
			grid[0] = i;
			grid[1] = j;
			grid_to_double(grid, ps, analysis);

			fprintf(out_psp, "%1.5f %1.5f %d\n", 
				ps[0], ps[1], 
				ps_matrix[grid[0]][grid[1]]);
		}
		fprintf(out_psp, "\n");
	}
	// close exit files
	fclose(out_psp);

	dealloc_2d_int(&ps_matrix, 
					analysis.grid_resolution);

	printf("Data written in output folder\n");

	return 0;
}
