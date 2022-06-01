#include "test_spin_orbit.h"

int basin_of_attraction_no_opt	(double *ref, dynsys system, 
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
	out_boa = fopen("output/basin_no_opt.dat", "w");
	out_ref = fopen("output/basin_ref_no_opt.dat", "w");

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
					evolve_basin(y, ref, &converged,
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

int basin_of_attraction_no_grid(double *ref, dynsys system, 
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
	int grid[2];
	double orb[4], orb_ini[4];
	init_orbital(orb_ini, e);
	int basin_counter = 0;

	bool converged;

	// open exit files
	out_boa = fopen("output/basin_no_grid.dat", "w");
	out_ref = fopen("output/basin_ref_no_grid.dat", "w");

	fprintf(out_ref, "%1.15e %1.15e\n", ref[0], ref[1]);		

	// loop over coordinate values
	for (int i = 0; i < analysis.grid_resolution; i++)
	{
		// print progress on coordinate
		printf("Calculating set %d of %d\n", 
					i + 1, analysis.grid_resolution);

		omp_set_dynamic(0);     // Explicitly disable dynamic teams
		omp_set_num_threads(12); // Use 4 threads for all consecutive parallel regions

		#pragma omp parallel private(y, rot_ini, grid, \
				orbit_fw_size, orbit_fw, converged)
		{

		#pragma omp for
			// loop over velocity values
			for (int j = 0; j < analysis.grid_resolution; j++)
			{
				// print progress on velocity
				printf("Calculating subset %d of %d\n", 
							j + 1, analysis.grid_resolution);

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
						#pragma omp critical
						{
							fprintf(out_boa, "%1.5e %1.5e %d\n",
								rot_ini[0], rot_ini[1], orbit_fw_size);
							basin_counter++;
							// for (int k = 0; k < orbit_fw_size; k++)
							// {
							// 	fprintf(out_boa, "%1.5e %1.5e %d\n", 
							// 		angle_mod(orbit_fw[k][0]),
							// 		orbit_fw[k][1],
							// 		orbit_fw_size - 1 - k);
							// }
						}

					}

					// free memory
					dealloc_2d_double(&orbit_fw, 
							analysis.number_of_cycles);
			
			}
		} // end pragma

		fprintf(out_boa, "\n");

		// new line on terminal
		printf("\n");
	}

	printf("Basin counter = %d\n", basin_counter);

	// close exit files
	fclose(out_boa);
	fclose(out_ref);

	printf("Data written in output folder\n");

	return 0;
}

int basin_of_attraction_no_omp	(double *ref, dynsys system, 
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
	out_boa = fopen("output/basin_no_omp.dat", "w");
	out_ref = fopen("output/basin_ref_no_omp.dat", "w");

	fprintf(out_ref, "%1.15e %1.15e\n", ref[0], ref[1]);		

	// loop over coordinate values
	for (int i = 0; i < analysis.grid_resolution; i++)
	{
		// print progress on coordinate
		printf("Calculating set %d of %d\n", 
					i + 1, analysis.grid_resolution);

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
						basin_counter++;

						for (int k = 1; k < orbit_fw_size; k++)
						{
							// basin[0] = angle_mod_pos(orbit_fw[k][0]);
							// basin[0] = fmod(orbit_fw[k][0], M_PI);
							basin[0] = angle_mod(orbit_fw[k][0]);
							basin[1] = orbit_fw[k][1];
							double_to_grid(grid, basin, analysis);
							if (grid[0] >= 0 && 
								grid[0] < analysis.grid_resolution && 
								grid[1] >= 0 && 
								grid[1] < analysis.grid_resolution)
							{
								if (control_matrix[grid[0]][grid[1]] == 0)
								{
									basin_matrix[grid[0]][grid[1]] = 1;
									control_matrix[grid[0]][grid[1]] = 1;
									time_matrix[grid[0]][grid[1]] = 
										(double)(orbit_fw_size - k);
									basin_counter++;
								}
							}
						}
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

int basin_of_attraction_no_opt_no_omp	(double *ref,
							dynsys system, anlsis analysis)
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
	out_boa = fopen("output/basin_no_opt_no_omp.dat", "w");
	out_ref = fopen("output/basin_ref_no_opt_no_omp.dat", "w");

	fprintf(out_ref, "%1.15e %1.15e\n", ref[0], ref[1]);		

	// loop over coordinate values
	for (int i = 0; i < analysis.grid_resolution; i++)
	{
		// print progress on coordinate
		printf("Calculating set %d of %d\n", 
					i + 1, analysis.grid_resolution);

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

int basin_of_attraction_no_grid_no_omp(double *ref,
					dynsys system, anlsis analysis)
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
	int grid[2];
	double orb[4], orb_ini[4];
	init_orbital(orb_ini, e);
	int basin_counter = 0;

	bool converged;

	// open exit files
	out_boa = fopen("output/basin_no_grid_no_omp.dat", "w");
	out_ref = fopen("output/basin_ref_no_grid_no_omp.dat", "w");

	fprintf(out_ref, "%1.15e %1.15e\n", ref[0], ref[1]);		

	// loop over coordinate values
	for (int i = 0; i < analysis.grid_resolution; i++)
	{
		// print progress on coordinate
		printf("Calculating set %d of %d\n", 
					i + 1, analysis.grid_resolution);

			// loop over velocity values
			for (int j = 0; j < analysis.grid_resolution; j++)
			{
				// print progress on velocity
				printf("Calculating subset %d of %d\n", 
							j + 1, analysis.grid_resolution);

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
						fprintf(out_boa, "%1.5e %1.5e\n",
							rot_ini[0], rot_ini[1]);
						basin_counter++;
					}

					// free memory
					dealloc_2d_double(&orbit_fw, 
							analysis.number_of_cycles);
			
			}

		// new line on terminal
		printf("\n");
	}

	printf("Basin counter = %d\n", basin_counter);

	// close exit files
	fclose(out_boa);
	fclose(out_ref);

	printf("Data written in output folder\n");

	return 0;
}

int trace_orbit_map_image	(double *ic, dynsys system,
							anlsis analysis)
{
	// create output folder if it does not exist
	struct stat st = {0};
	if (stat("output", &st) == -1) {
		mkdir("output", 0700);
	}

	// declare variables
	int orbit_size;
	double **orbit;
	double orb[4], orb_ini[4];

	for (int i = 0; i < 4; i++)
	{
		orb_ini[i] = ic[i+2];
	}

	int grid[2];
	double real_space[2];
	int **grid_matrix, **control_matrix;
	alloc_2d_int(&grid_matrix, 
		analysis.grid_resolution, analysis.grid_resolution);
	alloc_2d_int(&control_matrix, 
		analysis.grid_resolution, analysis.grid_resolution);

	for (int i = 0; i < analysis.grid_resolution; i++)
	{
		for (int j = 0; j < analysis.grid_resolution; j++)
		{
			grid_matrix[i][j] = 0;
			control_matrix[i][j] = 0;
		}
	}

	// open exit files
	FILE *out_orb, *out_orb_ic,
		 *out_orb_ang_mom_err, *out_vis_viva_err,
		 *out_grid;
	out_orb = fopen("output/test_orbit.dat", "w");
	out_orb_ic = fopen("output/test_orbit_ic.dat", "w");
	out_orb_ang_mom_err = 
		fopen("output/test_orbit_orbital_angular_momentum_error.dat", "w");
	out_vis_viva_err = fopen("output/test_orbit_vis_viva_error.dat", "w");
	out_grid = fopen("output/test_orbit_grid.dat", "w");

	// evolve system
	evolve_orbit(ic, analysis.cycle_period, 
				analysis.number_of_cycles, 
				&orbit, &orbit_size, system);

	// write orbit and constant error to file
	fprintf(out_orb_ic, "%1.15e %1.15e\n", 
			angle_mod(orbit[0][0]), orbit[0][1]);
	for (int i = 0; i < orbit_size; i++)
	{
		fprintf(out_orb, "%1.15e %1.15e\n", 
				angle_mod(orbit[i][0]), orbit[i][1]);

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

		real_space[0] = angle_mod(orbit[i][0]);
		real_space[1] = orbit[i][1];
		double_to_grid(grid, real_space, analysis);
		if (grid[0] >= 0 && 
			grid[0] < analysis.grid_resolution && 
			grid[1] >= 0 && 
			grid[1] < analysis.grid_resolution)
		{
			if (control_matrix[grid[0]][grid[1]] == 0)
			{
				grid_matrix[grid[0]][grid[1]] = 1;
				control_matrix[grid[0]][grid[1]] = 1;
			}
		}
	}

	printf("w = %1.10e\n", 
		angular_dist(orbit[orbit_size-1][0], orbit[0][0])
		/ (double)analysis.number_of_cycles);

	for (int i = 0; i < analysis.grid_resolution; i++)
	{
		for (int j = 0; j < analysis.grid_resolution; j++)
		{
			grid[0] = i;
			grid[1] = j;
			grid_to_double(grid, real_space, analysis);

			fprintf(out_grid, "%1.5f %1.5f %d\n", 
				real_space[0], real_space[1], 
				grid_matrix[grid[0]][grid[1]]);
		}
		fprintf(out_grid, "\n");
	}

	// free memory
	dealloc_2d_double(&orbit, analysis.number_of_cycles);

	dealloc_2d_int(&grid_matrix, 
					analysis.grid_resolution);

	dealloc_2d_int(&control_matrix, 
					analysis.grid_resolution);

	// close files
	fclose(out_orb);
	fclose(out_orb_ic);
	fclose(out_orb_ang_mom_err);
	fclose(out_vis_viva_err);
	fclose(out_grid);

	printf("Data written in output folder\n");

	return 0;
}