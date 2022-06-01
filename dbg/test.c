#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "dynamical_system.h"
#include "spin_orbit.h"
#include "pendulum.h"

#include "test_spin_orbit.h"

int main(int argc, char **argv)
{
	/******************** Start clock **********************/

	clock_t begin = clock(), end;
	#ifdef _OMP_H
		double begin_time_omp = omp_get_wtime(), end_time_omp;
	#endif

	/***************** Declared variables ******************/

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

	dynsys system;
	dynsys system_two_body = init_two_body(*params);
	dynsys system_rigid = init_rigid(*params);
	dynsys system_rigid_kepler = init_rigid_kepler(*params);
	dynsys system_linear = init_linear(*params);
	dynsys system_linear_average = init_linear_average(*params);
	
	dynsys system_pendulum = init_pendulum(*params);

	anlsis analysis;

	int grid[2];
	double orbital[4], basin[2];

	/*******************************************************/

	/////////////////////////////////////////////////////////
	/*				Dynamical systems functions				*/
	/////////////////////////////////////////////////////////

	// basin[0] = angle_mod_pos(2.0*M_PI);
	// basin[1] = 0.0;

	// analysis.grid_resolution = 10;
	// analysis.grid_coordinate_min = 0.0;
	// analysis.grid_coordinate_max = 2.0 * M_PI;
	// analysis.grid_velocity_min = 0.0;
	// analysis.grid_velocity_max = 3.0;

	// double_to_grid(grid, basin, analysis);

	// printf("i = %d j = %d\n", grid[0], grid[1]);

	// grid_to_double(grid, basin, analysis);

	// printf("x = %1.3f y = %1.3f\n", basin[0], basin[1]);

	// grid[0] = 9;
	// grid[1] = 0;

	// grid_to_double(grid, basin, analysis);

	// printf("x = %1.3f y = %1.3f\n", basin[0], basin[1]);

	// analysis.nc = 5, analysis.nv = 15; //nc = 11, nv = 30;
	// analysis.number_of_cycles = 1e3; //1e3
	// analysis.cycle_period = 2.0 * M_PI;
	// analysis.coordinate_min = 0.0;
	// analysis.coordinate_max = M_PI; // M_PI 2.0* M_PI
	// analysis.velocity_min = 0.0;
	// analysis.velocity_max = 3.0;
	// analysis.grid_resolution = 100;
	// analysis.grid_coordinate_min = -M_PI; // 0.0
	// analysis.grid_coordinate_max = M_PI; // 2.0*M_PI
	// analysis.grid_velocity_min = 0.0;
	// analysis.grid_velocity_max = 3.0;

	// FILE *out_test = fopen("output/test_grid.dat", "w");

	// for (int i = 0; i < analysis.grid_resolution; i++)
	// {
	// 	// int i = 0;
	// 	// double coordinate = analysis.grid_coordinate_min + 
	// 	// 	(double)(i) * (analysis.grid_coordinate_max - 
	// 	// 	analysis.grid_coordinate_min) / 
	// 	// 	(double)(analysis.grid_resolution - 1);
	// 	for (int j = 0; j < analysis.grid_resolution; j++)
	// 	{
	// 		// int j = 9;
	// 		// double velocity = analysis.grid_velocity_min + 
	// 		// 	(double)(j) * (analysis.grid_velocity_max - 
	// 		// 	analysis.grid_velocity_min) / 
	// 		// 	(double)(analysis.grid_resolution - 1);

	// 		grid[0] = i;
	// 		grid[1] = j;

	// 		grid_to_double(grid, basin, analysis);
	// 		double_to_grid(grid, basin, analysis);

	// 		// fprintf(out_test, 
	// 		// 	"coordinate = %1.3e velocity = %1.3e basin[0] = %1.3e basin[1] = %1.3e diff_c = %1.10e diff_v = %1.10e\ni = %d j = %d grid[0] = %d grid[1] = %d diff_i = %d diff_j = %d\n", 
	// 		// 	coordinate, velocity,
	// 		// 	basin[0],
	// 		// 	basin[1],
	// 		// 	coordinate - basin[0],
	// 		// 	velocity - basin[1],
	// 		// 	i, j,
	// 		// 	grid[0], grid[1],
	// 		// 	i - grid[0],
	// 		// 	j - grid[1]);

	// 		// fprintf(out_test, 
	// 		// 	"i = %d j = %d grid[0] = %d grid[1] = %d\n", 
	// 		// 	i, j,
	// 		// 	grid[0], grid[1]);

	// 		fprintf(out_test, 
	// 			"%d %d\n", 
	// 			i - grid[0], j - grid[1]);
	// 	}
	// }

	// fclose(out_test);

	// system = system_linear_average;

	// gamma = (.89 * .89) / 3.;
	// e = 0.1;
	// m_secondary = 0.0;
	// m_primary = 1.0 - m_secondary;
	// G = 1.0;
	// a = 1.0;
 	// K = 1e-3;

	// analysis.nc = 5, analysis.nv = 15;
	// analysis.number_of_cycles = 2e3;
	// analysis.cycle_period = 2.0 * M_PI;
	// analysis.coordinate_min = 0.0;
	// analysis.coordinate_max = M_PI; 
	// analysis.velocity_min = 0.0;
	// analysis.velocity_max = 3.0;
	// analysis.grid_resolution = 300;
	// analysis.grid_coordinate_min = -M_PI;
	// analysis.grid_coordinate_max = M_PI;
	// analysis.grid_velocity_min = 0.0;
	// analysis.grid_velocity_max = 3.0;

	// double ic[system.dim];
	// ic[0] = 0.2;
	// ic[1] = 0.1;
	// init_orbital(orbital, e);
	// for (int i = 0; i < 4; i++) ic[i+2] = orbital[i];

	// test_trace_orbit_map(ic, system, analysis);


	/////////////////////////////////////////////////////////
	/*		  System and orbital motion analysis    	   */
	/////////////////////////////////////////////////////////

	// int system_dimension;
	// int orbit_size;
	// double **orbit;
	// double t;

	// gamma = 0.01;
	// e = 0.0549; 					//Moon
	// m_secondary = 1.215e-2; 		//Moon
	// m_primary = 1.0 - m_secondary;
	// G = 1.0;

	// t = 0.0;

	// // number_of_cycles = 2e3 + 1;
	// // number_of_cycles = 2.5e4;
	// // number_of_cycles = 6e3;
	// // number_of_cycles = 1e3;
	// number_of_cycles = 6e3;

	// cycle_period = -2.0 * M_PI;
	// system_dimension = 6;
	// system = "rigid";

	// // cycle_period = 2.0 * M_PI;
	// // // cycle_period = M_PI / 100.0;
	// // system_dimension = 4;
	// // system = "two_body";

	// double perigee_over_semimajor_axis = 0.9451; //Moon

	// double y[system_dimension];
	// y[0] = 0.1;
	// y[1] = 0.1;
	// // y[2] = perigee_over_semimajor_axis;
	// y[3] = 0.0;
	// y[4] = 0.0;
	// // y[5] = sqrt(G * mass * (1.0 + e) / 
	// // 			perigee_over_semimajor_axis);

	// double a = 1.0;
	// y[2] = a * (1.0 - e * e) / (1.0 + e);
	
	// // printf("y[2]  = %1.15e\n", y[2]);

	// double f_e = atan2(y[3], y[2]);
	// // double y_dot = (e + cos(f_e))/sqrt(1.0-e*e);
	// double y_dot = (e + 1.0/sqrt(y[3]*y[3]/(y[2]*y[2])+1.0)) 
	// 				/sqrt(1.0-e*e);
	// // printf("y[5]  = %1.15e\ny_dot = %1.15e\n", y[5], y_dot);

	// y[5] = y_dot;

	// // y[0] = perigee_over_semimajor_axis;
	// // y[1] = 0.0;
	// // y[2] = 0.0;
	// // y[3] = sqrt(G * mass * (1.0 + e) / 
	// // 			perigee_over_semimajor_axis);


	// // evolve_cycle(y, &e, cycle_period, &t);
	// // printf("%1.10e %1.10e\n", y[0] , y[1]);

	// FILE *tst_ev_orbit
	// 		= fopen("output/test_evolve_orbit.dat", "w");
	// FILE *tst_orbit_const
	// 		= fopen("output/test_orbit_angular_moment.dat", "w");
	// FILE *tst_orbit_const_2
	// 		= fopen("output/test_orbit_vis_viva.dat", "w");
	// FILE *tst_orbit_true_anomaly
	// 		= fopen("output/test_orbit_true_anomaly.dat", "w");
	// FILE *tst_orbit_radius
	// 		= fopen("output/test_orbit_radius.dat", "w");

	// evolve_orbit(*params, y, cycle_period, number_of_cycles, 
	// 			&orbit, &orbit_size, system);

	// double orbital_motion[4], orbital_motion_ini[4];

	// // printf("%d\n", orbit_size);

	// for (int i = 0; i < orbit_size; i++)
	// {
	// 	// fprintf(tst_ev_orbit, "%1.15e %1.15e\n", 
	// 	// 	fmod(orbit[i][0], 2.0*M_PI), orbit[i][1]);
	// 	fprintf(tst_ev_orbit, "%1.15e %1.15e\n", 
	// 		angle_mod(orbit[i][0]), orbit[i][1]);
	// 	for (int j = 0; j < 4; j++)
	// 	{
	// 		orbital_motion[j] = orbit[i][j+2];
	// 		// orbital_motion[j] = orbit[i][j];
	// 	}
	// 	if (i == 0)
	// 	{
	// 		copy(orbital_motion_ini, orbital_motion, 4);
	// 	}
	// 	fprintf(tst_orbit_const, "%1.15e %1.15e\n", 
	// 		(double)i * cycle_period, 
	// 		fabs(angular_momentum_two_body(orbital_motion) - 
	// 		angular_momentum_two_body(orbital_motion_ini)));
	// 	double T = 2.0 * M_PI;
	// 	double a = 1.0;
	// 	fprintf(tst_orbit_const_2, "%1.15e %1.15e\n", 
	// 		(double)i * cycle_period, 
	// 		fabs(vis_viva_two_body(orbital_motion, T, a) - 
	// 		vis_viva_two_body(orbital_motion_ini, T, a)));
	// 	double f_e = 
	// 		atan2(orbital_motion[1], orbital_motion[0]);
	// 	double r = 
	// 		sqrt((orbital_motion[0] * orbital_motion[0]) + 
	// 			 (orbital_motion[1] * orbital_motion[1]));
	// 	fprintf(tst_orbit_true_anomaly, "%1.15e %1.15e\n", 
	// 			(double)i * cycle_period, f_e);
	// 	fprintf(tst_orbit_radius, "%1.15e %1.15e\n", 
	// 			(double)i * cycle_period, r);
	// }

	// printf("w_test = %1.10e\n", 
	// 	angular_dist(orbit[orbit_size-1][0], orbit[0][0])
	// 	/ (double)number_of_cycles);

	// dealloc_2d_double(&orbit, number_of_cycles);

    // fclose(tst_ev_orbit);
	// fclose(tst_orbit_const);
	// fclose(tst_orbit_const_2);
	// fclose(tst_orbit_true_anomaly);
	// fclose(tst_orbit_radius);

	/////////////////////////////////////////////////////////
	/*				   	Kepler equation	    			  */
	/////////////////////////////////////////////////////////

    // FILE *tst_kep_plane 
	// 		= fopen("output/test_kepler_plane.dat", "w");
    // FILE *tst_kep_time 
	// 		= fopen("output/test_kepler_time.dat", "w");
	// FILE *tst_kep_comparison 
	// 		= fopen("output/test_kepler_comparison.dat", "w");
	// FILE *tst_kep_true_anomaly
	// 		= fopen("output/test_kepler_true_anomaly.dat", "w");
	// FILE *tst_kep_vis_viva
	// 		= fopen("output/test_kepler_vis_viva.dat", "w");
	// FILE *tst_kep_ang_mom
	// 		= fopen("output/test_kepler_angular_mom.dat", "w");

    // e = 0.0549; //Moon
	// int count = 0;

	// double ini[4];

    // for (t = 0.0; t < 20.0 * M_PI; t += M_PI / 100.0)
    // {
    //     double u = kepler_equation(e,t);
    //     double r = 1.0 - e * cos(u);
    //     double f = 2.0 * atan(sqrt((1.0 + e)/(1.0 - e)) 
	// 			* tan(0.5 * u));
    //     double x = r * cos(f);
    //     double y = r * sin(f);
    //     fprintf(tst_kep_plane, "%1.15e %1.15e\n", x , y);
    //     fprintf(tst_kep_comparison, "%1.15e %1.15e\n", 
	// 			orbit[count][0] - x , orbit[count][1] - y);
	// 	fprintf(tst_kep_time, "%1.15e %1.15e\n", t , r);
	// 	fprintf(tst_kep_true_anomaly, "%1.15e %1.15e\n",
	// 			 t , f);
	// 	double vv[4];
	// 	vv[0] = x;
	// 	vv[1] = y;
	// 	vv[2] = -sin(f)/sqrt(1.0-e*e);
	// 	vv[3] = (e + cos(f))/sqrt(1.0-e*e);
	// 	if (count == 0)
	// 	{
	// 		copy(ini,vv,4);
	// 	}
	// 	double T = 2.0 * M_PI;
	// 	double a = 1.0;
	// 	fprintf(tst_kep_vis_viva, "%1.15e %1.15e\n", t, 
	// 		fabs(vis_viva_two_body(vv, T, a) - 
	// 			 vis_viva_two_body(ini, T, a)));
	// 	fprintf(tst_kep_ang_mom, "%1.15e %1.15e\n", t, 
	// 		fabs(angular_momentum_two_body(vv) - 
	// 			 angular_momentum_two_body(ini)));
	// 	count++;
	// }

    // fclose(tst_kep_plane);
    // fclose(tst_kep_time);
	// fclose(tst_kep_comparison);
	// fclose(tst_kep_true_anomaly);
	// fclose(tst_kep_vis_viva);
	// fclose(tst_kep_ang_mom);

	// dealloc_2d_double(&orbit, number_of_cycles);

	////////////////////////////////////////////////////////
	/*						 Pendulum   			   */
	////////////////////////////////////////////////////////

	// system = system_pendulum;

	// gamma = (.89 * .89) / 3.;
	// e = 0.1;
	// m_secondary = 0.0;
	// m_primary = 1.0 - m_secondary;
	// G = 1.0;
	// a = 1.0;
 	// K = 1e-6;

	// analysis.nc = 3, analysis.nv = 50; //nc = 3, nv = 50;
	// analysis.number_of_cycles = 1e3; //1e3
	// analysis.cycle_period = 1e-1;
	// analysis.coordinate_min = -M_PI; // M_PI
	// analysis.coordinate_max = M_PI; // M_PI 2.0* M_PI
	// analysis.velocity_min = -2.0;
	// analysis.velocity_max = 2.0;
	// analysis.grid_resolution = 1000;
	// analysis.grid_coordinate_min = -M_PI; // 0.0
	// analysis.grid_coordinate_max = M_PI; // 2.0*M_PI
	// analysis.grid_velocity_min = -2.0;
	// analysis.grid_velocity_max = 2.0;

	// // draw_phase_space_grid(system, analysis);

	// double x[2], ref[2];
	// ref[0] = 0.0;
	// ref[1] = 0.0;
	// FILE *test_dis = 
	// 	fopen("output/test_dist_from_ref.dat","w");

	// for (double theta = -M_PI; theta < M_PI; theta += 1e-3)
	// {
	// 	for (double theta_dot = -M_PI; theta_dot < M_PI; theta_dot += 1e-3)
	// 	{
	// 		x[0] = theta;
	// 		x[1] = theta_dot;
	// 		double d = dist_from_ref(x, ref);
	// 		fprintf(test_dis, "%1.5f %1.5f %1.5f\n", 
	// 				theta, theta_dot, d);
	// 	}
	// 	fprintf(test_dis, "\n");
	// }

	// fclose(test_dis);

	////////////////////////////////////////////////////////
	/*					 Test Spin Orbit	  			   */
	////////////////////////////////////////////////////////

	// system = system_linear_average;
  	// double ref[2], ref_union[2][2];
  	// int num_of_basins;

	// gamma = (.89 * .89) / 3.;
	// e = 0.1;
	// m_secondary = 0.0;
	// m_primary = 1.0 - m_secondary;
	// G = 1.0;
	// a = 1.0;
 	// K = 1e-2;

	// analysis.nc = 5, analysis.nv = 15; //nc = 11, nv = 30;
	// analysis.number_of_cycles = 1e3; //1e3
	// analysis.cycle_period = 2.0 * M_PI;
	// analysis.coordinate_min = 0.0;
	// analysis.coordinate_max = M_PI; // M_PI 2.0* M_PI
	// analysis.velocity_min = 0.0;
	// analysis.velocity_max = 3.0;
	// analysis.grid_resolution = 600;
	// analysis.grid_coordinate_min = -M_PI;
	// analysis.grid_coordinate_max = M_PI;
	// analysis.grid_velocity_min = 0.0;
	// analysis.grid_velocity_max = 3.0;

 	// ref[0] = M_PI; ref[1] = 0.551540;
  	// // ref[0] = 0.0; ref[1] = 0.551540;
  	// // basin_of_attraction (ref, system, analysis);
	// // basin_of_attraction_no_opt (ref, system, analysis);
	// basin_of_attraction_no_grid (ref, system, analysis);
	// // basin_of_attraction_no_omp (ref, system, analysis);
	// // basin_of_attraction_no_grid_no_omp (ref, system, analysis);

	// // num_of_basins = 2;
	// // ref_union[0][0] = 0.0; ref_union[0][1] = 0.551540;
	// // ref_union[1][0] = M_PI; ref_union[1][1] = 0.551540;
	// // union_basin_of_attraction (ref_union, num_of_basins, system, analysis);
  

	/******************** Stop clock ***********************/

	end = clock();
	double time_spent 
		= (double)(end - begin) / CLOCKS_PER_SEC;
	if (time_spent < 60.0)
	{
		printf("time spent = %1.2e seconds\n", 
				time_spent);
	}
	else if (time_spent < 3600.0)
	{
		printf("time spent = %1.2e minutes\n", 
				time_spent/60.0);
	}
	else
	{
		printf("time spent = %1.2e hours\n", 
				time_spent/3600.0);
	}
	
	#ifdef _OMP_H
		end_time_omp = omp_get_wtime();
		double time_spent_omp
			= (double)(end_time_omp - begin_time_omp);
		if (time_spent_omp < 60.0)
		{
			printf("time_spent_omp = %1.2e seconds\n", 
					time_spent_omp);
		}
		else if (time_spent_omp < 3600.0)
		{
			printf("time_spent_omp = %1.2e minutes\n", 
					time_spent_omp/60.0);
		}
		else
		{
			printf("time_spent_omp = %1.2e hours\n", 
					time_spent_omp/3600.0);
		}
	#endif


	/*******************************************************/

	return 0;
}
