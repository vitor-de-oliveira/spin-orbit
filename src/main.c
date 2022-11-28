#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "dynamical_system.h"
#include "periodic_orbit.h"
#include "spin_orbit.h"

int main(int argc, char **argv)
{
	/******************** Start clock **********************/

	clock_t begin_time = clock(), end_time;
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
	double T;				// system period

	double *params[8] = {&gamma,
						 &e,
						 &m_primary,
						 &m_secondary,
						 &G,
						 &a,
						 &K,
						 &T};

	dynsys system;
	dynsys system_two_body = init_two_body(*params);
	dynsys system_rigid = init_rigid(*params);
	dynsys system_rigid_kepler = init_rigid_kepler(*params);
	dynsys system_linear = init_linear(*params);
	dynsys system_linear_average = init_linear_average(*params);

	anlsis analysis;

	perorb po;

	double orbital[4];

	/********************* Some values **********************/

	/** Moon
	 * e = 0.0549
	 * m_secondary = 1.215e-2
	**/

	double e_moon = 0.0549;

	/** Hyperion (Wisdom 1987)
	 * e = 0.1 
	 * gamma = (.89 * .89) / 3.
	 * m_secondary = 0. (assuming for now)
	**/

	double e_hyperion = 0.1;
	double gamma_hyperion = ((.89 * .89) / 3.); //~0.264
	  
	/////////////////////////////////////////////////////////
	/*				   		   Orbit		   	           */
	/////////////////////////////////////////////////////////

	// // system = system_rigid;
	// // system = system_linear_average;
	// system = system_linear;
	// double ic[system.dim];

	// gamma = gamma_hyperion;
	// e = 0.186;
	// m_secondary = 0.;
	// m_primary = 1.0 - m_secondary;
	// G = 1.0;
	// a = 1.0;
	// K = 1e-1;
	// T = 2.0 * M_PI;

	// analysis.number_of_cycles = 2e3; //1e3 6e3
	// analysis.cycle_period = T; // 1e-3
	// analysis.evolve_box_size = 1e8;
	// analysis.evolve_basin_eps = 1e-1;

	// // // ic[0] = 0.0, ic[1] = 1000.;
	// // // // // near the 1:1 stable fp in the rigid case
	// // // // ic[0] = M_PI; ic[1] = 0.551537;
	// // // init_orbital(orbital, e);
	// // // for (int i = 0; i < 4; i++) ic[i+2] = orbital[i];

	// ic[0] = -0.237849, ic[1] = 0.569783;
	// init_orbital(orbital, e);
	// for (int i = 0; i < 4; i++) ic[i+2] = orbital[i];
	// orbit_map(ic, system, analysis);

	/////////////////////////////////////////////////////////
	/*				   	Periodic Orbit		   	           */
	/////////////////////////////////////////////////////////

	// system = system_rigid;
	// system = system_linear_average;
	system = system_linear;
	double ic_po[system.dim];

	gamma = gamma_hyperion;
	e = 1.313131313131292e-02;
	m_secondary = 0.;
	m_primary = 1.0 - m_secondary;
	G = 1.0;
	a = 1.0;
	K = 1e-1;
	T = 2.0 * M_PI;

	analysis.cycle_period = T;
	analysis.evolve_box_size = 1e8;

	// po.period = 1;
	// // po.seed[0] = 0.0; po.seed[1] = 1.0; // e = 0.0 SFP 1/1 resonance
	// // po.seed[0] = 0.0; po.seed[1] = 2.21320; // e = 0.01 SFP 2/1 resonance
	// // po.seed[0] = 0.0; po.seed[1] = 0.551537; // e = 0.1 SFP 1/1 resonance
	// // po.seed[0] = 0.0; po.seed[1] = 2.32185; // e = 0.1 SFP 2/1 resonance  
	// // po.seed[0] = -1.56892; po.seed[1] = 0.868688; // e = 0.1 period 2 SPO 1/2 resonance 
	// // po.seed[0] = 0.0; po.seed[1] = 1.87878; // e = 0.1 period 2 UPO 2/2 resonance
	// // po.seed[0] = -1.57310; po.seed[1] =  1.71059; // e = 0.1 UFP 2/1 resonance 
	// // po.seed[0] = -1.57246; po.seed[1] =  2.14877; // e = 0.1 UPO 5/2 resonance
	// // po.seed[0] = -1.94124; po.seed[1] =  1.46147; // e = 0.1 period 2 UPO resonance 3/2
	// // po.seed[0] = 1.35558; po.seed[1] =  1.08285; // e = 0.1 UPO 2/2
	// // po.seed[0] = 1.94124; po.seed[1] =  1.46147; // e = 0.1 period 2 UPO resonance 3/2
	// // po.seed[0] = 0.0; po.seed[1] =  2.72177; // e = 0.1 SPO 5/2
	// // po.seed[0] = -1.57079; po.seed[1] =  1.95929; // e = 0.1 SPO 9/4
	// // po.seed[0] = -0.0257629; po.seed[1] = 0.484803; // e = 0.140 period 3 UPO around 1/1 resonance
	// // po.seed[0] = 0.7; po.seed[1] = 1.1; // e = 0.2 ?
	// // po.seed[0] = 0.0; po.seed[1] = 0.380929; // e = 0.2 UFP 1/1 resonance
	// // po.seed[0] = 0.0; po.seed[1] = 0.722967; // e = 0.2 SFP 2/2 resonance
	// po.seed[0] = 3.823807668771570e-12; po.seed[1] = 2.222212216031596e+00; // test
	// periodic_orbit(&po, system, analysis);

	// draw_periodic_orbit_on_phase_space (po, system);
	// draw_periodic_orbit_on_phase_space_clean (po, system);

	// FILE 	*out_po_track_e;
	// char    filename[200];
	// int		po_counter, po_number = 1000;
	// double	po_tracker[po_number][3]; // e |lambda_1| |lambda_2|
	// double	e_max, e_min, e_base, e_step_up, e_step_down;
	
	// e_base = 0.1;
	// e_max = 0.2;
	// e_min = 0.0;
	// e_step_up = fabs(e_base - e_max) / ((double) (po_number - 1));
	// e_step_down = fabs(e_base - e_min) / ((double) (po_number - 1));

	// e = e_base;
	// po.period = 2;
	// po.seed[0] = -1.94124; po.seed[1] =  1.46147; // e = 0.1 period 2 UPO resonance 3/2
	
	// periodic_orbit(&po, system, analysis);
	
	// if (strcmp(system.name, "rigid") == 0)
	// {
	// 	sprintf(filename, "output/periodic_orbit/periodic_orbit_track_gamma_%1.3f_system_%s_resonance_%d_%d.dat", 
	// 		gamma, system.name, po.winding_number_numerator, po.winding_number_denominator);
	// }
	// else
	// {
	// 	sprintf(filename, "output/periodic_orbit/periodic_orbit_track_gamma_%1.3f_system_%s_K_%1.5f_resonance_%d_%d.dat", 
	// 		gamma, system.name, K, po.winding_number_numerator, po.winding_number_denominator);
	// }
	// out_po_track_e = fopen(filename, "w");

	// po_counter = 0;
	// for (e = e_base; e > e_min - e_step_down/2.; e -= e_step_down)
	// {

	// 	if ((po.seed[0] != po.seed[0]) ||
	// 		(po.seed[1] != po.seed[1]))
	// 	{
	// 		printf("Can't calculate PO due to bad seed\n");
	// 		goto jump;
	// 	}

	// 	periodic_orbit(&po, system, analysis);

	// 	jacobian_eigenvalues_magnitude(&po, system, analysis);

	// 	jump:;

	// 	po_tracker[po_counter][0] = e;
	// 	po_tracker[po_counter][1] = po.eigenvalues_absolute_value[0];
	// 	po_tracker[po_counter][2] = po.eigenvalues_absolute_value[1];

	// 	copy(po.seed, po.initial_condition, 2);

	// 	po_counter++;
	// }

	// for (int i = po_number - 1; i > 0; i--)
	// {
	// 	fprintf(out_po_track_e, "%1.5e %1.5e %1.5e\n", 
	// 		po_tracker[i][0], po_tracker[i][1], po_tracker[i][2]);
	// }

	// po.seed[0] = -1.94124; po.seed[1] =  1.46147; // e = 0.1 period 2 UPO resonance 3/2

	// po_counter = 0;
	// for (e = e_base; e < e_max + e_step_up/2.; e += e_step_up)
	// {

	// 	if ((po.seed[0] != po.seed[0]) ||
	// 		(po.seed[1] != po.seed[1]))
	// 	{
	// 		printf("Can't calculate PO due to bad seed\n");
	// 		goto jump2;
	// 	}

	// 	periodic_orbit(&po, system, analysis);

	// 	jacobian_eigenvalues_magnitude(&po, system, analysis);

	// 	jump2:;

	// 	po_tracker[po_counter][0] = e;
	// 	po_tracker[po_counter][1] = po.eigenvalues_absolute_value[0];
	// 	po_tracker[po_counter][2] = po.eigenvalues_absolute_value[1];

	// 	po.seed[0] = po.initial_condition[0];
	// 	po.seed[1] = po.initial_condition[1];

	// 	po_counter++;
	// }

	// for (int i = 0; i < po_number; i++)
	// {
	// 	fprintf(out_po_track_e, "%1.5e %1.5e %1.5e\n", 
	// 		po_tracker[i][0], po_tracker[i][1], po_tracker[i][2]);
	// }
	
	// fclose(out_po_track_e);

	FILE 	*out_bif;
	char    filename[200];
	int 	number_of_iter = 2000;		// 1000
	int		transient = 1900;			// 900
	int 	number_of_points = 10;		// 1000
	double 	x[system.dim], t;
	double	e_max, e_min, e_base;
	double	e_step_up, e_step_down;

	e_base = 0.186004;					// e_hyperion
	e_max = 0.186006;				// e_hyperion + 0.1;
	e_min = e_hyperion - 0.1;
	e_step_up = fabs(e_base - e_max) / ((double) number_of_points);
	e_step_down = fabs(e_base - e_min) / ((double) number_of_points);

	// e = e_base;
	// po.period = 1;
	// po.seed[0] = 0.0; po.seed[1] = 0.551537; // e = 0.1 SFP 1/1 resonance
	// periodic_orbit(&po, system, analysis);

	e = e_base;
	po.period = 2;
	po.seed[0] = -0.245487; po.seed[1] = 0.558370; // e = 0.1 SFP 1/1 resonance
	periodic_orbit(&po, system, analysis);

	sprintf(filename, "output/periodic_orbit/test_bifurcation_diagram_gamma_%1.3f_system_%s_K_%1.5f_resonance_%d_%d.dat", 
			gamma, system.name, K, po.winding_number_numerator, po.winding_number_denominator);
	out_bif = fopen(filename, "w");

	// copy (x, po.initial_condition, 2);
	// // init_orbital(orbital, e);
	// // for (int i = 0; i < 4; i++) x[i+2] = orbital[i];
	// for (e = e_base; e > e_min - e_step_down/2.; e -= e_step_down)
	// {
	// 	printf("e = %1.4f\n", e);
	// 	x[0] = angle_mod(x[0]);
	// 	init_orbital(orbital, e);
	// 	for (int i = 0; i < 4; i++) x[i+2] = orbital[i];
	// 	t = 0.0;
	// 	for (int i = 0; i < number_of_iter; i++)
	// 	{
	// 		evolve_cycle(x, &t, system, analysis);
	// 		if (i > transient)
	// 		{
	// 			fprintf(out_bif, "%1.4f %1.10e %1.10e\n", 
	// 				e, angle_mod(x[0]), x[1]);
	// 		}
	// 	}
	// }

	copy (x, po.initial_condition, 2);
	// init_orbital(orbital, e);
	// for (int i = 0; i < 4; i++) x[i+2] = orbital[i];
	for (e = e_base; e < e_max + e_step_up/2.; e += e_step_up)
	{
		printf("e = %1.4f\n", e);
		x[0] = angle_mod(x[0]);
		init_orbital(orbital, e);
		for (int i = 0; i < 4; i++) x[i+2] = orbital[i];
		t = 0.0;
		for (int i = 0; i < number_of_iter; i++)
		{
			evolve_cycle(x, &t, system, analysis);
			if (i > transient)
			{
				fprintf(out_bif, "%1.10e %1.10e %1.10e\n", 
					e, angle_mod(x[0]), x[1]);
			}
		}
	}

	fclose(out_bif);

	// analysis.number_of_cycles = 5e4;
	// analysis.cycle_period = 1e-3;

	// ic_po[0] = po.initial_condition[0];
	// ic_po[1] = po.initial_condition[1];
	// init_orbital(orbital, e);
	// for (int i = 0; i < 4; i++) ic_po[i+2] = orbital[i];
	// orbit_map(ic_po, system, analysis);

	// analysis.number_of_cycles = 2;
	// analysis.cycle_period = 2.0 * M_PI;
	// analysis.grid_resolution = 600;
	// analysis.grid_coordinate_min = -M_PI;
	// analysis.grid_coordinate_max = M_PI;
	// analysis.grid_velocity_min = 0.0;
	// analysis.grid_velocity_max = 3.0;
	// analysis.evolve_box_size = 1e6;

	// look_for_resonance (2, 3, system, analysis);

	// draw_orbit_map(system);

	// draw_orbit_on_phase_space(system);
	// draw_orbit_on_phase_space_latex(system);

	/////////////////////////////////////////////////////////
	/*				   		Time series		   	           */
	/////////////////////////////////////////////////////////

	// analysis.number_of_cycles = 4e3; //1e3 6e3
	// analysis.cycle_period = 2.0 * M_PI; // 1e-3
	// analysis.evolve_box_size = 1e8;
	// analysis.evolve_basin_eps = 1e-1;
	// analysis.number_of_time_series = 20;
	// analysis.time_series_delta = 10;

	// multiple_time_series(system, analysis);

	// draw_multiple_time_series(system);

	// multiple_time_series_delta_theta_dot(system, analysis);
	// draw_multiple_time_series_delta_theta_dot(system, analysis);
	// draw_multiple_time_series_delta_theta_dot_latex(system, analysis);

	// // multiple_time_series_delta_theta(system, analysis);
	// draw_multiple_time_series_delta_theta(system);

	// time_series(system, analysis);

	// // draw_time_series(system);
	// draw_time_series_latex(system);

	// analysis.number_of_cycles = 2e2;
	// for (e = 0.00; e < 0.205; e += 0.01)
	// {
	// 	time_series(system, analysis);
	// }

	// draw_time_series_union_e(system);
	// draw_time_series_union_e_latex(system);
	// draw_time_series_union_e_eps(system);

	// analysis.number_of_cycles = 1e3;
	// for (K = 0.01; K > 0.00095; K -= 0.001)
	// {
	// 	time_series(system, analysis);
	// }

	// draw_time_series_union_K(system);
	// draw_time_series_union_K_latex(system);

	/////////////////////////////////////////////////////////
	/*				   		Phase space		   	           */
	/////////////////////////////////////////////////////////

	// system = system_rigid;

	// gamma = gamma_hyperion;
	// e = 1.313131313131292e-02;
	// m_secondary = 0.0;
	// m_primary = 1.0 - m_secondary;
	// G = 1.0;
	// a = 1.0;

	// analysis.nc = 3, analysis.nv = 50; //nc = 3, nv = 50;
	// analysis.number_of_cycles = 2e3; //1e3
	// analysis.cycle_period = 2.0 * M_PI;
	// analysis.coordinate_min = 0.0; // M_PI
	// analysis.coordinate_max = M_PI; // M_PI 2.0* M_PI
	// analysis.velocity_min = 0.0;
	// analysis.velocity_max = 3.0;
	// analysis.evolve_box_size = 1e8;

	// phase_space(system, analysis);
	// // draw_phase_space(system);
	// draw_phase_space_latex(system);

	// gamma = gamma_hyperion;

	// for (e = 0.0; e < 0.405; e += 0.01)
	// {
	// 	// phase_space(system, analysis);
	// 	// draw_phase_space(system);
	// 	draw_phase_space_clean(system);
	// }

	// e = e_hyperion;

	// for (gamma = 0.0; 
	// 	 gamma < gamma_hyperion * 2. + gamma_hyperion / 40.;
	// 	 gamma += gamma_hyperion / 20.)
	// {
	// 	// phase_space(system, analysis);
	// 	// draw_phase_space(system);
	// 	draw_phase_space_clean(system);
	// }

	/////////////////////////////////////////////////////////
	/*				Basin of attraction		   	           */
	/////////////////////////////////////////////////////////

	// // system = system_linear_average;
	// system = system_linear;
	// int ref_period = 1;
	// double ref[ref_period][2];

	// gamma = gamma_hyperion;
	// e = e_hyperion;
	// m_secondary = 0.0;
	// m_primary = 1.0 - m_secondary;
	// G = 1.0;
	// a = 1.0;
 	// K = 1e-2;
	// T = 2.0 * M_PI;

	// analysis.number_of_cycles = 1e3; //1e3
	// analysis.cycle_period = T;
	// analysis.grid_resolution = 600;
	// analysis.grid_coordinate_min = -M_PI;
	// analysis.grid_coordinate_max = M_PI;
	// analysis.grid_velocity_min = 0.0;
	// analysis.grid_velocity_max = 3.0;
	// analysis.evolve_box_size = 1e8;
	// analysis.evolve_basin_time_tol = 100;
	// analysis.evolve_basin_eps = 1e-1;

	// po.period = 1;
	// po.seed[0] = 0.0; po.seed[1] = 0.551537; // e = 0.1 SFP 1/1 resonance
	// periodic_orbit(&po, system, analysis);
	// ref[0][0] = po.initial_condition[0]; ref[0][1] = po.initial_condition[1];
  	// // basin_of_attraction (ref, ref_period, system, analysis);
	// draw_basin_of_attraction (ref, ref_period, system, analysis);
	// // draw_basin_of_attraction_clean (ref, ref_period, system, analysis);

	/////////////////////////////////////////////////////////
	/*						Benchmark		   	           */
	/////////////////////////////////////////////////////////

	// linear_average_benchmark();

	/******************** Stop clock ***********************/

	end_time = clock();
	double time_spent 
		= (double)(end_time - begin_time) / CLOCKS_PER_SEC;
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