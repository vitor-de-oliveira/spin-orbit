#include "auxiliar_functions_gsl.h"

int set_driver(gsl_odeiv2_driver **d, 
			   gsl_odeiv2_system *sys, 
			   const gsl_odeiv2_step_type *T, 
			   double h, double h_max, double h_min,
			   double error_abs, double error_rel, 
			   char *control)
{
	*d = gsl_odeiv2_driver_alloc_y_new(sys, 
		T, h, error_abs, error_rel);
		
	if (strcmp(control, "adaptive") == 0)
	{
		gsl_odeiv2_driver_set_hmax(*d, h_max);
		gsl_odeiv2_driver_set_hmin(*d, h_min);
	}
	else if (strcmp(control, "fixed") != 0)
	{
		printf("Warning: invalid GSL error control\n");
		exit(2);
	}
	return 0;
}

int set_integrator(const gsl_odeiv2_step_type **T, 
				   char *integrator)
{
    if (strcmp(integrator, "rk2") == 0)
	{
		*T = gsl_odeiv2_step_rk2;
	}
	else if (strcmp(integrator, "rk4") == 0)
	{
		*T = gsl_odeiv2_step_rk4;
	}
	else if (strcmp(integrator, "rk45") == 0)
	{
		*T = gsl_odeiv2_step_rkf45;
	}
	else if (strcmp(integrator, "rkck") == 0)
	{
		*T = gsl_odeiv2_step_rkck;
	}
	else if (strcmp(integrator, "rk8pd") == 0)
	{
		*T = gsl_odeiv2_step_rk8pd;
	}
	else if (strcmp(integrator, "rk4imp") == 0)
	{
		*T = gsl_odeiv2_step_rk4imp;
	}
	else if (strcmp(integrator, "bsimp") == 0)
	{
		*T = gsl_odeiv2_step_bsimp;
	}
	else if (strcmp(integrator, "msadams") == 0)
	{
		*T = gsl_odeiv2_step_msadams;
	}
	else if (strcmp(integrator, "msbdf") == 0)
	{
		*T = gsl_odeiv2_step_msbdf;	
	}
	else
	{
		printf("Warning: invalid GSL integrator\n");
		exit(2);
	}
	return 0;
}

int set_system(gsl_odeiv2_system *sys, void *par, 
			   dynsys system)
{
	sys->function = system.field;
	sys->jacobian = system.jac;
	sys->dimension = system.dim;
	sys->params = par;
	return 0;
}