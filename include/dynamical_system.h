#ifndef DYN_SYS_H
#define DYN_SYS_H

#include <math.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include <kepler.h>

// field for rotational frame
int field_two_body(double t, const double y[], double f[], 
                   void *params);

// jacobian for rotational frame
int jacobian_two_body(double t, const double y[], 
            double *dfdy, double dfdt[], void *params);

// field for rotational frame
int field_rigid(double t, const double y[], double f[], 
            void *params);

// jacobian for rotational frame
int jacobian_rigid(double t, const double y[], double *dfdy, 
            double dfdt[], void *params);

// field for rotational frame
int field_rigid_kepler(double t, const double y[], 
            double f[], void *params);

// jacobian for rotational frame
int jacobian_rigid_kepler(double t, const double y[], 
            double *dfdy, double dfdt[], void *params);

// field for rotational frame
double angular_momentum_two_body(double y[4]);

double vis_viva_two_body(double y[4], double T, double a);

int init_orbital(double y[4], double e);

// Dynamical System
typedef struct DynSys{
    char *name;
    int (*field) (double t, const double y[], double dydt[], void *params);
    int (*jac) (double t, const double y[], double *dfdy, double dfdt[],
                   void *params);
    size_t dim;
    void *params;
} dynsys;

dynsys init_rigid(void *params);

dynsys init_rigid_kepler(void *params);

dynsys init_two_body(void *params);

// System analysis
typedef struct AnlSis{
    size_t nc, nv;
	int number_of_cycles;
	double cycle_period;
	double coordinate_min;
	double coordinate_max;
	double velocity_min;
	double velocity_max;
} anlsis;

#endif