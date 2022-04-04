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

#endif