#ifndef PEND_H
#define PEND_H

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_roots.h>

#include "dynamical_system.h"
#include "aux_vmo.h"

int field_pendulum  (double t,
                    const double y[],
                    double f[], 
                    void *params);

dynsys init_pendulum(void *params);

int draw_phase_space_grid   (dynsys system,
                            anlsis analysis);

#endif