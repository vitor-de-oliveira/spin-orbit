#ifndef TEST_SPIN_ORB_H
#define TEST_SPIN_ORB_H

#include "spin_orbit.h"

int basin_of_attraction_no_opt  (double *ref,
                                dynsys system,
                                anlsis analysis);                        

int basin_of_attraction_no_grid(double *ref,
                                dynsys system,
                                anlsis analysis);

int basin_of_attraction_no_omp  (double *ref,
                                dynsys system,
                                anlsis analysis);

int basin_of_attraction_no_opt_no_omp  (double *ref,
                                        dynsys system,
                                        anlsis analysis);


int basin_of_attraction_no_grid_no_omp  (double *ref,
                                        dynsys system,
                                        anlsis analysis);

int trace_orbit_map_image   (double *ic,
                            dynsys system,
                            anlsis analysis);

#endif
