/**
 * personal library created by VMO to assist
 * writting c programming codes
**/

#ifndef AUX_VMO_H
#define AUX_VMO_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define PBSTR "=================================================="
#define PBWIDTH 50

/*************************** makes life easier ********************************/

// defines m as a square matrix of size dim by dim
// with elelents copied from the array x of size dim squared
int array_to_matrix(double **m, double *x, int dim);

// defines x as a copy of y
// where x and y are arrays of dimension dim
int copy(double *x, double *y, int dim);

// defines an interval of x as a copy of the same interval of y
// the copied indices go from initial_value to final_value
int copy_parts(double *x, double *y, int initial_value, int final_value);

// prints array of doubles of size dim
int print_array(double *x, int dim);

// prints the progress of a calculation as a filling bar
void print_prog(double percentage);

// returns one if n is positive and minus one if n is negative
double sign(double x);

// divides the dim-dimensional line connecting
// x_initial to x_final in number pieces
// and storages them in x
int split_section(double **x, double *x_initial, double *x_final, int number, int dim);

// switches the i-th and j-th elements of an array x
int switch_array_element(double *x, int i, int j);

/****************************** mathematics ***********************************/

// returns a value from 0 to 2PI for a given angle
double angle_mod_pos(double x);

// returns a value from -PI to PI for a given angle
double angle_mod(double x);

// returns the distance between two angles modulo 2PI
double angular_dist(double x1, double x2);

// return the norm of an real array x of dimension dim
double array_norm(double *x, int dim);

// returns the euclidean distance between
// arrays x and y of dimension dim
double dist(double *x, double *y, int dim);

// returns the euclidean distance between
// a range of indicies of arrays x and y
double dist_partial(double *x, double *y, int index_1, int index_2);

// returns the greatest common divisor of two integers
int greatest_common_divisor (int n1, int n2);

// defines x as a two-dimensional real-valued identity matrix
// of size order by order in array form starting from index
// arr_first_index
int identity_matrix_array_form(double *x, int arr_first_index, int order);

// v = a * x + b * y
void linear_combination(double *v, double a, double *x, double b, double *y, int dim);

// defines x as a linear combination of y and z
// where x, y, z and lamb are arrays of dimension dim
int lin_comb_vec(double *x, double *lamb, double *y, double *z, int dim);

// mx = m * x
void square_matrix_product_vector(double *mx, double **m, double *x, int dim);

// returns the product of two two-dimensional
// real square matrices of size dim by dim
void square_matrix_product_2d(double **mp, double **m1, double **m2, int dim);

// returns the transpose of a real two-dimensional 
// square matrix of size dim by dim
void square_matrix_tranpose_2d(double **mt, double **m, int dim);

/**************************** memory handling *********************************/

// allocates memory for one-dimensional double array x of size n
int alloc_1d_double(double **x, int n);

// allocates memory for one-dimensional integer array x of size n
int alloc_1d_int(int **x, int n);

// allocates memory for two-dimensional double array x of size n by n
int alloc_2d_double(double ***x, int n, int m);

// allocates memory for two-dimensional integer array x of size n by n
int alloc_2d_int(int ***x, int n, int m);

// allocates memory for three-dimensional double array x of size n by n by n
int alloc_3d_double(double ****x, int n, int m, int p);

// deallocates memory for one-dimensional double array x of size n
int dealloc_1d_double(double **x);

// deallocates memory for one-dimensional integer array x of size n
int dealloc_1d_int(int **x);

// deallocates memory for two-dimensional double array x of size n by n
int dealloc_2d_double(double ***x, int n);

// deallocates memory for two-dimensional integer array x of size n by n
int dealloc_2d_int(int ***x, int n);

// deallocates memory for three-dimensional double array x of size n by n by n
int dealloc_3d_double(double ****x, int n, int m);

// reallocates memory for one-dimensional double array x of size n
int realloc_1d_double(double **x, int n);

// reallocates memory for two-dimensional double array x of size n by n
int realloc_2d_double(double ***x, int n, int m, int o);

#endif