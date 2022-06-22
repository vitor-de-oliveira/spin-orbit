#include "aux_vmo.h"

int array_to_matrix(double **m, double *x, int dim)
{
	int k = 0;
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			m[i][j] = x[k];
			k++;
		}
	}
	return 0;
}

int copy(double *x, double *y, int dim)
{
	for (int i = 0; i < dim; i++)
	{
		x[i] = y[i];
	}
	return 0;
}

int copy_parts(double *x, double *y, int initial_value, int final_value)
{
	for (int i = initial_value; i < final_value; i++)
	{
		x[i] = y[i];
	}
	return 0;
}

int print_array(double *x, int dim)
{
	for (int i = 0; i < dim; i++)
	{
		printf("x[%d] = %1.15e\n", i, x[i]);
	}
	return 0;
}

void print_prog(double percentage)
{
	int val = (int)(percentage * 100);
	int lpad = (int)(percentage * PBWIDTH);
	int rpad = PBWIDTH - lpad;
	printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
	fflush(stdout);
}

double sign(double x)
{
	if (x >= 0)
	{
		return 1.0;
	}
	else
	{
		return -1.0;
	}
}

int split_section(double **x, double *x_initial, double *x_final, int number, int dim)
{
	for (int i = 0; i < number; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			x[i][j] = (x_final[j] - x_initial[j]) / ((double)number - 1.0) * (double)i + x_initial[j];
		}
	}
	return 0;
}

int switch_array_element(double *x, int i, int j)
{
	double temp = x[i];
	x[i] = x[j];
	x[j] = temp;
	return 0;
}

double angle_mod_pos(double x)
{
	double result;
	result = fmod (x , 2.0 * M_PI);
	if (result < 0.0)
	{
		result += 2.0 * M_PI;
	}
	return result;
}

double angle_mod(double x)
{
	double result;
	result = fmod (x , 2.0 * M_PI);
	if (result > M_PI)
	{
		result -= 2.0 * M_PI;
	}
	else if (result < -M_PI)
	{
		result += 2.0 * M_PI;
	}
	return result;
}

double angular_dist(double x1, double x2)
{
	double result;
	double dif, dif_2pi;
	dif = angle_mod_pos(x1 - x2);
	dif_2pi = fabs(2.0*M_PI - dif);
	if (dif < dif_2pi)
	{
		result = dif;
	}
	else
	{
		result = dif_2pi;
	}
	return result;
}

double array_norm(double *x, int dim)
{
	double norm_sq = 0;
	for (int i = 0; i < dim; i++)
	{
		norm_sq += x[i] * x[i];
	}
	return sqrt(norm_sq);
}

double dist(double *x, double *y, int dim)
{
	double d_squared = 0.;
	for (int i = 0; i < dim; i++)
	{
		d_squared += (x[i] - y[i]) * (x[i] - y[i]);
	}
	return sqrt(d_squared);
}

double dist_partial(double *x, double *y, int index_1, int index_2)
{
	double d_squared = 0.;
	for (int i = index_1; i <= index_2; i++)
	{
		d_squared += (x[i] - y[i]) * (x[i] - y[i]);
	}
	return sqrt(d_squared);
}

int greatest_common_divisor (int n1, int n2)
{
	int gcd = 1;

    for(int i = 2; i <= n1 && i <= n2; i++)
    {
        // Checks if i is factor of both integers
        if((n1 % i == 0) && (n2 % i == 0))
		{
            gcd = i;
		}
	}

	if (gcd == 1)
	{
		printf("Warning: Greatest common divisor not found.\n");
	}

	return gcd;
}

int identity_matrix_array_form(double *x, int arr_first_index, int order)
{
	for (int i = arr_first_index; i < arr_first_index + order * order; i++)
	{
		x[i] = 0.0;
	}
	for (int i = arr_first_index; i < arr_first_index + order * order; i += order + 1)
	{
		x[i] = 1.0;
	}
	return 0;
}

void linear_combination(double *v, double a, double *x, double b, double *y, int dim)
{
	for (int i = 0; i < dim; i++)
	{
		v[i] = a * x[i] + b * y[i];
	}
}

int lin_comb_vec(double *x, double *lamb, double *y, double *z, int dim)
{
	for (int i = 0; i < dim; i++)
	{
		x[i] = lamb[i] * y[i] + z[i];
	}
	return 0;
}

void square_matrix_product_vector(double *mx, double **m, double *x, int dim)
{
	for (int i = 0; i < dim; i++)
	{
		mx[i] = 0.0;
		for (int j = 0; j < dim; j++)
		{
			mx[i] += m[i][j] * x[j];
		}
    }
}

void square_matrix_product_2d(double **mp, double **m1, double **m2, int dim)
{
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			mp[i][j] = 0.0;
			for (int k = 0; k < dim; k++)
			{
				mp[i][j] += m1[i][k] * m2[k][j];
			}
		}
	}
}

void square_matrix_tranpose_2d(double **mt, double **m, int dim)
{
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			mt[j][i] = m[i][j];
		}
	}
}

int alloc_1d_double(double **x, int n)
{
	*x = (double *)malloc(n * sizeof(double));
	return 0;
}

int alloc_2d_double(double ***x, int n, int m)
{
	*x = (double **)malloc(n * sizeof(double *));
	for (int i = 0; i < n; i++)
	{
		(*x)[i] = (double *)malloc(m * sizeof(double));
	}
	return 0;
}

int alloc_2d_int(int ***x, int n, int m)
{
	*x = (int **)malloc(n * sizeof(int *));
	for (int i = 0; i < n; i++)
	{
		(*x)[i] = (int *)malloc(m * sizeof(int));
	}
	return 0;
}

int alloc_1d_int(int **x, int n)
{
	*x = (int *)malloc(n * sizeof(int));
	return 0;
}

int alloc_3d_double(double ****x, int n, int m, int p)
{
	*x = (double ***)malloc(n * sizeof(double **));
	for (int i = 0; i < n; i++)
	{
		(*x)[i] = (double **)malloc(m * sizeof(double *));
		for (int j = 0; j < m; j++)
		{
			(*x)[i][j] = (double *)malloc(p * sizeof(double));
		}
	}
	return 0;
}

int dealloc_1d_double(double **x)
{
	free(*x);
	return 0;
}

int dealloc_1d_int(int **x)
{
	free(*x);
	return 0;
}

int dealloc_2d_double(double ***x, int n)
{
	for (int i = 0; i < n; i++)
	{
		free((*x)[i]);
	}
	free(*x);
	return 0;
}

int dealloc_2d_int(int ***x, int n)
{
	for (int i = 0; i < n; i++)
	{
		free((*x)[i]);
	}
	free(*x);
	return 0;
}

int dealloc_3d_double(double ****x, int n, int m)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			free((*x)[i][j]);
		}
		free((*x)[i]);
	}
	free(*x);
	return 0;
}

int realloc_1d_double(double **x, int n)
{
	*x = (double *)realloc(*x, n * sizeof(double));
	return 0;
}

int realloc_2d_double(double ***x, int n, int m, int o)
{
	*x = (double **)realloc(*x, n * sizeof(double *));
	for (int j = o; j < n; j++)
	{
		(*x)[j] = (double *)malloc(m * sizeof(double));
	}
	return 0;
}
