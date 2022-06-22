#include "periodic_orbit.h"

void jacobian_periodic_orbit(double **J,
                             perorb po,
                             dynsys system, 
                             anlsis analysis)
{
    double  *x, *x_plus, *x_minus;
    double  dx = 1e-7;

    // allocation on heap
    alloc_1d_double(&x, 2);
    alloc_1d_double(&x_plus, 2);
    alloc_1d_double(&x_minus, 2);

    // approximate 1st col. of jacobian matrix
    x_plus[0] = po.initial_condition[0] + 0.5 * dx;
    x_plus[1] = po.initial_condition[1];
    evolve_n_cycles (po.period, x_plus, system, analysis);

    x_minus[0] = po.initial_condition[0] - 0.5 * dx;
    x_minus[1] = po.initial_condition[1];
    evolve_n_cycles (po.period, x_minus, system, analysis);

    J[0][0] = (x_plus[0] - x_minus[0]) / dx;
    J[1][0] = (x_plus[1] - x_minus[1]) / dx;

    // approximate 2nd col. of jacobian matrix
    x_plus[0] = po.initial_condition[0];
    x_plus[1] = po.initial_condition[1] + 0.5 * dx;
    evolve_n_cycles (po.period, x_plus, system, analysis);

    x_minus[0] = po.initial_condition[0];
    x_minus[1] = po.initial_condition[1] - 0.5 * dx;
    evolve_n_cycles (po.period, x_minus, system, analysis);

    J[0][1] = (x_plus[0] - x_minus[0]) / dx;
    J[1][1] = (x_plus[1] - x_minus[1]) / dx;

    // free memory
    dealloc_1d_double(&x);
    dealloc_1d_double(&x_plus);
    dealloc_1d_double(&x_minus);
}

void minimization_step  (double *x,
                         double lamb,
                         perorb po,
                         dynsys system, 
                         anlsis analysis)
{
    double *y, *rhs, *dx, *ic;
    double **M, **Jt, **J;

    alloc_1d_double(&y, 2);
    alloc_1d_double(&dx, 2);
    alloc_1d_double(&rhs, 2);
    alloc_1d_double(&ic, 2);
    alloc_2d_double(&M, 2, 2);
    alloc_2d_double(&J, 2, 2);
    alloc_2d_double(&Jt, 2, 2);

    jacobian_periodic_orbit(J, po, system, analysis);

    // functional jacobian
    J[0][0] -= 1.0; J[1][1] -= 1.0;
    square_matrix_tranpose_2d(Jt, J, 2);
    square_matrix_product_2d(M, Jt, J, 2);
    M[0][0] *= (1.0+lamb); M[1][1] *= (1.0+lamb);

    // right hand side vector
    copy(ic, po.initial_condition, 2);
    evolve_n_cycles (po.period, po.initial_condition,
                    system, analysis);
    linear_combination(y, 1.0, po.initial_condition, -1.0, ic, 2);
    square_matrix_product_vector(rhs, Jt, y, 2);
    gauss_solve (M, dx, rhs, 2);

    // new point
    linear_combination(x, 1.0, po.initial_condition, 1.0, dx, 2);

    dealloc_1d_double(&y);
    dealloc_1d_double(&dx);
    dealloc_1d_double(&rhs);
    dealloc_1d_double(&ic);
    dealloc_2d_double(&M, 2);
    dealloc_2d_double(&J, 2);
    dealloc_2d_double(&Jt, 2);
}

int periodic_orbit  (double seed[2],
                     perorb po,
                     dynsys system, 
                     anlsis analysis)
{
    if (po.period > 100)
    {
        printf("Warning: cannot calculate\n");
        printf("a periodic orbit with this period.\n");
        exit(2);
    }

    int count, max_steps;
    double r, lamb, err, err1, err2, tol;
    double *x1, *x2;

    // Levenberg-Marquardt parameters
    r = 0.5;
    lamb = 0.01;

    // iteration parameters
    tol = 1e-14;
    max_steps = 100;

    // dynamical memory allocation 
    alloc_1d_double(&x1, 2);
    alloc_1d_double(&x2, 2);

    // error of the seed
    printf("seed for initial condition = (%1.15e  %1.15e)\n",
            seed[0], seed[1]);
    copy (po.initial_condition, seed, 2);
    evolve_n_cycles (po.period, po.initial_condition, system, 
                     analysis);
    err = dist(po.initial_condition, seed, 2);
    printf("seed error = |M^%d(seed)-seed| = %1.5e\n", 
            po.period, err);

    // fixed point calculation 
    printf("Starting periodic orbit calculation\n");

    copy (po.initial_condition, seed, 2);
    count = 0;
    while (err > tol && count < max_steps)
    {
        printf ("step = %d err = %1.5e lamb = %f\n",
                count, err, lamb);

        minimization_step(x1, lamb, po, system, analysis);
        err1 = dist(x1, po.initial_condition, 2);

        minimization_step(x2, r * lamb, po, system, analysis);
        err2 = dist(x2, po.initial_condition, 2);

        if (err1 < err && err1 < err2)
        {
            err = err1;
            copy (po.initial_condition, x1, 2);
        }
        else if (err2 < err && err2 < err1)
        {
            err = err2;
            copy (po.initial_condition, x2, 2);
            r *= lamb;
        }
        else if (err2 < err1)
        {
            lamb *= r;
        }
        else
        {
            lamb /= r;
        }
        count ++;
    }

    // if max count was reached
    if (count == max_steps)
    {
        printf("method didn't converge\n");
	    printf("max. count reached during fixed point search\n");
        printf("approx. initial condition = (%1.15e  %1.15e)\n", 
                po.initial_condition[0], po.initial_condition[1]);
        printf("error = |M^%d(po_ic)-po_ic| = %1.5e\n", 
                po.period, err);
    }
    else
    {
        printf("periodic orbit of period %d found after %d steps\n", po.period, count);
        printf("initial condition = (%1.15e  %1.15e)\n", 
                po.initial_condition[0], po.initial_condition[1]);
        printf("error = |M^%d(po_ic)-po_ic| = %1.5e\n",
                po.period, err);
    }

    // save values on perorb
    fills_periodic_orbit(po, system, analysis);

    dealloc_1d_double(&x1);
    dealloc_1d_double(&x2);

    return 0;
}

int fills_periodic_orbit(perorb po,
                         dynsys system, 
                         anlsis analysis)
{
    double t = 0;
    double y[2];

    copy(y, po.initial_condition, 2);

    for (int i = 0; i < po.period; i++)
    {
        copy(po.orbit[i], y, 2);
        evolve_cycle(y, &t, system, analysis);
    }

    return 0;
}

int evolve_n_cycles (int    n,
                     double *y,
                     dynsys system,
                     anlsis analysis)
{
    double t = 0;

    // loop over cycles
	for (int i = 0; i < n; i++)
	{
		evolve_cycle(y, &t, system, analysis);
	
		// check if orbit diverges
		for (int j = 0; j < system.dim; j++)
		{
			if (fabs(y[j]) > analysis.evolve_box_size)
			{
                printf("Warning: evolve n cycles failed\n");
				printf("box limit reached\n");
				printf("y[%d] = %1.10e\n", j, y[j]);
				exit(2);
			}
		}
	}

	return 0;
}

void gauss_solve (double **A, double *x, double *b, int dim)
{
    double temp;

    bool back_sub = true;

    double** Ab; // augmented matrix
    Ab = malloc (dim * sizeof(double*));

    for (int i = 0; i < dim; i++)
    {
        Ab[i] = malloc ((dim + 1) * sizeof(double));
        for (int j = 0; j < dim; j++)
        {
            Ab[i][j] = A[i][j];
        }
        Ab[i][dim] = b[i];
    }

    for (int k = 0; k < dim; k++) // move in the pivots
    {
        if (Ab[k][k] == 0.0) // if zero swaps with a nonzero
        {
            int l = k + 1;
            while (Ab[l][k] == 0 && l < dim)
            {
                l++;
            }
            if (l == dim)
            {
                back_sub = false;
                break;
            }
            for (int j = 0; j < dim + 1; j++)
            {
                temp = Ab[k][j];
                Ab[k][j] = Ab[l][j];
                Ab[l][j] = temp;
            }
        }
        for (int j = dim; j >= k; j--)
        {
            Ab[k][j] /= Ab[k][k];
        }
        for (int i = k + 1; i < dim; i++)
        {
            if (Ab[i][k] != 0.0)
            {
                for (int j = dim; j >= k; j--)
                {
                    Ab[i][j] = Ab[i][j] - Ab[i][k]*Ab[k][j];
                }
            }
        }
    }
    
    if (back_sub) // performs the back_substitution
    {
        for (int i = dim - 1; i >= 0; i--)
        {
            temp = Ab[i][dim];
            for (int j = i + 1; j < dim; j++)
            {
                temp -= Ab[i][j]*x[j];
            }
            x[i] = temp;
        }
    }
}