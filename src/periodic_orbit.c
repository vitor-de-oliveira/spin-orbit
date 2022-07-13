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
    po.evolve_n_cycles(x_plus, po.period, system, analysis);

    x_minus[0] = po.initial_condition[0] - 0.5 * dx;
    x_minus[1] = po.initial_condition[1];
    po.evolve_n_cycles(x_minus, po.period, system, analysis);

    J[0][0] = (x_plus[0] - x_minus[0]) / dx;
    J[1][0] = (x_plus[1] - x_minus[1]) / dx;

    // approximate 2nd col. of jacobian matrix
    x_plus[0] = po.initial_condition[0];
    x_plus[1] = po.initial_condition[1] + 0.5 * dx;
    po.evolve_n_cycles(x_plus, po.period, system, analysis);

    x_minus[0] = po.initial_condition[0];
    x_minus[1] = po.initial_condition[1] - 0.5 * dx;
    po.evolve_n_cycles(x_minus, po.period, system, analysis);

    J[0][1] = (x_plus[0] - x_minus[0]) / dx;
    J[1][1] = (x_plus[1] - x_minus[1]) / dx;

    // free memory
    dealloc_1d_double(&x);
    dealloc_1d_double(&x_plus);
    dealloc_1d_double(&x_minus);
}

void minimization_step  (double *x,
                         double *err,
                         double lamb,
                         perorb po,
                         dynsys system, 
                         anlsis analysis)
{
    double *y, *rhs, *dx, *ic;
    double **M, **M_inv, **Jt, **J;

    alloc_1d_double(&y, 2);
    alloc_1d_double(&dx, 2);
    alloc_1d_double(&rhs, 2);
    alloc_1d_double(&ic, 2);
    alloc_2d_double(&M, 2, 2);
    alloc_2d_double(&M_inv, 2, 2);
    alloc_2d_double(&J, 2, 2);
    alloc_2d_double(&Jt, 2, 2);

    jacobian_periodic_orbit(J, po, system, analysis);

    // functional jacobian
    J[0][0] -= 1.0; J[1][1] -= 1.0;
    square_matrix_transpose_2d(Jt, J, 2);
    square_matrix_product_2d(M, Jt, J, 2);
    M[0][0] *= (1.0+lamb); M[1][1] *= (1.0+lamb);

    // right hand side vector
    copy(ic, po.initial_condition, 2);
    po.evolve_n_cycles(ic, po.period, system, analysis);
    linear_combination(y, 1.0, po.initial_condition, -1.0, ic, 2);
    y[0] = angle_mod(y[0]);
    square_matrix_product_vector(rhs, Jt, y, 2);
    // gauss_solve (dx, M, rhs, 2);
    square_2d_matrix_inverse(M_inv, M);
    square_matrix_product_vector(dx, M_inv, rhs, 2);

    // new point
    linear_combination(x, 1.0, po.initial_condition, 1.0, dx, 2);
    x[0] = angle_mod(x[0]);

    // new error
    copy(ic, x, 2);
    po.evolve_n_cycles(ic, po.period, system, analysis);
    *err = po.dist_on_phase_space(ic, x);

    dealloc_1d_double(&y);
    dealloc_1d_double(&dx);
    dealloc_1d_double(&rhs);
    dealloc_1d_double(&ic);
    dealloc_2d_double(&M, 2);
    dealloc_2d_double(&M_inv, 2);
    dealloc_2d_double(&J, 2);
    dealloc_2d_double(&Jt, 2);
}

int calculate_periodic_orbit_ic(perorb *po,
                                dynsys system, 
                                anlsis analysis)
{
    int count, max_steps;
    double r, lamb, err, err1, err2, tol;
    double *x1, *x2;

    // Levenberg-Marquardt parameters
    r = 2.0;
    lamb = 0.01;

    // iteration parameters
    tol = 1e-10;
    max_steps = 1000;

    // dynamical memory allocation 
    alloc_1d_double(&x1, 2);
    alloc_1d_double(&x2, 2);

    // error of the seed
    printf("seed for initial condition = (%1.5e  %1.5e)\n",
            (*po).seed[0], (*po).seed[1]);
    copy ((*po).initial_condition, (*po).seed, 2);
    (*po).evolve_n_cycles((*po).initial_condition, (*po).period, 
                    system, analysis);
    err = (*po).dist_on_phase_space((*po).initial_condition, (*po).seed);
    printf("seed error = |M^%d(seed)-seed| = %1.5e\n", 
            (*po).period, err);

    // fixed point calculation 
    printf("Starting periodic orbit calculation\n");

    copy ((*po).initial_condition, (*po).seed, 2);
    count = 0;
    while (err > tol && count < max_steps)
    {
        printf ("step = %d err = %1.5e lamb = %f\n",
                count, err, lamb);

        minimization_step(x1, &err1, lamb, *po, system, analysis);

        minimization_step(x2, &err2, r * lamb, *po, system, analysis);

        printf("err1 = %1.5e err2 = %1.5e\n", err1, err2);

        if (err1 < err && err1 < err2)
        {
            err = err1;
            copy ((*po).initial_condition, x1, 2);
        }
        else if (err2 < err && err2 < err1)
        {
            err = err2;
            copy ((*po).initial_condition, x2, 2);
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
                (*po).initial_condition[0], (*po).initial_condition[1]);
        printf("error = |M^%d(po_ic)-po_ic| = %1.5e\n", 
                (*po).period, err);
        exit(2);
    }
    else
    {
        printf("periodic orbit of period %d found after %d steps\n", (*po).period, count);
        printf("initial condition = (%1.15e  %1.15e)\n", 
                (*po).initial_condition[0], (*po).initial_condition[1]);
        printf("error = |M^%d(po_ic)-po_ic| = %1.5e\n",
                (*po).period, err);
    }

    dealloc_1d_double(&x1);
    dealloc_1d_double(&x2);

    return 0;
}

void gauss_solve (double *x, double **A, double *b, int dim)
{
    double temp;

    bool back_sub = true;

    double** Ab; // augmented matrix
    Ab = malloc (dim * sizeof(double*));

    alloc_2d_double(&Ab, dim, dim + 1);

    for (int i = 0; i < dim; i++)
    {
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

    dealloc_2d_double(&Ab, dim);
}