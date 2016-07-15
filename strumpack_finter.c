#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <mpi.h>
#include <stdbool.h>
#include "StrumpackSparseSolver.h"

void strumpack_serial_finter(int *job, int *nrows, int *nnz,
                             int *irptr, int *jcptr, float complex *a,
                             float complex *rhs, float complex *sol,
                             int *ierr)
{
    const char *fcnm = "strumpack_serial_finter\0";
    static STRUMPACK_SparseSolver spss;
    STRUMPACK_RETURN_CODE info;
    int i, n;
    int lsymmetric = 1; // This matrix is structurally symmetric
    static bool linit = false;
    static bool lfactored = false;
    *ierr = 0;
    n = *nrows;
    // Strumpack initialization
    if (*job == 1)
    {
        if (linit)
        {
            printf("%s: Error already initialized!\n", fcnm);
            *ierr = 1;
            return;
        }
        STRUMPACK_init(&spss, 0, STRUMPACK_FLOATCOMPLEX, STRUMPACK_MT,
                       0, NULL, 1);
        // Don't use MC64 to reorder for numerical stability
        STRUMPACK_set_mc64job(spss, 0);
        // Control verbosity
        //STRUMPACK_set_verbose(spss, 0);
        // Set the accuracy 
        STRUMPACK_set_rctol(spss, 0.0001);  // default is 0.0001
        STRUMPACK_set_actol(spss, 1.e-10);  // default is 1.e-10 
        // Specify ordering for regular grids -> otherwise METIS
        STRUMPACK_set_reorder_method(spss, STRUMPACK_METIS); //GEOMETRIC);
        STRUMPACK_set_from_options(spss);
        // Renumber Fortran -> C
        for (i=0; i<n+1; i++){irptr[i] = irptr[i] - 1;}
        for (i=0; i<*nnz; i++){jcptr[i] = jcptr[i] - 1;}
        // Perform the reordering
        STRUMPACK_set_csr_matrix(spss, &n, irptr, jcptr, a, lsymmetric);
        info = STRUMPACK_reorder(spss);
        if (info != STRUMPACK_SUCCESS)
        {
            if (info == STRUMPACK_REORDERING_ERROR)
            {
                printf("%s: Error reordering matrix!\n", fcnm);
            }
            else
            {
                printf("%s: Strumpack error!\n", fcnm);
            }
            *ierr = 1;
        }
        // Renumber C -> Fortran
        for (i=0; i<n+1; i++){irptr[i] = irptr[i] + 1;}
        for (i=0; i<*nnz; i++){jcptr[i] = jcptr[i] + 1;}
        linit = true;
    }
    // Strumpack reset matrix and factorize 
    else if (*job == 2)
    {
        lfactored = false;
        if (!linit)
        {
            printf("%s: Error not initialized!\n", fcnm);
            *ierr = 1;
            return;
        }
        // Renumber Fortran -> C
        for (i=0; i<n+1; i++){irptr[i] = irptr[i] - 1;} 
        for (i=0; i<*nnz; i++){jcptr[i] = jcptr[i] - 1;}
        STRUMPACK_set_csr_matrix(spss, &n, irptr, jcptr, a, lsymmetric);
        // Renumber C -> Fortran
        for (i=0; i<n+1; i++){irptr[i] = irptr[i] + 1;} 
        for (i=0; i<*nnz; i++){jcptr[i] = jcptr[i] + 1;}
        info = STRUMPACK_factor(spss);
        if (info != STRUMPACK_SUCCESS)
        {
            printf("%s: Error factoring matrix!\n", fcnm);
            *ierr = 1;
            return;
        }
        lfactored = true;
    }
    else if (*job == 3)
    {
        if (!linit)
        {
            printf("%s: Error matrix not factorized!\n", fcnm);
            *ierr = 1;
            return;
        }
        int linit_guess = 0;
        info = STRUMPACK_solve(spss, rhs, sol, linit_guess);
        if (info != STRUMPACK_SUCCESS)
        {
            if (info == STRUMPACK_MATRIX_NOT_SET)
            {
                printf("%s: Error matrix not set\n", fcnm);
            }
            printf("%s: Error in solution phase1\n", fcnm);
            *ierr = 1;
            return;
        }
    }
    // Strumpack finalize
    else if (*job ==-1)
    {
        if (linit){STRUMPACK_destroy(&spss);}
        linit = false;
    }
    return;
}

int strumpack_parallel_finter( )
{

}
