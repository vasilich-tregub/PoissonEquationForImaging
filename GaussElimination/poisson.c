/* TERMS OF USE
 * This source code is subject to the terms of the MIT License.
 * Copyright(c) 2026 Vladimir Vasilich Tregub
*/
#include "poisson.h"
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>

int main(int argc, char* argv[]) {
    int ihor = 20, ivert = 20;
    int uhor = ihor + 1, uvert = ivert + 1; // grid dimensions including boundaries
    int ahor = ihor - 1, avert = ivert - 1; // grid dimensions excluding boundaries
    double* u = calloc(uvert * uhor, sizeof(double)); // solution to find
    double* f = calloc(uvert * uhor, sizeof(double)); // charge distribution for Poisson eq
    int rank = ahor * avert; // length of solution vector
    double* a = calloc(rank * (rank + 1), sizeof(double)); // matrix A to construct for algebraic eq

    // let we have charge ditribution within solution range
    for (int iy = 0; iy <= ivert; ++iy) {
        for (int ix = 0; ix <= ihor; ++ix)
        {
            f[iy * uhor + ix] = -2.0 / (((ix - ihor / 2.0) * (ix - ihor / 2.0) + 0.01) * ((iy - ivert / 2.0) * (iy - ivert / 2.0) + 0.01));
            u[iy * uhor + ix] = 0;
        }
    }

    // boundary conditions
    for (int ix = 0; ix < uhor; ++ix) {
        u[0 * uhor + ix] = 0;
        u[ivert * uhor + ix] = -ix * (ihor - ix) * (1.0 / ihor);
    }
    for (int iy = 0; iy < uvert; ++iy) {
        u[iy * uhor + 0] = 0;
        u[iy * uhor + ihor] = iy * (ivert - iy) * (1.0 / ivert);
    }

    /*for (int iy = 0; iy < uvert; ++iy) {
        for (int ix = 0; ix < uhor; ++ix) {
            printf("%-7.3f ", u[iy * uhor + ix]);
        }
        printf("\n");
    }
    printf("\n");*/
    /* fill in matrix a with zeroes */
    for (int j = 0; j < rank; ++j)
        for (int i = 0; i <= rank; ++i)
            a[j * rank + i] = 0;

    // right-hand side column (#rank+1) for augmented matrix, a Poisson term contribution
    for (int iy = 1; iy < ivert; ++iy)
        for (int ix = 1; ix < ihor - 1; ++ix)
            a[((iy - 1) * ahor + ix - 1) * (rank + 1) + rank] = (1.0 / ihor) * (1.0 / ivert) * f[iy * (ihor + 1) + ix];

    // fill in the augmented matrix, including boundary coditions
    // 
    // boundary conds for outer border horizontals (iy = 0; iy = ivert)
    for (int ix = 1; ix < ihor; ++ix) {
        a[(ix - 1) * (rank + 1) + rank] -= u[ix];
        a[((ivert - 2) * ahor + (ix - 1)) * (rank + 1) + rank] -= u[ivert * uhor + ix];
    }

    // boundary conds for outer border verticals (ix = 0; ix = ihor)
    for (int iy = 1; iy < ivert; ++iy) {
        a[(iy - 1) * ahor * (rank + 1) + rank] -= u[iy * uhor + 0];
        a[((iy - 1) * ahor + ihor - 2) * (rank + 1) + rank] -= u[iy * uhor + ihor];
    }

    a[0] += -4;
    for (int ix = 1; ix < ihor - 1; ++ix)
    {
        a[(ix - 1) * (rank + 1) + ix] += 1;
        a[ix * (rank + 1) + ix - 1] += 1;
        a[ix * (rank + 1) + ix] += -4;
    }
    for (int iy = 1; iy < ivert - 1; ++iy) {
            a[((iy - 1) * ahor) * (rank + 1) + iy * ahor] += 1;
            a[(iy * ahor) * (rank + 1) + (iy - 1) * ahor] += 1;
            a[(iy * ahor) * (rank + 1) + iy * ahor] += -4;
        for (int ix = 1; ix < ihor - 1; ++ix) {
            a[(iy * ahor + ix - 1) * (rank + 1) + iy * ahor + ix] += 1;
            a[(iy * ahor + ix) * (rank + 1) + iy * ahor + ix - 1] += 1;
            a[((iy - 1) * ahor + ix) * (rank + 1) + iy * ahor + ix] += 1;
            a[(iy * ahor + ix) * (rank + 1) + (iy - 1) * ahor + ix] += 1;
            a[(iy * ahor + ix) * (rank + 1) + iy * ahor + ix] += -4;
        }
    }

    /*for (int irow = 0; irow < rank; ++irow) {
        for (int icol = 0; icol <= rank; ++icol) {
            printf("%-6.3f ", a[irow * (rank + 1) + icol]);
        }
        printf("\n");
    }*/

    double* x_arr = (double*)malloc(rank * sizeof(double));
    for (int i = 0; i < rank; ++i)
        x_arr[i] = 0;

    clock_t startSolver = clock();
    solvelinsys(a, rank + 1, rank, x_arr);
    clock_t finishSolver = clock();

    /*for (int iy = 1; iy < ivert; ++iy)
        for (int ix = 1; ix < ihor; ++ix)
            u[iy * (ivert + 1) + ix] = x_arr[(iy - 1) * (ihor - 1) + (ix - 1)];*/

    for (int iy = 1; iy < uvert - 1; ++iy) {
        for (int ix = 1; ix < uhor - 1; ++ix) {
            u[iy * uhor + ix] = x_arr[(iy - 1) * ahor + ix - 1];
        }
    }

    for (int iy = 0; iy < uvert; ++iy) {
        for (int ix = 0; ix < uhor; ++ix) {
            printf("%-6.2f ", u[iy * uhor + ix]);
        }
        printf("\n");
    }

    printf("Solved in %d ms\n", finishSolver - startSolver);

    return 0;
}