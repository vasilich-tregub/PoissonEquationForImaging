/* TERMS OF USE
 * This source code is subject to the terms of the MIT License.
 * Copyright(c) 2026 Vladimir Vasilich Tregub
*/
// PoissonWithCholesky.c : Defines the entry point for the application.

#include "PoissonWithCholesky.h"
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>

int main(int argc, char* argv[]) {
    int ihor = 20, ivert = 20;
    double h = 1.0 / (ihor + 1); //step
    double xleft = 0, xright = xleft + ihor * h, ybottom = 0, ytop = ivert * h + ybottom; //domain boundaries
    double* u = calloc((ivert + 1) * (ihor + 1), sizeof(double)); // solution to find
    double* f = calloc((ivert + 1) * (ihor + 1), sizeof(double)); // charge distr for Poisson eq
    int rank = (ihor - 1) * (ivert - 1); // length of solution vector

    double* b = calloc(rank, sizeof(double));

    for (int iy = 0; iy <= ivert; ++iy) {
        for (int ix = 0; ix <= ihor; ++ix)
        {
            f[iy * (ihor + 1) + ix] = -2.0 / (((ix - ihor / 2.0) * (ix - ihor / 2.0) + 0.01) * ((iy - ivert / 2.0) * (iy - ivert / 2.0) + 0.01));
            u[iy * (ihor + 1) + ix] = 0;
        }
    }

    // boundary conditions
    for (int ix = 0; ix <= ihor; ++ix) {
        u[0 * (ihor + 1) + ix] = 0;
        u[ivert * (ihor + 1) + ix] = -ix * (ihor - ix) * (1.0 / ihor);
    }
    for (int iy = 0; iy <= ivert; ++iy) {
        u[iy * (ihor + 1) + 0] = 0;
        u[iy * (ihor + 1) + ihor] = iy* (ivert - iy)* (1.0 / ivert);
    }

    /*for (int iy = 0; iy <= ivert; ++iy) {
        for (int ix = 0; ix <= ihor; ++ix) {
            printf("%-6.2f ", u[iy * (ihor + 1) + ix]);
        }
        printf("\n");
    }
    printf("\n");*/

    // b, right hand vector
    for (int iy = 1; iy < ivert; iy++)
        for (int ix = 1; ix < ihor; ix++)
            b[(iy - 1) * (ihor - 1) + ix - 1] = (1.0 / ihor) * (1.0 / ivert) * f[iy * (ihor + 1) + ix];
    for (int ix = 1; ix < ihor; ix++)
    {
        b[(ix - 1)] -= u[0 * (ihor + 1) + ix];
        b[(ihor - 1) * (ivert - 2) + (ix - 1)] -= u[ivert * (ihor + 1) + ix];
    }
    for (int iy = 1; iy < ivert; iy++)
    {
        b[(ihor - 1) * (iy - 1)] -= u[iy * (ihor + 1) + 0];
        b[(ihor - 1) * iy - 1] -= u[iy * (ihor + 1) + ihor];
    }
    /*for (int i = 0; i < rank; i++)
        printf("b[%d] = %.3f\n", i, b[i]);
    printf("\n");*/

    int nz = 4 * (ihor - 1) - 2 + (ivert - 3) * (5 * (ihor - 1) - 2) + 4 * (ihor - 1) - 2; // top row + (height - 2) middle rows + bottom row contributions

    int* Ap = calloc(rank + 1, sizeof(int));
    int* Ai = calloc(nz, sizeof(int));
    double* Ax = calloc(nz, sizeof(double));

    int col = 0;
    int ix = 0;
    int nzcount = 0;
    Ap[col++] = nzcount;
    nzcount += 3;
    Ap[col++] = nzcount;
    Ai[ix] = 0;
    Ax[ix++] = -4.0;
    Ai[ix] = 1;
    Ax[ix++] = 1.0;
    Ai[ix] = ihor - 1;
    Ax[ix++] = 1.0;
    for (int i = 1; i < ihor - 2; ++i)
    {
        nzcount += 4;
        Ap[col++] = nzcount;
        Ai[ix] = i - 1;
        Ax[ix++] = 1.0;
        Ai[ix] = i;
        Ax[ix++] = -4.0;
        Ai[ix] = i + 1;
        Ax[ix++] = 1.0;
        Ai[ix] = i + ihor - 1;
        Ax[ix++] = 1.0;
    }
    nzcount += 3;
    Ap[col++] = nzcount;
    Ai[ix] = ihor - 3;
    Ax[ix++] = 1.0;
    Ai[ix] = ihor - 2;
    Ax[ix++] = -4.0;
    Ai[ix] = ihor - 2 + ihor - 1;
    Ax[ix++] = 1.0;

    for (int iv = 1; iv < ivert - 2; ++iv)
    {
        nzcount += 4;
        Ap[col++] = nzcount;
        Ai[ix] = (iv - 1) * (ihor - 1);
        Ax[ix++] = 1.0;
        Ai[ix] = iv * (ihor - 1);
        Ax[ix++] = -4.0;
        Ai[ix] = iv * (ihor - 1) + 1;
        Ax[ix++] = 1.0;
        Ai[ix] = iv * (ihor - 1) + ihor - 1;
        Ax[ix++] = 1.0;
        for (int i = 1; i < ihor - 2; ++i)
        {
            nzcount += 5;
            Ap[col++] = nzcount;
            Ai[ix] = (iv - 1) * (ihor - 1) + i;
            Ax[ix++] = 1.0;
            Ai[ix] = iv * (ihor - 1) + i - 1;
            Ax[ix++] = 1.0;
            Ai[ix] = iv * (ihor - 1) + i;
            Ax[ix++] = -4.0;
            Ai[ix] = iv * (ihor - 1) + i + 1;
            Ax[ix++] = 1.0;
            Ai[ix] = iv * (ihor - 1) + (ihor - 1) + i;
            Ax[ix++] = 1.0;
        }
        nzcount += 4;
        Ap[col++] = nzcount;
        Ai[ix] = (iv - 1) * (ihor - 1) + ihor - 2;
        Ax[ix++] = 1.0;
        Ai[ix] = iv * (ihor - 1) + ihor - 3;
        Ax[ix++] = 1.0;
        Ai[ix] = iv * (ihor - 1) + ihor - 2;
        Ax[ix++] = -4.0;
        Ai[ix] = iv * (ihor - 1) + (ihor - 1) + ihor - 2;
        Ax[ix++] = 1.0;
    }

    nzcount += 3;
    Ap[col++] = nzcount;
    Ai[ix] = (ivert - 2) * (ihor - 1) - (ihor - 1);
    Ax[ix++] = 1.0;
    Ai[ix] = (ivert - 2) * (ihor - 1);
    Ax[ix++] = -4.0;
    Ai[ix] = (ivert - 2) * (ihor - 1) + 1;
    Ax[ix++] = 1.0;
    for (int i = 1; i < ihor - 2; ++i)
    {
        nzcount += 4;
        Ap[col++] = nzcount;
        Ai[ix] = (ivert - 2) * (ihor - 1) - (ihor - 1) + i;
        Ax[ix++] = 1.0;
        Ai[ix] = (ivert - 2) * (ihor - 1) + i - 1;
        Ax[ix++] = 1.0;
        Ai[ix] = (ivert - 2) * (ihor - 1) + i;
        Ax[ix++] = -4.0;
        Ai[ix] = (ivert - 2) * (ihor - 1) + i + 1;
        Ax[ix++] = 1.0;
    }
    nzcount += 3;
    Ap[col++] = nzcount;
    Ai[ix] = (ivert - 2) * (ihor - 1) - 1;
    Ax[ix++] = 1.0;
    Ai[ix] = (ivert - 2) * (ihor - 1) + ihor - 3;
    Ax[ix++] = 1.0;
    Ai[ix] = (ivert - 2) * (ihor - 1) + ihor - 2;
    Ax[ix++] = -4.0;

    /*for (int i = 0; i < nzcount; ++i)
    {
        printf("%d: %d, %.0f; ", i, Ai[i], Ax[i]);
    }
    printf("\n");
    for (int i = 0; i < rank + 1; ++i)
    {
        printf("%d: %d; ", i, Ap[i]);
    }
    printf("\n");*/

    double* D = calloc(rank, sizeof(double));
    double* Y = calloc(rank, sizeof(double));
    int* Lp = calloc(rank + 1, sizeof(int));
    int* Parent = calloc(rank, sizeof(int));
    int* Lnz = calloc(rank, sizeof(int));
    int* Pattern = calloc(rank, sizeof(int));
    int d;

    clock_t startSolver = clock();
    ldl_symbolic(rank, Ap, Ai, Lp, Parent, Lnz);
    //printf("Nonzeros in L, excluding diagonal: %d\n", Lp[rank]) ;
    double* Lx = calloc(Lp[rank], sizeof(double));
    int* Li = calloc(Lp[rank], sizeof(int));

    d = ldl_numeric(rank, Ap, Ai, Ax, Lp, Parent, Lnz, Li, Lx, D, Y, Pattern);

    /*for (int i = 0 ; i < Lp[rank] ; i++)
        printf("Li[%d] = %d, Lx[%d] = %.2f;   ", i, Li[i], i, Lx[i]);
    printf("\n");
    printf("\n");
    for (int i = 0 ; i < rank+1 ; i++)
        printf("Lp[%d] = %d;   ", i, Lp[i]);
    printf("\n");
    printf("\n");
    for (int i = 0 ; i < rank ; i++)
        printf("D[%d] = %.2f;   ", i, D[i]);
    printf("\n");
    printf("\n");*/
    
    // copy a[,rank] to b
    if (d == rank)
    {
        ldl_lsolve(rank, b, Lp, Li, Lx);
        ldl_dsolve(rank, b, D);
        ldl_ltsolve(rank, b, Lp, Li, Lx);
        /*for (int i = 0 ; i < rank ; i++)
            printf("x[%d] = %.3f\n", i, b[i]);
        printf("\n");*/
    }
    else
    {
        printf("ldl_numeric failed, D (%d, %d) is zero", d, d);
    }

    for (int iy = 1; iy < ivert; iy++)
        for (int ix = 1; ix < ihor; ix++)
            u[iy * (ihor + 1) + ix] = b[(iy - 1) * (ihor - 1) + (ix - 1)];

    clock_t finishSolver = clock();

    for (int iy = 0; iy <= ivert; ++iy) {
        for (int ix = 0; ix <= ihor; ++ix) {
            printf("%-6.2f ", u[iy * (ihor + 1) + ix]);
        }
        printf("\n");
    }
    printf("\n");

    printf("Solved in %d ms\n", finishSolver - startSolver);
}