//------------------------------------------------------------------------------
// LDL/Source/ldl.c: sparse LDL' factorization
//------------------------------------------------------------------------------

// LDL, Copyright (c) 2005-2022 by Timothy A. Davis. All Rights Reserved.
// SPDX-License-Identifier: LGPL-2.1+

// A minified version of 
// https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/dev/LDL/Source/ldl.c
// See CREDIT.md of this folder

#include <malloc.h>

void ldl_symbolic(int n, int* Ap, int* Ai, int* Lp, int* Parent, int* Lnz)
{
    int i, k, p, p2;
    int* Flag = calloc(n, sizeof(int));
    for (k = 0; k < n; k++)
    {
        // L(k,:) pattern: all nodes reachable in etree from nz in A(0:k-1,k)
        Parent[k] = -1; // parent of k is not yet known
        Flag[k] = k; // mark node k as visited 
        Lnz[k] = 0; // count of nonzeros in column k of L 
        p2 = Ap[k + 1];
        for (p = Ap[k]; p < p2; p++)
        {
            // A (i,k) is nonzero (original or permuted A) 
            i = Ai[p];
            if (i < k)
            {
                // follow path from i to root of etree, stop at flagged node
                for (; Flag[i] != k; i = Parent[i])
                {
                    // find parent of i if not yet determined 
                    if (Parent[i] == -1)
                        Parent[i] = k;
                    Lnz[i]++; // L (k,i) is nonzero 
                    Flag[i] = k;
                }
            }
        }
    }
    // construct Lp index array from Lnz column counts 
    Lp[0] = 0;
    for (k = 0; k < n; k++)
    {
        Lp[k + 1] = Lp[k] + Lnz[k];
    }
}

// returns n if successful, k if D (k,k) is zero
int ldl_numeric(int n, int* Ap, int* Ai, double* Ax, int* Lp, int* Parent, int* Lnz, int* Li, double* Lx, double* D, double* Y, int* Pattern)
{
    double yi, l_ki;
    int i, k, p, p2, len, top;
    int* Flag = calloc(n, sizeof(int));
    for (k = 0; k < n; k++)
    {
        // compute nonzero Pattern of kth row of L, in topological order 
        Y[k] = 0.0; // Y(0:k) is now all zero 
        top = n; // stack for pattern is empty 
        Flag[k] = k; // mark node k as visited
        Lnz[k] = 0; // count of nonzeros in column k of L
        p2 = Ap[k + 1];
        for (p = Ap[k]; p < p2; p++)
        {
            i = Ai[p]; // get A(i,k)
            if (i <= k)
            {
                Y[i] += Ax[p]; // scatter A(i,k) into Y (sum duplicates)
                for (len = 0; Flag[i] != k; i = Parent[i])
                {
                    Pattern[len++] = i; // L(k,i) is nonzero
                    Flag[i] = k; // mark i as visited
                }
                while (len > 0)
                    Pattern[--top] = Pattern[--len];
            }
        }
        // compute numerical values kth row of L (a sparse triangular solve)
        D[k] = Y[k]; // get D(k,k) and clear Y(k) 
        Y[k] = 0.0;
        for (; top < n; top++)
        {
            i = Pattern[top]; // Pattern [top:n-1] is pattern of L(:,k) 
            yi = Y[i]; // get and clear Y(i) 
            Y[i] = 0.0;
            p2 = Lp[i] + Lnz[i];
            for (p = Lp[i]; p < p2; p++)
            {
                Y[Li[p]] -= Lx[p] * yi;
            }
            l_ki = yi / D[i]; // the nonzero entry L(k,i) 
            D[k] -= l_ki * yi;
            Li[p] = k; // store L(k,i) in column form of L
            Lx[p] = l_ki;
            Lnz[i]++; // increment count of nonzeros in col i 
        }
        if (D[k] == 0.0)
            return (k); // failure, D(k,k) is zero
    }
    return (n); // success, diagonal of D is all nonzero
}

void ldl_lsolve(int n, double* X, int* Lp, int* Li, double* Lx)
{
    int j, p, p2;
    for (j = 0; j < n; j++)
    {
        p2 = Lp[j + 1];
        for (p = Lp[j]; p < p2; p++)
        {
            X[Li[p]] -= Lx[p] * X[j];
        }
    }
}

void ldl_dsolve(int n, double* X, double* D)
{
    int j;
    for (j = 0; j < n; j++)
    {
        X[j] /= D[j];
    }
}

void ldl_ltsolve(int n, double* X, int* Lp, int* Li, double* Lx)
{
    int j, p, p2;
    for (j = n - 1; j >= 0; j--)
    {
        p2 = Lp[j + 1];
        for (p = Lp[j]; p < p2; p++)
        {
            X[j] -= Lx[p] * X[Li[p]];
        }
    }
}
