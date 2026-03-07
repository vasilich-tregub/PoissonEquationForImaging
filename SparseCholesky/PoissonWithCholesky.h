// euclideangraphics.h : Include file for standard system include files,
// or project specific include files.

#pragma once

//#include <iostream>
void ldl_symbolic(int n, int* Ap, int* Ai, int* Lp, int* Parent, int* Lnz);
int ldl_numeric(int n, int* Ap, int* Ai, double* Ax, int* Lp, int* Parent, int* Lnz, int* Li, double* Lx, double* D, double* Y, int* Pattern);
void ldl_lsolve(int n, double* X, int* Lp, int* Li, double* Lx);
void ldl_dsolve(int n, double* X, double* D);
void ldl_ltsolve(int n, double* X, int* Lp, int* Li, double* Lx);
