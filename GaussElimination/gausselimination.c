/* TERMS OF USE
 * This source code is subject to the terms of the MIT License.
 * Copyright(c) 2026 Vladimir Vasilich Tregub
*/
#include <math.h>
#include <stdio.h>
int gausselim(double* AB, int rows, int cols)
{
    // convert matrix to a triangular matrix using Gaussian elimination
    // with partial pivoting
    for (int i = 0; i < rows; i++)
    {
        double maxpival = fabs(AB[i * cols + i]);
        int jpiv = i;
        for (int j = i + 1; j < rows; j++)
            if (maxpival < fabs(AB[j * cols + i]))
            {
                maxpival = fabs(AB[j * cols + i]);
                jpiv = j;
            }

        if (jpiv != i)
            for (int k = i; k < cols; k++)
            {
                double tmp = AB[i * cols + k];
                AB[i * cols + k] = AB[jpiv, k];
                AB[jpiv * cols + k] = tmp;
            }

        if (AB[i * cols + i] != 0)
        {
            for (int j = i + 1; j < rows; j++)
            {
                double rowScale = AB[j * cols + i] / AB[i * cols + i];
                for (int k = 0; k < cols; k++)
                {
                    if (k <= i)
                        AB[j * cols + k] = 0;
                    else
                        AB[j * cols + k] -= AB[i * cols + k] * rowScale;
                }
            }
        }
    }

    return 0;
}

int solvelinsys(double* AB, int rows, int cols, double* x_arr)
{
    int res = gausselim(AB, cols, rows);
    /*printf("\n");
    for (int iy = 0; iy < cols; ++iy) {
        for (int ix = 0; ix < rows; ++ix) {
            printf("%-6.2f ", AB[iy * rows + ix]);
        }
        printf("\n");
    }*/
    // get the solution vector x by backsubstitution
    //double* x_arr = (double*) malloc(rows * sizeof(double));
    for (int i = cols - 1; i >= 0; i--) {
        double x = AB[i * rows + rows - 1];
        for (int j = rows - 2; j > i; j--)
            x -= AB[i * rows + j] * x_arr[j];
        //if (AB[i * cols + i] != 0)
        x_arr[i] = (x / AB[i * rows + i]);
    }
    return 0;
}
