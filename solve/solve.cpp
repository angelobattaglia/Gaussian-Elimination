#include <stdio.h>
#include <math.h>

double *solve_backward(double **U, double *b, double *x, int n)
{
    /*
 Risolve il sistema Ux=b con U matrice n x n triangolare superiore non singolare
 usando la backward substitution.
 */
    int i,
        j;
    double tmp;
    x[n - 1] = b[n - 1] / U[n - 1][n - 1];
    for (i = n - 2; i >= 0; i--)
    {
        tmp = 0.0;
        for (j = i + 1; j < n; j++)
            tmp += U[i][j] * x[j];
        x[i] = (b[i] - tmp) / U[i][i];
    }
    return x;
}

double *solve_forward(double **L, double *b, double *x, int n)
{
    /*
 Risolve il sistema Lx=b con L matrice n x n triangolare inferiore non singolare
 usando la forward substitution.
 */
    int i,
        j;
    double tmp;
    x[0] = b[0] / L[0][0];
    for (i = 1; i < n; i++)
    {
        tmp = 0.0;
        for (j = 0; j < i; j++)
            tmp += L[i][j] * x[j];
        x[i] = (b[i] - tmp) / L[i][i];
    }

    return x;
}
