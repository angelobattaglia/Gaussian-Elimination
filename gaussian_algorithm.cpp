#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double *solve_backward(double **, double *, double *, int);
double *solve_forward(double **, double *, double *, int);

int main()
{
    int ndim, i, j, k, ind_pivot;
    double **A = NULL;
    double *p_tmp, *x, *y, *b, *diag;
    int *pivot = NULL;
    double el_pivot, tmp;

    // Si legge da standard input la dimensione della matrice

    printf("Inserire la dimensione della matrice di Hilbert:\n");
    scanf("%d", &ndim);

    // Si alloca dinamicamente la matrice e il vettore pivot

    A = (double **)malloc(ndim * sizeof(double *));
    pivot = (int *)malloc((ndim - 1) * sizeof(double));
    x = (double *)malloc(ndim * sizeof(double));
    y = (double *)malloc(ndim * sizeof(double));
    b = (double *)malloc(ndim * sizeof(double));
    diag = (double *)malloc(ndim * sizeof(double));

    for (i = 0; i < ndim; i++)
    {
        A[i] = (double *)malloc(ndim * sizeof(double));
    }

    // Si stampa la matrice e si costruisce il vettore b in maniera tale che la soluzione
    //x sia un vettore con
    // tutti gli elementi uguali a 1
    printf("\nLa matrice di Hilbert e':\n\n");
    for (i = 0; i < ndim; i++)
    {
        b[i] = 0;
        for (j = 0; j < ndim; j++)
        {
            A[i][j] = 1.0 / (i + j + 1);
            b[i] += A[i][j];
            printf("%f ", A[i][j]);
        }
        printf("\n");
    }
    // Si effettua l'eliminazione gaussiana con pivoting parziale
    for (k = 0; k < ndim - 1; k++)
    {
        // Si ricerca l'elemento pivot
        el_pivot = fabs(A[k][k]);
        ind_pivot = k;
        for (i = k + 1; i < ndim; i++)
            if (fabs(A[i][k]) > el_pivot)
            {
                el_pivot = fabs(A[i][k]);
                ind_pivot = i;
            }
        pivot[k] = ind_pivot;
        // Si scambiano le righe se necessario
        if (pivot[k] != k)
        {
            p_tmp = A[k];
            A[k] = A[pivot[k]];
            A[pivot[k]] = p_tmp;
        }
        // Si procede con l'algoritmo di eliminazione Gaussiana
        for (i = k + 1; i < ndim; i++)
            if (A[k][k] != 0.0)
            {
                A[i][k] = A[i][k] / A[k][k]; // si salva il modificatore
                //direttamente nella parte triangolare inferiore della matrice
                for (j = k + 1; j < ndim; j++)
                    A[i][j] -= A[i][k] * A[k][j];
            }
            else
            {
                printf("\nErrore: la matrice inserita e' singolare\n");
                return -1;
            }
    }
    // Si ristampa la matrice e il vettore pivot dopo la modifica
    printf("\nLa matrice modificata e':\n\n");
    for (i = 0; i < ndim; i++)
    {
        for (j = 0; j < ndim; j++)
        {
            printf("%f ", A[i][j]);
        }
        printf("\n");
    }
    printf("\nIl vettore pivot degli scambi e':\n\n");
    for (i = 0; i < ndim - 1; i++)
        printf("%d ", pivot[i]);
    /*
 Si procede con la risoluzione del sistema lineare Ax=b sfruttando la fattorizzazione
PA=LU
 che porta a risolvere due sistemi triangolari Ly=Pb e Ux=y
 */
    // Si salva la diagonale di A, la si imposta uguale a tutti 1 e si permuta b secondo
    // le informazioni del vettore pivot
    for (i = 0; i < ndim; i++)
    {
        diag[i] = A[i][i];
        A[i][i] = 1.0;
        if (i != ndim - 1 && pivot[i] != i)
        {
            tmp = b[i];
            b[i] = b[pivot[i]];
            b[pivot[i]] = tmp;
        }
    }
    y = solve_forward(A, b, y, ndim);
    // Si reimposta la diagonale di A
    for (i = 0; i < ndim; i++)
    {
        A[i][i] = diag[i];
    }
    x = solve_backward(A, y, x, ndim);
    // Si stampa il vettore soluzione x (come vettore riga)
    printf("\nLa soluzione del sistema e' (rappresentata come vettore riga):\n");
    for (i = 0; i < ndim; i++)
    {
        printf("%f ", x[i]);
    }
    free(A);
    free(x);
    free(y);
    free(b);
    free(diag);
    free(pivot);
    return 0;
}
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
