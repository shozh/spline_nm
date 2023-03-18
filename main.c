#include <stdio.h>
#include <stdlib.h>

double runge(double x) {
    return 1 / (1 + x * x);
}

double* tridiagonal_matrix_algorithm(double* a, double* b, double* c, double* d, int n) {
    double* x = (double*)malloc(n * sizeof(double*));

    for (int i = 1; i < n; i++) {
        double coef = a[i] / b[i - 1];
        b[i] -= coef * c[i-1];
        d[i] -= coef * d[i-1];
    }

    for (int i = n - 1; i > 0; i--) {
        double coef = c[i-1]/b[i];
        d[i - 1] -= d[i] * coef;
    }

    for (int i = 0; i < n; i++)
        x[i] = d[i]/b[i];

    return x;
}

double* get_gamma(double* D, double* E, int N) {


}

int main() {

    int N = 7;

    double* D = (double*)malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++)
        D[i] = i;

    double* E = (double*)malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++)
        E[i] = runge(i);



    return 0;
}
