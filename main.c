#include <stdio.h>
#include <stdlib.h>
#include <iso646.h>
#include "tester.c"
#include <math.h>

double runge(double x) {
    return 1 / (1 + x * x);
}

double* tridiagonal_matrix_algorithm(double* a, double* b, double* c, double* d, int n) {
    double* x = (double*)malloc(n * sizeof(double));

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

double* get_gamma(const double* D, const double* E, int N) {
    double *h = (double*)malloc(N * sizeof(double));
    double *a = (double*)malloc((N - 2) * sizeof(double));
    double *b = (double*)malloc((N - 2) * sizeof(double));
    double *c = (double*)malloc((N - 2) * sizeof(double));
    double *d = (double*)malloc((N - 2) * sizeof(double));

    h[0] = 0;
    for (int i = 1; i < N; i++)
        h[i] = D[i] - D[i - 1];

    a[0] = 0;
    for (int i = 2; i < N - 1; i++)
        a[i - 1] = h[i]/6;

    for (int i = 1; i < N - 1; i++)
        b[i - 1] = (h[i] + h[i + 1])/3;

    for (int i = 1; i < N - 2; i++)
        c[i - 1] = h[i + 1]/6;
    c[N - 3] = 0;

    for (int i = 1; i < N - 1; i++)
        d[i - 1] = (E[i + 1] - E[i])/h[i + 1] - (E[i] - E[i - 1])/h[i];

    double* gamma = tridiagonal_matrix_algorithm(a, b, c, d, N - 2);
    double* gamma_extend = (double*)malloc(N * sizeof(double));
    gamma_extend[0] = 0;
    gamma_extend[N - 1] = 0;
    for (int i = 1; i < N - 1; i++)
        gamma_extend[i] = gamma[i-1];

    free(h);
    free(a);
    free(b);
    free(c);
    free(d);
    free(gamma);
    return gamma_extend;
}

double* get_coefs(int i, double* D, double* E, int N) {
    double* gamma = get_gamma(D, E, N);
    double *h = (double*)malloc(N * sizeof(double));
    h[0] = 0;
    for (int j = 1; j < N; j++)
        h[j] = D[j] - D[j - 1];

    double* coefs = (double*)malloc(4 * sizeof(double));
    coefs[3] = (gamma[i] - gamma[i-1])/(6*h[i]);
    coefs[2] = (gamma[i-1] * D[i] - gamma[i] * D[i - 1])/(2 * h[i]);
    coefs[1] = (gamma[i] * D[i - 1] * D[i - 1] - gamma[i - 1] * D[i] * D[i])/(2 * h[i]) + (E[i] - E[i - 1])/h[i] - h[i] * (gamma[i] - gamma[i-1]) / 6;
    coefs[0] = (gamma[i-1] * D[i]*D[i]*D[i] - gamma[i] * D[i-1]*D[i-1]*D[i-1])/(6*h[i]) + (E[i-1] * D[i] - E[i] * D[i - 1])/h[i] + h[i]*(gamma[i] * D[i - 1] - gamma[i - 1] * D[i]) / 6;
    free(gamma);
    free(h);
    return coefs;
}

double P(int i, double x, double* D, double* E, int N) {
    double* a = get_coefs(i, D, E, N);
    return a[3] * x * x * x + a[2] * x * x + a[1] * x + a[0];
}

double Spline(double x, double* D, double* E, int N) {
    for (int i = 1; i < N + 1; i++) {
        if ( (D[i - 1] <= x) and (x <= D[i]))
            return P(i, x, D, E, N);
    }
}

int main() {

    const int N = 7;

    double* D = (double*)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++)
        D[i] = i;

    double* E = (double*)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++)
        E[i] = runge(i);

    for (int i = 1; i < N; i++) {
        double* a = get_coefs(i, D, E, N);
        printf("P_%d(x) = %.3lfx^3 + %.3lfx^2 + %.3lfx + %.3lf, if %.3lf<=x<=%.3lf\n", i, a[3], a[2], a[1], a[0], D[i-1], D[i]);
    }

    while (1) {
        printf("Input x_0 within D (Input 0 to exit)\n");
        double x;
        scanf("%lf", &x);
        if (!x)
            break;
        printf("Runge(%.3lf) = %.3lf\n", x,  runge(x));
        printf("Spline(%.3lf) = %.3lf\n", x, Spline(x, D, E, N));
        printf("Delta = %.3lf\n", fabs(Spline(x, D, E, N) - runge(x)));
    }

    return 0;
}
