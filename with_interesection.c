#include <stdio.h>
#include <stdlib.h>
#include <iso646.h>
#include <math.h>
#include "tester.c"
#define Te 1

double h(double x) {
    return (1/(5 - tan(x)));
}

double t(double x) {
    return (1/(1+x));
}

double f(double x) {
    return 1/(1 + x*x);
}

double g(double x) {
    return 1/(1 + exp(-x));
}

double* solve_cubic_eq(double a, double b, double c, double d) {
    double f = ((3 * c / a) - (b*b/(a * a)))/3;
    double g = ((2 * b * b * b)/(a * a * a) - (9 * b * c)/(a * a) + (27 * d)/a)/27;
    double h = (g * g / 4) + (f * f * f)/27;
    double i = sqrt(g * g / 4 - h);
    double* ret = (double*)malloc(3 * sizeof(double));

    // all 3 roots are real
    if (h <= 0) {
        double j = cbrt(i);
        double k = acos(-g / (2 * i));
        double l = -j;
        double m = cos(k/3);
        double n = sqrt(3) * sin(k/3);
        double p = -(b/(3*a));
        ret[0] = 2 * j * cos(k/3)-(b/(3*a));
        ret[1] = l * (m + n) + p;
        ret[2] = l * (m - n) + p;
        return ret;
    }
    // only 1 root is real
    if (h > 0) {
        double R = -g/2 + sqrt(h);
        double S = cbrt(R);
        double T = -(g/2) - sqrt(h);
        double U = cbrt(T);
        double xx = (S + U) - (b/(3 * a));
        ret[0] = xx;
        ret[1] = xx;
        ret[2] = xx;
        return ret;
    }
    // when all 3 roots are real and equal
}

double* solve_sq_eq(double a, double b, double c) {
    double D = b * b - 4 * a * c;
    double* ret = (double*)malloc(2 * sizeof(double));
    if (D < 0) {
        ret[0] = 55.5;
        ret[1] = 55.5;
        return ret;
    }
    else {
        double x1 = (-b + sqrt(D))/(2 * a);
        double x2 = (-b - sqrt(D))/(2 * a);
        ret[0] = x1;
        ret[1] = x2;
        return ret;
    }
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

typedef struct {
    double begin;
    double end;
} segment;

typedef struct {
    double x;
} polynom;

double max(double x, double y) {
    return x > y ? x : y;
}

double min(double x, double y) {
    return x < y ? x : y;
}

int Program() {

    const int N = 7;
    const int M = 10;
    double *D1 = (double *) malloc(N * sizeof(double));
    double *E1 = (double *) malloc(N * sizeof(double));
    double *D2 = (double *) malloc(M * sizeof(double));
    double *E2 = (double *) malloc(M * sizeof(double));

    for (int i = 0, j = -2; i < N; i++, j++)
        D1[i] = j;

    for (int i = 0; i < N; i++)
        E1[i] = f(D1[i]);

    for (int i = 0; i < M; i++)
        D2[i] = i;

    for (int i = 0; i < M; i++)
        E2[i] = g(D2[i]);

    if ((D2[0] <= D1[N - 1]) or (D1[0] <= D2[N - 1])) {
        for (int i = 1; i < N; i++) {
            for (int j = 1; j < M; j++) {
                double left_border = max(D1[i - 1], D2[j - 1]);
                double right_border = min(D1[i], D2[j]);
                if (left_border <= right_border) {
                    double *a = get_coefs(i, D1, E1, N);
                    double *b = get_coefs(j, D2, E2, M);
                    double A = a[0] - b[0];
                    double B = a[1] - b[1];
                    double C = a[2] - b[2];
                    double D = a[3] - b[3];
                    if (A != 0) {
                        double *R = solve_cubic_eq(A, B, C, D);
                        printf("coefs: %lf, %lf, %lf, %lf\n", A, B, C, D);
                        printf("left: %lf, right: %lf\n", left_border, right_border);
                        printf("D1[i-1]: %lf, D1[i], %lf\n", D1[i-1], D1[i]);
                        printf("D2[i-1]: %lf, D2[i], %lf\n", D2[i-1], D2[i]);
                        print_roots(R, 3);
                        puts("");
                        for (int k = 0; k < 3; k++) {
                            if ((left_border <= R[k]) and (R[k] <= right_border)) {
                                print_array(a, 4);
                                print_array(b, 4);
                                printf("coefs: %lf, %lf, %lf, %lf\n", A, B, C, D);
                                printf("left: %lf, right: %lf\n", left_border, right_border);
                                printf("D1[i-1]: %lf, D1[i], %lf\n", D1[i-1], D1[i]);
                                printf("D2[i-1]: %lf, D2[i], %lf\n", D2[i-1], D2[i]);
                                printf("(%lf, %lf)\n", R[k], Spline(R[k], D1, E1, N));
                                printf("(%lf, %lf)", R[k], Spline(R[k], D2, E2, M));
                                return 0;
                            }
                        }
                    } else if (B != 0) {
                        double* R = solve_sq_eq(B, C, D);
                        if (R[0] == 55.5)
                            continue;
                        for (int k = 0; k < 2; k++) {
                            if (left_border <= R[k] and R[k] <= right_border) {
                                printf("(%lf, %lf)", R[k], Spline(R[k], D1, E1, N));
                                return 0;
                            }
                        }
                    } else if (C != 0) {
                        double v = C / D;
                        if (left_border <= v and v <= right_border) {
                            printf("(%lf, %lf)", v, Spline(v, D1, E1, N));
                            return 0;
                        }
                    }
                }
            }
        }
        printf("There doesn't exist any intersection between these two splines within borders");
    } else {
        printf("There doesn't exist any intersection between these two splines within borders");
    }
}

int main(void) {
    int ret = Program();
    //double* R = solve_sq_eq(1, 0, -16);
    //printf("%.2lf, %.2lf", R[0], R[1]);
   //  double* R = solve_cubic_eq(-0.970066, 1.234188, -0.458725, -0.014917);
   //printf("%.2lf, %.2lf, %.2lf", R[0], R[1], R[2]);

    return 0;
}