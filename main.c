#include <stdio.h>
#include <stdlib.h>
#include <iso646.h>
#include <math.h>
#define M_PI 3.141592653589793
#define M_2PI 2. * M_PI

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

double max(double x, double y) {
    return x > y ? x : y;
}

double min(double x, double y) {
    return x < y ? x : y;
}

void print_array(double* A, int N) {
    printf("coefs of array: ");
    for (int i = 0; i < N; i++)
        printf("%.3lf ", A[i]);
    printf("\n");
}

int solve_cubic_eq(double* x,double a,double b,double c) {
    double q,r,r2,q3;
    q=(a*a-3.*b)/9.; r=(a*(2.*a*a-9.*b)+27.*c)/54.;
    r2=r*r; q3=q*q*q;
    if(r2<q3) {
        double t=acos(r/sqrt(q3));
        a/=3.; q=-2.*sqrt(q);
        x[0]=q*cos(t/3.)-a;
        x[1]=q*cos((t+M_2PI)/3.)-a;
        x[2]=q*cos((t-M_2PI)/3.)-a;
        return(3);
    }
    else {
        double aa,bb;
        if(r<=0.) r=-r;
        aa=-pow(r+sqrt(r2-q3),1./3.);
        if(aa!=0.) bb=q/aa;
        else bb=0.;
        a/=3.; q=aa+bb; r=aa-bb;
        x[0]=q-a;
        x[1]=(-0.5)*q-a;
        x[2]=(sqrt(3.)*0.5)*fabs(r);
        if(x[2]==0.) return(2);
        return(1);
    }
}

int solve_sq_eq(double* ret, double a, double b, double c) {
    double D = b * b - 4 * a * c;
    if (D < 0)
        return (0);
    else {
        double x1 = (-b + sqrt(D))/(2 * a);
        double x2 = (-b - sqrt(D))/(2 * a);
        ret[0] = x1;
        ret[1] = x2;
        return (2);
    }
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

double intersection_of_two_splines(double* D1, double* E1, int N, double* D2, double* E2, int M) {
    if ((D2[0] <= D1[N - 1]) or (D1[0] <= D2[M - 1])) {
        for (int i = 1; i < N; i++) {
            for (int j = 1; j < M; j++) {
                double left_border = max(D1[i - 1], D2[j - 1]);
                double right_border = min(D1[i], D2[j]);
                if (left_border <= right_border) {
                    double *a = get_coefs(i, D1, E1, N);
                    double *b = get_coefs(j, D2, E2, M);
                    double A = a[3] - b[3];
                    double B = a[2] - b[2];
                    double C = a[1] - b[1];
                    double D = a[0] - b[0];
                    if (A != 0) {
                        double *R = (double*)malloc(3 * sizeof(double));
                        int nums = solve_cubic_eq(R, B/A, C/A, D/A);
                            for (int k = 0; k < nums; k++) {
                                if ((left_border <= R[k]) and (R[k] <= right_border)) {
                                    print_array(a, 4);
                                    print_array(b, 4);
                                    printf("%.3lf, %.3lf, %.3lf, %.3lf\n", A, B, C, D);
                                    printf("%.3lf, %.3lf\n", D1[i - 1], D1[i]);
                                    printf("%.3lf, %.3lf\n", D2[j - 1], D2[j]);
                                    return R[k];
                                }
                            }
                        free(R);
                    } else if (B != 0) {
                        double *R = (double*)malloc(2 * sizeof(double));
                        int nums = solve_sq_eq(R, B, C, D);
                        if (nums == 0)
                            continue;
                        for (int k = 0; k < 2; k++) {
                            if (left_border <= R[k] and R[k] <= right_border) {
                                return R[k];
                            }
                        }
                        free(R);
                    } else if (C != 0) {
                        double v = C / D;
                        if (left_border <= v and v <= right_border)
                            return v;
                    }
                }
            }
        }
    }
    return 404;
}