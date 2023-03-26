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

double* intersection_of_two_splines(int N, double* D1, double* E1, double* D2, double* E2) {
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
                        double R[3] = {0, 0, 0};
                        int nums = solve_cubic_eq(R, B/A, C/A, D/A);
                        if (nums == 3)
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
                        else if (nums == 1)
                            if ((left_border <= R[0]) and (R[0] <= right_border)) {
                                print_array(a, 4);
                                print_array(b, 4);
                                printf("coefs: %lf, %lf, %lf, %lf\n", A, B, C, D);
                                printf("left: %lf, right: %lf\n", left_border, right_border);
                                printf("D1[i-1]: %lf, D1[i], %lf\n", D1[i-1], D1[i]);
                                printf("D2[i-1]: %lf, D2[i], %lf\n", D2[i-1], D2[i]);
                                printf("(%lf, %lf)\n", R[0], Spline(R[0], D1, E1, N));
                                printf("(%lf, %lf)", R[0], Spline(R[0], D2, E2, M));
                            }

                        free(R);
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
