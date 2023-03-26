#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iso646.h>
#include "main.c"

double runge(double x) {
    return 1 / (1 + x * x);
}

double sigmoid(double x) {
    return 1/(1 + exp(-x));
}


int Runge_starter() {

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
        printf("Input number in [%.2lf, %.2lf] (Input 404 to exit)\n", D[0], D[N - 1]);
        double x;
        scanf("%lf", &x);
        if (x == 404)
            break;
        if (x > D[N - 1] or x < D[0]) {
            printf("%.3lf not in [%.2lf, %.2lf]\n", x, D[0], D[N - 1]);
            continue;
        }
        printf("Runge(%.3lf) = %.3lf\n", x,  runge(x));
        printf("Spline(%.3lf) = %.3lf\n", x, Spline(x, D, E, N));
        printf("Error = %.3lf\n", fabs(Spline(x, D, E, N) - runge(x)));
    }
    return 0;
}

int Runge_sigmoid_intersection() {

    const int N = 7;
    const int M = 10;
    double *D1 = (double *) malloc(N * sizeof(double));
    double *E1 = (double *) malloc(N * sizeof(double));
    double *D2 = (double *) malloc(M * sizeof(double));
    double *E2 = (double *) malloc(M * sizeof(double));


    for (int i = 0; i < N; i++)
        D1[i] = i;

    for (int i = 0; i < N; i++)
        E1[i] = runge(i);

    for (int j = 0; j < M; j++)
        D2[j] = j;

    for (int j = 0; j < M; j++)
        E2[j] = sigmoid(j);

    double x = intersection_of_two_splines(D1, E1, N, D2, E2, N);
    if (x == 404)
        printf("There doesn't exist any intersection between these two splines within borders");
    else {
        printf("First spline: (%lf, %lf)\n", x, Spline(x, D1, E1, N));
        printf("Second spline: (%lf, %lf)", x, Spline(x, D2, E2, M));
    }
    return 0;
}

int main() {
    Runge_sigmoid_intersection();
    //Runge_starter();
    return 0;
}


