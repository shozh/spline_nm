#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iso646.h>
#include "main.c"


double runge(double x) {
    return 1 / (1 + x * x);
}



void Runge_starter() {

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


