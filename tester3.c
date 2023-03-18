#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double* solve_cubic_eq(double a, double b, double c, double d) {
    double f = ((3 * c / a) - (b*b/(a * a)))/3;
    double g = ((2 * b * b * b)/(a * a * a) - (9 * b * c)/(a * a) + (27 * d)/a)/27;
    double h = (g * g / 4) + (f * f * f)/27;
    double i = pow(g * g / 4 - h, 0.5);
    printf("h %lf\n", h);
    printf("f %lf\n", f);
    printf("f %lf\n", g);

    if (h <= 0) {
        double j = pow(i, 0.333333);
        double k = acos(-g / (2 * i));
        double l = -j;
        double m = cos(k/3);
        double n = pow(3, 0.5) * sin(k/3);
        double p = -(b/(3*a));
        double* ret = (double*)malloc(3 * sizeof(double));

#ifdef Te
        printf("j: %.2lf\n", j);
        printf("k: %.2lf\n", k);
        printf("m: %.2lf\n", m);
        printf("n: %.2lf\n", n);
        printf("p: %.2lf\n", p);

#endif
        ret[0] = 2 * j * cos(k/3)-(b/(3*a));
        ret[1] = l * (m + n) + p;
        ret[2] = l * (m - n) + p;
        return ret;
    }
    if (h > 0) {
        double R = -g/2 + sqrt(h);
        double S = pow(R, 0.3333333);
        double T = -(g/2) - sqrt(h);
        double U = cbrt(T);
#ifdef Te
        printf("R: %.2lf\n", R);
        printf("S: %.2lf\n", S);
        printf("T: %.2lf\n", T);
        printf("U: %.2lf\n", U);
        printf("h: %.2lf\n", h);
#endif

        double* ret = (double*)malloc(3 * sizeof(double));
        double xx = (S + U) - (b/(3 * a));
        ret[0] = xx;
        ret[1] = xx;
        ret[2] = xx;
        //       printf("helllo\n");
        return ret;
    }
    if (f == 0 and g == 0 and h == 0) {

        double* ret = (double*)malloc(3 * sizeof(double));
        double x = -pow(d/a, (double)1/3);
        ret[0] = x;
        ret[1] = x;
        ret[2] = x;
        printf("helllo");
        return ret;
    }
}


int main(void) {


}