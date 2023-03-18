#include "tinyexpr.h"
#include <stdio.h>
#include <math.h>

const double eps = 1e-3;
int err;

int main() {
    double x;
    char input[100];
    scanf("%s", input);
    te_variable vars[] = {{"x", &x}};
    x = 2;
    do {
        double f = te_eval(te_compile(input, vars, 1, &err));
        x += eps;
        double epsf= te_eval(te_compile(input, vars, 1, &err));
        double derf = (epsf - f)/eps;
        x -= eps;
        x -= f/derf;
    } while (abs(te_eval(te_compile(input, vars, 1, &err))) > eps);
    printf("%.2f", x);
    return 0;
}