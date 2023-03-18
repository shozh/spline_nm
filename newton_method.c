#include "tinyexpr.h"
#include <stdio.h>
#include <math.h>
#include "matrix_lib.c"

const double eps = 1e-3;
int err;

int main() {
    int n;
    scanf("%d", &n);
    char input[n][100];
    for (int i = 0; i < n; i++) {
        scanf("%s", input[i]);
    }

    double x[n];
    te_variable vars[n];
    char xi[n][100];
    for (int i = 0; i < n; i++) {
        sprintf(xi[i], "x%d", i+1);
        vars[i] = (te_variable){xi[i], &x[i]};
    }

    #��������� �����������
    for (int i = 0; i < n; i++)
        x[i] = 0.5;

    while (1) {

        #�������� �������-������� �� ������� �������-�����������
        double f[n];
        for (int i = 0; i < n; i++)
            f[i] = te_eval(te_compile(input[i], vars, n, &err));

        #��������, ��� �������� ������-������� ������ eps
        int b = 1;
        for (int i = 0; i < n; i++)
            if (abs(f[i]) > eps) {
                b = 0;
                break;
            }

        if (b) {
            for (int i = 0; i < n; i++)
                printf("%.2f\n", x[i]);
            break;
        }

        #����� �������
        double** J;
        allocate_memory(&J, n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                x[j] += eps;
                double epsf = te_eval(te_compile(input[i], vars, n, &err));
                J[i][j] = (epsf - f[i]) / eps;
                x[j] -= eps;
            }
        }

        #�������� ����� �������
        double** iJ = get_transposed_matrix(get_matrix_of_cofactors(J, n), n, n);
        double det = get_determinant(J, n);
        multiply_by_number(&iJ, n, n, 1/det);

        #������� �� �������-�������
        double** F;
        allocate_memory(&F, n, 1);
        for (int i = 0; i < n; i++)
            J[i][0] = f[i];

        #��������� ��������� �������� ����� ������� � ������� �� �������-�������
        double** second;
        second = multiply_by_matrix(iJ, n, n, F, n, 1);

        #������������ ������� ������ �������-�����������
        for (int i = 0; i < n; i++) {
            x[i] -= second[i][0];
        }
    }
}
