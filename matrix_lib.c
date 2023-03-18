#include <stdio.h>
#include <stdlib.h>

typedef double Type;//Тип значений в матрице(лучше не менять)

void allocate_memory(Type *** matrix, int row, int col) {//Выделяет память для матрицы
    (*matrix) = (Type **)malloc(sizeof(Type *) * (row));
    for (int i = 0; i < row; i++) {
        (*matrix)[i] = (Type *)malloc(col * sizeof(Type));
    }
}

Type ** get_cofactor(Type ** matrix, int order, int row, int col) {//Алгебраическое дополнение
    Type ** ans;
    allocate_memory(&ans, order - 1, order - 1);
    for (int i = 0; i < order - 1; ++i) {
        for (int j = 0; j < order - 1; ++j) {
            int x, y;
            x = i + (i >= row);
            y = j + (j >= col);
            ans[i][j] = matrix[x][y];
        }
    }
    return ans;
}

Type get_determinant(Type ** matrix, int order) {//Определитель
    if (order == 1) {
        return matrix[0][0];
    }
    Type ans = 0;
    for (int i = 0; i < order; ++i) {
        int pm = ((i & 1) ? -1 : 1);
        ans += (Type)pm * matrix[0][i] * get_determinant(get_cofactor(matrix, order, 0, i), order - 1);
    }
    return ans;
}

Type ** get_matrix_of_cofactors(Type ** matrix, int order) {//Матрица алгебраических дополнений
    Type ** ans;
    allocate_memory(&ans, order, order);
    for (int i = 0; i < order; ++i) {
        for (int j = 0; j < order; ++j) {
            int pm = (((i + j) & 1) ? -1 : 1);
            ans[i][j] = (Type)pm * get_determinant(get_cofactor(matrix, order, i, j), order - 1);
        }
    }
    return ans;
}

Type ** get_transposed_matrix(Type ** matrix, int row, int col) {//Транспонированная матрица
    Type ** ans;
    allocate_memory(&ans, col, row);
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            ans[j][i] = matrix[i][j];
        }
    }
    return ans;
}

Type ** multiply_by_matrix(Type ** matrix_a, int row_a, int col_a,
                           Type ** matrix_b, int row_b, int col_b) {//Умножение двух матриц
    Type ** ans;
    allocate_memory(&ans, row_a, col_b);
    for (int i = 0; i < row_a; ++i) {
        for (int j = 0; j < col_b; ++j) {
            ans[i][j] = 0;
            for (int k = 0; k < col_a; ++k) {
                ans[i][j] += matrix_a[i][k] * matrix_b[k][j];
            }
        }
    }
    return ans;
}

void multiply_by_number(Type *** matrix, int row, int col, Type number) {//Умножение матрицы на число
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            (*matrix)[i][j] *= number;
        }
    }
}

Type ** addition_matrixes(Type ** matrix_a, Type ** matrix_b, int row, int col) {//Сложение двух матриц
    Type ** ans;
    allocate_memory(&ans, row, col);
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            ans[i][j] = matrix_a[i][j] + matrix_b[i][j];
        }
    }
    return ans;
}

void init_matrix(Type *** matrix, int * row, int * col) {//Считывание матрицы
    printf("Input two integers:\n");
    scanf("%d%d", row, col);
    printf("Input matrix with %d rows ans %d columns:\n", *row, *col);
    allocate_memory(matrix, *row, *col);
    for (int i = 0; i < (*row); ++i) {
        for (int j = 0; j < (*col); ++j) {
            scanf("%lf", &(*matrix)[i][j]);
        }
    }
}

void output_matrix(Type ** matrix, int row, int col) {//Вывод матрицы
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            if (j == 0) {
                if (row == 1)
                    printf("( ");
                else if (i == 0)
                    printf("/ ");
                else if (i == row - 1)
                    printf("\\ ");
                else
                    printf("| ");
            }
            printf("%10.4lf", matrix[i][j]);
            if (i == row - 1 && j == col - 1)
                printf("  ");
            else
                printf(", ");
            if (j == col - 1) {
                if (row == 1)
                    printf(")");
                else if (i == 0)
                    printf("\\");
                else if (i == row - 1)
                    printf("/");
                else
                    printf("|");
            }
        }
        puts("");
    }
}
