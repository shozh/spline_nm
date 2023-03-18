#include <stdio.h>
#include <stdlib.h>

float* f(float b) {
    float* c = (float*)malloc(3 * sizeof(float));
    c[0] = b-2;
    c[1] = b + 8;
    c[2] = ++b;
    return c;
}

int main(void) {

    float* a = f(3);
    printf("%f", a[2]);

    return 0;
}
