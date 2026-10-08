/* fmttest.c -- compare fE/fF/ld_* helpers against gfortran reference output */
#include "whamp.h"

int main(void)
{
    static const double vals[] = {
        8.9786, 2.8, 1.0e6, 0.01, 0.001, 1.0, 5.0, 0.0, 1.0640799e-1,
        0.999, 9.999, 3.34e-11, -0.10640799, -1.2345e-100, 1.5, 3.141592653589793,
        123456.789, 1e-4, 1e-16, 1.2345678901234567, -1.0, 1.2345678901234567e17,
        1e16, 1e17, 1e18, 1e20, 1e-5, 1e-6, 1e-10, 0.5, 0.05, 12.34, 123.4, -12.34,
        0.0022, 1.0640799e-01, 3.34e-11, 1.0e-3, 2.0, 3.0
    };
    size_t k;
    /* E descriptors used in the Fortran sources */
    for (k = 0; k < sizeof vals / sizeof vals[0]; k++) {
        putchar('<'); fE(12, 2, 0, vals[k]); putchar('>'); putchar('\n');
        putchar('<'); fE(13, 5, 0, vals[k]); putchar('>'); putchar('\n');
        putchar('<'); fE(12, 3, 0, vals[k]); putchar('>'); putchar('\n');
        putchar('<'); fE(11, 5, 1, vals[k]); putchar('>'); putchar('\n');
        putchar('<'); fE(12, 5, 1, vals[k]); putchar('>'); putchar('\n');
        putchar('<'); fE(8, 2, 0, vals[k]); putchar('>'); putchar('\n');
        putchar('<'); fE(9, 3, 0, vals[k]); putchar('>'); putchar('\n');
        putchar('<'); fE(12, 3, 1, vals[k]); putchar('>'); putchar('\n');
        /* F descriptors */
        putchar('<'); fF(11, 4, vals[k]); putchar('>'); putchar('\n');
        putchar('<'); fF(10, 4, vals[k]); putchar('>'); putchar('\n');
        putchar('<'); fF(9, 5, vals[k]); putchar('>'); putchar('\n');
        putchar('<'); fF(4, 2, vals[k]); putchar('>'); putchar('\n');
        putchar('<'); fF(5, 2, vals[k]); putchar('>'); putchar('\n');
        putchar('<'); fF(6, 3, vals[k]); putchar('>'); putchar('\n');
        putchar('<'); fF(10, 7, vals[k]); putchar('>'); putchar('\n');
        /* list directed */
        putchar('<'); ld_real(vals[k]); putchar('>'); putchar('\n');
    }
    return 0;
}
