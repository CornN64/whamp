#include "whamp.h"
int main(int argc, char **argv)
{
    double v = atof(argv[1]);
    int mode = atoi(argv[2]);
    printf("mode=%d\n", mode);
    fflush(stdout);
    switch (mode) {
    case 0: putchar('<'); fE(12, 2, 0, v); puts(">"); break;
    case 1: putchar('<'); fE(13, 5, 0, v); puts(">"); break;
    case 2: putchar('<'); fE(12, 3, 0, v); puts(">"); break;
    case 3: putchar('<'); fE(11, 5, 1, v); puts(">"); break;
    case 4: putchar('<'); fE(12, 5, 1, v); puts(">"); break;
    case 5: putchar('<'); fE(8, 2, 0, v); puts(">"); break;
    case 6: putchar('<'); fE(9, 3, 0, v); puts(">"); break;
    case 7: putchar('<'); fE(12, 3, 1, v); puts(">"); break;
    case 8: putchar('<'); fF(11, 4, v); puts(">"); break;
    case 9: putchar('<'); fF(10, 4, v); puts(">"); break;
    case 10: putchar('<'); fF(9, 5, v); puts(">"); break;
    case 11: putchar('<'); fF(4, 2, v); puts(">"); break;
    case 12: putchar('<'); fF(5, 2, v); puts(">"); break;
    case 13: putchar('<'); fF(6, 3, v); puts(">"); break;
    case 14: putchar('<'); fF(10, 7, v); puts(">"); break;
    case 15: putchar('<'); ld_real(v); puts(">"); break;
    }
    return 0;
}
