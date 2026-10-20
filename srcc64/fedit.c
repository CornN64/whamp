/*
 * fcompat.c -- implementation of the gfortran editing emulation declared in
 * fcompat.h.  See the header for the contract of each routine.
 */
#include "fedit.h"

#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

void femit(int w, const char *s)
{
    int len = (int)strlen(s);
    if (len <= w) {
        printf("%*s", w, s);
        return;
    }

    /* Overflow fallback: gfortran drops the leading zero of "0."/"-0.". */
    {
        const char *p = (*s == '-') ? s + 1 : s;
        if (p[0] == '0' && p[1] == '.') {
            char alt[FCOMPAT_BUF];
            int n = 0;
            if (*s == '-')
                alt[n++] = '-';
            strcpy(alt + n, p + 1);
            if ((int)strlen(alt) <= w) {
                printf("%*s", w, alt);
                return;
            }
        }
    }

    /* Still overflowing: the field is filled with asterisks. */
    for (int i = 0; i < w; i++)
        putchar('*');
}

/* Extract the decimal exponent of a string produced by printf("%e"). */
static int estring_exponent(const char *s, char **mantissaEnd)
{
    char *ep = strchr(s, 'e');
    if (!ep)
        ep = strchr(s, 'E');
    *mantissaEnd = ep;
    return atoi(ep + 1);
}

/* Format the "E±ee" (or bare "±eee" when three digits are needed) suffix
 * that gfortran appends to Ew.d edited values. */
static void append_exponent(char *out, int expo)
{
    sprintf(out, "%s%c%02d", (expo > 99 || expo < -99) ? "" : "E",
            expo < 0 ? '-' : '+', expo < 0 ? -expo : expo);
}

void fE(int w, int d, int p, double v)
{
    char tmp[FCOMPAT_BUF], out[FCOMPAT_BUF];
    int neg = 0;
    char *ep;

    /* gfortran edits IEEE special values as fixed strings right-justified
     * in the field: NaN, Infinity, -Infinity. */
    if (isnan(v)) {
        femit(w, "NaN");
        return;
    }
    if (isinf(v)) {
        femit(w, v < 0.0 ? "-Infinity" : "Infinity");
        return;
    }
    if (v < 0.0 || (v == 0.0 && signbit(v))) {
        neg = 1;
        v = -v;
    }

    if (v == 0.0) {
        int o = 0;
        if (neg)
            out[o++] = '-';
        out[o++] = '0';
        out[o++] = '.';
        for (int i = 0; i < d; i++)
            out[o++] = '0';
        strcpy(out + o, "E+00");
        femit(w, out);
        return;
    }

    if (p == 1) {
        int o = 0, ex;
        snprintf(tmp, sizeof tmp, "%.*e", d, v);
        ex = estring_exponent(tmp, &ep);
        if (neg)
            out[o++] = '-';
        memcpy(out + o, tmp, (size_t)(ep - tmp));
        o += (int)(ep - tmp);
        append_exponent(out + o, ex);
        femit(w, out);
        return;
    }

    /* p == 0: exactly d significant digits, printed as 0.d1..d E±ee */
    snprintf(tmp, sizeof tmp, "%.*e", d > 0 ? d - 1 : 0, v);
    {
        int expo = estring_exponent(tmp, &ep) + 1;
        char digits[FCOMPAT_BUF];
        int nd = 0, o = 0;
        for (const char *q = tmp; q != ep; q++)
            if (isdigit((unsigned char)*q))
                digits[nd++] = *q;
        if (neg)
            out[o++] = '-';
        out[o++] = '0';
        out[o++] = '.';
        for (int i = 0; i < d; i++)
            out[o++] = digits[i < nd ? i : nd - 1];
        append_exponent(out + o, expo);
        femit(w, out);
    }
}

void fF(int w, int d, double v)
{
    char out[FCOMPAT_BUF];

    if (isnan(v)) {
        femit(w, "NaN");
        return;
    }
    if (isinf(v)) {
        femit(w, v < 0.0 ? "-Infinity" : "Infinity");
        return;
    }
    snprintf(out, sizeof out, "%.*f", d, v);
    femit(w, out);
}

void ld_real_w(int w, double v)
{
    char out[64];
    int expo;

    if (v == 0.0 || !isfinite(v)) {
        snprintf(out, sizeof out, "%.*f", 16, v);
        printf("%*s%5s", w - 5, out, "");
        return;
    }

    /* Adjusted exponent E of value = 0.d1d2... x 10^E.  gfortran uses fixed
     * notation for 0 <= E <= 17 and exponential P=1 notation otherwise. */
    expo = (int)floor(log10(fabs(v))) + 1;
    if (expo >= 0 && expo <= 17) {
        char fixed[64];
        int decimals = 17 - expo;
        if (decimals > 0)
            snprintf(fixed, sizeof fixed, "%.*f", decimals, v);
        else
            snprintf(fixed, sizeof fixed, "%.0f.", v); /* Fortran keeps the point */
        printf("%*s%5s", w - 5, fixed, "");
    } else {
        char tmp[64], mant[64];
        char *ep;
        int ex;
        snprintf(tmp, sizeof tmp, "%.16e", v);
        ex = estring_exponent(tmp, &ep);
        *ep = '\0';
        snprintf(mant, sizeof mant, "%s", tmp);
        snprintf(out, sizeof out, "%sE%c%03d", mant, ex < 0 ? '-' : '+', ex < 0 ? -ex : ex);
        printf("%*s", w, out);
    }
}

void ld_real(double v) { ld_real_w(26, v); }
void ld_int(int v) { printf("%12d", v); }
void ld_logical(int v) { printf("%12s", v ? "T" : "F"); }

void ld_complex(cd z)
{
    printf("(");
    ld_real_w(25, creal(z));
    printf(",");
    ld_real_w(25, cimag(z));
    printf(")");
}

int fortran_read_line(char *buf, int len)
{
    int c, n = 0;
    while (n < len) {
        c = getchar();
        if (c == EOF) {
            if (n == 0) {
                buf[0] = '\0';
                return 1;
            }
            break;
        }
        if (c == '\n')
            break;
        buf[n++] = (char)c;
    }
    while (n < len)
        buf[n++] = ' ';
    buf[len] = '\0';
    return 0;
}
