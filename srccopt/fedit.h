/*
 * fcompat.h -- emulation of the gfortran output editing that WHAMP relies on.
 *
 * WHAMP's terminal protocol is part of its public interface: the numbers
 * printed by the interactive host are consumed by scripts and by the MATLAB
 * front end, so the C translation must reproduce gfortran's edited output
 * byte for byte.  The helpers here implement the subset of Fortran edit
 * descriptors used by the program:
 *
 *   femit()  - right-justified field with gfortran's overflow fallbacks
 *   fE()     - Ew.d with P scaling of 0 or 1
 *   fF()     - Fw.d
 *   ld_*()   - list-directed items (real, integer, logical, complex)
 *
 * All of them write to stdout.
 */
#ifndef WHAMP_FEDIT_H
#define WHAMP_FEDIT_H

#include <complex.h>
#include <stdio.h>

typedef double complex cd;

/* Scratch buffer size shared by the editors; comfortably larger than any
 * field the formats in use can produce. */
#define FCOMPAT_BUF 512

/* Emit an edited value right-justified in a field of width w, applying
 * gfortran's fallbacks: when the value does not fit, the leading zero of
 * "0."/"-0." is dropped; if it still does not fit the field is filled
 * with '*'. */
void femit(int w, const char *s);

/* Emulate Fortran Ew.d editing with P scaling of 0 or 1.
 * p=0 prints 0.d1d2...d E±ee (Fortran default), p=1 prints d.d1...d E±ee.
 * gfortran leaves out the 'E' when the exponent needs three digits.
 * Rounding is done by the C library on the exact binary value, which is
 * what gfortran does as well. */
void fE(int w, int d, int p, double v);

/* Emulate Fortran Fw.d editing (round-half-even on the binary value). */
void fF(int w, int d, double v);

/* Emulate Fortran list-directed editing of a real(8) item in field w.
 * See the comment in fcompat.c for the exact gfortran rules reproduced. */
void ld_real_w(int w, double v);

/* Fortran list-directed write of one real(8) item (field width 26). */
void ld_real(double v);
/* Fortran list-directed write of one integer item (field width 12). */
void ld_int(int v);
/* Fortran list-directed write of one logical item. */
void ld_logical(int v);
/* Fortran list-directed write of one complex(8) item. */
void ld_complex(cd z);

/* Emulate Fortran READ(*,'(A)'): read one line, pad with blanks to exactly
 * len characters (buf must hold len+1 bytes), and return 0 on success or
 * nonzero on end-of-file. */
int fortran_read_line(char *buf, int len);

#endif /* WHAMP_FEDIT_H */
