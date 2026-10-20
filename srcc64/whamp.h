/*
 * whamp.h -- shared definitions for the C translation of the WHAMP
 * Fortran sources in ../src.
 *
 * THIS IS THE ALL-DOUBLE VARIANT (srccopt64).
 *
 * The reference translation (../srccopt) reproduces gfortran bit-for-bit,
 * which requires emulating Fortran's single-precision (`real`) entities:
 * the unsuffixed literals and the `real` locals of typin.f90.  That
 * emulation is done with the F32() macro and C `float` variables.
 *
 * In this variant that emulation is removed: every quantity is a C
 * `double`, so no F32() macro and no `float` variables exist.  The
 * numerical results therefore differ slightly from gfortran (typically in
 * the 6th-8th significant digit of printed frequencies and in the tiny
 * imaginary parts).  Use ../srccopt when byte-exact output matters; use
 * this variant when uniform double precision matters more.
 *
 * Conventions used throughout the translation:
 *   - Fortran `real(kind(1.0d0))` (d2p) -> C `double`
 *   - Fortran `complex(kind(1.0d0))`    -> C `double complex` (typedef cd)
 *   - Fortran `real` (default)          -> C `double` (see note above)
 *   - Fortran arrays keep 1-based indexing: an array declared
 *     A(10) in Fortran becomes A[11] in C and is used with indices 1..10.
 *   - Fortran module variables become globals defined in comin.c,
 *     comcout.c and comoutput.c.
 *   - All Fortran output editing goes through the helpers in fedit.h.
 */
#ifndef WHAMP_H
#define WHAMP_H

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "constants.h"
#include "fedit.h"

/* Number of plasma components WHAMP keeps track of (Fortran arrays 1..10). */
#define NSPECIES 10

/* Fortran intrinsic helpers ------------------------------------------- */

/* Fortran INT(x): truncation toward zero. */
static inline double fint(double x) { return trunc(x); }

/* Fortran NINT(x): nearest integer, ties away from zero. */
static inline int fnint(double x)
{
    return (int)(x >= 0.0 ? floor(x + 0.5) : ceil(x - 0.5));
}

/* Fortran SIGN(A,B): magnitude of A, sign of B. */
#define FSIGN(A, B) copysign((A), (B))

/* Fortran CONJG(). */
#define conjg(z) conj(z)

/* Prototypes ------------------------------------------------------------ */

/* energy.c */
void ENERGY(cd U1, cd U3, cd U2, cd U12, cd U32);
/* rint.c */
void RINT(cd YY, double AL, cd *RC);
/* rtay.c */
void RTAY(cd Y, double AL, cd *RC);
/* rasy.c */
void RASY(cd Y, double AL, cd *RC);
/* ryla.c */
void RYLA(cd Y, double AL, cd *RC);
/* xsi.c */
void CHI(cd *XSI, int J, int IB, int KOL, int *IERR);
/* difu.c */
void DIFU(int KOL, int JMAX, int *IERR);
/* av.c */
void AV(void);
/* input.c */
void read_input_file(const char *FILENAME);
/* output.c */
void OUTPT(void);
void INOUT(void);
/* typin.c */
void TYPIN(int *NPL, int *KFS);
/* whamp_engine.c */
void whamp_engine(void);

#endif /* WHAMP_H */
