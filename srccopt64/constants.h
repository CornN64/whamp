/*
 * constants.h -- every physical, mathematical and algorithmic constant
 * used by the WHAMP C translation, in one place.
 *
 * The values come from the Fortran sources in ../src.  Where those sources
 * were internally inconsistent, this header keeps a single value and the
 * divergence is documented below.  Where a value is a deliberate tuning
 * parameter rather than a constant of nature, its origin is named.
 *
 * Two places where the original Fortran disagreed with itself:
 *
 *   PI      av.f90:36 used 3.1415926535897 (14 digits) while
 *           whamp_engine.f90:22, rasy.f90:8, rint.f90:51 and output.f90:44
 *           used 3.14159265358979 (15 digits).  This header keeps the
 *           longer value.  The short value only ever reached the POYN[]
 *           Poynting estimate in av.c, which the original computes and then
 *           discards, so the choice cannot change any output.
 *
 *   C_LIGHT av.f90 used 299.79 / 2.9979e8 while energy.f90:66 used
 *           299.792458.  This header keeps the exact SI value.  See the
 *           note on C_LIGHT for the numerical consequence.
 *
 * Include this header through whamp.h; do not include it directly.
 */
#ifndef WHAMP_CONSTANTS_H
#define WHAMP_CONSTANTS_H

/* Mathematical ---------------------------------------------------------- */

/*
 * PI, the longer of the two values used in the Fortran sources.
 * 15 significant digits; the nearest double is 3.1415926535897931160.
 */
#define PI 3.14159265358979

/* Physical -------------------------------------------------------------- */

/*
 * Speed of light in the units WHAMP normalises to: km/s, which is also
 * 1/c in the m/s-based normalisation energy.c uses.
 *
 * The exact SI value is 299792.458 m/s, so C_LIGHT is exact rather than a
 * rounded constant, and C_MS is the same speed expressed in m/s.
 *
 * energy.c scales the wave magnetic field with V = 1/C_LIGHT, and av.c then
 * forms the field-energy ratio as |C_LIGHT*B|^2/|E|^2.  Both must therefore
 * use the same C_LIGHT, which the Fortran did not do (av.f90 used 299.79
 * against energy.f90's 299.792458), leaving that ratio wrong by 1.6e-5.
 * Using one value fixes it.
 *
 * NOTE ON A NUMERICAL SENSITIVITY.  av.c computes the field energy as
 *
 *     ENDEN = 2*Re(ERG*conjg(ERG)/(ERG + conjg(ERG))) - ENFIELD
 *
 * where ENFIELD = 1 + |B|^2/|E|^2 * C_LIGHT^2.  The two terms are nearly
 * equal, so ENDEN keeps only the leading digits of their difference.  On
 * the Models/Ex1 case p0z.0022f.1 the terms are 320862.74 and 320847.20,
 * giving ENDEN = 15.54.  Replacing the Fortran 299.79 with the exact
 * 299.792458 moves ENFIELD by 1.6e-5 relative -- but ENDEN from 15.54 to
 * 10.28, a 34% change.  The printed enden= is therefore unreliable in
 * *both* spellings of this constant; the exact value is used because it is
 * correct, not because it is more stable.  Fixing the cancellation would
 * mean reformulating the energy partition, which is out of scope here.
 */
#define C_LIGHT 299.792458
#define C_MS (C_LIGHT * 1.0e6) /* exact: 299792458 m/s */

/*
 * Vacuum permittivity in F/m.  The Fortran used 8.8542e-12, a rounded
 * value (the exact SI value is 8.8541878128e-12); kept as-is because it is
 * part of the unit conversion the original calibrated against.
 */
#define EPSILON_0 8.8542e-12

/*
 * Proton-to-electron mass ratio, truncated to five digits by the Fortran
 * (whamp.f90:82, untyped so single precision there).  Used to turn the
 * per-species mass numbers ASS[] into mass-to-charge ratios.
 */
#define PROTON_MASS_RATIO 1836.1

/*
 * Inverse Debye length factor: DEK is the wavenumber (in 1/m) of a plasma
 * frequency of 1 rad/s, so that DEK * RN normalises lengths.  The Fortran
 * used 12405 (whamp.f90:99).
 */
#define DEBYE_KM 12405.

/*
 * Electron rest energy in keV, used for the relativistic correction
 * CV = TR*(1022 + TR)/(511 + TR)^2 of whamp.f90:104.  1022 is 2*511, so
 * only ELECTRON_REST_ENERGY_KEV is a constant.
 */
#define ELECTRON_REST_ENERGY_KEV 511.

/* Algorithmic ----------------------------------------------------------- */

/*
 * Convergence of the Newton iteration on the complex frequency: the root is
 * accepted once the correction is small relative to the current estimate
 * (FREQ_REL_TOL, whamp.f90:277) or absolutely small (FREQ_ABS_TOL,
 * whamp.f90:278).
 */
#define FREQ_REL_TOL 1.0e-6
#define FREQ_ABS_TOL 1.0e-6

/*
 * Fractional change in |D| that counts as a local minimum of the dispersion
 * function rather than convergence (whamp.f90:292).
 */
#define MINIMA_TOL 1.0e-5

/*
 * Initial step size for the Taylor expansion of the dispersion function:
 * the search radius is RTAY_SCALE * |dD/dY| (rtay.f90:27).
 */
#define RTAY_SCALE 1.0e8

/*
 * A |dD/dY| this small, relative to the previous derivative, means the
 * asymptotic expansion has stopped resolving the root (rasy.f90:37).
 */
#define ASYMP_TOL 1.0e-7

/*
 * Large sentinel used by the asymptotic expansion in place of the Fortran's
 * 1.E99, which rasy.f90:12-15 comments as "too big for S/370 hardware, set
 * to largest possible for IBM machines" before setting C = 7.2d35.  It is a
 * stand-in for infinity in the `C >= T` convergence test, so its exact size
 * does not matter as long as it is larger than any derivative reached.
 */
#define OVERFLOW_CEILING 7.2e35

/*
 * Sentinel that ends the P or Z sweep.  The Fortran used 1.D99, a value no
 * log10 wavenumber approaches; 1e99 is the double equivalent.
 */
#define SWEEP_END 1e99

/*
 * Coefficient converting the E x B product the code forms into a Poynting
 * flux in uW/m^2 for the reference wave <E^2> = 0.5 (mV/m)^2.  Written as
 * the Fortran wrote it (10/(4*PI*2)) so the rounding is unchanged.
 */
#define POYNTING_COEF (10.0 / 4.0 / PI / 2.0)

/* Shared array-indexing macros ----------------------------------------- */

/*
 * Element (I, K) of the 6x4 work array XSI that the Fortran declared
 * dimensioned (6,4) and C flattens column-major.  av.c and difu.c hold it
 * as a local and xsi.c receives it as a parameter, so the macro refers to
 * whichever XSI is in scope.
 */
#define XS(I, K) XSI[(size_t)((I) - 1) + (size_t)((K) - 1) * 6]

/*
 * The four elements of the RC array that RINT, RTAY and RASY fill in and
 * return through their cd *RC parameter: the dispersion function and its
 * derivatives, in the order the Fortran wrote them.
 */
#define RC11 RC[0]
#define RC21 RC[1]
#define RC12 RC[2]
#define RC22 RC[3]

#endif /* WHAMP_CONSTANTS_H */
