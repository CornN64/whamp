/*
 * constants.h -- the constants the branch tracker needs, in one place.
 *
 * This is the MINIMAL subset for the tracker (whamp_track.c) and the MATLAB
 * wrappers.  The engine keeps spelling its own constants inline, because this
 * variant reproduces the gfortran reference bit-for-bit and centralising them
 * here would invite a value change that breaks that property.
 *
 * Deliberately absent: PI.  av.c, output.c, rasy.c, rint.c and whamp_engine.c
 * each declare a file-scope `static const double PI`, and a macro of that name
 * would make those declarations ill-formed.  av.c's is the short 14-digit
 * value the Fortran used there; the rest use the 15-digit one.  POYNTING_COEF
 * below therefore spells the 15-digit value out rather than referring to a
 * macro, which keeps it equal to what whamp_engine.c computes.
 *
 * THE VALUES ARE SINGLE-ROUNDED WHERE THE ENGINE SINGLE-ROUNDS.  This variant
 * reproduces gfortran bit-for-bit, and the Fortran wrote 1836.1, 1.E-6 and 1e-5
 * as unsuffixed (single-precision) literals.  whamp_engine.c emulates that with
 * F32(); the tracker has to reach the same mass ratio and the same convergence
 * tolerances as the engine it shares its root map with, or the two disagree
 * about where the root is.  So these three macros expand through (float).
 * DEBYE_KM and ELECTRON_REST_ENERGY_KEV are exactly representable in single and
 * are written plain, matching how the engine writes them.
 *
 * The all-double tree ../srcc64 keeps the same macros unrounded, which is
 * what ../srcm/whamp_track.m uses.
 *
 * Include this header through whamp.h; do not include it directly.
 */
#ifndef WHAMP_CONSTANTS_H
#define WHAMP_CONSTANTS_H

/*
 * Proton-to-electron mass ratio, truncated to five digits by the Fortran
 * (whamp.f90:82, untyped so single precision there).  Used to turn the
 * per-species mass numbers ASS[] into mass-to-charge ratios.
 */
#define PROTON_MASS_RATIO ((double)(float)1836.1)

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

/*
 * Convergence of the Newton iteration on the complex frequency: the root is
 * accepted once the correction is small relative to the current estimate
 * (FREQ_REL_TOL, whamp.f90:277) or absolutely small (FREQ_ABS_TOL,
 * whamp.f90:278).
 */
#define FREQ_REL_TOL ((double)(float)1.0e-6)
#define FREQ_ABS_TOL ((double)(float)1.0e-6)

/*
 * Fractional change in |D| that counts as a local minimum of the dispersion
 * function rather than convergence (whamp.f90:292).
 */
#define MINIMA_TOL ((double)(float)1.0e-5)

/*
 * Coefficient converting the E x B product the code forms into a Poynting
 * flux in uW/m^2 for the reference wave <E^2> = 0.5 (mV/m)^2.  Written as
 * the Fortran wrote it (10/(4*PI*2)) so the rounding is unchanged; the
 * literal is the 15-digit PI that whamp_engine.c uses, so the tracker and
 * the engine report the same Poynting flux.
 */
#define POYNTING_COEF (10.0 / 4.0 / 3.14159265358979 / 2.0)

/*
 * Sentinel that ends the P or Z sweep.  The Fortran used 1.D99, a value no
 * log10 wavenumber approaches; 1e99 is the double equivalent.
 */
#define SWEEP_END 1e99

#endif /* WHAMP_CONSTANTS_H */
