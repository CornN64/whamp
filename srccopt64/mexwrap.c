/*
 * mexwrap.c -- MATLAB MEX wrapper with whamp.m's struct interface.
 *
 *     Output = mexwrap(PlasmaModel, InputParameters)
 *
 * It takes exactly the same PlasmaModel / InputParameters structures as
 * whamp.m and returns the same Output structure, but the dispersion solving
 * is done by the compiled C engine in this directory rather than by the
 * MATLAB implementation.  It is therefore a compiled, self-contained
 * equivalent of whamp.m: the struct parsing, the unit conversions, the sweep
 * axis handling and the output shaping all live here, and libwhamp.a is
 * linked in, so no MATLAB-side glue is needed to run it.
 *
 * Relationship to mexwhamp.c:
 *
 *   mexwhamp.c is the faithful translation of ../src/mexwhamp.F.  It speaks
 *   the engine's flat 13-input/21-output list and is what matlab/+whamp/run.m
 *   calls.  mexwrap.c is the higher-level interface: it speaks the structs and
 *   does in C what whamp.m does in MATLAB.  Both link the same engine.
 *
 * Semantics copied from whamp.m:
 *
 *   - Missing structure fields take whamp.m's defaults (default_plasma and
 *     default_parameters).  A missing field is filled in; a field that is
 *     present but malformed is an error.
 *   - Units: density cm^-3 -> m^-3 (x1e6), temperature eV -> keV (/1e3),
 *     gyrofrequency kHz = 0.0279928 * B(nT).
 *   - kperp/kpar are written [start step end]; the engine keeps the first and
 *     third numbers as the bounds and the second as the step, so the request
 *     becomes [PM1 PM2 PM3] = [start end step].  A bare scalar means one
 *     point and gets the engine's default step of 10.
 *   - Outputs are transposed to the established [kpar x kperp] shape, and
 *     every point the engine reports as non-converged is returned as NaN.
 *   - useLog interprets the axes logarithmically and the reported axes are
 *     converted back with 10.^axis.
 *   - Like whamp.m, every result field is returned complex.  The engine
 *     computes the Poynting flux, EB, the group-velocity and growth-rate
 *     components and u as real quantities, so their imaginary parts are zero.
 *
 * One deliberate difference from whamp.m: a zero axis step is rejected.  The
 * sweep advances the axis by adding the step, so a zero step never leaves the
 * first point and the loop never terminates.  whamp.m accepts a zero step
 * whenever start equals end and then hangs; a MEX that hangs MATLAB is not an
 * acceptable failure mode, so the case is reported instead.  Use a scalar
 * axis, or any non-zero step, for a single point.
 *
 * Build (from this directory, with MATLAB's mex on PATH):
 *
 *     make mexwrap
 *
 * which links libwhamp.a with mexwrap.c into mexwrap.mex<platform>.  The
 * matching matlab/mexwrap.m finds the result on the MATLAB path and calls it,
 * so only the engine directory needs to be added to the path.
 */
#include <math.h>
#include <string.h>

#include "mex.h"

#include "whamp.h"
#include "comin.h"
#include "comoutput.h"
#include "comtrack.h"
#include "whamp_track.h"

#define NSLOT 10 /* the engine keeps ten plasma component slots */

/* Refuse grids large enough to exhaust memory before the engine can report
 * anything; a runaway allocation kills MATLAB rather than the MEX. */
#define MAX_AXIS_POINTS 1e6

/* Electron gyrofrequency in kHz per nanotesla, as whamp.m and run.m use it. */
#define FCE_PER_NANOTESLA 0.0279928

/* Result fields in whamp.m's Output order, with the engine buffer behind
 * each one and whether that buffer holds C99 complex values. */
#define NRESULT 16
#define NFLAGS 3
#define NTRACKEXTRA 5 /* whamp_track_v2.m's five tracking arrays */

/* Output structure size: the four leading fields and the results are always
 * there; the tracker adds the five tracking arrays and `diagnostics`. */
#define NOUTFIELDS_ENGINE (4 + NRESULT + NFLAGS)
#define NOUTFIELDS_TRACK (NOUTFIELDS_ENGINE + NTRACKEXTRA + 1)

typedef struct {
    const char *name;
    const void *buffer;
    int isComplex;
} ResultField;

/* ---- scalar and vector input helpers ----------------------------------- *
 *
 * Any real numeric class is accepted and converted to double, because
 * callers routinely build these inputs with whatever class is at hand.
 */

/* Convert element `i` of the real numeric array `a` to double.  Returns 0 if
 * `a` is not a real numeric array or the element is not finite. */
static int numeric_element(const mxArray *a, size_t i, double *out)
{
    mxClassID id;
    const void *base;
    double v;

    if (!mxIsNumeric(a) || mxIsLogical(a) || mxIsComplex(a))
        return 0;

    if (mxIsDouble(a)) {
        v = ((const double *)mxGetData(a))[i];
    } else {
        id = mxGetClassID(a);
        base = mxGetData(a);
        switch (id) {
        case mxSINGLE_CLASS:
            v = ((const float *)base)[i];
            break;
        case mxINT32_CLASS:
            v = (double)((const int *)base)[i];
            break;
        case mxUINT32_CLASS:
            v = (double)((const unsigned int *)base)[i];
            break;
        case mxINT16_CLASS:
            v = (double)((const short *)base)[i];
            break;
        case mxUINT16_CLASS:
            v = (double)((const unsigned short *)base)[i];
            break;
        case mxINT8_CLASS:
            v = (double)((const signed char *)base)[i];
            break;
        case mxUINT8_CLASS:
            v = (double)((const unsigned char *)base)[i];
            break;
        default:
            return 0;
        }
    }
    if (!mxIsFinite(v))
        return 0;
    *out = v;
    return 1;
}

static int numeric_scalar(const mxArray *a, double *out)
{
    return mxGetNumberOfElements(a) == 1 ? numeric_element(a, 0, out) : 0;
}

/* Read one scalar field of a structure.  A missing or empty field yields
 * `dflt`; a field that is present but not a finite real scalar is an error,
 * which is whamp.m's rule: complete_* fills the gaps, then validates. */
static double field_scalar(const mxArray *s, mwIndex element, const char *field,
                           double dflt, const char *context)
{
    const mxArray *a = mxGetField(s, element, field);
    double v;

    if (a == NULL || mxIsEmpty(a))
        return dflt;
    if (!numeric_scalar(a, &v))
        mexErrMsgIdAndTxt("mexwrap:invalidInput",
                          "%s.%s must be a finite real scalar.", context,
                          field);
    return v;
}

/* ---- plasma model ------------------------------------------------------ */

/* whamp.m's default_plasma: 100 nT with an O+ beam on electrons, as
 * m, n, t, a, b, d, vd. */
static const double defaultSpecies[2][7] = {
    {16.0, 1.0, 10.0, 5.0, 0.0, 1.0, 1.0},
    {0.0, 1.0, 100.0, 1.0, 0.0, 1.0, 0.0}};

/* whamp.m's speciesDefaults, in the same field order. */
static const double speciesDefault[7] = {0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0};
static const char *const speciesField[7] = {"m", "n", "t", "a", "b", "d", "vd"};

/* Fill the engine's component arrays from PlasmaModel.Species.
 *
 * The species list may be a struct array or a cell array of structs, as
 * whamp.m's normalise_species allows; a missing or empty Species selects the
 * built-in two-species model.  Unused component slots are zero except delta,
 * which is one, matching whamp.m's species_arrays. */
static void load_species(const mxArray *plasma)
{
    const mxArray *list = mxGetField(plasma, 0, "Species");
    int useDefaults, nSpecies, j, f;

    useDefaults = (list == NULL || mxIsEmpty(list));
    if (useDefaults) {
        nSpecies = 2;
    } else if (mxIsCell(list)) {
        nSpecies = (int)mxGetNumberOfElements(list);
    } else if (mxIsStruct(list)) {
        nSpecies = (int)mxGetNumberOfElements(list);
    } else {
        mexErrMsgIdAndTxt("mexwrap:invalidSpecies",
                          "PlasmaModel.Species must be a structure or cell "
                          "array.");
        return; /* not reached */
    }

    if (nSpecies > NSLOT)
        mexErrMsgIdAndTxt("mexwrap:tooManySpecies",
                          "WHAMP supports at most %d species (%d given).",
                          NSLOT, nSpecies);

    for (j = 1; j <= NSLOT; j++) {
        ASS[j] = 0.0;
        DN[j] = 0.0;
        TA[j] = 0.0;
        DD[j] = 1.0;
        VD[j] = 0.0;
        AA[j][1] = 0.0;
        AA[j][2] = 0.0;
    }

    for (j = 0; j < nSpecies; j++) {
        /* For a struct array every species lives inside `list` at element j;
         * for a cell array the cell contents are themselves the structure, so
         * the element index is 0. */
        const mxArray *element = NULL;
        mwIndex elementIndex = 0;
        double value[7];

        if (useDefaults) {
            for (f = 0; f < 7; f++)
                value[f] = defaultSpecies[j][f];
        } else {
            if (mxIsCell(list)) {
                element = mxGetCell(list, (mwSize)j);
                if (element == NULL || !mxIsStruct(element) ||
                    mxGetNumberOfElements(element) != 1)
                    mexErrMsgIdAndTxt("mexwrap:invalidSpecies",
                                      "PlasmaModel.Species{%d} must be a "
                                      "scalar structure.",
                                      j + 1);
            } else {
                element = list;
                elementIndex = (mwIndex)j;
            }
            for (f = 0; f < 7; f++)
                value[f] = field_scalar(element, elementIndex, speciesField[f],
                                        speciesDefault[f], "Species");
        }

        ASS[j + 1] = value[0];
        DN[j + 1] = value[1] * 1e6; /* cm^-3 -> m^-3 */
        TA[j + 1] = value[2] / 1e3; /* eV -> keV */
        AA[j + 1][1] = value[3];
        AA[j + 1][2] = value[4];
        DD[j + 1] = value[5];
        VD[j + 1] = value[6];
    }
}

/* ---- sweep axes -------------------------------------------------------- */

/* Read InputParameters.<field> and reduce it to the engine's [PM1 PM2 PM3].
 * The caller writes [start step end]; a bare scalar means one point and keeps
 * the engine's default step of 10.  A missing field takes `dflt`, which is
 * already in engine order. */
static void axis_spec(const mxArray *params, const char *field,
                      const double *dflt, double *spec)
{
    const mxArray *a = mxGetField(params, 0, field);
    double v[3];
    int n, i;

    if (a == NULL || mxIsEmpty(a)) {
        spec[1] = dflt[0];
        spec[2] = dflt[1];
        spec[3] = dflt[2];
        return;
    }

    n = (int)mxGetNumberOfElements(a);
    if (n != 1 && n != 3)
        mexErrMsgIdAndTxt("mexwrap:invalidWaveNumber",
                          "InputParameters.%s must be a finite scalar or "
                          "[start step end].",
                          field);

    for (i = 0; i < n; i++) {
        if (!numeric_element(a, (size_t)i, &v[i]))
            mexErrMsgIdAndTxt("mexwrap:invalidWaveNumber",
                              "InputParameters.%s must be a finite scalar or "
                              "[start step end].",
                              field);
    }

    if (n == 1) {
        spec[1] = v[0];
        spec[2] = v[0];
        spec[3] = 10.0;
    } else {
        spec[1] = v[0];
        spec[2] = v[2];
        spec[3] = v[1];
    }

    /* A zero step makes the sweep advance by nothing, so it never reaches the
     * end of the axis and never terminates.  whamp.m accepts it when start
     * equals end; here it is reported instead. */
    if (spec[3] == 0.0)
        mexErrMsgIdAndTxt("mexwrap:invalidWaveNumberStep",
                          "The %s step must be non-zero. Use a scalar %s, or "
                          "any non-zero step, for a single point.",
                          field, field);
    if (fabs(spec[2] - spec[1]) / fabs(spec[3]) > MAX_AXIS_POINTS)
        mexErrMsgIdAndTxt("mexwrap:invalidWaveNumberRange",
                          "InputParameters.%s spans more than %g grid points.",
                          field, MAX_AXIS_POINTS);
}

/* ---- output helpers ---------------------------------------------------- */

/* The engine stores results as [kperp x kpar] in Fortran order; whamp.m
 * returns [kpar x kperp].  These helpers transpose on the way out and blank
 * the non-converged points, which whamp.m does with NaN. */
static mxArray *push_transposed_real(const double *buf, int nrow, int ncol,
                                     const int *missing)
{
    mxArray *a = mxCreateDoubleMatrix((mwSize)ncol, (mwSize)nrow, mxCOMPLEX);
    int i, j;

#if MX_HAS_INTERLEAVED_COMPLEX
    double *d = (double *)mxGetData(a);
    for (j = 1; j <= nrow; j++)
        for (i = 1; i <= ncol; i++) {
            size_t o = (size_t)(i - 1) + (size_t)(j - 1) * (size_t)ncol;
            int blank = F2D(missing, j, i);
            d[2 * o] = blank ? mxGetNaN() : F2D(buf, j, i);
            d[2 * o + 1] = blank ? mxGetNaN() : 0.0;
        }
#else
    /* Split storage: mxGetData is the real plane alone and mxGetImagData the
     * imaginary one, each nrow*ncol wide.  Writing the imaginary part through
     * the real pointer would run past the end of the buffer. */
    double *re = (double *)mxGetData(a);
    double *im = (double *)mxGetImagData(a);
    for (j = 1; j <= nrow; j++)
        for (i = 1; i <= ncol; i++) {
            size_t o = (size_t)(i - 1) + (size_t)(j - 1) * (size_t)ncol;
            int blank = F2D(missing, j, i);
            re[o] = blank ? mxGetNaN() : F2D(buf, j, i);
            im[o] = blank ? mxGetNaN() : 0.0;
        }
#endif
    return a;
}

/* A real double grid, transposed the same way.  Used for certWinding, which
 * whamp_track_v2.m keeps real and leaves NaN where certification did not run,
 * so there is nothing to blank. */
static mxArray *push_transposed_double(const double *buf, int nrow, int ncol)
{
    mxArray *a = mxCreateDoubleMatrix((mwSize)ncol, (mwSize)nrow, mxREAL);
    double *d = (double *)mxGetData(a);
    int i, j;

    for (j = 1; j <= nrow; j++)
        for (i = 1; i <= ncol; i++)
            d[(size_t)(i - 1) + (size_t)(j - 1) * (size_t)ncol] =
                F2D(buf, j, i);
    return a;
}

/* `missing` may be NULL, which leaves every point as computed; whamp.m blanks
 * only the result fields, not the tracker's `predicted`. */
static mxArray *push_transposed_complex(const cd *buf, int nrow, int ncol,
                                        const int *missing)
{
    mxArray *a = mxCreateDoubleMatrix((mwSize)ncol, (mwSize)nrow, mxCOMPLEX);
    int i, j;

#if MX_HAS_INTERLEAVED_COMPLEX
    double *d = (double *)mxGetData(a);
    for (j = 1; j <= nrow; j++)
        for (i = 1; i <= ncol; i++) {
            size_t o = (size_t)(i - 1) + (size_t)(j - 1) * (size_t)ncol;
            cd v = F2D(buf, j, i);
            int blank = missing != NULL && F2D(missing, j, i);
            d[2 * o] = blank ? mxGetNaN() : creal(v);
            d[2 * o + 1] = blank ? mxGetNaN() : cimag(v);
        }
#else
    double *re = (double *)mxGetData(a);
    double *im = (double *)mxGetImagData(a);
    for (j = 1; j <= nrow; j++)
        for (i = 1; i <= ncol; i++) {
            size_t o = (size_t)(i - 1) + (size_t)(j - 1) * (size_t)ncol;
            cd v = F2D(buf, j, i);
            int blank = missing != NULL && F2D(missing, j, i);
            re[o] = blank ? mxGetNaN() : creal(v);
            im[o] = blank ? mxGetNaN() : cimag(v);
        }
#endif
    return a;
}

static mxArray *push_transposed_int(const int *buf, int nrow, int ncol)
{
    mxArray *a = mxCreateNumericMatrix((mwSize)ncol, (mwSize)nrow,
                                       mxINT32_CLASS, mxREAL);
    int *d = (int *)mxGetData(a);
    int i, j;

    for (j = 1; j <= nrow; j++)
        for (i = 1; i <= ncol; i++)
            d[(size_t)(i - 1) + (size_t)(j - 1) * (size_t)ncol] =
                F2D(buf, j, i);
    return a;
}

/* The engine axis buffers are 1-based; whamp.m returns them as row vectors,
 * converted back from the logarithmic sweep coordinate when useLog is set. */
static mxArray *push_axis(const double *buf, int n, int useLog)
{
    mxArray *a = mxCreateDoubleMatrix(1, (mwSize)n, mxREAL);
    double *d = (double *)mxGetData(a);
    int i;

    for (i = 1; i <= n; i++)
        d[i - 1] = useLog ? pow(10.0, buf[i]) : buf[i];
    return a;
}

/* Read Options.<field> as a logical/numeric scalar for the tracker.  A
 * missing or empty field keeps `dflt`, which is how whamp_track_v2.m's
 * local_get treats them. */
static double option_scalar(const mxArray *options, const char *field,
                            double dflt)
{
    const mxArray *a;

    if (options == NULL)
        return dflt;
    a = mxGetField(options, 0, field);
    if (a == NULL || mxIsEmpty(a))
        return dflt;
    if (mxIsLogical(a))
        return (double)(*((const mxLogical *)mxGetData(a)) != 0);
    if (!numeric_scalar(a, &dflt))
        mexErrMsgIdAndTxt("mexwrap:invalidOption",
                          "Options.%s must be a finite real scalar.", field);
    return dflt;
}

/* Fill trackOptions from the optional third argument, mirroring
 * whamp_track_v2.m's local_options: defaults first, then any field the caller
 * supplied.  scanOffsets accepts up to six values; certRadius and certSamples
 * fall back to the defaults when empty, exactly as local_get does. */
static void load_track_options(const mxArray *options, double iterations)
{
    const mxArray *a;
    int i, n;

    track_options_default();

    if (options == NULL)
        return;
    if (!mxIsStruct(options) || mxGetNumberOfElements(options) != 1)
        mexErrMsgIdAndTxt("mexwrap:invalidOptions",
                          "Options must be a scalar structure.");

    trackOptions.trustFactor = option_scalar(options, "trustFactor", 3.0);
    trackOptions.jumpFactor = option_scalar(options, "jumpFactor", 6.0);
    trackOptions.useScan = option_scalar(options, "useScan", 1.0) != 0.0;
    trackOptions.certify = option_scalar(options, "certify", 1.0) != 0.0;
    trackOptions.certRadius = option_scalar(options, "certRadius", 0.0);
    trackOptions.certSamples =
        (int)option_scalar(options, "certSamples", 8.0);
    trackOptions.numScanOffsets = 0;

    a = mxGetField(options, 0, "scanOffsets");
    if (a != NULL && !mxIsEmpty(a)) {
        n = (int)mxGetNumberOfElements(a);
        if (n > 6)
            mexErrMsgIdAndTxt("mexwrap:invalidOption",
                              "Options.scanOffsets has at most six entries "
                              "(%d given).",
                              n);
        for (i = 0; i < n; i++) {
            if (!numeric_element(a, (size_t)i, &trackOptions.scanOffsets[i]))
                mexErrMsgIdAndTxt("mexwrap:invalidOption",
                                  "Options.scanOffsets must hold finite real "
                                  "values.");
        }
        trackOptions.numScanOffsets = n;
    } else {
        /* local_options' default is [0.5 -0.5 1 -1 2 -2] * 0.5, which is the
         * same list track_options_default() installs. */
        trackOptions.numScanOffsets = 6;
    }

    /* whamp_track_v2.m takes maxIterations from Options when given and from
     * InputParameters otherwise; the caller already resolved the latter. */
    if (mxGetField(options, 0, "maxIterations") != NULL &&
        !mxIsEmpty(mxGetField(options, 0, "maxIterations"))) {
        double v = option_scalar(options, "maxIterations", iterations);
        if (v < 1.0 || v != fint(v))
            mexErrMsgIdAndTxt("mexwrap:invalidMaxIterations",
                              "Options.maxIterations must be a positive "
                              "integer scalar.");
        maxIterations = (int)fint(v);
    }
}

/* Output.diagnostics, shaped as whamp_track_v2.m builds it.  `attempts` holds
 * the tracker's per-point record and `nAttempts` its length; the caller owns
 * both.  The `perPoint` cell of Version 1's dump has no C equivalent and is
 * returned as an empty cell. */
static mxArray *push_diagnostics(const TrackDiagnostics *d, int nAttempts)
{
    static const char *const counters[10] = {
        "clampedSteps", "rejected", "rescued", "flagged", "certified",
        "notCertified", "certUnknown", "converged", "pointsVisited",
        "meanAttempts"
    };
    static const char *names[] = {
        "attemptsPerPoint", "clampedSteps", "rejected", "rescued", "flagged",
        "certified",        "notCertified", "certUnknown", "converged",
        "pointsVisited",    "meanAttempts", "perPoint"
    };
    mxArray *s = mxCreateStructMatrix(1, 1, 12, names);
    mxArray *a;
    double *v;
    int k;

    a = mxCreateDoubleMatrix(1, (mwSize)(nAttempts > 0 ? nAttempts : 0),
                             mxREAL);
    if (nAttempts > 0) {
        v = (double *)mxGetData(a);
        for (k = 0; k < nAttempts; k++)
            v[k] = d->attemptsPerPoint[k];
    }
    mxSetField(s, 0, "attemptsPerPoint", a);

    /* The counters in the same order whamp_track_v2.m stores them. */
    {
        double values[10];
        values[0] = (double)d->clampedSteps;
        values[1] = (double)d->rejected;
        values[2] = (double)d->rescued;
        values[3] = (double)d->flagged;
        values[4] = (double)d->certified;
        values[5] = (double)d->notCertified;
        values[6] = (double)d->certUnknown;
        values[7] = (double)d->converged;
        values[8] = (double)d->pointsVisited;
        values[9] = d->meanAttempts;
        for (k = 0; k < 10; k++) {
            a = mxCreateDoubleScalar(values[k]);
            mxSetField(s, 0, counters[k], a);
        }
    }

    mxSetField(s, 0, "perPoint", mxCreateCellMatrix(0, 0));
    return s;
}

/* Release the engine's output buffers so they are not pinned between calls.
 * Setting the globals to NULL is required: ALLOC1/ALLOC2 free whatever the
 * pointer holds at the start of the next call. */
static void release_output_buffers(void)
{
    free(kperpOUT);
    free(kparOUT);
    free(fOUT);
    free(ExOUT);
    free(EyOUT);
    free(EzOUT);
    free(BxOUT);
    free(ByOUT);
    free(BzOUT);
    free(SxOUT);
    free(SyOUT);
    free(SzOUT);
    free(EBOUT);
    free(VGPOUT);
    free(VGZOUT);
    free(SGPOUT);
    free(SGZOUT);
    free(uOUT);
    free(flagSolutionFoundOUT);
    free(flagTooHeavilyDampedOUT);
    free(flagNoConvergenceOUT);

    free(predictedOUT);
    free(certWindingOUT);
    free(flagSuspectBranchOUT);
    free(flagNotCertifiedOUT);
    free(flagCertUnknownOUT);
    predictedOUT = NULL;
    certWindingOUT = NULL;
    flagSuspectBranchOUT = NULL;
    flagNotCertifiedOUT = NULL;
    flagCertUnknownOUT = NULL;

    kperpOUT = kparOUT = NULL;
    fOUT = ExOUT = EyOUT = EzOUT = BxOUT = ByOUT = BzOUT = SxOUT = SyOUT =
        SzOUT = NULL;
    EBOUT = VGPOUT = VGZOUT = SGPOUT = SGZOUT = uOUT = NULL;
    flagSolutionFoundOUT = flagTooHeavilyDampedOUT = flagNoConvergenceOUT =
        NULL;
}

/* ---- gateway ----------------------------------------------------------- */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const mxArray *plasma, *params, *options = NULL;
    const char *names[NOUTFIELDS_TRACK];
    ResultField fields[NRESULT];
    TrackDiagnostics diagnostics;
    mxArray *out;
    double pm[4], zm[4];
    double B, fstart, iterations, useLog, varyKzFirst;
    int useLogFlag;
    int useTracker;
    int kp, kz, k;

    /* whamp.m's default_parameters: one point at kperp = 0, kpar = 0.0022.
     * Stored already in the engine's [start end step] order. */
    static const double defaultKperp[3] = {0.0, 0.0, 10.0};
    static const double defaultKpar[3] = {0.0022, 0.0022, 10.0};

    if (nrhs < 2 || nrhs > 3)
        mexErrMsgIdAndTxt("mexwrap:invalidInput",
                          "mexwrap takes a PlasmaModel, an InputParameters "
                          "structure and optionally Options (%d inputs given).",
                          nrhs);
    if (nlhs > 1)
        mexErrMsgIdAndTxt("mexwrap:invalidOutput",
                          "mexwrap returns a single Output structure.");

    plasma = prhs[0];
    params = prhs[1];
    if (!mxIsStruct(plasma) || mxGetNumberOfElements(plasma) != 1)
        mexErrMsgIdAndTxt("mexwrap:invalidPlasmaModel",
                          "PlasmaModel must be a scalar structure.");
    if (!mxIsStruct(params) || mxGetNumberOfElements(params) != 1)
        mexErrMsgIdAndTxt("mexwrap:invalidInputParameters",
                          "InputParameters must be a scalar structure.");

    /* ---- plasma model ---- */
    B = field_scalar(plasma, 0, "B", 100.0, "PlasmaModel");
    if (B <= 0.0)
        mexErrMsgIdAndTxt("mexwrap:invalidMagneticField",
                          "PlasmaModel.B must be a positive finite scalar.");
    load_species(plasma);
    XC = FCE_PER_NANOTESLA * B;

    /* ---- run parameters ---- */
    fstart = field_scalar(params, 0, "fstart", 0.1, "InputParameters");
    if (fstart <= 0.0)
        mexErrMsgIdAndTxt("mexwrap:invalidStartFrequency",
                          "InputParameters.fstart must be positive.");
    iterations = field_scalar(params, 0, "maxIterations", 50.0,
                              "InputParameters");
    if (iterations < 1.0 || iterations != fint(iterations))
        mexErrMsgIdAndTxt("mexwrap:invalidMaxIterations",
                          "InputParameters.maxIterations must be a positive "
                          "integer scalar.");
    useLog = field_scalar(params, 0, "useLog", 0.0, "InputParameters");
    if (useLog != 0.0 && useLog != 1.0)
        mexErrMsgIdAndTxt("mexwrap:invalidLogOption",
                          "InputParameters.useLog must be 0 or 1.");
    varyKzFirst =
        field_scalar(params, 0, "varyKzFirst", 1.0, "InputParameters");

    axis_spec(params, "kperp", defaultKperp, pm);
    axis_spec(params, "kpar", defaultKpar, zm);

    PM[1] = pm[1];
    PM[2] = pm[2];
    PM[3] = pm[3];
    ZM[1] = zm[1];
    ZM[2] = zm[2];
    ZM[3] = zm[3];
    XOI = fstart;
    PZL = useLog;
    cycleZFirst = (int)varyKzFirst;
    maxIterations = (int)fint(iterations);

    /* OUTPT/INOUT write to stdout and are not reachable from whamp_engine(),
     * but keep the flag off so nothing can print into a MEX call. */
    printDebugInfo = 0;

    /* ---- solver selection ----
     *
     * whamp_track_v2.m is the reference the caller expects, so the branch
     * tracker is the default and the plain sweep is opt-in with
     * Options.solver = 'engine'.  The two differ where a sweep line's
     * first-order extrapolation lands in a neighbouring basin: the tracker
     * predicts by implicit differentiation against a map keyed by grid cell,
     * judges each root against the branch's own measured rate, and certifies
     * it with the argument principle.  On well-behaved lines the tracker's
     * first seed is the engine's own extrapolation tried with an infinite
     * trust radius, so it lands on the same root at the same cost. */
    useTracker = 1;
    if (nrhs == 3) {
        options = prhs[2];
        if (mxIsStruct(options) && mxGetNumberOfElements(options) == 1) {
            const mxArray *a = mxGetField(options, 0, "solver");
            if (a != NULL && !mxIsEmpty(a) && mxIsChar(a)) {
                char which[16];
                if (mxGetString(a, which, sizeof(which)) != 0)
                    mexErrMsgIdAndTxt("mexwrap:invalidOption",
                                      "Options.solver is too long.");
                if (strcmp(which, "engine") == 0)
                    useTracker = 0;
                else if (strcmp(which, "track") != 0)
                    mexErrMsgIdAndTxt("mexwrap:invalidOption",
                                      "Options.solver must be 'track' or "
                                      "'engine'.");
            }
        }
    }

    /* ---- solve ---- */
    if (useTracker) {
        load_track_options(options, iterations);
        memset(&diagnostics, 0, sizeof(diagnostics));
        whamp_track_v2(&diagnostics);
    } else {
        whamp_engine();
    }

    kp = kperpSize;
    kz = kparSize;
    useLogFlag = (useLog == 1.0);

    /* ---- Output structure, in whamp.m's field order ---- */
    fields[0].name = "f";
    fields[0].buffer = fOUT;
    fields[0].isComplex = 1;
    fields[1].name = "Ex";
    fields[1].buffer = ExOUT;
    fields[1].isComplex = 1;
    fields[2].name = "Ey";
    fields[2].buffer = EyOUT;
    fields[2].isComplex = 1;
    fields[3].name = "Ez";
    fields[3].buffer = EzOUT;
    fields[3].isComplex = 1;
    fields[4].name = "Bx";
    fields[4].buffer = BxOUT;
    fields[4].isComplex = 1;
    fields[5].name = "By";
    fields[5].buffer = ByOUT;
    fields[5].isComplex = 1;
    fields[6].name = "Bz";
    fields[6].buffer = BzOUT;
    fields[6].isComplex = 1;
    fields[7].name = "Sx";
    fields[7].buffer = SxOUT;
    fields[7].isComplex = 1;
    fields[8].name = "Sy";
    fields[8].buffer = SyOUT;
    fields[8].isComplex = 1;
    fields[9].name = "Sz";
    fields[9].buffer = SzOUT;
    fields[9].isComplex = 1;
    fields[10].name = "EB";
    fields[10].buffer = EBOUT;
    fields[10].isComplex = 0;
    fields[11].name = "VGP";
    fields[11].buffer = VGPOUT;
    fields[11].isComplex = 0;
    fields[12].name = "VGZ";
    fields[12].buffer = VGZOUT;
    fields[12].isComplex = 0;
    fields[13].name = "SGP";
    fields[13].buffer = SGPOUT;
    fields[13].isComplex = 0;
    fields[14].name = "SGZ";
    fields[14].buffer = SGZOUT;
    fields[14].isComplex = 0;
    fields[15].name = "u";
    fields[15].buffer = uOUT;
    fields[15].isComplex = 0;

    names[0] = "InputParameters";
    names[1] = "PlasmaModel";
    names[2] = "kperp";
    names[3] = "kpar";
    for (k = 0; k < NRESULT; k++)
        names[4 + k] = fields[k].name;
    names[4 + NRESULT] = "flagSolutionFound";
    names[5 + NRESULT] = "flagTooHeavilyDamped";
    names[6 + NRESULT] = "flagNoConvergence";
    if (useTracker) {
        names[7 + NRESULT] = "flagSuspectBranch";
        names[8 + NRESULT] = "certWinding";
        names[9 + NRESULT] = "flagNotCertified";
        names[10 + NRESULT] = "flagCertUnknown";
        names[11 + NRESULT] = "predicted";
        names[12 + NRESULT] = "diagnostics";
    }

    out = mxCreateStructMatrix(1, 1,
                               useTracker ? NOUTFIELDS_TRACK
                                          : NOUTFIELDS_ENGINE,
                               names);
    mxSetField(out, 0, "InputParameters", mxDuplicateArray(params));
    mxSetField(out, 0, "PlasmaModel", mxDuplicateArray(plasma));
    mxSetField(out, 0, "kperp", push_axis(kperpOUT, kp, useLogFlag));
    mxSetField(out, 0, "kpar", push_axis(kparOUT, kz, useLogFlag));

    for (k = 0; k < NRESULT; k++) {
        mxArray *value;
        if (fields[k].isComplex)
            value = push_transposed_complex((const cd *)fields[k].buffer, kp,
                                            kz, flagNoConvergenceOUT);
        else
            value = push_transposed_real((const double *)fields[k].buffer, kp,
                                         kz, flagNoConvergenceOUT);
        mxSetField(out, 0, fields[k].name, value);
    }

    mxSetField(out, 0, "flagSolutionFound",
                     push_transposed_int(flagSolutionFoundOUT, kp, kz));
    mxSetField(out, 0, "flagTooHeavilyDamped",
                     push_transposed_int(flagTooHeavilyDampedOUT, kp, kz));
    mxSetField(out, 0, "flagNoConvergence",
                     push_transposed_int(flagNoConvergenceOUT, kp, kz));

    if (useTracker) {
        /* whamp_track_v2.m's five tracking arrays.  Unlike the result fields
         * they keep whamp_track_v2.m's own classes: the three verdict flags
         * are int32, certWinding is real double, and predicted is complex
         * because the reference is a complex frequency.  None is blanked here
         * -- the tracker leaves them NaN where no verdict was reached. */
        mxSetField(out, 0, "flagSuspectBranch",
                   push_transposed_int(flagSuspectBranchOUT, kp, kz));
        mxSetField(out, 0, "certWinding",
                   push_transposed_double(certWindingOUT, kp, kz));
        mxSetField(out, 0, "flagNotCertified",
                   push_transposed_int(flagNotCertifiedOUT, kp, kz));
        mxSetField(out, 0, "flagCertUnknown",
                   push_transposed_int(flagCertUnknownOUT, kp, kz));
        mxSetField(out, 0, "predicted",
                   push_transposed_complex(predictedOUT, kp, kz, NULL));
        mxSetField(out, 0, "diagnostics",
                   push_diagnostics(&diagnostics, diagnostics.pointsVisited));
        free(diagnostics.attemptsPerPoint);
    }

    release_output_buffers();

    if (nlhs > 0)
        plhs[0] = out;
    else
        mxDestroyArray(out);
}
