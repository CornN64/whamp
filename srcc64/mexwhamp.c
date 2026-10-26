/*
 * mexwhamp.c -- MATLAB MEX gateway for the C translation of WHAMP.
 *
 * Translation of ../src/mexwhamp.F.  Same 13 inputs and 21 outputs in the
 * same order, so matlab/+whamp/run.m drives it unchanged.
 *
 * Differences from the Fortran gateway (all deliberate):
 *
 *   - The legacy mxCopy*Ptr routines are Fortran-only entry points; they are
 *     not exported from the C import library, so a C gateway cannot use them.
 *     The arrays are read and written directly through mxGetData(), with the
 *     complex layout selected by MX_HAS_INTERLEAVED_COMPLEX so this file
 *     builds both for the -R2017b (interleaved) and pre-R2018a (separate
 *     real/imaginary) APIs.
 *
 *   - Every input is validated.  The Fortran gateway copied m*n elements of
 *     whatever array was handed in into the module scalar it was filling, so
 *     a wrong size silently corrupted memory -- and, inside a MEX, takes
 *     MATLAB down with it.  Here a wrong size is a clean mexErrMsgIdAndTxt.
 *
 *   - AA is a Fortran (10,2) array, i.e. column-major, while the C
 *     translation stores it as AA[11][3].  The 10-by-2 MATLAB matrix is
 *     therefore unpacked element by element rather than copied wholesale.
 *
 *   - The engine allocates its output buffers with calloc() on every call
 *     (ALLOC1/ALLOC2 free whatever the pointer holds first).  Once the
 *     results are copied into MATLAB arrays the buffers are released and the
 *     globals set to NULL, so a long MATLAB session does not pin the grid
 *     from the previous call.
 *
 *   - The Fortran gateway required exactly 21 outputs; here fewer are
 *     allowed, which is the usual MEX convention and lets a caller grab just
 *     the frequencies.
 *
 * Build (from this directory, with MATLAB's mex on PATH):
 *
 *     make mex
 *
 * which links libwhamp.a with mexwhamp.c into mexwhamp.mex<platform>.  Copy
 * or move the result next to matlab/+whamp/run.m, i.e. into matlab/+whamp/.
 *
 * The engine interface (comin.h, comoutput.h) is identical across srcc/,
 * srcc64/ and srccopt64/, so this file drops into any of them unchanged.
 * srccopt/ is the natural home because it reproduces the Fortran results
 * bit-for-bit; see ../REPORT.md, section 2.  The newer mexwrap.c is the
 * recommended gateway and is built in all four C trees.
 */
#include <string.h>

#include "mex.h"

#include "whamp.h"
#include "comin.h"
#include "comcout.h"
#include "comoutput.h"

#define NINPUT  13
#define NOUTPUT 21
#define NSPEC   10 /* Fortran species arrays are (10) and (10,2) */

/* Refuse grids large enough to exhaust memory before the engine can report
 * anything; a runaway allocation kills MATLAB rather than the MEX. */
#define MAX_AXIS_POINTS 1e6

/* ---- input helpers ---------------------------------------------------- *
 *
 * Any real numeric class is accepted and converted to double, because
 * run.m passes int32 for maxIterations and callers routinely build the
 * species vectors with whatever class is at hand.  64-bit integers beyond
 * 2^53 lose precision, which is not meaningful for these inputs.
 */

/* Validate that `a` is a real numeric array with exactly `nel` elements and
 * return a pointer to its data, converted to double if necessary.  The
 * returned pointer is either the array's own storage (read-only) or a
 * pointer into `scratch`, which must hold at least `nel` doubles. */
static const double *real_in(const mxArray *a, const char *name, size_t nel,
                             double *scratch)
{
    mxClassID id;
    const void *base;
    size_t i;

    if (!mxIsNumeric(a) || mxIsLogical(a) || mxIsComplex(a))
        mexErrMsgIdAndTxt("WHAMP:inputType",
                          "Input '%s' must be a real numeric array, not '%s'.",
                          name, mxGetClassName(a));
    if ((size_t)mxGetNumberOfElements(a) != nel)
        mexErrMsgIdAndTxt("WHAMP:inputSize",
                          "Input '%s' must have exactly %lu elements (%lu given).",
                          name, (unsigned long)nel,
                          (unsigned long)mxGetNumberOfElements(a));

    if (mxIsDouble(a))
        return (const double *)mxGetData(a);

    id = mxGetClassID(a);
    base = mxGetData(a);
    for (i = 0; i < nel; i++) {
        switch (id) {
        case mxSINGLE_CLASS:
            scratch[i] = ((const float *)base)[i];
            break;
        case mxINT32_CLASS:
            scratch[i] = (double)((const int *)base)[i];
            break;
        case mxUINT32_CLASS:
            scratch[i] = (double)((const unsigned int *)base)[i];
            break;
        case mxINT16_CLASS:
            scratch[i] = (double)((const short *)base)[i];
            break;
        case mxUINT16_CLASS:
            scratch[i] = (double)((const unsigned short *)base)[i];
            break;
        case mxINT8_CLASS:
            scratch[i] = (double)((const signed char *)base)[i];
            break;
        case mxUINT8_CLASS:
            scratch[i] = (double)((const unsigned char *)base)[i];
            break;
        default:
            mexErrMsgIdAndTxt("WHAMP:inputType",
                              "Input '%s' has unsupported type '%s'.",
                              name, mxGetClassName(a));
            return NULL; /* not reached */
        }
    }
    return scratch;
}

static double scalar_in(const mxArray *a, const char *name)
{
    double scratch[1];
    if (mxGetNumberOfElements(a) != 1)
        mexErrMsgIdAndTxt("WHAMP:inputSize",
                          "Input '%s' must be a scalar (%lu elements given).",
                          name, (unsigned long)mxGetNumberOfElements(a));
    return real_in(a, name, 1, scratch)[0];
}

static int int_in(const mxArray *a, const char *name)
{
    return (int)fint(scalar_in(a, name));
}

/* Copy a length-n MATLAB vector into a 1-based Fortran-style array. */
static void vector_in(const mxArray *a, const char *name, double *dst, int n)
{
    double scratch[NSPEC];
    const double *src = real_in(a, name, (size_t)n, scratch);
    int i;
    for (i = 1; i <= n; i++)
        dst[i] = src[i - 1];
}

/* AA is Fortran AA(10,2) (column-major) but C AA[11][3]; unpack it. */
static void species_ab_in(const mxArray *a, const char *name)
{
    double scratch[2 * NSPEC];
    const double *src;
    int j;

    if (!mxIsNumeric(a) || mxIsLogical(a) || mxIsComplex(a))
        mexErrMsgIdAndTxt("WHAMP:inputType",
                          "Input '%s' must be a real numeric array, not '%s'.",
                          name, mxGetClassName(a));
    if (mxGetM(a) != (mwSize)NSPEC || mxGetN(a) != 2)
        mexErrMsgIdAndTxt("WHAMP:inputSize",
                          "Input '%s' must be %d-by-2 (%lu-by-%lu given).",
                          name, NSPEC, (unsigned long)mxGetM(a),
                          (unsigned long)mxGetN(a));

    src = real_in(a, name, 2 * NSPEC, scratch);
    for (j = 1; j <= NSPEC; j++) {
        AA[j][1] = src[(size_t)(j - 1)];
        AA[j][2] = src[(size_t)NSPEC + (size_t)(j - 1)];
    }
}

/* Validate an (min,max,step) axis.
 *
 * A zero step is rejected even when min == max: the engine divides the
 * offset of each output point by the step to find its slot
 * (save_output's indexKpar = 1 + NINT((ZLG - ZM(1)) / ZM(3))), so a zero
 * step divides by zero and writes the result out of bounds.  That is a
 * pre-existing bug in the Fortran too -- src/whamp_engine.f90 segfaults the
 * same way -- but crashing MATLAB is not an acceptable failure mode for a
 * MEX, so it is caught here instead.
 */
static void check_axis(const double *spec, const char *name)
{
    if (spec[3] == 0.0)
        mexErrMsgIdAndTxt("WHAMP:inputSize",
                          "Input '%s' has a zero step; use a step larger than "
                          "zero even for a single grid point.", name);
    if (fabs(spec[2] - spec[1]) / fabs(spec[3]) > MAX_AXIS_POINTS)
        mexErrMsgIdAndTxt("WHAMP:inputSize",
                          "Input '%s' spans more than %g grid points.",
                          name, MAX_AXIS_POINTS);
}

/* ---- output helpers --------------------------------------------------- */

/* Engine axis buffers are 1-based (element 0 unused); return a 1-by-n row. */
static mxArray *push_axis(const double *buf, int n)
{
    mxArray *a = mxCreateDoubleMatrix(1, (mwSize)n, mxREAL);
    if (n > 0)
        memcpy(mxGetData(a), buf + 1, (size_t)n * sizeof(double));
    return a;
}

/* Engine 2-D buffers are Fortran order, which is MATLAB's column-major
 * order, so the block can be copied straight across. */
static mxArray *push_real(const double *buf, int nrow, int ncol)
{
    mxArray *a = mxCreateDoubleMatrix((mwSize)nrow, (mwSize)ncol, mxREAL);
    size_t n = (size_t)nrow * (size_t)ncol;
    if (n > 0)
        memcpy(mxGetData(a), buf, n * sizeof(double));
    return a;
}

static mxArray *push_complex(const cd *buf, int nrow, int ncol)
{
    mxArray *a = mxCreateDoubleMatrix((mwSize)nrow, (mwSize)ncol, mxCOMPLEX);
    size_t n = (size_t)nrow * (size_t)ncol;
    size_t i;

#if MX_HAS_INTERLEAVED_COMPLEX
    double *d = (double *)mxGetData(a);
    for (i = 0; i < n; i++) {
        d[2 * i] = creal(buf[i]);
        d[2 * i + 1] = cimag(buf[i]);
    }
#else
    double *re = (double *)mxGetData(a);
    double *im = (double *)mxGetImagData(a);
    for (i = 0; i < n; i++) {
        re[i] = creal(buf[i]);
        im[i] = cimag(buf[i]);
    }
#endif
    return a;
}

static mxArray *push_int(const int *buf, int nrow, int ncol)
{
    mxArray *a = mxCreateNumericMatrix((mwSize)nrow, (mwSize)ncol,
                                       mxINT32_CLASS, mxREAL);
    size_t n = (size_t)nrow * (size_t)ncol;
    if (n > 0)
        memcpy(mxGetData(a), buf, n * sizeof(int));
    return a;
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

    kperpOUT = kparOUT = NULL;
    fOUT = ExOUT = EyOUT = EzOUT = BxOUT = ByOUT = BzOUT = SxOUT = SyOUT =
        SzOUT = NULL;
    EBOUT = VGPOUT = VGZOUT = SGPOUT = SGZOUT = uOUT = NULL;
    flagSolutionFoundOUT = flagTooHeavilyDampedOUT = flagNoConvergenceOUT = NULL;
}

/* ---- gateway ---------------------------------------------------------- */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int kp, kz;

    if (nrhs != NINPUT)
        mexErrMsgIdAndTxt("WHAMP:nInput",
                          "%d inputs required (%d given).", NINPUT, nrhs);
    if (nlhs > NOUTPUT)
        mexErrMsgIdAndTxt("WHAMP:nOutput",
                          "At most %d outputs allowed (%d requested).",
                          NOUTPUT, nlhs);

    /* ---- inputs, in the order matlab/+whamp/run.m passes them ---- */
    XC = scalar_in(prhs[0], "fceWHAMP");
    PZL = scalar_in(prhs[1], "pzlWHAMP");
    cycleZFirst = int_in(prhs[2], "zfirstWHAMP");
    vector_in(prhs[3], "nWHAMP", DN, NSPEC);
    vector_in(prhs[4], "tWHAMP", TA, NSPEC);
    vector_in(prhs[5], "dWHAMP", DD, NSPEC);
    species_ab_in(prhs[6], "aaWHAMP");
    vector_in(prhs[7], "assWHAMP", ASS, NSPEC);
    vector_in(prhs[8], "vdWHAMP", VD, NSPEC);
    vector_in(prhs[9], "pWHAMP", PM, 3);
    vector_in(prhs[10], "zWHAMP", ZM, 3);
    XOI = scalar_in(prhs[11], "fstartWHAMP");
    maxIterations = int_in(prhs[12], "maxIterationsWHAMP");

    check_axis(PM, "pWHAMP");
    check_axis(ZM, "zWHAMP");

    /* OUTPT/INOUT write to stdout and are not reachable from
     * whamp_engine(), but keep the flag off so nothing can print into a
     * MEX call. */
    printDebugInfo = 0;

    /* ---- solve ---- */
    whamp_engine();

    /* ---- outputs ---- */
    kp = kperpSize;
    kz = kparSize;

    if (nlhs > 0)
        plhs[0] = push_axis(kperpOUT, kp);
    if (nlhs > 1)
        plhs[1] = push_axis(kparOUT, kz);
    if (nlhs > 2)
        plhs[2] = push_complex(fOUT, kp, kz);
    if (nlhs > 3)
        plhs[3] = push_complex(ExOUT, kp, kz);
    if (nlhs > 4)
        plhs[4] = push_complex(EyOUT, kp, kz);
    if (nlhs > 5)
        plhs[5] = push_complex(EzOUT, kp, kz);
    if (nlhs > 6)
        plhs[6] = push_complex(BxOUT, kp, kz);
    if (nlhs > 7)
        plhs[7] = push_complex(ByOUT, kp, kz);
    if (nlhs > 8)
        plhs[8] = push_complex(BzOUT, kp, kz);
    if (nlhs > 9)
        plhs[9] = push_complex(SxOUT, kp, kz);
    if (nlhs > 10)
        plhs[10] = push_complex(SyOUT, kp, kz);
    if (nlhs > 11)
        plhs[11] = push_complex(SzOUT, kp, kz);
    if (nlhs > 12)
        plhs[12] = push_real(EBOUT, kp, kz);
    if (nlhs > 13)
        plhs[13] = push_real(VGPOUT, kp, kz);
    if (nlhs > 14)
        plhs[14] = push_real(VGZOUT, kp, kz);
    if (nlhs > 15)
        plhs[15] = push_real(SGPOUT, kp, kz);
    if (nlhs > 16)
        plhs[16] = push_real(SGZOUT, kp, kz);
    if (nlhs > 17)
        plhs[17] = push_real(uOUT, kp, kz);
    if (nlhs > 18)
        plhs[18] = push_int(flagSolutionFoundOUT, kp, kz);
    if (nlhs > 19)
        plhs[19] = push_int(flagTooHeavilyDampedOUT, kp, kz);
    if (nlhs > 20)
        plhs[20] = push_int(flagNoConvergenceOUT, kp, kz);

    release_output_buffers();
}
