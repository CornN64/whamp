/*
 * whamp_track_dump.c -- dump whamp_track_v2() output for the equivalence test.
 *
 * Runs one named case (see CASES below) and writes every output array in a
 * machine-readable form that srcm/whamp_track_v2_vs_c.m parses and compares
 * against the MATLAB implementation.  Cases and their plasma models mirror the
 * MATLAB side of that script exactly.
 *
 * Usage:  whamp_track_dump <caseName>
 */
#include "whamp_track.h"
#include "comin.h"
#include "comcout.h"
#include "comoutput.h"
#include "comtrack.h"

#include <string.h>

/* One species, in the same order PlasmaModel.Species uses. */
typedef struct {
    double m, n, t, a, b, d, vd;
} DumpSpecies;

typedef struct {
    const char *name;
    double xc;              /* = 0.0279928 * B(nT) */
    DumpSpecies sp[3];
    int nspecies;
    double pm1, pm2, pm3;   /* kperp  [start stop step] */
    double zm1, zm2, zm3;   /* kpar   [start stop step] */
    int cycleZFirst;
    int useLog;
    double xoi;
} DumpCase;

/* Models/Ex1: O+ beam on electrons, XC 2.8 kHz. */
static const DumpCase CASES[] = {
    { "line_kpar", 2.8,
      { {16,1,10,5,0,1,1}, {0,1,1,1,0,1,0} }, 2,
      0.0, 0.0, 10.0,   0.002, 0.008, 0.001, 1, 0, 0.1 },
    { "line_kperp", 2.8,
      { {16,1,10,5,0,1,1}, {0,1,1,1,0,1,0} }, 2,
      0.005, 0.010, 0.001,   0.0022, 0.0022, 10.0, 1, 0, 0.1 },
    { "line_kpar_pfirst", 2.8,
      { {16,1,10,5,0,1,1}, {0,1,1,1,0,1,0} }, 2,
      0.0, 0.0, 10.0,   0.002, 0.008, 0.001, 0, 0, 0.1 },
    { "grid2d_zfirst", 2.8,
      { {16,1,10,5,0,1,1}, {0,1,1,1,0,1,0} }, 2,
      0.005, 0.008, 0.001,   0.002, 0.006, 0.001, 1, 0, 0.1 },
    { "grid2d_pfirst", 2.8,
      { {16,1,10,5,0,1,1}, {0,1,1,1,0,1,0} }, 2,
      0.005, 0.008, 0.001,   0.002, 0.006, 0.001, 0, 0, 0.1 },
    { "grid2d_negstep", 2.8,
      { {16,1,10,5,0,1,1}, {0,1,1,1,0,1,0} }, 2,
      0.008, 0.005, -0.001,   0.006, 0.002, -0.001, 1, 0, 0.1 },
    { "grid2d_log", 2.8,
      { {16,1,10,5,0,1,1}, {0,1,1,1,0,1,0} }, 2,
      -2.301, -2.096, 0.1,   -2.699, -2.522, 0.1, 1, 1, 0.1 },
    { "ex2_grid", 0.008978,
      { {0,1,10,1,0,1,0}, {1,1,1,1,0,1,0} }, 2,
      0.0, 0.02, 0.01,   0.001, 0.005, 0.002, 1, 0, 0.1 },
    { "ex2_log", 0.008978,
      { {0,1,10,1,0,1,0}, {1,1,1,1,0,1,0} }, 2,
      -2.3, -1.8, 0.1,   -3.0, -2.7, 0.1, 1, 1, 0.1 },
    { "isotropic", 2.8,
      { {16,1,10,1,0,1,0}, {0,1,1,1,0,1,0} }, 2,
      0.0, 0.0, 10.0,   0.002, 0.006, 0.001, 1, 0, 0.1 },
    { "three_species", 2.8,
      { {16,1,10,5,0,1,1}, {1,0.5,1,1,0,1,0.5}, {0,1.5,1,1,0,1,0} }, 3,
      0.005, 0.008, 0.001,   0.002, 0.006, 0.001, 1, 0, 0.1 },
    { "start_high", 2.8,
      { {16,1,10,5,0,1,1}, {0,1,1,1,0,1,0} }, 2,
      0.0, 0.0, 10.0,   0.002, 0.008, 0.001, 1, 0, 0.9 },
    /* Whistlereanisotropy_.m: 50 nT electrons, kperp = 0.01, kpar 0.01:0.01:0.90.
     * The root is strongly complex along this line, so it is the case that
     * checks the seed is carried between points as a COMPLEX frequency. */
    { "line89_whistler", 1.39964,
      { {0,1,100,5,0,1,0} }, 1,
      0.01, 0.01, 10.0,   0.01, 0.90, 0.01, 1, 0, 0.006208543881328393 },
    /* Same plasma one line further out in kperp, where the sweep crosses the
     * boundary between "too heavily damped" and "no convergence" and the
     * Newton iteration raises the damping flag part-way through a solve. */
    { "damped_boundary", 1.39964,
      { {0,1,100,5,0,1,0} }, 1,
      0.25, 0.25, 10.0,   0.01, 0.90, 0.01, 1, 0, 0.006208543881328393 }
};

static void set_model(const DumpCase *c)
{
    int j, k;
    for (j = 1; j <= 10; ++j) {
        DN[j] = 0.0;
        TA[j] = 0.0;
        DD[j] = 0.0;
        VD[j] = 0.0;
        ASS[j] = 0.0;
        AA[j][1] = 0.0;
        AA[j][2] = 0.0;
    }
    for (k = 0; k < c->nspecies; ++k) {
        j = k + 1;
        ASS[j] = c->sp[k].m;
        DN[j] = c->sp[k].n * 1.0e6;      /* cm^-3 -> m^-3 */
        TA[j] = c->sp[k].t / 1.0e3;      /* eV    -> keV  */
        AA[j][1] = c->sp[k].a;
        AA[j][2] = c->sp[k].b;
        DD[j] = c->sp[k].d;
        VD[j] = c->sp[k].vd;
    }
    XC = c->xc;
    PM[1] = c->pm1;
    PM[2] = c->pm2;
    PM[3] = c->pm3;
    ZM[1] = c->zm1;
    ZM[2] = c->zm2;
    ZM[3] = c->zm3;
    PZL = c->useLog ? 1.0 : 0.0;
    cycleZFirst = c->cycleZFirst;
    XOI = c->xoi;
    maxIterations = 50;
}

static double axis_value(double raw, int useLog)
{
    return useLog ? pow(10.0, raw) : raw;
}

static void print_complex(cd value)
{
    printf(" %.17g %.17g", creal(value), cimag(value));
}

static void print_real(double value)
{
    printf(" %.17g", value);
}

int main(int argc, char **argv)
{
    TrackDiagnostics diagnostics;
    const DumpCase *c = NULL;
    int i, ip, iz;

    if (argc < 2) {
        fprintf(stderr, "usage: %s <caseName>\n", argv[0]);
        return 2;
    }
    for (i = 0; i < (int)(sizeof(CASES) / sizeof(CASES[0])); ++i) {
        if (strcmp(argv[1], CASES[i].name) == 0) {
            c = &CASES[i];
            break;
        }
    }
    if (!c) {
        fprintf(stderr, "unknown case '%s'\n", argv[1]);
        return 2;
    }

    set_model(c);
    track_options_default();
    whamp_track_v2(&diagnostics);

    printf("CASE %s\n", c->name);
    printf("SIZE %d %d\n", kperpSize, kparSize);
    for (ip = 1; ip <= kperpSize; ++ip)
        printf("KPERP %d %.17g\n", ip, axis_value(kperpOUT[ip], c->useLog));
    for (iz = 1; iz <= kparSize; ++iz)
        printf("KPAR %d %.17g\n", iz, axis_value(kparOUT[iz], c->useLog));

    for (iz = 1; iz <= kparSize; ++iz) {
        for (ip = 1; ip <= kperpSize; ++ip) {
            int cell = (ip - 1) + (iz - 1) * kperpSize;
            printf("CELL %d %d", ip, iz);
            print_complex(F2D(fOUT, ip, iz));
            print_complex(F2D(ExOUT, ip, iz));
            print_complex(F2D(EyOUT, ip, iz));
            print_complex(F2D(EzOUT, ip, iz));
            print_complex(F2D(BxOUT, ip, iz));
            print_complex(F2D(ByOUT, ip, iz));
            print_complex(F2D(BzOUT, ip, iz));
            print_complex(F2D(SxOUT, ip, iz));
            print_complex(F2D(SyOUT, ip, iz));
            print_complex(F2D(SzOUT, ip, iz));
            print_real(F2D(EBOUT, ip, iz));
            print_real(F2D(VGPOUT, ip, iz));
            print_real(F2D(VGZOUT, ip, iz));
            print_real(F2D(SGPOUT, ip, iz));
            print_real(F2D(SGZOUT, ip, iz));
            print_real(F2D(uOUT, ip, iz));
            printf(" %d %d %d %d %.17g %d %d",
                F2D(flagSolutionFoundOUT, ip, iz),
                F2D(flagTooHeavilyDampedOUT, ip, iz),
                F2D(flagNoConvergenceOUT, ip, iz),
                F2D(flagSuspectBranchOUT, ip, iz),
                certWindingOUT[cell],
                F2D(flagNotCertifiedOUT, ip, iz),
                F2D(flagCertUnknownOUT, ip, iz));
            print_complex(predictedOUT[cell]);
            printf("\n");
        }
    }

    printf("DIAG %d %d %d %d %d %d %d %d %d %.17g\n",
        diagnostics.pointsVisited, diagnostics.converged,
        diagnostics.rejected, diagnostics.clampedSteps,
        diagnostics.rescued, diagnostics.flagged, diagnostics.certified,
        diagnostics.notCertified, diagnostics.certUnknown,
        diagnostics.meanAttempts);
    return 0;
}
