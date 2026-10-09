/*
 * whamp_engine.c -- translation of SUBROUTINE WHAMP_ENGINE (src/whamp_engine.f90)
 *
 * The Fortran contained subroutines (plasma_setup, input_setup,
 * root_finding, allocate_output_matrices, save_output) access the host
 * variables through host association; they are translated as static
 * file-scope functions with the host variables as file-scope statics.
 */
#include "whamp.h"
#include "comin.h"
#include "comcout.h"
#include "comoutput.h"

/* Host variables of WHAMP_ENGINE (Fortran locals, re-created on each
 * call; zero-initialised here which is deterministic and matches the
 * Fortran behaviour for all defined uses). */
static int LoopI, IERR, IRK, J;
static int rootFindingConverged;
static int solutionIsTooHeavilyDamped;
static double REN[11]; /* particle mass expressed in masses of first particles */
static double RN;      /* mass of first particle in electron masses */
static double ADIR;    /* abs(D) */
static double DEK, DKP, DKZ, KV, PFQ, PLG, PLO, PO, PVO;
static double RED, TR, XA, XI, ZLG, ZLO, ZO, ZVO;
static cd XO, XVO, ddDX, OME, FPX, DOX, DOZ, DOP;
static double T[11], ST[11];
static double coef_poynt; /* coefficient used to simplify Poynting flux estimate */

static void plasma_setup(void)
{
    DEN = 0.;
    RED = 0.;
    for (J = 1; J <= 10; J++) {
        REN[J] = PROTON_MASS_RATIO * ASS[J];
        if (REN[J] == 0.)
            REN[J] = 1.;
        T[J] = TA[J] / TA[1];
        if (DN[J] == 0.)
            continue;
        JMA = J;
        RED = RED + DN[J] / REN[J];
        if (ASS[J] == 0.)
            DEN = DEN + DN[J];
    }

    RN = REN[1];
    /*                  ****  NORMALIZED TEMPERATURES AND VELOCITIES.  **** */
    for (J = 1; J <= JMA; J++) {
        REN[J] = REN[J] / RN;
        T[J] = T[J] * REN[J];
        ST[J] = sqrt(T[J]);
    }

    DEK = DEBYE_KM;
    PFQ = RED / DEK;
    PX = sqrt(PFQ);
    XA = XC / RN;
    TR = TA[1] / RN;
    CV = TR * (2. * ELECTRON_REST_ENERGY_KEV + TR) / ((ELECTRON_REST_ENERGY_KEV + TR) * (ELECTRON_REST_ENERGY_KEV + TR));
    CV = 1. / sqrt(CV);
    DEK = DEK * RN;
}

static void input_setup(void)
{
    /* empty in the Fortran source */
}

static void root_finding(void)
{
    cd CX; /* correction in Newton's iteration method */

    DIFU(2, JMA, &IERR);
    if (IERR != 0)
        solutionIsTooHeavilyDamped = 1;
    /*                  ****  START OF ITERATION.  **** */
    for (LoopI = 1; LoopI <= maxIterations; LoopI++) {
        ADIR = cabs(D);
        IRK = 1;
        CX = D / DX;
        /*CX=CX*X/(2*CX+X) ! finding zero of (w^2 D), faster convergence */
        for (;;) {
            X = X - CX;
            OME = (X * XA) * (X * XA);
            FPX = PFQ / OME;
            for (J = 1; J <= JMA; J++) {
                XP[J] = DN[J] / DEK / REN[J] / OME;
                XX[J] = X * REN[J];
            }
            DIFU(2, JMA, &IERR);
            if (IERR != 0) {
                solutionIsTooHeavilyDamped = 1;
                goto end_loop_iteration;
            }
            if (cabs(D) < ADIR) {
                if ((cabs(CX) <= FREQ_REL_TOL * cabs(X)) /* relative frequency precision */
                    || (cabs(CX) < FREQ_ABS_TOL)) {       /* absolute precision */
                    rootFindingConverged = 1;
                    if (LoopI >= 2) { /* at least 2 steps have been made */
                        goto end_loop_iteration;
                    } else {
                        goto cycle_loop_iteration;
                    }
                } else {
                    goto cycle_loop_iteration;
                }
            } else {
                X = X + CX;
                CX = CX / 2.;
                if (IRK > 3) { /* check for sitting at local minima */
                    if ((cabs(D) - ADIR) / ADIR < MINIMA_TOL) {
                        /*PRINT*,' Local minima!' */
                        rootFindingConverged = 0;
                        goto end_loop_iteration;
                    }
                }
                IRK = IRK + 1;
                if (IRK > maxIterations)
                    goto end_loop_iteration;
            }
        }
    cycle_loop_iteration:;
    }
end_loop_iteration:;
}

/* Number of grid points spanned by (start, stop, step), matching the Fortran
 * expression 1 + floor((max(start,stop) - min(start,stop)) * sign(1,step) /
 * step).  A single-valued axis (start == stop) yields one point. */
static int axis_size(double start, double stop, double step)
{
    double lo, hi;
    if (start == stop)
        return 1; /* single-valued axis: the step is irrelevant */
    lo = (start < stop) ? start : stop;
    hi = (start > stop) ? start : stop;
    return 1 + (int)floor((hi - lo) * copysign(1.0, step) / step);
}

/* Fill a freshly allocated 1-based axis with start + (i-1)*step. */
static void fill_axis(double *axis, int size, double start, double step)
{
    for (int i = 1; i <= size; i++)
        axis[i] = start + (i - 1.0) * step;
}

static void allocate_output_matrices(void)
{
    /* Free any previous allocation and (re)allocate the k-perp axis. */
    ALLOC1(kperpOUT, kperpSize, double);
    kperpSize = axis_size(PM[1], PM[2], PM[3]);
    ALLOC1(kperpOUT, kperpSize, double);
    fill_axis(kperpOUT, kperpSize, PM[1], PM[3]);

    ALLOC1(kparOUT, kparSize, double);
    kparSize = axis_size(ZM[1], ZM[2], ZM[3]);
    ALLOC1(kparOUT, kparSize, double);
    fill_axis(kparOUT, kparSize, ZM[1], ZM[3]);

    /* Two-dimensional (kperp x kpar) result arrays, Fortran order.  calloc
     * zeroes them, matching the Fortran "= 0" initialisation. */
    ALLOC2(fOUT, cd);
    ALLOC2(ExOUT, cd);
    ALLOC2(EyOUT, cd);
    ALLOC2(EzOUT, cd);
    ALLOC2(BxOUT, cd);
    ALLOC2(ByOUT, cd);
    ALLOC2(BzOUT, cd);
    ALLOC2(SxOUT, cd);
    ALLOC2(SyOUT, cd);
    ALLOC2(SzOUT, cd);
    ALLOC2(EBOUT, double);
    ALLOC2(VGPOUT, double);
    ALLOC2(VGZOUT, double);
    ALLOC2(SGPOUT, double);
    ALLOC2(SGZOUT, double);
    ALLOC2(uOUT, double);
    ALLOC2(flagSolutionFoundOUT, int);
    ALLOC2(flagTooHeavilyDampedOUT, int);
    ALLOC2(flagNoConvergenceOUT, int);
}

static void save_output(void)
{
    int indexKperp, indexKpar;
    if (PLG == PM[1]) {
        indexKperp = 1;
    } else {
        indexKperp = 1 + fnint((PLG - PM[1]) / PM[3]);
    }
    if (ZLG == ZM[1]) {
        indexKpar = 1;
    } else {
        indexKpar = 1 + fnint((ZLG - ZM[1]) / ZM[3]);
    }

    if (rootFindingConverged) {
        F2D(flagSolutionFoundOUT, indexKperp, indexKpar) = 1;
        F2D(fOUT, indexKperp, indexKpar) = X;
        F2D(ExOUT, indexKperp, indexKpar) = EFL[1];
        F2D(EyOUT, indexKperp, indexKpar) = EFL[2];
        F2D(EzOUT, indexKperp, indexKpar) = EFL[3];
        F2D(BxOUT, indexKperp, indexKpar) = BFL[1];
        F2D(ByOUT, indexKperp, indexKpar) = BFL[2];
        F2D(BzOUT, indexKperp, indexKpar) = BFL[3];
        /*     write Poynting vector uW/m^2 */
        coef_poynt = POYNTING_COEF;
        F2D(SxOUT, indexKperp, indexKpar) =
            creal(EFL[2] * conjg(BFL[3]) - EFL[3] * conjg(BFL[2])) * coef_poynt;
        F2D(SyOUT, indexKperp, indexKpar) =
            creal(EFL[3] * conjg(BFL[1]) - EFL[1] * conjg(BFL[3])) * coef_poynt;
        F2D(SzOUT, indexKperp, indexKpar) =
            creal(EFL[1] * conjg(BFL[2]) - EFL[2] * conjg(BFL[1])) * coef_poynt;
        F2D(EBOUT, indexKperp, indexKpar) = sqrt(
            (creal(EFL[1] * conjg(EFL[1]) + EFL[2] * conjg(EFL[2]) +
                   EFL[3] * conjg(EFL[3]))) /
            (creal(BFL[1] * conjg(BFL[1]) + BFL[2] * conjg(BFL[2]) +
                   BFL[3] * conjg(BFL[3]))));
        F2D(VGPOUT, indexKperp, indexKpar) = VG[1];
        F2D(VGZOUT, indexKperp, indexKpar) = VG[2];
        F2D(SGPOUT, indexKperp, indexKpar) = SG[1];
        F2D(SGZOUT, indexKperp, indexKpar) = SG[2];
        F2D(uOUT, indexKperp, indexKpar) = ENE;
    } else if (solutionIsTooHeavilyDamped) {
        F2D(flagTooHeavilyDampedOUT, indexKperp, indexKpar) = 1;
    } else {
        F2D(flagNoConvergenceOUT, indexKperp, indexKpar) = 1;
    }
}

void whamp_engine(void)
{
    IERR = 0;
    allocate_output_matrices();
    plasma_setup();
    if ((PM[1] == 0.0) && (ZM[1] == 0.0))
        return;
    input_setup();

    KV = 1;
    PLG = PM[1];
    if (PM[3] < 0.)
        PLG = PM[2];
    ZLG = ZM[1];
    if (ZM[3] < 0.)
        ZLG = ZM[2];
    if (PZL == 1.) {
        P = pow(10., PLG);
        Z = pow(10., ZLG);
    } else {
        P = PLG;
        Z = ZLG;
    }
    X = XOI;
    for (;;) { /* loop_z_p */
        OME = (X * XA) * (X * XA);
        FPX = PFQ / OME;
        for (J = 1; J <= JMA; J++) {
            XX[J] = X * REN[J];
            PP[J] = P * ST[J];
            ZZ[J] = Z * ST[J];
            XP[J] = DN[J] / DEK / REN[J] / OME;
        }

        solutionIsTooHeavilyDamped = 0; /* default not heavily damped */
        rootFindingConverged = 0;       /* default no convergence */
        root_finding();

        if (rootFindingConverged && !solutionIsTooHeavilyDamped) {
            /*                  ****  CONVERGENCE!  **** */
            DIFU(4, JMA, &IERR);

            XI = cimag(X);
            VG[1] = -creal(DP / DX);
            VG[2] = -creal(DZ / DX);
            RI = sqrt(P * P + Z * Z) * CV / X; /* refractive index */
            if (VG[1] != 0.)
                SG[1] = XI / VG[1];
            if (VG[2] != 0.)
                SG[2] = XI / VG[2];
            /*          ****  PRINT THE RESULTS.  **** */
            /*if (printDebugInfo)              CALL OUTPT */
            save_output();
            PO = P;
            ZO = Z;
            XO = X;
            if (KV != 0) {
                XVO = X;
                ZVO = Z;
                ZLO = ZLG;
                PVO = P;
                PLO = PLG;
                DOX = DX;
                DOZ = DZ;
                DOP = DP;
                KV = 0;
            }
        } else {
            save_output();
        }
        if (!rootFindingConverged) {
            /*if (printDebugInfo) PRINT 125,P,Z,X,LoopI,IRK */
            if (cycleZFirst == 1)
                PLG = SWEEP_END; /* end cycling in P */
            if (cycleZFirst == 2)
                ZLG = SWEEP_END; /* end cycling in Z */
        }
        if (solutionIsTooHeavilyDamped) {
            /*if (printDebugInfo) PRINT*,' TOO HEAVILY DAMPED!' */
            IERR = 0;
            /*if (printDebugInfo) CALL OUTPT */
            if (cycleZFirst == 1)
                PLG = SWEEP_END;
            if (cycleZFirst == 2)
                ZLG = SWEEP_END;
        }
        if (cycleZFirst == 0) { /* cycle first P */
            PLG = PLG + PM[3];
            /*                   ****  UPDATE P AND Z.  **** */
            if (PLG >= PM[1] && PLG <= PM[2]) {
                P = PLG + PZL * (pow(10., PLG) - PLG);
            } else {
                ZLG = ZLG + ZM[3];
                if (ZLG < ZM[1] || ZLG > ZM[2]) {
                    break; /* exit loop_z_p */
                }
                KV = 1;
                PLG = PLO;
                P = PVO;
                Z = ZLG + PZL * (pow(10., ZLG) - ZLG);
            }
        } else if (cycleZFirst == 1) { /* cycle first Z */
            ZLG = ZLG + ZM[3];
            if (ZLG >= ZM[1] && ZLG <= ZM[2]) {
                Z = ZLG + PZL * (pow(10., ZLG) - ZLG);
            } else {
                PLG = PLG + PM[3];
                if (PLG < PM[1] || PLG > PM[2])
                    break; /* exit loop_z_p */
                KV = 1;
                ZLG = ZLO;
                Z = ZVO;
                P = PLG + PZL * (pow(10., PLG) - PLG);
            }
        }
        /*                    ****  NEW START FREQUENCY.  **** */
        if (KV != 0) {
            DKP = P - PVO;
            DKZ = Z - ZVO;
            ddDX = (DKP * DOP + DKZ * DOZ) / DOX;
            X = XVO - ddDX;
        } else {
            DKP = P - PO;
            DKZ = Z - ZO;
            ddDX = (DKP * DP + DKZ * DZ) / DX;
            X = XO - ddDX;
        }
    }
    return;
}
