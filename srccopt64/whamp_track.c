/*
 * whamp_track.c -- Version 2 branch continuation for the C WHAMP engine.
 *
 * The tracker shares the dispersion kernel and output arrays with
 * whamp_engine(), but keeps its own root map and sweep history.
 */
#include "whamp_track.h"
#include "comin.h"
#include "comcout.h"
#include "comoutput.h"
#include "comtrack.h"

#include <float.h>
#include <limits.h>

typedef struct {
    int valid;
    double p, z;
    cd x, dp, dz, dx;
} TrackRoot;

typedef struct {
    int status;
    int clamped;
    cd x, dp, dz, dx;
} TrackSolution;

typedef struct {
    int converged, rejected, clamped, rescued, flagged;
    int certified, notCertified, certUnknown, pointsVisited;
    double attempts;
} TrackCounts;

TrackOptions trackOptions = {
    3.0, 6.0, {0.25, -0.25, 0.5, -0.5, 1.0, -1.0}, 6, 1, 1, 0.0, 8
};

static double trRen[11], trTemp[11], trSt[11];
static double trDek, trPfq, trXa;

static int tr_axis_size(double a, double b, double step)
{
    double lo, hi;
    if (a == b)
        return 1;
    lo = fmin(a, b);
    hi = fmax(a, b);
    return 1 + (int)floor((hi - lo) * copysign(1.0, step) / step);
}

static double tr_origin(double a, double b, double step)
{
    return step < 0.0 ? b : a;
}

static double tr_physical(double value)
{
    return PZL == 1.0 ? pow(10.0, value) : value;
}

static double tr_grid_value(double origin, double step, int index)
{
    return origin + (index - 1.0) * step;
}

static int tr_grid_index(double value, double origin, double step, int size)
{
    int index;
    if (step == 0.0)
        return 1;
    index = 1 + fnint((value - origin) / step);
    if (index < 1)
        return 1;
    if (index > size)
        return size;
    return index;
}

static void tr_plasma_setup(void)
{
    double red = 0.0, tr, rn;
    int j, last = 0;

    DEN = 0.0;
    for (j = 1; j <= 10; ++j) {
        trRen[j] = PROTON_MASS_RATIO * ASS[j];
        if (trRen[j] == 0.0)
            trRen[j] = 1.0;
        trTemp[j] = TA[j] / TA[1];
        if (DN[j] == 0.0)
            continue;
        last = j;
        red += DN[j] / trRen[j];
        if (ASS[j] == 0.0)
            DEN += DN[j];
    }
    JMA = last;
    rn = trRen[1];
    tr = TA[1] / rn;
    for (j = 1; j <= JMA; ++j) {
        trRen[j] /= rn;
        trTemp[j] *= trRen[j];
        trSt[j] = sqrt(trTemp[j]);
    }
    trDek = DEBYE_KM * rn;
    trPfq = red / DEBYE_KM;
    PX = sqrt(trPfq);
    trXa = XC / rn;
    CV = 1.0 / sqrt(tr * (2.0 * ELECTRON_REST_ENERGY_KEV + tr) /
                    ((ELECTRON_REST_ENERGY_KEV + tr) *
                     (ELECTRON_REST_ENERGY_KEV + tr)));
}

static void tr_set_state(double p, double z, cd x)
{
    cd omega = (x * trXa) * (x * trXa);
    int j;
    P = p;
    Z = z;
    X = x;
    for (j = 1; j <= JMA; ++j) {
        XX[j] = x * trRen[j];
        PP[j] = p * trSt[j];
        ZZ[j] = z * trSt[j];
        XP[j] = DN[j] / trDek / trRen[j] / omega;
    }
}

static int tr_evaluate(double p, double z, cd x, int order)
{
    int ierr = 0;
    tr_set_state(p, z, x);
    DIFU(order, JMA, &ierr);
    return ierr;
}

static void tr_store_solution(TrackSolution *sol, double p, double z)
{
    if (tr_evaluate(p, z, sol->x, 4) != 0) {
        sol->status = 2;
        return;
    }
    sol->x = X;
    sol->dx = DX;
    sol->dp = DP;
    sol->dz = DZ;
    sol->status = 1;
}

static int tr_limit_step(cd *correction, cd x, cd center, double radius)
{
    double budget, magnitude;
    if (!isfinite(radius))
        return 0;
    budget = radius - cabs(x - center);
    if (budget <= 0.0) {
        *correction = 0.0;
        return 2;
    }
    magnitude = cabs(*correction);
    if (magnitude > budget && magnitude > 0.0) {
        *correction *= budget / magnitude;
        return 1;
    }
    return 0;
}

static TrackSolution tr_solve(double p, double z, cd seed, double radius,
                              cd center)
{
    TrackSolution sol = {0, 0, seed, 0.0, 0.0, 0.0};
    cd correction;
    double oldMagnitude;
    int iteration, irk, ierr;

    if (tr_evaluate(p, z, seed, 2) != 0) {
        sol.status = 2;
        return sol;
    }
    for (iteration = 1; iteration <= maxIterations; ++iteration) {
        oldMagnitude = cabs(D);
        irk = 1;
        correction = D / DX;
        for (;;) {
            int limit = tr_limit_step(&correction, X, center, radius);
            int wasClamped = limit != 0;
            if (limit == 2) {
                ++sol.clamped;
                goto done;
            }
            if (wasClamped)
                ++sol.clamped;
            X -= correction;
            ierr = tr_evaluate(p, z, X, 2);
            /* A damping flag raised part-way through the iteration is not the
             * same statement as one raised by the very first evaluation: the
             * reference reports the latter as "too heavily damped" but the
             * former only as a search that failed, because its solver routine
             * resets the status to "no root found" whenever it returns
             * without one. */
            if (ierr != 0) {
                sol.status = 0;
                return sol;
            }
            if (cabs(D) < oldMagnitude) {
                if (!wasClamped &&
                    (cabs(correction) <= FREQ_REL_TOL * cabs(X) ||
                     cabs(correction) < FREQ_ABS_TOL)) {
                    if (iteration >= 2) {
                        sol.status = 1;
                        goto done;
                    }
                }
                break;
            }
            X += correction;
            correction *= 0.5;
            if (irk > 3 &&
                (cabs(D) - oldMagnitude) / oldMagnitude < MINIMA_TOL)
                goto done;
            ++irk;
            if (irk > maxIterations)
                goto done;
        }
    }
done:
    sol.x = X;
    return sol;
}

static void tr_save_root(int ip, int iz, const TrackSolution *sol)
{
    double poynting = POYNTING_COEF;
    double xi = cimag(sol->x);
    double b2, e2;

    X = sol->x;
    DX = sol->dx;
    DP = sol->dp;
    DZ = sol->dz;
    VG[1] = -creal(DP / DX);
    VG[2] = -creal(DZ / DX);
    if (VG[1] != 0.0)
        SG[1] = xi / VG[1];
    else
        SG[1] = 0.0;
    if (VG[2] != 0.0)
        SG[2] = xi / VG[2];
    else
        SG[2] = 0.0;

    F2D(flagSolutionFoundOUT, ip, iz) = 1;
    F2D(fOUT, ip, iz) = sol->x;
    F2D(ExOUT, ip, iz) = EFL[1];
    F2D(EyOUT, ip, iz) = EFL[2];
    F2D(EzOUT, ip, iz) = EFL[3];
    F2D(BxOUT, ip, iz) = BFL[1];
    F2D(ByOUT, ip, iz) = BFL[2];
    F2D(BzOUT, ip, iz) = BFL[3];
    F2D(SxOUT, ip, iz) =
        creal(EFL[2] * conj(BFL[3]) - EFL[3] * conj(BFL[2])) * poynting;
    F2D(SyOUT, ip, iz) =
        creal(EFL[3] * conj(BFL[1]) - EFL[1] * conj(BFL[3])) * poynting;
    F2D(SzOUT, ip, iz) =
        creal(EFL[1] * conj(BFL[2]) - EFL[2] * conj(BFL[1])) * poynting;
    e2 = creal(EFL[1] * conj(EFL[1]) + EFL[2] * conj(EFL[2]) +
               EFL[3] * conj(EFL[3]));
    b2 = creal(BFL[1] * conj(BFL[1]) + BFL[2] * conj(BFL[2]) +
               BFL[3] * conj(BFL[3]));
    F2D(EBOUT, ip, iz) = sqrt(e2 / b2);
    F2D(VGPOUT, ip, iz) = VG[1];
    F2D(VGZOUT, ip, iz) = VG[2];
    F2D(SGPOUT, ip, iz) = SG[1];
    F2D(SGZOUT, ip, iz) = SG[2];
    F2D(uOUT, ip, iz) = ENE;
}

static double tr_winding(double p, double z, cd center, double radius,
                         int samples)
{
    int refine, i;
    double winding = NAN;
    for (refine = 0; refine <= 2; ++refine) {
        cd first = 0.0, previous = 0.0;
        double phase = 0.0, worst = 0.0;
        int valid = 1;
        for (i = 0; i <= samples; ++i) {
            double theta = 2.0 * PI * i / samples;
            cd point = center + radius * (cos(theta) + I * sin(theta));
            cd value;
            if (tr_evaluate(p, z, point, 1) != 0) {
                valid = 0;
                break;
            }
            value = D;
            if (!isfinite(creal(value)) || !isfinite(cimag(value)) ||
                cabs(value) == 0.0) {
                valid = 0;
                break;
            }
            if (i == 0) {
                first = value;
            } else {
                double step = carg(value * conj(previous));
                phase += step;
                worst = fmax(worst, fabs(step));
            }
            previous = value;
        }
        if (!valid)
            return NAN;
        (void)first;
        winding = nearbyint((phase / (2.0 * PI)) * 2.0) / 2.0;
        if (worst <= PI / 2.0)
            return winding;
        samples *= 2;
    }
    return winding;
}

static int tr_cell(int ip, int iz)
{
    return ip - 1 + (iz - 1) * kperpSize;
}

static TrackRoot *tr_predecessor(TrackRoot *roots, int ip, int iz,
                                 int zFirst, int *axis, int *aip, int *aiz)
{
    int i;
    *axis = 0;
    if (zFirst) {
        if (iz > 1 && roots[tr_cell(ip, iz - 1)].valid) {
            *axis = 1; *aip = ip; *aiz = iz - 1;
            return &roots[tr_cell(ip, iz - 1)];
        }
        /* MATLAB walks back at most 64 cells along the outer axis. */
        for (i = ip - 1; i >= 1 && ip - i <= 64; --i) {
            if (roots[tr_cell(i, iz)].valid) {
                *axis = 2; *aip = i; *aiz = iz;
                return &roots[tr_cell(i, iz)];
            }
        }
    } else {
        if (ip > 1 && roots[tr_cell(ip - 1, iz)].valid) {
            *axis = 1; *aip = ip - 1; *aiz = iz;
            return &roots[tr_cell(ip - 1, iz)];
        }
        /* MATLAB walks back at most 64 cells along the outer axis. */
        for (i = iz - 1; i >= 1 && iz - i <= 64; --i) {
            if (roots[tr_cell(ip, i)].valid) {
                *axis = 2; *aip = ip; *aiz = i;
                return &roots[tr_cell(ip, i)];
            }
        }
    }
    return NULL;
}

static double tr_nearest_rate(const double *rates, int axis, int ip, int iz)
{
    int i, j, best = INT_MAX, found = 0;
    double rate = 0.0;
    if (axis < 1 || axis > 2)
        return NAN;
    for (j = 1; j <= kparSize; ++j) {
        for (i = 1; i <= kperpSize; ++i) {
            int distance;
            double value = rates[2 * tr_cell(i, j) + axis - 1];
            if (!isfinite(value))
                continue;
            distance = abs(i - ip) + abs(j - iz);
            if (distance < best) {
                best = distance;
                rate = value;
                found = 1;
            } else if (distance == best) {
                rate = fmax(rate, value);
            }
        }
    }
    return found ? rate : NAN;
}

static cd tr_predict(const TrackRoot *anchor, double p, double z,
                     int *informative)
{
    double dp = p - anchor->p, dz = z - anchor->z;
    if (anchor->dx == 0.0) {
        *informative = 0;
        return anchor->x;
    }
    *informative = !((dp != 0.0 && anchor->dp == 0.0) ||
                     (dz != 0.0 && anchor->dz == 0.0));
    return anchor->x - (dp * anchor->dp + dz * anchor->dz) / anchor->dx;
}

static void tr_record_solution(int cell, TrackRoot *roots,
                              const TrackSolution *sol, double p, double z)
{
    roots[cell].valid = 1;
    roots[cell].p = p;
    roots[cell].z = z;
    roots[cell].x = sol->x;
    roots[cell].dp = sol->dp;
    roots[cell].dz = sol->dz;
    roots[cell].dx = sol->dx;
}

static void tr_allocate_outputs(void)
{
    size_t n;
    int i;
    kperpSize = tr_axis_size(PM[1], PM[2], PM[3]);
    kparSize = tr_axis_size(ZM[1], ZM[2], ZM[3]);
    ALLOC1(kperpOUT, kperpSize, double);
    ALLOC1(kparOUT, kparSize, double);
    for (i = 1; i <= kperpSize; ++i)
        kperpOUT[i] = tr_grid_value(tr_origin(PM[1], PM[2], PM[3]), PM[3], i);
    for (i = 1; i <= kparSize; ++i)
        kparOUT[i] = tr_grid_value(tr_origin(ZM[1], ZM[2], ZM[3]), ZM[3], i);
    ALLOC2(fOUT, cd); ALLOC2(ExOUT, cd); ALLOC2(EyOUT, cd); ALLOC2(EzOUT, cd);
    ALLOC2(BxOUT, cd); ALLOC2(ByOUT, cd); ALLOC2(BzOUT, cd);
    ALLOC2(SxOUT, cd); ALLOC2(SyOUT, cd); ALLOC2(SzOUT, cd);
    ALLOC2(EBOUT, double); ALLOC2(VGPOUT, double); ALLOC2(VGZOUT, double);
    ALLOC2(SGPOUT, double); ALLOC2(SGZOUT, double); ALLOC2(uOUT, double);
    ALLOC2(flagSolutionFoundOUT, int);
    ALLOC2(flagTooHeavilyDampedOUT, int);
    ALLOC2(flagNoConvergenceOUT, int);
    track_allocate_output();
    n = (size_t)kperpSize * (size_t)kparSize;
    for (i = 0; (size_t)i < n; ++i)
        predictedOUT[i] = NAN + I * NAN;
}

void track_options_default(void)
{
    trackOptions.trustFactor = 3.0;
    trackOptions.jumpFactor = 6.0;
    trackOptions.scanOffsets[0] = 0.25;
    trackOptions.scanOffsets[1] = -0.25;
    trackOptions.scanOffsets[2] = 0.5;
    trackOptions.scanOffsets[3] = -0.5;
    trackOptions.scanOffsets[4] = 1.0;
    trackOptions.scanOffsets[5] = -1.0;
    trackOptions.numScanOffsets = 6;
    trackOptions.useScan = 1;
    trackOptions.certify = 1;
    trackOptions.certRadius = 0.0;
    trackOptions.certSamples = 8;
}

void whamp_track_v2(TrackDiagnostics *diagnostics)
{
    TrackRoot *roots = NULL, last = {0}, lineAnchor = {0};
    double *rates = NULL;
    TrackCounts counts = {0};
    double pOrigin, zOrigin, plg, zlg, p, z;
    cd xseed;
    int ip, iz, zFirst, lineFirst = 1, done = 0;
    int maxRateCells, i;

    if (trackOptions.certSamples < 4)
        trackOptions.certSamples = 8;
    if (trackOptions.numScanOffsets < 0)
        trackOptions.numScanOffsets = 0;
    if (trackOptions.numScanOffsets > 6)
        trackOptions.numScanOffsets = 6;
    if (maxIterations < 1)
        maxIterations = 1;

    tr_allocate_outputs();
    tr_plasma_setup();

    maxRateCells = kperpSize * kparSize;
    roots = (TrackRoot *)calloc((size_t)maxRateCells, sizeof(*roots));
    rates = (double *)malloc((size_t)maxRateCells * 2 * sizeof(*rates));
    if (diagnostics) {
        memset(diagnostics, 0, sizeof(*diagnostics));
        /* The sweep visits each grid cell at most once, so the grid bounds
         * the number of per-point records. */
        diagnostics->attemptsPerPoint =
            (double *)malloc((size_t)maxRateCells * sizeof(double));
    }
    if (!roots || !rates) {
        free(roots);
        free(rates);
        fprintf(stderr, "whamp_track_v2: unable to allocate sweep history\n");
        if (diagnostics) {
            free(diagnostics->attemptsPerPoint);
            memset(diagnostics, 0, sizeof(*diagnostics));
        }
        return;
    }
    for (i = 0; i < maxRateCells * 2; ++i)
        rates[i] = NAN;

    pOrigin = tr_origin(PM[1], PM[2], PM[3]);
    zOrigin = tr_origin(ZM[1], ZM[2], ZM[3]);
    plg = pOrigin;
    zlg = zOrigin;
    p = tr_physical(plg);
    z = tr_physical(zlg);
    xseed = XOI;
    zFirst = cycleZFirst != 0;

    while (!done) {
        TrackRoot *anchor;
        TrackSolution accepted = {0};
        cd xPred, reference;
        double step, judgeScale, searchScale, localRate, observed;
        int axis, aip = 1, aiz = 1, informative = 0, hasRef;
        int cell, attempts = 0, rejected = 0, flagged = 0, rescued = 0;
        int engineStatus = 0, foundClean = 0, seedCount = 0, seedIndex;
        int newLine = 0;
        cd seeds[9];

        ip = tr_grid_index(plg, pOrigin, PM[3], kperpSize);
        iz = tr_grid_index(zlg, zOrigin, ZM[3], kparSize);
        cell = tr_cell(ip, iz);
        anchor = tr_predecessor(roots, ip, iz, zFirst, &axis, &aip, &aiz);
        if (!anchor) {
            xPred = xseed;
            reference = xseed;
            hasRef = 0;
        } else {
            xPred = tr_predict(anchor, p, z, &informative);
            reference = informative ? xPred : anchor->x;
            hasRef = 1;
            step = fabs(p - anchor->p) + fabs(z - anchor->z);
        }
        if (!anchor)
            step = 0.0;

        localRate = tr_nearest_rate(rates, axis, aip, aiz);
        if (isfinite(localRate) && step > 0.0) {
            judgeScale = fmax(localRate * step,
                sqrt(DBL_EPSILON) * fmax(1.0, cabs(reference)));
        } else {
            judgeScale = sqrt(DBL_EPSILON) * fmax(1.0, cabs(reference));
        }
        if (!hasRef) {
            searchScale = INFINITY;
        } else {
            observed = isfinite(localRate) ? localRate * step : 0.0;
            searchScale = fmax(fmax(cabs(reference - anchor->x), observed),
                fmax(judgeScale, sqrt(DBL_EPSILON) * fmax(1.0, cabs(anchor->x))));
        }
        predictedOUT[cell] = reference;

        seeds[seedCount++] = xseed;
        if (anchor) {
            if (xPred != xseed)
                seeds[seedCount++] = xPred;
            if (anchor->x != xseed && anchor->x != xPred)
                seeds[seedCount++] = anchor->x;
            if (trackOptions.useScan && isfinite(searchScale)) {
                for (i = 0; i < trackOptions.numScanOffsets; ++i)
                    seeds[seedCount++] =
                        xPred + trackOptions.scanOffsets[i] * searchScale;
            }
        }

        ++counts.pointsVisited;
        for (seedIndex = 0; seedIndex < seedCount; ++seedIndex) {
            TrackSolution candidate;
            double radius = seedIndex == 0 ? INFINITY :
                trackOptions.trustFactor * searchScale;
            cd center = hasRef ? reference : seeds[seedIndex];
            if (!isfinite(creal(seeds[seedIndex])) ||
                !isfinite(cimag(seeds[seedIndex])) || creal(seeds[seedIndex]) <= 0.0)
                continue;
            ++attempts;
            candidate = tr_solve(p, z, seeds[seedIndex], radius, center);
            if (seedIndex == 0)
                engineStatus = candidate.status;
            counts.clamped += candidate.clamped;
            if (candidate.status != 1)
                continue;
            if (!isfinite(creal(reference)) || !isfinite(cimag(reference)) ||
                !isfinite(localRate) || step <= 0.0 ||
                cabs(candidate.x - reference) <= trackOptions.jumpFactor * judgeScale) {
                accepted = candidate;
                foundClean = 1;
                flagged = 0;
                rescued = seedIndex > 0;
                break;
            }
            ++rejected;
            if (accepted.status != 1) {
                accepted = candidate;
                flagged = 1;
            }
        }
        counts.attempts += attempts;
        counts.rejected += rejected;

        /* "Too heavily damped" is a physical statement about the plasma, not
         * a numerical one about the iteration, so it outranks any root a
         * fallback seed happened to find.  The reference returns immediately
         * on a CLEAN acceptance -- which therefore survives the override --
         * but a root kept only as a last resort against the continuity check
         * falls through and is reported as damped instead. */
        if (engineStatus == 2 && !foundClean)
            accepted.status = 2;

        if (diagnostics && diagnostics->attemptsPerPoint &&
            counts.pointsVisited <= maxRateCells)
            diagnostics->attemptsPerPoint[counts.pointsVisited - 1] =
                (double)attempts;

        if (accepted.status == 1) {
            double winding = NAN;
            tr_store_solution(&accepted, p, z);
            if (accepted.status == 1) {
                tr_save_root(ip, iz, &accepted);
                tr_record_solution(cell, roots, &accepted, p, z);
                ++counts.converged;
                counts.flagged += flagged && !foundClean;
                counts.rescued += rescued;
                if (trackOptions.certify) {
                    double radius = trackOptions.certRadius > 0.0
                        ? trackOptions.certRadius
                        : 1.0e-4 * fmax(1.0, cabs(accepted.x));
                    winding = tr_winding(p, z, accepted.x, radius,
                                         trackOptions.certSamples);
                    certWindingOUT[cell] = winding;
                    if (isnan(winding)) {
                        flagCertUnknownOUT[cell] = 1;
                        ++counts.certUnknown;
                    } else if (winding != 1.0) {
                        flagNotCertifiedOUT[cell] = 1;
                        ++counts.notCertified;
                    } else {
                        ++counts.certified;
                    }
                }
                if (!flagged && anchor && axis >= 1 && step > 0.0)
                    rates[2 * cell + axis - 1] =
                        cabs(accepted.x - anchor->x) / step;
                if (flagged && !foundClean)
                    flagSuspectBranchOUT[cell] = 1;

                last.valid = 1;
                last.p = p; last.z = z;
                last.x = accepted.x; last.dp = accepted.dp;
                last.dz = accepted.dz; last.dx = accepted.dx;
                if (lineFirst) {
                    lineAnchor = last;
                    lineFirst = 0;
                }
            }
        }
        if (accepted.status != 1 && (accepted.status == 2 || engineStatus == 2)) {
            F2D(flagTooHeavilyDampedOUT, ip, iz) = 1;
        } else if (accepted.status != 1) {
            F2D(flagNoConvergenceOUT, ip, iz) = 1;
        }

        if (accepted.status == 1) {
            xseed = accepted.x;
        } else if (last.valid) {
            int seedInfo;
            xseed = tr_predict(&last, p, z, &seedInfo);
        }

        if (accepted.status != 1) {
            if (zFirst)
                zlg = SWEEP_END;
            else
                plg = SWEEP_END;
        }
        if (!zFirst) {
            plg += PM[3];
            if (plg >= PM[1] && plg <= PM[2]) {
                p = tr_physical(plg);
            } else {
                zlg += ZM[3];
                if (zlg < ZM[1] || zlg > ZM[2]) {
                    done = 1;
                    continue;
                }
                lineFirst = 1;
                newLine = 1;
                plg = pOrigin;
                p = tr_physical(plg);
                z = tr_physical(zlg);
                if (lineAnchor.valid) {
                    int seedInfo;
                    xseed = tr_predict(&lineAnchor, p, z, &seedInfo);
                } else {
                    xseed = XOI;
                }
            }
        } else {
            zlg += ZM[3];
            if (zlg >= ZM[1] && zlg <= ZM[2]) {
                z = tr_physical(zlg);
            } else {
                plg += PM[3];
                if (plg < PM[1] || plg > PM[2]) {
                    done = 1;
                    continue;
                }
                lineFirst = 1;
                newLine = 1;
                zlg = zOrigin;
                z = tr_physical(zlg);
                p = tr_physical(plg);
                if (lineAnchor.valid) {
                    int seedInfo;
                    xseed = tr_predict(&lineAnchor, p, z, &seedInfo);
                } else {
                    xseed = XOI;
                }
            }
        }
        if (!done && accepted.status == 1 && !newLine) {
            int seedInfo;
            xseed = tr_predict(&last, p, z, &seedInfo);
        }
    }

    if (diagnostics) {
        diagnostics->pointsVisited = counts.pointsVisited;
        diagnostics->converged = counts.converged;
        diagnostics->rejected = counts.rejected;
        diagnostics->clampedSteps = counts.clamped;
        diagnostics->rescued = counts.rescued;
        diagnostics->flagged = counts.flagged;
        diagnostics->certified = counts.certified;
        diagnostics->notCertified = counts.notCertified;
        diagnostics->certUnknown = counts.certUnknown;
        diagnostics->meanAttempts = counts.pointsVisited
            ? counts.attempts / counts.pointsVisited : NAN;
    }
    free(roots);
    free(rates);
}
