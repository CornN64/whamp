/*
 * whamp_track.h -- branch-tracking sweep over the WHAMP dispersion relation.
 *
 * C translation of srcm/whamp_track_v2.m, the Version 2 branch tracker.  It
 * drives the same kernel the engine does (DIFU/CHI/ENERGY) and writes the same
 * comoutput arrays, plus the five tracking arrays in comtrack.h.
 *
 * Where whamp_engine() continues a sweep line with a first-order extrapolation
 * of the previous root and accepts whatever Newton lands on, whamp_track_v2()
 *
 *   1. predicts the next root by implicit differentiation of D(p, z, x(p,z)) = 0
 *          dx = -((p - p0) DP + (z - z0) DZ) / DX
 *      which is the same expression the engine already uses for its own seed,
 *      but evaluated against a map of accepted roots keyed by GRID CELL rather
 *      than against the single carried predictor -- so it survives a rejected
 *      point and it is correct across a two-dimensional sweep;
 *
 *   2. judges a converged root against the branch's own measured RATE
 *      |dx|/|dk| along the axis being stepped, read from the nearest grid cell
 *      that recorded one.  Nothing here is hand-tuned.  Version 1 used the
 *      median recent step floored by RELFLOOR = 0.25|x|, and that floor was
 *      wide enough to let a lucky seed land in a neighbouring basin;
 *
 *   3. certifies every accepted root with the argument principle, counting
 *      Z - P for D around a small contour.  A residual threshold cannot do
 *      this job: |D| spans some thirty orders of magnitude over the domain.
 *
 * The input interface is exactly whamp_engine()'s -- the comin module globals
 * (DN, TA, AA, DD, VD, ASS, PM, ZM, XC, PZL, cycleZFirst, XOI, maxIterations)
 * -- so a caller can swap one for the other.  Behavioural differences from
 * whamp_engine() are listed at whamp_track_v2().
 */
#ifndef WHAMP_TRACK_H
#define WHAMP_TRACK_H

#include "whamp.h"

/* Tunables.  The global is initialized to these defaults; callers can also
 * reset it explicitly with track_options_default().  Factors are relative to
 * measured quantities, not absolute frequency scales. */
typedef struct {
    double trustFactor;  /* Newton disc radius, in units of the search scale */
    double jumpFactor;   /* reject a root further than this many scales away */
    double scanOffsets[6];
    int numScanOffsets;  /* last-resort seeds: prediction + offset*scale */
    int useScan;
    int certify;         /* run the argument principle on accepted roots */
    double certRadius;   /* <= 0 selects 1e-4 * max(1, |x|) */
    int certSamples;     /* starting contour resolution; doubled if coarse */
} TrackOptions;

extern TrackOptions trackOptions;

void track_options_default(void);

/* Counts worth reporting after a sweep.  May be NULL. */
typedef struct {
    int pointsVisited;
    int converged;
    int rejected;    /* converged but off-branch, kept only as a last resort */
    int clampedSteps;
    int rescued;     /* needed a fallback seed */
    int flagged;     /* kept against the continuity check */
    int certified;
    int notCertified;
    int certUnknown;
    double meanAttempts;
} TrackDiagnostics;

/* Runs the sweep.  Reads the comin globals and fills the comoutput arrays
 * (fOUT, ExOUT, ..., flagNoConvergenceOUT, kperpOUT, kparOUT, kperpSize,
 * kparSize) plus the comtrack arrays.  Callers must not pre-allocate. */
void whamp_track_v2(TrackDiagnostics *diagnostics);

#endif /* WHAMP_TRACK_H */
