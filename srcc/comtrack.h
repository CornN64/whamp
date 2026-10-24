/* comtrack.h -- output arrays that whamp_track_v2() adds over whamp_engine().
 *
 * whamp_engine() writes the shared comoutput arrays (fOUT, ExOUT, ...,
 * flagNoConvergenceOUT).  Branch tracking needs five more, because it reports
 * not only whether a point converged but whether the converged value is the
 * branch and whether it is a root at all:
 *
 *   flagSuspectBranchOUT  converged, but further from the branch prediction
 *                         than the measured local rate allows.  Kept, because
 *                         refusing a converged root is worse than keeping a
 *                         doubtful one, but the reader must see it.
 *   certWindingOUT        (Z - P) for D around a small contour centred on the
 *                         reported root.  +1 is a simple zero and no pole.
 *                         NaN where certification was refused or did not run.
 *   flagNotCertifiedOUT   converged root the argument principle does NOT
 *                         report as a single simple zero -- "not a root", or
 *                         "sitting on a gyroharmonic pole".
 *   flagCertUnknownOUT    the contour could not be evaluated (the pole
 *                         expansion was refused somewhere on it).  Neither
 *                         trusted nor rejected.
 *   predictedOUT          the implicit-differentiation reference each point was
 *                         judged against, so every verdict can be audited
 *                         after the fact.
 *
 * Same 1-based Fortran-order storage and F2D accessor as comoutput.h; the
 * arrays are sized by the same kperpSize / kparSize.
 */
#ifndef COMTRACK_H
#define COMTRACK_H

#include "whamp.h"

extern cd *predictedOUT;
extern double *certWindingOUT;

extern int *flagSuspectBranchOUT;
extern int *flagNotCertifiedOUT;
extern int *flagCertUnknownOUT;

/* Allocate the five tracking arrays to the current kperpSize x kparSize.
 * certWindingOUT is filled with NaN, since zero would read as a verdict. */
void track_allocate_output(void);

#endif /* COMTRACK_H */
