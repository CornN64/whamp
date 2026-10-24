/* Small non-interactive smoke test for the Version 2 C tracker. */
#include "whamp_track.h"
#include "comin.h"
#include "comoutput.h"
#include "comtrack.h"

#include <assert.h>

int main(void)
{
    TrackDiagnostics diagnostics;
    cd engineRoot;
    cd zFirstRoots[6];
    int j, found = 0, certified = 0;

    DN[1] = DN[2] = 1.0e6;
    TA[1] = 0.01;
    TA[2] = 0.001;
    DD[1] = 1.0;
    AA[1][1] = 5.0;
    AA[1][2] = 0.1;
    AA[2][1] = 1.0;
    AA[2][2] = 0.1;
    ASS[1] = 16.0;
    VD[1] = 1.0;
    XC = 2.79928;
    PM[1] = 0.005;
    PM[2] = 0.006;
    PM[3] = 0.001;
    ZM[1] = 0.008;
    ZM[2] = 0.012;
    ZM[3] = 0.002;
    PZL = 0.0;
    cycleZFirst = 1;
    XOI = 0.1;
    maxIterations = 50;

    ZM[2] = ZM[1];
    whamp_engine();
    assert(flagSolutionFoundOUT[0]);
    engineRoot = fOUT[0];

    ZM[2] = 0.012;
    track_options_default();
    whamp_track_v2(&diagnostics);

    assert(diagnostics.pointsVisited == 6);
    assert(kperpSize == 2 && kparSize == 3);
    for (j = 0; j < kperpSize * kparSize; ++j) {
        if (flagSolutionFoundOUT[j])
            ++found;
        if (certWindingOUT[j] == 1.0)
            ++certified;
    }
    assert(found > 0);
    assert(certified > 0);
    assert(found == 6);
    for (j = 0; j < 6; ++j)
        zFirstRoots[j] = fOUT[j];
    assert(cabs(fOUT[0] - engineRoot) <=
           1.0e-12 * fmax(1.0, cabs(engineRoot)));
    assert(diagnostics.certified == certified);

    cycleZFirst = 0;
    whamp_track_v2(&diagnostics);
    assert(diagnostics.pointsVisited == 6);
    assert(kperpSize == 2 && kparSize == 3);
    for (j = 0; j < 6; ++j) {
        assert(flagSolutionFoundOUT[j]);
        assert(cabs(fOUT[j] - zFirstRoots[j]) <=
               1.0e-8 * fmax(1.0, cabs(zFirstRoots[j])));
    }
    assert(cabs(fOUT[0] - engineRoot) <=
           1.0e-12 * fmax(1.0, cabs(engineRoot)));
    return 0;
}
