/* comtrack.c -- storage for the branch-tracking output arrays. */
#include "comtrack.h"
#include "comoutput.h"

cd *predictedOUT = NULL;
double *certWindingOUT = NULL;

int *flagSuspectBranchOUT = NULL;
int *flagNotCertifiedOUT = NULL;
int *flagCertUnknownOUT = NULL;

void track_allocate_output(void)
{
    size_t n = (size_t)kperpSize * (size_t)kparSize;
    size_t i;

    if (predictedOUT) free(predictedOUT);
    predictedOUT = (cd *)calloc(n, sizeof(cd));

    if (certWindingOUT) free(certWindingOUT);
    certWindingOUT = (double *)malloc(n * sizeof(double));

    if (flagSuspectBranchOUT) free(flagSuspectBranchOUT);
    flagSuspectBranchOUT = (int *)calloc(n, sizeof(int));

    if (flagNotCertifiedOUT) free(flagNotCertifiedOUT);
    flagNotCertifiedOUT = (int *)calloc(n, sizeof(int));

    if (flagCertUnknownOUT) free(flagCertUnknownOUT);
    flagCertUnknownOUT = (int *)calloc(n, sizeof(int));

    /* Zero would read as the verdict "no zero and no pole inside", which is
     * the opposite of "not evaluated".  Seed with NaN so an uncertified cell
     * is visibly uncertified. */
    for (i = 0; i < n; i++)
        certWindingOUT[i] = NAN;
}
