/*
 * input.c -- translation of SUBROUTINE read_input_file (src/input.f90)
 *
 * The model file holds seven list-directed records of ten real values each
 * (density, temperature, delta, alpha1, alpha2, mass, drift) followed by two
 * scalar records (gyrofrequency and the log-scale flag).  As in the Fortran
 * original, reading stops at the first record that cannot be parsed.
 */
#include "whamp.h"
#include "comin.h"

#define NFIELD 10 /* values per vector record in the model file */

/* Read `count` list-directed reals into the array starting at `first`,
 * advancing `stride` bytes between elements (so that the columns of AA can
 * be read directly).  Returns 0 on success, nonzero on a read error or
 * premature end of file. */
static int read_reals(FILE *file, double *first, size_t stride, int count)
{
    for (int i = 0; i < count; i++) {
        double *item = (double *)((char *)first + (size_t)i * stride);
        if (fscanf(file, "%lf", item) != 1)
            return 1;
    }
    return 0;
}

/* Convenience wrapper for a contiguous record of NFIELD values stored in the
 * Fortran element range values[1..NFIELD]. */
static int read_record(FILE *file, double *values)
{
    return read_reals(file, values + 1, sizeof(double), NFIELD);
}

void read_input_file(const char *FILENAME)
{
    FILE *file;
    int readOk;

    if (FILENAME[0] == '\0') {
        printf(" ERROR model filename not given!\n");
        exit(0);
    }

    if (printDebugInfo)
        printf(" # read_input_file: file =  %s\n", FILENAME);

    file = fopen(FILENAME, "r");
    if (file == NULL) {
        printf(" READ_INPUT_FILE: ERROR IN OPEN, FILE =  %s\n", FILENAME);
        return;
    }

    /* The vector records are read in the order the Fortran namelist uses.
     * The two alpha columns live in one 2D array, so they are read with the
     * row stride of AA.  The && chain stops at the first failing record,
     * matching the Fortran control flow. */
    readOk = read_record(file, DN) == 0 && read_record(file, TA) == 0 &&
             read_record(file, DD) == 0 &&
             read_reals(file, &AA[1][1], sizeof(AA[0]), NFIELD) == 0 &&
             read_reals(file, &AA[1][2], sizeof(AA[0]), NFIELD) == 0 &&
             read_record(file, ASS) == 0 && read_record(file, VD) == 0 &&
             fscanf(file, "%lf", &XC) == 1 && fscanf(file, "%lf", &PZL) == 1;

    fclose(file);

    if (!readOk)
        printf(" READ_INPUT_FILE: ERROR IN READ, FILE =  %s\n", FILENAME);
}
