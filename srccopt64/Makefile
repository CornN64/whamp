#
# Makefile -- build the WHAMP C translation.
#
# Targets:
#   all                 whamp, libwhamp.a, whamp_engine_test
#                       and the whamp_track_test regression runner
#   whamp               interactive WHAMP program
#   libwhamp.a          engine library
#   whamp_engine_test   host program that drives the engine
#   whamp_track_test    non-interactive tracker regression test
#   clean
#

CC      = gcc
AR      = ar r
BINDIR  = .

# MATLAB's MEX compiler (must be on PATH; override with MEX=/path/to/mex).
MEX     = mex

CFLAGS  = -c -std=c11 -Wall -Wextra -O2

# Numerical kernel shared by both executables.
KERNEL_OBJS = energy.o rint.o rtay.o ryla.o rasy.o xsi.o difu.o av.o

WHAMP_OBJS       = $(KERNEL_OBJS) fedit.o output.o typin.o input.o \
                   comin.o comcout.o whamp.o
ENGINE_OBJS      = $(KERNEL_OBJS) fedit.o output.o typin.o input.o \
                   comin.o comcout.o comoutput.o whamp_engine.o \
                   comtrack.o whamp_track.o whamp_engine_test.o
LIB_OBJS         = whamp_engine.o whamp_track.o comin.o comoutput.o \
                   comtrack.o comcout.o difu.o \
                   xsi.o energy.o ryla.o rtay.o rasy.o rint.o
TRACK_TEST_OBJS  = $(KERNEL_OBJS) fedit.o output.o typin.o input.o \
                   comin.o comcout.o comoutput.o comtrack.o \
                   whamp_engine.o whamp_track.o whamp_track_test.o
TRACK_DUMP_OBJS  = $(KERNEL_OBJS) fedit.o output.o typin.o input.o \
                   comin.o comcout.o comoutput.o comtrack.o \
                   whamp_engine.o whamp_track.o whamp_track_dump.o

all: whamp libwhamp.a whamp_engine_test whamp_track_test

# Every source includes whamp.h, which pulls in constants.h and fedit.h, so
# a change to any of the three rebuilds the whole tree.
%.o: %.c whamp.h constants.h fedit.h
	$(CC) $(CFLAGS) $<

comtrack.o: comtrack.h comoutput.h
whamp_track.o: whamp_track.h comtrack.h comoutput.h comin.h comcout.h
whamp_track_test.o: whamp_track.h comtrack.h comoutput.h comin.h

whamp: $(WHAMP_OBJS)
	$(CC) $(WHAMP_OBJS) -lm -o $(BINDIR)/$@

whamp_engine_test: $(ENGINE_OBJS)
	$(CC) $(ENGINE_OBJS) -lm -o $(BINDIR)/$@

whamp_track_test: $(TRACK_TEST_OBJS)
	$(CC) $(TRACK_TEST_OBJS) -lm -o $(BINDIR)/$@

whamp_track_dump: $(TRACK_DUMP_OBJS)
	$(CC) $(TRACK_DUMP_OBJS) -lm -o $(BINDIR)/$@

libwhamp.a: $(LIB_OBJS)
	$(AR) $@ $(LIB_OBJS)

# MATLAB gateway (translation of ../src/mexwhamp.F).  Requires MATLAB's mex
# on PATH and produces mexwhamp.mexa64 / .mexw64 / .mexmaci64 depending on
# the host.  This is the all-double engine variant; ../srccopt is the variant
# that reproduces the Fortran reference bit-for-bit.
mex: libwhamp.a mexwhamp.c
	$(MEX) -O -R2017b mexwhamp.c libwhamp.a

# Struct-interface wrapper: mexwrap(PlasmaModel, InputParameters) exactly as
# whamp.m takes them, returning whamp.m's Output structure.  Links the engine
# in, so the MEX is self-contained; ../matlab/mexwrap.m only locates it.
mexwrap: libwhamp.a mexwrap.c
	$(MEX) -O -R2017b mexwrap.c libwhamp.a

clean:
	rm -f *.o libwhamp.a whamp whamp_engine_test whamp_track_test \
		whamp_engine_test.exe whamp_track_test.exe whamp.exe \
		mexwhamp.mexa64 mexwhamp.mexw64 mexwhamp.mexmaci64 \
		mexwrap.mexa64 mexwrap.mexw64 mexwrap.mexmaci64

.PHONY: all clean mex mexwrap
