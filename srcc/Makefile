#*
# Makefile     -- Makefile to build the WHAMP C translation.
#
# Mirrors src/Makefile: same three targets (whamp, libwhamp.a,
# whamp_engine_test) with the same object lists, plus the branch tracker
# and the two MATLAB gateways.
#
# Targets:
#   whamp               interactive WHAMP program
#   libwhamp.a          engine and tracker library
#   whamp_engine_test   host program that drives the engine
#   whamp_track_test    non-interactive tracker regression test
#   whamp_track_dump    tracker transcript for ../srcm/whamp_track_v2_vs_c.m
#   mex                 MATLAB gateway (needs mex on PATH)
#   mexwrap             MATLAB struct-interface wrapper (needs mex on PATH)
#
# This is the bit-for-bit variant: it emulates Fortran single-precision
# entities through F32(), so its printed output matches the gfortran
# reference.  The tracker uses the same kernel and so inherits that.
#
#===============================================================================

CC   = gcc
AR   = ar r
BINDIR = .

# MATLAB's MEX compiler (must be on PATH; override with MEX=/path/to/mex).
MEX = mex

CFLAGS = -c -Wall -O2

all: whamp libwhamp.a whamp_engine_test whamp_track_test

# General rule
%.o: %.c whamp.h constants.h
	$(CC) $(CFLAGS) $<


# Special dependencies
energy.o  output.o  difu.o: comcout.o
xsi.o     output.o  difu.o: comin.o
whamp_engine.o: comoutput.o
whamp_engine_test.o: comoutput.o
comtrack.o: comtrack.h comoutput.h
whamp_track.o: whamp_track.h comtrack.h comoutput.h comin.h comcout.h
whamp_track_test.o: whamp_track.h comtrack.h comoutput.h comin.h


OBJS =  energy.o rint.o\
	rtay.o xsi.o output.o ryla.o\
	rasy.o typin.o difu.o\
	whamp.o\
	input.o ctof.o av.o comin.o comcout.o

OBJSTEST = energy.o rint.o\
	rtay.o xsi.o output.o ryla.o\
	rasy.o typin.o difu.o\
	whamp_engine_test.o\
	whamp_engine.o comoutput.o comtrack.o whamp_track.o\
	input.o ctof.o av.o comin.o comcout.o

OBJSLIB = whamp_engine.o whamp_track.o comin.o comoutput.o comtrack.o\
	comcout.o difu.o xsi.o energy.o\
	ryla.o rtay.o rasy.o rint.o

TRACKTESTOBJS = energy.o rint.o\
	rtay.o xsi.o output.o ryla.o\
	rasy.o typin.o difu.o\
	whamp_track_test.o\
	whamp_engine.o comoutput.o comtrack.o whamp_track.o\
	input.o ctof.o av.o comin.o comcout.o

TRACKDUMPOBJS = energy.o rint.o\
	rtay.o xsi.o output.o ryla.o\
	rasy.o typin.o difu.o\
	whamp_track_dump.o\
	whamp_engine.o comoutput.o comtrack.o whamp_track.o\
	input.o ctof.o av.o comin.o comcout.o


whamp:  ${OBJS}
	${CC} ${OBJS} -lm -o ${BINDIR}/$@

whamp_engine_test: ${OBJSTEST}
	${CC} ${OBJSTEST} -lm -o ${BINDIR}/$@

whamp_track_test: ${TRACKTESTOBJS}
	${CC} ${TRACKTESTOBJS} -lm -o ${BINDIR}/$@

whamp_track_dump: ${TRACKDUMPOBJS}
	${CC} ${TRACKDUMPOBJS} -lm -o ${BINDIR}/$@

libwhamp.a: ${OBJSLIB}
	${AR} $@ ${OBJSLIB}


# Struct-interface wrapper: mexwrap(PlasmaModel, InputParameters) exactly as
# whamp.m takes them, returning whamp.m's Output structure.  Links the engine
# and tracker in, so the MEX is self-contained; ../matlab/mexwrap.m only
# locates it.  Requires MATLAB's mex on PATH.
mexwrap: libwhamp.a mexwrap.c
	${MEX} -O -R2017b mexwrap.c libwhamp.a


# Clean up everything
clean :
	rm -f *.o libwhamp.a whamp whamp_engine_test whamp_track_test\
	whamp_engine_test.exe whamp_track_test.exe whamp.exe\
	whamp_track_dump.exe\
	mexwrap.mexa64 mexwrap.mexw64 mexwrap.mexmaci64

.PHONY: all clean mexwrap
