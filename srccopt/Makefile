#
# Makefile -- build the WHAMP C translation.
#
# Targets:
#   all                 whamp, libwhamp.a, whamp_engine_test
#   whamp               interactive WHAMP program
#   libwhamp.a          engine library
#   whamp_engine_test   host program that drives the engine
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
                   whamp_engine_test.o
LIB_OBJS         = whamp_engine.o comin.o comoutput.o comcout.o difu.o \
                   xsi.o energy.o ryla.o rtay.o rasy.o rint.o

all: whamp libwhamp.a whamp_engine_test

%.o: %.c whamp.h
	$(CC) $(CFLAGS) $<

whamp: $(WHAMP_OBJS)
	$(CC) $(WHAMP_OBJS) -lm -o $(BINDIR)/$@

whamp_engine_test: $(ENGINE_OBJS)
	$(CC) $(ENGINE_OBJS) -lm -o $(BINDIR)/$@

libwhamp.a: $(LIB_OBJS)
	$(AR) $@ $(LIB_OBJS)

# MATLAB gateway (translation of ../src/mexwhamp.F).  Requires MATLAB's mex
# on PATH and produces mexwhamp.mexa64 / .mexw64 / .mexmaci64 depending on
# the host.  Move the result into ../matlab/+whamp/ where run.m calls it.
mex: libwhamp.a mexwhamp.c
	$(MEX) -O -R2017b mexwhamp.c libwhamp.a

clean:
	rm -f *.o libwhamp.a whamp whamp whamp_engine_test whamp_engine_test.exe whamp.exe \
		mexwhamp.mexa64 mexwhamp.mexw64 mexwhamp.mexmaci64

.PHONY: all clean mex
