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

# Every source includes whamp.h, which pulls in constants.h and fedit.h, so
# a change to any of the three rebuilds the whole tree.
%.o: %.c whamp.h constants.h fedit.h
	$(CC) $(CFLAGS) $<

whamp: $(WHAMP_OBJS)
	$(CC) $(WHAMP_OBJS) -lm -o $(BINDIR)/$@

whamp_engine_test: $(ENGINE_OBJS)
	$(CC) $(ENGINE_OBJS) -lm -o $(BINDIR)/$@

libwhamp.a: $(LIB_OBJS)
	$(AR) $@ $(LIB_OBJS)

clean:
	rm -f *.o libwhamp.a whamp whamp whamp_engine_test whamp_engine_test.exe whamp.exe

.PHONY: all clean
