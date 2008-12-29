##########################################
# Makefile for eccsvm                    #
# Copyright (C) 2000 Louis F. Rossi      #
# You may distribute this under the      #
# terms of the GNU Public License        #
# included in this package.              #
##########################################
TARGET		:= eflow
HOSTTYPE 	:= $(shell arch)

ifeq ($(mpi),on)
CC	= hcc
LN	= hcc
MPICFLAGS  = -DMULTIPROC
LAMLIB	= -lmpi
TARGET	:= $(addsuffix .mpi,$(TARGET))
else
CC	= gcc
LN	= gcc
LAMLIB	= 
endif

COPTFLAGS	= -O3 -Wall
CLIBS 		= -lm 

ifeq ($(HOSTTYPE),sun4)
COPTFLAGS	= -fast -xO5 -xcrossfile -xdepend -fsimple=2
CLIBS 		= -fast -xO5 -xcrossfile -xdepend  -fsimple=2 -lm
#COPTFLAGS 	= -O1 -Wall \
#		-fcaller-saves -fcse-follow-jumps -fcse-skip-blocks \
#              -felide-constructors \
#              -fexpensive-optimizations -ffast-math -ffloat-store \
#              -fforce-addr -fforce-mem -finline-functions \
#              -fkeep-inline-functions -fmemoize-lookups \
#              -fomit-frame-pointer -frerun-cse-after-loop \
#              -fschedule-insns -fschedule-insns2 \
#              -fstrength-reduce -fthread-jumps \
#              -funroll-loops
endif

ifeq ($(HOSTTYPE),i386-linux)
COPTFLAGS	= -O1 -Wall
CLIBS 		= -lm 
endif

ifeq ($(HOSTTYPE),macintosh)
COPTFLAGS	= -fast -Wall
CLIBS 		= -lm 
endif

ifeq ($(xantisymm),on) 
XANTISYMMCFLAGS  = -DXANTISYMM
TARGET	:= $(subst eflow,eflow.xanti,$(TARGET))
endif

#CFLAGS = -fast $(MPICFLAGS) $(XANTISYMMCFLAGS) $(TRACEFLAGS) $(SCHFLAG)
CFLAGS = $(COPTFLAGS) $(MPICFLAGS) $(XANTISYMMCFLAGS) $(TRACEFLAGS) $(SCHFLAG)

# CFLAGS	= -fast
# CPROF_FLAGS = -pg -g3

# CLIBS = -lm -lblas -llapack
# CLIBS = -lm -llapack -lblas -lfor

OBJECTS = deriv.o evolve.o files.o flip.o init.o \
	merge.o removal.o solve.o split.o veldev.o \
	split14.o split15.o split15asym.o \
	velocity.o partition.o mpsum.o mpcoeffs.o \
	reflect.o master-slave.o peer.o cache_resort.o \
	polynB.o potential.o
#	lamb-oseen.o

eflow: $(OBJECTS)
	$(LN) $(CPROF_FLAGS) -o $(TARGET) $(OBJECTS) $(CLIBS) $(LAMLIB)

$(OBJECTS): global.h multiproc.h 

.c.o:
	$(CC) -c $(CFLAGS) $(CPROF_FLAGS) $<

clean:
	rm -f *.o a.out core 
