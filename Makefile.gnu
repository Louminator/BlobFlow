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

COPTFLAGS	=  -Wall
CLIBS 		= -lm 

ifeq ($(HOSTTYPE),sun4)
COPTFLAGS	= -fast -xO5 -xcrossfile -xdepend -fsimple=2
CLIBS 		= -fast -xO5 -xcrossfile -xdepend  -fsimple=2 -lm
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

OBJECTS = deriv.o evolve.o files.o flip.o init.o \
	merge.o removal.o solve.o split.o veldev.o \
	split14.o split15.o split15asym.o \
	velocity.o partition.o mpsum.o mpcoeffs.o \
	reflect.o master-slave.o peer.o cache_resort.o \
	stream_fn.o \
	polynB.o


eflow: $(OBJECTS)
	$(LN) $(CPROF_FLAGS) -o $(TARGET) $(OBJECTS) $(CLIBS) $(LAMLIB)

$(OBJECTS): global.h multiproc.h 

.c.o:
	$(CC) -c $(CFLAGS) $(CPROF_FLAGS) $<

clean:
	rm -f *.o a.out core 
