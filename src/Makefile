##########################################
# Makefile for eccsvm                    #
# Copyright (C) 2004 Louis F. Rossi      #
# You may distribute this under the      #
# terms of the GNU Public License        #
# included in this package.              #
##########################################
#
# example: make eflow BINNAME=eflow


HOSTTYPE 	:= $(shell arch)

ifeq ($(mpi),on)
CC	= mpicc
LN	= mpicc
MPICFLAGS  = -DMULTIPROC
else
CC		= gcc
LN		= gcc
endif

# Set your optional C flags here.  It depends upon your compiler and 
# your architecture.  Caution, mpicc is built on top of a C compiler and
# you'll need to know which one.

# These are good with pgcc.
#COPTFLAGS	= -fast  
# These are good with gcc.
#COPTFLAGS	= -m64 -O3

# This one is an example of how to change the location of Rodrigo's
# data files at compile time.
#COPTFLAGS	= -fast  -D'DATA_ROOT="/home/rossi/ugh/"'
COPTFLAGS	= -O3
CWARNFLAGS	= -Wall

ifeq ($(ccsvm),on)
CCSVM   = -DCCSVM
endif

# Default libs.
CLIBS		= -lm -llapack -lblas

ifeq ($(HOSTTYPE),i686)
CLIBS 		= -lm -llapack -lblas
endif

ifeq ($(HOSTTYPE),i386-linux)
CLIBS 		= -lm -llapack -lblas
endif

ifeq ($(HOSTTYPE),x86_64)
# This is good with pgcc
#CLIBS 		=  -lm  -llapack -lblas -lpgftnrtl
# This works with gcc.
CLIBS 		=  -lm  -llapack -lblas
endif

ifeq ($(HOSTTYPE),macintosh)
CLIBS 		= -lm  -llapack  -lblas
endif

ifeq ($(HOSTTYPE),sun4)
COPTFLAGS	= -fast -xO5 -xcrossfile -xdepend -fsimple=2
CLIBS 		= -fast -xO5 -xcrossfile -xdepend  -fsimple=2 -lm
endif

ifeq ($(xantisymm),on) 
XANTISYMMCFLAGS  = -DXANTISYMM
endif

CFLAGS = $(CCSVM) $(COPTFLAGS) $(MPICFLAGS) $(XANTISYMMCFLAGS) $(CWARNFLAGS)

OBJECTS = deriv.o solve.o evolve.o files.o flip.o init.o \
	removal.o veldev.o deform.o \
	velocity.o partition.o mpsum.o mpcoeffs.o \
	reflect.o master-slave.o peer.o cache_resort.o \
	stream_fn.o potential.o boundary.o\
	polynB.o RHE_interp.o \
	eval_biot.o eval_doms.o find_coeffs.o read_data.o 

.PHONY: dummy

dummy:
	$(error You must specify a target)

eflow: $(OBJECTS)
	$(LN)  $(CPROF_FLAGS) -o $(BINNAME) $(OBJECTS) $(CLIBS)

$(OBJECTS): multiproc.h global_matrices.h global_min.h particle.h \
	biot-savart.h boundary.h field_interp.h linear_velocity.h potential.h

.c.o:
	$(CC) -c $(CFLAGS) $(CPROF_FLAGS) $<

clean:
	rm -f *.o a.out core 
