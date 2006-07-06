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
endif

ifeq ($(HOSTTYPE),i386-linux)
COPTFLAGS	= -O3 -Wall
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
	stream_fn.o potential.o boundary.o \
	polynB.o \
	biot.o biot_dom1.o biot_dom2.o biot_dom3.o Cheb.o \
	data_1a1p5_dom1.o  data_2a3_dom1.o  data_4a5_dom1.o  data_6a8_dom1.o \
	data_1a1p5_dom2.o  data_2a3_dom2.o  data_4a5_dom2.o  data_6a8_dom2.o \
	data_1a1p5_dom3.o  data_2a3_dom3.o  data_4a5_dom3.o  data_6a8_dom3.o \
	data_1p5a2_dom1.o  data_3a4_dom1.o  data_5a6_dom1.o  data_8a10_dom1.o \
	data_1p5a2_dom2.o  data_3a4_dom2.o  data_5a6_dom2.o  data_8a10_dom2.o \
	data_1p5a2_dom3.o  data_3a4_dom3.o  data_5a6_dom3.o  data_8a10_dom3.o

eflow: $(OBJECTS)
	$(LN) $(CPROF_FLAGS) -o $(TARGET) $(OBJECTS) $(CLIBS) $(LAMLIB)

$(OBJECTS): global.h multiproc.h biot_global.h

.c.o:
	$(CC) -c $(CFLAGS) $(CPROF_FLAGS) $<

clean:
	rm -f *.o a.out core 
