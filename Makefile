#==============================================================================#
# Makefile for eccsvm
#==============================================================================#

TARGET := eflow

ifeq ($(mpi),on)
TARGET := $(addsuffix .mpi,$(TARGET))
endif

ifeq ($(ccsvm),on)
TARGET := $(addsuffix .ccsvm,$(TARGET))
endif

ifeq ($(xantisymm),on) 
TARGET := $(subst eflow,eflow.xanti,$(TARGET))
endif

.PHONY: $(TARGET) clean

#==============================================================================#

$(TARGET):
	cd src && $(MAKE) TARGET=../$(TARGET) ../$(TARGET)

clean:
#	rm -f eflow.*
	cd src && $(MAKE) clean

#==============================================================================#
