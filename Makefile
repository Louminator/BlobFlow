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

.PHONY: $(TARGET) clean

#==============================================================================#

$(TARGET):
	cd src && $(MAKE) TARGET=../$(TARGET) ../$(TARGET)

clean:
#	rm -f eflow.*
	cd src && $(MAKE) clean

#==============================================================================#
