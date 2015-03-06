# Author: Varun Hiremath <vh63@cornell.edu>
# Date: Tue, 29 Nov 2011 20:21:26 -0600

# intel compiler is used for building, ensure that intel module is loaded:
# -> Lonestar: module load mkl (intel module is loaded by default)
# -> Kraken:   module load PrgEnv-intel
# Ranger use 'mpif90' wrapper; and on Kraken use 'ftn'

CC = icc
FC = ifort
ifeq ($(shell which mpif90 &> /dev/null; echo $$?), 0)
        FC = mpif90
	LINK_SYS_LIBS = -Bdynamic
endif
ifeq ($(shell which ftn &> /dev/null; echo $$?), 0)
        FC = ftn
	LINK_SYS_LIBS = -Bstatic
endif

LINK_LIBS = -Bstatic
# set mkl library path
#ifndef TACC_MKL_LIB
#        TACC_MKL_LIB = $(INTEL_PATH)/$(INTEL_MAJOR_VERSION)/$(INTEL_MINOR_VERSION)/mkl/lib/em64t/
#endif

# some compilations require ISAT modules
ISAT_MODS_PATH = ../../ISAT/isatab_mpi

# set mkl libraries for linking
MKL_LIBS = -mkl
#MKL_LIBS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lguide

# set debug/optimize modes
ifeq ($(BUILD_TYPE), debug)
	FCFLAGS = -fp-stack-check -traceback -g -debug all -debug-parameters all -fPIC -vec-report0 -DDEBUG -D_DEBUG -I$(ISAT_MODS_PATH)
else
	FCFLAGS = -O2 -fPIC -vec-report0 -I$(ISAT_MODS_PATH)
endif
CCFLAGS = -O2 -fPIC -vec-report0 -w

LN = ln -s
CP = cp
RM = rm -f
RMDIR = rm -rf
AR = ar crs
OBJ = o
MAKE=make BUILD_TYPE=$(BUILD_TYPE)
CREATE_LIBS = $(LIB_NAME).a $(LIB_NAME).so
