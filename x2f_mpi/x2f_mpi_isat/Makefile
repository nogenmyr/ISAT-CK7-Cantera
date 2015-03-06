# Author: Varun Hiremath <vh63@cornell.edu>
# Last Modified:  Thu,  8 Jul 2010 17:44:24 -0400

include ../build-files/vars.mk

LIB_NAME = libx2f_mpi_isat

# ------------- OBJECTS --------------------

include SOURCES.mk

OBJS=$(FILES:.f90=.$(OBJ))

# ------------- TARGETS --------------------

all: $(FILES) $(CREATE_LIBS)

include ../build-files/build-targets.mk

# ----------- DEPENDENCIES -----------------

include DEPENDS.mk
