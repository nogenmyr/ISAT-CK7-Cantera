# Author: Varun Hiremath <vh63@cornell.edu>
# Last Modified:  Thu,  8 Jul 2010 17:43:03 -0400

include ../build-files/vars.mk

LIB_NAME = libx2f_mpi

# ------------- OBJECTS --------------------

include SOURCES.mk

OBJS=$(FILESf90:.f90=.$(OBJ)) $(FILESc:.c=.$(OBJ)) 

# ------------- TARGETS --------------------

all: $(FILES) $(CREATE_LIBS)

include ../build-files/build-targets.mk

# ----------- DEPENDENCIES -----------------

include DEPENDS.mk
