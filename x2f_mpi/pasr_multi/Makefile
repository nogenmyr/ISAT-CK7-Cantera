# Author: Varun Hiremath <vh63@cornell.edu>
# Date: Thu, 29 Nov 2012 20:06:21 -0600

include ../build-files/vars.mk

X2F_LIBS_PATH = ../lib/
ISAT_LIBS_PATH = ../../ISAT/lib/

LDFLAGS= -L$(X2F_LIBS_PATH) -L$(ISAT_LIBS_PATH) -Xlinker --start-group -lisat7_mpi -lx2f_mpi -lx2f_mpi_isat $(MKL_LIBS)  -Xlinker --end-group

NAME = pasr_multi

# ------------- OBJECTS --------------------

include SOURCES.mk

OBJS = $(FILESf90:.f90=.$(OBJ)) $(FILESf:.f=.$(OBJ))

# ------------- TARGETS --------------------

all: $(FILES) $(NAME)

include ../build-files/build-targets.mk

# ----------- DEPENDENCIES -----------------

include DEPENDS.mk
