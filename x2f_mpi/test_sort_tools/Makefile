# Author: Varun Hiremath <vh63@cornell.edu>
# Date: Thu, 29 Nov 2012 20:06:10 -0600

include ../build-files/vars.mk

IFLAGS = -I../x2f_mpi/
X2F_LIBS_PATH = ../lib/
ISAT_LIBS_PATH = ../../ISAT/lib/
LDFLAGS= -L$(X2F_LIBS_PATH) -L$(ISAT_LIBS_PATH) -L$(TACC_MKL_LIB) -Xlinker --start-group -lx2f_mpi $(MKL_LIBS) -Xlinker --end-group

NAME = test_sort_tools

# ------------- OBJECTS --------------------

include SOURCES.mk

OBJS = $(FILES:.f90=.$(OBJ))

# ------------- TARGETS --------------------

all: $(FILES) $(NAME)

include ../build-files/build-targets.mk

# ----------- DEPENDENCIES -----------------

include DEPENDS.mk
