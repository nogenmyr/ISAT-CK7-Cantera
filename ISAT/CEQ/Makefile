# Author: Varun Hiremath <vh63@cornell.edu>
# Last Modified:  Thu, 12 Nov 2009 15:49:41 -0500

include ../build-files/vars.mk

LIB_NAME = libceq

# ------------- OBJECTS --------------------

include SOURCES.mk

OBJS=$(FILES:.f90=.$(OBJ))

# ------------- TARGETS --------------------

all: $(FILES) $(CREATE_LIBS)

include ../build-files/build-targets.mk

# ----------- DEPENDENCIES -----------------

include DEPENDS.mk
