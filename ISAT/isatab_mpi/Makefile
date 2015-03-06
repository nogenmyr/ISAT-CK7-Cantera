# Author: Varun Hiremath <vh63@cornell.edu>
# Last Modified:  Thu, 12 Nov 2009 16:42:03 -0500

include ../build-files/vars.mk

SELL_LIB_PATH=../sell_lib
IFLAGS = -I$(SELL_LIB_PATH)
LIB_NAME = libisatab_mpi

# ------------- OBJECTS --------------------

include SOURCES.mk

OBJS=$(FILES:.f90=.$(OBJ))

# ------------- TARGETS --------------------

all: $(FILES) BUILD-SELL $(CREATE_LIBS)

BUILD-SELL:
	(cd $(SELL_LIB_PATH); $(MAKE))
	touch $@

include ../build-files/build-targets.mk

clean::
	$(RM) BUILD-SELL

# ----------- DEPENDENCIES -----------------

include DEPENDS.mk
