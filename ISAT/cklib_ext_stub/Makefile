# Author: Varun Hiremath <vh63@cornell.edu>
# Last Modified: Tue Jan 20 09:37:26 CET 2015

include ../build-files/vars.mk

ISATAB_PATH=../../ISAT/isatab_ser
ISATCK_PATH=../../ISAT/isat-ck

IFLAGS = -I$(ISATCK_PATH) -I$(ISATAB_PATH) -I$(CANTERA_INCPATH) -I$(CANTERA_SRC)

LIB_NAME = libck_ext

LDFLAGS = -L$(CANTERA_LIBPATH) -lcantera_fortran -lcantera -lpthread -lstdc++ 

# ------------- OBJECTS --------------------

include SOURCES.mk

OBJS=$(FILESf90:.f90=.$(OBJ)) $(FILESf:.f=.$(OBJ))

# ------------- TARGETS --------------------

all: $(FILES) BUILD-ISATAB BUILD-ISATCK $(CREATE_LIBS)

BUILD-ISATAB:
	(cd $(ISATAB_PATH); $(MAKE))
	touch $@

BUILD-ISATCK:
	(cd $(ISATCK_PATH); $(MAKE))
	touch $@

include ../build-files/build-targets.mk

clean::
	$(RM) BUILD-ISATAB BUILD-ISATCK
