# Author: Varun Hiremath <vh63@cornell.edu>
# Last Modified: Sat,  3 Apr 2010 17:56:29 -0400

include ../build-files/vars.mk

CEQ_PATH=../CEQ
ISATAB_PATH=../isatab_ser
CANTERALIBPATH=../canteralib
ISATCK_EXT_PATH=../isat-ck-ext

IFLAGS = -I$(CEQ_PATH) -I$(ISATAB_PATH) -I$(ISATCK_EXT_PATH) -I$(CANTERALIBPATH) -I$(CANTERA_INCPATH)

LDFLAGS = -L/usr/local/lib -lcantera_fortran -lcantera -lpthread -lstdc++

LIB_NAME = libisatck

# ------------- OBJECTS --------------------

include SOURCES.mk

OBJS=$(FILESf90:.f90=.$(OBJ)) $(FILESf:.f=.$(OBJ))

# ------------- TARGETS --------------------

all: $(FILES) BUILD-CEQ BUILD-ISATAB BUILD-STREAMS-MOD $(CREATE_LIBS)

BUILD-CEQ:
	(cd $(CEQ_PATH); $(MAKE))
	touch $@

BUILD-ISATAB:
	(cd $(ISATAB_PATH); $(MAKE))
	touch $@

BUILD-STREAMS-MOD:
	(cd $(ISATCK_EXT_PATH); $(MAKE) STREAMS-MOD)
	touch $@

include ../build-files/build-targets.mk

clean::
	$(RM) BUILD-CEQ BUILD-ISATAB BUILD-STREAMS-MOD

# ----------- DEPENDENCIES -----------------

include DEPENDS.mk
