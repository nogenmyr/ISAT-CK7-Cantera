# Author: Varun Hiremath <vh63@cornell.edu>
# Last Modified: Sat,  3 Apr 2010 17:57:14 -0400

include ../build-files/vars.mk

ISATAB_PATH=../isatab_ser
ISATCK_PATH=../isat-ck
IFLAGS = -I$(ISATCK_PATH) -I$(ISATAB_PATH)

LIB_NAME = libisatck_ext

# ------------- OBJECTS --------------------

include SOURCES.mk

OBJS=$(FILESf90:.f90=.$(OBJ)) $(FILESf:.f=.$(OBJ))

# ------------- TARGETS --------------------

all: $(FILES) BUILD-ISATCK $(CREATE_LIBS)

STREAMS-MOD: streams_mod.$(OBJ)

BUILD-ISATCK:
	(cd $(ISATCK_PATH); $(MAKE))
	touch $@

include ../build-files/build-targets.mk

clean::
	$(RM) BUILD-ISATCK

# ----------- DEPENDENCIES -----------------

include DEPENDS.mk
