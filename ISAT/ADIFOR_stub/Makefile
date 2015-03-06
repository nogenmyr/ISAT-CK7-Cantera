# Author: Varun Hiremath <vh63@cornell.edu>
# Last Modified: Wed Nov 28 14:27:51 EST 2012 by Steve Lantz

include ../build-files/vars.mk

ISATAB_PATH=../isatab_ser
IFLAGS = -I$(ISATAB_PATH)

LIB_NAME = libadifor

# ------------- OBJECTS --------------------

include SOURCES.mk

OBJS=$(FILESf90:.f90=.$(OBJ)) $(FILESf:.f=.$(OBJ))

# ------------- TARGETS --------------------

all: $(FILES) BUILD-ISATAB $(CREATE_LIBS)

BUILD-ISATAB:
	(cd $(ISATAB_PATH); $(MAKE))
	touch $@

include ../build-files/build-targets.mk

clean::
	$(RM) BUILD-ISATAB
