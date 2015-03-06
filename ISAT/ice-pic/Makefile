# Author: Varun Hiremath <vh63@cornell.edu>
# Last Modified: Sat,  3 Apr 2010 18:02:59 -0400

include ../build-files/vars.mk

CEQ_PATH=../CEQ
ISATAB_PATH=../isatab_ser
ISATCK_PATH=../isat-ck
ISATCK_EXT_PATH=../isat-ck-ext
IFLAGS = -I$(CEQ_PATH) -I$(ISATAB_PATH) -I$(ISATCK_PATH) -I$(ISATCK_EXT_PATH)

LIB_NAME = libice_pic

# ------------- OBJECTS --------------------

include SOURCES.mk

OBJS=$(FILESf90:.f90=.$(OBJ)) $(FILESf:.f=.$(OBJ))

# ------------- TARGETS --------------------

all: $(FILES) BUILD-CEQ BUILD-ISATAB BUILD-ISATCK BUILD-ISATCK-EXT $(CREATE_LIBS)

BUILD-CEQ:
	(cd $(CEQ_PATH); $(MAKE))
	touch $@

BUILD-ISATAB:
	(cd $(ISATAB_PATH); $(MAKE))
	touch $@

BUILD-ISATCK:
	(cd $(ISATCK_PATH); $(MAKE))
	touch $@

BUILD-ISATCK-EXT:
	(cd $(ISATCK_EXT_PATH); $(MAKE))
	touch $@

include ../build-files/build-targets.mk

clean::
	$(RM) BUILD-CEQ BUILD-ISATAB BUILD-ISATCK BUILD-ISATCK-EXT

# ----------- DEPENDENCIES ------------------

include DEPENDS.mk
