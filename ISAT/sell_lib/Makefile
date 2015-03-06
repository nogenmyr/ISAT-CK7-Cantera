# Author: Varun Hiremath <vh63@cornell.edu>
# Last Modified:  Thu, 12 Nov 2009 16:32:31 -0500

include ../build-files/vars.mk

RNU_LIB_PATH=../rnu_lib
IFLAGS = -I$(RNU_LIB_PATH)
LIB_NAME = libsell

# ------------- OBJECTS --------------------

include SOURCES.mk

OBJS=$(FILES:.f90=.$(OBJ))

# ------------- TARGETS --------------------

all: $(FILES) BUILD-RNU $(CREATE_LIBS)

BUILD-RNU: 
	(cd $(RNU_LIB_PATH); make)
	touch $@

include ../build-files/build-targets.mk

clean::
	$(RM) BUILD-RNU

# ----------- DEPENDENCIES ------------------

include DEPENDS.mk
