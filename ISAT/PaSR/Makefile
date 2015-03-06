# Author: Varun Hiremath <vh63@cornell.edu>
# Last Modified: Wed Nov 28 10:18:09 EST 2012 by Steve Lantz

include ../build-files/vars.mk

ISATCK_PATH = ../isat-ck
ICEPIC_PATH = ../ice-pic
ISATAB_PATH = ../isatab_ser

IFLAGS = -I$(ISATAB_PATH) -I$(ISATCK_PATH) -I$(ICEPIC_PATH) -I$(CANTERA_INCPATH)
LIBS_PATH = $(DESTDIR)/lib
LDFLAGS= -L$(LIBS_PATH) -lisat7_ser -lck_ext $(MKL_LIST) -L$(CANTERA_LIBPATH) -lcantera_fortran -lcantera -lpthread -lstdc++ 

#$(CANTERA_LIBPATH)/libcantera_fortran.a

NAME = PaSR

# ------------- OBJECTS --------------------

include SOURCES.mk

OBJS = $(FILES:.f=.$(OBJ))

# ------------- TARGETS --------------------

all: $(FILES) $(NAME)

include ../build-files/build-targets.mk
