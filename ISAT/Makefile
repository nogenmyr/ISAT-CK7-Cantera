# Makefile to build all the projects on Linux
# Author: Varun Hiremath <vh63@cornell.edu>
# Last Modified: Tue Jan 20 09:37:26 CET 2015

PROJECTS=ADIFOR_stub CEQ canteralib cklib_ext_stub ell_lib isatab_mpi isatab_ser isat-ck isat-ck-ext ice-pic rnu_lib sell_lib # removed: 
DEP_LIBS=libadifor libct libck_ext liblapack libcantera_fortran

PROGRAMS= # removed PaSR isatab_test
# DEP_PROGRAMS=ckinterp  # Not needed, refer to ck2cti
ACTIVE= ice-pic  # removed: isat-ck PaSR

ISAT_SER=ADIFOR_stub CEQ canteralib cklib_ext_stub ell_lib isatab_ser isat-ck isat-ck-ext ice-pic rnu_lib sell_lib # removed:
ISAT_MPI=ADIFOR_stub CEQ canteralib cklib_ext_stub ell_lib isatab_mpi isat-ck isat-ck-ext ice-pic rnu_lib sell_lib # removed: 

DESTDIR=${PWD}

# Set build type
MAKE=make BUILD_TYPE=$(BUILD_TYPE)

include build-files/vars.mk

all: install
ifneq ($(BUILD_TYPE), debug)
	@$(RM) DEBUG-MODE
endif

debug:
	make BUILD_TYPE=debug
	@echo "THE CODE HAS BEEN BUILT UNDER THE DEBUG MODE" > DEBUG-MODE

build: build-libs build-isat-libs build-programs

build-libs:
	@echo -e "\n>> BUILDING LIBRARIES";
	@for project in $(PROJECTS); \
	do \
		echo -e "\n>>> BUILDING: $$project"; \
		cd $$project; $(MAKE); test "$$?" != "0" && exit 1; cd ..; \
	done

build-isat-libs: isat7_ser isat7_mpi

isat7_ser:
	@echo -e "\n>>> BUILDING ISAT7 SERIAL LIBRARY";
	@$(RMDIR) tmp_$@
	@mkdir -p tmp_$@
	@for project in $(ISAT_SER); \
		do (cp $$project/*.a tmp_$@;); done
	@for dep_lib in $(DEP_LIBS); \
	do \
		if [ -f ../ISAT_dependencies/lib/$$dep_lib.a ]; \
		then \
			cp -p ../ISAT_dependencies/lib/$$dep_lib.a tmp_$@; \
			echo "Using ../ISAT_dependencies/lib/$$dep_lib.a"; \
		fi; \
	done

	@for alib in tmp_$@/*.a; \
		do (cd tmp_$@; ar x ../$$alib;); done \

	@(cd tmp_$@; \
	$(AR) lib$@.a *.o; \
	$(FC) -shared -o lib$@.so *.o)
	@mkdir -p $(DESTDIR)/lib
	@cp tmp_$@/lib$@.* $(DESTDIR)/lib
	@$(RMDIR) tmp_$@

isat7_mpi:
	@echo -e "\n>>> BUILDING ISAT7 MPI LIBRARY";
	@$(RMDIR) tmp_$@
	@mkdir -p tmp_$@
	@for project in $(ISAT_MPI); \
		do (cp $$project/*.a tmp_$@;); done
	@for dep_lib in $(DEP_LIBS); \
	do \
		if [ -f ../ISAT_dependencies/lib/$$dep_lib.a ]; \
		then \
			cp -p ../ISAT_dependencies/lib/$$dep_lib.a tmp_$@; \
			echo "Using ../ISAT_dependencies/lib/$$dep_lib.a"; \
		fi; \
	done

	@for alib in tmp_$@/*.a; \
		do (cd tmp_$@; ar x ../$$alib;); done \

	@(cd tmp_$@; \
	$(AR) lib$@.a *.o; \
	$(FC) -shared -o lib$@.so *.o)
	@mkdir -p $(DESTDIR)/lib
	@cp tmp_$@/lib$@.* $(DESTDIR)/lib
	@$(RMDIR) tmp_$@

build-programs: install-libs clean
	@echo -e "\n>> BUILDING PROGRAMS";
	@for program in $(PROGRAMS); \
	do \
		echo -e "\n>>> BUILDING: $$program"; \
		cd $$program; $(MAKE) DESTDIR=$(DESTDIR); test "$$?" != "0" && exit 1; cd ..; \
	done

install: build install-libs install-programs

install-libs:
	@mkdir -p $(DESTDIR)/lib
	@echo -e "\n>>> INSTALLING LIBRARIES to $(DESTDIR)/lib";
	@for project in $(PROJECTS); \
		do (cd $$project; make install-libs DESTDIR=$(DESTDIR);); done
	@echo -e "\n>>> INSTALLING ISAT DEPENDENCIES to $(DESTDIR)/lib";
	@for dep_lib in $(DEP_LIBS); \
	do \
		if [ -f ../ISAT_dependencies/lib/$$dep_lib.a ]; \
		then \
			rm $(DESTDIR)/lib/$$dep_lib.*; \
			cp -p ../ISAT_dependencies/lib/$$dep_lib.a $(DESTDIR)/lib; \
		else \
			echo -e "\n    Warning: stub libs have been installed for $$dep_lib"; \
		fi; \
		if [ -f ../ISAT_dependencies/lib/$$dep_lib.so ]; \
		then \
			cp -p ../ISAT_dependencies/lib/$$dep_lib.so $(DESTDIR)/lib; \
		fi; \
	done
	@echo -e "\n>>> INSTALLING LIBRARIES DONE";
	@echo -e "\n>>> THE FOLLOWING LIBS HAVE BEEN GENERATED IN $(DESTDIR)/lib";
	@ls $(DESTDIR)/lib

install-programs:
	@mkdir -p $(DESTDIR)/bin
	@echo -e "\n>>> INSTALLING PROGRAMS to $(DESTDIR)/bin";
	@for program in $(PROGRAMS); \
		do (cd $$program; make install-programs DESTDIR=$(DESTDIR);); done
	@echo -e "\n>>> INSTALLING ISAT DEPENDENCIES to $(DESTDIR)/bin";
	@for dep_program in $(DEP_PROGRAMS); \
	do \
		if [ -f ../ISAT_dependencies/bin/$$dep_program ]; \
		then \
			cp -p ../ISAT_dependencies/bin/$$dep_program $(DESTDIR)/bin; \
		fi; \
	done
	@echo -e "\n>> INSTALLING PROGRAMS DONE";
	@echo -e "\n>>> THE FOLLOWING EXECUTABLES HAVE BEEN GENERATED IN $(DESTDIR)/bin";
	@ls $(DESTDIR)/bin
	@echo -e "\n"

clean: clean-programs

dist-clean: clean-libs clean-programs
	@$(RM) DEBUG-MODE

clean-libs:
	@for project in $(PROJECTS); \
	do \
		echo -e "\n>>> CLEANING: $$project"; \
		(cd $$project; make clean); \
	done

clean-programs:
	@for program in $(PROGRAMS); \
	do \
		echo -e "\n>>> CLEANING: $$program"; \
		(cd $$program; make clean); \
	done

clean-active:
	@for program in $(ACTIVE); \
	do \
		echo -e "\n>>> CLEANING: $$program"; \
		(cd $$program; make clean); \
	done

rebuild: clean build

.PHONY: clean
