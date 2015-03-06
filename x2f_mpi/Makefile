# Makefile to build all the projects on Linux
# Author: Varun Hiremath <vh63@cornell.edu>
# Last Modified:   Mon Jan 31 2011 17:05:34 -0400 by Steve Lantz

PROJECTS=x2f_mpi x2f_mpi_isat
PROGRAMS=test_x2f_mpi test_sort_tools pasr_multi

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

build: build-libs build-programs

build-libs:
	@echo -e "\n>> BUILDING LIBRARIES";
	@for project in $(PROJECTS); \
	do \
		echo -e "\n>>> BUILDING: $$project"; \
		cd $$project; $(MAKE); test "$$?" != "0" && exit 1; cd ..; \
	done

build-programs: install-libs clean
	@echo -e "\n>> BUILDING PROGRAMS";
	@for program in $(PROGRAMS); \
	do \
		echo -e "\n>>> BUILDING: $$program"; \
		cd $$program; $(MAKE); test "$$?" != "0" && exit 1; cd ..; \
	done

install: build install-libs install-programs

install-libs:
	@mkdir -p $(DESTDIR)/lib
	@for project in $(PROJECTS); \
		do (cd $$project; cp *.so $(DESTDIR)/lib; cp *.a $(DESTDIR)/lib;); done
	@echo -e "\n>>> BUILDING LIBRARIES DONE";
	@echo -e "\n>>> THE FOLLOWING LIBS HAVE BEEN GENERATED IN $(DESTDIR)/lib";
	@ls $(DESTDIR)/lib

install-programs:
	@mkdir -p $(DESTDIR)/bin
	@for program in $(PROGRAMS); \
		do (cd $$program; cp $$program $(DESTDIR)/bin;); done
	@echo -e "\n>> BUILDING PROGRAMS DONE";
	@echo -e "\n>>> THE FOLLOWING EXECUTABLES HAVE BEEN GENERATED IN $(DESTDIR)/bin";
	@ls $(DESTDIR)/bin
	@echo -e "\n"

clean: clean-programs

dist-clean: clean-libs clean-programs

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

rebuild: clean build

.PHONY: clean
