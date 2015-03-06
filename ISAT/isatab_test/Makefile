# Author: Varun Hiremath <vh63@cornell.edu>
# Date: Tue,  4 Sep 2012 16:56:19 -0500

all:
	@echo -e "\n>> BUILDING ISATAB SERIAL TEST CASE";
	make -f Makefile_ser
	@echo -e "\n>> BUILDING ISATAB MPI TEST CASE";
	make -f Makefile_mpi

install-programs:
	make -f Makefile_ser install-programs DESTDIR=$(DESTDIR)
	make -f Makefile_mpi install-programs DESTDIR=$(DESTDIR)

clean:
	make -f Makefile_ser clean
	make -f Makefile_mpi clean
