# Author: Varun Hiremath <vh63@cornell.edu>
# Last Modified:  Thu,  8 Jul 2010 17:43:24 -0400

# List of source files to build

FILESf90 = sort_tools.f90 x2f_check_atmpts.f90 x2f_loclev.f90		\
	x2f_mpi_bksize.f90 x2f_mpi.f90 x2fmpi_file_name.f90		\
	x2f_mpi_nproc.f90 x2f_pref.f90 x2f_strata.f90 x2f_uran.f90

FILESc = buildGraph.c glib.c match.c trymatch.c

FILESh = graphtypes.h 

FILES = $(FILESf90) $(FILESc) $(FILESh)
