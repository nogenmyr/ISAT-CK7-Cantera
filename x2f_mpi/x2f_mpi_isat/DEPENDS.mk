# Author: Varun Hiremath <vh63@cornell.edu>
# Date: Fri, 11 May 2012 15:49:30 -0500

# List of dependencies

ISAT_MP.$(OBJ): x2f_mpi_grpdata.$(OBJ) x2f_mpi_adapt.$(OBJ)
x2f_mpi_grpdata.$(OBJ): $(MPI_OBJ)
x2f_mpi_adapt.$(OBJ): $(MPI_OBJ)
x2f_rw_isatab.$(OBJ): x2f_mpi_grpdata.$(OBJ)
