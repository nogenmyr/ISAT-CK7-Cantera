# Author: Varun Hiremath <vh63@cornell.edu>
# Last Modified:  Thu,  8 Jul 2010 17:44:09 -0400

# List of dependencies

x2f_mpi.$(OBJ): $(MPI_OBJ) x2f_uran.$(OBJ) x2f_loclev.$(OBJ) x2f_strata.$(OBJ) x2f_pref.$(OBJ) sort_tools.$(OBJ) 
x2f_uran.$(OBJ): $(MPI_OBJ)
x2f_pref.$(OBJ): $(MPI_OBJ)
x2f_loclev.$(OBJ): $(MPI_OBJ)
x2f_strata.$(OBJ): $(MPI_OBJ)
x2f_mpi_nproc.$(OBJ): $(MPI_OBJ)
x2f_mpi_bksize.$(OBJ): $(MPI_OBJ)
