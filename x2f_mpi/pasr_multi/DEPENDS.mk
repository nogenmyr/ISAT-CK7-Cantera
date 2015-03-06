# Author: Varun Hiremath <vh63@cornell.edu>
# Last Modified: Fri, 15 Mar 2013 18:03:56 -0500 by Steve Lantz

# List of dependencies

#isat_rnu.$(OBJ): isat_abort_mpi.$(OBJ)
#isat_abort_mpi.$(OBJ): $(MPI_OBJ) isat_prec_dp.$(OBJ) 
#pasr_multi_subs.$(OBJ): isat_rnu.$(OBJ)
pasr_multi.$(OBJ): $(MPI_OBJ) pasr_multi_subs.$(OBJ) 

