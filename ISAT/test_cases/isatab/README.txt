# Author: Varun Hiremath <vh63@cornell.edu>
# Date: Thu, 29 Nov 2012 15:28:01 -0600

ISATAB: In Situ Adaptive Tabulation
===================================
Test cases are provided to run ISATAB in serial and parallel (MPI).

Executables
-----------
Use 'test_ser' and 'test_mpi' executables built using the isatab_test
source code provided with ISAT-CK7.

Serial Test Case
----------------
Copy the 'test_ser' executable into the 'test_case_ser' directory
provided, which contains three isat_#.nml files for the three test
cases explained in the test_ser.f90 source file (and listed below).

Parallel Test Case
------------------
Copy the 'test_mpi' executable into the 'test_case_mpi' directory
provided, which contains three isat_#.nml files for the three test
cases explained in the test_mpi.f90 source file (and listed below).

A sample job script 'job.sh' is provided to run the parallel test case
on 16 cores on the TACC Ranger cluster.

Test Cases
----------
As listed in the source files:

Table 1: 3D, linear, fixed domain
Table 2: 5D, non-linear, fixed domain (40Mbytes), ifull=1
Table 3: 5D, non-linear, shifting domain

Input Parameters
----------------
The following parameters can be set using the isat_#.nml input file:
etola : ISAT error tolerance
stomby: ISAT table size

Expected Output
---------------
The output files generated from the serial test case are included in
the 'test_case_ser_output' directory; and the output files from the
parallel test case are included in the 'test_case_mpi_output'
directory. 
