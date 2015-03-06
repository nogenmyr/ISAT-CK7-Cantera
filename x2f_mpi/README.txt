x2f_mpi package for parallel function evaluation.
Version 1.1, July 26, 2010
L. Lu, S.R. Lantz, S.B. Pope, Cornell University
(with acknowledgment to D. Rowinski and V. Hiremath)

All rights reserved.  This software may be used for non-commercial 
academic research.  It may not be re-distributed.

This package contains six sub-directories: x2f_mpi, x2f_mpi_isat, build-files,
               test_x2f_mpi, pasr_multi,and test_files_pasr_multi.

x2f_mpi      : contains all the source codes for building the x2f_mpi library.
               The overview of x2f_mpi library can be found in ISAT_MP.f90,
               which exists in both test_x2f_mpi and x2f_mpi_isat directory.

x2f_mpi_isat : contains source for building a library of interface routines
               (e.g., ISAT_MP) designed specially for users who are interested
               in computations of chemical kinetics using In Situ Adaptive 
               Tabulation (ISAT) as the "x2f" function.  Therefore the program
               uses not only the x2f_mpi library but also the ISAT library.

build-files  : contains lists of targets and variables needed during the build		   
               processes for Linux and Windows.

test_x2f_mpi : contains a small program for testing the correctness of the
               x2f_mpi library.

pasr_multi   : contains PaSR (Partially Stirred Reactor) source codes.

pasr_multi_test_files : has input files, instructions, and comparison output
                        to be used in conjunction with the pasr_multi test.

Besides the directories described above, there are 3 more files: a Makefile
for Linux systems; an nmake-based Makefile.windows for Windows systems; and
(for ISAT) the pdf file flowchart.pdf containing a diagram showing how all the 
blocks of code are linked -- PaSR, ISAT_MP, x2f_mpi, ISAT-CK, ISATAB etc. 

You must have a Fortran 90/95 compiler to build the x2f_mpi library and the
related test executables.  Pre-build setup steps are discussed below.  To
initiate a build on Windows systems, open up a command prompt and issue the
following command (make sure the Microsoft "nmake" utility is in your path):

    nmake -f Makefile.windows

On Linux systems, enter:

    make

To clean up the library and executables, issue the Windows command:

    nmake -f Makefile.windows clean

For Linux systems, use:

    make clean

Setup notes:

1/ Before compiling, you must modify the compiler name and/or library paths
   in the vars, vars.windows, and/or PATHS.windows files in the build-files
   directory. For the compiler name, mpif90 (if it exists) is usually safest.

2/ If your local (Windows) MPI implementation has an mpi.f90 or mpif.f90 file,
   set MPI_FILE to its name and MPI_PATH to its full path without the name.
   These variables are defined in PATHS.windows. If your MPI implementation
   has an mpif.h instead of mpi.f90, create a file called mpi.f90 in the 
   x2f_mpi directory with the following contents:

      module MPI
         include "mpif.h"
      end module MPI

   In this case MPI_FILE and MPI_PATH should be set to the name and full path
   of the above mpi.f90, while MPI_INC should point to the include directory
   containing the mpif.h.

3/ If you do not wish to use x2f_mpi with ISAT, remove "pasr" from the "all"
   dependency list in makeinstall (or makeinstall_Linux).

>>> the remaining setup steps pertain only to PaSR and ISAT >>>

4/ To run pasr_multi, you must provide the files chem.bin and streams.in.

5/ In the pasr_multi test, isat_1.nml and pasr.nml can be used to change ISAT 
   parameters and PaSR parameters, respectively.

6/ In parallel simulations, file names of chem.bin, streams.in and isat_1.nml
   need to be modified to the form chem_#_P.bin, streams_P.in and
   isat_#_P.nml. Here # denotes the table number and P is the MPI rank of the
   process. If there is only one table on each processor, the file names will
   look like chem_1_P.bin and isat_1_P.nml.
