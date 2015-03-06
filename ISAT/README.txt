# Note by Kalle: See ChangeLogCantera.txt!

# Author: Varun Hiremath <vh63@cornell.edu>
# Date: Fri, 30 Nov 2012 13:21:04 -0600

ISAT-CK7
========

ISAT-CK7 is a collection of the following projects:

Libraries
---------
CEQ
ell_lib
ice-pic
isatab_mpi
isatab_ser
isat-ck
isat-ck-ext
rnu_lib
sell_lib

Programs
--------
PaSR

Makefiles are provided to build all the projects on Linux and Windows.

Dependencies
------------

ISAT-CK7 requires the following external libraries for full
functionality:

ADIFOR
CHEMKIN (including ckinterp)
LAPACK

More information on how to link to these libraries is provided in
README_dependencies.txt

HOWTO BUILD
===========

Linux
-----

To build all the sub-projects (libraries and executables) on Linux
run the following command from project's root directory:

    $ make

This would build all the libraries and copy the built (shared *.so and
static *.a) libraries to ${DESTDIR}/lib and the executables (ckinterp
and PaSR) to ${DESTDIR}/bin.

The DESTDIR is by default set to ${PWD} on Linux but may be changed
by passing the desired DESTDIR while invoking `make` like this: 

    $ make DESTDIR=<install path>

If the DESTDIR is changed, the LIBS_PATH variable in PaSR/Makefile
must be set to the same DESTDIR/lib before building PaSR, since PaSR
needs to link to these libraries for building.

Windows
-------

To build all the sub-projects (libraries and executables) on Windows
run the following command from project's root directory:

    $ nmake -f Makefile.windows

This would build all the libraries and copy the built (static *.lib)
libraries to ${DESTDIR}/lib and the executables (ckinterp and PaSR) to
${DESTDIR}/bin.

The DESTDIR is by default set to the current working directory (.) but
may be changed by editing the Makefile.windows. With the default
settings, after the build all the libraries must be present in .\lib
and executables in .\bin

If the DESTDIR is changed, the LIBS_PATH variable in PaSR/Makefile
must be set to the same DESTDIR\lib before building PaSR, since PaSR
needs to link to these libraries for building.

MPI Build
---------

Build of MPI codes (isatab_mpi) on Linux uses intel's mpif90. The
Makefiles for Linux are configured to use mpif90. On the Lonestar
cluster, before calling make please ensure that you have loaded the
mkl module (the intel module should be loaded by default):

   $ module load mkl
   $ make DESTDIR=<install path>

To build isatab_mpi on windows, the path to mpi.f90 file needs to be
specified in build-files/MPI.windows before calling make.


TEST CASES
==========

The test_cases folder provides example test cases to demonstrate the
use of isatab and PaSR. Please read the README.txt files provided with
the test cases for more information.
