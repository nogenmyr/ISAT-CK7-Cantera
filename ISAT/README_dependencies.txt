Full ISAT-CK7 capability is available only when the provided source
code is linked against the libraries (listed below) that belong to
other copyright holders. If any of these third-party libraries is
unavailable to you, the code will instead be linked against a
corresponding "stub" library, provided along with the ISAT-CK7 source
code. The stub libraries will satisfy linker requirements, but they do
not provide true functionality. ISAT-CK7 will still work, but some of
its modes will merely return an error if you attempt to use them.


The ISAT-CK7 makefiles will look for a directory
../ISAT_dependencies/lib, where the path is relative to the main ISAT
directory. Libraries found in this directory will be used
preferentially over the stub libraries.

ADIFOR
------
ISAT-CK7 relies on ADIFOR to compute an analytical Jacobian.
If you have access to the ADIFOR source code, this mode can become
available. Note that it relies only on the subroutine ehbfdo.

To link ISAT-CK7 to ADIFOR, you need to provide a static library for
ADIFOR named "libadifor.a" and place it in ../ISAT_dependencies/lib/.

CANTERA
-------
This version of ISAT-CK7 relies on this software to compute chemical
reactions. This mode can become fully available if the fortran
interface library (libcantera_fortran.a) of Cantera is placed in 
../ISAT_dependencies/lib/. This may be a symbolic link.

CHEMKIN EXTENSIONS
------------------
The Turublence and Combustion Group (TCG) at Cornell has written
extensions to CHEMKIN to enhance its performance. This library is
completely optional, and the use of the "cklib_ext_stub" stub library
is adequate in this case.

LAPACK
------
This is the well-known linear algebra package from netlib.org. Many
implementations (e.g., ATLAS) are freely available. LAPACK can be
downloaded and built from source, in case it is not otherwise
available on your platform.

Currently, ISAT-CK7 is made to link with the Intel MKL LAPACK
libraries on TACC clusters. The include path and libraries are listed
in the file "ISAT/build-files/vars.mk". Update the variables
TACC_MKL_LIB and MKL_LIST to link to LAPACK on your system.



