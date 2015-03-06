INFORMATION
===========
NAME: CTLIB
DESCRIPTION: A library of gas-phase kinetics subroutines using Cantera

AUTHOR: Karl-Johan Nogenmyr

LICENSE:
 This is from the "public domain" version of Chemkin developed at
 Sandia. This code may be freely used within TCG. However, because of
 Sandia's license to Reaction Design, the code must not be
 redistributed.

BUILD
=====

PREREQUISITES
-------------
1/ ifort (intel fortran compiler, any version >= 9) compiler
2/ make (on Linux) and nmake (on windows) for building source code

Linux
-----
On linux, in the source directory run:
  $ make

This would build the following libraries:

1/ libct.so -> dynamic library
2/ libct.a  -> static library

Windows
-------
On windows, in the source directory run:
   $ nmake -f Makefile.windows

This would build the following library:

1/ ct.lib -> static library
