INFORMATION
===========
NAME: CKLIB
DESCRIPTION: A library of gas-phase kinetics subroutines for the
 collection of codes known as CHEMKIN
VERSION: 4.5
AUTHOR: Sandia National Laboratories
Homepage: http://www.ca.sandia.gov/chemkin/
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

1/ libck.so -> dynamic library
2/ libck.a  -> static library

Windows
-------
On windows, in the source directory run:
   $ nmake -f Makefile.windows

This would build the following library:

1/ ck.lib -> static library
