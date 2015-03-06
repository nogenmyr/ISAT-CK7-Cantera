# Author: Varun Hiremath <vh63@cornell.edu>
# Last Modified:  Thu, 12 Nov 2009 16:17:01 -0500

# List of source files to build

FILESf90 = ci_sens.f90 ci_stats.f90

FILESf = ci_1.f ci_2.f ci_6.f ci_ck.f ci_ck_flag.f ci_cksubs.f		\
	 ci_dasubs.f ci_dat.f ci_dat6.f ci_prec.f ci_rmap.f ci_subs.f	\
	 cireal.f cirmap1.f

FILES = $(FILESf) $(FILESf90)

# Left out: g_jacadf.f (not ported to Cantera)
