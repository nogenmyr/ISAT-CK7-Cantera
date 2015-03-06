# Author: Varun Hiremath <vh63@cornell.edu>
# Last Modified:  Thu, 12 Nov 2009 16:17:56 -0500

# List of dependencies

ci_1.$(OBJ): ci_dat.$(OBJ)
ci_2.$(OBJ): ci_dat.$(OBJ)
ci_6.$(OBJ): ci_ck.$(OBJ) ci_dat.$(OBJ) ci_rmap.$(OBJ)
ci_ck.$(OBJ): ci_dat.$(OBJ) ci_dat6.$(OBJ)
ci_cksubs.$(OBJ): ci_dat6.$(OBJ)
ci_dasubs.$(OBJ): ci_cksubs.$(OBJ) ci_dat6.$(OBJ)
ci_dat.$(OBJ): ci_prec.$(OBJ)
ci_dat6.$(OBJ): ci_prec.$(OBJ) ckstrt.h
ci_rmap.$(OBJ): ci_cksubs.$(OBJ) ci_stats.$(OBJ)
ci_sens.$(OBJ): ci_dat6.$(OBJ) ci_stats.$(OBJ)
ci_stats.$(OBJ): ci_dat.$(OBJ) ci_prec.$(OBJ)
ci_subs.$(OBJ): ci_1.$(OBJ) ci_2.$(OBJ) ci_6.$(OBJ) ci_cksubs.$(OBJ) ci_dat.$(OBJ) ci_dat6.$(OBJ) ci_prec.$(OBJ)
cireal.$(OBJ): ci_cksubs.$(OBJ)
cirmap1.$(OBJ): ci_dat.$(OBJ) ci_rmap.$(OBJ) ci_sens.$(OBJ)
