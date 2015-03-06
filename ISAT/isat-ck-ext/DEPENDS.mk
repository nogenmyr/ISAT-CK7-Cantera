# Author: Varun Hiremath <vh63@cornell.edu>
# Last Modified:  Fri, 12 Mar 2010 17:05:38 -0500

# List of dependencies

streams_mod.$(OBJ):
ci_subs_ext.$(OBJ): streams_mod.$(OBJ)
ci_ext_routines.$(OBJ): streams_mod.$(OBJ)
