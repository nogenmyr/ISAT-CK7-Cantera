# Author: Varun Hiremath <vh63@cornell.edu>
# Last Modified:  Thu, 12 Nov 2009 16:07:01 -0500

# List of dependencies

ceq_linprog.$(OBJ): linprog.$(OBJ)
ceq_solve.$(OBJ): ceq_types.$(OBJ)
ceq_state.$(OBJ): ceq_solve.$(OBJ)
ceq_system.$(OBJ): ceq_state.$(OBJ)
