# Author: Varun Hiremath <vh63@cornell.edu>
# Last Modified:  Thu, 12 Nov 2009 16:32:52 -0500

# List of dependencies

id_list.$(OBJ): 
sell_bt.$(OBJ): sell_m.$(OBJ)
sell_ebt.$(OBJ): sell_m.$(OBJ)
sell_ll.$(OBJ): sell_m.$(OBJ)
sell_m.$(OBJ): id_list.$(OBJ)

