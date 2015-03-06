# Author: Varun Hiremath <vh63@cornell.edu>
# Last Modified:  Thu, 12 Nov 2009 16:10:24 -0500

# List of dependencies

isat_cdf.$(OBJ): isat_abort_ser.$(OBJ) isat_prec.$(OBJ)
isat_cdf_action.$(OBJ): isat_cdf.$(OBJ) isat_types.$(OBJ)
isat_change_m.$(OBJ): isat_defaults.$(OBJ) isat_subs.$(OBJ) isat_types.$(OBJ)
isat_defaults.$(OBJ): isat_types.$(OBJ)
isat_init.$(OBJ): isat_cdf_action.$(OBJ) isat_change_m.$(OBJ) isat_defaults.$(OBJ) isat_io.$(OBJ) \
	isat_lu_m.$(OBJ) isat_subs.$(OBJ)
isat_io.$(OBJ): isat_cdf.$(OBJ) isat_subs.$(OBJ) isat_types.$(OBJ)
isat_kill.$(OBJ): isat_cdf.$(OBJ) isat_types.$(OBJ)
isat_lu.$(OBJ): isat_abort_ser.$(OBJ) isat_lu_m.$(OBJ)
isat_mpi_type.$(OBJ): isat_prec.$(OBJ)
isat_rnu.$(OBJ): isat_abort_ser.$(OBJ)
isat_subs.$(OBJ): isat_kill.$(OBJ) isat_rnu.$(OBJ)
isat_types.$(OBJ): isat_abort_ser.$(OBJ) isat_cdf.$(OBJ)
isat_val.$(OBJ): isat_abort_ser.$(OBJ) isat_rnu.$(OBJ)
isatab.$(OBJ): isat_cdf.$(OBJ) isat_init.$(OBJ) isat_io.$(OBJ) isat_kill.$(OBJ) isat_prec.$(OBJ) \
	isat_subs.$(OBJ) isat_val.$(OBJ)
