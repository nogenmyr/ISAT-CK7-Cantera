# Author: Varun Hiremath <vh63@cornell.edu>
# Last Modified:  Thu, 18 Feb 2010 19:17:28 -0500

# List of source files to build

FILESf90 = ci_chem_one_step.f90 ci_recon.f90 ci_cem_recon.f90		\
	   ci_ice_facet_newton.f90 ci_ice_facet_newton_lmap.f90		\
	   ci_ice_recon_pic.f90 ci_ice_recon_pic_lmap.f90		\
	   ci_ice_pic_bound.f90 ci_ice_pic_newton.f90			\
	   ci_ice_pic_newton_lmap.f90 ci_ice_pic_ray.f90		\
	   ci_ice_pim_tan.f90 ci_linrmap.f90 ci_ice_recon_bc.f90	\
	   ci_ice_recon_bc_lmap.f90 ci_ice_recon.f90 ci_test_grad.f90	\
	   ci_test_jacobian.f90 ci_test_linrmap.f90 ci_all_errors.f90	\
	   ci_cem_species_sr.f90 ci_cem_species_rm.f90			\
	   ci_cem_species_dr.f90 ci_rzmap.f90 ci_test_rzmap.f90		\
	   ci_ice_species_sr.f90 ci_ice_species_rm.f90 ci_utils.f90	\
	   ci_dr_subs.f90 ci_test_interface.f90 ci_test_dpt.f90		\
	   ci_dpt_dr.f90 ci_cem_species_rm_gali.f90			\
	   ci_test_rxn_map.f90 ci_test_skeletal.f90			\
	   ci_test_detailed.f90


FILESf = ci_cem_tan.f ci_ice_tan.f ci_ice_chem_map.f ci_ice_cksubs.f	\
	 ci_ice_min_norm.f cirmap_dr_ice_pic.f cirmap_dr_rcce.f		\
	 ddasac5.f gecp.f ci_ice_facet.f ci_8.f ci_dat8.f

FILES = $(FILESf) $(FILESf90)

# removed g_jacadf.f (not ported to Cantera) 
