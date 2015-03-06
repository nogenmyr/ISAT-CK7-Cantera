# Author: Varun Hiremath <vh63@cornell.edu>
# Last Modified:  Tue,  9 Mar 2010 19:52:59 -0500

# List of dependencies

ci_8.$(OBJ): ci_dat8.$(OBJ) ci_ice_cksubs.$(OBJ) ci_utils.$(OBJ)
ci_all_errors.$(OBJ): ci_cem_recon.$(OBJ) ci_dat8.$(OBJ) ci_ice_cksubs.$(OBJ) ci_utils.$(OBJ)
ci_cem_recon.$(OBJ): ci_dat8.$(OBJ) ci_utils.$(OBJ)
ci_cem_species_dr.$(OBJ): ci_cem_recon.$(OBJ) ci_dat8.$(OBJ) ci_ice_cksubs.$(OBJ) ci_utils.$(OBJ)
ci_cem_species_rm.$(OBJ): ci_cem_recon.$(OBJ) ci_dat8.$(OBJ) ci_ice_cksubs.$(OBJ) ci_utils.$(OBJ)
ci_cem_species_sr.$(OBJ): ci_cem_recon.$(OBJ) ci_dat8.$(OBJ) ci_ice_cksubs.$(OBJ) ci_utils.$(OBJ)
ci_cem_species_rm_gali.$(OBJ): ci_cem_recon.$(OBJ) ci_dat8.$(OBJ) ci_ice_cksubs.$(OBJ) ci_utils.$(OBJ)
ci_cem_tan.$(OBJ): ci_utils.$(OBJ)
ci_chem_one_step.$(OBJ): ci_utils.$(OBJ)
ci_dat8.$(OBJ): 
ci_dpt_dr.$(OBJ): ci_dat8.$(OBJ)
ci_dr_subs.$(OBJ): ci_8.$(OBJ) ci_dat8.$(OBJ) ci_dpt_dr.$(OBJ)
ci_ice_chem_map.$(OBJ): ci_dat8.$(OBJ) ci_ice_cksubs.$(OBJ) ci_utils.$(OBJ)
ci_ice_cksubs.$(OBJ): ci_dat8.$(OBJ) ci_dpt_dr.$(OBJ)
ci_ice_facet.$(OBJ): ci_dat8.$(OBJ)
ci_ice_facet_newton.$(OBJ): ci_cem_recon.$(OBJ) ci_dat8.$(OBJ) ci_ice_cksubs.$(OBJ) ci_utils.$(OBJ)
ci_ice_facet_newton_lmap.$(OBJ): ci_cem_recon.$(OBJ) ci_dat8.$(OBJ) ci_ice_cksubs.$(OBJ) ci_utils.$(OBJ)
ci_ice_min_norm.$(OBJ): ci_dat8.$(OBJ) ci_utils.$(OBJ)
ci_ice_pic_bound.$(OBJ): ci_dat8.$(OBJ) ci_utils.$(OBJ)
ci_ice_pic_newton.$(OBJ): ci_cem_recon.$(OBJ) ci_dat8.$(OBJ) ci_ice_cksubs.$(OBJ) ci_utils.$(OBJ)
ci_ice_pic_newton_lmap.$(OBJ): ci_cem_recon.$(OBJ) ci_dat8.$(OBJ) ci_ice_cksubs.$(OBJ) ci_utils.$(OBJ)
ci_ice_pic_ray.$(OBJ): ci_cem_recon.$(OBJ) ci_dat8.$(OBJ) ci_ice_cksubs.$(OBJ) ci_utils.$(OBJ)
ci_ice_pim_tan.$(OBJ): ci_dat8.$(OBJ) ci_utils.$(OBJ)
ci_ice_recon.$(OBJ): ci_cem_recon.$(OBJ) ci_dat8.$(OBJ) ci_utils.$(OBJ)
ci_ice_recon_bc.$(OBJ): ci_cem_recon.$(OBJ) ci_dat8.$(OBJ) ci_ice_cksubs.$(OBJ) ci_utils.$(OBJ)
ci_ice_recon_bc_lmap.$(OBJ): ci_cem_recon.$(OBJ) ci_dat8.$(OBJ) ci_ice_cksubs.$(OBJ) ci_utils.$(OBJ)
ci_ice_recon_pic.$(OBJ): ci_cem_recon.$(OBJ) ci_dat8.$(OBJ) ci_ice_cksubs.$(OBJ) ci_utils.$(OBJ)
ci_ice_recon_pic_lmap.$(OBJ): ci_cem_recon.$(OBJ) ci_dat8.$(OBJ) ci_ice_cksubs.$(OBJ) ci_utils.$(OBJ)
ci_ice_species_rm.$(OBJ): ci_8.$(OBJ) ci_cem_recon.$(OBJ) ci_dat8.$(OBJ) ci_ice_cksubs.$(OBJ) ci_utils.$(OBJ)
ci_ice_species_sr.$(OBJ): ci_8.$(OBJ) ci_cem_recon.$(OBJ) ci_dat8.$(OBJ) ci_ice_cksubs.$(OBJ) ci_utils.$(OBJ)
ci_ice_tan.$(OBJ): ci_ice_cksubs.$(OBJ) ci_utils.$(OBJ)
ci_linrmap.$(OBJ): ci_dat8.$(OBJ) ci_ice_cksubs.$(OBJ) ci_utils.$(OBJ)
ci_recon.$(OBJ): ci_8.$(OBJ) ci_cem_recon.$(OBJ) ci_dat8.$(OBJ)
ci_rzmap.$(OBJ): ci_cem_recon.$(OBJ) ci_dat8.$(OBJ) ci_utils.$(OBJ)
ci_test_dpt.$(OBJ): ci_dat8.$(OBJ) ci_cem_recon.$(OBJ)
ci_test_grad.$(OBJ): ci_8.$(OBJ) ci_utils.$(OBJ)
ci_test_interface.$(OBJ): ci_utils.$(OBJ)
ci_test_jacobian.$(OBJ): ci_cem_recon.$(OBJ) ci_dat8.$(OBJ) ci_utils.$(OBJ)
ci_test_linrmap.$(OBJ): ci_cem_recon.$(OBJ) ci_dat8.$(OBJ) ci_utils.$(OBJ)
ci_test_rzmap.$(OBJ): ci_8.$(OBJ) ci_cem_recon.$(OBJ) ci_dat8.$(OBJ) ci_ice_cksubs.$(OBJ)
ci_test_rxn_map.$(OBJ): ci_dat8.$(OBJ) ci_utils.$(OBJ) ci_ice_cksubs.$(OBJ)
ci_test_skeletal.$(OBJ): ci_dat8.$(OBJ) ci_utils.$(OBJ) ci_ice_cksubs.$(OBJ)
ci_test_detailed.$(OBJ): ci_dat8.$(OBJ) ci_utils.$(OBJ) ci_ice_cksubs.$(OBJ)
ci_utils.$(OBJ): 
cirmap_dr_ice_pic.$(OBJ): ci_dat8.$(OBJ) ci_ice_cksubs.$(OBJ) ci_utils.$(OBJ) ci_dpt_dr.$(OBJ)
cirmap_dr_rcce.$(OBJ): ci_cem_recon.$(OBJ) ci_dat8.$(OBJ) ci_utils.$(OBJ) ci_dpt_dr.$(OBJ)
ddasac5.$(OBJ): 
g_jacadf.$(OBJ): 
gecp.$(OBJ): 
