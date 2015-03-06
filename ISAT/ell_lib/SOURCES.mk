# Author: Varun Hiremath <vh63@cornell.edu>
# Last Modified:  Thu, 12 Nov 2009 16:26:17 -0500

# List of source files to build

FILESf90 = ell_aff_pr.f90 ell_bbt2chol.f90 ell_bbt2eig.f90		\
	ell_chol2eig.f90 ell_chol_det.f90 ell_cov_2d.f90		\
	ell_eig2chol.f90 ell_full2eig.f90 ell_full2low.f90		\
	ell_house.f90 ell_line_proj.f90 ell_low2chol.f90		\
	ell_pair_cover_cv.f90 ell_pair_cover_cv_ql.f90			\
	ell_pair_cover.f90 ell_pair_cover_it.f90			\
	ell_pair_cover_query.f90 ell_pair_cover_sph.f90			\
	ell_pair_separate.f90 ell_pair_shrink.f90 ell_pt_dist.f90	\
	ell_pt_hyper.f90 ell_pt_in.f90 ell_pt_modify.f90		\
	ell_pt_near_far.f90 ell_pt_shrink.f90 ell_pts_uncover.f90	\
	ell_radii.f90 ell_rad_lower.f90 ell_rad_upper.f90		\
	ellu_bbt2chol.f90 ellu_chol2eig.f90 ellu_pt_near_far.f90	\
	ellu_radii.f90 ellu_rad_upper.f90

FILESf = destsv.f dgqt.f

FILES = $(FILESf) $(FILESf90)

