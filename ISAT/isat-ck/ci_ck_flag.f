!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

! ci_ck:  subroutine ci_ckinit to call ctinit
! this version with iflag

	subroutine ci_ckinit( iflag, gas )

	use cantera
	implicit none
	integer      :: iflag
	type(phase_t)::  gas
	
	iflag = 0
	call ctinit( iflag, gas )

     	return
	end subroutine ci_ckinit

