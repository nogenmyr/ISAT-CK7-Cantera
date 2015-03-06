!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module isat_defaults

   use isat_types, only: l_info, l_rinfo, k_xf
   implicit none

!  set default values (and limits) on info and rinfo

   integer,    save, dimension(l_info)  :: info_def, info_min, info_max
   real(k_xf), save, dimension(l_rinfo) :: rinfo_def, rinfo_min
           
contains

!====================================================================== 

subroutine isat_info_defaults( def_set )

   implicit none
   integer, intent(in) :: def_set  !  default settings to be used

!----------  set defaults for info -----------------------------
   
   info_def     = 0
   info_min     = 0  ! minimum allowed value
   info_max     = 1  ! maximun allowed value: 
!                      set info_max=0 for non-negative; = -1 for unbounded
   info_def( 4) = 5
   info_def( 5) = 30
   info_def( 6) = 10
   info_def( 7) = -1
   info_def( 8) =  1
   info_def( 9) = -1

   info_def(19) = 100

   info_def(35) = 1
   info_def(36) = 4
   info_def(38) = 3
   
   info_def(39) = -1
   info_def(40) = -1
   info_def(41) = -1
   info_def(42) = -1
   info_def(43) = 1

   info_def(44) = 1000
   info_def(45) = 100
   
   info_def(49) = -3

   info_def(50) = 2
   info_def(51) = 2
   info_def(52) = 2

   info_min( 4)   = -3
   info_min( 5)   = -3
   info_min( 6:9) = -1
   info_min(13)   = -1

   info_min(35)   = -1
   info_min(36)   = -1
   info_min(38)   = -1
   
   info_min(39:42) = -1
   info_min(43)   = -2
   
   info_min(49)   = -3

   info_min(53:55)  = -1

   info_max(4:9) = -1
   info_max(11)  =  2
   info_max(12)  =  2
   info_max(13)  =  -1
   info_max(15)  =  3
   info_max(19)  =  0
   info_max(28)  =  0

   info_max(35)  =  3
   info_max(36)  =  4
   info_max(37)  =  3
   info_max(38)  =  3

   info_max(39:42)  =  -1
   info_max(43)  = 2

   info_max(44)  =  0
   info_max(45)  =  -1

   info_max(48)  = -1
   info_max(49)  = -1

   info_max(50)  =  8
   info_max(51)  =  8
   info_max(52)  =  5
   info_max(53:55) = -1
   info_max(60)  =  2
   info_max(61)  =  0
   info_max(62)  =  0
   info_max(66)  =  2
   info_max(70)  =  2

!--------- set defaults for rinfo  ------------------------------

    rinfo_def = 0.d0
    rinfo_min = 0.d0
    rinfo_def( 1) = 1.d-4
	rinfo_def( 3) = 1.d20
	rinfo_def( 4) = 1.d-1
	rinfo_def( 5) = 10.d0
	rinfo_def( 6) = 0.5d0
	rinfo_def( 7) = 2.0d0
	rinfo_def( 8) = 500.d0
	rinfo_def( 9) = 1.02d0
	rinfo_def(10) = 1.2d0
	rinfo_def(11) = 1.2d0

	rinfo_def(12) = 10.d0
	rinfo_def(13) = 1.2d0
	rinfo_def(14) = 1.d-1
	rinfo_def(15) = 2.d0
	rinfo_def(16) = 0.9d0
	rinfo_def(17) = 0.9d0
	rinfo_def(18) = 1.0d0
	rinfo_def(19) = 0.5d0
	rinfo_def(20) = 0.5d0
	rinfo_def(21) = 0.1d0

	rinfo_min(6)  = -huge(1.d0)
	rinfo_min(7)  = -huge(1.d0)
	rinfo_min(8)  = -huge(1)
	rinfo_min( 9) = 1.d0
	rinfo_min(10) = 1.d0
	rinfo_min(11) = 1.d0
	rinfo_min(12) = -huge(1.d0)
	rinfo_min(13) = 1.d0
	rinfo_min(15) = 1.d0
	rinfo_min(16) = -huge(1.d0)
	rinfo_min(21) = -huge(1.d0)
	
!------- special default over-rides

    if( def_set == 1 ) then  !  easy case: set for speed
       info_def(40)  = 10       !  degen_g
       rinfo_def(6)  = 0.2d0    ! ret_frac
       rinfo_def(7)  = 0.5d0    ! grow_frac
       rinfo_def(18) = 1.d0     ! slow_prog_thresh
       rinfo_def(19) = rinfo_def(6) ! ret_frac_slow
       rinfo_def(20) = rinfo_def(7) ! grow_frac_slow
    
    elseif( def_set == 2 ) then  !  hard case/long run: set to maximize use of memory
       info_def(40)  = 100      !  degen_g
       rinfo_def(6)  = 1.d0     ! ret_frac
       rinfo_def(7)  = 10.d0    ! grow_frac
       rinfo_def(18) = 2.d0     ! slow_prog_thresh
       rinfo_def(19) = rinfo_def(6) ! ret_frac_slow
       rinfo_def(20) = 1.d0     ! grow_frac_slow    
    elseif( def_set == 3 ) then  !  hard case/short run: limit time on growing
       info_def(37)  = 3        ! mode_eoi - no shrinking
       info_def(40)  = 10       ! degen_g
       rinfo_def(6)  = 1.d0     ! ret_frac
       rinfo_def(7)  = 1.d0     ! grow_frac
       rinfo_def(18) = 1.d0     ! slow_prog_thresh
       rinfo_def(19) = rinfo_def(6) ! ret_frac_slow
       rinfo_def(20) = rinfo_def(7) ! grow_frac_slow    
    endif
    
	return

end subroutine isat_info_defaults

end module isat_defaults
