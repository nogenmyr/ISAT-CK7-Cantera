!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module isat_lu_m

   implicit none
   integer, save :: lu = 60, luf = 60, lul = 99 

contains

   subroutine isat_lu_set( lufset, lulset )

   integer, intent(in) :: lufset, lulset

   luf = lufset
   lul = lulset

   return
   end subroutine isat_lu_set

end module isat_lu_m
