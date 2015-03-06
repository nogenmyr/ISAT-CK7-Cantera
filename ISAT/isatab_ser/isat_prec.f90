!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module isat_prec	!  define precision to be used in isat

   implicit none

!  the following integer parameters are defined:



!  k_xf		

   integer, parameter :: k_sp   = kind(1.0)   !  kind for single precision
   integer, parameter :: k_dp   = kind(1.0d0) !  kind for double precision
   integer, parameter :: k_xf   = k_dp        !  kind for x, f, etc. in isatab
   integer, parameter :: k_g    = k_sp        !  kind for g in isatab

end module isat_prec
