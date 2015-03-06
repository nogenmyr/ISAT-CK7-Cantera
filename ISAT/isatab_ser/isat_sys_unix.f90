!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

function isat_sys( command )

!  call system to execute unix command: return 0 if successful

   implicit none
   character(128), intent(in) :: command
   integer :: isat_sys, system

   isat_sys = system( command )

   return
end function isat_sys

subroutine isat_systype( systype )

!  return system type

   character(6), intent(out) :: systype
   character(6)              :: this_systype = 'UNIX  '

   systype = this_systype

   return
end subroutine isat_systype
