!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

real(kind(1.d0)) function ceq_norm( n, x )

! return the 2-norm of the double-precision n-vector x

   integer, intent(in)          :: n
   real(kind(1.d0)), intent(in) :: x(n)

   if( n <= 0 ) then
      ceq_norm = 0.d0
	  return
   endif

   ceq_norm = sqrt( sum( x(1:n) * x(1:n) ) )

   return
end function ceq_norm 
