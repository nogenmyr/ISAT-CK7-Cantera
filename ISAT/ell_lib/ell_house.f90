!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ell_house( n, x, v, beta )

!  Determine a Householder vector v based on x.
!  The orthogonal matrix P = I - beta * v * v^T has the
!  property: Px = |x| * e_1.

!  From Golub & Van Loan Algorithm 5.1.1

integer, intent(in)           :: n
real(kind(1.d0)), intent(in)  :: x(n)
real(kind(1.d0)), intent(out) :: v(n), beta

real(kind(1.d0)) :: sigma, mu

v    = 0.d0
beta = 0.d0
if( n <= 1 ) return

sigma  = sum( x(2:n)*x(2:n) )
v(1)   = 1.d0
v(2:n) = x(2:n)

if( sigma == 0.d0 ) then
   beta = 0.d0
else
   mu = sqrt( x(1)**2 + sigma )
   if( x(1) < 0.d0 ) then
      v(1) = x(1) - mu
   else
      v(1) = -sigma / ( x(1) + mu )
   endif

   beta = 2.*v(1)**2 / (sigma+v(1)**2)
   v    = v / v(1)
endif

return

end subroutine ell_house
