!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ell_pt_hyper( n, c, gg, p, xh, v )

!  Given an ellipsoid E (centered at c) and a point p, determine the
!  hyperplane H which is the perpendicular bisector of the line c-p 
!  in the transformed space in which E is the unit ball.

!  E is defined by:  (x-c)^T * G * G^T * (x-c) <= 1.  The array
!  gg contains the lower triangular n x n matrix G in packed format.

!  The hyperplane H is defined by v^T * ( x - xh  ) = 0, where v is
!  a unit vector.  The quantity s(x) = v^T * ( x - xh  ) is the signed
!  distance of the point x from the hyperplane: s(c) is negative, and 
!  s(p) is positive. 

!  S.B. Pope  6/4/2006
implicit none

integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: n
real(k_dp), intent(in)  :: c(n), gg((n*(n+1))/2), p(n)
real(k_dp), intent(out) :: xh(n), v(n)

real(k_dp) :: vsq

xh = 0.5d0 * ( c + p )
v  = p-c   
call dtpmv( 'L', 'T', 'N', n, gg, v, 1 )  !   = G * G^T * [p-c]
call dtpmv( 'L', 'N', 'N', n, gg, v, 1 ) 

vsq = sum( v * v )
if( vsq > 0.d0 ) then
   v = v / sqrt( vsq ) 
   return
endif

write(0,*)'ell_pt_hyper: vsq = 0.d0'
stop 

end subroutine ell_pt_hyper
