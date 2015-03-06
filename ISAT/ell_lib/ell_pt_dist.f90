!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ell_pt_dist( n, c, gg, p, s )

!  Given the ellipsoid E (centered at c) and the point p, return in s
!  the relative distance from c to the boundary of E along the ray c-p
!  (relative to the distance |p-c|).
!  Let b denote the intersection of the boundary of E and the ray c-p.  
!  Then s is:     s = |b-c| / |p-c|. 
!  
!  E is given by { x | norm(G^T * (x-c) ) <=1 ),  where G is an 
!  n x n lower triangular matrix.  The array gg contains 
!  the matrix G in packed format.

!  Errors: s = -1 is returned if p = c
!          s = -2 is returned if G is not positive definite

!  S.B. Pope  1/1/06

implicit none

integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: n
real(k_dp), intent(in)  :: c(n), gg((n*(n+1))/2), p(n)
real(k_dp), intent(out) :: s

real(k_dp) :: y(n), ysq

y = p-c   !   form  y = G^T * [p-c]
call dtpmv( 'L', 'T', 'N', n, gg, y, 1 ) 

ysq = sum(y*y) 

if( ysq > 0.d0 ) then
   s = 1.d0 / sqrt( ysq )   !  s successfully evaluated
elseif( sum( (p-c)*(p-c) ) == 0.d0 ) then
   s = -1.d0                !  p=c
else
   s = -2.d0                !  G is not positive definite
endif

return
end subroutine ell_pt_dist
