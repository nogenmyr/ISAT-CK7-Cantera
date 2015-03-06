!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ell_line_proj( n, c, gg, x0, v, smin, smax )

!  Given the ellipsoid E = { x | norm(G^T * (x-c) ) <=1 ), 
!  where G is an n x n lower triangular, and given the line
!  L = { x | x = x0 + v * s }, determine the line segment
!  (smin, smax) which is the orthogonal projection of E onto L.
!  The array gg contains the matrix G in packed format.

!  Method:  E = { x | x = c + G^{-T} u,  |u| <=1}
!           s = v^T * ( x - x0 ) / v^T * v
!             = v^T * ( c - x0 + G^{-T} u ) / v^T * v
!             = s0 + w^T * u,  s0 = v^T * ( c - x0 ) / v^T * v,
!                              w  = G^{-1} * v  / v^T * v
!           smax = s0 + |w|, smin = s0 - |w|

!  S.B. Pope  6/12/04

implicit none

integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: n
real(k_dp), intent(in)  :: c(n), gg((n*(n+1))/2), x0(n), v(n)
real(k_dp), intent(out) :: smin, smax

real(k_dp) :: vnormsq, w(n), s0, wnorm

vnormsq = sum(v*v)   !  v^T * v = |v|^2
if( vnormsq == 0.d0 ) then
   write(0,*)'ell_line_proj: |v| = 0'
   stop
endif

w  = v / vnormsq
s0 = sum( w * (c-x0) )  !  s0 = v^T ( c - x0 ) / |v|^2

call dtpsv( 'L', 'N', 'N', n, gg, w, 1 ) ! w  = G^{-1} * v  / v^T * v
wnorm = sqrt( sum( w * w ) )

smin = s0 - wnorm
smax = s0 + wnorm

return
end subroutine ell_line_proj