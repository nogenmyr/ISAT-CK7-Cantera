!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ell_pt_modify( n, c, gg, p )

!  Modify the ellipsoid E so that its boundary intersects the given
!  point p.  The ellipsoid E is given by: { x | norm(G^T * (x-c) ) <=1 ), 
!  where G is an n x n lower triangular matrix.  The array gg contains 
!  the matrix G in packed format.
!  Method: 
!  1/ transform to y-space in which E is the unit ball: (p-c) transforms to yp.
!  2/ generate the modified ellipsoid in y-space in terms of a rank-one 
!     modification of the identity: one principal axis is yp.
!  3/ transform back to x-space and form the modified G.

!  S.B. Pope  12/24/05

implicit none

integer, parameter        :: k_dp = kind(1.d0)
integer, intent(in)       :: n
real(k_dp), intent(in)    :: c(n), p(n)
real(k_dp), intent(inout) :: gg((n*(n+1))/2)

integer    :: i, j, k
real(k_dp) :: yp(n), ypnorm, v(n), q(n,n), g(n,n), gamma

!  Transform to y-space: y = G^T * [x-c]
yp = p-c   
call dtpmv( 'L', 'T', 'N', n, gg, yp, 1 )  !  yp = G^T * [p-c]

!  Transform to z-space:  z = Q^T * y,  where yp = Q * R
ypnorm = sqrt( sum(yp*yp) )

if( ypnorm == 0.d0 ) then  ! degenerate case p = c
   gg = 0.d0  !  set G = I * 1.d30
   k  = 0  
   do j = 1, n
      do i = j,n
         k = k + 1
	     if( i == j ) gg(k) = 1.d30
      end do
   end do
   return
endif

!  The modified ellipsoid is { y | y^T B^2 y <=1 }
!  B = I + gamma * yp * yp^2

gamma = (1.d0/ypnorm - 1.d0) / ypnorm**2

!  Transform back to x
g = 0.d0
k = 0  !  unpack g
do j = 1, n
   do i = j,n
      k = k + 1
	  g(i,j) = gg(k)
   end do
end do

! q = G*B = G + (gamma*G*yp)*yp^T

v = gamma * matmul( g, yp )

do i = 1, n
do j = 1, n
   q(i,j) = g(i,j) + v(i)*yp(j)
end do
end do


call ell_bbt2chol( n, q, gg )

return
end subroutine ell_pt_modify
