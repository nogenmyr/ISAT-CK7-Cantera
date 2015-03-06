!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ell_pair_shrink( n, gg1, gg2, n_shrink )

!  Given a pair of concentric ellipsoids, E1 and E2, shrink E2
!  (if necessary) so that it is contained in E1.
!  E1 is defined by:  x^T * G1 * G1^T * x <= 1, where the array gg1
!  contains G1 in packed format.  Similarly for E2. 

!  Method:
!  1/ transform to y-space in which E1 is the unit ball: y = G1^T * x
!  2/ in y-space, SVD of G2 is G2 = U S V^T
!  3/ modify S by S(i,i) = max( S(i,i), 1. ); number of modifications is n_shrink
!  4/ transform back

!  S.B. Pope  6/13/04

implicit none

integer, parameter        :: k_dp = kind(1.d0)
integer, intent(in)       :: n
real(k_dp), intent(in)    :: gg1((n*(n+1))/2)
real(k_dp), intent(inout) :: gg2((n*(n+1))/2)
integer, intent(out)      :: n_shrink

integer    :: i, j, k, info, lwork
real(k_dp) :: g1(n,n), g2(n,n), s(n), u(n,n), y(n), work(n*n+10*n), r, r_in

lwork    = n*n+10*n+10
n_shrink = 0

! unpack gg1 and gg2
k  = 0
g1 = 0.d0
g2 = 0.d0
do j = 1, n
   do i = j, n
      k = k + 1
	  g1(i,j) = gg1(k)
	  g2(i,j) = gg2(k)	  
   end do
end do

! transform to y-space:  y = G1^T * x;  G2y = G1^{-1} * G2
call dtrsm( 'L', 'L', 'N', 'N', n, n, 1.d0, g1, n, g2, n )

!  find point in G2y furthest from the origin

call ellu_radii( n, g2, r_in, r )

if( r <= 1.d0 ) return  !  no modification required (G2 contained in G1)

! SVD:  G2y = U S V^T    (Note: vt is not referenceed)
call dgesvd( 'A', 'N', n, n, g2, n, s, u, n, y, n, work, lwork, info )

if( info /= 0 ) then
   write(0,*)'ell_pair_shrink: SVD failed, info = ', info
   stop
endif

!  shrink principal axes as needed
do i = 1, n
   if( s(i) < 1.d0 ) then
      s(i) = 1.d0
	  n_shrink = n_shrink + 1
   endif
end do

!  form modified G2y * G2y^T = (U * S) * (U * S)^T 
do j = 1, n
   u(:,j) = u(:,j) * s(j)  ! U * S
end do

! Form Chholesky G:  G * G^T = G2y * G2y^T = (U * S) * (U * S)^T
call ellu_bbt2chol( n, u, g2 )

!  Zero upper triangle
do j = 2, n
   do i = 1, j-1
      g2(i,j) = 0.d0
   end do
end do

!  transform back to x-space:  G2 = G1 * G2y
call dtrmm( 'L', 'L', 'N', 'N', n, n, 1.d0, g1, n, g2, n )

k = 0
do j = 1, n   !  pack g2 into gg2  ensuring that the diagonal elements are non-nehative
   if( g2(j,j) >= 0.d0 ) then
      do i = j,n
         k = k + 1
	     gg2(k) = g2(i,j)
      end do
   else
      do i = j,n
         k = k + 1
	     gg2(k) = -g2(i,j)
      end do
   endif
end do

return
end subroutine ell_pair_shrink
