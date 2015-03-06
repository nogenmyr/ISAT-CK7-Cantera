!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ell_pair_cover_query( n, c1, gg1, c2, gg2, rmax )

!  Given a pair of  ellipsoids, E1 and E2, determine whether E1
!  covers E2.  rmax <=1.d0 indicates that E1 covers E2.
!  E1 is defined by:  (x-c1)^T * G1 * G1^T * (x-c1) <= 1, 
!  where the array gg1 contains G1 in packed format.  Similarly for E2. 

!  Method:
!  1/ transform to y-space in which E1 is the unit ball at the origin: 
!        y = G1^T * (x-c)
!  2/ in y-space, determine rmax = the distance from the origin to the 
!        furthest point in E2
!  3/ if rmax <=1.d0, then E1 covers E2.


!  S.B. Pope  10/27/04

implicit none

integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: n
real(k_dp), intent(in)  :: gg1((n*(n+1))/2), gg2((n*(n+1))/2), c1(n), c2(n)
real(k_dp), intent(out) :: rmax 

integer    :: i, j, k
real(k_dp) :: g1(n,n), g2(n,n), c(n), origin(n), r(n)

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

c = matmul( (c2-c1), g1 )

!  find point in G2y furthest from the origin

origin = 0.d0
call ellu_pt_near_far( 2, n, c, g2, origin, r )

rmax = sqrt( sum(r*r) )

return
end subroutine ell_pair_cover_query
