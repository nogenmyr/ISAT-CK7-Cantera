!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ell_pair_cover_cv( n, c1, gg1, c2, gg2, c, gg )

!  Given a pair of ellipsoids, E1 and E2, determine a third ellipsoid
!  E which covers E1 and E2.
!  E is defined by:  (x-c)^T * G * G^T * (x-c) <= 1, where the array gg
!  contains G in packed format.  Similarly for E1 and E2. 
!  A = G * G^T, A1 = G1 * G1^T, A2 = G2 * G2^T

!  Method: c = (c1+c2)/2,  A = a * ( A1^{-1} + A2^{-1} + (c1-c2) * (c1-c2)^T /4 )^{-1},
!                          where a is as small as possible.

!  Note: this is a heuristic: E does not have minimal volume.

!  S.B. Pope  6/13/04,  5/14/07 (cv_max added)

implicit none

integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: n
real(k_dp), intent(in)  :: c1(n), gg1((n*(n+1))/2), c2(n), gg2((n*(n+1))/2)
real(k_dp), intent(out) :: c(n),  gg((n*(n+1))/2)

real(k_dp), parameter   :: cv_add = 1.d-10  !  frac. amount to add to diag. of CV to avoid num. sing.
integer    :: i, j, k, info
real(k_dp) :: g1(n,n), g2(n,n), cv1(n,n), cv2(n,n), dc(n), &
              orig(n), c3(n), xf(n), rsq1, rsq2, cv_max

c    = 0.5d0 * ( c1 + c2 )
dc   = 0.5d0 * ( c1 - c2 )
orig = 0.d0

! unpack gg1 and gg2 into g1 and g2
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

cv1 = g1
cv2 = g2

!  form covariances CV1 = (G1 * G1^T)^{-1}  and CV2
call dpotri( 'L', n, cv1, n, info )

if( info /= 0 ) then
   write(0,*)'ell_pair_cover_cv(1) dpotri failed for cv1: info = ', info
   stop
endif

call dpotri( 'L', n, cv2, n, info )

if( info /= 0 ) then
   write(0,*)'ell_pair_cover_cv(2) dpotri failed for cv2: info = ', info
   stop
endif

!  form CV = CV1 + CV2 + dc * dc^T

cv_max = 0.d0
do j = 1, n
   do i = j, n
      cv1(i,j) = cv1(i,j) + cv2(i,j) + dc(i) * dc(j)
   end do
   cv_max = max( cv1(j,j), cv_max )
end do

cv_max = cv_max * cv_add  !  amount to add to diagonal to avoid numerical singularity
do j = 1, n
   cv1(j,j) = cv1(j,j) + cv_max
end do  

!  Cholesky decomp of CV
call dpotrf( 'L', n, cv1, n, info )

if( info /= 0 ) then
   write(0,*)'ell_pair_cover_cv(3) dpotrf failed for cv: info = ', info
   stop
endif

!  Inverse of CV
call dpotri( 'L', n, cv1, n, info )

if( info /= 0 ) then
   write(0,*)'ell_pair_cover_cv(4) dpotri failed: info = ', info
   stop
endif

!  G:  CV^{-1) = G * G^T
call dpotrf( 'L', n, cv1, n, info )

if( info /= 0 ) then
   write(0,*)'ell_pair_cover_cv(5) dpotrf failed for CV^{-1}: info = ', info
   stop
endif

!  Transform so that E is the unit ball, and find the furthest point on E1
c3 = c1 - c   ! c3 = G^T (c1 - c )
call dtrmv( 'L', 'T', 'N', n, cv1, n, c3, 1 )

! g1 = G^{-1} * G1
call dtrsm( 'L', 'L', 'N', 'N', n, n, 1.d0, cv1, n, g1, n )

call ellu_pt_near_far( 2, n, c3, g1, orig, xf )  !  xf is furthest point 
rsq1 = sum( xf*xf )

!  Transform so that E is the unit ball, and find the furthest point on E2
c3 = c2 - c   ! c3 = G^T (c2 - c )
call dtrmv( 'L', 'T', 'N', n, cv1, n, c3, 1 )

! g2 = G^{-1} * G2
call dtrsm( 'L', 'L', 'N', 'N', n, n, 1.d0, cv1, n, g2, n )

call ellu_pt_near_far( 2, n, c3, g2, orig, xf )

rsq2 = sum( xf*xf )

rsq1 = max(rsq1,rsq2)
if( rsq1 > 0.d0 ) rsq1 = 1.d0 / sqrt( rsq1 )
k  = 0  !  pack re-scaled G into gg
do j = 1, n
   do i = j, n
      k = k + 1
	  gg(k) = rsq1 * cv1(i,j) 
   end do
end do

return
end subroutine ell_pair_cover_cv
