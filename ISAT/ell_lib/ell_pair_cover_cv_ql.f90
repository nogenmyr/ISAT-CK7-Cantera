!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ell_pair_cover_cv_ql( n, c1, gg1, c2, gg2, c, gg )

!  Given a pair of ellipsoids, E1 and E2, determine a third ellipsoid
!  E which covers E1 and E2.
!  E is defined by:  (x-c)^T * G * G^T * (x-c) <= 1, where the array gg
!  contains G in packed format.  Similarly for E1 and E2. 

!  Method:
!  c = (c1+c2)/2
!  With A = G * G^T, A1 = G1 * G1^T, A2 = G2 * G2^T, d = (c1-c2)/2, A is obtained from
!  A^{-1} = a * ( A1^{-1} + A2^{-1} + d * d^T )
!         = a * B^T * B,  where B = [G_1^{-1}; G_2^{-1}; d^T]
!         = a * L^T * L,  where B = Q * [0; L] is the QL factorization of B
!         = G^{-T} * G,   where a is determined so that E just covers E1 and E2.
!  Thus, G = L^{-1} * sqrt(a)

!  Note: this is a heuristic: E does not have minimal volume.

!  S.B. Pope  12/2/2007

implicit none

integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: n
real(k_dp), intent(in)  :: c1(n), gg1((n*(n+1))/2), c2(n), gg2((n*(n+1))/2)
real(k_dp), intent(out) :: c(n),  gg((n*(n+1))/2)

integer    :: i, j, k, info, lwork, m
real(k_dp) :: g1(n,n), g2(n,n), cv1(n,n), cv2(n,n), dc(n), &
              orig(n), c3(n), xf(n), rsq1, rsq2, &
              b(n+n+1,n), work(3*n*n), eye(n,n)
              
lwork = 3*n*n
m     = n+n+1

eye   = 0.d0 !  set identity
do i = 1, n
   eye(i,i) = 1.d0
end do

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

!  form inverses G1^{-1}, G2^{-1}
cv1 = eye
cv2 = eye
call dtrsm( 'L', 'L', 'N', 'N', n, n, 1.d0, g1, n, cv1, n )
call dtrsm( 'L', 'L', 'N', 'N', n, n, 1.d0, g2, n, cv2, n )

!  assemble B=[ G1^{-1}; G2^{-1}; d^T ]
b(1:n,1:n)     = cv1
b(n+1:n+n,1:n) = cv2
b(m,1:n)       = dc

!  B = Q [0; L]
call dgeqlf( m, n, b, m, c3, work, lwork, info )

if( info /= 0 ) then
   write(0,*)'ell_pair_cover_cv_ql dgeqlf failed: info = ', info
   stop
endif

!  extract L from the QL factorization: cv2 = L
cv2 = 0.d0
do i = 1, n
   do j = 1, i
      cv2(i,j) = b(i+n+1,j)
   end do
end do

!  change sign of rows as needed so that diagonals are positive
do j = 1, n
   if( cv2(j,j) < 0.d0 ) cv2(j,1:j) = -cv2(j,1:j)
end do

!  cv1=L_0 = L^{-1} = G / sqrt(a)
cv1 = eye
call dtrsm( 'L', 'L', 'N', 'N', n, n, 1.d0, cv2, n, cv1, n )

!  Determine a such that E just covers E1 and E2

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
end subroutine ell_pair_cover_cv_ql
