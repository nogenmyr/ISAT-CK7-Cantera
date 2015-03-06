!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ell_aff_pr( n, c, gg, m, d, T, cp, ggp )

!  Given an ellipsoid E (in R^n), and an m-dimensional affine space A, 
!  determine the ellipsoid P(E) which is the orthogonal projection of E onto A,
!  for  1 <= m <= n.

!  E is defined by:  || G^T * (x-c) || <= 1, where the array gg
!  contains the lower Cholesky triangle G in packed format. 

!  A is defined by { x | x = d + T * s }, where T is an n x m orthogonal
!  matrix, and s is any m-vector.

!  The orthogonal projection of the general point x is given by:
!  P(x) = { x | d + T * s,  s = T^T * ( x - d ) }.

!  The orthogonal projection of E onto A is given by:
!  P(E) = { x | d + T * s,  || GP^T * ( s - cp ) || <= 1 },  where the
!  array  ggp contains GP in packed format.  

!  Method:
!    cp = T^T * ( c - d )
!    F  = G^{-1} * T = U * [ S ; 0 ] * V^T
!    B  = V * S^{-1}
!    GP * GP^T = B * B^T

implicit none

integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: n, m
real(k_dp), intent(in)  :: c(n), gg((n*(n+1))/2), d(n), T(n,m)
real(k_dp), intent(out) :: cp(m),  ggp((m*(m+1))/2)

integer, parameter ::  ludiag = 0
integer    :: i, j, k, info, lwork
real(k_dp) :: G(n,n), F(n,m), S(m), B(m,m), VT(m,m), work(5*n*n+20*n)

!  check dimensions
if( n < 1 ) then
   write(ludiag,*)'ell_aff_pr: n < 1, n = ', n
   stop
endif

lwork = 5*n*n+20*n  ! amount of work space for SVD

if( m > n ) then
   write(ludiag,*)'ell_aff_pr: m > n, [m n] = ', m, n
   stop
elseif( m < 1 ) then
   write(ludiag,*)'ell_aff_pr: m < 1, m = ', m
   stop
endif

!  form cp = T^T * ( c - d )

cp = matmul( (c-d) , T )

!  unpack G

G = 0.d0

k  = 0
do j = 1, n
   do i = j, n
      k = k + 1
	  G(i,j) = gg(k)
   end do
end do

!  form F = G^{-1} * T
F = T
call dtrsm('L','L','N','N',n,m,1.d0,G,n,F,n)

!  form S and V^T in:  F = U * [ S ; 0 ] * V^T

call dgesvd('N','A',n,m,F,n,S,B,m,VT,m,work,lwork,info)

if( info /= 0 ) then
   write(ludiag,*)'ell_aff_pr: SVD failed, info = ', info
   stop
elseif( S(m) <= 0.d0 ) then
   write(ludiag,*)'ell_aff_pr: zero singular value, S = ', S
   stop
endif

!  form B = V * S^{-1}

do j = 1, m
   B(1:m,j) = VT(j,1:m) / S(j)
end do

!  form GP:  GP * GP^T = B * B^T

call ell_bbt2chol( m, B, ggp )

return
end subroutine ell_aff_pr 