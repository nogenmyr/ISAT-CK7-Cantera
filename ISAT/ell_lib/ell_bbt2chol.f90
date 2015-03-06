!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ell_bbt2chol( n, B, g )

!  Given an n x n matrix B, compute the lower  
!  Cholesky triangle  G:   G * G^T = B * B^T

!  The Cholesky triangle G is returned in g in packed format.

!  Method: B = L * Q;  so that B * B^T = L * L^T; and so G = L.

implicit none
integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: n
real(k_dp), intent(in)  :: B(n,n)
real(k_dp), intent(out) :: g((n*(n+1))/2)

integer    :: i, j, k, info, lwork
real(k_dp) :: gf(n,n), tau(n), work(10+n*(10+n))

!  perform LQ factorization: B = G*Q
lwork = 10+n*(10+n)
gf = B
call dgelqf( n, n, gf, n, tau, work, lwork, info )

if( info /=0 ) then
   write(0,*)'ell_bbt2chol: dgelqf failed, info = ', info
   stop
endif 

!  put in packed format, ensuring that diagonal elements are non-negative
k = 0
do j = 1, n
   if( gf(j,j) >= 0.d0 ) then
      do i = j,n
         k = k + 1
	     g(k) = gf(i,j)
      end do
   else
      do i = j,n
         k = k + 1
	     g(k) = -gf(i,j)
      end do
   endif
end do

return
end subroutine ell_bbt2chol
