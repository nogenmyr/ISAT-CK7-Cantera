!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ellu_bbt2chol( n, B, G )

!  Given an n x n matrix B, compute the lower  
!  Cholesky triangle  G:   G * G^T = B * B^T

!  The Cholesky triangle G is returned in full format.

!  Method: B = L * Q;  so that B * B^T = L * L^T; and so G = L.
!  Note: the upper triangle of G is NOT set to zero.

implicit none
integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: n
real(k_dp), intent(in)  :: B(n,n)
real(k_dp), intent(out) :: G(n,n)

integer    :: j, info, lwork
real(k_dp) :: tau(n), work(10+n*(10+n))

!  perform LQ factorization: B = G*Q
lwork = 10+n*(10+n)
G = B
call dgelqf( n, n, G, n, tau, work, lwork, info )

if( info /=0 ) then
   write(0,*)'ell_bbt2cholf: dgelqf failed, info = ', info
   stop
endif 

!  ensure that diagonal elements are non-negative
do j = 1, n
   if( G(j,j) < 0.d0 ) G(j:n,j) = -G(j:n,j)
end do

return
end subroutine ellu_bbt2chol
