!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ell_bbt2eig( n, B,  U, lam )

!  The n x n PSD matrix A is given by A = B * B^T.
!  The SVD of B is:  B = U * S * V^T, so that the
!  eigendecomposition of A is: A = U * S^2 * U^T.
!  Given B, return U and lam=diag(S).

implicit none
integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: n
real(k_dp), intent(in)  :: B(n,n)
real(k_dp), intent(out) :: U(n,n), lam(n)


real(k_dp) :: vt(n,n), work(5*n*n+20*n)
integer    :: lwork, info

lwork = 5*n*n+20*n

vt = B
call dgesvd( 'A', 'N', n, n, vt, n, lam, u, n, vt, n, work, lwork, info)

if( info /= 0 ) then
   write(0,*)'ell_bbt2eig: info, lwork, work(1) = ', info, lwork, work(1)
   stop
endif

return
end subroutine ell_bbt2eig

