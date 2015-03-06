!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ell_eig2chol( n, u, lam, g )

!  The n x n matrix A has the eigendecomposition A = u * lam^2 * u^T
!  and the Cholesky decomposition A = G * G^T.  Given u and lam, this
!  routine returns G in packed format.

!  Method:
!     A = u * lam^2 * u^T = B * B^T, where B = u * lam
!     B = G * Q  ( LQ factorization)

implicit none
integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: n
real(k_dp), intent(in)  :: u(n,n), lam(n)
real(k_dp), intent(out) :: g((n*(n+1))/2)

real(k_dp) :: B(n,n)
integer    :: j

!  form B
do j = 1, n  !  B = u * lam
   B(:,j) = u(:,j) * lam(j)
end do

call ell_bbt2chol( n, B, g )

return
end subroutine ell_eig2chol
