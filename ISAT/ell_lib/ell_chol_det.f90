!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ell_chol_det( n, g, det )

!  g contains (in packed format) the lower triangular n x n matrix G.
!  The determinant of G, det, is returned.
!  If G represents an ellipsoid E, then det is the inverse of the 
!  content (volume) of E.

implicit none
integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: n
real(k_dp), intent(in)  :: g((n*(n+1))/2)
real(k_dp), intent(out) :: det

integer    :: j, k

det = 1.d0

k  = 1
do j = 1, n
   det = det * g(k)
   k = k + n + 1 - j
end do

return
end subroutine ell_chol_det
