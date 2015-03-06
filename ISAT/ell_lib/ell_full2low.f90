!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ell_full2low( n, A, L )

!  Given the n x n PSD matrix L, return in L the lower
!  triangle of A in packed format.

implicit none
integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: n
real(k_dp), intent(in)  :: A(n,n)
real(k_dp), intent(out) :: L((n*(n+1))/2)

integer    :: i, j, k

k = 0
do j = 1, n
   do i = j,n
      k = k + 1
	  L(k) = A(i,j)
   end do
end do

return
end subroutine ell_full2low
