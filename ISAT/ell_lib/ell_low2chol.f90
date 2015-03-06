!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ell_low2chol( n, l, g )

!  l contains the lower triangle L of the n x n PSD matrix A,
!  in packed format.
!  The Cholesky lower triangle G,  A = G * G^T is returned in g
!  in packed format in g.

implicit none
integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: n
real(k_dp), intent(in)  :: l((n*(n+1))/2)
real(k_dp), intent(out) :: g((n*(n+1))/2)

integer :: info

g = l
call dpptrf( 'L', n, g, info ) 

if( info /= 0 ) then
   write(0,*)'ell_low2chol: info = ', info
   stop
endif

return
end subroutine ell_low2chol
