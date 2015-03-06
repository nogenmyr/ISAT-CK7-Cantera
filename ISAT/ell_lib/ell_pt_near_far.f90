!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ell_pt_near_far( knf, n, c, gg, p, xnf )

!  Determine the point xnf in the specified ellipsoid E that is
!  closest (for knf=1) or furthest (for knf=2) from the specified point p.  
!  The ellipsoid E is given by:
!  { x | norm(G^T * (x-c) ) <=1 ), where G is an n x n lower triangular
!  matrix.  The array gg contains the matrix G in packed format.

!  See ellu_pt_near_far for details

!  S.B. Pope  12/1/2006

implicit none

integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: knf, n
real(k_dp), intent(in)  :: c(n), gg((n*(n+1))/2), p(n)
real(k_dp), intent(out) :: xnf(n)

integer    :: i, j, k
real(k_dp) :: g(n,n)

!  unpack gg into g
k = 0
g = 0.d0
do j = 1, n
   do i = j, n
      k = k + 1
	  g(i,j) = gg(k)
   end do
end do

call ellu_pt_near_far( knf, n, c, g, p, xnf )

return
end subroutine ell_pt_near_far
