!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ceq_reorder(nr,nc,x,order,y)
!----------------
  implicit none
  integer, intent(in) :: nr, nc, order(nr)
  real(kind(1.d0)), intent(in) :: x(nr,nc)
  real(kind(1.d0)), intent(out) :: y(nr,nc)
!----------------
! re-order rows of input matrix
! input:
!	nr		- number of rows in x
!	nc		- number of columns in x
!   x       - matrix to be re-ordered
!   order   - pointer: j-th row of reordered matrix is row order(j) of x
! output:
!   y       - reordered matrix

!  S. B. Pope 5/25/02

integer :: j

do j=1,nr
    y(j,:)=x(order(j),:)
end do

return
end subroutine ceq_reorder
