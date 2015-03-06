!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ceq_lss(nb,nx,A,b,x,info)

!-------------
  implicit none
  integer,          intent(in)  :: nb, nx
  real(kind(1.d0)), intent(in)  :: A(nb,nx), b(nb)
  real(kind(1.d0)), intent(out) :: x(nx)
  integer,          intent(out) :: info
!--------------

!  Determine the least-squares/minimum-norm solution x to
!  the linear equation A x = b.
!  
!  Input:
!	nb	- number of rows in b
!	nx	- number of rows in A and x
!	A	- the nb x nx matrix A
!	b	- the nb-vector b
!  Output:
!	x	- the solution nx-vector
!	info=0 for successful solution

!	S.B. Pope 10/2/02

integer :: lwork, rank
real(kind(1.d0)) :: tol=1.d-9, aa(nb,nx), bb(nb+nx), sv(nb+nx), &
  work( 4*(nb+nx+1)*(nb+nx+1) )

lwork= size(work)
aa=A
bb=0.d0
bb(1:nb)=b

call dgelss(nb,nx,1,aa(1:nb,1:nx),nb,bb(1:nb+nx),nb+nx,sv(1:nb+nx), &
            tol,rank,work(1:lwork),lwork,info)
x=bb(1:nx)

end subroutine ceq_lss