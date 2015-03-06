!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ceq_linprog(nx,nb,f,A,b,xm,iret)

!--------------------

  use linprog

  implicit none

  integer, intent(in) :: nx, nb

  real(kind(1.d0)), intent(in) :: f(nx), A(nb,nx), b(nb)

  integer, intent(out) :: iret

  real(kind(1.d0)), intent(out) :: xm(nx)

!--------------------

!  Determine x=xm which minimizes g=f'*x subject to

!  x(i)>=0 and A*x=b, where A has full rank.



! Input:

!	nx	- number of components of x

!	nb	- number of components of b

!	f	- nx-vector f

!	A	- nx x nb matrix A

! Output:

!	xm	- solution	

!	iret= 0 if solution is found 

!	iret=-1 if g is unbounded

!	iret=-2 if there is no feasible solution

!	iret=-3 if A is rank deficient



!  S.B. Pope 10/1/06



real(kind(1.d0)) :: eps=1.d-9, ale(1,1), age(1,1), ble(1), bge(1)



call lp( nx, 0, 0, nb, ale, age, A, ble, bge, b, f, xm, iret, toler=eps )
if( iret<0 ) then
   !XXX write(0,*)'ceq_linprog, iret = ', iret  ! SBP XXX
endif

          

end subroutine ceq_linprog