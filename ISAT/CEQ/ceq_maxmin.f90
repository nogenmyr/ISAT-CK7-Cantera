!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ceq_maxmin(nz,nc,B,c,zm,zmin,iret)
!----------------------
implicit none
integer, intent(in) :: nz, nc
real(kind(1.d0)), intent(in) :: B(nz,nc), c(nc)
integer, intent(out) :: iret
real(kind(1.d0)), intent(out) :: zm(nz), zmin
!-----------------------

!  Determine the max-min composition.
!  Find zm which maximizes zmin=min_i(z(i)) subject to B'*z=c.

!  iret<0 indicates failure

!  Method: 
!   Initially assume zmin>=0.
!   Define: x=[ (z-zmin)' zmin]'
!   Maximize zmin (i.e., minimize -x(n)) subject to x(i)>=0
!       and B'*z=c, which is equivalent to A*x=c,
!       where A=[ B' -sum(B)'].
!   If feasible solution not found, zmin<0, and
!   Define: x=[ (z-zmin)' -zmin]'
!   Maximize zmin (i.e., minimize x(n)) subject to x(i)>=0
!       and B'*z=c, which is equivalent to A*x=c,
!       where A=[ B' sum(B)'].

!   S.B. Pope 10/1/02

integer :: nx, tries,  j, jj
real(kind(1.d0)) :: A(nc,nz+1), bsum(nc), f(nz+1), x(nz+1)


nx=nz+1
bsum=sum(B,dim=1)
f=0
!  first assume zmin>0, x=[z'-zmin zmin]'
f(nx)=-1.d0  ! minimize -zmin
A(1:nc,1:nz)=transpose(B)
A(1:nc,nx)=bsum

call ceq_linprog(nx,nc,f,A,c,x,iret)

if( iret == 0 ) then  !  success, zmin>=0
  zmin=x(nx)
  zm=x(1:nz)+zmin
  return
elseif( iret /= -2 ) then  !  failure
  return
endif
!  zmin<0, re-define x=[z'-zmin -zmin]'
f(nx)=1.d0  ! minimize -zmin
A(1:nc,nx)=-bsum

call ceq_linprog(nx,nc,f,A,c,x,iret)

if( iret /= 0 ) return  ! failure

zmin=-x(nx)
zm=x(1:nz)+zmin

return
end subroutine ceq_maxmin