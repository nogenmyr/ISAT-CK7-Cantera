!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ceq_ming(nz,nc,B,c,g,zg,iret)
!--------------------
  implicit none
  integer, intent(in) :: nz, nc
  real(kind(1.d0)), intent(in) :: B(nz,nc), c(nc), g(nz)
  integer, intent(out) :: iret
  real(kind(1.d0)), intent(out) :: zg(nz)
!--------------------
!  return zg, the value of z which minimizes g'z, 
!  subject to B'*z=c, z(i)>=0.
!  Exit flag from linprog: iret<=0 for failure.

!  S.B. Pope 10/1/02

integer :: i, iftest  
real(kind(1.d0)) :: BT(nc,nz), errc(nc), errmax, gmin

BT=transpose(B)
call ceq_linprog(nz,nc,g,BT,c,zg,iret)

do i=1,nz  ! guard against small negative values due to round-off
  zg(i)=max(zg(i),0.d0)
end do

iftest=0	! for testing only
if( iftest == 0  ) return
!  test satisfaction of constraints and value of g'z
errc=matmul(zg,B)-c
errmax=max( maxval(errc), -minval(errc) )
gmin=dot_product(zg,g)
write(0,*)'ceq_ming: iret, errmax, gmin = ', iret, errmax, gmin

return
end subroutine ceq_ming