!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ceq_input(ns,ne,nth,nsamp,ndi,Ein,thermo,DI)

!  Read in data
!----------------
implicit none
integer, intent(in) :: ns,ne,nth,nsamp,ndi
real(kind(1.d0)), intent(out) :: Ein(ns,ne), thermo(ns,nth), DI(nsamp,ndi)
!-----------------
integer :: i

open( 10, file='element.txt')
do i=1,ns
  read(10,*,err=100,end=101)Ein(i,:)
end do
close(10)

open( 10, file='ck_thermo.txt')
do i=1,ns
  read(10,*,err=102,end=103)thermo(i,:)
end do
close(10)

open( 10, file='pasrp.op')
do i=1,nsamp
  read(10,*,err=104,end=105)DI(i,:)
end do
close(10)

return

100		write(0,*)'Error reading element.txt'
return
102		write(0,*)'Error reading ck_thermo.txt'
return
104		write(0,*)'Error reading pasrp.op'
return

101		write(0,*)'Hit end reading element.txt'
return
103		write(0,*)'Hit end reading ck_thermo.txt'
return
105		write(0,*)'Hit end reading pasrp.op'
return

end subroutine ceq_input