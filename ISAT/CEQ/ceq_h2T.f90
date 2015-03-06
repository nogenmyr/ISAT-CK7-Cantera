!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ceq_h2T(ns,z,hin,T_low,T_high,thermo,T,iret)
!-------------------
  implicit none
  integer, parameter :: nth=15
  integer, intent(in) :: ns
  real(kind(1.d0)), intent(in)  :: z(ns), hin, T_low, T_high, thermo(ns,nth)
  real(kind(1.d0)), intent(out) :: T
  integer,          intent(out) :: iret
!-------------------
!  determine temperature given enthalpy
! input:
!	ns		- number of species
!   z       - moles of species
!   hin     - enthalpy/R (K) = z'*h
!   T_low   - lower bound on temperature range
!   T_high  _ upper bound on temperature range
!   thermo  - thermo data
! output:
!   T    - temperature (K)
!   iret =  0  for success
!        = -1  for T > T_high
!        = -2  for T < T_low
!        = -3  for failure

! Notes:  if the temperature is outside the range [T_low T_high]
!   then T is returned as the closest of these bounds.
!   If iteration fails, T is returned as T=-1.

! S. B. Pope 9/26/02

integer :: itmax,it
real(kind(1.d0)) :: T_tol, T0, hort(ns), hor, h_a, T_a, h_b, T_b, dT, &
	cpor(ns), hh, cpp

itmax=100      ! maximum number of Newton iterations (usually only 3 required)
T_tol=1.d-6    ! error tolerance
T0=1500.d0     ! initial guess
iret=0         ! anticipate success

!  determine if T>T0 and bracket T in [T_a T_b]
call ceq_h(ns,T0,thermo,hort)
hor=dot_product(z, hort)*T0

if( hin>hor ) then	! T > T0=T_a
    h_a=hor
    T_a=T0
    call ceq_h(ns,T_high,thermo,hort)
    h_b=dot_product(z, hort)*T_high
    if( hin>=h_b ) then
        T=T_high   ! T > T_high (return T=T_high)
		iret = -1
        return
    endif
    T_b=T_high
else
    h_b=hor	! T < T0=T_b
    T_b=T0
    call ceq_h(ns,T_low,thermo,hort)
    h_a=dot_product(z, hort)*T_low
    if( hin<=h_a ) then
        T=T_low    ! T < T_low (return T=T_low)
		iret = -2
        return
    endif
    T_a=T_low
endif

!  estimate of T based on linear interpolation
T=T_a+(hin-h_a)*(T_b-T_a)/(h_b-h_a)
    
!  Newton iterations
do it=1,itmax
    call ceq_cp(ns,T,thermo,cpor)
    call ceq_h(ns,T,thermo,hort)
    hh=dot_product(z, hort)*T
    cpp=dot_product(z, cpor)
    dT=(hin-hh)/cpp
    T=T+dT
    if( abs(dT) < T_tol ) return  ! success
end do

! failure
iret = -3
T=-1.d0

end subroutine ceq_h2T