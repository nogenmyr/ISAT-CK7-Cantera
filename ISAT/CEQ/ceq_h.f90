!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ceq_h(ns,T,thermo,hort)
!----------------
  implicit none
  integer, parameter :: ncof=7
  integer, intent(in) :: ns
  real(kind(1.d0)), intent(in) :: T, thermo(ns,2*ncof+1)
  real(kind(1.d0)), intent(out) :: hort(ns)
!-------------
!  return normalized enthalpies at temperature T
! input:
!	ns	   - number of species
!   T      - temperature (K)
!   thermo - thermo data for all species
! output:
!   hort   - h_j/(RT)  - normalized enthalpies

! S. B. Pope 9/26/02

real(kind(1.d0)) :: th(6), Tpnm1
integer :: k, n

th(1)=1.d0  ! coefficient multipliers for enthalpy
th(6)=1./T
Tpnm1=1.d0
do n=2,5
    Tpnm1=Tpnm1*T     ! =T.^(n-1)
    th(n)=Tpnm1/float(n)     ! =T.^(n-1) ./ n
end do

do k=1,ns
    if( T<thermo(k,1) ) then
        hort(k)=dot_product( thermo(k,2:7), th )  ! coefficients in lower temperature range
    else
        hort(k)=dot_product( thermo(k,9:14), th ) ! coefficients in upper temperature range
    endif
end do

end subroutine ceq_h