!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ceq_cp(ns,T,thermo,cpor)
!-------------
  implicit none
  integer, parameter :: ncof=7
  integer, intent(in) :: ns
  real(kind(1.d0)), intent(in) :: T, thermo(ns,2*ncof+1)
  real(kind(1.d0)), intent(out) :: cpor(ns)
!-------------

!  return normalized Cp's at temperature T
! input:
!	ns	   - number of species
!   T      - temperature (K)
!   thermo - thermo data for all species
! output:
!   cpor   - Cp_j/R    - normalized constant-pressure specific heats

! S. B. Pope 9/26/02

real(kind(1.d0)) :: tc(5)
integer :: k, n

cpor=0.d0


tc(1)=1.d0  ! coefficient multipliers for specific heats
do n=2,5
    tc(n)=T*tc(n-1)   ! =T.^(n-1) 
end do

do k=1,ns
    if( T < thermo(k,1) ) then
        cpor(k)=dot_product(thermo(k,2:6), tc)  ! coefficients in lower temperature range
    else
        cpor(k)=dot_product(thermo(k,9:13), tc) ! coefficients in upper temperature range
    endif
end do

end subroutine ceq_cp