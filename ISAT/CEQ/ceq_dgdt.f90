!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ceq_dgdT(ns,T,thermo,dgdT)
!----------------
  implicit none
  integer, parameter :: ncof=7
  integer, intent(in) :: ns
  real(kind(1.d0)), intent(in) :: T, thermo(ns,2*ncof+1)
  real(kind(1.d0)), intent(out) :: dgdT(ns)
!-------------
!  return d/dT of the normalized Gibbs functions at temperature T
! input:
!   T      - temperature (K)
!   thermo - thermo data for all species
! output:
!   dgdT   - d/dT (g_j/(RT) )

! S. B. Pope 7/1/03

real(kind(1.d0)) :: tc(ncof), th(ncof), ts(ncof), tg(ncof)
integer :: k, n

tc=0.d0  ! coefficient multipliers for specific heats
th=0.d0  ! coefficient multipliers for enthalpy
ts=0.d0  ! coefficient multipliers for entropy
tc(1)=1.d0/T
th(1)=0.d0
th(6)=-1.d0/T**2
ts(1)=1.d0/T
ts(7)=0.d0

do n=2,5
    tc(n)=T*tc(n-1)                  ! =T.^(n-2)
    th(n)=tc(n)*float(n-1)/float(n)  ! =T.^(n-2) * (n-1) / n
    ts(n)=tc(n)                      ! =T.^(n-2) 
end do
tg=th-ts

do k=1,ns
    if( T<thermo(k,1) ) then
        dgdT(k)=dot_product( thermo(k,2:8), tg)  ! coefficients in lower temperature range
    else
        dgdT(k)=dot_product( thermo(k,9:15), tg) ! coefficients in upper temperature range
    endif
end do

end subroutine ceq_dgdT