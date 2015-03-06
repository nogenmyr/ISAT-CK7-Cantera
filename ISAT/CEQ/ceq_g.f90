!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ceq_g(ns,T,p,thermo,gort)
!----------------
  implicit none
  integer, parameter :: ncof=7
  integer, intent(in) :: ns
  real(kind(1.d0)), intent(in) :: T, p, thermo(ns,2*ncof+1)
  real(kind(1.d0)), intent(out) :: gort(ns)
!-------------
!  return normalized Gibbs functions at temperature T
! input:
!   T      - temperature (K)
!   p      - pressure (atm)
!   thermo - thermo data for all species
! output:
!   gort   - g_j/(RT)  - normalized Gibbs functions

! S. B. Pope 9/26/02

real(kind(1.d0)) :: tc(ncof), th(ncof), ts(ncof), tg(ncof)
integer :: k, n

if( ns <= 0 ) return

tc=0.d0  ! coefficient multipliers for specific heats
th=0.d0  ! coefficient multipliers for enthalpy
ts=0.d0  ! coefficient multipliers for entropy
tc(1)=1.d0
th(1)=1.d0
th(6)=1./T
ts(1)=log(T)
ts(7)=1.d0
do n=2,5
    tc(n)=T*tc(n-1)   ! =T.^(n-1)
    th(n)=tc(n)/float(n)     ! =T.^(n-1) ./ n
    ts(n)=tc(n)/float((n-1)) ! =T.^(n-1) ./ (n-1)
end do
tg=th-ts

do k=1,ns
    if( T<thermo(k,1) ) then
        gort(k)=dot_product( thermo(k,2:8), tg)  ! coefficients in lower temperature range
    else
        gort(k)=dot_product( thermo(k,9:15), tg) ! coefficients in upper temperature range
    endif
end do

gort=gort+log(p)

end subroutine ceq_g