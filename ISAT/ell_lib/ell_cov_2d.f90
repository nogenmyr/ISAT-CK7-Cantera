!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ell_cov_2d( p, k_ell, L )

!  In 2D, given the point p (|p|<1, |p_1|>0), three ellipses are defined:
!  E_1 is the unit circle, shrunk in the p direction to intersect p 
!  E_2 is the unit circle, shrunk in the x_1 direction to intersect p, and
!  E_3 is the minimum-area ellipse covering E_1 and E_2.
!  The 2 x 2 lower trianglular matrix L is returned, 
!  where E_{k_ell} = { x | |L^T * x| <1 }.

!  S.B. Pope  12/16/05

implicit none

integer, parameter      :: k_dp = kind(1.d0)
real(k_dp), intent(in)  :: p(2)
integer,    intent(in)  :: k_ell
real(k_dp), intent(out) :: L(2,2)

real(k_dp), parameter   :: tol = 1.d-12

real(k_dp) :: psq, chi, theta, a, b, c, phi, sp, cp, sig1sq, sig2sq, sig1, &
              sig2, alpha, d11, d21, d12, d22, ca, sa


L = 0.d0  ! set L to unit circle
L(1,1) = 1.d0
L(2,2) = 1.d0

!  check that input is valid
if( abs(p(1))==0.d0 ) then
   write(0,*) 'ell_cov_2d: p_1 is zero'
   return
endif

if( k_ell < 1  .or.  k_ell > 3 ) then
   write(0,*) 'ell_cov_2d: invalid k_ell = ', k_ell
   stop
endif

psq = sum(p*p)
if( psq > 1.d0+tol ) then  !  treat boundary and exterior
   write(0,*) 'ell_cov_2d: |p| > 1, ', psq-1.d0
   return 
elseif( psq >= 1.d0-tol ) then
   return
endif

chi   = abs(p(1))/sqrt(1-p(2)**2)  ! intersection of E_2 and x_1 axis

if( k_ell == 2 ) then  !  short-cut for E_2
   L(1,1) = 1.d0 / chi
   return
endif

theta = (1-psq)/psq**2             ! E_1={x | x^T*(I+theta*p*p^T)*x <1 }

!  transform so that E_1 is the unit circle: (x,y)->(x/chi,y)
a = (1.d0+theta*p(1)**2)*chi**2     ! B=[a cc b] is the matrix describing E_1
b = 1.d0 + theta*p(2)**2
c = theta*chi*p(1)*p(2)

!  eigendecomposition of B
phi = .5d0*atan2(-2.d0*c,(a-b))
sp  = sin(phi)
cp  = cos(phi)

sig1sq = a*cp*cp-2.d0*c*sp*cp+b*sp*sp
sig2sq = a*sp*sp+2.d0*c*sp*cp+b*cp*cp

!  check that sigs bracket unity
if( min(sig1sq,sig2sq) >1.d0 + tol .or. max(sig1sq,sig2sq) <1.d0 - tol ) then 
    write(0,*) 'ell_cov_2d: sig1sq-1, sig2sq-1= ' ,sig1sq-1.d0, sig2sq-1.d0
endif

! set principal axes for E_1, E_2 or E_3
!  no change needed for E_1

if( k_ell == 2 ) then
    sig1sq=1.d0  ! set to E_2
	sig2sq=1.d0 
    
elseif( k_ell ==3 ) then
    sig1sq = min(sig1sq, 1.d0)  ! E_3
    sig2sq = min(sig2sq, 1.d0)
endif

!  the matrix for E_3 is: D*D^T,  D=C*U*Sig = [q r s t]
sig1 = sqrt(sig1sq)
sig2 = sqrt(sig2sq)
d11  = sig1*cp/chi
d12  = sig2*sp/chi
d21  =-sig1*sp
d22  = sig2*cp
!!  D = [d11 d12 d21 d22]
!  transform back
alpha = atan2(-d12, d11)
ca = cos(alpha)
sa = sin(alpha)
!!  QT = [ca sa; -sa ca]
!! L  = D*QT
L(1,1) = ca*d11 - sa*d12
L(2,1) = ca*d21 - sa*d22
L(2,2) = sa*d21 + ca*d22

return

end subroutine ell_cov_2d


