!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ellu_rad_upper( n, g, r )

!  Determine the radius r of the ball covering the ellipsoid E given
!  by { x | norm(G^T * x) <=1 ), where G is an n x n lower triangular
!  matrix.  The array g contains the matrix G (unpacked).

!  S.B. Pope  6/12/04

implicit none

integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: n
real(k_dp), intent(in)  :: g(n,n)
real(k_dp), intent(out) :: r

integer    :: itmax = 100  !  max. iterations in dgqt
real(k_dp) :: atol = 1.d-4, rtol = 1.d-4  !  tolerances for dgqt

integer    :: info
real(k_dp) :: gi(n,n), alpha, delta, par, f
real(k_dp) :: a(n,n), b(n), x(n), z(n), wa1(n), wa2(n)

!  Method:  E = { y | y^T * y <= 1 } where y = G^T * x.
!  Now  x^T * x = y^T * G^{-1} * G^{-T} * y.
!  Thus, with f(y) = -2 * y^T * G^{-1} * G^{-T} * y, 
!  r = sqrt(-f(y_m) ), where y_m is the minimizer of f, subject to |y|<=1.

gi = g
call dtrtri( 'L', 'N', n, gi, n, info )  !  gi = G^{-1}

if( info /= 0 ) then
   write(0,*)'ellu_rad_upper: dtrtri failed, info = ', info
   stop
endif

a = gi

! B = alpha * B * op(A) = -2 * G^{-1} * G^{-T}

alpha = -2.d0
call dtrmm ( 'R', 'L', 'T', 'N', n, n, alpha, gi, n, a, n )

delta = 1.d0
par   = 0.d0
b     = 0.d0
!  minimize f(y) = (1/2) * y^T * A * y + b^T * y,  subject to |y|<= delta
!                = -2 * y^T * G^{-1} * G^{-T} * y

call dgqt(n,a,n,b,delta,rtol,atol,itmax,par,f,x,info,z,wa1,wa2)

if( info >= 3 ) then
   write(0,*)'ellu_rad_upper: dgqt incomplete convergence, info = ', info
endif

if( f > 0.d0 ) then
   write(0,*)'ellu_rad_upper: f positive, f = ', f
   stop
endif

r = sqrt( -f ) 

return
end subroutine ellu_rad_upper
