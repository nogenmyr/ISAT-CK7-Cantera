!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ell_rad_upper( n, gg, r )

!  Determine the radius r of the ball covering the ellipsoid E given
!  by { x | norm(G^T * x) <=1 ), where G is an n x n lower triangular
!  matrix.  The array gg contains the matrix G in packed format.

!  S.B. Pope  6/12/04

implicit none

integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: n
real(k_dp), intent(in)  :: gg((n*(n+1))/2)
real(k_dp), intent(out) :: r

integer    :: itmax = 100  !  max. iterations in dgqt
real(k_dp) :: atol = 1.d-4, rtol = 1.d-4  !  tolerances for dgqt

integer    :: info, i, j, k
real(k_dp) :: gi((n*(n+1))/2), g(n,n), alpha, delta, par, f
real(k_dp) :: a(n,n), b(n), x(n), z(n), wa1(n), wa2(n)

!  Method:  E = { y | y^T * y <= 1 } where y = G^T * x.
!  Now  x^T * x = y^T * G^{-1} * G^{-T} * y.
!  Thus, with f(y) = -2 * y^T * G^{-1} * G^{-T} * y, 
!  r = sqrt(-f(y_m) ), where y_m is the minimizer of f, subject to |y|<=1.

gi = gg
call dtptri( 'L', 'N', n, gi, info )  !  G^{-1} in packed format

!  unpack gi into g
k = 0
g = 0.d0
do j = 1, n
   do i = j, n
      k = k + 1
	  g(i,j) = gi(k)
   end do
end do
a = g

! B = alpha * B * op(A) = -2 * G^{-1} * G^{-T}

alpha = -2.d0
call dtrmm ( 'R', 'L', 'T', 'N', n, n, alpha, g, n, a, n )

delta = 1.d0
par   = 0.d0
b     = 0.d0
!  minimize f(y) = (1/2) * y^T * A * y + b^T * y,  subject to |y|<= delta
!                = -2 * y^T * G^{-1} * G^{-T} * y

call dgqt(n,a,n,b,delta,rtol,atol,itmax,par,f,x,info,z,wa1,wa2)

if( info >= 3 ) then
   write(0,*)'ell_rad_upper: dgqt incomplete convergence, info = ', info
endif

if( f > 0.d0 ) then
   write(0,*)'ell_rad_upper: f positive, f = ', f
   stop
endif

r = sqrt( -f ) 

return
end subroutine ell_rad_upper
