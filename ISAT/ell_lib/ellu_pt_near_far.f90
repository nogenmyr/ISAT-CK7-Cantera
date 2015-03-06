!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ellu_pt_near_far( knf, n, c, g, p, xnf )

!  Determine the point xnf in the specified ellipsoid E that is
!  closest (for knf=1) or furthest (for knf=2) from the specified point p.  
!  The ellipsoid E is given by:
!  { x | norm(G^T * (x-c) ) <=1 ), where G is an n x n lower triangular
!  matrix.  The array g contains the matrix G (unpacked format).

!  S.B. Pope  6/12/04, 12/1/2006

implicit none

integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: knf, n
real(k_dp), intent(in)  :: c(n), g(n,n), p(n)
real(k_dp), intent(out) :: xnf(n)

integer    :: itmax = 100  !  max. iterations in dgqt
real(k_dp) :: atol = 1.d-12, rtol = 1.d-4  !  tolerances for dgqt

integer    :: info, i, j
real(k_dp) :: gi(n,n), delta, par, f, sig_min, sig_max, r, dist
real(k_dp) :: a(n,n), b(n), ynf(n), z(n), wa1(n), wa2(n), sign

!  Method:  E = { y | y^T * y <= 1 } where y = G^T * (x-c),  x = c + G^{-T} * y
!  Now  (x-p)^T * (x-p) = ( [c-p]^T + y^T * G^{-1} ) * ( G^{-T} * y + [c-p] ).
!  Thus in y-space, the nearest point ynf minimizes
!     fn(y) =  1/2 * y^T * G^{-1} * G^{-T} * y + [c-p]^T * G^{-T} * y, 
!           =  1/2 y^T * An * y + bn^T * y,  where
!        An = G^{-1} * G^{-T},  bn = G^{-1} * [c-p]
!  Similarly the furthest point ynf minimizes
!     ff(y) =  -1/2 * y^T * G^{-1} * G^{-T} * y - [c-p]^T * G^{-T} * y, 
!           =   1/2 y^T * Af * y + bf^T * y,  where
!        Af =  -G^{-1} * G^{-T},  bf = -G^{-1} * [c-p]

!  xnf = G^{-T} * ynf + c

if( knf == 1 ) then
   sign = 1.d0   ! for nearest point
elseif( knf == 2 ) then
   sign = -1.d0  ! for furthest point
else
   write(0,*)'ell_pt_near_far: bad value of knf = ', knf
   stop
endif

!-----------  special treatment for E being a spheroid
!  examine diagonal
sig_min =  huge(1.d0)
sig_max = -huge(1.d0)

do j = 1, n
   sig_min = min( abs(g(j,j)), sig_min )
   sig_max = max( abs(g(j,j)), sig_max )
end do

if( 1.d0 - sig_min/sig_max < rtol ) then
! diagonal is constant (within rtol); check off-diagonals
   sig_max = -huge(1.d0)
   do i = 2, n
      do j = 1, i-1
         sig_max = max( abs(g(i,j)), sig_max )
      end do
   end do
   if( n==1 ) sig_max = 0.d0  !  no off-diagonals for n=1 
   
   if( sig_max/sig_min < rtol ) then
!  E is a spheroid of radius 1/sig_min
      r    = 1.d0 / sig_min
      xnf  = c - p
      dist = sqrt( sum(xnf*xnf) )
      
      if( knf == 1 ) then  !  nearest point
         if( dist <= r ) then
            xnf = p
         else
            xnf = c - r * xnf / dist
         endif
         
      else  !  furthest point
      
         if( dist > 0.d0 ) then
            xnf = c + r * xnf / dist
         else
            xnf    = c
            xnf(1) = xnf(1) + r  !  arbitrarily chose 1-direction
         endif
      
      endif
      return  !  all done for E being a spheroid
      
   endif
endif

!---------- treat E not being a spheroid

gi = g
call dtrtri( 'L', 'N', n, gi, n, info )  !  G^{-1} 

b = sign * (c-p)
call dtrmv( 'L', 'N', 'N', n, gi, n, b, 1 )  !  b = G^{-1} * [c-p]

a = 0.d0  !  a = lower triangle of gi 


do j = 1, n
   a(j:n,j) = gi(j:n,j)
end do

! B = alpha * B * op(A) = G^{-1} * G^{-T}

call dtrmm ( 'R', 'L', 'T', 'N', n, n, sign, gi, n, a, n )

delta = 1.d0
par   = 0.d0

!  minimize f(y) = (1/2) * y^T * A * y + b^T * y,  subject to |y|<= delta

call dgqt(n,a,n,b,delta,rtol,atol,itmax,par,f,ynf,info,z,wa1,wa2)

if( info >= 3 ) then
   write(0,*)'ellu_pt_near_far: dgqt incomplete convergence, info = ', info
   write(0,*)'n = ', n, '  a, b ='
   write(0,'((1p,3e13.4))') a
   write(0,'((1p,3e13.4))') b
endif

if( knf == 1 ) then
   xnf = ynf
else
   xnf = ynf / sqrt( sum(ynf*ynf) )  !  small correction to enforce unit length
endif

call dtrmv( 'L', 'T', 'N', n, gi, n, xnf, 1 )  !  xnf = G^{-T} * ynf + c
xnf = xnf + c

return
end subroutine ellu_pt_near_far
