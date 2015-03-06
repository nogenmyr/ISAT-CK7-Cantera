!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ell_pair_separate( n, c1, gg1, c2, gg2, qual, max_it, xh, v, intersect )

!  Determine if two non-concentric ellipsoids, E1 and E2, intersect.
!  If they do not intersect, determine a separating hyperplane.
!  If they do intersect, determine the hyperplane which is the 
!     perpendicular bisector of the line of centers.

!  E1 is defined by:  (x-c1)^T * G1 * G1^T * (x-c1) <= 1.  The array
!  gg1 contains the lower triangular n x n matrix G1 in packed format.
!  Similarly for E2.  

!  The hyperplane is defined by v^T * ( x - xh  ) = 0, where v is
!  a unit vector.  The quantity s(x) = v^T * ( x - xh  ) is the signed
!  distance of the point x from the hyperplane: s(c1) is negative, and 
!  s(c2) is positive.  If E1 and E2 do not intersect, then s(x) is negative
!  for all points x in E1, and positive for all points x in E2.  

!  Method:
!
!Phase I
!  1/ transform to y-space in which E1 is the unit ball: y = G1^T * ( x - c1 )
!  2/ find the point y2 in E2 that is closest to the origin
!  3/ if |y2| <= 1, then E1 and E2 intersect; otherwise...
!  4/ define y1 = y2/|y2|;  yh = (y1+y2)/2;  vy = y2-y1; then vy^T * ( y - yh )
!     is a separating hyperplane
!  5/ transform back to x-space to obtain xh and v (and x1 and x2)
!  6/ all done if qual <=0 or max_it <=0
!
!Phase II
!  The objective of phase II is to improve the quality of the separating hyperplane.
!  The quality q ( q <= 1 ) is defined by: 
!     q= (distance between supporting hyperplanes)/|x1-x2|.
!  Starting from the value of x2 obtained in phase I, iterations are performed to:
!  7/ determine the point x1 in E1 closest to x2
!  8/ determine the point x2 in E2 closest to x1
!  9/ xh = (x1+x2)/2;  v = (x2-x1)/|x2-x1|
!  These operations (7-9) are performed iteratively until q >= qual, or the number
!  of iterations exceeds max_it.  The values of xh and v returned correspond to the
!  hyperplane of greatest separation encountered during (or before) the iterations
!  (which usually occurs on the final iteration). 
!
!  Notes:
!  1/ For 7 and 8, the same algorithm is used as in ell_pt_near_far, but there are 
!     efficiency gains in re-coding it here.
!  2/ A quadratic programming routine could be used to solve the whole problem.
!  3/ Convergence is slow if E1 and E2 nearly intersect; but then the separating
!     hyperplane obtained from phase I may be adequate.
!  4/ Aitken extrapolation could be used to speed up convergence.

!  S.B. Pope  6/13/04, 12/28/04

implicit none

integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: n, max_it
real(k_dp), intent(in)  :: c1(n), gg1((n*(n+1))/2), c2(n), gg2((n*(n+1))/2), qual
real(k_dp), intent(out) :: xh(n), v(n)
logical, intent(out)    :: intersect


integer    :: itmax = 100  !  max. iterations in dgqt
integer    :: lu_diag = -1 !  logical unit for diagnostics (set < 0 to suppress)
real(k_dp) :: atol = 1.d-6, rtol = 1.d-6  !  tolerances for dgqt

integer    :: info, i, j, k, iter
real(k_dp) :: g1(n,n), g2(n,n), c2y(n), y2(n), y2norm, x1(n), x2(n), g2y(n,n), &
              dist, smin1, smax1, smin2, smax2, q, qual_tol, sep, sep_max, &
			  v_best(n), xh_best(n), gi(n,n), delta, par1, par2, f, &
              a1(n,n), a2(n,n), b(n), z(n), wa1(n), wa2(n)

! unpack gg1 and gg2  -----------------------------------------------------------
k  = 0
g1 = 0.d0
g2 = 0.d0
do j = 1, n
   do i = j, n
      k = k + 1
	  g1(i,j) = gg1(k)
	  g2(i,j) = gg2(k)	  
   end do
end do

! transform to y-space:  y = G1^T * (x-c1);  G2y = G1^{-1} * G2
g2y = g2
call dtrsm( 'L', 'L', 'N', 'N', n, n, 1.d0, g1, n, g2y, n )

c2y = c2 - c1  !  c2y = G1^T * ( c2 - c1 )
call dtrmv( 'L', 'T', 'N', n, g1, n, c2y, 1 )

!  find point y2 in E2 closest to the origin
v = 0.d0  ! origin
call ellu_pt_near_far( 1, n, c2y, g2y, v, y2 ) 

y2norm = sum( y2*y2 )  !  if y2 in unit ball then E1 and E2 intersect

if( y2norm <= 1.d0 ) then
   intersect = .true.           ! E1 and E2 intersect: all done
   xh   = 0.5d0 * (c1 + c2)     ! mid-point between centers
   v    = c2 - c1
   dist = sum( v*v )
   if( dist > 0.d0 ) then
      v = v / sqrt( dist )      ! unit vector from c1 to c2
   else
      v = 0.d0                  ! E1 and E2 are concentric
	  v(1) = 1.d0               ! set v = e_1
   endif

   return         
endif

intersect = .false.

!  hyperplane in y-space is:  vy^T * ( y - yh ) = 0
y2norm = sqrt( y2norm )
v      = y2  
            
!  transform back to x-space
call dtrsv( 'L', 'T', 'N', n, g1, n, y2, 1 )
x2 = y2        + c1  ! x2 = G1^{-T} * y2 + c1
x1 = y2/y2norm + c1
xh = 0.5d0 * ( x1 + x2 )

call dtrmv( 'L', 'N', 'N', n, g1, n, v, 1 )  ! v = G1 * vy
v     = v / sqrt( sum(v*v) )  !  make unit vector

!-----------  end of phase I  -------------------------------------------

if( max_it <=0  .or.  qual <= 0.d0 ) then
   return  !  accept existing separating hyperplane
else
   qual_tol = min( qual, 1.d0-1.d-6 )  !  qual must be less than unity
endif

!---  form matrices A1 and A2 used in minimization ----

gi = g1  !  start formation of a1
call dtrtri( 'L', 'N', n, gi, n, info )  !  gi = G1^{-1} 

a1 = 0.d0  !  a1 = lower triangle of gi 
do j = 1, n
   a1(j:n,j) = gi(j:n,j)
end do

call dtrmm ( 'R', 'L', 'T', 'N', n, n, 1.d0, gi, n, a1, n )  ! a1 =  = G1^{-1} * G1^{-T}

gi = g2  !  start formation of a2
call dtrtri( 'L', 'N', n, gi, n, info )  !  gi = G2^{-1} 

a2 = 0.d0  !  a2 = lower triangle of gi 
do j = 1, n
    a2(j:n,j) = gi(j:n,j)
end do

call dtrmm ( 'R', 'L', 'T', 'N', n, n, 1.d0, gi, n, a2, n )  ! a2 =  = G2^{-1} * G2^{-T}

!---  initialization prior to iteration
par1  = 0.d0
par2  = 0.d0
delta = 1.d0

v_best  = v
xh_best = xh

call ell_line_proj( n, c1, gg1, xh, v, smin1, smax1 )
call ell_line_proj( n, c2, gg2, xh, v, smin2, smax2 )
sep_max = smin2 - smax1  !  distance between supporting hyperplanes

do iter = 1, max_it  !------------------  iterations  --------------------------

!  determine point x1 in E1 closest to x2 ----------
   b = c1 - x2
   call dtrsv( 'L', 'N', 'N', n, g1, n, b, 1 )  !  b = G1^{-1} * [c1-x2]

   call dgqt(n,a1,n,b,delta,rtol,atol,itmax,par1,f,x1,info,z,wa1,wa2)

   if( info >= 3 ) then
      write(0,*)'ell_pair_separate: dgqt incomplete convergence, info = ', info
   endif

   call dtrsv( 'L', 'T', 'N', n, g1, n, x1, 1 )  !  x1 = G1^{-T} * y1 + c1
   x1 = x1 + c1

!  determine point x2 in E2 closest to x1 ----------

   b = c2 - x1
   call dtrsv( 'L', 'N', 'N', n, g2, n, b, 1 )  !  b = G2^{-1} * [c2-x1]

   call dgqt(n,a2,n,b,delta,rtol,atol,itmax,par2,f,x2,info,z,wa1,wa2)

   if( info >= 3 ) then
      write(0,*)'ell_pair_separate: dgqt incomplete convergence, info = ', info
   endif

   call dtrsv( 'L', 'T', 'N', n, g2, n, x2, 1 )  !  x2 = G2^{-T} * y2 + c2
   x2 = x2 + c2

!  determine hyperplane, dist=distance(x1,x2), sep=separation distance between
!  supporting hyperplanes, and q = quality = sep/dist         -------------

   v    = x2 - x1
   dist = sqrt( sum( v*v ) )  !  distance between "close" points
   v    = v / dist            !  unit vector from x1 to x2
   xh   = 0.5d0 * ( x1 + x2 ) !  mid-point between x1 and x2

   call ell_line_proj( n, c1, gg1, xh, v, smin1, smax1 )
   call ell_line_proj( n, c2, gg2, xh, v, smin2, smax2 )

   sep = smin2 - smax1  !  distance between supporting hyperplanes
   q   = sep / dist     !  quality of separating hyperplane

   if( lu_diag >= 0 ) write(lu_diag,'(a,i4,1p,10e13.4)')  &
      'iter, sep, dist, q = ', iter, sep, dist, q

   if( sep > sep_max ) then
      sep_max = sep  !  store xh and v with largest separation
	  xh_best = xh
	  v_best  = v
   endif

   if( q >= qual_tol ) exit

end do

xh = xh_best  !  return xh and v with largest separation
v  =  v_best  

return
end subroutine ell_pair_separate
