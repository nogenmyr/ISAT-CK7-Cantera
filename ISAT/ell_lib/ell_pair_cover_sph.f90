!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ell_pair_cover_sph( n, c1, gg1, c2, gg2, shrink, c, gg )

!  Given a pair of ellipsoids, E1 and E2, determine a third ellipsoid
!  E (specifically a spheroid) which covers E1 and E2.
!  E is defined by:  (x-c)^T * G * G^T * (x-c) <= 1, where the array gg
!  contains G in packed format.  Similarly for E1 and E2. 
!  A = G * G^T, A1 = G1 * G1^T, A2 = G2 * G2^T

!  Method:
!     Determine ball which covers bounding balls of E1 and E2.
!     Optionally shrink ball as much as possible (if shrink true).

!  S.B. Pope  7/3/2006, 11/25/2006

implicit none

integer,    parameter   :: k_dp = kind(1.d0)
integer,    intent(in)  :: n
real(k_dp), intent(in)  :: c1(n), gg1((n*(n+1))/2), c2(n), gg2((n*(n+1))/2)
logical,    intent(in)  :: shrink
real(k_dp), intent(out) :: c(n),  gg((n*(n+1))/2)

real(k_dp), parameter   :: extend = 1.d-6   !  small amount to extend radius to ensure covering  

integer    :: i, j, k
real(k_dp) :: d, r1, r2, r, smin, smax, dc(n), zero(n), xfar(n)
              
d	= sqrt( sum( (c2-c1)*(c2-c1) ) )  ! distance between centers
call ell_radii( n, gg1, r, r1 )       ! radius of ball covering E1
call ell_radii( n, gg2, r, r2 )       ! radius of ball covering E2

!  Consider line of centers, with distance s measured from c1 towards c2.
!  Thus c2 corresponds to s = d.

smin = min( -r1, -r2 + d )  !  intersections of E with line
smax = max(  r1,  r2 + d )

c = c1  !  set c to correspond to (smin+smax)/2
if( d > 0.d0 ) then
   c = c1 + 0.5d0*(smin+smax)*(c2-c1)/d
endif

if( shrink ) then
   zero = 0.d0  !  find furthest point in E1 from c
   dc   = c1 - c
   call ell_pt_near_far( 2, n, dc, gg1, zero, xfar )
   r1 = sum( xfar*xfar)
   
   dc   = c2 - c !  find furthest point in E2 from c
   call ell_pt_near_far( 2, n, dc, gg2, zero, xfar )
   r2 = sum( xfar*xfar)
   
   r  = sqrt( max(r1,r2) )  ! set r to max. distance from c  
else
   r = 0.5d0*(smax-smin)
endif

r = (1.d0 - extend) / r  !  diagonal components of G = 1/r

k  = 0  !  G into gg
gg = 0.d0
do j = 1, n
   do i = j, n
      k = k + 1
	  if( i==j ) gg(k) = r
   end do
end do

return
end subroutine ell_pair_cover_sph
