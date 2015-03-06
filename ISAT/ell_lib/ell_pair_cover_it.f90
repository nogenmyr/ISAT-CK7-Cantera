!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ell_pair_cover_it( n, c1, gg1, c2, gg2, cc, gg )

!  Given a pair of ellipsoids, E1 and E2, determine a third ellipsoid
!  E which covers E1 and E2.
!  E is defined by:  (x-cc)^T * G * G^T * (x-cc) <= 1, where the array gg
!  contains G in packed format.  Similarly for E1 and E2. 

implicit none

integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: n
real(k_dp), intent(in)  :: c1(n), gg1((n*(n+1))/2), c2(n), gg2((n*(n+1))/2)
real(k_dp), intent(out) :: cc(n),  gg((n*(n+1))/2)

integer, parameter ::  ludiag = 0

integer    :: i, j, k, info, lwork, iter, niter

real(k_dp) :: B1x(n,n), B2x(n,n), B2y(n,n), B4z(n,n), B12z(n,n,2), U(n,n), &
              Vt(n,n), temp(n,n), work(5*n*n+20*n), &
              g12(((n+1)*n)/2,2), gc(((n+1)*n)/2), &
			  dc(n), c12(n,2), sig(n), sigh(n), sigt(n), sigb(n),  &
              w(n), hv(n), tv(n), xfar(n,2), ctry(n), &
			  r1, r2, cdist,  spt(4), smin, smax, beta, re(2), rmax_best, &
			  tt, rho1b0, rho1b1, rho2b0, rho2b1, fb0, fb1, eta, err_last, &
			  f, err, eta_best, errtol, dfdeta, f_last, eta_last, r2_min, r2_max

errtol   = 1.d-2  ! error tolerance in root-finding
niter    = 10     ! maximum number of iterations

lwork = 5*n*n+20*n  ! amount of work space for SVD

!----- STAGE 1: given E1 and E2 ----------------------

! unpack gg1 and gg2 into B1x and B2x

B1x = 0.d0
B2x = 0.d0

k  = 0
do j = 1, n
   do i = j, n
      k = k + 1
	  B1x(i,j) = gg1(k)
	  B2x(i,j) = gg2(k)	  
   end do
end do

!------ STAGE 2: shift so that E1 is centered at the origin ------

dc = c2 - c1
cc = dc

!------ STAGE 3: transform (to y-space) so that E1 is the unit ball -------

!  B2y = B1x^{-1} * B2x
B2y = B2x
call dtrsm('L','L','N','N',n,n,1.d0,B1x,n,B2y,n)

! SVD:  B2y = U sig VT
temp(:,:) = B2y(:,:)
call dgesvd( 'A', 'A', n, n, temp, n, sig, U, n, Vt, n, work, lwork, info)

if( info /= 0 ) then
   write(ludiag,*)'ell_pair_cover_it: info, lwork, work(1) = ', info, lwork, work(1)
   stop
endif

!  c2y = B1x^T * (c2x-c1x)
cc = matmul( cc, B1x )

!  determine whether E1 covers E2 -------------------------

r2_min = 1.d0 / maxval(sig)    !  radius of E2's inscribed ball
r2_max = 1.d0 / minval(sig)    !  radius of E2's bounding ball
cdist  = sqrt( sum( cc*cc ) )  !  distance between centers

if( cdist + r2_max <= 1.d0 ) then  
   cc = c1   !  E1 covers E2; return E = E1
   gg = gg1
   return
elseif( r2_max <= 1.d0  .and.  cdist + r2_min <= 1.d0 ) then
!  E1 may cover E2:  determine distance from center of E1 (i.e. the origin)
!  to the furthest point in E2.
   call ell_bbt2chol( n, B2y, gc )
   ctry = 0.d0
   call ell_pt_near_far( 2, n, cc, gc, ctry, tv )
   
   if( sum( tv*tv ) <= 1.d0 ) then
      cc = c1   !  E1 covers E2; return E = E1
      gg = gg1
      return
   endif
endif

!  determine whether E2 covers E1 --------------------------

if( r2_min >= cdist + 1.d0 ) then  
   cc = c2   !  E2 covers E1; return E = E2
   gg = gg2
   return
elseif( r2_min >= 1.d0  .and.  r2_max >= cdist + 1.d0 ) then
!  E2 may cover E2:  determine distance from center of E1 (i.e. the origin)
!  to the furthest point in E2.
!  Transform E2 to unit ball at the origin.  Determine transformed E1...

   hv = -matmul( cc, U ) * sig  !  center of transformed E1

   do i = 1, n
      temp(i,:) = U(:,i)/sig(i)  !  B matrix for transformed E1
   end do

   call ell_bbt2chol( n, temp, gc )
   ctry = 0.d0
   call ell_pt_near_far( 2, n, hv, gc, ctry, tv )
   
   if( sum( tv*tv ) <= 1.d0 ) then
      cc = c2   !  E2 covers E2; return E = E2
      gg = gg2
      return
   endif
endif

!  E1 and E2 do not cover each other

!  define E3 which covers E1 and E2^0
do i = 1, n
   sigh(i) = min( sig(i), 1.d0 )
   sigt(i) = 1.d0 / max( sig(i), 1.d0 )
end do

!------ STAGE 4: transform so that E3 is the unit ball --------------------

cc = sigh * matmul( cc, U ) 
cdist = sqrt( sum( cc*cc ) )

r1 = maxval( sigh )  !  radius of ball bounding E1
r2 = maxval( sigt )  !  radius of ball bounding E2

!---- treat simple degenerate case in which c1=c2.
!  G*G^T = B*B^T, B=B1*U*Sigh
if( cdist == 0.d0 ) then
   B2y = matmul( B1x, U )
   do j = 1, n
      B2x(:,j) = B2y(:,j) * sigh(j)
   end do
   call ell_bbt2chol( n, B2x, gg )
   cc = c1
   return
endif
!-------------------------------------------------

!  determine unit vector w and distance cdist between the centers
w     = cc / cdist
spt   = (/ -r1, r1, cdist-r2, cdist+r2 /)  ! intersections on line and balls
smin  = minval( spt )
smax  = maxval( spt )

!  define E4: E3 stretched in w-direction to intersect extrema

call ell_house( n, w, hv, beta )  ! Householder vector

sigb    = 1.d0  !  \bar{Sigma}
sigb(1) = 2.d0 / ( smax - smin )

!  B_4 = P sigb = (I-beta*hv*hv^T) * sigb

tv    = beta * hv
tv(1) = tv(1) * sigb(1)

do i = 1, n
   do j = 1, n
      B4z(i,j) = -hv(i)*tv(j)
   end do
   B4z(i,i) = B4z(i,i) + sigb(i)
end do 

!------ STAGE 5: transform so that E4 is the unit ball --------------------

do i = 1, n
   tv    = 0.d0
   tv(i) = 1.d0  !  identity
   do j = 1, n
      B12z(i,j,1) = (tv(j) - beta * hv(i) * hv(j)) / ( sigb(i) * sigh(j) )
      B12z(i,j,2) = (tv(j) - beta * hv(i) * hv(j)) / ( sigb(i) * sigt(j) )
   end do
end do

call ell_bbt2chol( n, B12z(:,:,1), g12(:,1) )  ! Cholesky factorization of B(1:2,5)
call ell_bbt2chol( n, B12z(:,:,2), g12(:,2) )

cc    = cc - beta * hv * sum( hv * cc )
cc(1) = cc(1) * sigb(1) 

!--------------- root finding -----------------------------
!  eta = fractional distance along line of centers E1-E2
!  rho1(eta) = furthest distance to point in E1
!  rho2(eta) = furthest distance to point in E2
!  f(eta)    = rho1(eta) - rho2(eta)

tt     = 0.5d0 * ( cdist + r1 + r2 )  !  compression

rho1b0 = r1  !  estimates of f(0) and f(1)
rho1b1 = sqrt( max( ((cdist+r1)/tt)**2 , (cdist/tt)**2+r1**2 ) ) 

rho2b0 = sqrt( max( ((cdist+r2)/tt)**2 , (cdist/tt)**2+r2**2 ) ) 
rho2b1 = r2

fb0 = rho1b0 - rho2b0;
fb1 = rho1b1 - rho2b1;

eta = -fb0 / (fb1-fb0)  !  initial guess

tv = cc

err_last = huge(1.d0)

c12(:,1) = 0.d0
c12(:,2) = cc

do iter = 1, niter

!  evaluate rho1/2(eta)

   ctry = eta * cc
   do k = 1, 2
      call ell_pt_near_far( 2, n, c12(:,k), g12(:,k), ctry, xfar(:,k) )
      tv = xfar(:,k) - ctry
      re(k) = sqrt( sum(tv*tv) )
   end do

   f = re(1) - re(2)
   err = abs(f)

   if( err < err_last ) then  ! best solution so far
      err_last  = err
	  eta_best  = eta
	  rmax_best = maxval( re )
   else
      exit  ! diverging -- give up
   endif

   if( err <= errtol ) exit  !  converged
      
   if( iter == 1 ) then  !  estimate df/d(eta)
      dfdeta = fb1 - fb0
   else
      dfdeta = (f - f_last)/(eta - eta_last)
   endif

   eta_last = eta
   f_last   = f

   eta = eta -f / dfdeta  !  value for next iteration

   if( eta < -1.d0  .or.  eta >= 2.d0 ) exit
end do

!-------------------------------------------------------------------------------

!define E5 to be expanded E4


!------ STAGE 6: transform back to original x-space --------------------

do i = 1, n
   temp(:,i) = U(:,i) / sigh(i)
end do

call dtrsm('L','L','T','N',n,n,1.d0,B1x,n,temp,n)
cc = c1 + eta_best * dc

do i = 1, n
   temp(i,:) = sigh(i) * B4z(i,:)
end do

temp = matmul( U, temp )
temp = matmul( B1x, temp ) / rmax_best

call ell_bbt2chol( n, temp, gg )

return

end subroutine ell_pair_cover_it
