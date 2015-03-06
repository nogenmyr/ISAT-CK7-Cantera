!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ell_pt_shrink( n, c, gg, p, k_ell, info )

!  Given the ellipse E = { x | (x-c)^T * G * G^T * (x-c) <= 1 }, 
!  where G is lower triangular, and the point p which is covered by E, we define:
! 1/ Ev: the maximum volume ellipse covered by E with p on its boundary
! 2/ En: the ellipse of maximal near content covered by E with p on its
! boundary
! 3/ Ec: the minimum volume ellipse covering Ev and En.
! For k=|k_ell|=1,2 or 3,  G is modified to correspond to Ev, En or Ec, respectively.
!  G (in packed format) is contained in gg.
!
! info = 0 - p is in the interior of E (intended case): G is returned corresponding to Ev, En or Ec.
! info = 1 - p is on the boundary of E (within tol): G is returned unchanged
! info = 2 - p is not covered by E: G is returned unchanged

!  S.B. Pope  12/26/05

implicit none

integer,    parameter     :: k_dp = kind(1.d0)
integer,    intent(in)    :: n, k_ell
real(k_dp), intent(in)    :: c(n), p(n)
real(k_dp), intent(inout) :: gg( (n*(n+1))/2 )
integer,    intent(out)   :: info

real(k_dp), parameter     :: tol = 1.d-10

integer    :: j, k, lwork
real(k_dp) :: U(n,n), sig(n), ph(n), phsqi, wt(n), pt(n), R(n,2), tau(n), work(n+10), &
              Lb(2,2), Lch(n,n), B(n,n), pch(n), dpsq
logical    :: in

lwork = n+10
info  = 0

!  check input
k = iabs( k_ell )
if( k < 1  .or.  k > 3 ) then
   write(0,*) 'ell_pt_shrink: invalid k_ell = ', k_ell
   stop
endif

if( k_ell == 1 ) then  !  efficient determination of Ev (avoids SVD)
   call ell_pt_in( n, c, gg, p, in )
   if( .not.in ) then
      info = 2
	  return
   endif

   call ell_pt_modify( n, c, gg, p )
   return
endif

call ell_chol2eig( n, gg, U, sig )  !  SVD  G = U * diag(sig)

ph  = matmul( p-c , U )  ! ph   = U'*(p-c) 
pt  = sig*ph
dpsq     = sum( pt*pt )

if( dpsq > 1.d0 + tol ) then
   info = 2  !  p is not covered by E
   return
elseif( dpsq >= 1.d0 - tol ) then
   info = 1  !  p is on the boundary of E
   return  
elseif( dpsq == 0.d0 ) then
   write(0,*) 'ell_pt_shrink_cons: p = c'
   stop
endif

if( n==1 ) then !  special case, n=1
   gg(1) = 1.d0 / abs(p(1)-c(1))
   return
endif

phsqi= 1.d0/sum(ph*ph)
do j=1,n 
   wt(j)=max(0.d0, phsqi-sig(j)**2) 
end do

wt  = wt*ph/sig

R(:,1) = wt
R(:,2) = pt

call dgeqrf( n, 2, R, n, tau, work, lwork, info )  ! [Q,R] = qr([wt pt])

if( info /= 0 ) then
   write(0,*) 'ell_shrink: dgeqrf failed, info = ', info
   stop
endif

pch = R(:,2)  !  pch   = Q'*pt

call ell_cov_2d( pch(1:2), k, Lb )

Lch = 0.d0
do j = 1, n
   Lch(j,j) = 1.d0
end do

Lch(1:2,1:2) = Lb(1:2,1:2)

call dormqr( 'L', 'N', n, n, 2, R, n, tau, Lch, n, work, lwork, info )  ! Q*Lch

if( info /= 0 ) then
   write(0,*) 'ell_shrink_cons: dormqr failed, info = ', info
   stop
endif

do j = 1, n
  Lch(j,:) = sig(j) * Lch(j,:)  ! Sig*Q*Lch
end do

B = matmul( U, Lch )  !  B = U*Sig*Q*Lch

call dgelqf( n, n, B, n, tau, work, lwork, info )  !  L * Q = B

if( info /= 0 ) then
   write(0,*) 'ell_shrink_cons: dgelqf failed, info = ', info
   stop
endif

call ell_full2low( n, B, gg )

end subroutine ell_pt_shrink