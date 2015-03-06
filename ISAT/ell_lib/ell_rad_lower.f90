!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ell_rad_lower( n, gg, r )

!  Determine the radius r of the ball inscribed in the ellipsoid E given
!  by { x | norm(G^T * x) <=1 ), where G is an n x n lower triangular
!  matrix.  The array gg contains the matrix G in packed format.

!  Method:  r = 1 / sqrt(lam_max)  where lam_max is the largest
!           eigenvalue of G * G^T

!  S.B. Pope 6/12/04

implicit none

integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: n
real(k_dp), intent(in)  :: gg((n*(n+1))/2)
real(k_dp), intent(out) :: r

integer    :: info, i, j, k, nev, lwork, iwork(5*n), ifail(n)
real(k_dp) :: g(n,n), a(n,n), dum(1), ev(n), dumz(n,1), work(n*(n+8))

lwork = n*(n+8)

!  unpack gg into g
k = 0
g = 0.d0
do j = 1, n
   do i = j, n
      k = k + 1
	  g(i,j) = gg(k)
   end do
end do

!  form a = A = G * G^T
a = g
call dtrmm ( 'R', 'L', 'T', 'N', n, n, 1.d0, g, n, a, n )

!  find first (smallest) eigenvalues of A

call dsyevx( 'N', 'I', 'L', n, a, n, dum, dum, n, n, 0.d0, nev, ev, dumz, n, &
             work, lwork, iwork, ifail, info )

if( info /= 0 ) then  !  use SVD instead
   call ell_chol2eig( n, gg, a, work )
   r = 1.d0 / work(1)
   return
endif

if( ev(1) <= 0.d0 ) then
   write(0,*)'ell_rad_lower: lam_min <=0 ', ev(1)
   stop
endif

r = 1.d0 / sqrt( ev(1) )

return
end subroutine ell_rad_lower
