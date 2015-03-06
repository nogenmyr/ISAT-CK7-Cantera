!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ellu_radii( n, g, r_in, r_out )

!  g contains  the lower triangular matrix G 
!  which defines the ellipsoid E = {x | |G^T * x | <=1 }.
!  This routine returns in  r_in  and  r_out  the radii of the 
!  inscribed and circumscribed balls.  These are the inverses of the
!  larges and smallest singular vaues of G.

implicit none
integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: n
real(k_dp), intent(in)  :: g(n,n)
real(k_dp), intent(out) :: r_in, r_out

real(k_dp)   :: a(n,n), sv(n), u(n,n), vt(n,n), work(10*n*n+20*n)
integer      :: lwork, info
character(1) :: jobu, jobvt

lwork = 10*n*n+20*n

if( n > 1 ) then
!  form singular values
   a = g
   jobu  = 'A'  !XXX 'N' should be sufficient, but leads to unexplained error
   jobvt = 'N'

   call dgesvd( jobu, jobvt, n, n, a, n, sv, u, n, vt, n, work, lwork, info )

   if( info /= 0 ) then
      write(0,*)'ellu_radii: info, n, lwork, work(1) = ', info, n, lwork, work(1)
      write(0,*)'sv = '
      write(0,'(1p,5e13.4)') sv
      stop
   endif

elseif( n ==1 ) then

   sv(1) = g(1,1)

else
   write(0,*)'ellu_radii: invalid n = ', n
   stop
endif

if( sv(n) <= 0.d0 ) then
   write(0,*)'ellu_radii: bad singular value ', sv
   stop
endif

r_in  = 1.d0 / sv(1)
r_out = 1.d0 / sv(n)

return

return
end subroutine ellu_radii
