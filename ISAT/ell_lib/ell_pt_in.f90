!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ell_pt_in( n, c, gg, p, in )

!  Return  in=.true. if the point p is in the ellipsoid E.
!  E is given by { x | norm(G^T * (x-c) ) <=1 ),  where G is an 
!  n x n lower triangular matrix.  The array gg contains 
!  the matrix G in packed format.

!  Note: select below method to be used.

!  S.B. Pope  6/12/04

implicit none

integer, parameter     :: k_dp = kind(1.d0)
integer, intent(in)    :: n
real(k_dp), intent(in) :: c(n), gg((n*(n+1))/2), p(n)
logical, intent(out)   :: in

integer    :: j, lm1, kst, method=2
real(k_dp) :: y(n), ysq1, ysq, yj

in = .false.  !  assume p outside E until found otherwise

! Method 1:  directly evaluate  y = G^T * [p-c] and test y^T * y <= 1.
!            (Simple, but not optimal.)

if( method == 1 ) then

!  Transform to y-space: y = G^T * [x-c]
   y = p-c   
   call dtpmv( 'L', 'T', 'N', n, gg, y, 1 )  !  yp = G^T * [p-c]
   ysq1 = sum(y*y) 
   if( ysq1 <= 1.d0 ) in = .true.

else

!  Method 2:  successively evaluate components of y = G^T * [p-c]; accumulate sum; test
!             (More involved, but possibly quicker)

   y   = p-c   
   ysq = 0.d0
   kst = (n*(n+1))/2
   do j = n, 1, -1  !  loop backwards over columns of G
      lm1 = n-j     !  length of column vector - 1
      yj  = sum( gg(kst:kst+lm1) * y(j:n) )
      ysq = ysq + yj*yj

	  if( ysq > 1.d0 ) return

      kst = kst - lm1 - 2
   end do

   in = .true.

endif

return
end subroutine ell_pt_in
