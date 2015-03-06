!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ell_pair_cover( algorithm, n, c1, gg1, c2, gg2, cc, gg )

!  Given a pair of ellipsoids, E1 and E2, determine a third ellipsoid
!  E which covers E1 and E2.
!  E is defined by:  (x-cc)^T * G * G^T * (x-cc) <= 1, where the array gg
!  contains G in packed format.  Similarly for E1 and E2. 

!  Algorithms to compute E are:
!     algorithm = 1 - spheroid (no shrinking)
!     algorithm = 2 - covariance (QL implementation)
!     algorithm = 3 - iterative
!     algorithm = 4 - spheroid (shrinking)
!     algorithm = 5 - covariance (Cholesky implementation)

implicit none

integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: algorithm, n
real(k_dp), intent(in)  :: c1(n), gg1((n*(n+1))/2), c2(n), gg2((n*(n+1))/2)
real(k_dp), intent(out) :: cc(n),  gg((n*(n+1))/2)

   if( algorithm == 1 ) then
      call ell_pair_cover_sph( n, c1, gg1, c2, gg2, .false., cc, gg )

   elseif( algorithm == 2 ) then
      call ell_pair_cover_cv_ql( n, c1, gg1, c2, gg2, cc, gg )

   elseif( algorithm == 3 ) then
      call ell_pair_cover_it( n, c1, gg1, c2, gg2, cc, gg )

   elseif( algorithm == 4 ) then
      call ell_pair_cover_sph( n, c1, gg1, c2, gg2, .true., cc, gg )
      
   elseif( algorithm == 5 ) then
      call ell_pair_cover_cv( n, c1, gg1, c2, gg2, cc, gg )     
   else
      write(0,*)'ell_pair_cover: called for invalid algorithm = ', algorithm
	  stop
   endif

   return

end subroutine ell_pair_cover
