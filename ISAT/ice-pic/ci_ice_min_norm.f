!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!


      subroutine ci_ice_min_norm( m, n, A, b, z, tol,
     1                                x, r_norm, vec, sv_rat )
     
! Given the linear system A * y = b (where A is m x n, m<=n), 
! return an approximate solution x, such that
! A * x = b + r, |r| = tol, and |x/z| is minimized, where z
! are positive, specified scaling factors, and tol is a non-negative tolerance.
! Also returned are:
!   r_norm = |r| - the norm of the residual 
!   vec          - an n-vector, which (for m<n) is in the null space of A,
!   sv_rat       - the ratio of singular values of the scaled A-matrix, S( min(m,n) ) / S(1).
! If |b| < tol, x=0 is returned (in which case |r| < tol).
! If A is ill-conditioned, the residual may exceed tol.

      use isat_abort_m
      use ci_dat8, only: k_dp
      use ci_utils

      implicit none
	
      integer, intent(in)     :: m, n
      real(k_dp), intent(in)  :: A(m,n), b(m), z(n), tol
      real(k_dp), intent(out) :: x(n), r_norm, vec(n), sv_rat
      
! Local arrays
      real(k_dp) :: work(m*n+10*(m+n)), U(m,m), VT(n,n), 
     1             Sig(n), A_t(m,n), bh(m), xt(m), b_norm,
     2             r(n), S(m), SI(m), S0
      integer    :: i, info, its
      
!  Check input
       if( minval(z) <= 0.d0 ) then
          call isat_abort('ci_ice_min_norm', 1,
     1            mess='z must be positive: min(z)= ', rsv = minval(z) )
       endif
     
        if( tol < 0.d0 ) call isat_abort('ci_ice__min_norm', 2,
     1            mess='tol must be non-negative: tol= ', rsv = tol )
     
        if( n < m ) call isat_abort('ci_ice__min_norm', 3,
     1            mess='n less than m: m, n = ', ivar = (/m, n/) ) 
          
!  Form scaled matrix A_t
       do i = 1, n
          A_t(:,i) = A(:,i) * z(i)
       end do
       
!    Form SVD of matrix:  A_t = U*S*VT      

       call dgesvd( 'A', 'A', m, n, A_t, m, Sig, U, m, VT, n,
     $                   work, m*n+10*(m+n), info )
      
      if (info /=0 ) call isat_abort('ci_ice_min_norm', 4, mess=
     1     'SVD failed, info = ', isv = info )
     
       if ( Sig(1) <= 0.d0 ) call isat_abort('ci_ice_min_norm', 5, 
     1     mess='bad singular value, Sig(1) = ', rsv = Sig(1) )
     
       S      = Sig(1:m)
       sv_rat = S(m)/S(1)
       vec    = VT(n,:) * z
       vec    = vec / norm( vec )
       
       bh     = matmul( b, U )
       b_norm = norm( bh )
       
       if( b_norm <= tol ) then
          x      = 0.d0
          r_norm = b_norm
          return
       endif
       
       call find_S0
          
       SI = S/(S0**2 + S**2)     ! modified inverse of S
       xt = SI * bh
       x  = matmul( xt, VT(1:m,1:n) ) * z ! solution
       
       if( .false. ) then ! optional output
          r  = b - matmul(A,x)
          b_norm = norm( r )  ! r_norm
          write(0,'(20x,i4,1p,10e13.4)') its, b_norm, b_norm/tol-1.d0,
     1        r_norm/tol-1.d0
       endif
       
      return
      
      contains  !------------------------------------------------------
      
      subroutine find_S0
      
 !  Determine S0 which yields r_norm = tol.
      
      integer,    parameter :: its_max  = 40
      integer,    save      :: its_most = 0
      real(k_dp), parameter :: sv_min   = 1.d-12
      
      real(k_dp) :: rup, sup, rlow, slow, S02, drdS02, S0_up
      
      !  minimum allowed value of S0
      S0 = sv_min * S(1)
      r_norm = norm( bh * S0**2/(S0**2 + S**2) )
      if( r_norm >= tol ) return !  accept S0

      S0 = sqrt( S(1)*S(m)*tol / b_norm )  ! initial guess
      ! upper bound
      S0_up = S(m) * sqrt( tol / (m*maxval( abs(bh) ) - tol ) )
      S0 = min( S0, S0_up )
      r_norm = norm( bh * S0**2/(S0**2 + S**2) )
      
      if( r_norm > tol ) then  !  bracket root
         rup = r_norm
         sup = S0
         do while ( r_norm > tol )  !  bound below
            S0 = S0 / 16.d0
            r_norm = norm( bh * S0**2/(S0**2 + S**2) )
         end do
         rlow = r_norm
         slow = S0     
      else
         rlow = r_norm
         slow = S0
         do while ( r_norm <= tol )  !  bound above
            S0 = S0 * 16.d0
            r_norm = norm( bh * S0**2/(S0**2 + S**2) )
         end do
         rup = r_norm
         sup = S0   
      endif
      
      its = 0  !  iterate to convergence
      do while ( abs(r_norm-tol) > 1.e-6 * tol )
         its    = its + 1
         drdS02 = sum( (bh*S0*S)**2/(S0**2 + S**2)**3 )/r_norm
         S02    = S0**2 + (tol-r_norm)/drdS02
         
         if( S02 > slow**2  .and.  S02 < sup**2 ) then
            S0 = sqrt(S02)  !  Newton applied to S0**2
         else
            S0 = .5d0*(sup+slow)  !  bisection
         endif
         
         r_norm = norm( bh * S0**2/(S0**2 + S**2) )
         
         if( r_norm > tol ) then
            rup = r_norm
            sup = S0
         else
            rlow = r_norm
            slow = S0
         endif
         
         if( its > its_max ) then
            call isat_abort('ci_ice_min_norm', 6, 
     1         mess= 'convergence failure in find_S0' )
         endif
 
      end do
      
      if( .false.  .and.  its > its_most ) then !optional diagnostic o/p
         its_most = its
         write(0,'(a,i4,1p,20e10.2)') 'find_S0 ', its, tol, 
     1        r_norm/tol-1.d0, S(m)/S(1), S0/S(1)
      endif
            
      return
      end subroutine find_S0
      
      end subroutine ci_ice_min_norm
