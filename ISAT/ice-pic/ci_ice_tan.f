!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

      subroutine ci_ice_tan(z_g, T_g, z_ICE, T_ICE, h, p,       
     1     A_ODE, B, k_facet, r_normal, A_ICE, info)
     
! This subroutine calculates the tangent space of the ICE manifold.
! The generating point may, or may not, be on the boundary.
   
!--------------------------------------------------------------------------
! Input:

!   z_g(ns)    : generating point composition
!   T_g        : generating point temperature 
!   z_ICE(ns)  : ICE manifold point
!   T_ICE      : ICE point temperature       
!   h          : total enthalpy 
!   p          : pressure
!   A_ODE (ns+1, ns+1): sensitivity matrix from ODE integration   
!   B          : constraint matrix: (ns*nr) B=[CS  CEu]  
!   k_facet    : = 0, interior
!                > 0, index of facet
!   r_normal   : normal on facet k

!  Output : 
!         A_ICE, (ns+1)* (nrc+2): Tangent vectors => d{z, T}^ICE/d{r,h,p}
!                  A_ICE(1:ns, 1:nrc  ) : dz/dr
!                  A_ICE(1:ns, nrc+1  ) : dz/dh
!                  A_ICE(1:ns, nrc+2  ) : dz/dp
!                  A_ICE(ns+1, 1:nrc+2) : dT/d{r,h,p}
!         info = 1  successful
!               -1  singular    

      use ci_ice_cksubs
      use ci_utils
      use ci_dat6, only: gascon, thermo_ns
      use ci_stats

      implicit none

      integer, intent(out) :: info
      integer, intent(in)  :: k_facet
      real(kind(1.d0)), intent(in) :: z_g(ns), T_g, z_ICE(ns), T_ICE,
     1     h, p, A_ODE(ns+1,ns+1), B(ns, nrc), r_normal(nrc)

      real(kind(1.d0)), intent(out) :: A_ICE(ns+1, nrc+2)
      
!---------local variables --------------
      integer          :: i, ipiv(nrc+1)

!     S.n < Sn_tol, then S is deemed to be parallel to the facet
      real(kind(1.d0)), parameter :: Sn_tol = 1.d-3

      real(kind(1.d0)) :: A_CEM(ns+1,nrc+2), C(ns, nrc),  Cn(ns),
     1     Sigma(nrc), U(ns,nrc), VT(1,1), work(ns*ns+10*(ns+nrc)),
     2     UTB(nrc,nrc), T_T(nrc,ns), S_g(ns), Sr_g(nrc), Sn

      real(kind(1.d0)):: hort(ns), h_unnorm(ns), R_gas_const, cpor(ns),
     1     cp(ns), cp_mix, tmp(ns), hs_n, r_g(nrc)

      call routine_start(i_ci_ice_tan)
      
!     determine A_CEM - tangent vectors to the CEM at z_g	 
      A_CEM = 0.d0
      call ci_cem_tan(z_g(1:ns), h, T_g, p, thermo_ns, ns,
     1     nrc, B, A_CEM, info) 

!---- the generating point is in the interior ------------

!     form C = A_ODE * A_CEM
      C = matmul( A_ODE(1:ns,1:ns), A_CEM(1:ns,1:nrc) )

!     obtain S_g: reaction rate at z_g
      call ciS( z_g, T_g, p, S_g )

!     Compute BBT*S.n
      Sr_g = matmul(BBT, S_g)
      Sn = dot_product( -r_normal, Sr_g )/norm(Sr_g)

!     Replace C(:,k_facet) with S only if S is not parallel to facet
      if( k_facet /= 0 .and. Sn > Sn_tol ) then
!---- the generating point is on the boundary ------------

!     write r_g for gradient test at the boundary
         r_g = matmul(transpose(B),z_g)

!     form modified matrix C using facet normal n
         Cn = matmul(C,r_normal) ! Cn = A*T*n
         S_g = matmul(A_ODE(1:ns,1:ns), S_g)/norm(S_g)  ! A*Sg/|Sg|
         S_g = S_g/norm(S_g)
         do i = 1,nrc
            C(:,i) = C(:,i) - ( Cn - S_g )*r_normal(i)
         enddo

      endif

!     compute SVD of C = U * Sigma * VT
      call dgesvd( 'S', 'N', ns, nrc, C, ns, Sigma, 
     1     U, ns, VT, nrc, work, ns*ns+10*(ns+nrc), info )
      
      if( info /= 0 ) call isat_abort( 'ci_ice_ICE_Tan', 1, 
     1     mess = 'dgesvd failed: info = ', isv = info )

!     solve using lu factorisation
!     solve: (U^T * B) * A_ICE^T = U^T
      UTB = matmul( transpose( U ), B )
      T_T = transpose( U )

      call dgesv( nrc, ns, UTB, nrc, ipiv, T_T, nrc, info )
      
      if( info /= 0 ) call isat_abort( 'ci_ice_ICE_Tan', 2, mess =
     1     'dgesv failed: info = ', isv = info )
      
      A_ICE(1:ns,1:nrc) = transpose( T_T )

! Enthalpy
      tmp = matmul(A_ODE(1:ns,1:ns), A_CEM(1:ns,nrc+1))
     1     + A_ODE(1:ns,ns+1)
      A_ICE(1:ns,nrc+1) = tmp - matmul(A_ICE(1:ns,1:nrc), 
     1     matmul(transpose(B(1:ns,1:nrc)), tmp))

! Pressure
! NOTE: neglecting dR/dp.  A_ODE(1:ns,ns+2) not present.
      tmp = matmul(A_ODE(1:ns,1:ns), A_CEM(1:ns,nrc+2))
      A_ICE(1:ns,nrc+2) = tmp - matmul(A_ICE(1:ns,1:nrc), 
     1     matmul(transpose(B(1:ns,1:nrc)), tmp))

! Last row  dT/d{r,h,p}
      R_gas_const =  gascon
      call ceq_h( ns, T_ICE, thermo_ns, hort )
      h_unnorm = hort * R_gas_const * T_ICE

      call ceq_cp( ns, T_ICE, thermo_ns, cpor )
      cp = cpor * R_gas_const
      cp_mix= dot_product(cp, z_ICE)

!     dT/dr
      A_ICE(ns+1,1:nrc) = -1.d0/cp_mix * 
     1     matmul(h_unnorm,A_ICE(1:ns,1:nrc))

!     dT/dh
      A_ICE(ns+1,nrc+1) = 1.d0/cp_mix * 
     1     (1.d0 - dot_product(h_unnorm, A_ICE(1:ns,nrc+1)))

!     dT/dp
      A_ICE(ns+1,nrc+2) = -1.d0/cp_mix * 
     1     dot_product(h_unnorm, A_ICE(1:ns,nrc+2))

      call routine_stop(i_ci_ice_tan)

      return
      end subroutine ci_ice_tan
