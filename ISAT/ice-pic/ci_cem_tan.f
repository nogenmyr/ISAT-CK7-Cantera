!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

      subroutine ci_cem_tan(z, h, T, p, thermo, ns, nrc, B,
     1     A_CEM, info)

!  This subroutine calculates the tangent space of the CE Manifold at fixed enthalpy
!  and pressure. 
!  Note: This subroutine can also be used to compute the tangent space of the CE 
!  manifold at the edge face.  In that case, ns is the number of non-zero species 
!  at the boundary and thermo are the corresponding thermo data for only those species.
!
!  Input:
!         z:     species specific moles on the CEM  
!         h:     enthalpy 
!         T:     temperature on CEM
!         p:     pressure
!         thermo, ns*15: thermo data
!         ns:    number of species in the whole system
!         nrc:   number of constraints (enthalpy is not included)
!         B:     constraint matrix: (ns*nrc) (enthalpy is not included)

!  Output : 
!         A_CEM, (ns+1)* (nrc+2): Tangent vectors => d{z, T}^CE/d{r,h,p}
!                  A_CEM(1:ns, 1:nrc  ) : dz/dr
!                  A_CEM(1:ns, nrc+1  ) : dz/dh
!                  A_CEM(1:ns, nrc+2  ) : dz/dp
!                  A_CEM(ns+1, 1:nrc+2) : dT/d{r,h,p}
!         info = 1  successful
!               -1  singular    

      use ci_utils
      use isat_abort_m 
      use ci_dat6, only: k_dp, gascon
      use ci_stats

      implicit none
      
      real(k_dp), intent(in)  :: z(ns), h, p, B(ns,nrc),
     1     T, thermo(ns, 15) 
      integer, intent(in)     :: ns, nrc
       
      real(k_dp), intent(out) :: A_CEM(ns+1, nrc+2)
      integer, intent(out)    :: info
      
!     Local variables --------------
      integer          :: i, j, ipiv(nrc+1), info_svd, info_sv

      real(kind(1.d0)) :: zlim=1.d-30 !  lower limit on z

      real(kind(1.d0)) :: zf(ns), r(nrc), rf(nrc), M(ns,nrc), rtr, 
     1     zBr(ns), U(ns,nrc), Sigma(nrc), VT(1,1), 
     2     UTB(nrc,nrc), T_T(nrc,ns), work(ns*ns+10*(ns+nrc))

      real(kind(1.d0)) :: dzdr(ns,nrc), dzdp(ns), dzdT(ns), zTh, Nb

      real(kind(1.d0)) :: hort(ns), h_unnorm(ns), R_gas_const, cpor(ns),
     1     cp(ns), cp_mix, cp_corr, tmp2(ns,ns), tmp(ns)

      
      call routine_start(i_ci_cem_tan)

!     kb_exact not implemented yet
!     NOTE: Exact derivatives at the boundary not imlemented yet.
     
!     Initialization       
      A_CEM = 0.d0
      info    = 1
      r = matmul(transpose(B), z)

!     Filter out zero species
      do i = 1, ns
         zf(i) = max( z(i), zlim )
      enddo
      rf = matmul(transpose(B),zf)

!     Compute the tangent space of CEM at fixed T and p

!     Compute dzdr
!     Form M = M/r = (z - Z*B*r/|r|^2) + ZB/r

      rtr = dot_product(r,r)
      zBr = matmul(diagmul(z,B),r)/rtr
      do i=1,ns
         do j =1,nrc
            M(i,j) = z(i) - zBr(i) + (zf(i)/rf(j))*B(i,j)
         enddo
      enddo

!     compute SVD of M = U * Sigma * VT
      call dgesvd( 'S', 'N', ns, nrc, M, ns, Sigma, 
     1     U, ns, VT, nrc, work, ns*ns+10*(ns+nrc), info_svd )

      if( info_svd /= 0 ) call isat_abort( 'ci_cem_tan', 2, mess =
     1     'dgesvd failed: info = ', isv = info_svd )


!     solve using lu factorisation
!     solve: (U^T * B) * dzdr^T = U^T
      UTB = matmul( transpose(U), B )
      T_T = transpose( U )
      
      call dgesv( nrc, ns, UTB, nrc, ipiv, T_T, nrc, info_sv )
      
      if( info_sv /= 0 ) call isat_abort( 'ci_cem_tan', 3, mess =
     1     'dgesv failed: info = ', isv = info_sv )

      dzdr = transpose(T_T)

!     Compute dzdT
      call ceq_h(ns, T, thermo, hort)
      zTh = dot_product(z,hort)
      tmp = (zTh*z + diagmul(z,hort) - zTh*zBr)/T
      dzdT = tmp - matmul(dzdr, matmul(transpose(B), tmp))

!     Compute dzdp
      Nb = sum(z)
      tmp = (Nb*zBr - (Nb + 1)*z)/p
      dzdp = tmp - matmul(dzdr, matmul(transpose(B), tmp))

!     compute the tangent space of CEM at fixed h and p.
!     unnormalized h
      R_gas_const =  gascon
      h_unnorm = hort * R_gas_const * T

      call ceq_cp(ns,T,thermo,cpor)

      cp = cpor *R_gas_const
      cp_mix= dot_product(cp, z)
      cp_corr = dot_product(h_unnorm, dzdT)
       
      do i=1, ns
         do j=1,ns
            tmp2(i,j)= -1.d0/(cp_mix +cp_corr) * dzdT(i)* h_unnorm(j)
         end do
         tmp2(i,i)=tmp2(i,i) + 1.d0
      end do

!     dz/dr       
      A_CEM(1:ns,1:nrc) = matmul(tmp2, dzdr)

!     dz/dh             
      A_CEM(1:ns,nrc+1) = 1.d0/(cp_mix +cp_corr) * dzdT
       
!     dz/dp
      A_CEM(1:ns,nrc+2) = matmul(tmp2, dzdp)

!     dT/dr
      A_CEM(ns+1,1:nrc) = -1.d0/cp_mix * 
     1     matmul(h_unnorm,A_CEM(1:ns,1:nrc))

!     dT/dh
      A_CEM(ns+1,nrc+1) = 1.d0/cp_mix * 
     1     (1.d0 - dot_product(h_unnorm, A_CEM(1:ns,nrc+1)))

!     dT/dp
      A_CEM(ns+1,nrc+2) = -1.d0/cp_mix * 
     1     dot_product(h_unnorm, A_CEM(1:ns,nrc+2))

      call routine_stop(i_ci_cem_tan)

      return
      end subroutine ci_cem_tan
