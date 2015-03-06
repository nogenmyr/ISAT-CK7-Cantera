!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

      subroutine ci_ice_chem_map( map, tau, press, phi, h, TRphi,
     1                            Rphi, A, iflag_ODE)
        
! This subroutine performs the numerical integration of the ODE system.

! Input: 
!    map =1, sensitivity is to be returned
!    map =0, no sensitivities required
!    tau,    integration time
!    press:  pressure (cmgs)
!    phi :   initial species specific moles
!    h   :   enthalpy (total)
!    TRphi:  tempeature after reaction

! Output:
!    Rphi:   species specific moles after reaction
!    A:      sensitivity matrix
!    iflag_ODE:  =0 successful integration
!                >0 unsucessful ODE itegration
     
      use ci_dat
      use ci_dat8
      use ci_ice_cksubs
      use ci_stats
      use ci_utils

      implicit none
	
      integer, intent(in)    :: map
      real(k_dp), intent(in) :: tau, press, phi(ns), h
      real(k_dp), intent(out):: TRphi, Rphi(ns), A(ns+1,ns+1)
      integer, intent(out)   :: iflag_ODE

!local arrays
      integer    :: check, info, i, j
      real(k_dp) :: y(ns), T, z0(ns+1), zt(ns+1), amwA(ns)

! ------- Update ci_stats
      if(map==0) then
         call routine_start(i0_ci_ice_chem_map)
      else
         call routine_start(i1_ci_ice_chem_map)
      endif

! Initialization
       iflag_ODE=0
       TRphi=0.d0
       Rphi=0.d0
       A=0.d0
       
       check = 1
       call temphz( h, phi, t, check, info )
        
        if ( info /= 0 ) then  ! temperature out of range
           iflag_ODE = 1
          return
        endif  
           
        z0(1:ns) = phi
        z0(ns+1) = h
        
        call ci_chem_one_step( map, phi, T, press, tau, zt(1:ns),
     1                        TRphi, A, iflag_ODE )
     
 !  call DDASSAC only if ci_chem_one_step fails
        if( iflag_ODE /= 0 ) then
           call ci_chem_ddassac(map, z0,press, tau, TRphi, zt, A, 
     1          iflag_ODE)
        
           if (iflag_ODE >0 ) then
             ! failure
              if(map==0) then
                 call routine_failed(i0_ci_ice_chem_map)
              else
                 call routine_failed(i1_ci_ice_chem_map)
              endif
              return
           endif 
           
        endif

       Rphi(1:ns) = zt(1:ns)

       !XXX perform corrections to satisfy normalization
       Rphi(1:ns) = Rphi(1:ns) / dot_product( Rphi(1:ns), amolwt )
       amwA = amolwt - matmul( amolwt, A(1:ns,1:ns) )
       do i = 1, ns
          do j = 1, ns
             A(i,j) = A(i,j) + amolwt(i) * amwA(j) / norm( amolwt ) **2
          end do
       end do
       !XXX end of normalization
      
       ! Update ci_stats
       if(map==0) then
          call routine_stop(i0_ci_ice_chem_map)
       else
          call routine_stop(i1_ci_ice_chem_map)
       endif
        
       end subroutine ci_ice_chem_map
