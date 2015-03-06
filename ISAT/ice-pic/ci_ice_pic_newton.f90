!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ci_ice_pic_newton( r, h, p, r0, tau0, r_tol, r_acc, &
                              rg, tau, zg, Tg, rR, zR, TR, &
                              dzds, drds, dtauds, dzTRds, niters, &
                              dr_norm, resid, sv_rat, angles, info )
                              
! Perform Newton iterations to determine the pre-image point  zg
!    and its related properties.  
! Also determined are the PIC tangent vectors, drds and dtauds.

!  Input:
!    r:           reduced composition at reconstruction point {zr,zue}
!    h:           enthalpy (cgs)
!    p:           pressure (cgs)
!    r0:          estimated of reduced composition of generating point, rg
!    tau0:        estimated integration time
!    r_tol:       error tolerance, |rR - r| < r_tol
!    r_acc:       acceptable error in |rR - r|

! Output: 
!    rg:          reduced composition at pre-image point {zr,zue}
!    tau:         integration time
!    zg:          composition of pre-image point
!    Tg:          temperature at zg
!    zR:          mapping from zg
!    rR:          reduced composition at zR
!    TR:          temperature at zR
!    dzds:        PIC tangent vector dz/ds
!    drds:        PIC tangent vector dr/ds
!    dtauds:      PIC tangent vector d(tau)/ds
!    dzTRds:      d/ds[zR;TR]
!    niters:      number fo Newton iterations (counted as ODE integrations)
!    dr_norm:     residual |rR-r| 
!    resid:       residual predicted by min_norm
!    sv_rat:      ratio of singular values of matrix M
!    angles:      principal angles (theta_SB, theta_SC, theta_BC, theta_BP, theta_CP)
!                 (requires diagnostics = .true. - set below)

!    info      =  0, success with |rR-r| < r_tol
!              =  1, success with |rR-r| < r_acc
!              <  0, failure
!              = -1, atten < atten_min .and. residual not decreased
!              = -3, not converged in iter_max
!              = -4, failure in CEQ, T < T_low
!              = -5, failure in ci_ice_chem_map
!              = -6, hitting the boundary

!  Setting lu_PN=-1 (below) produces the file PN.op which can be post-processed using F_Newt.m

use ci_dat
use ci_dat8
use ci_utils
use ci_ice_cksubs
use ci_cem_recon
use ci_stats
implicit none
	
real(k_dp), intent(in)  :: r(nrc), h, p, r0(nrc), tau0, r_tol, r_acc

integer,    intent(out) :: niters, info
real(k_dp), intent(out) :: rg(nrc), tau, zg(ns), Tg, zR(ns), rR(nrc), TR, &
                           dzds(ns), drds(nrc), dzTRds(ns+1), dtauds, dr_norm, &
                            resid, sv_rat, angles(2+3*nrc)

! Parameters affecting performance
logical,    parameter :: diagnostics = .false.
integer,    parameter :: iter_min  = 0      ! min. no of iterations (>=0)
real(k_dp), parameter :: atten_min = 0.5d0  ! min. attenuation of Newton step
real(k_dp), parameter :: mn_fac    = 0.d0   ! min. norm tol: mn_tol = mn_fac * r_tol

! local arrays
integer, save :: lu_PN = -1  ! set < -1 for no diag. o/p
integer :: op_max = 1000 ! max. no. of calls for diag. o/p

integer   :: i, iter, iflag_CEQ, iflag_ODE, iflag_CEM, k_min, rankA, rankB, npim
               
real(k_dp), save :: calls = 0.d0
real(k_dp):: dr(nrc), t_scale, mn_tol, alpha, r_ref, &
             M(nrc, nrc+1), z(ns+1), Sg(ns), &
             stats(20), CEM_tan(ns+1,nrc+2), T_CE(ns,nrc), &
             S_norm, Sg0(ns), S_R(ns), atten, drt_scl(nrc+1), &
             r_neg_tol, vec(nrc+1), alpha_k, &
             frac_step, A_ODE(ns+1,ns+1), dr_dtau(nrc+1), &
             SR(ns), SRA(ns), SRA_err, S_rat
             
! Quantities saved from last iteration             
real(k_dp) :: rg_last(nrc), tau_last, zg_last(ns), Tg_last, zR_last(ns), &
              rR_last(nrc), TR_last, dzds_last(ns), drds_last(nrc), &
              dtauds_last, dr_norm_last, resid_last, vec_n, dr_dtau0(nrc+1), &
              gamma, beta_r, beta_v, &
              theta_SB(1), theta_SC(1), theta_BC(nrc), theta_BP(nrc), &
              theta_CP(nrc), svrat, T_PIM(ns,ns+1-nrc)

!stats
call routine_start(i_ci_ice_pic_newton)
              
calls  = calls + 1.d0

! initializing the returning variables
rg = 0.d0; tau = 0.d0; zg = 0.d0; Tg = 0.d0; zR = 0.d0; rR = 0.d0; TR = 0.d0
dzds = 0.d0; drds = 0.d0; dtauds = 0.d0; dr_norm = 0.d0; resid = 0.d0; angles = 0.d0
niters = 0
      
t_scale =  0.d0  ! initialize quantities that may be o/p
r_ref   =  0.d0
S_norm  =  0.d0
sv_rat  =  0.d0
resid   = -1.d0
dr_norm = -1.d0
k_min   = -2
Tg      = -1.d0
  
!==============================================================================

!  Starting point of Newton iterations: r0 = {zr,zue}.
      
mn_tol = mn_fac * r_tol

info      = -1   !  anticipate failure
iflag_ODE = 1
alpha     = 1.d0
atten     = 1.d0
dr_norm   = huge(0.d0)     
rg        = r0  !  initial estimates
tau       = max( tau0, 0.d0 )
      
call put_last  !  values on last accepted step
      
do iter = 1, iter_max  !  start of Newton iterations =============

!  Obtain zg = z^CE(rg)

   call ci_ceq( 0, rg, h, p, .false., z, Tg,  &
                   zg, Tg, stats, iflag_CEQ )
                   
!  ci_ice_ceq may fail  due to low/high temperature 
   if ( iflag_CEQ < 0 )  then 
      if( iflag_CEQ <= -18 ) then  !  T out of allowed range
         dice(16) = dice(16) + 1.d0
         !XXX write(0,*)'ICE_PIC: CEQ failure', iflag_CEQ, Tg !XXXX
         info = -4
         exit  !  failure due to T_CE out of range
      else
         call isat_abort('ci_ice_pic_newton', 1, &
               mess='CEQ failure, info = ', isv=iflag_CEQ )
      endif
   endif 
   
!  Note: if a perturbation is performed in CEQ, then B^T * zg may not equal rg.
!  However, rg should not be reset, so that repeat call with returned value of 
!  rg yields identical results
   
!  get CEM tangent vectors
   call ci_cem_tan(zg, h, Tg, p, thermo_ns, ns, nrc, &
                       BB, CEM_tan,  iflag_CEM )  
     
   if( iflag_CEM /= 1 ) call isat_abort('ci_ice_pic_newton', 2, &
                     mess='CEMTan failed' )
     
   T_CE = CEM_tan(1:ns,1:nrc)    
      
!  get Sg = S(zg)
   call ciS( zg, Tg, p, Sg )
   S_norm  = norm( Sg )
   
   if( S_norm <= 0.d0 ) call isat_abort('ci_ice_pic_newton', 3, mess='S_norm = 0' )
      
   if( iter == 1 ) Sg0 = Sg      
      
! integrating the ODE for R and A
    iflag_ODE = 0;
        
    call ci_ice_chem_map( 1, tau, p, zg, h, TR, &
                          zR, A_ODE, iflag_ODE ) 
       
    niters   = niters + 1   
    dice(25) = dice(25) + 1.d0 ! count of ODE integrations
    if ( iflag_ODE > 0 ) then
      dice(17)  = dice(17) + 1.d0  !  count of failures
      !write(0,*)'ICE_PIC: ODE failure', iflag_ODE, tau !XXXX
      info = -5
      exit	   
    endif  
       
!  residual
    rR      = matmul( BBT, zR )
    dr      = r - rR
    dr_norm = norm( dr )        ! residual
         
!  make decisions on if and how to continue iterations  
     
    if( dr_norm <= r_tol .and. iter > iter_min ) then
       info = 0   ! success: Newton has converged
       ! But need to continue to obtain tangent vector based on current solution.

    elseif( iter == iter_max ) then
       if( dr_norm <= r_acc ) then
          info = 1   ! success: Newton has converged acceptably
       else
          info = -3
          exit  !  failure: not converged in iter_max iterations
       endif
          
    elseif( dr_norm >= 0.9d0 * dr_norm_last ) then
       !  residual has not decreased
       
       if( dr_norm_last <= r_acc ) then
          !  take previous step, which was acceptable
          call get_last
          info = 1
          exit
       endif

       ! reject step: attenuate previous step and continue
       ! (Note: if here, iter > 1,  dr_dtau defined on previous step)
        atten  = 0.25d0 * atten
        if( atten < atten_min ) then
           info = -1
           exit !  failure
        endif
        call diagnostic_op( -1 )
        call get_last  
        call r_step( rg, tau, dr_dtau, atten )
        cycle
    endif
       
!  Newton step accepted: values from this accepted step stored below, 
!                       once drds etc. have been evaluated

    atten    = min( 2.d0 * atten, 1.d0 ) !attenuation factor for step

!  Form Newton matrix
    M(1:nrc,1:nrc) = matmul( BBT, matmul( A_ODE(1:ns,1:ns), T_CE ) )
       
    M(:,nrc+1) = matmul( BBT, matmul( A_ODE(1:ns,1:ns), Sg ) ) 

!  specify scaling factor for 
    r_ref   = maxval( rg ) 
    t_scale = r_ref / S_norm
       
    drt_scl        = r_ref  
    drt_scl(nrc+1) = t_scale 
       
!  Definitions: rg_new = rg + drg = r0 + Delta_rg, and similarly for tau.
!  Obtain minimum norm solution for Delta_rg.

    dr_dtau0(1:nrc) = rg - r0
    dr_dtau0(nrc+1) = tau - tau0
    dr = dr + matmul( M, dr_dtau0 )
   
!  Solve for Delta_rg
    call ci_ice_min_norm( nrc, nrc+1, M, dr, drt_scl, mn_tol, dr_dtau, resid, vec, sv_rat)
    
!  Express solution as drg
    call re_norm( dr_dtau(1:nrc), amolwt_n, beta_r )
    dr_dtau = dr_dtau - dr_dtau0
    
!  Normalize vec by r-components
   call re_norm( vec(1:nrc), amolwt_n, beta_v )
   
   vec_n   = norm( vec(1:nrc) )
   if( vec_n == 0.d0 ) call isat_abort('ci_ice_pic_newton', 4, &
                     mess='vec_n = 0' )
                     
!  Take the solution  rg_new=rg+drg+gamma*vec  which minimizes |rg_new - r0|
    vec_n   = norm( vec(1:nrc) )
    gamma   = -dot_product( rg+dr_dtau(1:nrc)-r0, vec(1:nrc) ) / vec_n
    dr_dtau = dr_dtau + gamma * vec / vec_n
    
    drds   = vec(1:nrc) / vec_n
    dtauds = vec(nrc+1) / vec_n 
    dzds   = matmul( T_CE, drds )
    dzTRds = matmul( A_ODE(1:ns+1,1:ns), dzds(1:ns) + Sg(1:ns)*dtauds )
      
    if( norm(dr_dtau) == 0.d0 ) info = 0  !  convergence - zero increment
    if( info == 0 ) exit  !  convergence
                  
!  Determine alpha,  such that rg + alpha*drg >= -r_neg_tol,
!                    where -r_neg_tol is a small negative number 
!  (i.e., it takes alpha Newton steps to reach the (extended) boundary)

    alpha     = huge(0.d0)
    r_neg_tol = min( 1.d-5 * r_ref, 1.d-2 * norm( dr_dtau(1:nrc) ) )
       
    do i = 1, nrc
       if( dr_dtau(i) < 0.d0 ) then
          alpha_k = -(rg(i)+r_neg_tol) / dr_dtau(i)
          alpha   = min( alpha, alpha_k )
       endif  
    end do
    
    if( dr_dtau(nrc+1) < 0.d0 )  then
       alpha_k = -tau / dr_dtau(nrc+1)
       alpha   = min( alpha, alpha_k )
     endif
              
!  Decide on fraction of Newton step -- or quit 
     frac_step = atten * min( 0.99d0 * alpha, 1.d0 )
     
     if( alpha >= 1.d0  .or.  &
         (1.d0-frac_step)*dr_norm < 0.9d0 * r_tol  .or. &
         (iter == 1  .and.  alpha > 0.5d0 ) ) then
         !  take step
     else
     
         if( dr_norm < r_acc ) then
            info = 1  !  accept this step
            exit      !  success
            
         elseif( dr_norm_last < r_acc ) then
            info = 1  !  accept previous step 
            call get_last
            exit      !  success
         endif

         info = -6
         exit  !  failure likely - hitting boundary
     endif
     
!  Store results from previous, accepted step
    call put_last  
    
!  Take Newton step
    call r_step( rg, tau, dr_dtau, frac_step )
    
    call diagnostic_op( 0 )
        
end do  !==========end of Newton iterations ==================
      
call diagnostic_op( info )

!stats
call routine_stop(i_ci_ice_pic_newton)

if( info /= 0 ) return
if( .not.diagnostics ) return

!  For diagnostics, compare: SR, A*Sg
   call ciS( zR, TR, p, SR )
   SRA = matmul( A_ODE(1:ns,1:ns), Sg )
   SRA_err = norm( SRA - SR ) / norm(SR)
   S_rat   = norm( SR ) / norm( Sg )
   
   write(1,'(1p,100e13.4)') Tg, TR, SRA_err, S_rat

!  For diagnostic purposes, examine angles between sub-spaces

   npim = ns+1-nrc

   call ci_ice_pim_tan( A_ODE, Sg, T_PIM, svrat )

   theta_BP = 0.d0
   theta_CP = 0.d0
   call principal_angles( ns, 1,   nrc,  Sg,   BB,    rankA, rankB, theta_SB )
   call principal_angles( ns, 1,   nrc,  Sg,   T_CE,  rankA, rankB, theta_SC )
   call principal_angles( ns, nrc, nrc,  BB,   T_CE,  rankA, rankB, theta_BC )
   call principal_angles( ns, nrc, npim, BB,   T_PIM, rankA, rankB, theta_BP )
   call principal_angles( ns, nrc, npim, T_CE, T_PIM, rankA, rankB, theta_CP )

   angles = (/ theta_SB, theta_SC, theta_BC, theta_BP, theta_CP /)

return
       
contains  !-------------------------------------------------

subroutine put_last
!  save results from last iteration

rg_last      = rg
tau_last     = tau 
zg_last      = zg 
Tg_last      = Tg 
zR_last      = zR 
rR_last      = rR 
TR_last      = TR
dzds_last    = dzds 
drds_last    = drds 
dtauds_last  = dtauds 
dr_norm_last = dr_norm 
resid_last   = resid

return
end subroutine put_last

subroutine get_last  !------------------------------------
!  get results from last iteration

rg      = rg_last
tau     = tau_last 
zg      = zg_last 
Tg      = Tg_last 
zR      = zR_last 
rR      = rR_last 
TR      = TR_last
dzds    = dzds_last 
drds    = drds_last 
dtauds  = dtauds_last 
dr_norm = dr_norm_last 
resid   = resid_last

return
end subroutine get_last
   
subroutine r_step( rg, tau, dr_dtau, atten )  !------------
       
!  increment rg and tau based on attenuated step dr_dtau

real(k_dp), intent(inout) :: rg(nrc), tau, dr_dtau(nrc+1)
real(k_dp), intent(in)    :: atten
       
rg       = rg  + dr_dtau(1:nrc) * atten
tau      = tau + dr_dtau(nrc+1) * atten
tau      = max( tau, 0.d0 )
       
do i = 1, nrc
    rg(i) = max( rg(i), 0.d0 )
end do

rg = rg / dot_product( rg, amolwt_n )
       
return
end subroutine r_step
       
subroutine diagnostic_op( flag )  !--------------------------
    
!  diagnostic output
     
integer, intent(in) :: flag
   
if( lu_PN == -1 ) then
   call isat_lu( lu_PN )
   open( lu_PN, file='PN.op' )
endif
       
if( lu_PN < 0  .or.  calls > op_max ) return
        
call ciS( zR, TR, p, S_R )
                 
write(lu_PN,'(4i6,1p,14e13.4)') nint(calls), iter, flag, 0, &
                 dr_norm, sv_rat, alpha, tau, Tg, TR, &
                  r_ref, S_norm, norm_dot(Sg,Sg0), &
                  norm_dot(Sg,S_R), atten, r_tol, r_acc, resid
                      
return
        
end subroutine diagnostic_op

end subroutine ci_ice_pic_newton
