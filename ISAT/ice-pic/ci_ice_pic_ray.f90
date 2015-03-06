!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ci_ice_pic_Ray( r, h, p, r0, tau0, zg0, Tg0, zR0, TR0, drds0, dtauds0, sb, ds0, npts)
 
! Diagnostic routine to evaluate quantities used in ci_ice_pic_newton along a ray.
! Output on ice_pic_ray.op, post-process with pic_ray.m

! The ray is rg(s) = r0 + drds0 * s,  taug(s) = tau0 + dtauds0 * s, for 0 <= s <= sb.
                              
!  Input:
!    r:           reduced composition at reconstruction point {zr,zue}
!    h:           enthalpy (cgs)
!    p:           pressure (cgs)
!    r0:          composition of initial generating point, rg(0)
!    tau0:        integration time of initial generating point, taug(0)
!    zg0:
!    Tg0
!    zR0
!    TR0
!    drds0:       tangent vector along ray
!    dtauds0:     d(tau)/ds along ray
!    sb:          upper limit on s
!    ds0:         = 0.d0 - use linear spacing;  
!                 >0 - use log spacing with first two points being {0, ds0/sb} 
!    npts:        number of points to generate


use ci_dat
use ci_dat8
use ci_utils
use ci_ice_cksubs
use ci_cem_recon
implicit none
	
real(k_dp), intent(in)  :: r(nrc), h, p, r0(nrc), tau0, zg0(ns), Tg0, zR0(ns), &
                           TR0, drds0(nrc), dtauds0, sb, ds0
integer,    intent(in)  :: npts

!//////////////////////////////////////////////////////////////////////////////////////
real(k_dp) :: rg(nrc), tau, zg(ns), Tg, zR(ns), rR(nrc), TR, s(npts), amp, &
                           dzds(ns), drds(nrc), dzTRds(ns+1), dtauds, dr_norm, &
                            resid, sv_rat, angles(2+3*nrc)


integer   :: i, ipt, iflag_CEQ, iflag_ODE, iflag_CEM, k_min, rankA, rankB, npim, lu
               
real(k_dp):: dr(nrc), t_scale, mn_tol, t_ref, alpha, r_ref, &
             M(nrc, nrc+1), z(ns+1), Sg(ns), &
             stats(20), CEM_tan(ns+1,nrc+2), T_CE(ns,nrc), &
             S_norm, Sg0(ns), S_R(ns), atten, drt_scl(nrc+1), &
             r_neg_tol, TR_min=huge(0.d0), vec(nrc+1), alpha_k, &
             frac_step, A_ODE(ns+1,ns+1), dr_dtau(nrc+1), zero_r, &
             zero_rr, sva(ns), svm(nrc), SR(ns), SRA(ns), SRA_err, S_rat, &
             diffs(4)
             
             
! Quantities saved from last iteration             
real(k_dp) :: rg_last(nrc), tau_last, zg_last(ns), Tg_last, zR_last(ns), &
              rR_last(nrc), TR_last, dzds_last(ns), drds_last(nrc), &
              dtauds_last, dr_norm_last, resid_last, vec_n, dr_dtau0(nrc+1), &
              gamma, beta_r, beta_v, beta_drR, beta_min=0.d0, &
              theta_SB(1), theta_SC(1), theta_BC(nrc), theta_BP(nrc), &
              theta_CP(nrc), svrat, T_PIM(ns,ns+1-nrc), angle(2)
              
!/////////////////////////////////////////////////////////////////////////////////

call isat_lu(lu)
open( lu, file = 'ice_pic_ray.op')
 
!  specify s(1:npts)
if( ds0 == 0.d0 ) then
   do i = 1, npts
      s(i) = (i-1)*sb/float(npts-1)
   end do
   
elseif( ds0 > 0.d0 ) then
   s(1) = 0.d0
   amp  = ds0**(-1.d0/(npts-1))
   do i = 2, npts
      s(i) = ds0 * amp**(i-1)
   end do
   s = s *(sb/s(npts))

else
   call isat_abort('ci_ice_pic_Ray',1,mess='bad ds0 = ', rsv=ds0 )
endif

drds_last = drds0
      
do ipt = 1, npts  !  loop over points
   rg  = r0   + drds0   * s(ipt)
   tau = tau0 + dtauds0 * s(ipt)

!  Obtain zg = z^CE(rg)
   call ci_ceq( 0, rg, h, p, .false., z, Tg,  &
                   zg, Tg, stats, iflag_CEQ )
                   
   if ( iflag_CEQ < 0 )  exit

   rg = matmul( BBT, zg ) ! re-compute rg (because it may have been perturbed in CEQ) 
   diffs(1) = abs(Tg-Tg0)/Tg0
   diffs(2) = norm( zg - zg0 )/norm(zg0)
   
!  get CEM tangent vectors
   call ci_cem_tan(zg, h, Tg, p, thermo_ns, ns, nrc, &
                      BB, CEM_tan,  iflag_CEM )  
     
   if( iflag_CEM /= 1 ) exit
     
   T_CE = CEM_tan(1:ns,1:nrc)    
      
!  get Sg = S(zg)
   call ciS( zg, Tg, p, Sg )
   S_norm = norm(Sg)
            
! integrating the ODE for R and A
    iflag_ODE = 0;
        
    call ci_ice_chem_map( 1, tau, p, zg, h, TR, &
                          zR, A_ODE, iflag_ODE ) 
       
    if ( iflag_ODE > 0 ) exit
    
    diffs(3) = abs(TR-TR0)/TR0
    diffs(4) = norm( zR - zR0 )/norm(zR0)
    
!  residual
    rR      = matmul( BBT, zR )
    dr      = r - rR
    dr_norm = norm( dr )        ! residual
    

!  Form Newton matrix
    M(1:nrc,1:nrc) = matmul( BBT, matmul( A_ODE(1:ns,1:ns), T_CE ) )
       
    M(:,nrc+1) = matmul( BBT, matmul( A_ODE(1:ns,1:ns), Sg ) ) 
    
!  specify scaling factor for 
    r_ref   = maxval( rg ) 
    t_scale = r_ref / max( S_norm, 1.d-20*tau )  
    drt_scl        = r_ref
    drt_scl(nrc+1) = t_scale 
       
!  Solve:  M * dr_dtau = dr;  M * vec = 0
    mn_tol = 0.d0
    call ci_ice_min_norm( nrc, nrc+1, M, dr, drt_scl, mn_tol, dr_dtau, resid, vec, sv_rat)
    
!  Normalize vec by r-components
   call re_norm( vec(1:nrc), amolwt_n, beta_v )
   
   vec_n   = norm( vec(1:nrc) )
   if( vec_n == 0.d0 ) call isat_abort('ci_ice_pic_Ray', 2, &
                     mess='vec_n = 0' )
                     
   if( dot_product( vec(1:nrc), drds_last ) < 0.d0 ) vec = -vec
   drds = vec(1:nrc) / vec_n
   call principal_angles( nrc, 1, 1, drds, drds0,     rankA, rankB, angle(1) )
   call principal_angles( nrc, 1, 1, drds, drds_last, rankA, rankB, angle(2) )
   drds_last = drds
        
   write(lu,'(i5,1p,100e25.15)') ipt, s(ipt), tau, Tg, TR, dr_norm, resid, sv_rat, angle(1:2), &
                                 diffs(1:4), norm(dr_dtau(1:nrc)),  dr_dtau(nrc+1)               
end do   

end subroutine ci_ice_pic_Ray
