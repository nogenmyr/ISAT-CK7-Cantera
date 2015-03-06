!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ci_ice_recon_pic( r, h, p, r_tol, r_acc, Sg_min, npt, np_acc, npmx,  &
                       np, ipt_best, z_CE, rg, zg, Tg, taug, rR, zR, TR, &
                       rg_ex, zR_ex, TR_ex, Sgn, SRn, kfi, sbi, r_normal, &
                       n_iters, pic_flag )

! Given a reduced composition (r,h,p), determine npts point on the PIC from r.
    
!  Input:
!    r:           reduced composition at reconstruction point {zr,zue}
!    h:           enthalpy (cgs)
!    p:           pressure (cgs)
!    r_tol:       tolerance on the accuracy of the reconstruction
!                 |rR - r| <= r_tol
!    r_acc:       acceptance tolerance:  on iteration failure, accept last (or
!                 previous) estimate if |rR - r| <= r_acc
!    Sg_min:      minimun rate of change (kmol/kg/s).
!                 PIC is terminated if |dz/dt| < Sg_min
!    npt:         target number of points to be generated on the PIC (>=2)
!    np_acc:      accept "very attractive (r_acc) point" if np >= np_acc (see note below)
!    npmx:        maximum number of points to be generated on the PIC (>=2)

! Output:
!    np           number of PI points generated
!    ipt_best     =np if boundary is reached; otherwise index of most attractive PI point
!    z_CE         CEQ point
!    rg           reduced composition of the PI point 
!    zg           full composition of the PI point 
!    Tg           temperature of the PI point 
!    taug         time from the PI point 
!    rR           reduced composition of the mapped point 
!    zR           full composition of the mapped point 
!    TR           temperature of the mapped point 
!    rg_ex        extrapolation (forward on PIC) for boundary generating point
!    zR_ex        extrapolation (backward on PIC) for ICE composition, zR
!    TR_ex        extrapolation (backward on PIC) for ICE temperature, TR
!    Sgn          norm of rate of change at generating point
!    SRn          norm of rate of change at mapped point
!    kfi          index of predicted facet
!    sbi          distance to predicted facet
!    n_iters      total number of Newton iterations
!    pic_flag     >=0 for success, < 0 for failure
!                 =  0 - last PIC point on boundary (normal case)
!                 =  1 - last PIC point very close to boundary (controlled by tol_rceb)
!                 =  2 - very attractive (r_tol)
!                 =  3 - very attractive (r_acc) and becoming inert:  |Sg|<1d-3*|SR|
!                 =  4 - very attractive (r_acc) point found, but subsequent step collapse
!                 =  5 - very attractive (r_acc) point accepted, since np = npmx
!                 =  6 - very attractive (r_acc) point accepted, since np > np_acc
!                 =  7 - inert: |dz/dt| < Sg_min at the initial (CE) point
!                 =  8 - PIC terminated because inert mixture encountered (|dz/dt| < Sg_min)

!                 = -1 - failure to obtain initial CE point
!                 = -2 - npmx points generated without reaching boundary
!                 = -3 - collapse of step size
!                 = -4 - excessive perturbation of CE point
!                 = -5 - ci_ice_pic_bound, max-min both failed

!  Note: a point is deemed "very attractive (r_t)" if |zR-zR_extrapolated| < r_t.

use ci_dat
use ci_dat8
use ci_utils
use ci_ice_cksubs
use ci_cem_recon
use ci_stats
implicit none
	
integer, intent(in)     :: npt, np_acc, npmx
real(k_dp), intent(in)  :: r(nrc), h, p, r_tol, r_acc, Sg_min
      
integer,    intent(out) :: np, ipt_best, kfi(npmx), n_iters, pic_flag
real(k_dp), intent(out) :: rg(nrc,npmx), zg(ns,npmx), Tg(npmx), taug(npmx), &
                           rR(nrc,npmx), zR(ns,npmx), TR(npmx), z_CE(ns), &
                           rg_ex(nrc,npmx), zR_ex(ns,npmx), TR_ex(npmx), &
                           Sgn(npmx), SRn(npmx), sbi(npmx), r_normal(nrc)
    
!  NOTE: very small steps may be needed in pathological cases in order for the turning angle 
!        between successive segments of the PIC to be small.  Hence step_min is et very small.            
real(k_dp), parameter :: dsm_fac  = 1.d-12 !  smallest PIC initial step allowed is dsm_fac * norm(r)
real(k_dp), parameter :: step_min = 1.d-15 !  smallest fraction of PIC step allowed
real(k_dp), parameter :: rcl_fac  = 1.d-12 !  attempt solution on facet if closer than rcl_fac*norm(r)
real(k_dp), parameter :: cos_min  = 0.9d0  !  max allowed cosine between successive tangent vectors
real(k_dp), parameter :: tau_fac  = 2.d0   !  factor by which tau is allowed to change (>1)
                               
logical, parameter :: diagnostic = .false.                           
integer, save      :: lu_PIC = -1       ! set < -1 for no diag. o/p
integer, save      :: lu_PIC2, lu_ang
integer            :: op_max = 1000     ! max. no. of calls for diag. o/p 
integer, parameter :: calls_stop = -1   ! call on which to o/p statistics and stop (-1 for never)
real(k_dp), save   :: calls = 0.d0, cpu_start, cpu_end
                         
integer    :: ipt, ipl, niters, kfa, isub, icount(11), &
              iflag, indic_rpos_o(nrc), indic_zpos_o(ns), kf(npmx), iflag_CEM, &
              rankA, rankB, npim, i, info   
real(k_dp) :: zdum(ns), Tdum, stats(20), z(ns+1), dzdt(ns+1), dens, Sg(ns), &
              drds(nrc), dtauds, step_size, drds_dirn, sb, ds_ref, s, r0(nrc), &
              tau0, drds_n(nrc), dtauds_n, A_ODE(ns+1,ns+1), dr_norm, resid, &
              norm_r, dzds(ns), rc(nrc), r_acc_int, r_tol_f, &
              r_acc_f, z_max_min, r_tol_n, dzTRds(ns+1), &
              CEM_tan(ns+1,nrc+2), T_CE(ns,nrc), BBTT(nrc,nrc), rhs(nrc), &
              dr_step, dr_diff, dr_step0, dvec, dvecp, dvecq, dvecr, sv_rat, sv_rat0, &
              theta_SB(1), theta_SC(1), theta_BC(nrc), theta_BP(nrc), &
              theta_CP(nrc), svrat, T_PIM(ns,ns+1-nrc), angles(2+3*nrc), dzR_ex, &
              dzR_ex_best
              
if( calls == 0.d0 ) then
   icount = 0
   call isat_lu( lu_ang )
   open( lu_ang, file = 'ice_pic_angles.op' )
   call cpu_time( cpu_start )
   
elseif( calls_stop /= -1 .and. nint(calls) >= calls_stop ) then !  output statistics and stop
   call cpu_time( cpu_end )
   write(0,*)' '
   write(0,'(a,i5)')      'number of calls   = ', icount(1)
   write(0,'(a,1p,e13.2)')'points            = ', icount(2)/float(icount(1))
   write(0,'(a,1p,e13.2)')'subs              = ', icount(3)/float(icount(1))
   write(0,*)' '
   write(0,'(a,1p,e13.2)')'internal attempts = ', icount(4)/float(icount(1))
   write(0,'(a,1p,e13.2)')'internal success% = ', icount(5)*100.d0/float(icount(4))
   write(0,'(a,1p,e13.2)')'iters on success  = ', icount(6)/float(icount(5))
   write(0,'(a,1p,e13.2)')'iters on failure  = ', icount(7)/float(icount(4)-icount(5))
   write(0,*)' '
   write(0,'(a,1p,e13.2)')'facet attempts    = ', icount(4)/float(icount(1))
   write(0,'(a,1p,e13.2)')'facet success%    = ', icount(9)*100.d0/float(icount(8))
   write(0,'(a,1p,e13.2)')'iters on success  = ', icount(10)/float(icount(9))
   write(0,'(a,1p,e13.2)')'iters on failure  = ', icount(11)/float(icount(8)-icount(9))
   write(0,*)' '
   write(0,'(a,1p,e13.2)')'CPU time (secs)   = ', cpu_end - cpu_start
   write(0,'(a,1p,e13.2)')'CPU ( mu sec/call)= ', (cpu_end - cpu_start)*1.d6/calls
   write(0,*)'ci_ice_pic: stopping, calls = ', nint(calls)
   stop 
endif 

if( npmx < 2 )  call isat_abort('ci_ice_pic', 1, mess = &
                                'too few points, npmx = ', isv = npmx )
calls = calls + 1.d0
icount(1) = icount(1) + 1

r_tol_n = r_tol
r_tol_f = r_tol
r_acc_f = r_acc

Tdum  = 0.d0
zdum  = 0.d0
kf    = 0
n_iters   = 0
step_size = 1.d0

ipt       =  1 ! set quantities used in diagnostic_op
isub      =  0
kfa       = -1
s         = -1.d0
sb        = -1.d0
sv_rat    = -1.d0
ipt_best  =  1
dzR_ex_best = huge(0.d0)
                                
!  treat initial point at the feasible end ----------------------------

call ci_ceq( 0, r, h, p, .false., zdum, Tdum, & 
                      zg(:,1), Tg(1), stats, info )

z_CE = zg(:,1)
                      
if ( info < 0 )  then 
   write(0,*)'ci_ice_pic: initial CE failure, info = ', info
   pic_flag = -1
   call diagnostic_op( -2 )
   return
endif

if( stats(7) == 0.d0 ) then     ! check for perturbation in CEQ
   rc = r
else
   rc = matmul( BBT, zg(:,1) )  !  corrected r if it has been perturbed in CEQ
   dr_norm = norm( rc - r )
   if( dr_norm > r_acc ) then
      write(150,*) ' FAILED -4: dr_norm = ', dr_norm
      write(150,*) ' FAILED -4: norm = ', dot_product(r, amolwt_n)
      call print_var( r, "r", 150 )
      call print_var( rc, "rc", 150 )
      call print_var( z_CE, "z_CE", 150 )
      pic_flag = -4
      return
   endif
endif

! the properties of the first generated point are: zg(:,1), Tg(:,1) and... 

ipt     = 1
rg(:,1) = rc 
taug(1) = 0.d0
rR(:,1) = rc
zR(:,1) = zg(:,1)
TR(1)   = Tg(1)
norm_r  = norm( rc )

!  get Sg = S(zg)
z(1:ns) = zg(:,1)
z(ns+1) = Tg(1)
call cidzdt( z, p, dzdt, dens )
Sg = dzdt(1:ns) 

Sgn(1) = norm( Sg )
SRn(1) = Sg(1) 

!  get CEM tangent vectors
call ci_cem_tan( zg(:,1), h, Tg(1), p, thermo_ns, ns, nrc, &
                    BB, CEM_tan,  iflag_CEM )  
     
if( iflag_CEM /= 1 ) call isat_abort('ci_ice_pic', 2, &
                     mess='CEMTan failed' )
     
T_CE = CEM_tan(1:ns,1:nrc)    

! initial PIC tangent vector
!         (BBT*T_CE) * drds = - BBT*Sg * dtauds
! solve:  (BBT*T_CE) * v    = - BBT*Sg

BBTT =  matmul( BBT, T_CE )
rhs  = -matmul( BBT, Sg )
dzds =  1.d0  ! scale factors
call ci_ice_min_norm( nrc, nrc, BBTT, rhs, dzds, 0.d0, drds, Tdum, zdum, sv_rat0 )

dtauds = 1.d0 / norm( drds )  ! d(taug)/ds
drds   = drds * dtauds        ! d(rg)/ds
dzds   = matmul( T_CE, drds ) ! d(zg)/ds

dzTRds(1:ns) = dzds + Sg * dtauds  ! d(zR)/ds
dzTRds(ns+1) = dot_product( CEM_tan(ns+1,1:nrc), drds ) + dzdt(ns+1) * dtauds  ! d(TR)/ds

!//////////////////////////////////////////////////////////////////////////////////////
!  For diagnostic purposes, examine angles between sub-spaces

if( diagnostic ) then

   npim = ns+1-nrc
   A_ODE = 0.d0
   do i = 1, ns+1
      A_ODE(i,i) = 1.d0
   end do

   call ci_ice_pim_tan( A_ODE, Sg, T_PIM, svrat )

   call principal_angles( ns, 1,   nrc,  Sg,   BB,    rankA, rankB, theta_SB )
   call principal_angles( ns, 1,   nrc,  Sg,   T_CE,  rankA, rankB, theta_SC )
   call principal_angles( ns, nrc, nrc,  BB,   T_CE,  rankA, rankB, theta_BC )
   call principal_angles( ns, nrc, npim, BB,   T_PIM, rankA, rankB, theta_BP )
   call principal_angles( ns, nrc, npim, T_CE, T_PIM, rankA, rankB, theta_CP )

   write(0,*)' '
   write(0,'(a, 1p,10e13.4)') 'ICE_PIC - initial angle = ', theta_SB(1), theta_SC(1)
   write(0,'((1p,6e12.2))')theta_BC
   write(0,*)' '
   write(0,'((1p,6e12.2))')theta_BP
   write(0,*)' '
   write(0,'((1p,6e12.2))')theta_CP
   write(0,*)' '
endif

!//////////////////////////////////////////////////////////////////////////////////////

! Test for inert initial point
if( Sgn(1) < Sg_min ) then
   pic_flag   = 7
   np         = 1
   rg_ex(:,1) = rg(:,1)
   zR_ex(:,1) = zR(:,1)
   kfi        = 0
   sbi        = 0.d0
   return
endif

r_acc_int  = r_tol ! acceptable error at internal pre-image points
pic_flag   = -2    ! assume failure until success achieved

points: do ipt = 2, npmx  !  loop over pre-image points ------------------
   icount(2) = icount(2) + 1
   
   ipl = ipt-1 !  previous point
   np  = ipt   !  may be changed to ipl below   

! at the start of this loop, the following are defined at the previous PI point, ipl:
! rg, zg, Tg, taug, rR, zR, TR, drds, dzds and dtauds

   !  simple determination of direction (not foolproof)
   !  note: dtauds returned positive; no change for first point
   
   if( ipt > 2 ) then
      drds_dirn = dot_product( drds, rg(:,ipt-1)-rg(:,ipt-2) ) / &
                               norm( rg(:,ipt-1)-rg(:,ipt-2) )
      if( drds_dirn < 0.d0 ) then
         drds   = -drds
         dtauds = -dtauds
         dzds   = -dzds
         dzTRds = -dzTRds
      endif
   endif
   
   ! extrapolate to find boundary point:  rg_ex = rg + sb * drds
   call ci_ice_pic_bound( rg(:,ipl), drds, &
                       sb, kf(ipl), indic_rpos_o, indic_zpos_o, r_normal, iflag )
                       
   kfi(ipl) = kf(ipl)
   sbi(ipl) = sb
                          
   if( iflag < 0 ) then
      !  failure can be caused by rg(:,ipl) being very close to boundary
      !  find max-min composition to test
      call ceq_maxmin( ns, nrc, BB, rg(:,ipl), zR(:,ipt), z_max_min, iflag )
      
      if( iflag < 0 ) call isat_abort('ci_ice_pic', 2, mess = &
                      'ceq_maxmin failed, iflag = ', isv = iflag ) 
                      
      if( z_max_min < tol_rceb ) then
         np       = ipl
         pic_flag = 1
         exit points  !  success  ! XXX note: r_normal not computed
      endif
      
      !call isat_abort('ci_ice_pic', 2, mess = &
      !  'ci_ice_pic_bound failed, iflag, z_max_min = ', isv = iflag, rsv = z_max_min ) 
      pic_flag = -5
      return
   endif
      
!  boundary point successfully identified    
   rg_ex(:,ipl) = rg(:,ipl) + sb * drds
   zR_ex(:,ipl) = zR(:,ipl) + sb * dzTRds(1:ns)
   TR_ex(ipl)   = TR(ipl)   + sb * dzTRds(ns+1)
  
   if( sb < tol_rceb ) then  ! point (ipl) is essentially on the boundary
      np       = ipl
      pic_flag = 1
      exit points   !  success  ! XXX note: r_normal not computed
   endif  
       
    !  set reference step length
   ds_ref    = max( sb / max(1,npt+1-ipt) , dsm_fac * norm_r ) 

   isub      = 0  
   sub_steps: do  !  loop over decreasing step size until Newton converges ------
      icount(3) = icount(3) + 1
   
      isub = isub + 1

    !  set step length for this attempt
      s   = step_size * ds_ref
      
      if( s > 0.9d0 * sb ) then
         s   = sb
         kfa = kf(ipl)  ! indicate on facet kfa
         step_size = s / ds_ref
      else
         kfa = 0        ! indicate in the interior
      endif
                
    !  initial estimate of next PI point
      
      r0   = rg(:,ipl) + drds   * s
      tau0 = taug(ipl) + dtauds * s
      tau0 = max( tau0, 0.d0 )
      
      if( diagnostic  .and. ipt== -18 ) then 
         !  diagnostic examination of properties along ray
         call ci_ice_pic_Ray( r, h, p, rg(:,ipl), taug(ipl), zg(:,ipl), Tg(ipl), &
                              zR(:,ipl), TR(ipl), drds, dtauds, sb, 1.d-20, 20 )
      endif
         
 !  attempt to determine PI point 
                 
      if( kfa == 0 ) then  ! interior      
         call ci_ice_pic_newton( r, h, p, r0, tau0, r_tol_n, r_acc_int, &
                                 rg(:,ipt), taug(ipt), zg(:,ipt), Tg(ipt),  &
                                 rR(:,ipt), zR(:,ipt), TR(ipt), &
                                 dzds, drds_n, dtauds_n, dzTRds, niters, &
                                 dr_norm, resid, sv_rat, angles, info )

         if(diagnostic) write(30,'(a,i3,i4, 4i3,1p,10e10.2)') 'PN (full): ', nint(calls), ipt, isub, &
             kfa, info, niters, dr_norm, s, sb, Tg(ipt), taug(ipt), sv_rat 
                                 
         icount(4) = icount(4) + 1
         if( info >= 0 ) then
            icount(5) = icount(5) + 1
            icount(6) = icount(6) + niters
            if( diagnostic ) write(lu_ang,'(i10,i4,1p,1000e20.10)') nint(calls), ipt, angles
         else
            icount(7) = icount(7) + niters
            call routine_failed( i_ci_ice_pic_newton )
         endif
                              
      else                ! on  facet                        
         call ci_ice_facet_newton( kfa, r, h, p, tau0, r0, r_normal, r_tol_f, r_acc_f,  &
                                   taug(ipt), zg(:,ipt), Tg(ipt), &
                                   zR(:,ipt), TR(ipt), A_ODE, niters, &
                                   dr_norm, resid, info )

         if(diagnostic) write(20,'(a,i3,i4, 4i3,1p,10e10.2)') 'FN (full): ', nint(calls), ipt, isub, &
             kfa, info, niters, dr_norm, s, sb, Tg(ipt), taug(ipt), sv_rat 
                                   
         icount(8) = icount(8) + 1
         if( info >= 0 ) then
            icount(9) = icount(9) + 1
            icount(10) = icount(10) + niters
         else
            icount(11) = icount(11) + niters
            call routine_failed( i_ci_ice_facet_newton )
         endif
                                   
      endif
      
      if( diagnostic ) then 
         write(0,'(a,i3,i4, 4i3,1p,10e10.2)') 'ICE_PIC: ', nint(calls), ipt, isub, &
             kfa, info, niters, dr_norm, s, sb, Tg(ipt), sv_rat 
      endif
             
      rg(:,ipt)  = matmul( BBT, zg(:,ipt) )
      rR(:,ipt)  = matmul( BBT, zR(:,ipt) )
      n_iters    = n_iters    + niters
                                     
      if( kfa == 0  .and.  info >= 0 ) then      !  success
      
         !  terminate PIC if composition is essentially inert
         if( S_norm( zg(:,ipt), Tg(ipt) )  < Sg_min ) then
            pic_flag = 8
            np       = ipt
            exit points  
         endif
      
         dr_step  = norm( rg(:,ipt)-rg(:,ipl) )  ! magnitude of step in r
         dvec  = dot_product( drds, drds_n ) / ( norm(drds) * norm(drds_n) ) ! successive tangent vectors
                   
         if( diagnostic ) then
            dr_step0 = norm( r0 - rg(:,ipl) )       ! magnitude of Euler step in r
            dr_diff  = norm( r0 - rg(:,ipt) )       ! difference in step compared to Euler
            
            dvecp = dot_product( drds,   rg(:,ipt)-rg(:,ipl) ) / ( norm(drds)   * norm(rg(:,ipt)-rg(:,ipl)) )
            dvecq = dot_product( drds_n, rg(:,ipt)-rg(:,ipl) ) / ( norm(drds_n) * norm(rg(:,ipt)-rg(:,ipl)) )
            dvecr = 0.d0
            if( ipt > 2 ) dvecr = dot_product( rg(:,ipl)-rg(:,ipl-1), rg(:,ipt)-rg(:,ipl) ) / &
                              ( norm(rg(:,ipl)-rg(:,ipl-1)) * norm(rg(:,ipt)-rg(:,ipl)) )
         
            write(99,'(4i4,1p,100e13.4)') nint(calls), ipt, isub, niters, step_size, dr_step/dr_step0, &
                                       dr_diff/dr_step, dvec, dvecp, dvecq, dvecr, sv_rat 
         endif
         
         if( sv_rat > 1.d-8 ) then
            ! check that change in direction of tangent vector is acceptably small
            if( abs(dvec)< cos_min ) then 
               if( diagnostic ) then 
                  write(0,'(a,1p,10e13.4)')'ICE_PIC dr: ', s, dr_step, dr_diff, dr_diff/dr_step, &
                                           dot_product( drds, r0 - rg(:,ipt) ) 
               endif
               info = -1  ! rejecting step
            endif
            
         else  !  multi-dimensional null space; drds_n may not be reliable; use secant instead
            drds_n      =  ( rg(:,ipt)   - rg(:,ipl)   ) / dr_step
            dtauds_n     = ( taug(ipt)   - taug(ipl)   ) / dr_step
            dzds         = ( zg(:,ipt)   - zg(:,ipl)   ) / dr_step
            dzTRds(1:ns) = ( zR(:,ipt)   - zR(:,ipl)   ) / dr_step
            dzTRds(ns+1) = ( TR(ipt)     - TR(ipl)     ) / dr_step
         endif
      endif
      
      !  reject step if tau has changed too much
      if( ipt > 2  .and. ( taug(ipt) > taug(ipl) * tau_fac  .or.  taug(ipt) < taug(ipl) / tau_fac ) ) then
         info = -1
         if( diagnostic ) write(0,'(a,1p,10e13.4)')'rejecting due to tau change', taug(ipt)/taug(ipl)
      endif
      
      if( info >= 0 ) then  !  PI point successfully generated
      
         if( kfa == 0 ) then     ! check to see if zg is on the boundary
         
            if( minval( zg(:,ipt) ) == 0.d0 ) then     !  zg on boundary
               if( minval( rg(:,ipt) ) == 0.d0 ) then  !  regular intersection
                  kf(ipt:ipt) = minloc( rg(:,ipt) )
               else
                  kf(ipt) = nrc + 1          ! indicate irregular intersection
               endif
               kfa     = kf(ipt)
            endif
            
         endif
      
         if( kfa > 0 ) then      !  all done: boundary generating point found
            np      = ipt
            kf(ipt) = kfa
            call diagnostic_op( 0 )
            np       = ipt
            pic_flag = 0
            exit points  ! success     
         endif
         
         !  boundary not reached: check other PIC termination criteria
         
         !  estimated change in |zR| between this PI point and the boundary
         dzR_ex = sb * norm( dzTRds(1:ns) ) 
         
         if( dzR_ex < dzR_ex_best ) then  !  most attractive (r_acc) point encountered 
            dzR_ex_best = dzR_ex
            ipt_best    = ipt
         endif
         
         if( dzR_ex < r_tol ) then
             pic_flag = 2  !  very attractive (r_tol)
             np       = ipt
             exit points
             
         elseif( dzR_ex < r_acc ) then  !  very attractive (r_acc)
             
             if( S_norm( zg(:,ipt), Tg(ipt) ) < 1.d-3 * S_norm( zR(:,ipt), TR(ipt) ) ) then
                pic_flag = 3  ! acceptable success: almost inert
                np       = ipt
                exit points
             endif
          endif
          
          if(  dzR_ex_best < r_acc ) then  ! acceptable point found
                   
             if( ipt == npmx ) then
                pic_flag = 5  !  acceptable success: max. points
                np       = ipt
                exit points
                
             elseif( ipt >= np_acc ) then
                pic_flag = 6  !  acceptable success: ipt >= np_acc
                np       = ipt
                exit points
             endif
         endif
         
         drds      = drds_n    !  accept step
         dtauds    = dtauds_n
         
         call diagnostic_op( 2 )
         step_size = min( 1.d0, 2.d0 * step_size )
         r_acc_int = max( r_tol_n, 1.1d0*dr_norm )  !  acceptance tolerance for next point
         r_acc_int = min( r_acc_int, r_acc )
         
         exit sub_steps  !  proceed to next PI point
         
      else  ! failure
      
         call diagnostic_op( -1 )
         step_size = 0.25d0 * step_size  !  cut step size
         
         if( step_size < step_min ) then ! step size too small
            if( diagnostic ) write(0,*)'ci_ice_pic: substep size collapsed, = ', step_size
            
            if( dzR_ex_best < r_acc ) then
               np = ipt_best
               pic_flag = 4
            else
               np   = ipl
               pic_flag = -3  ! failure due to step size collapse
            endif
            
            call diagnostic_op( -3 )
            exit points
         endif
      endif
      
   end do sub_steps  !--------------------------------------------------------
   
end do points  !-------------------------------------------------------------

if( pic_flag == -2 ) then  ! npmx points successfully generated
   np = npmx ! maximum number of points without reaching boundary
endif

!   finish up

if( np > ipl ) then
   rg_ex(:,np) = rg(:,np) 
   zR_ex(:,np) = zR(:,np)
   TR_ex(np)   = TR(np)
   sbi(np)     = 0.d0
   kfi(np)     = kfi(np-1)
endif

!  set Sgn and SRn

do ipt = 1, np
   z(1:ns)  = zg(:,ipt)
   z(ns+1)  = Tg(ipt)
   call cidzdt( z, p, dzdt, dens )
   Sgn(ipt) = norm( dzdt(1:ns) ) 

   z(1:ns)  = zR(:,ipt)
   z(ns+1)  = TR(ipt)
   call cidzdt( z, p, dzdt, dens )
   SRn(ipt) = norm( dzdt(1:ns) ) 
end do

return

contains  !-----------------------------------------------------------------

function S_norm( zz, TT )
!  evaluate norm of rate of change
real(k_dp) :: zz(ns), TT, S_norm

   z(1:ns)  = zz(1:ns)
   z(ns+1)  = TT
   call cidzdt( z, p, dzdt, dens )
   S_norm   = norm( dzdt(1:ns) ) 
   
return
end function S_norm

subroutine diagnostic_op( flag )  !--------------------------
    
!  diagnostic output
     
integer, intent(in) :: flag

! flag = -3   failure due to step size too small
! flag = -2   failure of initial CEQ
! flag = -1   Newton failed, reduce step size

! flag =  0   success, Newton converged on boundary
! flag =  1   success, reconstruction point on boundary
! flag =  2   success, Newton converged at interior PI point

if( .not. diagnostic ) return
   
if( lu_PIC == -1 ) then
   call isat_lu( lu_PIC )
   open( lu_PIC, file='PIC.op' )   
   call isat_lu( lu_PIC2 )
   open( lu_PIC2, file='PIC2.op' )
endif
       
if( lu_PIC < 0  .or.  calls > op_max ) return
        
write(lu_PIC,'(7i4,1p,14e10.2)') nint(calls), ipt, isub, flag, info, kfa, niters, &
                                    step_size, s, sb, r_tol_n, dr_norm, resid, sv_rat, sv_rat0
                                    
if( flag == 0 ) write(lu_PIC2,'(i8,4i4)') nint(calls), ipt, kfa, niters

return
        
end subroutine diagnostic_op

subroutine ray_resid( r, h, p, r0, tau0, drds, dtauds, sb )

!  evaluate and print the residual along the ray: 
!     rg(s) = r0 + s * drds;   taug = tau0 + s * dtauds;   0 <=s <= sb.

real(k_dp), intent(in) :: r(nrc), h, p, r0(nrc), tau0, drds(nrc), dtauds, sb

integer    :: iflag_CEQ, iflag_ODE, luray
real(k_dp) :: s, ds, rg(nrc), taug, z(ns), Tg, zg(ns), stats(30), TR, zR(ns), &
              A_ODE(ns+1,ns+1), rR(nrc), res, rR0(nrc)

call isat_lu(luray)
open( luray, file='ray.op')

write(4,'(1p,100e20.11)') tau0, dtauds
write(4,'(1p,100e20.11)') r0
write(4,'(1p,100e20.11)') drds

ds = 1.d-7 * sb
s  = 0.d0
do
   
   rg   = r0   + s * drds
   taug = tau0 + s * dtauds
   taug = max( taug, 0.d0 )
   
   call ci_ceq( 0, rg, h, p, .false., z, Tg,  &
                   zg, Tg, stats, iflag_CEQ )
   if( iflag_CEQ < 0 ) exit
   
   call ci_ice_chem_map( 1, taug, p, zg, h, TR, zR, A_ODE, iflag_ODE ) 
   if( iflag_ODE > 0 ) exit
   
   rR = matmul( BBT, zR ) 
   
   if( s == 0.d0 ) then
      rR0 = rR
   else
      res = norm( rR0-rR )
      write(luray,'(1p,100e20.10)') s, res, taug, Tg, TR
   endif
   
   s  = s + ds
   ds = sqrt(2.d0) * ds
end do

return

end subroutine ray_resid

end subroutine ci_ice_recon_pic
