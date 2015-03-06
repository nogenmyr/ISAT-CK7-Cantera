!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ci_ice_facet_newton( kfa, r, h, p, tau0, r0, &  
                                r_normal0, r_tol, r_acc, tau, zg, Tg, zR, &
                                TR, A_g, iter, dr_norm, resid,  iter_flag )
                                    
! Perform Newton iterations to determine manifold generating point, zg.
! Also determined are the mapping, zR and the mapping gradient A_g.
! The initial guess r0 is a point on the boundary facet kfa on which the
! outward normal is r_normal0.

!  Input:
!    kfa:         index of facet
!    r:           reduced composition at reconstruction point {zr,zue}
!    h:           enthalpy (cgs)
!    p:           pressure (cgs)
!    tau0:        estimated integration time
!    r0:          estimated of reduced composition of generating point, rg
!    r_normal0:   outward unit normal on facet
!    r_tol:       tolerance on |rR-r|
!    r_acc:       acceptable error |rR-r|

! Output: 
!    tau:         integration time
!    zg:          composition of boundary generating point
!    Tg:          temperature at zg
!    zR:          mapping from zg
!    TR:          temperature at zR
!    A_g:         mapping gradient
!    iter:        number of Newton iterations
!    dr_norm:     residual |rR-r| (before last Newton)
!    resid:       residual predicted by min_norm

!    iter_flag = 0, success with |rR-r| < r_tol
!              = 1, success with |rR-r| < r_acc
!              < 0 failure
!              = -1, atten < atten_min .and. residual not decreased
!              = -2, atten < atten_min .and. no intersection point
!              = -3, not converged in iter_max
!              = -4, failure in CEQ, T < T_low
!              = -5, failure in ci_ice_chem_map

! diagnostic flags:
!              = -10, reject step: attenuate previous step and continue
!              = -11, end of do CEQ_TRY
!              = -12, end of do NEWTON

  use ci_dat
  use ci_dat8
  use ci_utils
  use ci_ice_cksubs
  use ci_cem_recon
  use ci_stats

  implicit none

  real(k_dp), parameter :: atten_min = 0.1d0  ! min. atten of Newton step
  real(k_dp), parameter :: tau_fac   = 1.d0   ! limit on tau increment
  real(k_dp), parameter :: mn_fac    = 0.1d0  ! min. norm tol: mn_tol = mn_fac * r_tol

  integer, intent(in)     :: kfa
  real(k_dp), intent(in)  :: r(nrc), h, p, tau0, r0(nrc), &
       r_normal0(nrc), r_tol, r_acc

  integer,    intent(out) :: iter, iter_flag
  real(k_dp), intent(out) :: tau, zg(ns), Tg, zR(ns), TR, &
       A_g(ns+1,ns+1), dr_norm, resid

  ! local arrays
  integer, save :: lu_NF = -1  ! set < -1 for no diag. o/p
  integer :: op_max = 1000 ! max. no. of calls for diag. o/p

  integer   :: info, iflag, indic_zpos(ns), kfa_g, &
       iflag_CEQ, iflag_ODE, indic_rpos(nrc)

  real(k_dp), save :: calls = 0.d0
  real(k_dp):: dr_norm_last, rg(nrc), dr(nrc), t_scale, &
       M(nrc, nrc+1), dr_dtau(nrc+1), z(ns+1), &
       sv_rat, stats(20), CEM_tan(ns+1,nrc+2), T_CE(ns,nrc), &
       mn_tol, drg(nrc), alpha, r_ref, beta, &
       S_norm, S_g0(ns), S_R(ns), atten, S_g(ns), &
       TR_min=huge(0.d0), tau_last, vec(nrc+1), &         
       drt_scl(nrc+1), dtau, beta_lim=0.1d0, &
       rgp(nrc), sb, r_normal(nrc), zg_last(ns), Tg_last, &
       n_dot_v, a_vec, dr_norm_p, rg_last(nrc), &
       tau_inc, tau_lim, drga(nrc)

  ! stats
  call routine_start(i_ci_ice_facet_newton)

  calls = calls + 1.d0

  !  properties of initial point
  rg       = r0  
  tau      = max( tau0, 0.d0 )
  kfa_g    =  kfa
  r_normal =  r_normal0

  ! initialization     
  mn_tol    =  mn_fac * r_tol
  iter_flag = -1                 !  anticipate failure
  iflag_ODE =  1
  atten     =  1.d0

  ! initialize quantities that may be o/p      
  t_scale  =  0.d0  
  r_ref    =  0.d0
  S_norm   =  0.d0
  sv_rat   =  0.d0
  resid    = -1.d0
  dr_norm  = huge(0.d0)
  iter     =  0

  ! determine CE point, zg  
  call ci_ceq( 0, rg, h, p, .false., z, Tg, &
       zg, Tg, stats, iflag_CEQ )

  if ( iflag_CEQ < 0 )  then 
     dice(16) = dice(16) + 1.d0
     iter_flag = -4  !  failure due to T_CE out of range
     call diagnostic_op( iter_flag )
     call routine_stop(i_ci_ice_facet_newton)
     return
  endif

  call put_last  !  save info from initial point

  !==============================================================================

  !  Starting point of Newton iterations: rg = {zr,ze} on facet kfa.

  NEWTON : do iter = 1, iter_max  !  start of Newton iterations =============
     ! At the start of each iteration, the current estimate of the boundary
     ! generating point is: zg, tau.  The previous estimate is: zg_last, tau_last

     ! recompute rg which may have been perturbed in CEQ       
     rg = matmul( BBT, zg )

     !  get CEM tangent vectors
     call ci_cem_tan(zg, h, Tg, p, thermo_ns, ns, nrc, &
          BB, CEM_tan,  info)  

     if( info /= 1 ) call isat_abort('ci_ice_Facet_newton', 2, &
          mess='CEMTan failed' )

     T_CE = CEM_tan(1:ns,1:nrc)    

     !  get S(z^g)
     call ciS( zg, Tg, p, S_g )
     S_norm = norm( S_g )
     if( iter == 1 ) S_g0 = S_g 

     if( S_norm <= 0.d0 ) call isat_abort('ci_ice_Facet_newton', 2, &
          mess='S_norm = 0' )     

     ! integrating the ODE for R and A

     iflag_ODE = 0;
     call ci_ice_chem_map( 1, tau, p, zg, h, TR, &
          zR, A_g, iflag_ODE ) 

     dice(25) = dice(25) + 1.d0 ! count of ODE integrations
     if ( iflag_ODE > 0 ) then
        dice(17)  = dice(17) + 1.d0  !  count of failures
        iter_flag = -5
        exit NEWTON  !failure
     endif

     !  residual
     dr      = r - matmul( BBT, zR )
     dr_norm = norm( dr )        ! residual

     !  make decisions on if and how to continue iterations  

     if( dr_norm <= r_tol ) then
        iter_flag = 0   ! success: Newton has converged
        exit NEWTON

     elseif( iter == iter_max ) then
        if( dr_norm <= r_acc ) then 
           iter_flag = 1 !  acceptable error...
        else
           iter_flag = -3 !  ...otherwise, failure: not converged in iter_max its
        endif
        exit NEWTON  

     elseif( dr_norm >= 0.9d0 * dr_norm_last ) then
        !  residual has not decreased adequately
        if( dr_norm <= r_acc ) then
           iter_flag = 1 !  acceptable error...
           exit NEWTON 
        endif

        ! reject step: attenuate previous step and continue
        ! (Note: if here, iter > 1)
        call diagnostic_op( -10 )
        call get_last
        atten = 0.25d0 * atten 
        if( atten < atten_min ) then 
           iter_flag = -1
           exit NEWTON ! failure
        endif
        cycle NEWTON
     endif

     !  Newton step accepted: store values

     rg = matmul( BBT, zg )
     call put_last

     atten    = min( 2.d0 * atten, 1.d0 ) !attenuation factor for step

     !  Form Newton matrix
     M(1:nrc,1:nrc) = matmul( BBT, matmul( A_g(1:ns,1:ns), T_CE ) )
     M(:,nrc+1)     = matmul( BBT, matmul( A_g(1:ns,1:ns), S_g ) ) 

     !  specify scaling factors
     r_ref          = maxval( rg )
     drt_scl(1:nrc) = r_ref    
     drt_scl(nrc+1) = r_ref / S_norm

     ! solve for [dr; dtau]

     call ci_ice_min_norm( nrc, nrc+1, M, dr, drt_scl, mn_tol, &
          dr_dtau, resid, vec, sv_rat )

     !  Make correction to drg to satisfy the normalization condition
     !  (Note that vec automatically satisfies this condition.)

     drg = dr_dtau(1:nrc)  !  change in reduced composition
     call re_norm( drg, amolwt_n, beta=beta )

     if( abs(beta)  > beta_lim ) then
        ! write(0,*)'ci_ice_pic_Facet: beta= ', beta !XXX
        ! re-solve exactly
        call ci_ice_min_norm( nrc, nrc+1, M, dr, drt_scl, 0.d0, &
             dr_dtau, resid, vec, sv_rat )
     endif

     drg  = dr_dtau(1:nrc)  !  change in reduced composition
     dtau = dr_dtau(nrc+1)  !  change in tau

     !  Take Newton step on the existing facet without checking realizability 
     !  Provisionally accept step if (a) CEQ succeeds and (b) the predicted residual decreases

     CEQ_TRY : do  !  do just once, but possibly exit
        n_dot_v = dot_product( r_normal, vec(1:nrc) )
        if( abs( n_dot_v ) < 1.d-3 * norm( vec(1:nrc) ) ) exit CEQ_TRY

        ! rg = rg + atten * ( I - v*n^T/(n^T*v) ) * drg

        drga  =  atten * drg 
        a_vec = -dot_product( r_normal, drga ) / n_dot_v 
        rg    =  rg + drga + a_vec * vec(1:nrc)    ! new point

        tau_inc = atten*dtau + a_vec * vec(nrc+1)  ! increment in tau

        tau_lim = tau_fac* max( tau_last, tau0 )   ! limit on increase in tau
        tau_inc = min( tau_inc, tau_lim )
        tau_inc = max( tau_inc, -tau )
        tau     = tau + tau_inc                    !  new tau

        call ci_ceq( 0, rg, h, p, .false., z, Tg, &
             zg, Tg, stats, iflag_CEQ )

        if( iflag_CEQ < 0 ) exit CEQ_TRY  !  reject point

        ! if CEQ has perturbed composition, check that residual is predicted to decrease
        if( stats(7) /= 0.d0 ) then 
           dr_dtau(1:nrc) = rg - rg_last
           dr_dtau(nrc+1) = tau - tau_last
           drga = matmul( M, dr_dtau )
           dr_norm_p = norm( drga )
           if( dr_norm_p >= dr_norm ) exit CEQ_TRY   !  reject point
        endif

        call diagnostic_op( -11 )
        cycle NEWTON  ! provisionally accept step

     end do CEQ_TRY

     !  acceptable point not found on current facet 
     !  find boundary point: rg = rg + atten * drg * s * vec   

     do  ! decrease attenuation factor as needed to find intersection
        if( atten < atten_min ) then
           iter_flag = -2
           exit NEWTON  !  failure
        endif

        rgp = rg_last + atten * drg 

        call ci_ice_pic_bound( rgp, vec(1:nrc), sb, kfa_g, indic_rpos, &
             indic_zpos, r_normal, iflag ) 

        if( iflag >= 0 ) then
           exit  !  success:  intersection found

        elseif( iflag >= -2 ) then
           ! no intersection of line with realizable region
           atten = 0.25d0 * atten  !  attenuate and try again

        else
           call isat_abort('ci_ice_Facet_newton', 3, &
                mess='ci_ice_pic_bound failed', isv=iflag )
        endif
     end do

     !  new boundary point (on facet kfa_g with normal r_normal)    
     rg = rgp + sb * vec(1:nrc)                 ! new point 

     tau_inc = atten * dtau + sb * vec(nrc+1)   
     tau_lim = tau_fac * max( tau_last, tau0 )
     tau_inc = min( tau_inc, tau_lim )
     tau_inc = max( tau_inc, -tau )
     tau     = tau + tau_inc                    !  new tau

     !  Obtain z^g = z^CE(rg)

     call ci_ceq( 0, rg, h, p, .false., z, Tg, &
          zg, Tg, stats, iflag_CEQ )

     !  ci_ice_ceq may fail  due to low/high temperature 

     if ( iflag_CEQ < 0 )  then 
        dice(16) = dice(16) + 1.d0
        iter_flag = -4
        exit NEWTON  !  failure due to T_CE out of range
     endif

     call diagnostic_op( -12 )

  end do NEWTON  !==========end of Newton iterations ==================

  if( iter_flag>=0 .and. TR < TR_min ) then
     TR_min = TR
     !XXX write(0,'(a,i5,1pe13.4)') 'Facet, low T ', nint(calls), TR_min
  endif

  call diagnostic_op( iter_flag )

  !stats
  call routine_stop(i_ci_ice_facet_newton)

  return

contains  !-------------------------------------------------

  subroutine put_last
    !  save results from last iteration

    rg_last      = rg
    tau_last     = tau 
    zg_last      = zg 
    Tg_last      = Tg 
    dr_norm_last = dr_norm

    return
  end subroutine put_last

  subroutine get_last  !------------------------------------
    !  get results from last iteration

    rg      = rg_last
    tau     = tau_last 
    zg      = zg_last 
    Tg      = Tg_last 
    dr_norm = dr_norm_last

    return
  end subroutine get_last


  subroutine diagnostic_op( flag )  !--------------------------

    !  diagnostic output

    integer, intent(in) :: flag

    if( lu_NF == -1 ) then
       call isat_lu( lu_NF )
       open( lu_NF, file='NF.op' )
    endif

    if( lu_NF < 0  .or.  calls > op_max ) return

    call ciS( zR, TR, p, S_R )

    write(lu_NF,'(4i6,1p,14e13.4)') nint(calls), iter, flag, kfa_g, &
         dr_norm, sv_rat, alpha, tau, Tg, TR, &
         r_ref, S_norm, norm_dot(S_g,S_g0), &
         norm_dot(S_g,S_R), atten, r_tol, r_acc, resid

    return

  end subroutine diagnostic_op

end subroutine ci_ice_facet_newton
