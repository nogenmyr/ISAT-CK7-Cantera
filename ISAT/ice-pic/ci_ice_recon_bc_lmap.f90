!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ci_ice_recon_bc_lmap( r, h, p, r_tol, r_acc, z_CE, &
     tau, rg, zg, Tg, rR, zR, TR, A_g, &
     iter, dr_norm, resid, kfacet, r_normal,  bc_flag )

  ! Determine ICE point using boundary continuation.
  ! ICE point properties are: zR, TR, rR=B^T*zR
  ! Boundary generating point properties are:  tau, zg, Tg, rg=B^T*zg
  ! The initial guess of the boundary generating point is {tau0,r0}.

  !  Input:
  !    r:           reduced composition at reconstruction point {zr,zue}
  !    h:           enthalpy (cgs)
  !    p:           pressure (cgs)
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

  !    bc_flag   = 0, success with |rR-r| < r_tol
  !              = 1, success with |rR-r| < r_acc
  !              = 2, r0 close to the boundary 
  !              < 0 failure
  !              = -1, niters > iter_max
  !              = -2, CEQ failure
  !              = -3, Too much perturbation
  !              = -4, ci_ice_chem_map failed
  !              = -5, S is zero
  !              = -6, atten < atten_min
  !              = -7, v_dot_n is zero
  !              = -8, ci_ice_pic_bound failed

  ! Definitions:
  !   Last accepted point:  rg, zg, Tg, tau,  zR,  TR,  rR
  !   Current (test) point: rc, zc, Tc, tauc, zRc, TRc, rRc

  use ci_dat
  use ci_dat8
  use ci_utils
  use ci_sens
  use ci_ice_cksubs
  use ci_cem_recon
  use ci_stats

  implicit none

  real(k_dp), parameter :: tau_inc_max = 3.d0    ! maximum fractional increase in tau
  real(k_dp), parameter :: tau_dec_min = 0.25d0   ! minimum fractional decrease in tau
  real(k_dp), parameter :: atten_min   = 1.d-20  ! minimum attenuation of step
  real(k_dp), parameter :: mn_fac      = 0.0d0   ! min. norm tol: mn_tol = mn_fac * r_tol

  real(k_dp), intent(in)  :: r(nrc), h, p, r_tol, r_acc

  integer,    intent(out) :: iter, bc_flag
  real(k_dp), intent(out) :: z_CE(ns), tau, rg(nrc), zg(ns), Tg, rR(nrc), zR(ns), TR, &
       A_g(ns+1,ns+1), dr_norm, resid, r_normal(nrc)

  ! local arrays
  integer   :: niters, max_iters
  integer   :: info, iflag, iflag_lmap, iflag_ODE

  real(k_dp), save :: calls = 0.d0

  real(k_dp):: dr(nrc), &
       M(nrc, nrc+1), dr_dtau(nrc+1), z(ns+1), dzdt(ns+1), &
       sv_rat, stats(20), CEM_tan(ns+1,nrc+2), T_CE(ns,nrc), &
       mn_tol, alpha, r_ref, S_norm, atten, S_g(ns), dens, &
       vec(nrc+1), drt_scl(nrc+1), sb

  integer, save :: lu_bc     
  logical :: diagnostic = .false.
  logical :: choose_c, success, RESET           
  integer   :: kfacet, kfacet_new, iflag_cem, indic_rpos_o(nrc), indic_zpos_o(ns), kfacet_c
  real(k_dp):: z0(ns), T0, r0(nrc), S0(ns), S0n, drR(nrc), ra(nrc), taua, rg_new(nrc), tau_old, &
       tau_new, zg_new(ns), zg_pert(ns), Tg_new, zR_new(ns), TR_new, v_dot_n, vdn_old, dr_dot_n, &
       rR_new(nrc), A_g_new(ns+1,ns+1), dr_norm_new, dr_norm_old, zdum(ns), Tdum, BBTT(nrc,nrc), &
       rhs(nrc), dzds(ns), drds(nrc), sv_rat0, dtauds, r_normal_c(nrc), &
       z_max_min, r_normal_new(nrc), drR_new(nrc), drR_lin(nrc), &
       drR_diff(nrc), err_lin, rgb(nrc), rgc(nrc), disb, disc, sc, v_dot_nc, &
       dr_norm_red, dr_dot_vec, dist_gc, n_dot_nc, n_dot_Sg, nc_dot_Sc, &
       zgc(ns), Tgc, tauc, S_c(ns), dzgc(ns), dzgc_l(ns), dzgc_n, dzgc_e, Adum(ns+1,ns+1)

  ! linear mapping variables
  integer  :: i, points
  real(k_dp) :: rmap1(nrc), zmap(ns), rmap2(nrc), dir(nrc), dzmap(ns), dtmap, Tmap, J(ns,ns), &
       C(ns+1, ns+1), eCdt(ns+1,ns+1), maxel

  ! global variables for linear mapping
  real(k_dp) :: zlin(ns), tlin, dzlin(ns), dtlin, Tzlin

  if( nint(calls) == 0 .and. diagnostic) then
     call isat_lu( lu_bc )
     open( lu_bc, file = 'ice_bc.op' )
  endif
  
  calls = calls + 1.d0
  max_iters = 20
  RESET = .true.

  !-----  obtain all needed information at the feasible end of the PIC -----  

  !  obtain z0 = z^CE(r) 

  call ci_ceq( 0, r, h, p, .false., zdum, Tdum, & 
       z0, T0, stats, info )

  z_CE = z0

  if ( info < 0 )  then 
     if(diagnostic) write(0,*)'ci_ice_recon_bc: initial CE failure, info = ', info
     bc_flag = -2  ! failure
     return
  endif

  if( stats(7) == 0.d0 ) then     ! check for perturbation in CEQ
     r0 = r
  else
     r0 = matmul( BBT, z0 )  !  corrected r if it has been perturbed in CEQ
     dr_norm = norm( r0 - r )
     if( dr_norm > r_tol ) then
        bc_flag = -2
        return
     endif
  endif

  !  get S0 = S(z0)
  call ciS( z0, T0, p, S0 )
  S0n = norm( S0 )

  if( S0n <= 0.d0 ) then
     bc_flag = -3
     return
  endif

  !  get CEM tangent vectors
  call ci_cem_tan( z0, h, T0, p, thermo_ns, ns, nrc, &
       BB, CEM_tan,  iflag_CEM )  

  if( iflag_CEM /= 1 ) call isat_abort('ci_ice_recon_bc', 1, &
       mess='CEMTan failed' )

  T_CE = CEM_tan(1:ns,1:nrc)    

  ! initial PIC tangent vector
  !         (BBT*T_CE) * drds = - BBT*Sg * dtauds
  ! solve:  (BBT*T_CE) * v    = - BBT*Sg

  BBTT =  matmul( BBT, T_CE )
  rhs  = -matmul( BBT, S0 )
  dzds =  1.d0  ! scale factors
  call ci_ice_min_norm( nrc, nrc, BBTT, rhs, dzds, 0.d0, drds, Tdum, zdum, sv_rat0 )

  dtauds = 1.d0 / norm( drds )  ! d(tau)/ds
  drds   = drds * dtauds        ! dr/ds
  dzds   = matmul( T_CE, drds ) ! dz/ds

  !----- find initial boundary point and obtain needed information ---------

  ! extrapolate to find boundary point:  rg = r0 + sb * drds
  call ci_ice_pic_bound( r0, drds, sb, kfacet, indic_rpos_o, indic_zpos_o, r_normal, iflag )
  
  if( iflag < 0 ) then
     !  failure can be caused by r0 being very close to boundary
     !  find max-min composition to test
     call ceq_maxmin( ns, nrc, BB, r0, zg, z_max_min, iflag )
     
     if( iflag < 0 ) call isat_abort('ci_ice_recon_bc', 2, mess = &
          'ceq_maxmin failed, iflag = ', isv = iflag ) 
     
     if( z_max_min < tol_rceb ) then
        bc_flag = 2  !  z0 is essentially on boundary
        tau = 0.d0 
        zg  = z0
        rg  = matmul( BBT, zg )
        Tg  = T0
        zR  = z0
        rR  = matmul( BBT, zR )
        TR  = T0
        A_g = eye(ns+1)
        return  !  success  ! XXX note: r_normal not computed
     endif
     
     !call isat_abort('ci_ice_recon_bc', 3, mess = &
     !     'ci_ice_pic_bound failed, iflag, z_max_min = ', isv = iflag, rsv = z_max_min ) 
     bc_flag = -8  ! failure
     return

  endif
  
  !  boundary point successfully identified   
  
  rg   = r0 + sb * drds
  tau  = norm(rg - r0)/norm(matmul(BBT, S0))
  if( tau > 10*dtc ) tau = 10*dtc ! XXX VH

  call realizable( rg )

  call compute_values

  if( .not. success ) then
     if(diagnostic) call print_var(rg, "rg", 0)
     if(diagnostic) write(0,*) 'ci_ice_recon_bc_lmap: failed at initial rg', tau
     bc_flag = -2
     return
  endif

  !-------------------------------------------------------------------------------------------------

  !  At the start of each iteration, the following are defined:
  !     boundary generating point: rg, tau, zg, Tg, r_normal
  !     reaction mapping:          rR, zR, TR, A_g
  !     step attenuation factor:   atten

  atten = 1.d0
  iter  = 0
  niters = 0

  steps: do 

     iter = iter + 1
     niters = niters + 1

     !write(0,*) 'niters = ', niters, dr_norm, tau

     if( niters > iter_max/2 ) then
        bc_flag = -1
        return
     endif

     ! compute actual norm using ci_ice_chem_map 0
     if(dr_norm <= r_acc) then
        call compute_current_norm
     endif

     !  test for termination
     if( dr_norm <= r_tol ) then
        bc_flag = 0  !  boundary generating point found 
        return
     endif

     if( atten < atten_min ) then
        bc_flag = -6
        return
     endif

     ! determine step for mapped point
     drR = atten * (r0 - rR)

     ! set up and solve linearized ICE equations to determine rg_new and tau_new 

     !  get S(z^g)
     call ciS( zg, Tg, p, S_g )
     S_norm = norm( S_g )

     if( S_norm <= 0.d0 ) call isat_abort('ci_ice_recon_bc', 7, mess='S_norm = 0' )     

     call ci_cem_tan( zg, h, Tg, p, thermo_ns, ns, nrc, &
          BB, CEM_tan,  info)  
     
     if( info /= 1 ) then
        call print_var( zg, "CEM Tan at zg failed", 0 )
        call isat_abort('ci_ice_recon_bc', 1, &
             mess='CEMTan failed' )
     endif

     T_CE = CEM_tan(1:ns,1:nrc)
     
     !  Form Newton matrix
     M(1:nrc,1:nrc) = matmul( BBT, matmul( A_g(1:ns,1:ns), T_CE ) )
     M(:,nrc+1)     = matmul( BBT, matmul( A_g(1:ns,1:ns), S_g ) ) 

     !  specify scaling factors
     r_ref          = maxval( rg )
     drt_scl(1:nrc) = r_ref    
     drt_scl(nrc+1) = r_ref / S_norm

     ! solve for [dr; dtau]
     mn_tol = mn_fac * r_tol
     call ci_ice_min_norm( nrc, nrc+1, M, drR, drt_scl, mn_tol, &
          dr_dtau, resid, vec, sv_rat )

     ! make vec^T r_normal >=0
     v_dot_n = dot_product( vec(1:nrc), r_normal )                    

     if( v_dot_n == 0.d0 ) then
        bc_flag = -7
        return
     endif

     if( v_dot_n < 0.d0 ) then
        vec     = -vec
        v_dot_n = -v_dot_n
     endif

     ! find the feasible point on the affine space containing the facet
     ! ra = rg + dr_dtau(1:nrc) + alpha * vec(1:nrc)  

     reduce_step: do ! possibly decrease step-size until change in tau is acceptable

        if( atten < atten_min ) then
           bc_flag = -6
           return
        endif

        dr_dot_n = dot_product( dr_dtau(1:nrc), r_normal ) 
        alpha    = -dr_dot_n / v_dot_n
        ra       = rg +  dr_dtau(1:nrc) + alpha * vec(1:nrc)                 
        taua     = tau + dr_dtau(nrc+1) + alpha * vec(nrc+1)
        
        ! find rg_new - feasible point on boundary                           
        call ci_ice_pic_bound( ra, vec(1:nrc), sb, kfacet_new, indic_rpos_o, indic_zpos_o, r_normal_new, iflag )

        if( iflag < 0 ) then

           ! Force a new intersection
           rg_new = rg + dr_dtau(1:nrc)
           call normalized( rg_new )

           dir = (rg_new - r0)/norm(rg_new - r0)

           call ci_ice_pic_bound( r0, dir, sb, kfacet_new, indic_rpos_o, indic_zpos_o, r_normal_new, iflag )

           if( iflag < 0 ) then
              atten = 0.25d0 * atten
              dr_dtau = 0.25d0 * dr_dtau
              if(diagnostic) write(0,*) 'second bound failed'
              cycle reduce_step
           endif

           ! new rg and tau
           rg_new = r0 + sb*dir
           tau_new = norm(rg_new - r0)/norm(matmul(BBT, S0))
           if( tau_new > tau_inc_max*tau ) tau_new = tau_inc_max*tau ! XXX VH

           call compute_new_values
           
           if(.not. success) then
              atten = 0.25d0 * atten
              dr_dtau = 0.5d0 * dr_dtau
              if(diagnostic) write(0,*) 'compute_values failed' , atten, tau
              cycle steps
           endif

           call copy_new_values

           ! reset attenuation and iter
           atten = 1.d0
           iter = 0
           cycle steps
           
        endif

        rg_new  = ra  +  sb * vec(1:nrc)                 
        tau_new = norm(rg_new - r0)/norm(matmul(BBT, S0))
        if(tau_new > tau_inc_max*tau) tau_new = tau_inc_max*tau

        ! compute new values for valid tau_new
        if( tau_new >= 0.d0 ) then
           call compute_new_values
           if(.not. success) then
              dr_dtau = 0.25d0 * dr_dtau
              atten = 0.25d0 * atten
              if(diagnostic) write(0,*) 'compute_new_values: failed' , atten, tau
              cycle steps
           endif
           exit reduce_step
        endif

        ! attenuate
        dr_dtau = 0.25d0 * dr_dtau
        atten = 0.25d0 * atten
     end do reduce_step

     !  for diagnostic purposes, compare predicted and observed changes in rR
     dr_dtau(1:nrc) = rg_new - rg
     dr_dtau(nrc+1) = tau_new - tau
     drR_lin  = matmul( M, dr_dtau )
     drR_new  = matmul( BBT, zR_new - zR )
     drR_diff = drR_new - drR_lin
     err_lin  = norm( drR_diff ) / max( norm(drR_lin), norm(drR_new) )

     if( dr_norm_new < dr_norm * ( 1.d0 - 0.1d-30 * atten ) .or. kfacet_new /= kfacet ) then

        ! accept step if facet changed irrespective of change in residue 
        call copy_new_values
        dr_norm_red = (1-dr_norm_new/dr_norm)/atten
        atten    = min( 2.d0*atten, 1.d0 )
     
     else
     
        if( dr_norm <= r_acc ) then
           bc_flag = 1  !  acceptable boundary generating point found 
           return
        endif

        atten    = 0.25d0 * atten
        dr_norm_red = dr_norm - dr_norm_new
     
     endif

     !write(0,'(a,4i4,1p,15e11.2)')'BC: ', nint(calls), iter, kfacet, kfacet_new, atten, dr_norm, dr_norm_new, dr_norm_old

  end do steps

contains

  subroutine compute_values
    success = .true.

    !  zg = z^CE(rg)
    call ci_ceq( 0, rg, h, p, .false., zdum, Tdum, & 
         zg, Tg, stats, info )

    if ( info < 0 ) then
       if(diagnostic) write(0,*)' CEQ failed'
       success = .false.
       return
    endif

    !  reaction mapping zR
    if( RESET ) then
       zlin = zg
       tlin = tau
       Tzlin = Tg
       dzlin = 0.d0
       dtlin = 0.d0
       RESET = .false.
    else
       dzlin = zg - zlin
       dtlin = tau - tlin
    endif
    call ci_linrmap( zlin, tlin, dzlin, dtlin, Tzlin, h, p, zR, TR, A_g, iflag_lmap ) 

    if ( iflag_lmap < 0 ) then
       call routine_failed(i_ci_linrmap)
       if(diagnostic) write(0,*)' linrmap failed', iflag_lmap
       success = .false.
       return
    endif

    ! reduced composition at mapped point
    rR = matmul( BBT, zR )
    dr_norm = norm( r0 - rR ) 

  end subroutine compute_values


  subroutine compute_new_values
    success = .true.

    !  new boundary point rg_new determined  
    call ci_ceq( 0, rg_new, h, p, .false., zdum, Tdum, & 
         zg_new, Tg_new, stats, info )

    if ( info < 0 ) then
       success = .false.
       return
    endif

    !  reaction mapping zR
    if( RESET ) then
       zlin = zg_new
       tlin = tau_new
       Tzlin = Tg_new
       dzlin = 0.d0
       dtlin = 0.d0
       RESET = .false.
    else
       dzlin = zg_new - zlin
       dtlin = tau_new - tlin
    endif
    call ci_linrmap( zlin, tlin, dzlin, dtlin, Tzlin, h, p, zR_new, TR_new, A_g_new, iflag_lmap ) 

    if ( iflag_lmap < 0 ) then
       call routine_failed(i_ci_linrmap)
       success = .false.
       return
    endif

    ! reduced composition at mapped point
    rR_new = matmul( BBT, zR_new )
    dr_norm_new = norm( r0 - rR_new ) 

  end subroutine compute_new_values

  subroutine compute_current_norm
    
    call ci_ice_chem_map( 0, tau, p, zg, h, TR, zR, Adum, iflag_ODE ) 
    
    if ( iflag_ODE > 0 ) then
       if(diagnostic) write(0,*) 'ci_ice_chem_map 0 failed '
       bc_flag = -4
       return
    endif

    ! reduced composition at mapped point
    rR = matmul( BBT, zR )
    dr_norm = norm( r0 - rR ) 

  end subroutine compute_current_norm

  subroutine copy_new_values

    tau      = tau_new
    rg       = rg_new
    zg       = zg_new
    Tg       = Tg_new
    rR       = rR_new
    zR       = zR_new
    TR       = TR_new
    A_g      = A_g_new
    dr_norm  = dr_norm_new
    kfacet   = kfacet_new
    r_normal = r_normal_new

  end subroutine copy_new_values

  subroutine realizable( r )

    real(k_dp), intent(inout) :: r(nrc)
    integer :: i

    do i = 1, nrc
       if(r(i) < 0.d0) r(i) = 0.d0
    enddo

    call normalized( r )

  end subroutine realizable

  subroutine normalized( r )

    real(k_dp), intent(inout) :: r(nrc)

    r = r / dot_product( r, amolwt_n )

  end subroutine normalized

  subroutine find_tau( rin, rout, tol, tau0, tau )
    
    implicit none
    real(k_dp), intent(in)  :: rin(nrc), rout(nrc), tol, tau0
    real(k_dp), intent(out) :: tau

    integer :: info, iter = 0
    real(k_dp) :: tc, zCE(ns), TCE, error, eCdt(ns+1,ns+1), J(ns,ns)
    real(k_dp) :: zmap(ns), rmap(nrc), C(ns+1, ns+1), S(ns), maxel

    iter = 0
    !  new boundary point rg_new determined  
    call ci_ceq( 0, rin, h, p, .false., zdum, Tdum, & 
         zCE, TCE, stats, info )
    
    call jacobian( tau0, ns, zCE, p, TCE, J )
    call ciS( zCE, TCE, p, S )
    C = 0.d0
    C(1:ns, 1:ns) = J
    C(1:ns, ns+1) = S

    tc = tau0
    call expm( ns+1, C*tc, eCdt, maxel )
    zmap = zCE + eCdt(1:ns, ns+1)
    rmap = matmul( BBT, zmap)
    error = norm(rmap - rout)
    iter = iter + 1
    
    bisect: do
       if( error <= tol ) exit bisect
       if( iter > 5 ) exit bisect
       tc = tc * norm(rout - rin)/norm(rmap - rin)
       call expm( ns+1, C*tc, eCdt, maxel )
       zmap = zCE + eCdt(1:ns, ns+1)
       rmap = matmul( BBT, zmap )
       error = norm(rmap - rout)
       iter = iter + 1
    enddo bisect

    tau = tc

  end subroutine find_tau
end subroutine ci_ice_recon_bc_lmap
