!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ci_ice_recon( r, h, p, z_CE, rg, zg, Tg, tau, r_ICE, z_ICE, T_ICE, A_ODE, &
     r_normal, kfacet, n_iters, n_integs, method, info )

  ! Determine z_ICE point on the ICE Manifold giver r, h and p.
  ! ICE point properties are: z_ICE, T_ICE, r_ICE = B^T * z_ICE
  ! Boundary generating point properties are:  tau, zg, Tg, rg=B^T*zg

  ! Input:
  !    r:           reduced composition at reconstruction point {zr,zue}
  !    h:           enthalpy (cgs)
  !    p:           pressure (cgs)

  ! Output:
  !    z_CE:        CEQ point
  !    rg:          reduced generating point
  !    zg:          generating point
  !    Tg:          temperature at the generating point
  !    tau:         integration time from zg to z_ICE
  !    r_ICE:       reduced mapping point 
  !    z_ICE:       mapped point on ICE Manifold
  !    T_ICE:       temperature at T_ICE
  !    A_ODE:       sensitivity matrix
  !    kfacet:      index of predicted facet
  !    n_iters:     no. of iterations
  !    n_integs:    no. of ODE integrations with map = 1
  !    method:      routine that succeeded
  !          = 0      special cases
  !           info = 0, z_CE close to equilibrium point
  !           info = 1, z_CE close to the boundary
  !          = 1      ci_ice_recon_bc_lmap
  !          = 2      ci_ice_recon_bc
  !          = 3      ci_ice_pic_lmap
  !          = 4      ci_ice_pic
  !    info: >= 0    successful outcome from the method used;
  !                 look at the individual routines for more info.
  !          = -1   all methods failed

  use ci_dat
  use ci_dat8
  use ci_stats
  use ci_utils
  use ci_cem_recon
  implicit none

  real(k_dp), intent(in)  :: r(nrc), h, p
  
  integer, intent(out)    :: kfacet, n_iters, n_integs, method, info
  real(k_dp), intent(out) :: z_CE(ns), rg(nrc), zg(ns), Tg, tau, r_ICE(nrc), z_ICE(ns), & 
       T_ICE, A_ODE(ns+1, ns+1), r_normal(nrc)
  
  ! local variables
  real(k_dp), parameter :: r_tol     = 1.d-8
  real(k_dp), parameter :: r_acc     = 1.d-6
  real(k_dp), parameter :: rmin_tol  = 1.d-20 ! r <= rmin_tol => near boundary assumption
  
  ! special check variables
  integer :: ceq_flag, minl(1)
  real(k_dp) :: rmin, stats(20)

  ! ci_ice_recon_bc variables (same for lmap also)
  integer    :: bc_flag, iters_bc, kfc, kfa
  real(k_dp) :: tau_bc, rg_bc(nrc), zg_bc(ns), Tg_bc, rR_bc(nrc), zR_bc(ns), TR_bc, &
       A_g_bc(ns+1,ns+1), dr_norm, resid, dr_norm_bc

  ! ci_ice_pic variables (same for lmap also)
  integer, parameter    :: npmx = 100
  integer    :: npt, np_acc, np, ipt_best, iters_pic, info_pic, kfi(npmx), iflag
  real(k_dp) :: rg_pic(nrc,npmx), zg_pic(ns,npmx), Tg_pic(npmx), &
       taug(npmx), rR_pic(nrc,npmx), zR_pic(ns,npmx), TR_pic(npmx), &
       rg_ex(nrc,npmx), zR_ex(ns,npmx), TR_ex(npmx), &
       Sgn(npmx), SRn(npmx), sbi(npmx), zz(ns), Sg_min

  ! ------------------------ special checks ------------------------
  method = 0
  ! 0/ z_CE near equilibirum not yet checked XXX

  ! 1/ z_CE close to boundary
  rmin = minval(r(1:nrc))
  if( rmin <= rmin_tol ) then
     minl = minloc( r(1:nrc) )
     kfacet = minl(1)
     call ci_ceq( 0, r, h, p, .false., z_CE, Tg, z_CE, Tg, stats, ceq_flag )
     zg = z_CE
     rg = matmul( BBT, zg )
     z_ICE = z_CE
     r_ICE = matmul( BBT, z_ICE )
     T_ICE = Tg
     tau = 0.d0
     minl = minloc( r(1:nrc) )
     kfacet = minl(1)
     A_ODE  = eye(ns+1)
     r_normal = 0.d0
     r_normal(kfacet) = -1.d0
     info = 1
     return
  endif

  ! ------------------------ special checks ------------------------

  ! n_iters and n_integs not implemented yet
  n_iters = -1
  n_integs = -1

  ! ------------------   ci_ice_recon_bc_lmap   -----------------------------
  method = 1

  call routine_start(i_ci_ice_recon_bc_lmap)
  
  call ci_ice_recon_bc_lmap( r, h, p, r_tol, r_acc, z_CE, tau_bc, rg_bc, zg_bc, Tg_bc, rR_bc, & 
       zR_bc, TR_bc, A_g_bc, iters_bc, dr_norm_bc, resid, kfc, r_normal, bc_flag )

  call routine_stop(i_ci_ice_recon_bc_lmap)

  if( bc_flag >= 0 ) then
     rg = rg_bc
     zg = zg_bc
     Tg = Tg_bc
     tau = tau_bc
     kfacet = kfc
     info = bc_flag

     call ci_ice_chem_map( 1, tau, p, zg, h, T_ICE, z_ICE, A_ODE, iflag ) 
     if( iflag > 0 ) then
        call isat_abort( 'ci_ice_Map', 2, mess='ci_ICE_chem_map failed = ', rsv=zg(ns) )
     endif
     r_ICE = matmul( BBT, z_ICE )
     return
  else
     call routine_failed(i_ci_ice_recon_bc_lmap)
  endif

  ! ------------------ END ci_ice_recon_bc_lmap -----------------------------

  ! ------------------   ci_ice_recon_bc   -----------------------------
  method = 2
  call routine_start(i_ci_ice_recon_bc)
  
  call ci_ice_recon_bc( r, h, p, r_tol, r_acc, z_CE, tau_bc, rg_bc, zg_bc, Tg_bc, rR_bc, & 
       zR_bc, TR_bc, A_g_bc, iters_bc, dr_norm_bc, resid, kfc, r_normal, bc_flag )

  call routine_stop(i_ci_ice_recon_bc)

  if( bc_flag >= 0 ) then
     rg = rg_bc
     zg = zg_bc
     Tg = Tg_bc
     r_ICE = rR_bc
     z_ICE = zR_bc
     T_ICE = TR_bc
     tau = tau_bc
     kfacet = kfc
     info = bc_flag
     A_ODE = A_g_bc
     return
  else
     call routine_failed(i_ci_ice_recon_bc)
  endif

  ! ------------------ END ci_ice_recon_bc -----------------------------
  
  ! ------------------   ci_ice_pic_lmap  -----------------------
  method = 3
  npt    = 2
  np_acc = 4
  Sg_min = 1.d-2

  call routine_start(i_ci_ice_recon_pic_lmap)

  call ci_ice_recon_pic_lmap( r, h, p, r_tol, r_acc, Sg_min, npt, np_acc, &
       npmx, np, ipt_best, z_CE, rg_pic, zg_pic, Tg_pic, taug, rR_pic, zR_pic, & 
       TR_pic, rg_ex, zR_ex, TR_ex, Sgn, SRn, kfi, sbi, r_normal, iters_pic, info_pic )

  call routine_stop(i_ci_ice_recon_pic_lmap)

  if( info_pic >= 0 ) then
     rg = rg_pic(1:nrc, np)
     zg = zg_pic(1:ns, np)
     Tg = Tg_pic(np)
     tau = taug(np)
     kfacet = kfi(np)
     info = info_pic
     
     call ci_ice_chem_map( 1, tau, p, zg, h, T_ICE, z_ICE, A_ODE, iflag ) 
     if( iflag > 0 ) then
        call isat_abort( 'ci_ice_Map', 3, mess='ci_ICE_chem_map failed = ', rsv=zg(ns) )
     endif

     r_ICE = matmul( BBT, z_ICE )
     return
  else
     call routine_failed(i_ci_ice_recon_pic_lmap)
  endif

  ! ------------------ END ci_ice_pic_lmap -----------------------


  ! ------------------   ci_ice_pic  -----------------------
  method = 4
  npt    = 2
  np_acc = 4
  Sg_min = 1.d-2

  call routine_start(i_ci_ice_recon_pic)

  call ci_ice_recon_pic( r, h, p, r_tol, r_acc, Sg_min, npt, np_acc, &
       npmx, np, ipt_best, z_CE, rg_pic, zg_pic, Tg_pic, taug, rR_pic, zR_pic, & 
       TR_pic, rg_ex, zR_ex, TR_ex, Sgn, SRn, kfi, sbi, r_normal, iters_pic, info_pic )

  call routine_stop(i_ci_ice_recon_pic)

  if( info_pic >= 0 ) then
     rg = rg_pic(1:nrc, np)
     zg = zg_pic(1:ns, np)
     Tg = Tg_pic(np)
     tau = taug(np)
     kfacet = kfi(np)
     info = info_pic
     
     call ci_ice_chem_map( 1, tau, p, zg, h, T_ICE, z_ICE, A_ODE, iflag ) 
     if( iflag > 0 ) then
        call isat_abort( 'ci_ice_Map', 4, mess='ci_ICE_chem_map failed = ', rsv=zg(ns) )
     endif
     r_ICE = matmul( BBT, z_ICE )
     return
  else
     call routine_failed(i_ci_ice_recon_pic)
  endif

  ! ------------------ END ci_ice_pic -----------------------

  ! all methods failed ( return info from last method )
  call routine_failed( i_ci_ice_recon )
  info = info_pic
  return
  
end subroutine ci_ice_recon
