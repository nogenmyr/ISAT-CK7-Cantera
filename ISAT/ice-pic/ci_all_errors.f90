!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ci_all_errors( nx, nf )
  ! This subroutine computes various species reconstruction and
  ! dimension reduction errors by comparing the results obtained from
  ! the Full Chemsitry (FC) (mode = 7) case with the Dimension
  ! Reduction (DR) methods: RCCE (mode = 8) and ICE-PIC (mode = 9)

  ! Input:
  !   nx - number of components of x
  !   nf - number of components of f

  use ci_dat,  only: modeci
  use ci_dat6
  use ci_dat8
  use ci_utils
  use ci_cksubs
  use ci_ice_cksubs
  use ci_rmap
  use ci_cem_recon

  implicit none
  integer, intent(in) :: nx, nf
  

  ! local variables
  integer          :: i, s, lu_in, nr_max = 2500
  integer          :: nx_FC, nf_FC, need(3), iusr(1)
  real(kind(1.d0)) :: dt, press

  ! full mode
  real(kind(1.d0)) :: hin, hsin, hsout, T_FC, S_FC(ns)
  real(kind(1.d0)), allocatable :: x_FC(:), f_FC(:)

  ! for reduced modes 8 and 9
  real(kind(1.d0)) :: xe(nx), fe(nf), fa(nf), dfdx(nf,nx), hvar(2)
  real(kind(1.d0)) :: hsin_n, hsout_n, rin(nrc), zr(nsp1,nsp4)
  
  ! RCCE variables
  integer :: check = 0, info, lu_ceq, lu_ce
  real(kind(1.d0)) :: hsin_CE, hsout_CE
  real(kind(1.d0)) :: z_CE(ns), x_CE(nsp1), f_CE(nsp2), T_CE
  real(kind(1.d0)) :: rout_CE(nrc), stats(20), S_CE(ns)


  ! ICE-PIC variables
  integer :: k_facet, n_iters, n_integs, method, info_ice, lu_ice
  real(kind(1.d0)) :: hsin_ICE, hsout_ICE, tau_ICE, S_ICE(ns)
  real(kind(1.d0)) :: z_ICE(ns), x_ICE(nsp1), f_ICE(nsp2), rout_ICE(nrc)
  real(kind(1.d0)) :: r_e(nrc), r_g(nrc), z_e_CE(ns), z_e_ICE(ns), z_g(ns)
  real(kind(1.d0)) :: A_ICE_ODE(ns+1, ns+1), r_normal(nrc), T_e_ICE, T_g

  ! Errors
  real(kind(1.d0)) :: nerr_ce, nerr_ice
  real(kind(1.d0)) :: sr_err_ce, sr_err_ice, zmap_err_ce, zmap_err_ice
  real(kind(1.d0)) :: dzdt_err_ce, dzdt_err_ice

  ! set nx_Fc and nf_FC for the FC, mode = 7
  nx_FC = nsp1
  nf_FC = nsp1
  
  ! allocate arrays for FC
  allocate( x_FC(nx_FC) )
  allocate( f_FC(nf_FC) )

  call isat_lu( lu_in )
  open( lu_in, file='fc_xe_fe.op' )

  if(modeci == 8) then
     call isat_lu( lu_ce )
     open( lu_ce, file='rcce_xe_fe.op' )
  else if( modeci == 9 ) then
     call isat_lu( lu_ice )
     open( lu_ice, file='icepic_xe_fe.op' )
  endif
  
  ! NOTE: Everywhere assumed const. press and dt
  dt = dtc
  press = prc

  do s = 1, nr_max
     read(lu_in,*, end=100, err=100) x_FC, f_FC
     hsin = x_FC(nsp1)
     hsout = f_FC(nsp1)
     
     call hs2h( hsin, href, x_FC(1:ns), ns, hin )
     call temphz( hin, x_FC(1:ns), T_FC, check, info )

     call ciS( x_FC(1:ns), T_FC, press, S_FC )

     ! DR: mode = 8 and 9
     ! For given nrs, form the reduced representation, r.
     rin(1:nrc) = matmul(BBT(1:nrc, 1:ns), x_FC(1:ns))
     call h2hs_n( hin, href_n, rin(1:nrc), nrc, hsin_n )

     if( modeci == 8 ) then
        ! RCCE reconstruction
        z_CE = 0.d0
        T_CE = 0.d0
        call ci_ceq( 0, rin, hin, press, .false., z_CE, T_CE, &
             z_CE, T_CE, stats, info, lu_ceq )

        x_CE(1:ns) = z_CE(1:ns)
        x_CE(nsp1) = T_CE
        
        ! RCCE: reaction rate vector error
        call ciS( z_CE, T_CE, press, S_CE )
        dzdt_err_ce = norm(S_CE - S_FC)
        
        ! RCCE: normalized error
        nerr_ce = 2*norm(x_CE(1:ns) - x_FC(1:ns))/(norm(x_CE(1:ns)) + norm(x_FC(1:ns)))
        
        ! RCCE: species reconstrcution error
        sr_err_ce = norm(x_CE(1:ns) - x_FC(1:ns))
        
        call rmap2( 0, nsp1, nsp3, x_CE, press, dt, zr )

        do i = 1, ns       !  enforce realizability 
           zr(i,1)  = max( zr(i,1), 0.d0 )
        end do

        call h2hs( hin, href, z_CE(1:ns), ns, hsin_CE )
        
        x_CE(nsp1) = hsin_CE
        
        f_CE(1:ns) = zr(1:ns,1)
        f_CE(nsp2) = zr(nsp1,1)
        
        ! RCCE: species mapping error
        zmap_err_ce = norm(f_CE(1:ns) - f_FC(1:ns))

        call h2hs( hin, href, f_CE(1:ns), ns, hsout_CE )
        
        f_CE(nsp1) = hsout_CE
        
        rout_CE(1:nrc) = matmul(BBT(1:nrc, 1:ns), f_CE(1:ns))

        ! write z(0) and z(dt) values
        write(lu_ce,'(1p,500e26.17e3)') x_CE(1:nsp1), f_CE(1:nsp1)

     else if( modeci == 9 ) then
        ! ICE-PIC reconstruction
        call routine_start(i_ci_ice_recon)

        call ci_ice_recon( rin, hin, press, z_e_CE, r_g, z_g, T_g, tau_ICE,  &
             r_e, z_e_ICE, T_e_ICE, A_ICE_ODE, r_normal, k_facet,  &
             n_iters, n_integs, method, info_ice )
        
        call routine_stop(i_ci_ice_recon)

        x_ICE(1:ns) = z_e_ICE( 1:ns )
        x_ICE(nsp1) = T_e_ICE
        
        ! ICE-PIC: reaction rate vector error
        call ciS( z_e_ICE, T_e_ICE, press, S_ICE )
        dzdt_err_ice = norm(S_ICE - S_FC)
        
        ! ICE-PIC: normalized error
        nerr_ice = 2*norm(x_ICE(1:ns) - x_FC(1:ns))/(norm(x_ICE(1:ns)) + norm(x_FC(1:ns)))

        ! ICE-PIC: species reconstrcution error
        sr_err_ice = norm(x_ICE(1:ns) - x_FC(1:ns))
        
        call rmap2( 0, nsp1, nsp3, x_ICE, press, dt, zr )

        do i = 1, ns       !  enforce realizability 
           zr(i,1)  = max( zr(i,1), 0.d0 )
        end do

        call h2hs( hin, href, z_e_ICE(1:ns), ns, hsin_ICE )
        
        x_ICE(nsp1) = hsin_ICE
        
        f_ICE(1:ns) = zr(1:ns,1)
        f_ICE(nsp2) = zr(nsp1,1)
        
        ! ICE-PIC: species mapping error
        zmap_err_ice = norm(f_ICE(1:ns) - f_FC(1:ns))
        
        call h2hs( hin, href, f_ICE(1:ns), ns, hsout_ICE )
        
        f_ICE(nsp1) = hsout_ICE
        
        rout_ICE(1:nrc) = matmul(BBT(1:nrc, 1:ns), f_ICE(1:ns))

        ! write z(0) and z(dt) values
        write(lu_ice,'(1p,500e26.17e3)') x_ICE(1:nsp1), f_ICE(1:nsp1)
     endif
  enddo

100 continue
end subroutine ci_all_errors
