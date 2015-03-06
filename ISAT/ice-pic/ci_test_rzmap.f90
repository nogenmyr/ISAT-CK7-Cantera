!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ci_test_rzmap

  use ci_8
  use ci_dat8
  use ci_cem_recon
  use ci_rmap
  use ci_ice_cksubs

  implicit none

  ! local variables
  integer, parameter  :: nrmax = 100
  integer, parameter  :: ndt = 20
  real(kind(1.d0))    :: dtmin = 1.d-10
  real(kind(1.d0))    :: dtmax

  integer    :: ni, nt, lui, luo1, luo2, luo3, luo4, luo5, check, info
  real(k_dp) :: amp, dt(ndt), stats(20), tau, dz(ns)
  real(k_dp) :: p, hsn, h, hs, xa(nsp1), xb(nsp1), r(nrc), r_lmap(nrc), r_z_lmap(nrc), r_zCE_lmap(nrc)
  real(k_dp) :: z_CE(ns), T_CE, zCE_map(ns), z_map(ns), z_lmap(ns), z(ns), r_zCE_map(nrc), r_z_map(nrc)
  real(k_dp) :: zr(nsp1,nsp4), A_ODE(ns+1,ns+1), TR, Tz
  real(k_dp) :: rz_lmap_err(ndt), rzCE_lmap_err(ndt), r_zCE_lmap_err(ndt), r_lmap_err(ndt)

  call isat_lu( lui )
  call isat_lu( luo1 )
  call isat_lu( luo2 )
  call isat_lu( luo3 )
  call isat_lu( luo4 )
  call isat_lu( luo5 )

  open( lui, file='fc_xe_fe.op' )
  open( luo1, file='r_lmap_err.op' )
  open( luo2, file='r_zCE_lmap_err.op' )
  open( luo3, file='rz_lmap_err.op' )
  open( luo4, file='rzCE_lmap_err.op' )
  open( luo5, file='rzmap_errs.op' )

  ! compute dt values
  dt(1) = dtmin
  dtmax = 2*dtc
  amp  = (dtmax/dtmin)**(1.d0/(ndt-1))
  do nt = 2, ndt
     dt(nt) = dt(nt-1)*amp
  enddo

  write(luo1,'(1p,200e26.17e3)') dt
  write(luo2,'(1p,200e26.17e3)') dt
  write(luo3,'(1p,200e26.17e3)') dt
  write(luo4,'(1p,200e26.17e3)') dt

  ! set p to prc
  p = prc
  tau = 0.d0
  dz = 0.d0
  check = 0
  do ni = 1, nrmax
     read( lui,*,end=100,err=100 ) xa, xb
     z( 1:ns ) = xa( 1:ns )
     hs = xa( ns+1 )
     call hs2h( hs, href, z(1:ns), ns, h )
     call temphz( h, z(1:ns), Tz, check, info )

     r(1:nrc) = matmul( BBT, z(1:ns) )

     ! constrained-equilibrium composition, z_CE
     call ci_ceq( 0, r, h, p, .false., z_CE, T_CE, z_CE, T_CE, stats, info, luout )

     rz_lmap_err = 0.d0
     rzCE_lmap_err = 0.d0
     r_lmap_err = 0.d0
     r_zCE_lmap_err = 0.d0

     do nt = 1, ndt
        ! z_map: Full ODE soln starting from z
        call ci_ice_chem_map( 0, dt(nt), p, z, h, TR, z_map, A_ODE, info )
        r_z_map(1:nrc) = matmul(BBT, z_map)

        ! zCE_map: Full ODE soln starting from z_CE
        call ci_ice_chem_map( 0, dt(nt), p, z_CE, h, TR, zCE_map, A_ODE, info )
        r_zCE_map(1:nrc) = matmul(BBT, zCE_map)

        ! z_lmap: linearized zmap
        call ci_linrmap( z, tau, dz, dt(nt), Tz, h, p, z_lmap, TR, A_ODE, info )
        r_z_lmap(1:nrc) = matmul(BBT, z_lmap)

        write(0,*) 'ni, nt =', ni, nt

        ! linearized solutions
        call ci_rzmap(r, dt(nt), h, p, r_lmap, r_zCE_lmap)
        
        ! linearization errors
        rz_lmap_err(nt) = norm(r_z_lmap - r_z_map)
        rzCE_lmap_err(nt) = norm(r_zCE_lmap - r_zCE_map)

        ! error in r(dt) values
        r_lmap_err(nt) = norm(r_lmap - r_z_map)
        r_zCE_lmap_err(nt) = norm(r_zCE_lmap - r_z_map)

        if(r_lmap_err(nt) < 1) then
           write(luo5,'(1p,200e26.17e3)') rz_lmap_err(nt), rzCE_lmap_err(nt), r_lmap_err(nt), r_zCE_lmap_err(nt)
        endif
     enddo
     if(maxval(r_lmap_err) < 1) then
        write(luo1,'(1p,200e26.17e3)') r_lmap_err
        write(luo2,'(1p,200e26.17e3)') r_zCE_lmap_err
        write(luo3,'(1p,200e26.17e3)') rz_lmap_err
        write(luo4,'(1p,200e26.17e3)') rzCE_lmap_err
     endif
  enddo

100 continue
end subroutine ci_test_rzmap
