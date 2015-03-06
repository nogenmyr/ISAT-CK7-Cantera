!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ci_test_linrmap( nx, nf, kr_out )
  ! Test the linearized reaction mapping function

  use ci_dat8
  use ci_cksubs
  use ci_utils
  use ci_cem_recon

  implicit none
  
  integer, intent(in) :: nx, nf, kr_out

  ! local variables
  integer, parameter  :: np=60        ! number of values of s (ns >=1)
  real(kind(1.d0))    :: tau=1.d-7    ! some value of tau
  real(kind(1.d0))    :: t0=1.d-14    ! initial value of dt
  real(kind(1.d0))    :: tmax=1.d-7   ! maximum value of dt
  real(kind(1.d0))    :: s0=1.d-6     ! initial value of s (s0<=s<=smax)
  real(kind(1.d0))    :: smax=1.d0    ! maximum value of s (s0<=s<=smax)
  integer             :: nr_min=1, nr_max = 20 ! first and last record to be treated
  
  character(80) :: pfname = '', nfname = '' 
  integer       :: methodp = -1, methodn = -1 ! used only in ci_linrmap_all
  integer       :: i, k, lui, luep, luen, lust, ludz, luJ, luT, logz, logt, luout = 0
  integer       :: index_pos(ns), iref, iflag_ODE, info, ipar(ns+2), ieform, ires

  real(k_dp)    :: incz, ampz, inct, ampt, s, dt, ndt, sref, xa(nx), xb(nx), zc(ns), yc(ns), ss(np), tt(np), taupdt
  real(k_dp)    :: p, hsin_n, ha, hb, za_CE(ns), zb_CE(ns), Ta_CE, Tb_CE, Tc, CBg(1,1), errp(np), errn(np)
  real(k_dp)    :: zc_r(ns), Tc_r, Tmap, J(ns, ns), A(ns+1,ns+1), zc_linr(ns), dz(ns), ceq_stats(20), rpar(ns+2)

  call isat_lu( lui )
  call isat_lu( luep )
  call isat_lu( luen )
  call isat_lu( lust )
  call isat_lu( ludz )

  if(methodp >= 0) then
     write(pfname,'(a,i0,a)') 'lmap_ep',methodp,'.op'
  else
     write(pfname,'(a)') 'lmap_ep.op'
  endif

  if(methodn >= 0) then
     write(nfname,'(a,i0,a)') 'lmap_en',methodn,'.op'
  else
     write(nfname,'(a)') 'lmap_en.op'
  endif

  open( lui, file='badds.op' )
  open( luep, file=pfname )
  open( luen, file=nfname )
  open( lust, file='lmap_st.op' )
  open( ludz, file='lmap_dz.op' )

  ! set p to prc
  p = prc

  ! log scale or linear scale
  logz = 1
  logt = 1
  
  if(logz==1) then
     ampz  = (smax/s0)**(1.d0/(np-1))
     do i = 1, np
        ss(i)  = s0 * ampz**(i-1)
     enddo
  else
     incz = (smax - s0)/(np - 1)
     do i = 1, np
        ss(i)  = s0 + incz*(i-1)
     enddo
  endif

  if(logt==1) then
     ampt  = (tmax/t0)**(1.d0/(np-1))
     do i = 1, np
        tt(i)  = t0 * ampt**(i-1)
     enddo
  else
     inct = (tmax - t0)/(np - 1)
     do i = 1, np
        tt(i)  = t0 + inct*(i-1)
     enddo
  endif

  ! store s and dt values in lmap_st.op
  write( lust,'(1p,200e25.16)' ) ss
  write( lust,'(1p,200e25.16)' ) tt
  write( lust,'(1p,200e25.16)' ) tau - tt

  do nr = 1, nr_max
     read(lui,*,end=100,err=100) xa, xb
     if( nr < nr_min ) cycle

     index_pos = (/(i, i=1,ns)/)

     hsin_n = xa(nrs+ne+1)
     call hs2h_n( hsin_n, href_n, xa(1:nrc), nrc, ha)

     ! compute za_CE = CE(xa)
     call ci_ceq( 0, xa(1:nrc), ha, p, .false., za_CE, Ta_CE, & 
          za_CE, Ta_CE, ceq_stats, info, luout )

     hsin_n = xb(nrs+ne+1)
     call hs2h_n( hsin_n, href_n, xb(1:nrc), nrc, hb )

     ! compute zb_CE = CE(xb)
     call ci_ceq( 0, xb(1:nrc), hb, p, .false., zb_CE, Tb_CE, & 
          zb_CE, Tb_CE, ceq_stats, info, luout )
     
     ! write delta_z to lmap_dz.op file
     write( ludz, '(1p,e25.16)' ) norm(zb_CE - za_CE)

     ! vary dt and dz and store error in lmap_e.op
     do k = 1, np
        dt = tt(k)
        
        do i = 1, np
           s = ss(i)
           zc = za_CE + s*(zb_CE - za_CE)
           call phi2y( zc, amolwt, ns, yc )
           call temphy( ha, yc, Tc )
           dz = zc - za_CE

           ! --> test positive dt <--
           ! linearized reaction mapping using mapping at za_CE
           ! call ci_linrmap_all( za_CE, tau, dz, dt, Ta_CE, ha, p, zc_linr, methodp, methodn )
           call ci_linrmap( za_CE, tau, dz, dt, Ta_CE, ha, p, zc_linr, Tmap, A, info )
           
           ! actual reaction mapping
           call ci_ice_chem_map( 0, tau+dt, p, zc, ha, Tc_r, zc_r, A, iflag_ODE )

           errp(i) = norm(zc_linr - zc_r)

           ! --> test negative dt <--
           ! linearized reaction mapping using mapping at za_CE
           ! call ci_linrmap_all( za_CE, tau, dz, -dt, Ta_CE, ha, p, zc_linr, methodp, methodn )
           ndt = -dt
           call ci_linrmap( za_CE, tau, dz, ndt, Ta_CE, ha, p, zc_linr, Tmap, A, info )

           ! actual reaction mapping
           taupdt = tau - dt
           if( taupdt < 0.d0 ) then
              ! write(0, *) 'Warning: Setting (tau - dt) less than zero = ',taupdt, ' to zero.'
              taupdt = 0.d0
           endif
           call ci_ice_chem_map( 0, taupdt, p, zc, ha, Tc_r, zc_r, A, iflag_ODE )

           errn(i) = norm(zc_linr - zc_r)

        enddo

        write( luep,'(1p,200e25.16)' ) errp
        write( luen,'(1p,200e25.16)' ) errn

     enddo
  enddo

100 continue
  return

end subroutine ci_test_linrmap
