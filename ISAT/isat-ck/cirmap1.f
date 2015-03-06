!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

        subroutine cirmap1( need, nx, x, nf, nh, iusr, rusr,
     1	                  f, dfdx, hvar )

!  Determine the reaction mapping (f), its gradient (dfdx), and additional 
!  variables (hvar) as a function of the ISAT tabulation parameters (x).

!  ISATAB routine usrfgh  28/4/00

!  definitions:
!   x = phi(0), hs(0), p, dt  for mode_pdt = 1, general case,   nx=nsp3
!   x = phi(0), hs(0), p      for mode_pdt = 2, const dt,       nx=nsp2
!   x = phi(0), hs(0), dt     for mode_pdt = 3, const p,        nx=nsp2
!   x = phi(0), hs(0)         for mode_pdt = 4, const p and dt, nx=nsp1

!   f  = phi(dt), hs(dt), T(dt);   nf = nsp2

!   hvar - none

!  input:
!       need(1)	= 1 - f is needed (otherwise need(1) = 0 )
!       need(2)	= 1 - df/dx is needed
!       need(3)	= 1 - hvar is needed
!	nx	- number of parameters
!	x	- parameters 
!	nf	- number of dependent variables (nf=nsp2)
!	nh	- number of additional variables (nh=0)
!  output:
!	f	- variables
!	dfdx	- df/dx
!	hvar	- additional variables - none
!  notes:
!	Calls  rmap2  to determine dz/dv, where:
!	   z = {phi(dt), T(dt)},   v = {phi(0), T(0), p, dt}
!	The subroutine  rmapt  can be used to test the correct operation
!	   of cirmap1 and rmap2.
!	Set ichdas=1 to monitor ddasac accuracy (file ddasaca.op) on calls
!	   with need(2)=1.  Set ichdas=2 to monitor on all calls - expensive!
! 	   The columns in the file ddasaca.op consists of:
!	      number of the call to cirmap1 on which record error occurred
!	      max. abs. error in species
!	      max. abs. error in temperature
!	      temperature
!	      time step

	use ci_dat
	use ci_rmap
	use ci_sens, only: sens_init, sens_end, sens_match, sens_eval,
     1                   jac_dd, sens_store
	implicit none

	integer, intent(in)     :: need(3), nx, nf, nh, iusr(1)
	real(k_dp), intent(in)  :: x(nx), rusr(1)
	real(k_dp), intent(out) :: f(nf), dfdx(nf,nx)
	real(k_dp)              :: hvar(1)

	real(k_dp) :: cp0, cp1, temp, press, dens, h, hh, hhs, phi(ns),
     1	              y(ns), zin(nsp1), hs, enth0(ns), enth1(ns),
     2	              dt, vsmall

	real(k_dp) :: dfdy(nsp2,nsp3), dzdy(nsp1,nsp3), dzdv(nsp1,nsp3),
     1	              zy(nsp1,nsp4)
     	real(k_dp) :: zr(nsp1,nsp4), jac(nsp1,nsp1), pd(nsp1*nsp1)

	real(k_dp), save :: termax=0.d0, permax=0.d0, calls=0.d0

	integer    :: i, map, ismall, mode_sens, jac_use, jac_tot
        integer    :: check = 0, info
	logical    :: call_rmap2, match
	real(k_dp) :: sens_max

! fc_xe_fe.op output
        integer, save :: lu_out = -2 ! set lu_out = -1 to enable output
        integer, save :: nout_max = 2500, nout = 0

!----------------------  start of execution  -------------------------------

	calls = calls + 1.d0
	if( need(1) == 0  .and.  need(2) == 0 ) return

!------ check input  and  determine  press  and  dt-----------------------

	if( nf /= nsp2 ) call isat_abort('cirmap1', 1,
     1	      mess = 'bad value of nf', ivar = (/ nf, nsp2 /) )
	if( nh /= 0 ) call isat_abort('cirmap1', 2,
     1	      mess = 'bad value of nh', isv = nh )

	if( mode_pdt == 1 ) then   !  variable  p  and  dt
	   if( nx /= nsp3 ) call bad_nx( nx, nsp3 )
     	   press = x(nsp2)
	   dt    = x(nsp3)
	elseif( mode_pdt == 2 ) then   !  variable  p, const  dt
	   if( nx /= nsp2 ) call bad_nx( nx, nsp2 )
     	   press = x(nsp2)
	   dt    = dtc
	elseif( mode_pdt == 3 ) then   !  const  p,  variable  dt
	   if( nx /= nsp2 ) call bad_nx( nx, nsp2 )
	   press = prc
     	   dt    = x(nsp2)
	elseif( mode_pdt == 4 ) then   !  const  p,  and  dt
	   if( nx /= nsp1 ) call bad_nx( nx, nsp1 )
	   press = prc
     	   dt    = dtc
	else
	   call isat_abort('cirmap1', 4, mess = 'bad value of mode_pdt',
     1	             isv = mode_pdt )
     	endif

!--------  enforce realizability  ----------------------------------

	if( kreal == 1 ) then
	   do i = 1, ns
	      phi(i)  = max( x(i), 0.d0 )
	   end do
	else
	   phi(1:ns)  = x(1:ns)
	   call cireal( 2, treal, phi, ismall, vsmall )
	endif

!------  determine mass fractions and temperature  ----------------------

	hs = x(nsp1)

	call phi2y( phi, amolwt, ns, y )
	call hs2h( hs, href, phi, ns, h )
	call temphz( h, phi, temp, check, info )

	zin(1:ns) = phi(1:ns)
	zin(nsp1) = temp

!---------------  initial Cp and enthalpies  ----------------------

	call ctcpbs( temp, y, cp0, gas )
	call cthml( temp, enth0, gas )
	enth0 = enth0 - href

!--------------  return to 50 with mode_sens=1 if sens_eval fails for mode_sens >=2
	mode_sens = m_sens
50	continue

!---------  determine  map  and  call_rmap2  ----------------------

	if( need(2) == 0 ) then  ! gradient not needed
	   map = 0
	   call_rmap2 = .true.
	else                     !  gradient needed
	   if( mode_sens == 1 ) then  !  sensitivities from ODE integration
	      map = 1
	      call_rmap2 = .true.
	   else                       !  sensitivities from Jacobians
	      call sens_match( nsp1, zin, press, dt, match, 
     1                       zr(1:nsp1,1), dzdv(1:nsp1,nsp3) )

	      if( match ) then    
	         map = 0
	         call_rmap2 = .false.
	      else
	         write(lu_err,*)'cirmap1: match failed'
	         map = 1
	         call_rmap2 = .true.
	      endif
	   endif
	endif

!------  call rmap2_test if so required

	if( map == 1  .and.  nint(rmap2t) /= 0 ) 
     1	   call rmap2_test( map, nsp1, nsp3, zin, press, dt, zr )

!-------------  call rmap2 -------------------------------------

	if( call_rmap2 ) then
	   if( mode_sens >= 2 ) 
     1       call sens_init( nsp1, zin, press, dt, 
     2                       mode_sens, njacs, sens_lim )

	   call rmap2( map, nsp1, nsp3, zin, press, dt, zr )

	   if( mode_sens >= 2 ) 
     1      call sens_end( nsp1, zr(1:nsp1,1), zr(1:nsp1,nsp4) )

	endif

!------------  check accuracy of ode integration  -------------

	if( ichdas == 2  .or. 
     1	   (ichdas == 1  .and.  map == 1) ) call acc_check

!----  successful integration:  transform results  ---------------

	f(1:ns) = zr(1:ns,1)

!----  enforce realizability
        do i = 1, ns
           f(i)  = max( f(i), 0.d0 )
        end do

	call phi2y( f, amolwt, ns, y )
	temp    = zr(nsp1,1)
	call cthbms( temp, y, hh, gas )
	call h2hs( hh, href, f, ns, hhs )
	call ctrhoy( press, temp, y, dens, gas )
	f(nsp1) = hhs
	f(nsp2) = temp
      
        if( lu_out == -1 ) then
           call isat_lu( lu_out )
           open( lu_out, file = 'fc_xe_fe.op', form = 'formatted' )
        endif

        ! Write phi(0), hs(0) and phi(dt), hs(dt) to lu_out
        if( lu_out > 0 .and. call_rmap2 .and. nout < nout_max ) then
           write(lu_out,'(1p,500e26.17)') x(1:nsp1), f(1:nsp1)
           nout = nout + 1
        endif

	if( need(2) == 0 ) return

!-----  extract d(z)/d(v):  z = {phi,T},  v = {phi(0), T(0), p, dt}

	if( map == 1 ) then  !  calculated by DDASAC
	   dzdv(1:nsp1,1:nsp3) = zr(1:nsp1,2:nsp4)

	else  !  calculated from Jacobians

	   call jac_dd( dt, nsp1, zr(1:nsp1,1), 
     1                press, atolc, rtolc, jac )  ! form and store J(dt)
	   pd = -reshape(jac,(/nsp1*nsp1/)) 
	   call sens_store(nsp1, dt, zr(1:nsp1,1), 0.d0, pd, nsp1*nsp1)

	   call sens_eval( nsp1, dzdv(1:nsp1,1:nsp1), 
     1                   jac_use, jac_tot, sens_max )

         if( lu_sens >= 0 ) write(14,'(2i4,1pe15.6)') jac_use, 
     1                         jac_tot, sens_max

! if computed sensitivity is too large, reject and calculate directly by DDASAC
	   if( sens_max > sens_lim ) then
	      mode_sens = 1
	      go to 50
	   endif

	   dzdv(1:nsp1,nsp2) = 0.d0	!  dz/dp set to zero
	endif

!-----  transform to d(z)/d(y):  y = {phi(0), hs(0), p, dt}  -------

	dzdy(1:nsp1,nsp1) = dzdv(1:nsp1,nsp1) / cp0     ! T
	dzdy(1:nsp1,nsp2:nsp3) = dzdv(1:nsp1,nsp2:nsp3) ! p, dt
	do i = 1, ns                                    ! phi
	   dzdy(1:nsp1,i) = dzdv(1:nsp1,i) - enth0(i) * dzdy(1:nsp1,nsp1)
	end do

!-----  transform to d(f)/d(y): f = {phi, hs, T}

	call ctcpbs( temp, y, cp1, gas )
	call cthml( temp, enth1, gas )
	enth1 = enth1 - href

	dfdy(1:ns,1:nsp3) = dzdy(1:ns,1:nsp3)              ! phi
	dfdy(nsp2,1:nsp3) = dzdy(nsp1,1:nsp3)              ! T
	dfdy(nsp1,1:nsp3) = dzdy(nsp1,1:nsp3) * cp1 + 
     1	        matmul( enth1(1:ns), dzdy(1:ns,1:nsp3) )   ! hs

!------  transform to d(f)/d(x): x = y possibly with p and/or dt omitted

	if( mode_pdt == 1 ) then
	   dfdx(1:nsp2,1:nsp3) = dfdy(1:nsp2,1:nsp3)   ! phi, hs, p, dt
	elseif( mode_pdt == 2 ) then
	   dfdx(1:nsp2,1:nsp2) = dfdy(1:nsp2,1:nsp2)   ! phi, hs, p
	elseif( mode_pdt == 3 ) then
	   dfdx(1:nsp2,1:nsp1) = dfdy(1:nsp2,1:nsp1)   ! phi, hs,
	   dfdx(1:nsp2,nsp2)   = dfdy(1:nsp2,nsp3)     !          dt
	else
	   dfdx(1:nsp2,1:nsp1) = dfdy(1:nsp2,1:nsp1)   ! phi, hs,
	endif

!===============  internal subroutines   =====================================

	contains

	subroutine bad_nx( n1, n2 )

	integer :: n1, n2

	call isat_abort('cirmap1', 3, mess = 'bad value of nx',
     1         ivar = (/ n1, n2 /) )

	end subroutine bad_nx

!-----------------------------------------------------------------------------

	subroutine acc_check

	integer, save :: ludas=0, iwarn=0
	real(k_dp)    :: atolx, rtolx,  phierr, terr 
	character(30) :: blank, head, tail, name

!  on first call open lu=ludas

	if( ludas == 0 .and. ddalog > 0 ) then
	   blank = repeat(' ',30)!  generate file name:  ddasaca_P.op
	   head  = blank
	   head  = 'ddasaca'
	   tail  = blank   
	   tail  = 'op'
	   call isat_file_name( head, -1, idproc, tail, name )

	   call isat_lu(ludas)
	   open( ludas, file = name )
	endif

	atolx = atolc
	rtolx = rtolc
	atolc = atolc * errfac
	rtolc = rtolc * errfac

	call rmap2( 0, nsp1, nsp3, zin, press, dt, zy )

	atolc   = atolx
	rtolc   = rtolx

	phierr  = 0.
	do i = 1, ns
	   phierr  = max( phierr, abs( zr(i,1) - zy(i,1) ) )
	end do
	terr    = abs( zr(nsp1,1) - zy(nsp1,1) )

	if( phierr > permax  .or.  terr > termax ) then
	   termax  = max( termax, terr )
	   permax  = max( permax, phierr )

           if(ddalog > 0) then
              write(ludas,601) calls, permax, termax, zin(nsp1), dt
              call isat_flush(ludas)
           endif
601	   format(1p,10e11.4)
	   if( iwarn == 0  .and. 
     1	        ( permax > pewarn  .or. termax > tewarn )) then
		 iwarn = 1
		 write(lu_err,*)'cirmap1: warning, inaccurate ddasac ',
     1	                   'integration: see file ddasaca.op'
	   endif
	endif

	end subroutine acc_check

	end subroutine cirmap1
