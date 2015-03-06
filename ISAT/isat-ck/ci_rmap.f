!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

	module ci_rmap

!  contains: rmap2, rmap1_test, rmap2_test, rmdump, ddasac0

	use ci_dat, only : ddalog
	use ci_cksubs
	use ci_stats

	contains

!============================================================================

        subroutine rmap2( map, nvar, npar, z0, press, dt, z )

!  Integrate the chemical kinetics equations to determine the reaction
!  mapping and its derivatives.

!  given z0 (initial value of {phi_0,T_0}), p = press and dt, determine
!	z = {phi,T} and (for map=1) its derivatives with respect
!      	to phi_0, T_0, p and dt.

!  input:
!	map  = 0  do not compute derivatives
!	     = 1  compute derivatives
!	nvar - number of dependent variables ( nvar = ns + 1 )
!	npar - number of parameters ( npar = ns + 3 )
!	z0   - initial value of z
!	press- pressure
!  output:
!	z    - z(:,1) = final value of z
!	     - z(:,1+j) = derivative of z w.r.t. j-th parameter
!  notes:
!   1/	ddasac is used to integrate the ode's and to determine the
!	sensitivities w.r.t. initial conditions and pressure.  The
!	sensitivity w.r.t. dt is simply the time derivative at the
!	final state.
!   2/	if  const_pr  is true, constant pressure is assumed.
!   3/	ode dependent variables = z   (1:nvar=ns+1)
!     	ode parameters          = z_0 (1:nvar=ns+1), p (ns+2=npar-1)
!	remaining parameter is    dt (ns+3=npar)
!       (number of ode parameters nspar = npar-1)
!   4/  diagnostic output is on: ddasacd.op
!       warning    output is on: ddasacw.op
!   5/  in case of complete failure, the calling arguments are dumped to
!	the file  rmdump.op.  This file can be renamed rmdump.ip and then
!	the program remap run to investigate the failure.
!

	implicit real(k_dp) (a-h, o-z), integer (i-n)

	integer, intent(in)     :: map, nvar, npar
	real(k_dp), intent(in)  :: z0(nvar), press, dt
	real(k_dp), intent(out) :: z(nvar,npar+1)

        external cires5, ciesub5, cijac5, cibsub5 !, cijac5adf


        integer    :: info(18), ipar(ns+3)
        real(k_dp) :: zp(ns+1,ns+4)
	real(k_dp) :: rtol((ns+1)*(ns+4)), atol((ns+1)*(ns+4)),
     1    rpar(ns+3), zstart(ns+1), varef(ns+1), paref(ns+3),
     1	  dzdt(ns+1)

	integer, save                 :: liwdas, lrwdas
     	integer, allocatable, save    :: iwork(:)
     	real(k_dp), allocatable, save :: rwork(:)

	character(30) :: blank, head, tail, name

        data ntry, kdiag, k1step, kover/ 5, 1, 1, 1 /

        data ifhmax, ifjac, nonneg, knewt, ksener/ 1, 1, 0, 0, 0 /
        data icall, mcall / 0, 0 /
        data lud, luw /  0, 0 /
	data iwarn /0/
	data ifst/0/

! ------- Update ci_stats
	if(map==0) then
	   call routine_start(i0_rmap2)
	else
	   call routine_start(i1_rmap2)
	endif

!=======================  preliminaries  ==========================

!---------  count calls in units and millions  --------------------

        icall = icall + 1
        if( icall .ge. 1000000 ) then
           mcall = min( mcall + 1, 1000000 )
           icall = 0
        endif

!---------  on first call allocate work arrays  --------------------

	if( ifst == 0 ) then
	   ifst = 1
	   liwdas = 20 + (ns+1) * (ns+2)
	   lrwdas = 40+10*(ns+1)*(ns+2)+(ns+1)**2
	   allocate( iwork(liwdas) )
	   allocate( rwork(lrwdas) )
	endif

!--------- open file for ddasac diagnostic output  -----------------

        if( lud == 0 .and. ddalog == 1) then
           call isat_lu( lud )
           call isat_lu( luw )

	     blank = repeat(' ',30) !  generate file name:  ddasacd_P.op
	     head  = blank
	     head  = 'ddasacd'
	     tail  = blank
	     tail  = 'op'
	     call isat_file_name( head, -1, idproc, tail, name )
	     open( lud , file = name )

	     head  = blank
	     head  = 'ddasacw'
	     call isat_file_name( head, -1, idproc, tail, name )
           open( luw , file = name )
       endif

!------  set initial condition   -------------------------------

	nover = 0
	ifxh0 = 0
	tleft = dt
	h0set = 1.1 * dt

!===========  ddasac restarts (with original i.c.) starts here  ======

10	continue
	do i   = 1, nvar
	   z(i,1)    = max( z0(i), 0.d0 )
	   zstart(i) = z(i,1)
	end do

!--------  set initial sensitivities: identity for sensitivities w.r.t.
!    the initial state; zero w.r.t. pressure       ---------------------

	nspar = npar - 1

        if( map == 1 ) then
           do i   = 1, nvar
              do j   = 2, nspar + 1
                 z(i,j)   = 0.d0
	      end do
              z(i,i+1) = 1.d0
	   end do

           do i = 1, nvar
              ipar(i) = i
	      rpar(i) = zstart(i)
	   end do
        endif

!---------  set pressure (set to rpar(nspar) in cires5 etc.)

	if( nspar /= nvar+1 ) then
	   write(lu_err,*)
     1	      'rmap2: indexes badly set for pressure in cires5'
	   write(lu_err,*)'rmap2: nspar, nvar+1 = ', nspar, nvar+1
	   call isat_abort('ci_rmap',1)
	endif

	ipar(nspar) = nspar
	rpar(nspar) = press

!========  ddasac restarts (with existing state as i.c.) start here ======

!------------  prepare for call to ddasac  ------------------------------

35	continue
	iwork(1:18) = 0
	info(1:18)  = 0

	rwork(1:44) = 0.d0

        if( kdiag == 3 ) k1step = 1

        info(2)  = 1
        info(3)  = k1step

	ifjac   = m_jac

	if (radiation) then
         ifjac  = 0
	endif

! VH: 12/08/2010: set ifjac =0, if user_rate = .true., becuase g_cidzdt
! is not yet implemented for 'usrate' routine
	if (user_rate) then
	   ifjac = 0
	endif

	info(5)  = ifjac

! if radiation is included, the Jacobian is computed by dividend difference

        info(7)  = ifhmax
        info(10) = nonneg
        info(11) = 1
        info(12) = map * nspar
        info(14) = 1
        info(15) = ksener
        info(18) = knewt

        tstart   = 0.
        tout     = tleft

	rwork(2) = tout
        rwork(3) = 1.d0

!  fix initial time step size h0?

	if( ifxh0 == 1 ) then
	   info(8)  = 1
	   rwork(3) = h0set
	   ifxh0    = 0
	endif

!  set tolerances

	varef(1:ns) = phiref

	varef(ns+1) = tempref

	atol(1:nvar) = atolc * varef(1:nvar)
	rtol(1:nvar) = rtolc

!  additional info for sensitivities ( map = 1 )

        if( map == 1 ) then
	   paref(1:nvar) = varef(1:nvar)

	   paref(nspar) = pressref

	   ij   = nvar
           do j = 1, nspar
	   do i = 1, nvar
	      ij   = ij + 1
              atol(ij) = atols * varef(i) / paref(j)
              rtol(ij) = rtols
	   end do
	   end do
        endif

!------  loop over tries at integration (500 ddasac steps) -----------

	itry = 0
100	continue
	itry = itry + 1

!-- on repeat steps call ddasac0 to write out input to ddasac5

	if( nover == 1  .and. ddalog == 1) then

          if (ifjac == 0) then
	    call ddasac0( tstart, tout, nvar, z, zp, rtol, atol,
     1        info, rwork, lrwdas, iwork, liwdas, rpar, ipar,
     2        idid, luw, icall, cires5, ciesub5, cijac5, cibsub5 )
	   else

	     call ddasac0adf( tstart, tout, nvar, z, zp, rtol, atol,
     1        info, rwork, lrwdas, iwork, liwdas, rpar, ipar,
     2        idid, luw, icall, cires5, ciesub5, cibsub5 )!, cijac5adf
	    endif
	endif

!------- regular ddasac substeps start here ---------------------

110	continue

      if (ifjac == 0 ) then
	call ddasac5( tstart, tout, nvar, z, zp, rtol, atol,
     1        info, rwork, lrwdas, iwork, liwdas, rpar, ipar,
     2        idid, luw, 0, cires5, ciesub5, cijac5, cibsub5 )
      else
        call ddasac5( tstart, tout, nvar, z, zp, rtol, atol,
     1        info, rwork, lrwdas, iwork, liwdas, rpar, ipar,
     2        idid, luw, 0, cires5, ciesub5, cibsub5 )!, cijac5adf
       endif
	nsteps = iwork(11)
	nfuncs = iwork(12)
	tlast  = rwork(4)

	zmin     = z(1,1)
	do i = 1, ns
	   zmin     = min( zmin, z(i,1) )
	end do

!---------------  diagnostic information  (ddasac.op)  ------------

	if( kdiag == 3  .or.
     1	  ( kdiag == 2  .and.  idid >= 2 )  .or.
     2	  ( kdiag == 1  .and.  idid .lt. 0 )  .or.
     3	  ( kdiag == 0  .and.  idid .lt. -2 ) .or.
     4	    nover == 1 ) then

	   iordnx = iwork(7)
	   iordls = iwork(8)
	   njacs  = iwork(13)
	   nerrfl = iwork(14)
	   ncnvfl = iwork(15)

	   dtnext = rwork(3)
	   dtlast = rwork(7)

	   if(ddalog > 0) then
	      write(lud, 600) icall, itry, idid, iordnx, iordls,
     1	      nsteps, nfuncs, njacs, nerrfl, ncnvfl, tstart,
     1	      tout, tlast, dtnext, dtlast, zmin
	      call isat_flush( lud )
	   endif
600	   format(i9, 9i5, 1p,20e13.4)

!--------------  diagnostic output on repeat steps  (ddasacd.op) -----

	   if( nover == 1 .and. ddalog > 0 ) then
	      write(luw,*) 'repeated step: icall=', icall
	      write(luw,610) nsteps, tlast, (z(i,1),i=1,nvar)
	      call isat_flush( luw )
610	      format(i5,1p,100e13.4)
	    endif
	endif

!=============  determine outcome of integration and act accordingly =========

!-----  success

	if( idid >= 2 ) then
	   if( nover == 1 .and. ddalog > 0 ) then
	      write(luw,*) 'success after failure: ',
     1	                       'icall=', icall
	      call isat_flush(luw)
	   endif
	   go to 195
	endif

!-----  no problems, but more steps required

	if( idid >= 0  .and.  zmin >= -tolneg ) go to 110

!-----  negative variable; set to zero and restart ddasac at current time

	if( tlast .lt. tleft  .and.  zmin .lt. -tolneg ) then
	   iwarn = iwarn + 1
	   tleft = tleft - tlast
	   ifxh0 = 0

	   if( iwarn <= 100 .and. ddalog > 0 ) then
              write(luw,*)'rmap2: neg. z '
	      write(luw, 620) iwarn, icall, map, nsteps
	      write(luw, 621) zmin, dt, tleft, tlast
620           format(3x,'iwarn,icall,map,nsteps= ',i4,i9,i4,i4)
621           format(3x,'zmin,dt,tleft,tlast   = ',1p,10e13.4)
	      call isat_flush(luw)
	   endif

	   do i = 1, ns
	      z(i,1)   = max( z(i,1), 0.d0 )
	   end do
	   go to 35
	endif

!-----  keep trying if tolerances reset

	if( idid == -2 ) then
	   if( kdiag .ne. 0 .and. ddalog > 0 ) then
	      write(luw,*)'rmap2: warning, tolerances ',
     1	                     'reset, icall, itry=', icall, itry
	      call isat_flush(luw)
	   endif
	   go to 110
	endif

!-----  failure to find initial step size

	   if( idid == -33  .and.  nsteps == 0  .and.
     1	       nfuncs >= 100  .and.  ifxh0 == 0  ) then

	      call cidzdt( z, p, dzdt, densx )
	      h0set  = 1.1 * tout
	      ifxh0  = 1

	      do i = 1, nvar
	         if( z(i,1) + dzdt(i) * h0set .lt. 0. )
     1	            h0set = - z(i,1) / dzdt(i)
     	      end do

	      if(ddalog > 0) write(luw,*)'rmap2: icall=', icall,
     1	                     'fixing h0: h0, dt=', h0set, dt
	      if(ddalog > 0) call isat_flush(luw)
	      go to 10
	   endif

!-----  500 steps without reaching end: try again?

	if( idid == -1  .and.  itry .lt. ntry ) then
	   info(1) = 1
	   if(kdiag.ne.0 .and. ddalog>0)
     1        write(luw,*)'rmap2: try another 500,',
     1	                        'icall, itry=',icall, itry
	   if(ddalog>0) call isat_flush(luw)
	   go to 100
	endif

!=================  irrecoverable failure:  repeat step with diagnostics
!    turned on, call rmdump, then quit  ================================

	if( kover == 1 .and.  nover == 0 ) then
	   if(ddalog>0) then
	      write(luw,*)
     1	        'rmap2: failure; repeating step with diagnostics'
	      write(luw,*)'rmap2: icall, idid=',icall, idid
	      call isat_flush(luw)
	   endif

	   nover    = 1
	   k1step0  = k1step
	   k1step   = 1
	   kdiag0   = kdiag
	   kdiag    = 3
	   ifxh0    = 0
	   tleft    = dt

	   z(1:nvar,1)   = zstart(1:nvar)

	   go to 10
	endif

!-----  stop
	if(ddalog>0) then
	   write(luw,*)'rmdump called from rmap2: ',
     1	            'icall=', icall
	   call rmdump( map, nvar, npar, z0, press, dt )

	   write(lu_err,*)'rmap2: ddasac failure: see file ddasacw.op',
     1	             ' icall=', icall
	   write(luw,*)'rmap2: ddasac failed, idid=', idid
	   write(luw,*)'rmap2: icall, map, itry=', icall, map, itry
	   write(luw,900)tstart,tout
	   write(luw,900)zstart
	   write(luw,900)(z(i,1),i=1,nvar)
 900	   format((1p,5e13.4))
	endif
	call isat_abort('ci_rmap',2)

!----------------------  successful integration  ---------------

!  set sensitivity w.r.t. dt

195	continue
	z(1:nvar,npar+1) = zp(1:nvar,1)

! ------- Update ci_stats
	if(map==0) then
	   call routine_stop(i0_rmap2)
	else
	   call routine_stop(i1_rmap2)
	endif

	return
	end subroutine rmap2

!============================================================================

	subroutine rmap1_test( idtab, nx, x, nf, nh, nhd, usrfgh,
     1	             iusr, rusr, info, rinfo, fa, ga, ha, stats )

!  Routine to test the derivative matrix df/dx produced by rmap1.
!  The values of df/dx from rmap1 are compared to those obtained
!  from a finite-difference approximation.

!  This routine is called from cirnx7 if  rmap1t /= 0.
!  The max abs normalized difference (diffmax) is output on rmt1.op
!       if it exceeds rmt1_op.
!  Full output is generated for the call  icall = rmapt1.

	implicit none

	integer, intent(in)    :: idtab, nx, nf, nh, nhd, info(20)
	integer                :: iusr(*)

	real(k_dp), intent(in) :: x(nx), rinfo(20+nx+nf)
	real(k_dp)             :: fa(nf), ha(nhd), stats(100)
	real(k_dp), target     :: ga(nf,nx)
	real(k_dp)             :: rusr(*)

!---------------------------------------------------------------------------

	interface
	subroutine usrfgh( need, nx, x, nf, nh, iusr, rusr, fa, ga, ha )
        integer :: need(3),  nx, nf, nh, iusr(*)
	  real(kind(1.d0)) :: x(nx), rusr(*), fa(nf), ga(nf,nx), ha(*)
	end subroutine usrfgh
	end interface

!---------------------------------------------------------------------------

	integer :: need(3), i, j
	real(k_dp) :: x0(nx), xp(nx), dx, dfdx0(nf,nx),
     1	  f0(nf), fp(nf), dfdxfd(nf,nx), fnorm(nf),
     2	  xnorm(nx), dfdx0sc(nf,nx), dfdxfdsc(nf,nx),
     3    diff, diffmax, diffsc, xm(nx), fm(nf), dfdxp(nf,nx), h0(nh)
	integer, save :: lu = 0, icall = 0
	real(k_dp), save :: aeps=1.d-6, reps=1.d-3

!  open file rmt1.op for output

	if( lu == 0 ) then
	   call isat_lu(lu)
	   open( lu, file = 'rmt1.op', action='write', status='replace')
	endif

	icall = icall + 1

!  obtain f0 and dfdx0 from usrfgh

	x0        = x
	need(1:2) = 1
	need(3)   = 0

	call usrfgh( need, nx, x0, nf, nh, iusr, rusr, f0, dfdx0, h0 )

!  obtain dfdxfd

	need(2) = 0
	do i = 1, nx

	xp  = x0
	xm  = x0

	xp(i) = x0(i) * ( 1.d0 + reps ) + aeps
	xm(i) = 2*x0(i) - xp(i)
	if( xm(i) <= 0.d0 ) xm(i) = x0(i)
	dx    = xp(i) - xm(i)

	call usrfgh( need, nx, xp, nf, nh, iusr, rusr, fp, dfdxp, h0 )
	call usrfgh( need, nx, xm, nf, nh, iusr, rusr, fm, dfdxp, h0 )

	do j = 1, nf
	   dfdxfd(j,i) = (fp(j)-fm(j)) / dx
	end do
	end do

!  form scaled dfdx

	xnorm = rinfo(21:20+nx)
	fnorm = rinfo(20+nx+1:20+nx+nf)

	diffmax = -1.d0
	do i = 1, nx
	do j = 1, nf
	   dfdx0sc(j,i) = dfdx0(j,i)*xnorm(i)/fnorm(j)
	   dfdxfdsc(j,i) = dfdxfd(j,i)*xnorm(i)/fnorm(j)
	   diffmax = max( diffmax, abs(dfdx0sc(j,i)-dfdxfdsc(j,i)) /
     1	                           max(1.d0, abs(dfdx0sc(j,i)) ) )
	end do
	end do

!--- output

	if( diffmax >= rmt1_op ) then
	   write(lu,600) icall, diffmax
600	   format(i8,1p,e13.4)
	   call isat_flush(lu)
	endif

	if(  nint(rmap1t) /= icall ) return
	write(lu,*)' '
	write(lu,*)'Output from rmap1_test for icall= ', icall
	write(lu,*)' '
	write(lu,*)' i        x0(i)      f0(i)   f0(i)-x0(i)'
	write(lu,*)' '
	do i = 1, min(nx,nf)
	   write(lu,500) i, x0(i), f0(i), f0(i)-x0(i)
	end do
500	format(i4,1p,3e13.4)

!  write full ouput if diffmax > rmt1_op

	if( diffmax < rmt1_op ) return

	do j = 1, nf
	write(lu,*)' '
	write(lu,605) j, f0(j), f0(j)/fnorm(j)
605	format('j, f0, f0norm = ', i4, 1p,2e13.4)
	write(lu,*)' '
	write(lu,*)'   dfdx0        dfdxfd      diff     i    ',
     1	           '  dfdx0sc      dfdxfdsc    diffsc'

	do i = 1, nx
	diff = dfdx0(j,i) - dfdxfd(j,i)
	diffsc = diff*xnorm(i)/fnorm(j)
	write(lu,610) dfdx0(j,i), dfdxfd(j,i), diff, i,
     1	  dfdx0sc(j,i), dfdxfdsc(j,i), diffsc
610	format(1p,3e12.3,i4,1p,3e12.4)
	end do
	end do

	call isat_flush(lu)

	stop
	end subroutine rmap1_test

!============================================================================

	subroutine rmap2_test( map, nz, nv, zin, press, dt, zr )

!  Routine to test the derivative matrix dz/dv produced by rmap2.
!  The values of dz/dv from rmap2 are compared to those obtained
!  from a finite-difference approximation.
!  z = {phi(dt), T(dt)}
!  v = {phi(0), T(0), press, dt}

!  This routine is called from cirmap1 if  rmap2t /= 0.
!  The max abs normalized difference (diffmax) is output on rmt2.op
!       if it exceeds rmt2_op.
!  Full output is generated for the call  icall = rmapt2.

	implicit none

	integer, intent(in)    :: map, nz, nv
	real(k_dp), intent(in) :: zin(nz), press, dt
	real(k_dp), target     :: zr(nsp1,nsp4)

	integer :: mapx, i, j
	real(k_dp) :: v0(nsp1), p0, dt0, vp(nsp1), pp, dtp, dv,
     1	  z0(nsp1), zp(nsp1), dzdvfd(nsp1,nsp3), znorm(nsp1),
     2	  vnorm(nsp3), dzdv0sc(nsp1,nsp3), dzdvfdsc(nsp1,nsp3),
     3    diff, diffmax, diffsc, vm(nsp1), pm, dtm, zm(nsp1)
	integer, save :: lu = 0, icall = 0
	real(k_dp), pointer :: dzdv0(:,:)
	real(k_dp), save :: aepsphi=1.d-6, repsphi=1.d-3,
     1	                    aepstemp=0.d0, repstemp=1.d-3,
     2	                    aepspress=0.d0, repspress=1.d-3,
     3	                    aepsdt=0.d0, repsdt=1.d-3

!  open file rmt2.op for output

	if( lu == 0 ) then
	   call isat_lu(lu)
	   open( lu, file = 'rmt2.op', action='write', status='replace')
	endif

	icall = icall + 1

!  obtain z0 and dzdv0 from rmap2

	v0        = zin
	p0        = press
	dt0       = dt
	mapx      = 1

	call rmap2( mapx, nsp1, nsp3, v0, p0, dt0, zr )

	z0 = zr(1:nsp1,1)
	dzdv0 => zr(1:nsp1,2:nsp4)

!  obtain dzdvfd

	mapx = 0
	do i = 1, nsp3

	pp  = p0
	dtp = dt0
	vp  = v0

	pm  = p0
	dtm = dt0
	vm  = v0

	if( i <= ns ) then
	   vp(i) = v0(i) * ( 1.d0 + repsphi ) + aepsphi
	   vm(i) = 2*v0(i) - vp(i)
	   if( vm(i) <= 0.d0 ) vm(i) = v0(i)
	   dv    = vp(i) - vm(i)
	elseif( i == nsp1 ) then
	   vp(i) = v0(i) * ( 1.d0 + repstemp ) + aepstemp
	   vm(i) = 2*v0(i) - vp(i)
	   if( vm(i) <= 0.d0 ) vm(i) = v0(i)
	   dv    = vp(i) - vm(i)
	elseif( i == nsp2 ) then
	   pp  = p0 * ( 1.d0 + repspress ) + aepspress
	   pm  = 2*p0 - pp
	   dv  = pp - pm
	elseif( i == nsp3 ) then
	   dtp = dt0 * ( 1.d0 + repsdt ) + aepsdt
	   dtm = 2*dt0 - dtp
	   dv  = dtp - dtm
	endif

	prc = pp
	dtc = dtp
	call rmap2( mapx, nsp1, nsp3, vp, pp, dtp, zp )

	prc = pm
	dtc = dtm
	call rmap2( mapx, nsp1, nsp3, vm, pm, dtm, zm )
	dtc = dt0
	prc = p0

	do j = 1, nsp1
	   dzdvfd(j,i) = (zp(j)-zm(j)) / dv
	end do
	end do

!  form scaled dzdv

	if( dtref == 0.d0 ) dtref = dt
	znorm(1:ns) = phiref
	znorm(nsp1) = tempref

	vnorm(1:ns) = phiref
	vnorm(nsp1) = tempref
	vnorm(nsp2) = pressref
	vnorm(nsp3) = dtref

	diffmax = -1.d0
	do i = 1, nsp3
	do j = 1, nsp1
	   dzdv0sc(j,i) = dzdv0(j,i)*vnorm(i)/znorm(j)
	   dzdvfdsc(j,i) = dzdvfd(j,i)*vnorm(i)/znorm(j)
	   diffmax = max( diffmax, abs(dzdv0sc(j,i)-dzdvfdsc(j,i)) /
     1	                           max(1.d0, abs(dzdv0sc(j,i)) ) )
	end do
	end do

!--- output

	if( diffmax >= rmt2_op ) then
	   write(lu,600) icall, diffmax
600	   format(i8,1p,e13.4)
	   call isat_flush(lu)
	endif

	if(  nint(rmap2t) /= icall ) return
	write(lu,*)' '
	write(lu,*)'Output from rmap2_test for icall= ', icall
	write(lu,*)' '
	write(lu,*)' i        v0(i)      z0(i)   z0(i)-v0(i)'
	write(lu,*)' '
	do i = 1, nsp1
	   write(lu,500) i, v0(i), z0(i), z0(i)-v0(i)
	end do
	write(lu,500) nsp2, p0
	write(lu,500) nsp3, dt0
500	format(i4,1p,3e13.4)


!  write full ouput if diffmax > rmt2_op

	if( diffmax < rmt2_op ) return

	do j = 1, nsp1
	write(lu,*)' '
	write(lu,605) j, z0(j), z0(j)/znorm(j)
605	format('j, z0, z0norm = ', i4, 1p,2e13.4)
	write(lu,*)' '
	write(lu,*)'   dzdv0        dzdvfd      diff     i    ',
     1	           '  dzdv0sc      dzdvfdsc    diffsc'

	do i = 1, nsp3
	diff = dzdv0(j,i) - dzdvfd(j,i)
	diffsc = diff*vnorm(i)/znorm(j)
	write(lu,610) dzdv0(j,i), dzdvfd(j,i), diff, i,
     1	  dzdv0sc(j,i), dzdvfdsc(j,i), diffsc
610	format(1p,3e12.3,i4,1p,3e12.4)
	end do
	end do

	call isat_flush(lu)

	stop
	end subroutine rmap2_test

!============================================================================

        subroutine rmdump( map, nvar, npar, z0, press, dt )

!  routine to dump input parameters to rmap2
!
	implicit none

	integer, intent(in)   :: map, nvar, npar
	real(k_dp), intent(in):: z0(nvar), press, dt

	integer, save :: lu=0
	integer       :: i

	if( lu == 0 ) then
	   call isat_lu(lu)
	   open( lu, file='rmdump.op', form = 'unformatted' )
	endif

	write(lu) map, nvar, npar
	write(lu) press, dt
	write(lu) (z0(i),i=1,nvar)

	return
	end subroutine rmdump

!============================================================================

       subroutine ddasac0adf( tstart, tout, nvar, z, zp, rtol, atol,
     1        info, rwork, lrwdas, iwork, liwdas, rpar, ipar,
     2        idid, luw, icall, cires5, ciesub5, cibsub5 )!, cijac5adf

!  routine to write out (on luw) the calling arguments to ddasac5

	implicit none

	integer, intent(in)  :: liwdas, nvar, info(18), lrwdas,
     1	   iwork(liwdas), ipar(:), luw, icall

     	real(k_dp), intent(in)  :: tstart, tout, z(*), zp(*), rtol(:),
     1	   atol(:),  rwork(lrwdas), rpar(:)

     	integer    :: idid, i, neq, n1, n2, nspar, ieform, ires
	real(k_dp) :: rtmin, rtmax, atmin, atmax

         external cires5, ciesub5, cijac5, cibsub5   !,cijac5adf
	if( luw < 0 ) call isat_abort('ddasac0',1,mess=
     1                       'negative luw = ', isv=luw)

	write(luw,*)' '
	write(luw,*)'ddasac0: start of output of calling arguments'
	write(luw,*)'for rmap2 icall=', icall
	write(luw,*)' '

	write(luw,*)'tstart = ', tstart
	write(luw,*)'tout   = ', tout
	write(luw,*)'nvar   = ', nvar
	write(luw,*)'lrwdas = ', lrwdas
	write(luw,*)'liwdas = ', liwdas
	write(luw,*)'luw    = ', luw

	write(luw,*)' '
	write(luw,*)'info:'
	write(luw,*)' 1: other than initial call?   ', info(1)
	write(luw,*)' 2: tolerances as arrays?      ', info(2)
	write(luw,*)' 3: one step at a time?        ', info(3)
	write(luw,*)' 4: force stop at tout?        ', info(4)
	write(luw,*)' 5: jacobian provided?         ', info(5)
	write(luw,*)' 6: banded jacobian?           ', info(6)
	write(luw,*)' 7: max dt specified?          ', info(7)
	write(luw,*)' 8: sign of initial step?      ', info(8)
	write(luw,*)' 9: maxord specified?          ', info(9)
	write(luw,*)'10: non-negativity?            ', info(10)
	write(luw,*)'11: consistent zp provided?    ', info(11)
	write(luw,*)'12: no. of sensitivities       ', info(12)
	write(luw,*)'13: non-unit matrix?           ', info(13)
	write(luw,*)'14: FD for B matrix?           ', info(14)
	write(luw,*)'15: error inc. sensitivities?  ', info(15)
	write(luw,*)'16: euclidean norm?            ', info(16)
	write(luw,*)'17: do not cross tout?         ', info(17)
	write(luw,*)'18: temp1 treatment            ', info(18)

	write(luw,*)' '
	write(luw,*)'initial conditions, rtol, atol'

	do i = 1, nvar
	   write(luw,1) i, z(i), rtol(i), atol(i)
1	   format(i4,1p,10e16.7)
	end do

	write(luw,*)' '
	write(luw,*)'initial species compared to rxnip'

	do i = 1, nvar-1
	   write(luw,1) i, z(i), rxnip(i)
	end do

	nspar = info(12)
	if( nspar > 0 ) then
	   write(luw,*)' '
	   write(luw,*)'ipar, rpar'
	   do i = 1, nspar
	      write(luw,3) i, ipar(i), rpar(i)
3	      format(2i4,1p,10e16.7)
	   end do

	   n1 = nvar+1
	   n2 = n1 + nvar*nspar - 1
	   write(luw,*)' '
	   write(luw,*)'initial sensitivities: n1, n2=', n1, n2
	   write(luw,2)(z(i),i=n1,n2)
2	   format((1p,8e10.2))
	endif

!  check range of tolerances

	neq   = nvar*(1+nspar)
	rtmin =  1.d30
	atmin =  1.d30
	rtmax = -1.d30
	atmax = -1.d30
	do i = 1, neq
	   rtmin    = min( rtmin, rtol(i) )
	   atmin    = min( atmin, atol(i) )
	   rtmax    = max( rtmax, rtol(i) )
	   atmax    = max( atmax, atol(i) )
	end do

	write(luw,*)' '
	write(luw,*)'rtmin, rtmax, atmin, atmax = '
	write(luw,4) rtmin, rtmax, atmin, atmax
4	format(1p,10e13.4)

	write(luw,*)' '
	if( info(4) == 1 ) write(luw,*) 'rwork(1)=tstop=',rwork(1)
	if( info(7) == 1 ) write(luw,*) 'rwork(2)=hmax =',rwork(2)
	write(luw,*) 'rwork(3)=h0   =',rwork(3)

	write(luw,*)' '
	write(luw,*)'about to call cires5 '
	ieform = 0
	call cires5( tout, nvar, z, zp, rpar, ipar, ieform, ires )
	write(luw,*)'completed call to cires5 '
	write(luw,*)'i, z(i), dzdt(i)'

	do i = 1, nvar
	   write(luw,1) i, z(i), zp(i)
	end do

	write(luw,*)' '
	write(luw,*)'ddasac0: end of output of calling arguments'
	write(luw,*)' '

	call isat_flush(luw)

	return
	end subroutine ddasac0adf

!============================================================================

	subroutine ddasac0( tstart, tout, nvar, z, zp, rtol, atol,
     1        info, rwork, lrwdas, iwork, liwdas, rpar, ipar,
     2        idid, luw, icall, cires5, ciesub5, cijac5, cibsub5 )

!  routine to write out (on luw) the calling arguments to ddasac5

	implicit none

	integer, intent(in)  :: liwdas, nvar, info(18), lrwdas,
     1	   iwork(liwdas), ipar(:), luw, icall

     	real(k_dp), intent(in)  :: tstart, tout, z(*), zp(*), rtol(:),
     1	   atol(:),  rwork(lrwdas), rpar(:)

     	integer    :: idid, i, neq, n1, n2, nspar, ieform, ires
	real(k_dp) :: rtmin, rtmax, atmin, atmax

       external cires5, ciesub5, cijac5, cibsub5
 	if( luw < 0 ) call isat_abort('ddasac0',1,mess=
     1                       'negative luw = ', isv=luw)

	write(luw,*)' '
	write(luw,*)'ddasac0: start of output of calling arguments'
	write(luw,*)'for rmap2 icall=', icall
	write(luw,*)' '

	write(luw,*)'tstart = ', tstart
	write(luw,*)'tout   = ', tout
	write(luw,*)'nvar   = ', nvar
	write(luw,*)'lrwdas = ', lrwdas
	write(luw,*)'liwdas = ', liwdas
	write(luw,*)'luw    = ', luw

	write(luw,*)' '
	write(luw,*)'info:'
	write(luw,*)' 1: other than initial call?   ', info(1)
	write(luw,*)' 2: tolerances as arrays?      ', info(2)
	write(luw,*)' 3: one step at a time?        ', info(3)
	write(luw,*)' 4: force stop at tout?        ', info(4)
	write(luw,*)' 5: jacobian provided?         ', info(5)
	write(luw,*)' 6: banded jacobian?           ', info(6)
	write(luw,*)' 7: max dt specified?          ', info(7)
	write(luw,*)' 8: sign of initial step?      ', info(8)
	write(luw,*)' 9: maxord specified?          ', info(9)
	write(luw,*)'10: non-negativity?            ', info(10)
	write(luw,*)'11: consistent zp provided?    ', info(11)
	write(luw,*)'12: no. of sensitivities       ', info(12)
	write(luw,*)'13: non-unit matrix?           ', info(13)
	write(luw,*)'14: FD for B matrix?           ', info(14)
	write(luw,*)'15: error inc. sensitivities?  ', info(15)
	write(luw,*)'16: euclidean norm?            ', info(16)
	write(luw,*)'17: do not cross tout?         ', info(17)
	write(luw,*)'18: temp1 treatment            ', info(18)

	write(luw,*)' '
	write(luw,*)'initial conditions, rtol, atol'

	do i = 1, nvar
	   write(luw,1) i, z(i), rtol(i), atol(i)
1	   format(i4,1p,10e16.7)
	end do

	write(luw,*)' '
	write(luw,*)'initial species compared to rxnip'

	do i = 1, nvar-1
	   write(luw,1) i, z(i), rxnip(i)
	end do

	nspar = info(12)
	if( nspar > 0 ) then
	   write(luw,*)' '
	   write(luw,*)'ipar, rpar'
	   do i = 1, nspar
	      write(luw,3) i, ipar(i), rpar(i)
3	      format(2i4,1p,10e16.7)
	   end do

	   n1 = nvar+1
	   n2 = n1 + nvar*nspar - 1
	   write(luw,*)' '
	   write(luw,*)'initial sensitivities: n1, n2=', n1, n2
	   write(luw,2)(z(i),i=n1,n2)
2	   format((1p,8e10.2))
	endif

!  check range of tolerances

	neq   = nvar*(1+nspar)
	rtmin =  1.d30
	atmin =  1.d30
	rtmax = -1.d30
	atmax = -1.d30
	do i = 1, neq
	   rtmin    = min( rtmin, rtol(i) )
	   atmin    = min( atmin, atol(i) )
	   rtmax    = max( rtmax, rtol(i) )
	   atmax    = max( atmax, atol(i) )
	end do

	write(luw,*)' '
	write(luw,*)'rtmin, rtmax, atmin, atmax = '
	write(luw,4) rtmin, rtmax, atmin, atmax
4	format(1p,10e13.4)

	write(luw,*)' '
	if( info(4) == 1 ) write(luw,*) 'rwork(1)=tstop=',rwork(1)
	if( info(7) == 1 ) write(luw,*) 'rwork(2)=hmax =',rwork(2)
	write(luw,*) 'rwork(3)=h0   =',rwork(3)

	write(luw,*)' '
	write(luw,*)'about to call cires5 '
	ieform = 0
	call cires5( tout, nvar, z, zp, rpar, ipar, ieform, ires )
	write(luw,*)'completed call to cires5 '
	write(luw,*)'i, z(i), dzdt(i)'

	do i = 1, nvar
	   write(luw,1) i, z(i), zp(i)
	end do

	write(luw,*)' '
	write(luw,*)'ddasac0: end of output of calling arguments'
	write(luw,*)' '

	call isat_flush(luw)

	return
	end subroutine ddasac0

!============================================================================

	end module ci_rmap
