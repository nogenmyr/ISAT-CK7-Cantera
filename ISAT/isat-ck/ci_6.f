!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

	module ci_6

	use ci_dat
	use ci_dat6
	use ci_ck  
	use isat_val, only: isat_lic, coderet
	use ci_rmap
	use streams_mod

	contains

!=========================================================================

        subroutine cinit6

!  chemistry interface initialization routine for  modeci = 6 and 7
!     which use Cantera.

!  files:
!       streams.in - input file of stream information.
!BEGEXTRACT

!      Specification of the file  streams.in  for modeci = 6 and 7

!         1st record - modeci
!         2nd record - nstr, nsin  
!                      nstr - number of streams (integer)
!                      nsin - number of non-zero species (integer)
!         subsequent - k, p, T, c(1), c(2),...,c(nsin)
!                      k = 1 if stream is as stated (integer)
!                      k = 2 if stream is an equilibrium mixture with the
!                            same elemental composition, pressure and
!                            enthalpy as that stated.
!                      p    - pressure in atm  (real)
!                      T    - temperature in K (real)
!                      c(i) - composition in relative volume/mole units (real)


!ENDEXTRACT

!  comment: 
!        c(i) does not need to be normalized
!       
	implicit none

	integer :: retcode

	if(op_rank) write(luout,*)' '

	if( myrank == 0 ) then
	   call isat_lic( luout, retcode )
	   if( retcode /= coderet ) call isat_abort('cinit6',0, 
     1	       mess='Invalid license code' )
	endif

!----------------  set defaults, and read changes from ci.nml

!  set all parameters to zero if they have not been set in ciparam
	if( .not.ciparam_called ) then
	   ci_info_n  = 0
	   ci_rinfo_n = 0.d0
	   info_n     = 0
	   rinfo_n    = 0.d0
	endif

	call ci_info_set

!--------------  read number of streams and non-zero species  ----------------

	if( luin > 0 ) then
	   if( strm_format == 1 ) then
	      read(luin,*,err=140,end=142) nstr, nsin
	      if( nstr < 0 ) call isat_abort( 'cinit6',3,
     1	                  mess='bad value of  nstr', isv = nstr )
	      if(op_rank) write(luout,*)'number of streams, nstr = ', nstr
	   else
	      nstr = StrmDat%nStream
	   endif
     	else
	   nstr = 0
	   nsin = 0
	endif
     	   

!--------------  initialize Cantera, and determine number of species, etc

	!call cickin

!---------------  perform all other initialization

	call ci_init6

	return

!-------------  error conditions

140	call isat_abort('cinit6',1, mess='error reading nstr and nsin' )
142	call isat_abort('cinit6',2, mess='hit end reading nstr and nsin' )

	end subroutine cinit6

!=========================================================================

! subroutines ci_info_set and cickin have been moved to ci_ck.f

!=========================================================================

	subroutine ci_init6

!  routine to get basic data from Cantera and to initialize for modeci = 6 and 7

!  NOTE:  reference values are obtained from: ci.nml; or if not set
!         there, from streams; or if not from there, set here.

      use ceq_state_m
	implicit none

	logical       :: kerr
	integer       :: rv(ns,nr), ev(ns,ne), kstr(nstr),
     1	                kpress, ist, is, info
	real(k_dp)    ::  phiavg(ns), pstr(nstr),
     1	                tstr(nstr), ystr(ns,nstr), hstr(nstr), 
     2	                rhostr(nstr), hs, tempcp, phi_eq(ns)

	nc    = nsp1
	nfull = nsp4

	call ci_alloc
!  obtain  reaction vectors

	call ctnu( ns, rv, gas )
	ev = dev

!-----  define names of full composition variables

	cmpsym(1:ns) = snames(1:ns)
	cmpsym(ns+1) = symb_dens
	cmpsym(ns+2) = symb_temp
	cmpsym(ns+3) = symb_press
	cmpsym(ns+4) = symb_enth

!  print out species info

	if(op_rank) call op_species

	if( nstr == 0 ) go to 700

!----- treatment with streams -------------

	if ( strm_format == 1 ) then
	   if( nsin < 1  .or.  nsin > ns ) call isat_abort( 'ci_init6', 1,
     1	         mess = 'bad value of nsin = ', isv = nsin )

	   ccstrm = 0.d0

	   do ist = 1, nstr
	      read(luin,*,err=144,end=146)kstr(ist),pstr(ist),tstr(ist),
     1	                   (ccstrm(is,ist),is=1,nsin)
     	   end do
	else
	   kstr = StrmDat%kstrm
	   pstr = StrmDat%Cvalue(ns+1,:)
	   tstr = StrmDat%Cvalue(ns+2,:)
	   do ist = 1,nstr
	      if(StrmDat%kmole(ist) == 0) then
		 call y2phi( StrmDat%Cvalue(1:ns,ist), amolwt, ns, 
     1		       ccstrm(1:ns,ist) )
	      else
		 ccstrm(1:ns,ist) = StrmDat%Cvalue(1:ns,ist)
	      endif
	   enddo
	endif
!-----  derived stream properties  -----------------------------------------


!  normalize to form specific mole numbers, enthalpy and density
! (note that hstr is enthalpy, not sensible enthalpy)

	kpress = 0
	do ist    = 1, nstr

	   pstr(ist)     = pstr(ist) * patm
	   if( pstr(ist) .ne. pstr(1) ) kpress = max( 1, kpress )

	   call phi2y( ccstrm(1:ns,ist), amolwt, ns, ystr(:,ist) )
	   call y2phi( ystr(:,ist), amolwt, ns, ccstrm(1:ns,ist) )
	   call cthbms( tstr(ist), ystr(:,ist), hstr(ist),gas )

!  determine equilibrium if required: kstr=2 for fixed (p,h); kstr>=3 for fixed (p,T)

	   if( kstr(ist) >= 2 ) then
	      if(op_rank) then
	      write(luout,*)' '
	      write(luout,'(a,a,i3,a,i2)')
     1                 'Starting equilibrium calculation ',
     1	               'for stream',ist, ', kstr = ', kstr(ist)
	      endif

	      if( kstr(ist) == 2 ) then
		 call ceq_state( sys_uc, N=ccstrm(1:ns,ist),
     1                      p_cgs=pstr(ist), N_h=ccstrm(1:ns,ist),
     2                      T_h=tstr(ist), N_eq=phi_eq, T_eq=tstr(ist),
     3                      info=info )
	      else
		 call ceq_state( sys_uc, N=ccstrm(1:ns,ist),
     1                      p_cgs=pstr(ist), T=tstr(ist),
     2                      N_eq=phi_eq, info=info )
	      endif

	      if(op_rank) then
              write(luout,*)' '
	      write(luout,'(a,f8.2,a)')'Equilibrium temperature =',
     1                  tstr(ist),' (K)'
	      endif
	      call phi2y( phi_eq, amolwt, ns, ystr(:,ist) )
	      call y2phi( ystr(:,ist), amolwt, ns, ccstrm(1:ns,ist) )
	      call cthbms( tstr(ist), ystr(:,ist),hstr(ist),gas)
	   endif

	   call ctrhoy( pstr(ist), tstr(ist), ystr(:,ist),
     1		rhostr(ist), gas )

!  set sensible enthalpy

	   call h2hs( hstr(ist), href, ccstrm(1:ns,ist), ns, hs )
	   ccstrm(nc,ist)= hs

	   dptstr(1,ist) = rhostr(ist)
	   dptstr(2,ist) = pstr(ist)
	   dptstr(3,ist) = tstr(ist)
	end do

!  obtain  tempref and pressref from streams

	if( tempref  == 0.d0 ) tempref  = maxval( tstr )
	if( pressref == 0.d0 ) pressref = maxval( pstr )

!  check assumption of constant pressure

	if( const_pr  .and.  kpress /= 0 ) then
	   if(op_rank) then
	      write(luout,*)' '
	      write(luout,*)'===WARNING===: const_pr = .true., but',
     1	                 ' streams have different pressures.'
	      write(luout,*)'Using variable pressure'
	   endif
	   const_pr = .false.
	endif

!  print stream info
	if(op_rank) call op_streams

!---------  end of streams considerations

700	continue

!----------- define reference values
!  For isatab:  x = phi(0), hs(0), press,    dt,   and these are scaled by:
!                   phiref, hsref, pressref, dtref
!               f = phi(dt), hs(dt), temp(dt),   which are scaled by
!                   phiref,  hsfref, tempref
!  The reference values can be set in ci.nml

	if( phiref   == 0.d0 ) phiref   = 0.1d0
	if( tempref  == 0.d0 ) tempref  = 1000.d0
	if( pressref == 0.d0 ) pressref = patm

!  define reference specific heat at  tempref   with
!      mole fraction = 1/ns for all species -------

	phiavg = 1.d0/float(ns)

	tempcp = min( tempref, tbadhi )
	call ctcpbs( tempcp, phiavg, cpmref, gas)

	if( hsref == 0.d0 ) hsref = cpmref * tempref

	if(op_rank) then
	   write(luout,*)' '
	   write(luout,*)'Reference values: '
	   write(luout,600) phiref
	   write(luout,601) tempref
	   write(luout,602) pressref
	   write(luout,603) hsref
	endif

	if( hsfref == 0.d0 ) then
	   hsfref = hsref
	else
	   if(op_rank) write(luout,604) hsfref
	endif

	if( dtref /= 0.d0 ) then
	   if(op_rank) write(luout,605) dtref
	endif
	   
600	format('Species:     phiref   = ', 1p,e12.4)
601	format('Temperature: tempref  = ', 1p,e12.4)
602	format('Pressure:    pressref = ', 1p,e12.4)
603	format('Sens. enth.: hsref    = ', 1p,e12.4)
604	format('Sens. enth.: hsfref   = ', 1p,e12.4)
605	format('Time interv: dtref    = ', 1p,e12.4)

!---------------------  set mode_pdt and nxx

	if(op_rank) write(luout,*)' '

	if( .not.const_pr ) then
	   if(op_rank) write(luout,*)'Pressure is taken to be variable.'
	   prc = -huge(1.d0)
	   if( .not.const_dt ) then
	      mode_pdt = 1
	      nxx      = nsp3
	   else
	      mode_pdt = 2
	      nxx      = nsp2
	   endif
	else
	   if( nstr > 0 ) then
	      prc = pstr(1)
	   else
	      prc = pressref
	   endif 

	   if( .not.const_dt ) then
	      mode_pdt = 3
	      nxx      = nsp2
	   else
	      mode_pdt = 4
	      nxx      = nsp1
	   endif
	   if(op_rank) write(luout,710) prc
710	   format('Pressure is taken to be constant, = ', 1p,e13.4)
	endif

	if(op_rank) write(luout,*)' '

	if( const_dt ) then
	   if(op_rank) write(luout,*)'Time step is taken to be constant.'
	else
	   if(op_rank) write(luout,*)'Time step is taken to be variable.'
	endif

	if(op_rank) write(luout,*)' '

	if( radiation ) then
	   if(op_rank) write(luout,*)'Enthalpy equation includes ',
     1	                 'radiative heat loss.'
	else
	   if(op_rank) write(luout,*)'Enthalpy equation does not ',
     1	                 'include radiative heat loss.'
	endif

	if(op_rank) write(luout,*)' '
	if(op_rank) write(luout,*)
     1     'Number of independent ISAT variables, nx = ', nxx


	if(op_rank) call isat_flush(luout)

	return

!---------  error conditions

140	call isat_abort(  'ci_init6', 2, mess = 'error reading nstr' )
142	call isat_abort(  'ci_init6', 3, mess = 'hit end reading nstr' )
144	call isat_abort(  'ci_init6', 4, mess =
     1	   ' error reading stream data. ist=', isv = ist )
146	call isat_abort(  'ci_init6', 5, mess =
     1	   ' hit end reading stream data. ist=', isv = ist )

	return

!------------  internal subroutines  ------------------------------------

	contains

	subroutine op_species !---  print information about species

	integer :: is, ie, ir

	write(luout,*)' '
	write(luout,*)'Species, mol. wt., and their elemental composition'
	write(luout,*)' '
	do is = 1, ns
	   write(luout,600)is,snames(is),amolwt(is),(ev(is,ie),ie=1,ne)
	end do
	call isat_flush(luout)

	if( iorv == 1 ) then
	   write(luout,*)' '
	   write(luout,*)'Reaction vectors (transposed) '
	   write(luout,*)'(set  iorv=0  in  ci.nml  to suppress ',
     1	                 'this output.)'
	   write(luout,*)' '
	   write(luout,610)(is, is=0, ns)
	   do ir = 1, nr
	      write(luout,610)ir,(rv(is,ir),is=1,ns)
	   end do

	   write(luout,*)' '
	   call isat_flush(luout)
600	   format(1x,i4,2x,a16,f8.4,20i3)
610	   format(1x,25i3,/,(4x,24i3))
	endif

	end subroutine op_species 

	subroutine op_streams !---  print information about streams

	integer :: is

	write(luout,*)' '
	write(luout,*)'Output from ci_init6: stream information ',
     1	              '(Chemkin units)'
	write(luout,*)' '
	write(luout,600)(pstr(ist),ist=1,nstr)
	write(luout,610)(tstr(ist),ist=1,nstr)
	write(luout,620)(rhostr(ist),ist=1,nstr)
	write(luout,630)(hstr(ist),ist=1,nstr)
	write(luout,635)(ccstrm(nc,ist),ist=1,nstr)
	write(luout,*)' '
	write(luout,*)'Specific mole numbers'
	write(luout,*)' '
	do is = 1, ns
	   write(luout,640)is,snames(is),(ccstrm(is,ist),ist=1,nstr)
	end do
	write(luout,*)' '
	write(luout,*)'Mass fractions '
	write(luout,*)' '
	do is = 1, ns
	   write(luout,640)is,snames(is),(ystr(is,ist),ist=1,nstr)
	end do

600	format('Press.  ',1p,20e12.4)
610	format('Temp    ',1p,20e12.4)
620	format('Dens.   ',1p,20e12.4)
630	format('Enth    ',1p,20e12.4)
635	format('Sen.Enth',1p,20e12.4)
640	format(i4,2x,a16,2x,     1p,20e12.4)

	end subroutine op_streams

	end subroutine ci_init6

!=========================================================================

	subroutine cicmp6( cc, mode, comp )

!  given compact representation, return composition:

!  input:
!	cc    - compact representation
!	mode  = 1 - express species as mole fractions
!	      = 2 - express species as mass fractions
!	      = 3 - express species as specific mole numbers
!	comp(1) - pressure in Chemkin units
!  output:
!	comp - composition consisting of:
!	comp(1:ns) = species concentrations
!	comp(ns+1) = density (g/cc)
!	comp(ns+2) = temperature (K)
!	comp(ns+3) = pressure (atm)
!	comp(ns+4) = enthalpy (ergs/mol)

	implicit none

	real(k_dp), intent(in)    :: cc(nsp1)
	real(k_dp), intent(inout) :: comp(nsp4)
	integer, intent(in)       :: mode

	real(k_dp) :: hs, h, y(ns), temp, press, dens
	integer    :: check = 0, info

	press      = comp(1)
	hs         = cc(nsp1)

!  evaluate and set thermodynamic variables

	call hs2h( hs, href, cc(1:ns), ns, h )
	call phi2y( cc(1:ns), amolwt, ns, y )
	call temphz( h, cc(1:ns), temp, check, info )
	call ctrhoy( press, temp, y, dens, gas )

	comp(ns+1) = dens
	comp(ns+2) = temp
	comp(ns+3) = press
	comp(ns+4) = h

!  set species and make conversions

	if( mode == 1 ) then
	   comp(1:ns) = cc(1:ns) / sum( cc(1:ns) )
	elseif( mode == 2 ) then
	   comp(1:ns) = y
	elseif( mode == 3 ) then
	   comp(1:ns) = cc(1:ns)
	endif

	return
	end subroutine cicmp6

!=========================================================================

	subroutine cirxn6( t, c0, ct, dpt )

!  chemistry interface routine to return composition c(t) resulting from
!       reaction for a time t from the initial composition c(0).  Also
!       returned are density, pressure and temperature.
!       This version for  modeci = 6  - direct integration.

!  input:
!       t     - time, duration of reaction (real)
!       c0    - initial composition vector (real)
!       dpt(2) - pressure
!  output:
!       ct    - final composition vector (real)
!       dpt   - density, pressure and temperature (real)

	implicit none

	real(k_dp), intent(in)    :: t, c0(nsp1)
	real(k_dp), intent(out)   :: ct(nsp1)
	real(k_dp), intent(inout) :: dpt(3)

	logical, save    :: first_call=.true.
	integer    :: need(3), iusr(1)
	real(k_dp) :: press, temp, x(nxx), rusr(1), f(nsp2),
     1	              dfdx(nsp2,nxx), hvar(1)

!  on first call set  dtc;  check t on subsequent calls

	if( first_call ) then
	   first_call = .false.
	   dtc = t
	else
	   if( const_dt ) then
	      if( t /= dtc ) call isat_abort( 'cirxn6', 1,
     1 	        mess='time step not constant: dt, dtc = ',
     2	        rvar=(/ t, dtc /) )
     	   endif
	endif

	need(1) = 1
	need(2) = 0
	need(3) = 0

!  assemble x

	x(1:nsp1) = c0
	press     = dpt(2)

	if( mode_pdt == 1 ) then
	   x(nsp2) = press
	   x(nsp3) = t
	elseif( mode_pdt == 2 ) then
	   x(nsp2) = press
	elseif( mode_pdt == 3 ) then
	   x(nsp2) = t
	endif

!  perform direct integration to obtain f = {phi, hs, T}

        call cirmap1( need, nxx, x, nsp2, 0, iusr, rusr,
     1	            f, dfdx, hvar )

!  extract ct, dpt

	ct     = f(1:nsp1)
	temp   = f(nsp2)
	dpt(1) = press / ( gascon * temp * sum( f(1:ns) ) )
	dpt(3) = temp

	return
	end subroutine cirxn6

!=========================================================================

	subroutine cirxn7( t, c0, ct, dpt )

!  chemistry interface routine to return composition c(t) resulting from
!       reaction for a time t from the initial composition c(0).  Also
!       returned are density, pressure and temperature.
!       This version for new version of ISATAB 4/28/00

!  input:
!       t     - time, duration of reaction (real)
!       c0    - initial composition vector (phi, hs) (real)
!	dpt(2) - pressure
!  output:
!       ct    - final composition vector (real)
!       dpt   - density, pressure and temperature (real)

	implicit none

	real(k_dp), intent(in)    :: t, c0(nsp1)
	real(k_dp), intent(inout) :: dpt(3)
	real(k_dp), intent(out)   :: ct(nsp1)

	real(k_dp), allocatable, save :: rinfo(:)
	integer, save    :: ifst=0, idtab, info(100), nx, nf, nh, mode

      real(k_dp), save :: rusr(1)
	real(k_dp) :: x(nxx), f(nsp2), press, phi(ns), y(ns),
     1	              dfdx(nsp2,nxx), hvar(1), stats(100),
     2	              vsmall, temp, hs, h, hlo, hhi

	integer    :: iusr(1), ismall
	integer    :: tcheck = 0, tinfo

! set lu_rec >=0 to output diagnostics on record temperatures to lu=lu_rec
	integer            :: lu_rec = -1 !XXX laniu 08-22-07
	real(k_dp), save   :: icall=0.d0, thib=0.d0, thia=0.d0, 
     1                      tlob=1.d30, tloa=1.d30, tb, tin
	
	external cirmap1

!--------------  start of execution  -------------------------------

	press = dpt(2)

	if( ifst /= 0 ) go to 60
	ifst = 1

!===============  initialization on first call =======================

!----------  check press and dt,  set dtc  -----------------------

	if( press <= 0.d0 ) call isat_abort('cirxn7', 1,
     1	   mess='bad pressure', rsv = press )

	if( t <= 0.d0 ) call isat_abort('cirxn7', 2,
     1	   mess='bad dt', rsv = t )

     	if( const_dt ) then
	   dtc = t
	else
	   dtc = -1.d0
	endif

	if( dtref <= 0.d0 ) dtref = t   !  SBP  8/20/04

!----------  check nxx set correctly

	if( mode_pdt == 1 ) then
	   if( nxx /= nsp3 ) call isat_abort('cirxn7', 3,
     1	       mess='bad nxx: nxx, nsp3, mode_pdt=',
     2	       ivar = (/ nxx, nsp3, mode_pdt/) )
	elseif( mode_pdt == 2  .or.   mode_pdt == 3 ) then
	   if( nxx /= nsp2 ) call isat_abort('cirxn7', 4,
     1	       mess='bad nxx: nxx, nsp2, mode_pdt=',
     2	       ivar = (/ nxx, nsp2, mode_pdt/) )
	else
	   if( nxx /= nsp1 ) call isat_abort('cirxn7', 5,
     1	       mess='bad nxx: nxx, nsp1, mode_pdt=',
     2	       ivar = (/ nxx, nsp1, mode_pdt/) )
     	endif

!-----------  set isat controls  ----------------------------------

	idtab  = 1
	mode   = 0
	nx     = nxx
	nf     = nsp2
	nh     = 0

	info     = info_n
	info( 1) = 1            !  iscale
	! VH - info(12) should be set using ciparam
	! info(12) = 2            !  isatop
	info(13) = lu_err       !  error output
	info(66) = ci_info_n(15)          !  mpi_uniq
	if( rmap1t /= 0.d0 ) info(22) = 1 ! return statistics

	allocate( rinfo(50+nxx+nsp2) )
	rinfo = 0.d0
	rinfo(1:50) = rinfo_n

!  scaling for x

	rinfo(51:50+ns) = phiref
	rinfo(50+nsp1)  = hsref

	if( mode_pdt == 1 ) then
	   rinfo(50+nsp2)  = pressref 
	   !SBP 4/26/02:  
	   rinfo(50+nsp3)  = dtref
	elseif( mode_pdt == 2 ) then
	   rinfo(50+nsp2)  = pressref 
	elseif( mode_pdt == 3 ) then
		!SBP 4/26/02:
	   rinfo(50+nsp2)  = dtref
	endif

	allocate( xref(nx) )
	xref = rinfo(51:50+nx)
!  scaling for f

	rinfo(50+nx+1:50+nx+ns) = phiref 

	rinfo(50+nx+nsp1) = hsfref
	rinfo(50+nx+nsp2) = tempref

      rusr(1)       = 1.d0
	isatab_called = .true.

!==========  end of initialization  ============================

60	continue

!-----------  check constant press and dt  ----------------------

	if( const_pr ) then
	   if( press /= prc ) call isat_abort('cirxn7', 6,
     1	      mess='pressure not constant', rvar=(/ press, prc/) )
     	endif

	if( const_dt ) then
	   if( t /= dtc ) call isat_abort('cirxn7', 7,
     1	      mess='dt not constant', rvar=(/ t, dtc/) )
     	endif

!-----------  set parameters for isat ---------------------------

!  x = { phi(1:ns), hs, p, t }  (at time 0),  mode_pdt = 1
!  x = { phi(1:ns), hs, p }     (at time 0),  mode_pdt = 2
!  x = { phi(1:ns), hs, t }     (at time 0),  mode_pdt = 3
!  x = { phi(1:ns), hs }        (at time 0),  mode_pdt = 4

!  f = { phi(1:ns), hs, T }     (at time t)

	x(1:nsp1) = c0

	if( mode_pdt == 1 ) then
	   x(nsp2) = press
	   x(nsp3) = t
	elseif( mode_pdt == 2 ) then
	   x(nsp2) = press
	elseif( mode_pdt == 3 ) then
	   x(nsp2) = t
	endif
	
!-----------  enforce realizability

	call cireal( 2, treal, x(1:ns), ismall, vsmall )

!--- call isatab  =================

! crashing in isatab
	call isatab( idtab, mode, nx, x, nf, nh, 1, cirmap1,
     1	   iusr, rusr, info, rinfo, f, dfdx, hvar, stats )

!--- test for non-retrieve
	if( nint( stats(9) ) == 8 ) then
	   dpt(1) = -1.d0  !  signal that f has not been returned
	   return
	endif

!--- test  cirmap1  if so required
	if( rmap1t /= 0 ) then
	   if( nint( stats(9) ) == 5  .and.  rmap1t /= 0 )
     1	     call rmap1_test( idtab, nx, x, nf, nh, 1, cirmap1,
     2	     iusr, rusr, info, rinfo, f, dfdx, hvar, stats )
     	endif
     	
!  enforce realizability

      phi  = f(1:ns)  !  store uncorrected species....
      temp = f(nsp2)  !  ...and temperature
      tb   = temp

	call cireal( 2, treal, f(1:ns), ismall, vsmall )
	
	if( ismall > 0  .and.  kreal_h == 1 ) then  ! adjust hs and T to conserve enthalpy
		 hs = f(nsp1)
	   call hs2h( hs, href, phi, ns, h )  !  enthalpy prior to correction
	   call phi2y( f(1:ns), amolwt, ns, y )  ! corrected mass fractions	   
	   
	   if( kreal_t > 0 ) then  ! check that T is in range
	      call cthbms( tbadlo*1.01, y, hlo, gas )
	      call cthbms( tbadhi*0.99, y, hhi, gas )
	   
	      if( h < hlo  .or.  h > hhi ) then            !laniu
	         if( kreal_t == 2 ) call isat_abort('cirxn7', 8,
     1	          mess='temperature out of range after realiz. corr.',
     2            rvar=(/ h, hlo, hhi/) )
     
	         if( h < hlo ) then  !  make adjustment
	            h    = hlo
	            temp =1.01* tbadlo
	         else
	            h    = hhi
	            temp = 0.99*tbadhi
	         endif 
	      endif
	   endif 
	      
	   call temphz( h, f(1:ns), temp, tcheck, tinfo )  !  adjusted T and hs
	   call h2hs( h, href, f(1:ns), ns, hs )
	   f(nsp1) = hs
	endif
	
	ct = f(1:nsp1)

!  set density

	dpt(1) = press / ( gascon * temp * sum( f(1:ns) ) )
	dpt(3) = temp
	
	if( lu_rec < 0 ) return  !  all done if diagnostics not needed
	
!  treat record temperatures
      icall = icall + 1.d0
      if( tb   >= tlob  .and.  tb   <= thib  .and.  
     1    temp >= tloa  .and.  temp <= thia )  return  ! no record
     
	tlob = min( tlob, tb )
	thib = max( thib, tb )
	tloa = min( tloa, temp )
	thia = max( thia, temp )

	call hs2h( x(nsp1), href, x(1:ns), ns, h ) 
	call phi2y( x(1:ns), amolwt, ns, y )  
	call temphz( h, x(1:ns), tin, tcheck, tinfo )  
      
         write(lu_rec,'(1p,e20.12,1p,9e14.5)') icall, tlob, tloa, 
     1                              thib, thia, tin, tb, temp

	return
	end subroutine cirxn7

!=========================================================================

	end module ci_6
