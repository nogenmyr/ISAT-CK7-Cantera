!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

	program pasr

!  PASR - partially-stirred reactor, using pairwise-mixing model.

!--------------------------------------------------

!  S.B. Pope  10/7/2006

!  BUG FIX: 04/06/2012
!  Prior to this date, the residence time was exactly twice that specified

!  PaSR test cases set up for testing ISAT-CK 5.0
!  16-species skeletal mechanism for methane.

!  Input files:
!	streams.in
!	pasr.nml  (optional)
!	pasr.inp  (optional)
!
!  Output:
!	ci.op		     - check for correctness of test
!	pasr_op.txt    - (generated if screen=.false. in pasr.nml)
!	pasrm.op	     - postprocess with pasr_m.m
!	isat_1.op	     - postprocess with isat5.m
!     isat_err_1.cdf - postprocess with cdf_err.m
!     isat_1.leaf    - postprocess with leaf.m
!
!--------------------------------------------------

!  The reactor consists of  n  particles.  Particles flow in and out
!  at a rate determined by the residence time  tres.
!  There are  nstr  inflowing streams, with relative mass flow
!  rates  flstrm.  The initial condition for all particles
!  is set to the composition of the first stream.   The particles
!  react; and they mix on a timescale  tmix  with partners which are
!  randomly re-assigned on a timescale  tpair.
!  The pressure varies sinusoidally in time between prmin and prmax
!  with a period of prtime.

!  physical parameters:
!	tres	- residence time (s)
!	tmix	- mixing time scale (s)
!	tpair	- pairing time scale (s)
!	anres	- length of run (number of residence times)
!	flstrm  - relative flow rates of streams
!	prmin	- minimum pressure (atm.)
!	prmax	- maximum pressure (atm.)
!	prtime	- period of pressure variation (s)

!  The largest time step  dtmax  is determined by the specified constant  cdtrp.
!  The smallest time step  dtmin  is smaller than  dtmax  by the factor
!  dtfac (i.e., dtmin = dtmax * dtfac).  The time step used is random,
!  uniformly distributed between dtmin and dtmax.
!  The largest time step  dtsmax  for mixing and reaction substeps
!  is determined by the specified constant  cdtmix.

!  numerical parameters:
!	n	    - number of particles (must be even)
!	cdtrp	- time-step control (must be less than 0.25)
!	dtsfac	- min-to-max substep ratio
!	cdtmix	- sub-step control, based on tmix (less than 0.5)

!  input files:
!	streams.in - required by ISAT-CK (provides value of nstr)
!	pasr.nml   - namelist file which can be used to change the values
!		     of the physical and numerical parameters

!  output files
!	pasrs.op - properties of a single sample particle
!	pasrm.op - ensemble means

!  control
!	cpu_lim	 - limit on CPU time (in seconds)
!	call_cisave - if .true., subroutine cisave is called at the end to
!                     checkpoint the ISAT table.

!------------------------------------------------------------------------

	use isat_rnu
	use ci_dat,  only: modeci
	use ci_dat8, only: nrs, autou, href
	use ci_cksubs, only: hs2h
	use ci_stats

	implicit none

	real, allocatable :: f(:,:), fin(:,:), flstrm(:), fmean(:),
     1	    work(:), temp(:), dens(:), dfin(:), cc(:), dfm(:), ct(:)

     	real :: tres, tmix, tpair, prmax, prmin, prtime, cdtrp, cdtmix,
	1        dtfac, anres, etola, dxsmax, stomby, cpu_lim, patm, pi,
     1        dtmax, dtmin, dtave, dpt(3), flstrmd(3), cpu_secs,
     1        cpu_start, dtsub, tend, anstep, anpsub, tempmaxx, prout,
     1        prpar, presso, dt, t, tp, ftemp,  prave, pramp,
     1        pnow, press, dpress, tempmean, tempmin, tempmax,
     1        queries, out_inc, out_next, rinfo_val1, rinfo_val2,
     1        thresh, ecc_freq,
     1        T_min = huge(0.)
      character(10), allocatable :: cname(:)
	logical :: exist, test_acc = .false., call_cisave = .false.,
     1           screen=.true. , defer=.false., leaf_op=.false.,
     2           part_in=.false., part_out=.false., test_grad=.false.,
     3           test_jac=.false., test_lmap=.false., find_errs=.false.,
     4           sp_order=.false., test_subs=.false., full_op=.true.,
     5           test_map=.false., test_skel=.false., gen_map=.false.

	integer :: n, iskip, iblock, ichin, ichout, isatop, idites,
     1           ci_info(20), info(100), nrec, isat_vers, lu_op, icheck,
     1           irnsed, lu, ncomp, nfull, nstr, ns, lus, luf,
     1           lum, lut, lur, nsub, ist, i, istep, nprout, npick, ii,
     1           istrm, inflow, nprpar, ifst, k, isub, nfulio, pre_set,
     1           info_i, rinfo_i1, rinfo_i2, info_val, lu_dist, nx,
     1           kr_out, lug, imap = 1

	real(kind(1.e0)) :: ci_rinfo(20), rinfo(50), stats(100)
	real(kind(1.e0)) :: cpu_1, cpu_2 !XXX

	real(kind(1.d0)), allocatable :: dcc(:), dct(:)
	real(kind(1.d0)) :: ddt, ddpt(3), h, hs, tmp
	integer :: itcheck = 0, itinfo

! random seed - make compiler independent
	integer, allocatable :: iseed(:)
	integer :: isize

!  dens. temp update variables
	real(kind(1.d0)) :: etol, stats_dpt(10)
	logical :: success

!  input parameters can be changed through namelist
	namelist / pasrnml / n, tres, tmix, tpair, prmax, prmin,
     1      prtime, cdtrp, cdtmix, dtfac, anres, iskip, iblock,
     2	    cpu_lim, call_cisave, isat_vers, etola,
     3      ichin, ichout, isatop, idites, icheck, stomby, sp_order,
     4      screen, test_grad, test_jac, test_lmap, find_errs,
     5      test_subs, full_op, test_map, test_skel, gen_map, kr_out

! initlal values
	success = .false.

!  default parameters for CH4-air combustion
	data n/ 100 /
	data isat_vers/ 50 /  ! version of ISATAB
	data tres, tmix, tpair/ 10.e-3, 1.0e-3, 1.0e-3 /
	data flstrmd/ 0.05, 0.85, 0.1 /
    	data prmax, prmin, prtime/ 1.0, 1.0, 1.0e-3 /
	data cdtrp, cdtmix, dtfac/ 0.1, 0.040, 1.0 /
	data anres/ 3 /

	data iskip, iblock/ 1, 100 /
	data out_inc/ 1.02 /	!  frequency of ISAT performance output

!  ISATAB parameters
	data ichin, ichout/ 0, 0 /
	data isatop, idites, icheck, pre_set/ 1, 0, 0, 0 /
	data etola, stomby, ecc_freq/ 1.e-4, 50., -2.d0 /
	data cpu_lim/ 1.e30 /

	data patm/ 1.01325e6 /
	data irnsed/ 1 /
	data kr_out/ 1 /


!----------------  initialization  ---------------------------------------
!--- read namelist (if the file pasr.nml exists)

	inquire( file = 'pasr.nml', exist = exist )
	if( exist ) then
	   call isat_lu( lu )
	   open( lu, file = 'pasr.nml', position = 'rewind',
     1              action = 'read' )
	   read( lu, nml = pasrnml )
	   close( lu )
	endif

	if( screen ) then
	   lu_op = 0
	else
	   call isat_lu( lu_op )
         open( lu_op, file='pasr_op.txt')
	endif

!--- initialize chemistry interface
	ci_info  = 0
	ci_rinfo = 0.d0

        ci_info(1)  = 1   ! const pr
	if( prmax /= prmin ) ci_info(1) = 0

	ci_info(2)  = 1   ! const dt
	if( dtfac /= 1.d0 ) ci_info(2) = 0

	ci_info(7)  = 1   ! ickcorr - correct thermodata?
	ci_info(8)  = 2   ! m_sens
	ci_info(11) = 2   ! exp_m
	ci_info(14) = 1   ! op_params
	ci_info(16) = 0   ! lu_err
	ci_info(17) = 0   ! kreal_h
	ci_info(18) = 0   ! kreal_t
	ci_info(19) = 1	  ! m_jac
	ci_info(20) = 1	  ! m_ckwyp

!--- initialize ISATAB parameters
	info     = 0
	rinfo    = 0.d0
	info(10) = ichin
	info(11) = ichout
	info(12) = isatop
	info(15) = pre_set
	info(20) = 1  !  stats0
	info(28) = idites
	info(60) = icheck
	info(61) = 0  !  isat_vers

	rinfo(1)  = etola
	rinfo(8)  = stomby
	rinfo(21) = ecc_freq

	if( .false. ) then    !XXX emulate ISAT4
	   info(37)  =  3    ! mode_eoi
	   info(39)  = -1    ! degen_r
	   info(41)  = -1    ! degen_s
	   rinfo(5)  = 1d2   ! eoi_lim
	   rinfo(16) = -1.d0 ! ad_theta
	endif

      if( defer ) info(29) = 1
	if( idites > 0  .or.  ecc_freq < 1.d0 ) then
	   info(70) = 1  !  cdf_err
	   test_acc = .true.
	endif

!--- read from pasr.inp
!  This file can be used (for parametric testing) to change one value
!  of info, and/or up to two values of rinfo.
!  (A value set here will be over-ridden by a value set in isat_1.nml.)

	inquire( file = 'pasr.inp', exist = exist )
	if( exist ) then
	   call isat_lu( lu )
	   open( lu, file = 'pasr.inp', position = 'rewind',
     1         action = 'read' )
	   read( lu, * ) info_i, info_val, rinfo_i1, rinfo_val1,
     1                                   rinfo_i2, rinfo_val2
	   close( lu )

	   if( info_i > 0 ) then
	      info( info_i ) = info_val
	      write(lu_op,'(a,i2,a,i6)') 'Value of  info(', info_i,
     1                                 ') set to ', info_val
	   endif

	   if( rinfo_i1 > 0 ) then
	      rinfo( rinfo_i1 ) = rinfo_val1
	      write(lu_op,'(a,i2,a,1p,e13.4)') 'Value of rinfo(',
     1                           rinfo_i1,') set to ', rinfo_val1
	   endif

	   if( rinfo_i2 > 0 ) then
	      rinfo( rinfo_i2 ) = rinfo_val2
	      write(lu_op,'(a,i2,a,1p,e13.4)') 'Value of rinfo(',
     1                           rinfo_i2,') set to ', rinfo_val2
	   endif
	endif

!--- perform initialization

	call sciparam( ci_info, ci_rinfo, info, rinfo )

	call ciinit( ncomp, nfull, nstr )

! compute ns using nfull
	if(modeci > 5) then
	   ! nfull = ns + 4 (dens, p, T, h)
	   ns = nfull - 4
	else
	   ns = nfull
	endif

!--- allocate arrays that do not depend on n

	allocate( fin(ncomp,nstr) )
	allocate( flstrm(nstr) )
	allocate( fmean(ncomp) )
	allocate( dfin(ncomp) )
	allocate( cc(ncomp) )
	allocate( dfm(ncomp+4) )
	allocate( ct(ncomp) )
	allocate( cname(ncomp+4) )

	allocate( dcc(ncomp) )
	allocate( dct(ncomp) )

!--- set default flow rates

	flstrm = 1./float(nstr)
	if( nstr >= 3 ) flstrm(1:3) =  flstrmd
	if( nstr >= 4 ) flstrm(4:nstr) = 0.

	flstrm = flstrm / sum(flstrm)

!--- allocate arrays that depend on n

	allocate( f(n,ncomp) )
	allocate( temp(n) )
	allocate( dens(n) )
	allocate( work(n) )

!--- open files, etc.

	pi     = 4. * atan( 1. )
        call rnused( irnsed )

!--- Set random seed in a compiler independent way
!--- 03/09/2012, changed by Varun
	call random_seed(size = isize)
	allocate(iseed(isize))

	iseed = 12345 * (/ (i - 1, i = 1, isize) /)
	call random_seed( put=iseed )
!---end random seed

        call isat_lu( lus )
        if(full_op) call isat_lu( luf )
        call isat_lu( lum )
        call isat_lu( lut )
        call isat_lu( lur )
	call isat_lu( lug )
        call isat_lu( lu_dist )
        open ( lus, file ='pasrs.op' )
        if(full_op) open ( luf, file ='pasrfull.op' )
	if(gen_map) open ( lug, file ='rxn_map.op' )
        open ( lum, file ='pasrm.op' )
	open ( lut, file ='cputime.op', position='append' )
	open ( lur, file ='runtime.op' )
        open ( lu_dist, file ='use_dist.txt' )

!--- derived parameters

        dtmax  = min( cdtrp, 0.25) * min( tres, tpair )
	dtmin  = min( dtfac, 1.0 ) * dtmax
	dtave  = 0.5 * ( dtmax + dtmin )

        dtsub  = min( cdtmix, 0.5) * tmix
        nsub   = 1 + ifix( dtmax / dtsub )

        tend   = anres * tres
        anstep = tend / dtave
        anpsub = anstep * n * nsub

!--- write values of parameters

        write(lu_op,6 ) ns
        write(lu_op,8 ) nstr
        write(lu_op,9 ) flstrm
        write(lu_op,10) n
        write(lu_op,12) tres
	write(lu_op, *) '04/06/2012: Prior to this date, due to a bug',
     1	                ' tres was twice the specified tres'
        write(lu_op,14) tmix
        write(lu_op,16) tpair
        write(lu_op,17) dtmin, dtmax
        write(lu_op,18) anstep
        write(lu_op,20) nsub
        write(lu_op,22) anpsub
        write(lu_op,24) prmax
        write(lu_op,26) prmin
        write(lu_op,28) prtime

6       format(' Number of species                 =', i8 )
8       format(' Number of streams                 =', i8 )
9       format(' Flow rates of streams             =', ((5f8.4)/) )
10      format(' PASR starting: number of particles=', i8 )
12      format(' Residence time, tres (seconds)    =', 1p,e13.4)
14      format(' Mixing time, tmix (seconds)       =', 1p,e13.4)
16      format(' Pairing time, tpair (seconds)     =', 1p,e13.4)
17      format(' Min. and max. time step (seconds) =', 1p,2e13.4)
18      format(' Number of steps                   =', 1p,e13.4)
20      format(' Number of sub-steps per step      =', i8 )
22      format(' Number of particle-sub-steps      =', 1p,e13.4)
24      format(' Maximum pressure (atm.)           =', 1p,e13.4)
26      format(' Minumum pressure (atm.)           =', 1p,e13.4)
28      format(' Period of pressure variation (sec)=', 1p,e13.4)

!-----  get stream compositions

      do ist = 1, nstr
         call scistrm( ist, ncomp, dfin, dpt )

	   fin(:,ist) = dfin(:)
	end do

!-----  set initial condition

      do i = 1, n
	   f(i,:)   = fin(:,1)
	end do

	if( part_in ) then  !  read particle properties
	   call isat_lu( lu )
	   open( lu, file='particle_f.dat', form='unformatted')
	   read(lu) f(1:n,1:ncomp)
	   close(lu)
	   write(lu_op,*)'Particle data read from file'
	endif

	tempmaxx = 0.
      prout    = 0.
      prpar    = 0.
	presso   = prmin * patm

!----  final initializations

      write(0,*)'PaSR starting'
      write(0,*)'...'
      istep  = 0
	dt     = dtave
      t      = -dt
	tp     = -dt / ( nsub * prtime )
	out_next = 1.
	call cpu_time( cpu_start )

!========= optionally, perform some tests and stop =========================
	if( test_grad .or. test_jac .or. test_lmap .or. find_errs
     1  .or. test_subs .or. sp_order .or. test_map .or. test_skel ) then

! Procedure:
! 1/ run pasr with test_grad=.false. and n_add_op>1 in isatab, to generate samples of x in adds.op
! 2/ run pasr with test_grad=.true. and n_add_op<0 in isatab to measure errors
! 3/ optionally, reset kr_out (e.g., based on observed large error) and re-run (as in 2)
! 4/ post-process with ci_test_grad.m
! Note: test_grad and kr_out can be set in pasr.nml; n_op_adds can be set in isat_1.nml.

	   dtsub = dt/nsub
	   press = prmax*patm
	   cc = fin(:,1)
	   dpt(2) = press
	   call scirxn( dtsub, ncomp, cc, ct, dpt )  !  make call to complete initialization

	   nx = ncomp
	   if( dtfac /= 1. )    nx = nx + 1  ! variable time step
	   if( prmin /= prmax ) nx = nx + 1  ! variable pressure

	   if(test_grad) then
              ! gradient test
	      call ci_test_grad( nx, ncomp+1, kr_out )
	   else if(test_lmap) then
	      ! linrmap test
	      call ci_test_linrmap( nx, ncomp+1, kr_out )
	   else if(test_jac) then
	      ! jacobian test
              call ci_test_jacobian( nx, ncomp+1, kr_out )
	   else if(test_map) then
	      ! rxn mapping errors
	      call ci_test_rxn_map( press, dtsub )
	   else if(test_skel) then
	      ! rxn mapping errors
	      call ci_test_skeletal( 110, press, dtsub )
	   else if(find_errs) then
	      ! compute errors
	      call ci_all_errors( nx, ncomp+1 )
	      !call ci_test_dpt
	      !call ci_test_rzmap
	   else if(sp_order) then
	      ! compute species ordering
	      call ci_cem_species_rm
	      call ci_cem_species_rm_gali
	   else if(test_subs) then
	      dcc = cc
	      ddpt = dpt
	      ddt = dtsub
	      call ci_test_interface( ddt, nx, dcc, dct, ddpt )
	   endif

	   call cpu_time( cpu_secs )
	   cpu_secs = cpu_secs - cpu_start

!       append to cputime.op: mode, nrs and cpu time secs, mins, hours
	   if(modeci == 6 .or. modeci == 7) then
	      write(lut,'(2i8,1p,3e13.4)') modeci, ns, cpu_secs,
     1        cpu_secs/60.0, cpu_secs/3600.0
	   else if(modeci == 8 .or. modeci == 9) then
	      write(lut,'(2i8,1p,3e13.4)') modeci, nrs, cpu_secs,
     1	      cpu_secs/60.0, cpu_secs/3600.0
	   endif

!       Print ci_stats
	   call print_ci_stats(luout=0)

	   stop
	endif
!====================================================================================

!===============  start of time steps  =================================

150   istep  = mod( istep + 1, 1000000000 )

!--- set the time step

      if( istep == 1 ) then
         dt = dtmax ! ...so that dtref=dtmax
	elseif( mod(istep,2) == 1 ) then
	   dt = 2.*dtave - dt
	else
	   if( dtmin == dtmax ) then
	      dt = dtmin
	   else
	      call rnu( dt )
	      dt     = dtmin + dt * ( dtmax - dtmin )
	   endif
	endif

	t = t + dt

!-------  outflow/ inflow -------------------------------------------

! BUG FIX: 04/06/2012
! Prior to this date, the residence time was exactly twice that specified
!       prout  = prout + 0.5 * n * dt / tres
! This is now changed to:
!       prout  = prout + n * dt / tres

	prout  = prout + n * dt / tres
        nprout = prout
        prout  = prout - nprout

        if( nprout .eq. 0 ) go to 250

!  select 2 * nprout pairs of particles at random; put at end of f array

        npick  = 2 * nprout
        call pickpr( npick, n, ncomp, n, f )

!  set alternate particles to inflowing properties

        i         = n + 2
        do ii = 1, nprout
           i         = i - 2
           istrm     = inflow( nstr, flstrm )
	   f(i,:)    = fin(:,istrm)
	end do

!-------  pairing  ------------------------------------------------

250     prpar  = prpar + 0.5 * n * dt / tpair
        nprpar = prpar

        if( nprpar <= 1 ) go to 290

        prpar  = prpar - nprpar

!  select nprpar pairs of particles at random; put at end of f array

        call pickpr( nprpar, n, ncomp, n, f )

!  rotate particle properties

        ifst = n - 2 * nprpar + 2
        do k = 1, ncomp
           ftemp    = f(n,k)

           do i = n, ifst+2, -2
	      f(i,k)   = f(i-2,k)
	   end do

	   f(ifst,k)= ftemp
	end do

!--------  loop over mix-react-mix substeps

290	dtsub    = dt / nsub

        do isub = 1, nsub

!-----   sub-step of mixing

           call mix( dtsub, tmix, n, ncomp, n, f )

!--- set the pressure

	   tp     = tp + dtsub / prtime
	   tp     = mod( tp, 1. )
	   prave  = 0.5 * ( prmax + prmin )
	   pramp  = 0.5 * ( prmax - prmin )
	   pnow   = prave - pramp * cos( 2. * pi * tp )
	   press  = pnow * patm
	   dpress = press - presso
	   presso = press

!-----  sub-step of reaction
           do i = 1, n
	      cc(:)    = f(i,:)
	      dpt(2)   = press

	      call scirxn( dtsub, ncomp, cc, ct, dpt ) !  call to ISAT-CK

              ! write rxn mappings to a file
	      if( gen_map .and. modeci == 6 ) then
		 if( mod(istep,3) == 0 .and. isub == 1 ) then
		    dcc = cc
		    hs = cc(ns+1)
		    call hs2h( hs, href, dcc(1:ns), ns, h )
		    call temphz( h, dcc(1:ns), tmp, itcheck, itinfo )
		    write(lug, '(i8,1p1000e20.11)') imap, cc, ct, h, tmp, dpt
		    imap = imap + 1
		 endif
	      endif

	      dens(i)  = dpt(1)
	      temp(i)  = dpt(3)

	      if( temp(i) < T_min ) then !XXXXXXXXXX
		 T_min = temp(i)
!       write(0,'(a,3i5,1pe13.4)') 'PaSR, low T ', istep,
!     1                                     isub, i, T_min
	      endif

!----- increment sensible enthalpy

	      ct(ncomp) = ct(ncomp) + dpress / dpt(1)
	      f(i,:)   = ct(:)
	     end do


	end do

!========================  end of particle evolution ==============

	if( defer ) then  !  check for deferred maintenance and perform
         call scisat( 6, info, rinfo, stats )  ! obtain stats
	   if( stats(21) > 0.e0 ) call scisat(17, info, rinfo, stats)
	   if( stats(22) > 0.e0 ) call scisat(9 , info, rinfo, stats)
	   if( stats(23) > 0.e0 ) call scisat(12, info, rinfo, stats)
	endif

!----- other output
      if( mod(istep,iskip) .ne. 0 ) go to 500

      if( mod(istep,iblock) .eq. 0  .and.  iblock < huge(1)/2 ) then
           iblock = 2 * iblock
           iskip  = 2 * iskip
      endif

!---  properties of particle 1
	write(lus,405) t, dens(1), press/patm, temp(1), f(1,ncomp),
     1               (f(1,k),k=1,ncomp-1)

!---- output compositions of all the particles
	if(full_op) then
	   do i=1,n
	      write(luf,405) t, dens(i), press/patm, temp(i), f(i,ncomp),
     1	             (f(i,k),k=1,ncomp-1)
	   enddo
	endif

405	format(1p,200e16.7)
	call isat_flush(lus)

!---  evaluate and output mean composition

      call means( n, ncomp, n, f, fmean )

	tempmean = sum( temp ) / float(n)
	tempmin  = temp(1)
	tempmax  = temp(1)
	do i     = 2, n
	   tempmin = min( tempmin, temp(i) )
	   tempmax = max( tempmax, temp(i) )
	end do
	tempmaxx = max( tempmax, tempmaxx )

      nfulio = min(13, ncomp)
      write(lum,600) istep, nprout, nprpar, t,
     1	             (fmean(k), k =1, nfulio), f(1,1), f(2,1),
     2	             tempmean, tempmin, tempmax, tempmaxx
600   format(i9,2i5,1p,30e13.4)
      call isat_flush(lum)

!----  end of output

500   call cpu_time( cpu_secs )
      cpu_secs = cpu_secs - cpu_start

!--- occasional output of ISAT performance (set out_inc < 0 to suppress)
      if( istep > out_next  .and.  out_inc > 0. ) then
         out_next = ceiling( istep * out_inc )
	 queries  = float(n)*float(nsub)*float(istep)
         write(lu_op,'(a,1p,10e10.2)')'q, mu_s, sec, min, hour ',
     1         queries, 1e6*cpu_secs/queries, cpu_secs, cpu_secs/60. ,
     2         cpu_secs/3600.
! print run time stats to runtime.op
	 write(lur,'(1p,2e25.14)') queries, cpu_secs
	endif

	if(mod(istep,50) == 0 .and. autou == 0 .and. .not.success) then
	   etol = 0.01
	   call ci_dens_temp_update(etol, success, stats_dpt)
	   if( success .and. stats_dpt(1) > etol ) then
	      write(0,*) ">>> density and temperature values reset."
	      write(0,*) "old error = ",stats_dpt(1)
	      write(0,*) "new error = ",stats_dpt(2)
	   endif
	endif

!--- test isat_change

	if( .false.  .and.  istep == 1 ) then  !XXXXX
	! rebuild EBTs
	   info  = -12345
	   rinfo = -12345.
	   info(52) = 4 ! pair_cover
	   info(60) = 0 ! icheck
	   call sciparam( ci_info, ci_rinfo, info, rinfo )

	   call cpu_time( cpu_1 )
	   info(81) = 1  ! mode for rebuilding
	   call  scisat( 19, info, rinfo, stats )

	   call cpu_time( cpu_2 )
	   cpu_2 = cpu_2 - cpu_1
	   call scisat( 6, info, rinfo, stats )
	   write(0,'(a,1pe13.4,a,1pe13.4)')'EBTs rebuilt: CPU secs = ',
     1             cpu_2, '   mu sec/leaf = ', 1.e6*cpu_2/stats(12)

       !  redefine affine spece
     	   info  = -12345
	   rinfo = -12345.
	   info(48) = 3            ! aff_na
	   rinfo(15) = 1.d0+1.d-9  ! aff_ratio
	   call sciparam( ci_info, ci_rinfo, info, rinfo )

	   call cpu_time( cpu_1 )
	   call  scisat( 17, info, rinfo, stats )
	   call cpu_time( cpu_2 )
	   cpu_2 = cpu_2 - cpu_1
	   write(0,'(a,1pe13.4,a,1pe13.4)')
     1            'Affine redefined: CPU secs = ', cpu_2

	endif

	if( istep== -1 ) then  !  re-set CDF  XXXX
	   info(80) = 2
	   info(81) = 10000
	   info(84) = 100
	   rinfo(41) = 1.d-9
	   rinfo(44) = 1.d-2
	   call  scisat( 25, info, rinfo, stats )
	endif

      if( istep==2**nint(log(float(istep))/log(2.)) ) then  !  test isatab_mfu_dist  XXXX
         rinfo(40) = 1.25
         thresh    = 1.5
         call scisat( 28, info, rinfo, stats )
         do i = 1, 100
            write(lu_dist,'(2i10,1p,2e20.10)') istep, i, stats(i),thresh
            if( thresh > 1.d8 ) then
               thresh = rinfo(40) * thresh
            else
               thresh = max( thresh+1.d0, 0.5d0+nint(rinfo(40)*thresh) )
            endif
         end do
         call isat_flush(lu_dist)

      endif

	if( t < tend  .and.  cpu_secs < cpu_lim ) go to 150

!============== end of time steps  ===============================================

      if( part_out ) then  !  checkpoint particle properties
	   call isat_lu( lu )
	   open( lu, file='particle_f.dat', form='unformatted')
	   write(lu) f(1:n,1:ncomp)
	   close(lu)
	   write(lu_op,*)'Particle data written to file'
	endif

	if( test_acc ) then  !  output CDF of error
	   info(80) = 3
	   call  scisat( 25, info, rinfo, stats )
	endif

	nrec = -2
      if( call_cisave ) call cisave(nrec)

      write(0,*)'pasr: stopped, t, cpu_secs, nrec =',
     1	           t, cpu_secs, nrec

      write(lu_op,*)'pasr: stopped, t, cpu_secs, nrec =',
     1	           t, cpu_secs, nrec

! get MFU information
	rinfo(41) = 0.9
	call scisat( 10, info, rinfo, stats )
	write(lu_op,*)' '
	write(lu_op,'(a,1p,e13.4,i9)')'MFU: fraction, index = ',
     1          rinfo(41), nint(stats(95))

! get retrieve error information
	if( test_acc ) then
	   write(lu_op,*)' '
	   write(lu_op,*)'Retrieve errors'
	   write(lu_op,*)'Fraction    error'

         do i = 1, 5
	      rinfo(41) = 1. - 0.1**i
	      call scisat( 11, info, rinfo, stats )
	      write(lu_op,'(1p,2e13.4)') rinfo(41), stats(96)
	   end do
	endif

!  print x and f ranges
	info(81) = 1
	info(82) = 0
	call scisat( 23, info, rinfo, stats )

!  output leaf properties
	if( leaf_op ) then
	   info(81:83) = 0
	   info(83)    = 1
	   call scisat( 13, info, rinfo, stats)
	endif

!  get stats and complete writing of output
	info(81) = 2
	call scisat( 6, info, rinfo, stats )

	write(lu_op,*)' '
	write(lu_op,'(a,1p,e13.4)')'queries  = ', stats(1)
	write(lu_op,'(a,1p,e13.4)')'prim_ret = ', stats(2)
	write(lu_op,'(a,1p,e13.4)')'sec_ret  = ', stats(3)
	write(lu_op,'(a,1p,e13.4)')'grow     = ', stats(4)
	write(lu_op,'(a,1p,e13.4)')'add      = ', stats(5)
	write(lu_op,'(a,1p,e13.4)')'replace  = ', stats(6)
	write(lu_op,'(a,1p,e13.4)')'dir eval = ', stats(7)
	write(lu_op,*)
	write(lu_op,'(a,1p,e13.4)')'CPU (mu s)/query = ',
     1                   1.e6*cpu_secs/stats(1)

!   append to cputime.op: mode, nrs and cpu time (mu s/query), total time
	if(modeci == 6 .or. modeci == 7) then
	   write(lut,'(2i8,1p,10e13.4)') modeci, ns,
     1       1.e6*cpu_secs/stats(1), cpu_secs,
     1       cpu_secs/60.0, cpu_secs/3600.0
	else if(modeci == 8 .or. modeci == 9) then
	   write(lut,'(2i8,1p,10e13.4)') modeci, nrs,
     1       1.e6*cpu_secs/stats(1), cpu_secs,
     1       cpu_secs/60.0, cpu_secs/3600.0
	endif

! Print ci_stats
	call print_ci_stats(luout=0)

	stop
	end program pasr
