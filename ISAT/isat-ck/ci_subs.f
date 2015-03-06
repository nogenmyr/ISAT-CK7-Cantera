!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

!  This file contains all of the Chemistry Interface routines, namely:

!=====================================================================
! The following comments are to precede the listing formed by extract

!BEGEXTRACT
! The calling sequences for the double-precision Chemistry Interface
! subroutines are given in the following order.  The names of the
! single-precision versions are shown to the right.

!   double precision                             single precision

!   ciinit( ncv, nfullv, nstrms )                          ciinit

!   cistrm( istr, ncv, c, dpt )                            scistrm

!   cicomp( ncv, c, krep, nfullv, comp, cname )            scicomp

!   cirxn( t, ncv, c0, ct, dpt )                           scirxn

!   ciconv( ncv, jd,   jp,   je,   js,                     sciconv
!                din,  pin,  ein,  spin,
!                kd,   kp,   ke,   ks,
!                dout, pout, eout, spout )

!   cisave( nrec )                                         cisave

!   dtchem( t, kspec, nspec, spec0, press, kht, ht, 
!          modecp, modeci, spect, temp, dens )             sdtchem

!   ciparam( ci_info, ci_rinfo, info, rinfo )              sciparam

!   cisat( mode, info, rinfo, stats )                      scisat

! VH 02/23/2010: The following new routines have been added in
! isat-ck-ext/ci_subs_ext.f

!   cistrm_ext(cstrm, ncv, c, dpt)                         scistrm_ext

!   cisize( nsv, nrsv, nev )                               cisize

!   cicomp_sr( ncv, c, krep, nfullv, comp, cname )         scicomp_sr

!   ci_mix_frac( mode, es_names, nc, z, workspace, mf )    sci_mix_frac

!ENDEXTRACT
!=====================================================================
!BEGEXTRACT

        subroutine ciinit( ncv, nfullv, nstrms )

!  chemistry interface initialization routine.

!  input:      - none
!  output:
!      ncv     - number of composition variables (integer)
!      nfullv  - number of items in the full representation (integer)
!      nstrms  - number of streams (integer)

!  files:
!      streams.in - input file of stream information.
!              the first record of streams.in contains the value of the
!              integer variable  modeci  which determines the mode of
!              thermochemistry to be used.  The specification of the
!              subsequent records of  streams.in  depends on  modeci.
!        modeci = 1 - inert, constant-density flow
!        modeci = 2 - mixture fraction formulation (f)
!        modeci = 3 - reaction progress variable formulation (c)
!        modeci =[4]- (f,c) formulation
!        modeci =[5]- general user-supplied thermochemistry
!        modeci = 6 - Cantera/direct integration
!        modeci = 7 - Cantera/ISAT
!        modeci = 8 - Cantera/RCCE/ISAT
!        modeci = 9 - Cantera/RCCE/ICE-PIC
!      ci.op - output file

!ENDEXTRACT
!      ci.op   - diagnostic output: any limitations of the implementation
!                must be given as warnings on the output file  ci.op.

	use ci_1
	use ci_2
	use ci_6
        use ci_ck
        use ci_dat
        use ci_dat6
        use streams_mod
	implicit none

	integer, intent(out) ::  ncv, nfullv, nstrms 
	logical :: exists, error
	integer :: modeci_in, ncomps
	character(30) :: blank, head, tail, name
        external cinit_dr
        character(16), allocatable :: comp_names(:)
        character(len=100) :: mess
!  return if already called

	if( initialized ) then
	   ncv    = nc
	   nfullv = nfull
	   nstrms = nstr
	   return
	endif

	initialized = .true.
	ciinit_called = .true.

!  determine mpi status; set idproc for file names

	call isat_mpi_rank( myrank, nprocs )  

	if( nprocs == 1 ) then
	   idproc = -1
	elseif( .not.ciparam_called ) then
	   idproc = myrank
	else
	   if( ci_info_n(15) == 0 ) then
	      idproc = myrank
	   else
	      idproc = -1
	   endif
	endif

!  set op_rank based on sing_io
        if( sing_io == 1 ) then
           idproc = -1
           if ( myrank == 0 ) then
              op_rank = .true.
           else
              op_rank = .false.
           endif
        endif


!  set logical unit number lu_err for error output

	blank = repeat(' ',30)!  generate file name:  ci_P.err
	head  = blank
	head  = 'ci'

	if( ciparam_called  .and.  ci_info_n(16) >=0 ) then
	   lu_err =  ci_info_n(16)  !  use specifed lu
	else
!  open file ci_P.err  on  lu_err  for error output
	   tail = blank   
	   tail = 'err'
	   call isat_file_name( head, -1, idproc, tail, name )

	   call isat_lu( lu_err )
	   open( lu_err, file = name, status='replace', 
     1         action='write', err=5 )
	endif

!  open output file (luout, ci_P.op) and write header

        luout = 0
        if(op_rank) then
           tail = blank   
           tail = 'op'
           call isat_file_name( head, -1, idproc, tail, name )
           
           call isat_lu( luout )
           open( luout, file = name, status='replace',
     1	      action='write', err=10 )
           call write_head
        endif

!  open input file (luin, streams.in) if it exists

	blank = repeat(' ',30)!  generate file name:  ci_P.err
	head  = blank  !  generate file name:  streams_P.in
	head  = 'streams'
	tail = blank   
	tail = 'in'
        call isat_file_name( head, -1, idproc, tail, name )

	inquire( file = name, exist = exists )

	if( .not.exists ) then

           if( modeci /= 6 .and. modeci /= 7 
     1          .and. modeci /= 8 .and. modeci /=9) call
     1 isat_abort('ciinit',1,'File streams.in required but not found',
     1 chv=name )

          if(op_rank) then
           write(luout,*)' '
           write(luout,*)'Initialization from (s)dtchem without streams'
          endif
          luin = -1
	else
	   call isat_lu( luin )
	   open( luin,  file = name, action='read', err=20 )

           ! find the streams format
           strm_format = streams_format(luin, error)
           if(error) then
              call streams_ErrMsg(error=mess)
              call isat_abort('ciinit: error reading streams.in',1,
     1             mess=trim(mess))
           endif
           rewind(luin)

           ! read modeci
           if( strm_format == 1) then
              read(luin, *, end=100, err=110)  modeci_in
           else
              modeci_in = streams_modeci(luin, error)
           endif
           if(error) then
              call streams_ErrMsg(error=mess)
              call isat_abort('ciinit: error reading streams.in',1,
     1             mess=trim(mess))
           endif
          


	   if(modeci/=6 .and. modeci/=7 .and. modeci/=8 .and. 
     1            modeci/=9 ) then 
	      modeci = modeci_in
	   else
             if(op_rank) then
              write(luout,*)' '
              write(luout,*)'Initialization from (s)dtchem with streams'
             endif

	      if( modeci /= modeci_in ) then
                if(op_rank) then
	         write(luout,*)' ' 
		 write(luout,*)'WARNING: value of modeci in argument ',
     1	                       'to (s)dtchem = ', modeci
		 write(luout,*)'         value of modeci set in ',
     1	                       'streams.in        = ', modeci_in
                endif
     	      endif
     	   endif

	   if(op_rank) write(luout,*)' '
	   if(op_rank) write(luout,*)'modeci = ', modeci
	endif

!  initialize Cantera
        if (modeci > 5) call cickin

!  initialization for new streams' format
        if( strm_format == 2) then
           if(modeci == 2 .or. modeci == 3) then
              ncomps = 1
              allocate(comp_names(ncomps))
              comp_names(1) = symb_f
              call streams_init(luin, error, comp_names, ncomps)
           elseif(modeci > 5) then
              ncomps = ns + 2
              allocate(comp_names(ncomps))
              comp_names(1:ns) = snames(1:ns)
              comp_names(ns+1) = symb_press
              comp_names(ns+2) = symb_temp
              call streams_init(luin, error, comp_names, ncomps)
           else
              call streams_init(luin, error)
           endif
           ! check for error
           if(error) then
              call streams_ErrMsg(error=mess)
              call isat_abort('ciinit: error reading streams.in',1,
     1             mess=trim(mess))
           endif
           if(.not.streams_status()) then
              call streams_ErrMsg(error=mess)
              call isat_abort('ciinit: error reading streams.in',2,
     1             mess=trim(mess))
           endif
        endif

!  call appropriate cinit (in which nc, nfull and nstr are set,
!  and in which arrays are allocated)

	if( modeci == 1 ) then
	   call cinit1
	elseif( modeci == 2  .or.  modeci == 3 ) then
	   call cinit2
	elseif( modeci == 6  .or.  modeci == 7 ) then
	   call cinit6
	elseif( modeci == 8  .or.  modeci == 9 ) then 
	   call cinit_dr 
	else
	   call isat_abort('ciinit',3, mess=
     1	                     ' bad modeci=', isv = modeci )
	endif

	ncv    = nc
	nfullv = nfull
	nstrms = nstr

!  close files

	if(op_rank) write(luout,*)' '
	if(op_rank) write(luout,*)
     1	 '======  ISAT-CK initialization completed  ====== '

	if( luin >0 ) close( luin )
	if(op_rank) close( luout )

	return gas%thermo_id

!--- error conditions

5	call isat_abort('ciinit', 1, mess=
     1	           ' error opening file  ci.err' )
10	call isat_abort('ciinit', 1, mess=
     1	           ' error opening file  ci.op' )
20	call isat_abort('ciinit', 2, mess=
     1	           ' error opening file  streams.in' )
100	call isat_abort('ciinit', 3, mess=
     1	           ' hit end of file trying to read  modeci' )
110	call isat_abort('ciinit', 4, mess= 
     1	           ' error trying to read  modeci' )

!-----  internal subroutine

	contains

	subroutine write_head

	integer :: v(8)

	write(luout,*)' '
	write(luout,180) version
180	format('=====  ISAT-CK Version: ', a10,
     1   ' initialization starting  ====== ' )

	write(luout,*)' '
	call date_and_time( values=v )

	write(luout,200) v(2),v(3),v(1), v(5),v(6),v(7)

200	format('Date (MM/DD/YEAR): ', i2,'/',i2,'/',i4,
     1       ' Time: ',i2,':',i2,':',i2,'.' )

	return
	end subroutine write_head

	end subroutine ciinit

!=====================================================================

!BEGEXTRACT

       subroutine cistrm( istr, ncv, c, dpt )

!  chemistry interface routine to return composition of specified stream.

!  input:
!      istrm - index of stream ( 1 <= istrm <= nstr ) (integer)
!      ncv   - number of composition variables (integer)
!  output:
!      c     - composition vector for stream istrm (length ncv, double)
!      dpt   - density, pressure and temperature of stream istrm
!	       (length 3, double)

!ENDEXTRACT

	use ci_dat
	implicit none

	integer, intent(in)     :: istr, ncv
	real(k_dp), intent(out) :: c(ncv), dpt(3)

!  check istr

	if( istr < 1  .or.  istr > nstr ) call isat_abort('cistrm',
     1	      1, mess='bad value of istr; istr, nstr=', 
     2	      ivar = (/ istr, nstr /) )

	c   = ccstrm(:,istr)
	dpt = dptstr(:,istr)

	return
	end subroutine cistrm

!=====================================================================
!BEGEXTRACT

        subroutine cicomp( ncv, c, krep, nfullv, comp, cname )

! chemistry interface routine to return full composition comp,
!	corresponding to composition c.

! input:
!      ncv     - number of composition variables (integer)
!              - If ncv=0, then only cname is returned: other input ignored.
!      c       - composition vector (length ncv, double)
!      krep    - type of representation required
!          = 1 - express species as mole fractions
!          = 2 - express species as mass fractions
!          = 3 - express species as specific mole numbers
!      nfullv  - number of full composition variables ( = ncv + 3 )
!      comp(1) - pressure in Chemkin units	
!      comp(2) - temperature (modeci = 8 or 9, only)	
!      comp(3) - density     (modeci = 8 or 9, only)
! output:
!      comp    - full composition vector (length nfullv, double) (see below)
!      cname   - names of composition variables (nfullv, character*(*))

!  for  modeci = 6 or 7, ncv = ns + 1; nfullv = ncv + 3
!       modeci = 8 or 9, ncv = nrs + ne + 1; nfullv = ncv + 3, where
!              ns  = total no. of species
!              nrs = no. of repr. species
!              ne  = no. of elements
!	comp(i)     = species(i), i=1, ncv-1
!       comp(ncv)   = density (Chemkin units)
!       comp(ncv+1) = temperature (K)
!       comp(ncv+2) = pressure (Chemkin units)
!       comp(ncv+3) = enthalpy (Chemkin units)

!ENDEXTRACT

        use ci_dat
	use ci_2
	use ci_6
	implicit none

	integer, intent(in)       :: ncv, krep, nfullv
	real(k_dp), intent(in)    :: c(ncv)
	real(k_dp), intent(inout) :: comp(nfullv)
	character*(*), intent(out):: cname(nfullv)
        external cicmp_dr

	cname = cmpsym

	if( ncv == 0 ) return

	if( modeci == 1 ) then
           ! set comp equal to density
	   comp = dptstr(1,1)
	   return
	elseif( modeci == 2  .or.  modeci == 3 ) then
	   call cicmp2( c, comp )
	elseif( modeci == 6  .or.  modeci == 7 ) then
	   call cicmp6( c, krep, comp )
	elseif( modeci == 8  .or.  modeci == 9 ) then 
	   call cicmp_dr( c, krep, comp )  	   
	else
	   call isat_abort( 'cicomp', 2,
     1	       mess = ' not coded for modeci=', isv=modeci )
	endif

	return
	end subroutine cicomp

!=====================================================================

!BEGEXTRACT

	subroutine cirxn( t, ncv, c0, ct, dpt )

!  chemistry interface routine to return composition c(t) resulting from
!      reaction for a time t from the initial composition c(0).  Also
!      returned are density, pressure and temperature.

!  input:
!      t     - time (seconds), duration of reaction (double)
!      ncv   - number of composition variables (integer)
!      c0    - initial composition vector (length ncv, double)
!      dpt(2)- pressure in Chemkin units (required for modeci=6,7)
!  output:
!      ct    - final composition vector (length ncv, double)
!      dpt   - density, pressure and temperature (length 3, double)

!ENDEXTRACT

	use ci_1
	use ci_2
	use ci_6
	implicit none

	integer, intent(in)       :: ncv
	real(k_dp), intent(in)    :: t, c0(ncv)
	real(k_dp), intent(out)   :: ct(ncv)
	real(k_dp), intent(inout) :: dpt(3)
	integer I
        external cirxn_dr
        
	if( modeci == 1 ) then
	   call cirxn1( t, c0, ct, dpt )
	elseif( modeci == 2  .or.  modeci == 3 ) then
	   call cirxn2( t, c0, ct, dpt )
	elseif( modeci == 6 ) then
	   call cirxn6( t, c0, ct, dpt )
	elseif( modeci == 7 ) then
	   call cirxn7( t, c0, ct, dpt )
	elseif( modeci == 8 .or. modeci == 9) then 
	   call cirxn_dr( t, c0, ct, dpt )      
	else
	   call isat_abort( 'cirxn', 1, mess=' not coded for modeci=',
     1	                    isv = modeci )
	endif

	return
	end subroutine cirxn

!=====================================================================
!BEGEXTRACT

        subroutine ciconv( ncv, jd,   jp,   je,   js,
     1                          din,  pin,  ein,  spin,
     2                          kd,   kp,   ke,   ks,
     3                          dout, pout, eout, spout )

!  chemistry interface routine to perform conversion of thermochemical
!     variables for modeci = 6 or 7.

!  input:
!	ncv  - number of composition variables (=ns+1) (integer)
!	jd   - units of input density (integer):
!	     =1, CGS (g/cm^3)
!	     =2, SI (kg/m^3)
!	jp   - units of input pressure (integer):
!	     =1, standard atmospheres
!	     =2, CGS (dyne/cm^2)
!	     =3, SI (Pa)
!	je   - type of input energy variable (integer):
!	     =1, temperature, T (K)
!	     =2, sensible enthalpy, h_s (CGS) ergs/g
!	     =3, enthalpy, h (CGS) ergs/g
!	     =4, enthalpy, h (SI) J/kg
!	js   - type of species variable (integer):
!	     =1, mole fraction, X
!	     =2, mass fraction, Y
!	     =3, specific mole number, Z (mole/g)
!	din  - input density (double)
!	pin  - input pressure (double)
!	ein  - input energy variable (double)
!	spin - input species variables (length ns=ncv-1, double)
!	kd   - units of output density (integer) (=1 or 2, as jd)
!	kp   - units of output pressure (integer) (=1,2 or 3 as jp)
!	ke   - type of output energy variable (integer) (=1,2,3 or 4, as je)
!	ks   - type of species variable (integer) (=1,2, or 3, as js)
!	dout - output density (double)
!	pout - output pressure (double)
!	eout - output energy variable (double)
!	spout- output species variables (length ns, double)

!ENDEXTRACT

	use ci_dat
	use ci_cksubs

	implicit none

	integer, intent(in) :: ncv, jd, jp, je, js, kd, kp, ke, ks
	real(k_dp), intent(in)  :: din, pin, ein, spin(ns)
	real(k_dp), intent(out) :: dout, pout, eout, spout(ns)
        integer :: check = 0, info ! used in temphz
	real(k_dp) :: press, hs, h, phi(ns), y(ns)

!  check that Cantera has been initialized; check modeci and ncv

	if( .not.initialized ) call isat_abort( 'ciconv', 0,
     1 	        mess='Cantera has not been initialized' )

	if( modeci /= 6  .and. modeci /= 7 ) call isat_abort(
     1	       'ciconv', 1, mess='invalid modeci=', isv=modeci )

	if( ncv /= ns+1 ) call isat_abort(
     1	       'ciconv', 2, mess='ncv /= ns+1', isv=ncv )

!  convert density

	if( jd == 1 ) then
	   if( kd == 1 ) then
	      dout = din
	   elseif( kd == 2 ) then
	      dout = din * 1.d3
	   else
	      call isat_abort( 'ciconv', 3,
     1 	        mess='invalid jd= ', isv=jd )
	   endif
	elseif( jd == 2 ) then
	   if( kd == 1 ) then
	      dout = din * 1.d-3
	   elseif( kd == 2 ) then
	      dout = din 
	   else
	      call isat_abort( 'ciconv', 4,
     1 	        mess='invalid jd= ', isv=jd )
	   endif
	else
	   call isat_abort( 'ciconv', 5,
     1 	     mess='invalid jd= ', isv=jd )
     	endif

!  get pressure in Chemkin units

	if( jp == 1 ) then
	   press = pin * patm
	elseif( jp == 2 ) then
	   press = pin
	elseif( jp == 3 ) then
	   press = pin * 10.d0
	else
	   call isat_abort( 'ciconv', 6,
     1 	        mess='invalid jp= ', isv=jp )
     	endif

! convert pressure if necessary

	if( kp == jp ) then
	   pout = pin
	elseif( kp == 1 ) then
	   pout = press / patm
	elseif( kp == 2 ) then
	   pout = press
	elseif( kp ==3 ) then
	   pout = press * 0.1d0
	else
	   call isat_abort( 'ciconv', 7,
     1 	        mess='invalid kp= ', isv=kp )
     	endif

!  get specific mole number

	if( js == 1 ) then
	   call m2phi( spin, amolwt, ns, phi )
	elseif( js == 2 ) then
	   call y2phi( spin, amolwt, ns, phi )
	elseif( js == 3 ) then
	   phi = spin
	else
	   call isat_abort( 'ciconv', 8,
     1 	        mess='invalid js= ', isv=js )
     	endif

!  convert species if necessary

	if( ks == js ) then
	   spout = spin
	elseif( ks == 1 ) then
	   call phi2m( phi, ns, spout )
	elseif( ks == 2 ) then
	   call phi2y( phi, amolwt, ns, spout )
	elseif( ks == 3 ) then
	   spout = phi
	else
	   call isat_abort( 'ciconv', 9,
     1 	        mess='invalid ks= ', isv=ks )
     	endif

!  quick return if no energy conversion required

	if( ke == je ) then
	   eout = ein
	   return
	endif

!  obtain hs

	if( je == 1 ) then
	   call phi2y( phi, amolwt, ns, y )
	   call cthbms( ein, y, h, gas )
	   call h2hs( h, href, phi, ns, hs )
	elseif( je == 2 ) then
	   hs = ein
	elseif( je == 3 ) then
	   call h2hs( ein, href, phi, ns, hs )
	elseif( je == 4 ) then
	   call h2hs( ein, href, phi, ns, hs )
	   hs = hs * 1.d4
	else
	   call isat_abort( 'ciconv', 10,
     1 	        mess='invalid je= ', isv=je )
     	endif

!  convert from hs to required energy variable

	if( ke == 1 ) then
	   call hs2h( hs, href, phi, ns, h )
	   call temphz( h, phi, eout, check, info )
	elseif( ke == 2 ) then
	   eout = hs
	elseif( ke == 3 ) then
	   call hs2h( hs, href, phi, ns, eout )
	elseif( ke == 4 ) then
	   call hs2h( hs, href, phi, ns, eout )
	   eout = eout * 1.d-4
	else
	   call isat_abort( 'ciconv', 11,
     1 	        mess='invalid ke= ', isv=ke )
     	endif

	return
	end subroutine ciconv

!=====================================================================
!BEGEXTRACT

	subroutine cisave( nrec )

!   Force checkpointing of ISATAB table

! inout:  none

! output:
!	nrec - number of records in table (integer)

!ENDEXTRACT

	use ci_dat
	implicit none

	integer, intent(out) :: nrec

	real(k_dp) :: rinfo(50), stats(100)
	integer    :: info(100)

	if( modeci < 7 ) then  ! no ISAT table to save for modeci<7
	   nrec = -1
	   return
	endif

	call cisat( 12, info, rinfo, stats )
	call cisat( 6,  info, rinfo, stats )

     	nrec = stats(12)

	return
	end subroutine cisave

!=====================================================================

!BEGEXTRACT

	subroutine ciparam( ci_info, ci_rinfo, info, rinfo ) 

!  Specify parameters used in ISAT-CK for modeci=6 and 7, and those passed to ISATAB.
!  If called before  ciinit  then the specifed values of ci_info and ci_rinfo are used in 
!  the initialization of ISAT-CK, and the values of info and rinfo are used in the
!  initialization of ISATAB.
!  If called after the ISAT table has been initialiazed, then the ISATAB parameters are 
!  changed according to info and rinfo (and ci_info and ci_rinfo are not referenced).

! Input:
!    ci_info  - integer array (dimensioned 20)  containing parameters used in ISAT-CK - details below
!    ci_rinfo - real array    (dimensioned 20)  containing parameters used in ISAT-CK - details below
!    info     - integer array (dimensioned 100) containing parameters used in ISATAB 
!    rinfo    - real array    (dimensioned 50)  containing parameters used in ISATAB 

!    The details of info and rinfo are given in ISATAB.  All defaults are zero.

! ci_info:
!    Default values are shown in [ ] and are used if ci_info(i) is set to 0.

! ci_info(1)  - con_pr  = 1 for constant pressure                            [0]
! ci_info(2)  - con_dt  = 1 for constant time step                           [0]
! ci_info(3)  - us_rate = 1 for user-supplied reaction rate                  [0]
! ci_info(4)  - radiate = 1 for radiation                                    [0]
! ci_info(5)  - kreal   - realizability treatment                            [2]
! ci_info(6)  - ichdas  = 1, 2 perform accuracy checking on ddasac           [1]
! ci_info(7)  - ickcorr = 1 to correct thermo coeffs                         [0]
! ci_info(8)  - m_sens  - mode of determining sensitivities                  [2]
! ci_info(9)  - njacs   - number of Jacobians to store for sensitivities     [41]
! ci_info(10) - q_pade  - order of Pade approximation in expm                [5]
! ci_info(11) - exp_m   - algorithm for expm; =1 - Pade, =2 - Najeld & Havel [2]
! ci_info(12) - lu_sens > 0 diagnostic o/p from ci_sens on unit lu_sens      [-1]
! ci_info(13) - iorv    = 1 to output reaction vectors                       [0]
! ci_info(14) - op_param= 1 to output values of all parameters               [0]
! ci_info(15) - mpi_uniq= 0 for unique file names when using MPI             [0]
! ci_info(16) - lu_err >= 0, logical unit number for error output;
!                       < 0, error output goes to file ci_P.err              [-1]
! ci_info(17) - kreal_h = 1 to adjust hs and T (to conserve enthalpy)        [0]
!                           following a realizability correction
!                       = 0 make no adjustment
! ci_info(18) - kreal_t = 0 do not check temperature after realiz. correc.   [0]
!                       = 1 check and as necessary adjust T
!                       = 2 check and abort if T not in [tbadlo, tbadhi]

! ci_info(19) - m_jac = 0  divided difference to obtain Jacobin matrix       [0]
!                     = 1  analytical Jacobian by adifor
!                      if  radiation =true, m_jac is 0 regardless the value of m_jac 
!                      if  us_rate   =true, m_jac is 0 regardless the value of m_jac
! ci_info(20) - m_ckwyp = 0 use standard Cantera rates canteralib.f90        [0]
!                       = 1 use ckwyp_ext.f
!                        Kalle changed def.val from 1 to 0; How can we set it to
!                        0 if it defaults to 1 if 0?? Line 53 in ci_ck.f


! ci_rinfo:
!    Default values are shown in [ ] and are used if ci_rinfo(i) is set to 0.

! ci_rinfo(1)  = atolc    - absolute error tolerance for DDASAC            [1.d-8]
! ci_rinfo(2)  = rtolc    - relative error tolerance for DDASAC            [1.d-9]
! ci_rinfo(3)  = atols    - absolute error tolerance for sensitivities     [1.d-2]
! ci_rinfo(4)  = rtols    - relative error tolerance for sensitivities     [1.d-2]
! ci_rinfo(5)  = phiref   - reference species             (internally set if zero)
! ci_rinfo(6)  = tempref  - reference temperature         (internally set if zero)
! ci_rinfo(7)  = pressref - reference pressure            (internally set if zero)
! ci_rinfo(8)  = hsref    - reference sens. enth.(x)      (internally set if zero)
! ci_rinfo(9)  = hsfref   - reference sens. enth.(f)      (internally set if zero)
! ci_rinfo(10) = dtref    - reference time interval       (internally set if zero)
! ci_rinfo(11) = tbadlo   - lower bound on allowed temperature             [250.]
! ci_rinfo(12) = tbadhi   - upper bound on allowed temperature             [3000.]
! ci_rinfo(13) = sens_lim - threshold above which m_sens=1 is used         [2.]

!ENDEXTRACT
      use ci_dat
	use ci_dat6
	implicit none

	integer, intent(in)    :: ci_info(20), info(100)
	real(k_dp), intent(in) :: ci_rinfo(20), rinfo(50)

	integer    :: iusr(1)
	real(k_dp) :: x(1), rusr(1), f(1), dfdx(1,1), hvar(1), stats(100)

	external cirmap1
	external cirmap_dr_rcce, cirmap_dr_ice_pic 
! if called prior to ciinit, store values

	if( .not.ciinit_called ) then
	   ciparam_called = .true.
	   ci_info_n  = ci_info
	   ci_rinfo_n = ci_rinfo
	   info_n     = info
	   rinfo_n    = rinfo

	elseif( isatab_called ) then
!  change ISATAB parameters
       if (modeci == 6 .or. modeci == 7) then 
	   call isatab( 1, 1, 1, x, 1, 0, 1, cirmap1,
     1	   iusr, rusr, info, rinfo, f, dfdx, hvar, stats )
       elseif(modeci == 8)    then 
          ! 5/21/10, VH: set nh = 1, hvar(1) = T_DR(dt)
         call isatab( 1, 1, 1, x, 1, 1, 1, cirmap_dr_rcce,
     1	    iusr, rusr, info, rinfo, f, dfdx, hvar, stats )
       elseif(modeci == 9)    then 
          ! 5/21/10, VH: set nh = 1, hvar(1) = T_DR(dt)
         call isatab( 1, 1, 1, x, 1, 1, 1, cirmap_dr_ice_pic,
     1	    iusr, rusr, info, rinfo, f, dfdx, hvar, stats )
       endif

	endif
	return

	end subroutine ciparam

!=====================================================================

!BEGEXTRACT

	subroutine cisat( mode, info, rinfo, stats ) 

!    Make a special call to ISATAB

! Input:
!    mode   - special ISATAB mode
!    info   - integer array (dimensioned 100) containing parameters used in ISATAB 
!    rinfo  - real array    (dimensioned 50)  containing parameters used in ISATAB 
! Output:
!    stats  - ISATAB statistics

!    The details of mode, info, rinfo and stats are given in ISATAB.

!ENDEXTRACT
      use ci_dat 
	use ci_dat6
	implicit none

	integer,    intent(in)  :: mode, info(100)
	real(k_dp), intent(in)  :: rinfo(50)
	real(k_dp), intent(out) :: stats(100)

	integer    :: iusr(1)
	real(k_dp) :: x(1), rusr(1), f(1), dfdx(1,1), hvar(1)

	external cirmap1
	external cirmap_dr_rcce, cirmap_dr_ice_pic 

	stats = 0.d0
      if (modeci == 6 .or. modeci == 7) then 
	  call isatab( 1, mode, 1, x, 1, 0, 1, cirmap1,
     1	   iusr, rusr, info, rinfo, f, dfdx, hvar, stats )
     
      elseif(modeci == 8)   then 
        call isatab( 1, mode, 1, x, 1, 0, 1, cirmap_dr_rcce,
     1	   iusr, rusr, info, rinfo, f, dfdx, hvar, stats )
     
      elseif(modeci == 9)   then 
        call isatab( 1, mode, 1, x, 1, 0, 1, cirmap_dr_ice_pic,
     1	   iusr, rusr, info, rinfo, f, dfdx, hvar, stats )
      endif
      
	return
	end subroutine cisat

!=====================================================================

!BEGEXTRACT

       subroutine dtchem( t, kspec, nspec, spec0, press, kht, ht,
     1                    modecp, modeit, spect, temp, dens )

!  ISAT-CK routine to determine thermochemical composition after
!  isobaric, reaction for a specified time interval.

!  input:
!      t       - time interval (seconds)
!      kspec   - representation of species
!              = 1 - mole fractions
!              = 2 - mass fractions
!              = 3 - specific mole numbers
!      nspec   - number of species
!      spec0   - initial species vector (length nspec)
!      press   - pressure (Chemkin units)
!      kht     - type of second thermodynamic variable
!              = 1 - enthalpy (Chemkin units)
!              = 2 - temperature (K)
!      ht      - initial value of second thermodynamic variable
!      modecp  - not used (retained for backward compatibility)
!      modeit  = 6 - direct integration
!              = 7 - ISAT
!  output:
!      ht      - final value of second thermodynamic variable
!      spect   - final species vector (length nspec)
!      temp    - final temperature
!      dens    - final density (Chemkin units)
!  notes:
!       1/ All reals are double precision.
!       2/ Units are those used in Chemkin.
!       3/ If this routine is called, then ciinit should not be called.
!       4/ The file  streams.in  is used if it exists, but the value of
!          modeci  specified in  streams.in is ignored: the value is taken
!          from the argument modeit.
!       5/ The enthalpy changes solely due to radiative heat loss (if at all).

!ENDEXTRACT

!      6/ it is assumed that there is a Chemkin II interpreter link
!         file named  chem.bin, which describes the thermochemistry.

	use ci_dat
	use ci_cksubs

	implicit none

	integer, intent(in)       :: kspec, nspec, kht, modecp, modeit
	real(k_dp), intent(in)    :: t, spec0(nspec), press
	real(k_dp), intent(inout) :: ht
	real(k_dp), intent(out)   :: spect(nspec), temp, dens

	integer, save :: nspecp
	logical, save :: first_call=.true.
	integer       :: ncx, nfullx, nstrmx, ismall
	real(k_dp)    :: cf0(nspec+1), y(nspec), h0, hs, vsmall,
     1	                 dpt(3), ct(nspec+1)

!---------  first call, initialize  ---------------------------

	if( first_call ) then
	   first_call = .false.

!  check that ciinit has not already been called

	   if( initialized ) call isat_abort( 'dtchem', 1,
     !	       mess = 'ciinit should not be called before dtchem' )

!  check modeit and set modeci

     	   if( modeit == 6  .or.  modeit == 7 ) then
	      modeci   = modeit
	      pressref = press
	   else
	       call isat_abort( 'dtchem', 2, mess = 'bad value of modeit',
     1	                        isv = modeit )
     	   endif

	   call ciinit( ncx, nfullx, nstrmx )
	   nspecp = nsp1
	endif

!-------------  check consistency of nspec

	if( ns /= nspec ) call isat_abort( 'dtchem', 3,
     1	    mess = 'Mismatch in number of species: ns, nspec=',
     2	    ivar = (/ ns, nspec /) )

!-------  convert to specific mole number and sensible enthalpy --------

	if( kspec == 1 ) then
	   call m2phi( spec0, amolwt, nspec, cf0 )
	elseif( kspec == 2 ) then
	   call y2phi( spec0, amolwt, nspec, cf0 )
	elseif( kspec == 3 ) then
	   cf0(1:nspec) = spec0
	else
	   call isat_abort( 'dtchem', 4, mess='bad value of kspec',
     1	                    isv = kspec )
	endif

	call cireal( 2, treal, cf0(1:nspec), ismall, vsmall )

!  form enthalpy, h0

	if( kht == 1 ) then
	   h0 = ht
	elseif( kht == 2 ) then
	   call phi2y( cf0(1:nspec), amolwt, nspec, y )
	   call cthbms( ht, y, h0, gas )
        else
	   call isat_abort( 'dtchem', 5, mess='bad value of kht',
     1	                    isv = kht )
	endif

!  form sensible enthalpy, hs

	call h2hs( h0, href, cf0(1:nspec), nspec, hs )
	cf0(nspecp) = hs

	dpt(2) = press
	call cirxn( t, nspecp, cf0, ct, dpt )

!-----  convert species back

	if( kspec == 1 ) then
	   call phi2m( ct(1:nspec), nspec, spect )
	elseif( kspec == 2 ) then
	   call phi2y( ct(1:nspec), amolwt, nspec, spect ) 
	else 
	   spect = ct(1:nspec)
	endif

!------  set ht

	if( kht == 1 ) then
	   call hs2h( ct(nsp1), href, ct(1:nspec), ns, ht )
	else
	   ht = dpt(3)
	endif

!-----  set  dens and temp

	dens     = dpt(1)
	temp     = dpt(3)

	return
	end subroutine dtchem

!=====================================================================

	subroutine scistrm( istr, ncv, sc, sdpt )

	use ci_dat

	implicit none

	integer, intent(in)     :: istr, ncv
	real(k_sp), intent(out) :: sc(ncv), sdpt(3)
	real(k_dp)              ::  c(ncv),  dpt(3)

	call cistrm( istr, ncv, c, dpt )

	sc   = c
	sdpt = dpt

	return
	end subroutine scistrm

!=====================================================================

	subroutine scicomp( ncv, sc, krep, nfullv, scomp, cname )

	use ci_dat
	implicit none

	integer, intent(in)        :: ncv, krep, nfullv
	real(k_sp), intent(in)     :: sc(ncv)
	real(k_sp), intent(inout)  :: scomp(nfullv)
	character*(*), intent(out) :: cname(nfullv)

	real(k_dp) :: c(ncv), comp(nfullv)

	if( modeci == 6 .or. modeci == 7 .or. modeci == 8 
     1              .or. modeci ==9 )
     1	   comp(1) = scomp(1)
	c = sc
	call cicomp( ncv, c, krep, nfullv, comp, cname ) 

	scomp = comp

	return
	end subroutine scicomp

!=====================================================================

	subroutine scirxn( st, ncv, sc0, sct, sdpt )

	use ci_dat
	implicit none

	integer, intent(in)       ::  ncv
	real(k_sp), intent(in)    ::  st, sc0(ncv)
	real(k_sp), intent(out)   ::  sct(ncv)
	real(k_sp), intent(inout) ::  sdpt(3)
	real(k_dp)                :: t, c0(ncv), ct(ncv), dpt(3)

	t      = st
	c0     = sc0
 	if( modeci == 6  .or.  modeci == 7 .or. modeci ==8 
     1          .or. modeci ==9)
     1	  dpt(2) = sdpt(2) 
	call cirxn( t, ncv, c0, ct, dpt )  

	sct  = ct
	sdpt = dpt

	return
	end subroutine scirxn

!=====================================================================

	subroutine sciconv( ncv, jd,    jp,     je,    js,
     1	                         sdin,  s_pin,  sein,  sspin,
     2                           kd,    kp,     ke,    ks,
     3	                         sdout, s_pout, seout, sspout )

	use ci_dat6

	implicit none

	integer, intent(in)        :: ncv, jd, jp, je, js, kd, kp, ke, ks
	real(k_sp), intent(in)     :: sdin, s_pin, sein, sspin(ns)
	real(k_sp), intent(out)    :: sdout, s_pout, seout, sspout(ns)

	real(k_dp) :: din, pin, ein, spin(ns)
	real(k_dp) :: dout, pout, eout, spout(ns)

	din  = sdin
	pin  = s_pin
	ein  = sein
	spin = sspin

	call ciconv( ncv, jd, jp, je, js, din, pin,  ein,  spin,
     1               kd, kp, ke, ks, dout, pout, eout, spout )

	sdout  = dout
     	s_pout = pout
	seout  = eout
	sspout = spout

	return
	end subroutine sciconv

!=====================================================================

	subroutine sciparam( ci_info, ci_sinfo, info, sinfo ) 

	use ci_dat6
	implicit none

	integer,    intent(in) :: ci_info(20), info(100)
	real(k_sp), intent(in) :: ci_sinfo(20), sinfo(50)

	real(k_dp) :: ci_rinfo(20), rinfo(50)

	ci_rinfo = ci_sinfo
	rinfo    = sinfo

	call ciparam( ci_info, ci_rinfo, info, rinfo )

	return
	end subroutine sciparam

!=====================================================================

	subroutine scisat( mode, info, srinfo, sstats ) 

	use ci_prec
	implicit none

	integer,    intent(in)  :: mode, info(100)
	real(k_sp), intent(in)  :: srinfo(50)
	real(k_sp), intent(out) :: sstats(100)

	real(k_dp) :: rinfo(50), stats(100)

	rinfo = srinfo
	call cisat( mode, info, rinfo, stats )
	sstats = stats

	return
	end subroutine scisat

!=====================================================================

       subroutine sdtchem( st, kspec, nspec, sspec0, spress, kht, sht,
     1	                   modecp, modeit, sspect, stemp, sdens )

	use ci_dat
	use ci_dat6

	implicit none

	integer, intent(in)       :: kspec, nspec, kht, modecp, modeit
	real(k_sp), intent(in)    :: st, sspec0(nspec), spress
	real(k_sp), intent(inout) :: sht
	real(k_sp), intent(out)   :: sspect(nspec), stemp, sdens

	real(k_dp)                :: t, spec0(nspec), press
	real(k_dp)                :: ht
	real(k_dp)                :: spect(nspec), temp, dens

!----  convert to double precision

	t     = st
	spec0 = sspec0
	press = spress
	ht    = sht

!---  call dtchem

       call dtchem( t, kspec, nspec, spec0, press, kht, ht,
     1	                   modecp, modeit, spect, temp, dens )

!---  convert back

	sht    = ht
	sspect = spect
	stemp  = temp
	sdens  = dens

	return
	end subroutine sdtchem

!=====================================================================
