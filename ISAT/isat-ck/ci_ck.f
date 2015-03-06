!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

	module ci_ck

      use ci_dat
      use ci_dat6

	integer             :: con_pr, con_dt, us_rate, radiate, op_param,
     1                       mpi_uniq
	integer,    private :: ci_info_def(20), ci_info_min(20), 
     1                       ci_info_max(20)
	real(k_dp), private :: ci_rinfo_def(20)

	contains


!=========================================================================
	subroutine ci_info_set

	logical :: exists
! lines 1-3 contain, in order, the 15 components of ci_info
! lines 4-5 contain, in order, the 13 components of ci_rinfo
	namelist / cinml / 
     1   con_pr, con_dt, us_rate, radiate,
     2   kreal, ichdas, ickcorr, m_sens, njacs, q_pade, exp_m,
     3   lu_sens, iorv, op_param, mpi_uniq, kreal_h, kreal_t,
     4   atolc, rtolc, atols, rtols, phiref, tempref, pressref,
     5   hsref, hsfref, dtref, tbadlo, tbadhi, sens_lim, m_jac,
     6   m_ckwyp

	character(30) :: blank, head, tail, name

! set defaults and max and min values of ci_info

	ci_info_def = 0
	ci_info_max = 1
	ci_info_min = 0

	ci_info_def(5)  = 2
	ci_info_def(6)  = 1
	ci_info_def(8)  = 2
	ci_info_def(9)  = 41
	ci_info_def(10) = 5
	ci_info_def(11) = 2
	ci_info_def(12) = -1
	ci_info_def(16) = -1
	ci_info_def(20) = 0  ! originally 1 here

	ci_info_min(12) = -1
	ci_info_min(16) = -1

	ci_info_max(5)  = 2
	ci_info_max(6)  = 2
	ci_info_max(8)  = 3
	ci_info_max(9)  = -1
	ci_info_max(10) = -1
	ci_info_max(11) = 2
	ci_info_max(12) = -1
	ci_info_max(16) = -1
	ci_info_max(18) = 2

! set local values and info_n based on info_n and defaults

	call ci_iset( con_pr  , 1 )
	call ci_iset( con_dt  , 2 )
	call ci_iset( us_rate , 3 )
	call ci_iset( radiate , 4 )
	call ci_iset( kreal   , 5 )
	call ci_iset( ichdas  , 6 )
	call ci_iset( ickcorr , 7 )
	call ci_iset( m_sens  , 8 )
	call ci_iset( njacs   , 9 )
	call ci_iset( q_pade  , 10 )
	call ci_iset( exp_m   , 11 )
	call ci_iset( lu_sens , 12 )
	call ci_iset( iorv    , 13 )
	call ci_iset( op_param, 14 )  
	call ci_iset( mpi_uniq, 15 )  
	call ci_iset( kreal_h , 17 )
	call ci_iset( kreal_t , 18 )
        call ci_iset( m_jac ,   19 )
        call ci_iset( m_ckwyp , 20 )

! set defaults for ci_rinfo.  All values must be non-negative.

	ci_rinfo_def     = 0.d0
	ci_rinfo_def(1)  = 1.d-8
	ci_rinfo_def(2)  = 1.d-9
	ci_rinfo_def(3)  = 1.d-2
	ci_rinfo_def(4)  = 1.d-2
	ci_rinfo_def(11) = 250.
	ci_rinfo_def(12) = 3000.
	ci_rinfo_def(13) = 2.

	call ci_rset( atolc    ,1  )
	call ci_rset( rtolc    ,2  )
	call ci_rset( atols    ,3  )
	call ci_rset( rtols    ,4  )
	call ci_rset( phiref   ,5  )
	call ci_rset( tempref  ,6  )
	call ci_rset( pressref ,7  )
	call ci_rset( hsref    ,8  )
	call ci_rset( hsfref   ,9  )
	call ci_rset( dtref    ,10 )
	call ci_rset( tbadlo   ,11 )
	call ci_rset( tbadhi   ,12 )
	call ci_rset( sens_lim ,13 )

!-----------  read changes from namelist  -----------------------------

	blank = repeat(' ',30)!  generate file name:  ci_P.nml
	head  = blank
	head  = 'ci'
	tail  = blank   
	tail  = 'nml'
	call isat_file_name( head, -1, idproc, tail, name )
	inquire( file = name, exist = exists )

        if(op_rank) write(luout,*) ' '
	if( exists ) then
	   call isat_lu( lu )
	   open( lu, file = name, position = 'rewind',
     1	         action = 'read' )
	   read( lu, nml = cinml )
	   if(op_rank) 
     1         write(luout,*) 'Variables have been read from ', name
	else
	   if(op_rank) 
     1         write(luout,*) 'No variables have been read from ', name
	endif

!-----------  write to  luout  if( op_all or op_lev1)  -----------------

	if( op_param == 1 .and. op_rank ) then
	   write(luout,*) ' '
	   write(luout,*)'Values of CI constants: set in ci_info_set,',
     1	                 ' modified through ci.nml'
     	   write(luout, nml = cinml ) 
	endif

	if( con_pr == 0 ) then  ! constant pressure?
	   const_pr = .false.
	else
	   const_pr = .true.
	endif

	if( con_dt == 0 ) then  ! constant dt ?
	   const_dt = .false.
	else
	   const_dt = .true.
	endif

	if( us_rate == 0 ) then  ! user-supplied reaction rate?
	   user_rate = .false.
	else
	   user_rate = .true.
	endif

	if( radiate == 0 ) then  !  radiation?
	   radiation = .false.
	else
	   radiation = .true.
	endif

!----set other constants which cannot be changed by the user

!--------  quantities used in  cirmap1 

	treal  = 0.d0	! tolerance on realizability
	errfac = 0.1d0	! factor by which tolerances changed for acc. check
	pewarn = 1.d-3	! threshold for species-error warning
	tewarn = 1.d-1	! threshold for temperature-error warning

!--------  quantities used in  rmap2

	tolneg = 1.d-8	! tolerance for negative species values

!--------  quantities used in  rmap1_test and rmap2_test

	rmap1t  = 0.d0	! controls call to rmap1_test
	rmt1_op = 0.d0	! threshold for output
	rmap2t  = 0.d0	! controls call to rmap2_test
	rmt2_op = 0.d0	! threshold for output

!--------  quantities used in  temphy

	tlow   = 270.	! lower bound on initial guess for temperature
	thigh  = 2500.	! upper bound on initial guess for temperature
	temtol = 1.d-6	! tolerance for temperature evaluation 

	contains  !--------------------------------------------------------------

	subroutine ci_iset( iset, i )

! Set the value of iset to ci_info_n(i), or to ci_info_def(i) if ci_info_n(i)=0. 
! Reset ci_info_n(i) = iset.
! If ci_info_max(i)>0 , check that ci_info_n(i) <= ci_info_max(i).
! If ci_info_min(i)>=0, check that ci_info_n(i) >= ci_info_min(i).

	implicit none
	integer, intent(in)  :: i
	integer, intent(out) :: iset

	integer :: idef, imin, imax

	iset = ci_info_n(i)
	idef = ci_info_def(i)

	if( iset == 0 ) then
	   iset = idef
	   ci_info_n(i) = iset
	   return
	endif
	
	imin = ci_info_min(i)
	imax = ci_info_max(i)

	if( imax > 0 ) then
         if( iset > imax ) call isat_abort('ci_iset', 1,
     1           mess = 'invalid value of, ci_info(i): i = ', isv = i )
	endif

	if( imin >=0 ) then
         if( iset < imin ) call isat_abort('isat_iset', 2, 
     1           mess = 'invalid value of, ci_info(i): i = ', isv = i ) 
	endif
      
	return
	end subroutine ci_iset

	subroutine ci_rset( rset, i )

! Set the value of rset to ci_rinfo_n(i), or to ci_rinfo_def(i) if ci_rinfo_n(i)=0. 
! Reset ci_rinfo_n(i) = rset.
! Check that rest is non-negative.

	implicit none
	integer,    intent(in)  :: i
	real(k_dp), intent(out) :: rset

	rset = ci_rinfo_n(i)

	if( rset < 0.d0 ) call isat_abort('ci_rset', 1,
     1        mess = 'negative value of, ci_rinfo(i): i = ', isv = i,
     2        rsv = rset )

	if( rset == 0.d0 ) rset = ci_rinfo_def(i)
	ci_rinfo_n(i) = rset
      
	return
	end subroutine ci_rset

	end subroutine ci_info_set

!=========================================================================

	subroutine cickin

!  routine to initialize cantera

!  comments
!	The file mech.xml must have been produced by cantera previously.
!	Different versions of  ckinit  either have, or do not have, an
!       integer  iflag  as the last argument.  If the version being used
!       has this argument, use ci_ck_flag.f, otherwise use ci_ck_noflag.f
!	(Make the approriate setting in the makefile.)
!
!     This routine is called for modeci = 6,7,8,9.  As much coomon information,
!     especially from Cantera, should be intialized here.

      use ceq_system
	implicit none
	      
	integer       :: iflag, nfit, CS1(1)
	real(k_dp)    :: chmax, ruc, Bg1(1,1)
	logical       :: exists, kerr
	character(30) :: blank, head, tail, name
	
	integer,    allocatable :: ev(:,:), evt(:,:)
	real(k_dp), allocatable :: cpref(:)

	if(op_rank) then
	   write(luout, *) ' '
	   write(luout, *) 'cickin: Initializing Cantera '
	endif

	call ci_ckinit( iflag, gas )

     	if( iflag /= 0 ) call isat_abort('cickin', 2,
     1	          mess='Error in ckinit: iflag=', isv=iflag )

	if(op_rank) then
	   write(luout, *) ' '
	   write(luout, *) 'cickin: Cantera initialization successful.'
	endif

!DR-08 Perform a correction to the thermo data. It is necessary for accurate sensitivity matrix A 
! Cantera does this already
	
	call ctindx( ne, ns, nr, nfit, gas )

	if(op_rank) then
	   write(luout,*)' '
	   write(luout,*)'Number of elements  = ', ne
	   write(luout,*)'Number of species   = ', ns
	   write(luout,*)'Number of reactions = ', nr
	endif

	if( nr == 0 .and. .not.user_rate .and. op_rank ) then
	   write(luout,*)' '
	   write(luout,*)
     1	   '*******WARNING: no reaction rates specified',
     2	   '  *****************************'
     	endif
     	
 !  set variables in ci_dat6

	call ci6_alloc

	allocate( ev(ns,ne) )
	allocate( evt(ne,ns) )
	allocate( cpref(ns) )
	allocate( enames(ne) )
	allocate( snames(ns) )

!  obtain molecular weights, element vector and reaction vectors
	
	call ctwt (amolwt, gas )
	call ctncf( ne, evt, gas )
	call cthml( temp_std, href, gas)
	call ctcpml( temp_std, cpref, gas)
	call ctrp( gascon, ruc, patm )

!  obtain species and element names
	call ctsyms( 6, snames, kerr, gas )
	call ctsyme( 6, enames, kerr, gas )

	href  = href - cpref * temp_std

	ev  = transpose( evt )
	dev = ev
	devt = transpose(dev)

!  get thermo data for all species   
	call ct_thermo( ns, thermo_ns(1:ns,1:15), gas )
	
!  initialize CEQ system for unconstrained equilibrium calculations
      call ceq_sys_init( ns, ne, 0, 0, dev, CS1, Bg1, thermo_ns, 
     1                   luout, 1, sys_uc, iflag )

     	if( iflag /= 0 ) call isat_abort('cickin', 3,
     1	          mess='Error in ceq_sys_init: iflag=', isv=iflag )	
     
	return 

	end subroutine cickin
	
	end module ci_ck
