!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

	module ci_dat

!  data used by CI routines.

	use ci_prec
	use isat_abort_m

	implicit none
	save

!  quantities defined in the Chemistry Interface

	integer :: nc, nfull, nstr, modeci=0 
	integer :: strm_format ! streams.in format, set in ciinit


! VH - 07/15/2011 - additional variables to control I/O
	integer :: sing_io = 1 ! single input/output file (0 to disable)
	integer :: ddalog  = 0 ! disable some of the ddasac logging
	logical :: op_rank = .true. ! enable output on this rank

	real(k_dp), allocatable :: ccstrm(:,:), dptstr(:,:),
     1	     cflstr(:,:)

	character(16), allocatable :: cmpsym(:)
	character(16), allocatable :: enames(:), snames(:)

!  property names
	character(16), parameter :: symb_dens  = 'RHO             '
	character(16), parameter :: symb_temp  = 'T               '
	character(16), parameter :: symb_press = 'P               '
	character(16), parameter :: symb_enth  = 'ENTH            '
	character(16), parameter :: symb_f     = 'f               '
	character(16), parameter :: symb_elem  = 'E               '

!  quantities internal to CI routines

	integer       :: luin, luout 
	logical       :: initialized = .false.
	character(10) :: version = '7.0       '

!  The input and output logical units (luin => streams.in  and 
!  luout => ci.op) are open in ciinit only.

	contains

!==========================================================================

	subroutine ci_alloc

!  allocate arrays: nc, nfull, nstr  must be defined.

	if( .not.initialized ) call isat_abort( 'ci_alloc', 1, 
     1	                       mess = 'must call ciinit first' )
     
	if(op_rank) then
	   write(luout,*)' ' 
	   write(luout,10) nc, nfull, nstr
	   write(luout,*)' ' 
	endif
	call isat_flush(luout)
10	format('CI arrays to be allocated with: nc=', i4,
     1          ', nfull=', i4, ', nstr=', i4 )

	allocate( ccstrm(nc,nstr) )
	allocate( dptstr(3,nstr) )
	allocate( cflstr(nfull,nstr) )
	allocate( cmpsym(nfull) )

	end subroutine ci_alloc

	end module ci_dat
