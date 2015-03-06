!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

	module ci_2

	use ci_dat
	use streams_mod
	implicit none

	real(k_dp), save :: trxn3
	real(k_dp), allocatable, save :: svol(:)

	contains 

!======================================================================

	subroutine cinit2

!  chemistry interface initialization routine for 
!	modeci = 2,  mixture fraction formulation, and
!       modeci = 3,  reaction progress variable formulation.

!BEGEXTRACT

!      Specification of the file  streams.in  for modeci = 2 and 3

!	  1st record: modeci
!	  2nd record: nstr  - number of streams
!	  3rd record: nfull - number of variables in full representation
! modeci=3 only: 
!	  4th record: trxn3 - reaction time scale
!	  next nstr records:  f, dens, p, T, c(1), c(2),...,c(nfull)
!	  optionally, next nfull records:
!	          names of variables in full representation (character*10)
!  notes:
!	1/ f denotes mixture fraction (modeci=2) and reaction progress
!	   variable (modeci=3).
!	2/ the values of f of the streams must be strictly increasing.
!	3/ properties are linearly interpolated in f.
!	4/ the specific volume (1/dens) is linearly interpolated, not dens.
!	5/ for modeci=3 the reaction is:
!		db/dt = -b/trxn3,  where b = ( fmax - f ) / ( fmax - fmin )

!ENDEXTRACT

	implicit none
	integer :: istr, ifull, i

!--- set nc=1; read  nstr  and  nfull; and allocate arrays

	nc       = 1
	if( strm_format == 1 ) then
	   read( luin, *, end=100, err=110) nstr
	else
	   nstr = StrmDat%nTable
	endif

	if( nstr < 2 ) call isat_abort( 'cinit2', 3,
     1	   mess= 'must have at least 2 streams; nstr=', isv = nstr )

	if( strm_format == 1) then
	   read( luin, *, end=130, err=140) nfull
	else
	   nfull = StrmDat%nSp
	endif

	if( nfull < 0 )  call isat_abort( 'cinit2', 6, 
     1	   mess= ' negative nfull= ', isv = nfull )

	call ci_alloc
	allocate( svol(nstr) )

!----------- for modeci=3, read trxn3 

	if( modeci == 3 ) then
	   if( strm_format == 1 ) then
	      read( luin, *, end=160, err=170) trxn3
	   else
	      trxn3 = StrmDat%timescale
	   endif

	   if( trxn3 < 0.d0 )  call isat_abort( 'cinit2', 9,
     1	     mess=' trxn3 must be strictly positive', rsv = trxn3 )
	endif

!---------  loop over streams to read info
	do istr = 1, nstr
	   if ( strm_format == 1 ) then
	      read( luin, *, end=210, err=220 ) ccstrm(1,istr),
     1	          dptstr(:,istr), cflstr(:,istr)
	   else
	      ccstrm(1,istr) = StrmDat%Table(1,istr)
	      dptstr(:,istr) = StrmDat%Table(2:4,istr)
	      cflstr(:,istr) = StrmDat%Table(5:,istr)
	   endif

!  check density and store specific volume
	   if( dptstr(1,istr) < 0.d0 ) call isat_abort( 'cinit2', 12,
     1	         mess='  density must be positive', isv = istr,
     2	         rsv =  dptstr(1,istr) )

	   svol(istr)  = 1. / dptstr(1,istr)

!  check that f is strictly increasing

	   if( istr > 1 ) then
	      if( ccstrm(1,istr) <= ccstrm(1,istr-1) )
     1	         call isat_abort( 'cinit2', 13,
     2            mess='mixture fraction is not strictly increasing',
     3	          rvar = ccstrm(1,1:istr) )
	   endif

	end do

!----------  loop over variables in full reresentation to read names
	if( strm_format == 1 ) then
	   do ifull = 1, nfull
	      read( luin, 265, end=270, err=270 ) cmpsym(ifull)
265	      format(a16)
	      cycle
270	      cmpsym(ifull) = 'not set   '
	   end do
	else
	   do ifull = 1, nfull
	      cmpsym(ifull) = StrmDat%Spnames(ifull)
	   enddo
	endif


!-------  output on ci.op

	write(luout,*)'cinit2: CI initialized for modeci=',modeci
	if( modeci == 3 ) write(luout,*)'cinit2: trxn3=', trxn3
	write(luout,*)' '
	write(luout,600) ccstrm(1,:)
	write(luout,610) dptstr(1,:)
	write(luout,620) dptstr(2,:)
	write(luout,630) dptstr(3,:)

	write(luout,*)' '
	write(luout,*)'full representation'
	write(luout,*)' '

	do i = 1, nfull
	   write(luout,640) i, cmpsym(i), cflstr(i,:)
	end do

	return

600	format(' f          ', 1p,20e13.4 )
610	format(' density    ', 1p,20e13.4 )
620	format(' pressure   ', 1p,20e13.4 )
630	format(' temperature', 1p,20e13.4 )
640	format(i3, 2x, a16, 1p,20e13.4 )

!-------  error conditions

100     call isat_abort( 'cinit2', 1, mess=
     1        'hit end of file trying to read  nstr' )
110     call isat_abort( 'cinit2', 2, mess=
     1        'error trying to read  nstr' )
130	 call isat_abort( 'cinit2', 4,
     1	    mess='hit end of file trying to read  nfull' )
140	 call isat_abort( 'cinit2', 5,
     1	    mess='error trying to read  nfull' )
160	call isat_abort( 'cinit2', 7,
     1	       mess=' hit end of file trying to read  trxn3' )
170	call isat_abort( 'cinit2', 8,
     1	     mess='  error trying to read  trxn3' )
210	call isat_abort( 'cinit2', 10,
     1	    mess='hit end of file trying to read istr=', isv = istr )
220	call isat_abort( 'cinit2', 11,
     1	    mess='error trying to read istr=', isv = istr )

	end subroutine cinit2

!======================================================================

	subroutine cicmp2( cc, comp )

!  chemistry interface routine to return full composition comp,
!	corresponding to composition cc.
!	This version for  modeci = 2  - mixture fraction formulation, and
!	for modeci = 3 - reaction progress variable formulation.
!	Compositions assumed to be piecewise linear functions.

!  input:
!	cc	- composition vector (double)
!  output:
!	comp	- full composition vector (double)

	implicit none

	real(k_dp), intent(in)  :: cc(nc)
	real(k_dp), intent(out) :: comp(nfull)
	real(k_dp) :: ctt, eta
	integer    :: j, jp

!  interpolate for full compositions

	ctt = cc(1)
	if( ctt <= ccstrm(1,1) ) then

	   comp = cflstr(:,1)

	elseif( ctt .ge. ccstrm(1,nstr) ) then

	   comp = cflstr(:,nstr)

	else

	   jp = 1
200	   j  = jp
	   jp = j + 1
	   if( ctt > ccstrm(1,jp) ) go to 200

	   eta  = ( ctt - ccstrm(1,j) ) /
     1	          ( ccstrm(1,jp) - ccstrm(1,j) )

	   comp = cflstr(:,j) +
     1	            eta * ( cflstr(:,jp) - cflstr(:,j) )
	endif

	return
	end subroutine cicmp2

!======================================================================

	subroutine cirxn2( t, c0, ct, dpt )

!  chemistry interface routine to return composition c(t) resulting from
!	reaction for a time t from the initial composition c(0).  Also
!	returned are density, pressure and temperature.
!	This version for  modeci = 2  - mixture fraction formulation, and
!	for modeci = 3 - reaction progress variable formulation.

!  input:
!	t     - time, duration of reaction (double)
!	c0    - initial composition vector (double)
!  output:
!	ct    - final composition vector (double)
!       dpt   - density, pressure and temperature (double)

	implicit none

	real(k_dp), intent(in)  :: t, c0(nc)
	real(k_dp), intent(out) :: ct(nc), dpt(3)
	real(k_dp) :: ctt, b

!  no change in mixture fraction

	if( modeci == 2 ) then
	   ct(1) = c0(1)
	else

!  new reaction progress variable

	   b     = ( ccstrm(1,nstr) - c0(1) ) /
     1	           ( ccstrm(1,nstr) - ccstrm(1,1) )
	   b     = b * exp( -t / trxn3 )
	   ct(1) = ccstrm(1,nstr) - b * ( ccstrm(1,nstr) - ccstrm(1,1) )
	endif

	ctt = ct(1)

	call ci_dpt2(ctt, dpt)
	return
	end subroutine cirxn2

!======================================================================

	subroutine ci_dpt2( ctt, dpt )
!       interpolate for specific volume, p and T

	real(k_dp), intent(in)  :: ctt
	real(k_dp), intent(out) :: dpt(3)

	integer :: j, jp
	real(k_dp) :: eta

	if( ctt <= ccstrm(1,1) ) then

	   dpt = dptstr(:,1)

	elseif( ctt >= ccstrm(1,nstr) ) then

	   dpt = dptstr(:,nstr)

	else
	   jp = 1
200	   j  = jp
	   jp = j + 1
	   if( ctt > ccstrm(1,jp) ) go to 200

	   eta    = ( ctt - ccstrm(1,j) ) /
     1	            ( ccstrm(1,jp) - ccstrm(1,j) )

	   dpt(1) = 1. /  ( svol(j) +
     1	            eta * ( svol(jp) - svol(j) ) )

	   dpt(2) = dptstr(2,j) +
     1	            eta * ( dptstr(2,jp) - dptstr(2,j) )

	   dpt(3) = dptstr(3,j) +
     1	            eta * ( dptstr(3,jp) - dptstr(3,j) )
	endif

	end subroutine ci_dpt2

	end module ci_2
