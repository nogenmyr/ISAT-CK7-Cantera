!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

	module ci_1

	use ci_dat
	use streams_mod

	contains

!==========================================================================

	subroutine cinit1

!  chemistry interface initialization routine for  modeci = 1,
!  inert constant-density flow.


!BEGEXTRACT

!      Specification of the file  streams.in  for modeci = 1

!	  1st record: modeci
!	  2nd record: density

!ENDEXTRACT

	implicit none
	real(k_dp) dens

!  read density from streams.in
	if( strm_format == 1) then
	   read( luin, *, end=100, err=110) dens
	else
	   dens = StrmDat%density
	endif

!  initialize

	nc          = 1
	nfull       = 1
	nstr        = 1

	call ci_alloc

	dptstr(1,1) = dens
	dptstr(2,1) = 1.
	dptstr(3,1) = 300.
	cmpsym      = 'const dens'

	write( luout, 200 ) dens
200	format('cinit1: modeci=1 initialized; dens =', 1pe13.4)

	return
	
100	call isat_abort( 'cinit1', 1, mess=
     1	      'hit end of file trying to read  dens' )

110	call isat_abort( 'cinit1', 2, mess=
     1	      'error trying to read  dens' )

	end subroutine cinit1

!==========================================================================

	subroutine cirxn1( t, c0, ct, dpt )

!  chemistry interface routine to return composition c(t) resulting from
!	reaction for a time t from the initial composition c(0).  Also
!	returned are density, pressure and temperature.
!	This version for  modeci = 1  - inert constant-density flow.

!  input:
!	t     - time, duration of reaction (double)
!	c0    - initial composition vector (double)
!  output:
!	ct    - final composition vector (double)
!       dpt   - density, pressure and temperature (double)

	implicit none

	real(k_dp), intent(in)  :: t, c0(nc)
	real(k_dp), intent(out) :: ct(nc), dpt(3)

	ct(1)  = c0(1)

	dpt(1) = dptstr(1,1)
	dpt(2) = dptstr(2,1)
	dpt(3) = dptstr(3,1)

	return
	end subroutine cirxn1

!===========================================================================

	end module ci_1
