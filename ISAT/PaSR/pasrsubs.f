!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

        subroutine means( n, ncomp, ldf, f, fmean )

!  sum for ensemble means

!  input:
!       n       - number of particles
!       ncomp   - number of compositions
!	ldf     - leading dimension of f
!       f       - property array

!  output:
!       fmean   - ensemble mean

	implicit none
	integer, intent(in)           ::  n, ncomp, ldf
        real(kind(1.e0)), intent(in)  :: f(ldf,ncomp)
        real(kind(1.e0)), intent(out) :: fmean(ncomp)

	integer :: k, i
        real(kind(1.d0)) :: sum

        do k = 1, ncomp

           sum  = 0.
           do i = 1, n
	      sum      = sum + f(i,k)
	   end do

	   fmean(k) = sum  / n
	end do

        return
        end subroutine means

!============================================================================

        subroutine mix( dt, tmix, n, ncomp, ldf, f )

!  advance particle properties by mixing for time dt.
!  For odd i, particle i relaxes to particle i+1 on the time scale tmix.

!  input:
!       dt      - time interval
!       tmix	- reaction time scale
!       n	- number of particles
!       ncomp	- number of composition variables
!       ldf	- leading dimension of f
!	f	- composition

!  output:
!	f	- new composition

	implicit none
	integer, intent(in)             :: n, ncomp, ldf
	real(kind(1.e0)), intent(in)    :: dt, tmix
	real(kind(1.e0)), intent(inout) :: f(ldf,ncomp)

	integer          :: k, i, j
        real(kind(1.e0)) :: decay, dfj
        real(kind(1.d0)) :: exparg

!  check that n is even

        if( mod(n,2) /= 0 ) then
           write(0,*)'mix: called with odd n = ', n
           stop
        endif

        exparg   = 2. * dt / tmix
        decay    = 0.5 * ( 1. - exp( -exparg ) )

        do k = 1, ncomp
        do i = 1, n-1, 2
           j        = i + 1
           dfj      = decay * ( f(i,k) - f(j,k) )
           f(i,k)   = f(i,k) - dfj
	   f(j,k)   = f(j,k) + dfj
	end do
	end do

        return
        end subroutine mix

!============================================================================

        subroutine pickpr( npick, n, ncomp, ldf, f )

!  randomly select pairs and place them at the end of the array.
!  with probability 1/2, particles in a pair are commuted.

!  input:
!	npick	- number of pairs to be picked
!	n	- number of particles
!	ncomp	- number of compositions
!	ldf	- leading dimension of f
!	f	- array of particle compositions

!  output:
!	f	- reordered array

!  work
!	work(n)

	use isat_rnu

	implicit none

	integer, intent(in) ::  npick, n, ncomp, ldf
        real(kind(1.e0)), intent(inout) :: f(ldf,ncomp)

	integer :: nh, ipick, iplast, ip, i, j, k, jl, il, jt
        real(kind(1.e0)) :: work(n), temp

!  check that n is even

        if( mod(n,2) /= 0 ) then
           write(0,*)'pickpr: called with odd n=', n
           stop
        endif

        nh = n / 2

!  check npick

        if( npick == 0  .or.  npick == nh ) then
           return
        elseif( npick < 0  .or.  npick  >  nh ) then
           write(0,*)'pickpr: called with bad npick=', npick
           stop
        endif

        call rnu( work )

!  loop over pairs

        do ipick = 1, npick

!  select at random pair (i,j=i+1)

           iplast = nh + 1 - ipick
           ip     = 1 + work(ipick) * iplast
           j      = 2 * ip
           i      = j - 1

           jl     = 2 * iplast
           il     = jl - 1

!  commute pair or not at random

           if( work(ipick+npick) > 0.5 ) then
              jt = jl
              jl = il
              il = jt
           endif

!  commute (i,j) and (il,jl)

           do k = 1, ncomp
              temp     = f(i ,k)
              f(i,k)   = f(il,k)
              f(il,k)  = temp

              temp     = f(j ,k)
              f(j,k)   = f(jl,k)
	      f(jl,k)  = temp
	   end do
	end do

        return
        end subroutine pickpr

!============================================================================

        integer function inflow( nstr, flstrm )

!  return index of stream of next inflowing particle

!  input:
!	nstr  - number of streams
!	flstr - relative flow rates of streams

	implicit none
	integer, intent(in) :: nstr
	real(kind(1.e0)), intent(in) :: flstrm(nstr)

	integer :: i, infl
	integer, save :: ifst = 0
	real(kind(1.e0)) :: flsum, xmax
	real(kind(1.e0)), allocatable, save :: xstrm(:)

!  initialize

        if( ifst == 0 ) then
           ifst     = 1
	   allocate( xstrm(nstr) )
	   xstrm = 0.
        endif

!  determine inflow

        flsum = 0.
        xmax  = 0.

        do i = 1, nstr
           xstrm(i) = xstrm(i) + flstrm(i)
           flsum    = flsum + flstrm(i)

           if( xstrm(i) > xmax ) then
              infl  = i
              xmax  = xstrm(i)
           endif
	end do

!  check flsum

        if( flsum <= 0. ) then
           write(0,*)'inflow: flsum=', flsum
           stop
        endif

        xstrm(infl) = xstrm(infl) - flsum
        inflow      = infl

        return
        end function inflow
