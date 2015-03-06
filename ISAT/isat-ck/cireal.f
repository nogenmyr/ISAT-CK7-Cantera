!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

	subroutine cireal( mode, tol, cc, ismall, vsmall )

!  routine to check and enforce realizability.
!  method based on:
!  (a) transforming from  cc  to  phi
!  (b) setting negative phi(i) to zero
!  (c) restoring conserved variables by solving the undetermined system using 
!      Gaussian elimination with column pivoting.

!  input:
!	mode   = 0 - do nothing - return
!	       = 1 - check realizability; do not enforce
!	       = 2 - check and enforce realizability
!	tol    - tolerance: realizable if  phi(i) > -tol * phiref
!	cc     - compact representaion of composition to be checked
!  output:
!	cc     - corrected composition (mode=2 only)
!	ismall - 0 if realizable; otherwise index of non-realizable
!		   species with greatest negative phi
!	vsmall - 0. if realizable; otherwise value of greatest negative
!	            non-realizable phi
!  data:
!	itmax  - maximum number of iterations
!	elmin  - if element specific mole number is negative, reset to elmin
!	       

	use ci_cksubs

	implicit none

	integer, intent(in)       :: mode
	integer, intent(out)      :: ismall
	real(k_dp), intent(in)    :: tol
	real(k_dp), intent(inout) :: cc(ns)
	real(k_dp), intent(out)   :: vsmall

	integer    :: ipiv(ns), i, jsmall, is, iter, nrow, ie, info
	real(k_dp) :: cf(ns), a(ne,ns),
     1            dlnc(ns), cf0(ns), cft(ns), rhs(ns), 
     2            cf00(ns), xi(ne),
     3	          cfmin0, tolphi, ximin, cfmin, dlnmin, atten, fac
     	integer, save :: itmax = 5

	ismall    = 0
	vsmall    = 0.
	if( mode == 0 ) return

!===  obtain species specific mole numbers and check realizability  ====

	tolphi = -tol * phiref
!	tolphi = 0.d0

	cfmin0 = 0.
	do i   = 1, ns
	   cf(i) = cc(i)
	   if( cf(i) < cfmin0 ) then
	      cfmin0 = cf(i)
	      jsmall = i
	   endif
	end do

! SBP .gt. changed to .ge. 1/23/98

	if( cfmin0 >= tolphi ) return

	vsmall    = -cfmin0
	ismall    = jsmall

!-----  non-realizable, but do not enforce

	if( mode == 1 ) return

!======  non=realizable, enforce  =======================================

!  store cf0 and cf00 for diagnostics

	do i  = 1, ns
	   cf0(i)  = cf(i)
	   cf00(i) = cf(i)
	end do

!-----  check that element mole numbers are non-negative  --------------

	ximin     = huge(1.d0)

	xi    = 0.d0
	do ie = 1, ne
	   do is = 1, ns
	      xi(ie)    = xi(ie)+ cf(is) * dev(is,ie)
	   end do   
	   ximin     = min( xi(ie) , ximin )
	end do

!  if an element specific mole number is negative, set specific mole numbers
!  to zero for all species containing that element

	if( ximin < 0.d0 ) then
	   do ie = 1, ne
	      if( xi(ie) < 0.d0 ) then
	         do is = 1, ns
	            if( dev(is,ie) > 0.d0 ) cf(is) = 0.d0
		 end do
	      endif
	   end do

!  re-normalize cf

	   call phinm( ns, cf, amolwt )

!  re-check realizability, and reset cf0

	   cfmin0     = huge(1.d0)
	   do is      = 1, ns
	      cf0(is) = cf(is)
	      cfmin0  = min( cfmin0 , cf(is) )
	   end do

	   if( cfmin0 >= tolphi ) go to 300
	endif

!===== start of iterations over corrections to species  ====================

	iter      = 0
200	cfmin     = 0.

!  form cft = max( cf, 0 )

	do is = 1, ns
	   if( cf(is) >= 0.d0 ) then
	      cft(is)  = cf(is)
	   else
	      cfmin    = min( cfmin, cf(is) )
	      cft(is)  = 0.d0
	   endif
	end do
! 
!  give up if cfmin exceeds cfmin0

	if( cfmin < cfmin0 ) then
	   write(lu_err,*)'cireal: diverging: iter, cfmin0, cfmin=',
     1	              iter, cfmin0, cfmin
	   go to 400
	endif

!  solve underdetermined system for species corrections giving the correct
!  change in element specific mole number

	nrow  = 0
	do ie = 1, ne

	   if( xi(ie) <= 0. ) cycle

	   nrow       = nrow + 1
	   rhs(nrow)  = 0.d0

	   do is  = 1, ns
	      rhs(nrow)  = rhs(nrow) + dev(is,ie) * ( cf0(is) - cft(is) )
	      a(nrow,is) = dev(is,ie) * cft(is)
	   end do

	end do

!  solve for dcf / cft

	call gecp( nrow, ns, a, ne, rhs, dlnc, info, ipiv )

	if( info /= 0 ) then
	   write(lu_err,*)'cireal: zero pivot, iter=', iter
	   write(lu_err,*)'cireal: cfmin, cfmin0, tolphi = ',
     1	                   cfmin, cfmin0, tolphi
	   go to 400
	endif

!  determine attenuation factor needed for positive specific mole numbers

	dlnmin    = huge(1.d0)
	do is     = 1, ns
	   dlnmin = min( dlnmin, dlnc(is) )
	end do

	if( dlnmin > -1.d0 ) then
	   atten = 1.d0
	else
	   atten = -1.d0 / dlnmin
	endif

!  perform correction

	cfmin = 0.d0
	do is = 1, ns
	   if( dlnc(is) == 0.d0 ) then
	      cf(is) = cft(is)
	   else
	      fac    = max( 1.d0 + dlnc(is) * atten , 0.d0 )
	      cf(is) = cft(is) * fac
	      cfmin  = min( cfmin, cf(is) )
	   endif
	end do

!  test for realizability now being satisfied

	if( cfmin >= tolphi ) go to 300

	if( iter > itmax ) then
	   write(lu_err,*)'cireal: too many iterations, iter=', iter
	   go to 400
	endif

	iter     = iter + 1
	go to 200

!=====  convert back to compact representation ======================

300	continue
	cc = cf

	return

!======  failure - print diagnostics  ===============================

400	write(lu_err,*)'cireal: element specific mole numbers'
	write(lu_err,*)' '

	do ie = 1, ne
	   write(lu_err,420)ie, xi(ie)
	end do
420	format(i4,1p,10e13.4)

	write(lu_err,*)' '
	write(lu_err,*)'cireal: species specific mole numbers'
	write(lu_err,*)' '

	do is = 1, ns
	   write(lu_err,420)is, cf00(is), cf0(is), cf(is), cft(is)
	end do

	call isat_abort('cireal',1)
	end
