!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

	module ci_cksubs

!  This module conatins various subroutines that evaluate thermochemical
!  properties.  Some call Cantera routines.

	use ci_dat6
	use isat_abort_m

	implicit none

	contains

!==========================================================================

	subroutine cidzdt( z, p, dzdt, dens )

!  evaluate time-rate-of-change of  z(t)

!  input:
!	z    - {phi,T}
!	p    - pressure (distinct from press in common)
!  output:
!	dzdt - dz(t) / dt
!	dens - density
! SBP modified 2/5/98 to ensure realizability
! includes radiation

	implicit none

	real(k_dp), intent(in)  :: z(nsp1), p
	real(k_dp), intent(out) :: dzdt(nsp1), dens

	real(k_dp) :: y(ns), hh(ns), t, cp, tsum, divqr
	integer    :: is

!  form mass fractions (y)
      y(1:ns) = z(1:ns) * amolwt(1:ns)

!  evaluate: dens, cp, dxdt, hh

	t     = z(nsp1)
	call ctrhoy( p, t, y, dens, gas )
	call ctcpbs( t, y, cp, gas )
	if( user_rate ) then
	   call usrate( ns, p, t, y, dzdt )
	elseif( m_ckwyp == 0 ) then
	   call ctwyp( p, t, y, dzdt, gas )
	elseif( m_ckwyp == 1 ) then
	   call ckwyp_ext( p, t, y, dzdt, gas )
	else
	   call isat_abort('cidzdt', 1, mess=' bad value of m_ckwyp',
     1                    isv=m_ckwyp )
	endif

	call cthms( t, hh, gas )

!  form rates for species and sums

	tsum  = 0.
	do is = 1, ns

!  ensure realizable

	   if( z(is) <= 0.d0 ) dzdt(is) = max( dzdt(is), 0.d0 )

	   dzdt(is) = dzdt(is) / dens
	   tsum    = tsum + hh(is) * dzdt(is) * amolwt(is)
	end do

	dzdt(nsp1) = -tsum / cp

!  radiation

	if( radiation ) then
	   call cirad( p, t, z, divqr )
	   dzdt(nsp1) = dzdt(nsp1) - divqr / ( dens * cp )
	endif

	return
	end subroutine cidzdt

!==========================================================================

	subroutine ciS( z, T, p, S )

!       evaluate S = dz/dt

!  input:
!	z - Species concentration
!	T - Temperature
!	p - Pressure

!  output:
!	S - dz(t)/dt

	implicit none

	real(k_dp), intent(in)  :: z(ns), T, p
	real(k_dp), intent(out) :: S(ns)

	! local variables
	real(k_dp) :: zf(nsp1), dzdt(nsp1), dens

	! form zf = {z, T}
	zf(1:ns) = z(1:ns)
	zf(nsp1) = T

	call cidzdt(zf, p, dzdt, dens)

	S(1:ns) = dzdt(1:ns)

	return
	end subroutine ciS

!==========================================================================

	subroutine cirad( p, t, z, divqr )

!  return radiative heat loss rate

!  variables and units:
!	p    	pressure		dyne/cm**2 ( = 0.1 Pa )
!	t	temperature		K
!	z	specific mole number	kg mol / kg
!	pp	partial pressure	kPa
!	divqr	divergence of heat flux	ergs/(cm**3 sec)
!	sb	Stefan-Boltzmann const	W/(m**2 K**4)
!	ap	Absorption coefficient	1/(m kPa)

!BEGEXTRACT

!   Specification of the file rad.in

!  1st line:		nrad	number of species tabulated
!  2nd line:		tback	background temperature (K)
!  next  nrad  lines:	symrad(j), augrad(j)   (a16, 4x, 1pe13.4)
!                       symrad - Chemkin symbol of radiating species
!                       augrad - augmentation factor (=1.0 for normal usage)
!                       If augrad(j)=0, then the correspond species j
!                       will not be considered in the radiation calculation
!  next ntemp lines:	t(i), (apc(i,j),j=1,nrad) 
!                       t(i) -- temperature (K)
!                       apc(i,j) -- Plank Mean Absorption Coeff.  1/(m kPa) 

!  (Radiation can be supressed by setting nrad=0)

!ENDEXTRACT

	implicit none

	real(k_dp), intent(in)  :: p, t, z(nsp1)
	real(k_dp), intent(out) :: divqr

	character(16)    :: symb(ns), symrad(ns), cdummy
	real(k_dp)       :: augrad(ns), apin(ns), pp(ns), rdummy,
     1	                    tlast, pkpa, zsum, eta, apis, ap
	real(k_dp), save :: tback
	real(k_dp), allocatable, save :: apc(:,:), trad(:), aug(:)
	integer          :: ip(ns), j, k, i, ii, nrad, idummy,
     1	                    luin, luout, ifound, is

	integer, save    :: nradsp=1, nwarn=10, iwarn=0, ntemp
	logical, save    :: initialized=.false.
	logical          :: kerr, exists 
	character(30)    :: blank, head, tail, name
	real(k_dp), save :: sb = 5.669d-08

!=============  quick return if no radiating species  ======================

	divqr = 0.
	if( nradsp == 0 ) return

!============  initialization on first call  ===============================

	if( initialized ) go to 300
	initialized = .true.
	allocate( aug(ns) )

	blank = repeat(' ',30)!  generate file name:  rad_P.out
	head  = blank
	head  = 'rad'
	tail = blank   
	tail = 'out'
	call isat_file_name( head, -1, idproc, tail, name )

	call isat_lu( luout )
	open( luout, file = name, action='write', status='replace' )

!  no radiation assumed if the file  rad.in  does not exist

	head  = blank !  generate file name:  rad_P.in
	head  = 'rad'
	tail = blank   
	tail = 'in'
	call isat_file_name( head, -1, idproc, tail, name )

	inquire( file = name, exist = exists )
	if( .not.exists ) then
	   write(luout,*)'radiation file not found: ', name
	   write(luout,*)'no radiation assumed'
	   nradsp = 0
	   return
	endif

	call isat_lu( luin )
	open( luin,  file = name, action='read' )

!  read and check number of radiating species, nrad

	read(luin, *, end=100, err=110) nrad
	if( nrad == 0 ) then
	   nradsp = 0
	   return
	elseif( nrad < 0 ) then
	   call isat_abort('cirad', 0, mess=' nrad must be positive',
     1	                   isv = nrad )
	elseif( nrad > ns ) then
	   call isat_abort('cirad', 1, mess=' nrad must be less than ns',
     1	                   isv = nrad )
	endif

!  read and check background temperature, tback

	read(luin, *, end=120, err=130) tback
	if( tback < 0.d0 ) call isat_abort('cirad', 2,
     1	                     mess=' tback must be positive', rsv=tback )

!  read  symrad(k), augrad(k)

	do k = 1, nrad
	   read(luin, 30, end=140, err=150) symrad(k), augrad(k)
	end do
30	format(a16, 4x, 1pe13.4)

!  determine pointer from rad species to chemkin species

	call ctsyms( 6, symb, kerr, gas )

	do k = 1, nrad
	   inner_loop: do j = 1, ns
	      if( symrad(k) == symb(j) ) then
	         ip(k) = j
	         exit inner_loop
	      endif
	      ip(k)    = 0
	   end do inner_loop
	end do

	nradsp  = 0
	write( luout, *)'radiating species and augmentation factors'
	do k = 1, nrad
	   if( ip(k) /= 0 ) then
	      write(luout,600) symrad(k), augrad(k)
600	      format( a16, 1pe13.4)
	      if( augrad(k) /= 0. ) nradsp = nradsp + 1
	   endif
	end do

	ifound = 0
	do k = 1, nrad
	   if( ip(k) == 0 ) then
	      if( ifound == 0 ) then
	         write( luout, *)'radiating species not found in chemkin'
	         ifound = 1
	      endif
	      write(luout,600) symrad(k), augrad(k)
	   endif
	end do

!  set augmentation factors

	aug = 0.

	do k = 1, nrad
	   if( ip(k) /= 0 ) aug( ip(k) ) = augrad( k )
	end do

!  count number of table entries

	ntemp = 0
	do
	   read(luin, *, end = 65, err = 165 ) rdummy
	   ntemp = ntemp + 1
	end do
65	continue
	if( ntemp < 2 ) call isat_abort('cirad', 3,
     1	     mess='table must have at least 2 entries', isv=ntemp )

!  allocate arrays for table

	allocate( apc(ntemp,ns) )
	allocate( trad(ntemp) )
	apc = 0.

!  re-position file
	
	rewind( luin )
	read( luin, *, end=210, err=210 ) idummy
	read( luin, *, end=210, err=210 ) rdummy
	do k = 1, nrad
	   read( luin, *, end=210, err=210 ) cdummy
	end do

!  read in table

	tlast  = -1.

	do i = 1, ntemp
	   read(luin, *, end = 160, err = 170 ) trad(i),
     1	       (apin(k),k=1,nrad)

	   if( trad(i) <= tlast ) call isat_abort('cirad', 4,  
     1	    mess='temperature must be strictly increasing', rsv=trad(i) )
	   tlast = trad(i)

	   do k = 1, nrad
	      if( ip(k) /= 0 ) apc(i,ip(k) ) = apin(k)
	   end do
	end do

	write(luout,*)'radiation table read:'
	write(luout,620) trad(1)
	write(luout,622) trad(ntemp)
	write(luout,624) ntemp
620	format(' minimum temperature = ', 1pe13.4 )
622	format(' maximum temperature = ', 1pe13.4 )
624	format(' number of entries   = ', i4 )

	close( luin )
	close( luout )

!=================  subsequent calls ===============================

300	continue

!  determine partial pressures (in kPa)

	pkpa = p * 1.e-4

	zsum = 0.
	do is = 1, ns
	   zsum      = zsum + z(is)
	end do

	do is = 1, ns
	   pp(is)    = pkpa * z(is) / zsum
	end do

!  determine interpolation coefficients

	if( t <= trad(1) ) then
	   iwarn = iwarn + 1
	   if( iwarn < nwarn ) write(lu_err,*)'cirad: low temp',t
	   ii  = 1
	   eta = 0.
	   go to 340
	elseif( t >= trad(ntemp) ) then
	   iwarn = iwarn + 1
	   if( iwarn < nwarn ) write(lu_err,*)'cirad: high temp',t
	   ii  = ntemp-1
	   eta = 1.
	   go to 340
	endif

	do i = 2, ntemp
	   if( t <= trad(i) ) then
	      ii = i - 1
	      eta = ( t - trad(ii) ) / ( trad(i) - trad(ii) )
	      go to 340
	   endif
	end do
	call isat_abort('cirad', 5,  mess='should not be here' )

!  interpolate and sum for mean absorption coefficient

340	ap  = 0.
	do is = 1, ns
	   if( aug(is) /= 0. ) then
	      apis = apc(ii,is) + eta * ( apc(ii+1,is) - apc(ii,is) )
	      ap   = ap + aug(is) * pp(is) * apis
	   endif
	end do

!  evaluate   divqr  in SI units
	divqr = 4. * sb * ap * ( t**4 - tback**4 )

!  convert to CGS
	divqr = divqr * 10.0

	return

!==============  read errors  =======================================

100	call isat_abort('cirad', 11, mess='hit end reading nrad')
110	call isat_abort('cirad', 12, mess='error reading nrad')
120	call isat_abort('cirad', 13, mess='hit end reading tback')
130	call isat_abort('cirad', 14, mess='error reading tback')
140	call isat_abort('cirad', 15, mess='hit end reading symrad', isv=k)
150	call isat_abort('cirad', 16, mess='error reading symrad', isv=k)
160	call isat_abort('cirad', 17, mess='hit end reading temperature')  
165	call isat_abort('cirad', 18, mess='error reading temperature')  
170	call isat_abort('cirad', 19, mess='error reading temperature')  
210	call isat_abort('cirad', 20, mess='error repositioning file ')  

	end subroutine cirad
!==========================================================================

      subroutine hs2h( hs, href, phi, ns, h )

!  routine to determine enthalpy from sensible enthalpy

!  input:
!	hs   - sensible enthalpy
!	href - reference enthalpy of species
!	phi  - specific mole numbers
!	ns   - number of species
!  output:
!	h    - enthalpy

	implicit none

	integer, parameter :: k_dp = kind(1.d0)
	integer, intent(in)     :: ns
	real(k_dp), intent(in)  :: hs, href(ns), phi(ns)
	real(k_dp), intent(out) :: h

      h = hs + sum( href(1:ns)*phi(1:ns) )
      return  

	return 
	end subroutine hs2h

!==========================================================================

	subroutine hs2h_n( hs, href_n, phi, nrc, h )

!Laniu  routine to determine enthalpy from sensible enthalpy based on nominal sensible of formation

!  input:
!	hs   - sensible enthalpy
!	href - reference enthalpy of species
!	phi  - specific moles of reduced compositions
!	nrc  - number of reduced compositions (nrc=nrs+ne)
!  output:
!	h    - enthalpy

	implicit none

	integer, parameter :: k_dp = kind(1.d0)
	integer, intent(in)     :: nrc
	real(k_dp), intent(in)  :: hs, href_n(nrc), phi(nrc)
	real(k_dp), intent(out) :: h

      h = hs + sum( href_n(1:nrc)*phi(1:nrc) )
      return  

	return 
	end subroutine hs2h_n

!==========================================================================

	subroutine h2hs( h, href, phi, ns, hs )

!  routine to determine sensible enthalpy from enthalpy

!  input:
!	h    - enthalpy
!       href - reference enthalpy of species
!	phi  - specific mole numbers
!	ns   - number of species
!  output:
!	hs   - sensible enthalpy

	implicit none

	integer, parameter :: k_dp = kind(1.d0)
	integer, intent(in)     :: ns
	real(k_dp), intent(in)  :: h, href(ns), phi(ns)
	real(k_dp), intent(out) :: hs

      hs = h - sum( href(1:ns)*phi(1:ns) )
      return  

	return
	end subroutine h2hs
	
	
!==========================================================================

	subroutine h2hs_n( h, href_n, phi, nrc, hs )

! Laniu routine to determine approximate sensible enthalpy contained in the represented specis 
!  from enthalpy 

!  input:
!	h      - enthalpy
!     href_n - reference enthalpy of represented species and elements (nominal species)
!	phi    - specific mole numbers of represented species and elements
!	nrc    - number of nominal species nrc=nrs +ne
!  output:
!	hs   - sensible enthalpy

	implicit none

	integer, parameter :: k_dp = kind(1.d0)
	integer, intent(in)     :: nrc
	real(k_dp), intent(in)  :: h, href_n(nrc), phi(nrc)
	real(k_dp), intent(out) :: hs

      hs = h - sum( href_n(1:nrc)*phi(1:nrc) )

	return
	end subroutine h2hs_n

!==========================================================================

	subroutine phinm( ns, phi, amolwt )

!  normalize the specific mole number phi

!  input:
!	ns	- number of species
!	phi	- specific mole numbers
!	amolwt	- molecular weights
!  output:
!	phi	- normalized specific mole numbers
	use isat_abort_m
	implicit none
	integer, parameter :: k_dp = kind(1.d0)
	integer, intent(in)       :: ns
	real(k_dp), intent(inout) :: phi(ns)
	real(k_dp), intent(in)    :: amolwt(ns)

	real(k_dp) :: sum_y

!  form sum for normalization

      sum_y = sum( phi(1:ns)*amolwt(1:ns) )
      
!  normalize

	if( sum_y > 0.d0 ) then
	   phi   = phi / sum_y
	   return
	else
	   write(lu_err,*)'phinm: non-positive sum=',sum_y
	   write(lu_err,*)phi
	   write(lu_err,*)amolwt
	   call isat_abort( 'phinm', 1, mess='non-positive sum' )
	endif

	end subroutine phinm

!==========================================================================

	subroutine phi2y( phi, amolwt, ns, y )

!  routine to determine mass fractions from specific mole numbers

!  input:
!	phi    - specific mole numbers
!	amolwt - molecular weights
!	ns     - number of species
!  output:
!	y     - species mass fractions

!  comment:  specific mole numbers do not need to be correctly normalized
	use isat_abort_m
	implicit none

	integer, parameter :: k_dp = kind(1.d0)
	integer :: i
	integer, intent(in)     :: ns
	real(k_dp), intent(in)  :: phi(ns), amolwt(ns)
	real(k_dp), intent(out) :: y(ns)

	real(k_dp) :: sum_y

	y(1:ns) = phi(1:ns) * amolwt(1:ns)
	sum_y   = sum( y(1:ns) )

	if( sum_y <= 0.d0 ) then 
	   call isat_abort( 'phi2y', 1, 
     1	           mess='non-positive species sum')
	endif
	y = y / sum_y

	return
	end subroutine phi2y

!==========================================================================

	subroutine y2phi( y, amolwt, ns, phi )

!  routine to determine specific mole numbers from mass fractions

!  input:
!	y      - species mass fractions
!	amolwt - molecular weights
!	ns     - number of species
!  output:
!	phi    - specific mole numbers

!  comment:  mass fractions do not need to be correctly normalized

	use isat_abort_m
	implicit none

	integer, parameter :: k_dp = kind(1.d0)
	integer, intent(in)     :: ns
	real(k_dp), intent(in)  :: y(ns), amolwt(ns)
	real(k_dp), intent(out) :: phi(ns)

	real(k_dp) :: sum_y

	sum_y = sum( y(1:ns) )

	if( sum_y <= 0.d0 ) 
     1	    call isat_abort( 'y2phi', 1, mess='negative sum')

	phi(1:ns) = y(1:ns) / ( sum_y * amolwt(1:ns) )

	return
	end subroutine y2phi

!==========================================================================

	subroutine phi2m( phi, ns, y )

!  routine to determine mole fractions from specific mole numbers

!  input:
!	phi    - specific mole numbers
!	ns     - number of species
!  output:
!	y      - mole fractions

!  comment:  specific mole numbers do not need to be correctly normalized

	use isat_abort_m
	implicit none

	integer, parameter      :: k_dp = kind(1.d0)
	integer, intent(in)     :: ns
	real(k_dp), intent(in)  :: phi(ns)
	real(k_dp), intent(out) :: y(ns)

	real(k_dp) :: sum_y

	sum_y = sum( phi(1:ns) )
	
	if( sum_y <= 0.d0 ) 
     1	    call isat_abort( 'phi2m',1, mess='negative sum')

	y = phi / sum_y

	return
	end subroutine phi2m

!========================================================================

	subroutine m2phi( y, amolwt, ns, phi )

!  routine to determine specific mole numbers from mole fractions

!  input:
!	y      - mole fraction
!	amolwt - molecular weights
!	ns     - number of species
!  output:
!	phi    - specific mole number

!  comment:  mole fractions do not need to be correctly normalized

	use isat_abort_m
	implicit none

	integer, parameter      :: k_dp = kind(1.d0)
	integer, intent(in)     :: ns
	real(k_dp), intent(in)  :: y(ns), amolwt(ns)
	real(k_dp), intent(out) :: phi(ns)

	real(k_dp) :: sum_y

	sum_y = sum( y(1:ns) * amolwt(1:ns) ) 

	if( sum_y <= 0.d0 ) 
     1      call isat_abort( 'm2phi',1, mess='negative sum')

	phi = y(1:ns) / sum_y

	return
	end subroutine m2phi

!==========================================================================

	subroutine temphy( h, y, t )

!  routine to determine temperature, given enthalpy and mass fractions

!  input:
!	h  - specific enthalpy of gas mixture
!	y  - species mass fractions
!  output:
!	temperature (K)
!  notes:
!	program stops if computed temperature is outside range
!	(tbadlo, tbadhi) which is set in data.
!     fix_t0 introduced 1/6/09 to eliminate history effects

!  routines called:
!	from chemkin library

	implicit none

	real(k_dp), intent(in)  :: h, y(ns)
	real(k_dp), intent(out) :: t

	real(k_dp), save :: tlast=1000.d0
	integer, save    :: nitmax = 20

	real(k_dp) :: ht, cp, dtemp
	integer    :: niter, last
	logical    :: fix_t0 = .true.

!  initial guess

      if( fix_t0 ) then
         t = 0.5 * ( tlow + thigh )
      else
	   if( tlast > tlow  .and.  tlast < thigh ) then
	      t = tlast  
	   else
	      t = 0.5 * ( tlow + thigh )
	   endif
	endif

	niter = 0

!  start newton iterations

	last  = 0
100	niter = niter + 1

!  determine h(t) and cp(t)

	call cthbms( t, y, ht, gas )
	call ctcpbs( t, y, cp, gas )

!  temperature increment and new temperature

	dtemp = ( h - ht ) / cp
	t     = t + dtemp

!  test for convergence

	tlast = t

	if( last .eq. 1 ) then
	   if( t > tbadlo  .and.  t <= tbadhi ) return
	   call temphy2( h, y, t )
	endif

	if( abs(dtemp) .lt. temtol ) last = 1
	if( niter .lt. nitmax ) go to 100

!  Newton's method has failed to converge -- use regula falsi

	call temphy2( h, y, t )

	return
	end subroutine temphy

!==========================================================================

	subroutine temphy2( h, y, t )

!  Routine to determine temperature, given enthalpy and mass fractions.
!  Robust version used when  temphy  fails.  Uses regula falsi.

!  input:
!	h  - specific enthalpy of gas mixture
!	y  - species mass fractions
!  output:
!	temperature (K)
!  notes:
!	program stops if computed temperature is outside range
!	(tbadlo, tbadhi) which is set in data.

!  routines called:
!	from chemkin library

	use isat_abort_m
	implicit none

	real(k_dp), intent(in)  :: h, y(ns)
	real(k_dp), intent(out) :: t

	integer, save :: nitmax=100
	integer       :: iter, i
	real(k_dp)    :: tlo, thi, hlo, hhi, cpinv, hm, terr

!  evaluate enthalpy at tbadlo and tbadhi

	tlo = tbadlo
	thi = tbadhi
	call cthbms( tlo, y, hlo, gas )
	call cthbms( thi, y, hhi, gas )

	if( h < hlo  .or.  h > hhi ) then
	   write(lu_err,*)'temphy2: temperature out of range'
	   write(lu_err,*)'temphy2: tbadlo, tbadhi', tbadlo, tbadhi
	   write(lu_err,*)'temphy2: h, hlo, hhi', h, hlo, hhi
	   write(lu_err,*)'species mass fractions'
	   do i = 1, ns
	      write(lu_err,600) i, y(i)
600	      format(1x,i5,1p,e13.4)
	   end do

	   call isat_abort( 'temphy2',1, mess='h outside range')
	endif

!  find temperature using regula falsi

	iter = 0
100	iter = iter + 1

	cpinv = (thi-tlo) / (hhi-hlo)
	t     = tlo + (h - hlo) * cpinv
	call cthbms( t, y, hm, gas )

!  test for convergence

	terr = abs(hm-h) * cpinv
	if( terr < temtol ) return

!  test for too many iterations

	if( iter > nitmax ) call isat_abort( 'temphy2', 2, 
     1                   mess = 'too many iterations')

!  replace hi or low

	if( hm > h ) then
	   hhi = hm
	   thi = t
	else
	   hlo = hm
	   tlo = t
	endif

	go to 100

	end subroutine temphy2
	 
!==========================================================================

	end module ci_cksubs
