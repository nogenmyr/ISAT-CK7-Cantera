!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!


!  contains:  cires5, ciesub5, cijac5, cibsub5

!============================================================================

	subroutine cires5( time, n, z, dzdt, rpar, ipar,
     1                     ieform, ires )

!  subroutine to evaluate residual for 1995 ddasac
!	delta = z - dz/dt

!  input:
!	time   - not used
!	n	- number of equations (must be nsp1)
!	z      - {phi,T}
!	rpar   - not used, except rpar(n+1) = press
!	ipar   - not used, except ipar(n+1) = n+1
!       ieform - not used
!  output:
!	dzdt   - dz(t)/dt
!	ires   = 0 for successful evaluation
!  note:
!	ipar(n+1) must be set to n+1; and rpar(n+1) must be the pressure

	use ci_cksubs

	implicit none

	integer, intent(in)     :: n, ipar(n+1), ieform
	integer, intent(out)    :: ires
	real(k_dp), intent(in)  :: time, z(n), rpar(n+1)
	real(k_dp), intent(out) :: dzdt(n)
	real(k_dp) :: p, dens

	p    = rpar( ipar(n+1) )
	call cidzdt( z, p, dzdt, dens )
	ires = 0

	return
	end subroutine cires5

!============================================================================

	subroutine ciesub5

!  dummy
	use isat_abort_m
	call isat_abort('ciesub5',1,mess='called in error')

	end subroutine ciesub5

!============================================================================

	subroutine cijac5

!  dummy
	use isat_abort_m
	call isat_abort('cijac5',1,mess='called in error')

	end subroutine cijac5 


!============================================================================
! Removed along with ice-pic g_jacadf.f (not ported to Cantera)
!      subroutine cijac5adf( time, n, z, PD, rpar, ipar, ieform, ires )

!   Evaluate PD(i,j)=partial f_i/partial z_j analytically.
!   Here the  Jacobian is a square n*n matrix and calculated using ADIFOR.
!   Note that the matrix calculated by ADIFOR is the transpose of PD.

!	use ci_cksubs
!	implicit none

!	integer, intent(in)     :: n, ipar(n+1), ieform
!	integer, intent(out)    :: ires
!	real(k_dp), intent(in)  :: time, z(n), rpar(n+1)
!	real(k_dp), intent(out) :: PD(n,n)
      
!	  
!	real(k_dp) :: g_z(n,n), g_dzdt(n,n), dzdt(n), zp(n)
!	integer    :: i

!  set g_z to the identity

!	g_z = 0.d0
!	do i = 1, n
!	   g_z(i,i) = 1.d0
!	   zp(i)    = max( z(i), 0.d0 ) ! SBP 9/4/2010 - enforce z >= 0
!	end do

!	call g_cires5( n, time, n, zp, g_z, n, dzdt, g_dzdt, n,
!	1             rpar, ipar, ieform, ires)
!	
!	PD  = transpose( g_dzdt )
!	ires = 1
       
!	return
!	end subroutine cijac5adf

!============================================================================
! Removed along with ice-pic g_jacadf.f (not ported to Cantera)
!	subroutine jacobian( tau, n, z, p, T, J )

!       Returns Jacobian J = dS/dz at constant enthalpy
!	
!	use ci_cksubs

!	implicit none
!	
!	integer, intent(in)     :: n
!	real(k_dp), intent(in)  :: tau, z(n) , p, T
!	real(k_dp), intent(out) :: J(n, n)

!	! local variables
!	integer    :: i, k, ires, ipar(n+2), ieform
!	real(k_dp) :: rpar(n+2), zf(n + 1), y(n)
!	real(k_dp) :: Jf(n+1, n+1), dSdh(ns)
!        real(k_dp) :: cp, hi(n)

!	! zf = {z, T}
!	zf(1:n) = z(1:n)
!	zf(n+1) = T

!	! compute Jf at constant temperature
!        rpar(n+2) = p
!        ipar(n+2) = n+2
!        call cijac5adf( tau, n+1, zf, Jf, rpar, ipar, ieform, ires )

!        J = Jf(1:n, 1:n)
!	
        ! convert to const enthalpy
        ! J|h = J|T - dS/dh * h^T
        ! dS/dh = 1/Cp * dS/dT

        ! compute h_unnorm and cp
!	call phi2y( z(1:n), amolwt, n, y )
!	call ctcpbs( T, y, cp, gas )
!        call cthml( T, hi, gas )

!        ! dS/dh
!        dSdh(1:n) = Jf(1:n, n+1)/cp
!        
!        ! J at const. enthalpy
!        do i=1,n
!           do k=1,n
!              J(i,k)= J(i,k) - dSdh(i)*hi(k)
!           end do
!        enddo

!	return
!	end subroutine jacobian

!============================================================================



	subroutine cibsub5( t, nvar, z, dfdpj, j, rpar, ipar,
     1	                    ieform, ires )

!  dassac 1995 routine to specify df/dp(j).
!  df/dp(j) = 0 returned, corresponding to p(j) being the j-th parameter

!  input:
!	t	- time (not used)
!	nvar	- number of state variables
!	z	- value of state variables
!	j	- index of parameter
!	rpar	- not used, except rpar(n+1) = press
!	ipar	- not used, except ipar(n+1) = n+1
!	ieform	- not used
!  output:
!	dfdpj	- df(i)/dp(j)
!	ires	= 1 for successful evaluation
!		= 0 to force ddasac to compute df(i)/dp(j) numerically

      use ci_dat6
	implicit none

	integer, intent(in)     :: nvar, j, ipar(nvar+1), ieform
	integer, intent(out)    :: ires
	real(k_dp), intent(in)  :: t, z(nvar), rpar(nvar+1)
	real(k_dp), intent(out) :: dfdpj(nvar)
	
	dfdpj = 0.d0

	ires  = 1

!  for variable pressure, use finite difference for sensitivity with
!  respect to pressure

	if( .not.const_pr  .and.  j == nvar+1 ) ires = 0

	return
	end subroutine cibsub5

!============================================================================

