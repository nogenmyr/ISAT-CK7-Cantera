!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module ci_sens

!  Module for calculating sensitivity coefficients from the Jacobians
!  evaluated by DDASAC.

!  S.B. Pope 6/9/2003

!  DDASAC can be used to solve ODE's either with or without sensitivity
!  coefficients.  In the latter case, the CPU time required is more than
!  in the former case, typically by an order of magnitude.  The idea of
!  the method implemented here (suggested by Rodney Fox) is to calculate
!  (approximate) sensitivity coefficients based on the Jacobians that 
!  DDASAC evaluates while solving the ODE's (without sensitivity 
!  coefficients).

!  There are 3 modes.  The mode to be used is specified in mode_sens in the
!  call to sens_init.
!     mode = 1 - sensitivity coefficients are calculated by DDASAC
!     mode = 2 - sensitivity coefficients are estimated as the product of
!                exp( J_k dt_k ), where dt_k is the k-th time interval, 
!                and J_k is the Jacobian in that time interval.
!     mode = 3 - sensitivity coefficients are estimated as exp( J(t_end) t_end),
!                i.e., based on the Jacobian at the end of the time interval.

!  Operation:

!  Phase I: storing the Jacobians (and other information about the ODE solution)
!           This is performed (for mode_sens >= 2 ) for every ODE solve.

!  1/ sens_init is called prior to using DDASAC to solve the ODE without 
!     sensitivity coefficients
!  2/ sens_store is called by DDASAC to store the Jacobian each time one is evaluated
!  3/ sens_end is called after DDASAC has completed the integration

!  Phase II: estimating the sensitivity coefficients
!           This is performed (for mode_sens >= 2 ) only when sensitivity 
!           coefficients are to be calculated

!  4/ jac_dd is called to evaluate the Jacobian at the end of the time interval, J(t_end)
!  5/ sens_store is called to store J(t_end)
!  6/ sens_match is called to check that the stored information corresponds to the
!     specified initial conditions and time interval
!  7/ sens_eval is called to evaluate the sensitivity coefficients.

!  Limitations:
!  The sensitivities with respect to pressure are not evaluated.

  use ci_dat6, only: q_pade, exp_m
  use isat_abort_m
  use ci_stats

  implicit none

  real(kind(1.d0)), private, parameter :: tscale   = 1.d5  
  logical, parameter                   :: diagnostics = .false.

  logical, private, save :: initialized = .false., init_sens = .false., &
                            ended = .false., last_set
  integer, private, save :: n, mode, njacs, kjac, jac_tot
  real(kind(1.d0)), private, save :: sens_lim, press, t_end, t_last
  real(kind(1.d0)), allocatable, private, save :: u0(:), u1(:), up1(:), ts(:), cjs(:), &
    pds(:,:)

contains

  subroutine sens_init( nin, uin, pressin, tend, mode_sens, njacin, sens_limin )  !---------

! initialize sensitivity matrices: call before ddacac

! Input:
!  nin       - number of ODE dependent variables
!  uin       - initial conditions
!  pressin   - pressure
!  tend      - end time of the integration
!  mode_sens - mode of determining sensitivity coefficients 
!            = 2 - product of exp( J dt )
!            = 3 - exp( J(tend) tend )
!  njacin    - maximum number of Jacobians to be stored (even integer >=4)
!              (not referenced for mode_sens=3)
!  sens_limin= upper acceptable limit of diagonal sensitivities

  integer, intent(in)          :: nin, mode_sens, njacin
  real(kind(1.d0)), intent(in) :: uin(nin), pressin, tend, sens_limin

! if nin or mode_sens or njacin have changed, deallocate storage
  if( initialized  .and.  ( nin /=n  .or.  mode_sens /= mode .or. &
      ( mode == 2 .and. njacs /= (max(njacin,3)/2)*2+1 ) ) )then
     deallocate(u0,u1,up1,ts,cjs,pds)
	 initialized = .false.
  endif

! on first call, or if parameters have changed, allocate storage
  if( .not.initialized ) then 
     n    = nin
	 mode = mode_sens
	 sens_lim = sens_limin

! set njacs
	 if( mode == 3 ) then
	    njacs = 1
	 else
		njacs = (max(njacin,3)/2)*2 + 1  !  make njacs odd >=3
	 endif

	 allocate( u0(n) )
	 allocate( u1(n) )
	 allocate( up1(n) )
	 allocate( ts(njacs) )
	 allocate( cjs(njacs) )
	 allocate( pds(n*n,njacs) )
	 initialized = .true.
  endif

!  store initial conditions

  u0(1:n)  = uin(1:n)
  press    = pressin
  t_end    = tend

  t_last     = 0.d0
  kjac       = 0  !  index of last value set
  last_set   = .false.
  jac_tot    = 0

  init_sens  = .true.
  ended      = .false.

  return
  end subroutine sens_init

  subroutine sens_end( nin, uin, upin )  !-----------------------------

! store u and du/dt at the end of the integration

! Input:
!  nin  - number of ODE dependent variables
!  uin  - dependent variables at the end of the step
!  upin - du/dt at the end of the step

  integer, intent(in)          :: nin
  real(kind(1.d0)), intent(in) :: uin(nin), upin(nin)

  if( .not.init_sens ) return

  if( nin /= n ) then
     write(lu_err,*)'sens_end: mismatch, n, nin = ', n, nin
	 call isat_abort('sens_end',1)
  endif

  u1(1:n)  = uin(1:n)
  up1(1:n) = upin(1:n)
  ended    = .true.

  return
  end subroutine sens_end

  subroutine sens_match( nin, uin, pressin, tend, match, uout, upout )  !-------

!  determine if the input conditions match those stored, and if so return u1 and up1

! Input:
!  nin       - number of ODE dependent variables
!  uin       - initial conditions
!  pressin   - pressure
!  tend      - end time of the integration

! Output:
!  match     = .true. if conditions match
!  uout      = u1  = u(t_end)
!  upout     = up1 = du/dt at t_end

  integer, intent(in)           :: nin
  real(kind(1.d0)), intent(in)  :: uin(nin), pressin, tend

  logical, intent(out)          :: match
  real(kind(1.d0)), intent(out) :: uout(nin), upout(nin)

  integer :: i

  match = .false.

  if( .not.ended ) return
  if( nin /= n  .or.  pressin /= press  .or.  tend /= t_end ) return

  do i = 1, n-1
     if( uin(i) /= u0(i) ) return
  end do

  if( abs(uin(n)-u0(n)) > 1.d-6 ) return  !  allow small error in T

  match      = .true.
  uout(1:n)  = u1(1:n)
  upout(1:n) = up1(1:n)

  return
  end subroutine sens_match


  subroutine sens_store( nin, tin, uin, cjin, pdin, lpd )  !------------

!  store Jacobian information:  call from ddasac

! input:
!  nin - number of ODE dependent variables
!  tin - current time
!  uin - dependent variables
!  cjin- DDASAC variable cj
!  pdin- pdin( i+(j-1)*n ) = cj * delta(i,j) - J(i,j),  
!        J(i,j)=partial(f_i)/partial(u_j).
!  lpd - size of pdin, i.e., n*n

  integer, intent(in)          :: nin, lpd
  real(kind(1.d0)), intent(in) :: tin, uin(nin), cjin, pdin(lpd)

  integer :: j, jj

  if( .not.init_sens ) return

  if( nin /= n ) then
     write(lu_err,*)'sens_store: mismatch, n, nin = ', n, nin
	 call isat_abort('sens_store',1)
  endif

  jac_tot = jac_tot + 1

! VH -- BUG FIX, 10/24/2012 
! sens_store occasionally called twice with tin = t_end 
  if( kjac > 0 ) then
     ! return if jacobian at tin already stored
     if( tin == ts(kjac) ) then
        if( .false. ) print *, 'jac at tin already stored', tin
        return
     endif
  endif
! ---- VH

!  write Jacobian - for diagnostics only
  if( .false. ) call jac_write( nin, tin, uin, cjin, pdin, lpd )

!  for mode == 3, only store final Jacobian
  if( mode ==3  .and.  tin /= t_end ) return  

! ignore Jacobian for t outside [t_last, t_beyond=ts(njacs)]
  if( tin <  t_last  .or.  tin > t_end )    return

! if necessary, compress storage by discarding even entries
  if( mode == 2  .and. kjac == njacs ) then 
     do j = 3, njacs-2, 2
	   jj = (j+1)/2
	   ts(jj)        = ts(j)
	   cjs(jj)       = cjs(j)
	   pds(1:lpd,jj) = pds(1:lpd,j)
	 end do
	 kjac = jj
  endif

  kjac = kjac + 1

! ensure kjac <= njacs (size limit)
  if( kjac > njacs ) then
     write(lu_err,*) 'sens_store: error kjac > njacs'
     call isat_abort('sens_store', 1)
  endif

  ts(kjac)        = tin
  cjs(kjac)       = cjin
  pds(1:lpd,kjac) = pdin(1:lpd)

  t_last   = tin
  last_set = .true.

  return
  end subroutine sens_store

  subroutine sens_eval( nin, a, jac_use, jac_total, sens_max )  !-----------------------

!  evaluate sensitivity coefficients from Jacobians

! Input:
!  nin     - number of independent variables in ODE
! Output:
!  a       - sensitivity coefficients du(t)/du(0)
!  jac_use - number of Jacobians used
!  jac_tot - total number of Jacobians evaluated
!  sens_max- maximum diagonal sensitivity

  integer, intent(in)           :: nin
  integer, intent(out)          :: jac_use, jac_total
  real(kind(1.d0)), intent(out) :: a(nin,nin), sens_max

  integer :: i, k, ksub
  integer, parameter :: nsub = 1
  real(kind(1.d0)) :: jac(nin,nin), jacp(nin,nin),  &
                      jdt(nin,nin), expjdt(nin,nin), frac, dtsub, maxel

  if( .not.ended ) call isat_abort('sens_eval',1, mess = 'sens_end has not been called' )

  if( nin /= n ) then
     write(lu_err,*)'sens_inc: mismatch, n, nin = ', n, nin
	 call isat_abort('sens_eval',2)
  endif

  if( .not.last_set ) call isat_abort('sens_eval',3, mess = '.not.last_set' )

  sens_max = 10.d0 * sens_lim  !  set to indicate failure

  if( mode == 3 ) then
     
	 call jac_get( 1, n, jac )
	 jac_use = 1
	 jdt = jac * t_end

	 call expm( n, jdt, a, maxel )

	 if( maxel < 0.d0 ) return

	 sens_max = max_diag( n, a )

  else  !  mode = 2
     jac_use = kjac
	 call jac_get( 1, n, jac )

     if( ts(1) == 0.d0 ) then  ! first interval
        a = 0.d0  !  set a = I
	    do i = 1, n
	       a(i,i) = 1.d0
	    end do
	 else  ! set a = exp( J(1) t(1) )
		jdt = jac*ts(1)
		call expm( n, jdt, a, maxel )
		if( maxel < 0.d0 ) return
	 endif

	 do k = 2, kjac  !  loop over intervals
	   call jac_get( k, n, jacp )

	   dtsub = (ts(k)-ts(k-1)) / nsub
	   do ksub = 1, nsub  !  sub-stepping
	      sens_max = 10.d0 * sens_lim  !  set to indicate failure
	      frac = (ksub-0.5d0) / float(nsub)
	      jdt = dtsub * ( (1.d0-frac)*jac + frac * jacp)

	      call expm( n, jdt, expjdt, maxel )
		  if( maxel < 0.d0 ) return

	      a   = matmul( expjdt, a )

		  sens_max = max_diag( n, a )
		  if( sens_max > sens_lim ) return  !  failure
	   end do

	   jac = jacp
	 end do

  endif

!  remove temperature scaling

  a(n,1:n) = a(n,1:n)*tscale
  a(1:n,n) = a(1:n,n)/tscale

  jac_total = jac_tot
   
  return

  contains  ! internal routine

  real(kind(1.d0)) function max_diag( n, a )  
!  return the maximum diagonal element of the n x n matrix a

  implicit none
  integer, intent(in) :: n
  real(kind(1.d0)), intent(in)  :: a(n,n)
  integer :: i

  max_diag = 0.d0
  do i = 1, n
     max_diag = max( max_diag, a(i,i) )
  end do

  return
  end function max_diag

  end subroutine sens_eval

  subroutine jac_get( k, nin, jac )  !------------------------------

!  recover scaled Jacobian from stored information

! input:
!  k   - index in store ( 1<= k <= njacs)
!  nin - number of ODE dependent variables
! Output
!  jac - Jacobian 

  integer, intent(in)           :: k, nin
  real(kind(1.d0)), intent(out) :: jac(nin,nin)

  integer :: i, j, l

  if( .not.initialized ) call isat_abort('jac_get',1,mess= &
                     'initialization has not been performed' )

  if( k <=0  .or.  k > kjac ) call isat_abort('jac_get',2,mess= &
                     'cannot retrieve, k= ', isv=k )

  l = 0
  do j = 1, n
    do i = 1, n
	  l = l + 1
	  jac(i,j) = -pds(l,k)
	end do
	jac(j,j) = jac(j,j) + cjs(k)
  end do

! scale temperature

  jac(n,1:n) = jac(n,1:n)/tscale
  jac(1:n,n) = jac(1:n,n)*tscale

  return
  end subroutine jac_get


  subroutine jac_write( nin, tin, uin, cjin, pdin, lpd )  !------------

!  write Jacobian information

! input:
!  nin - number of ODE dependent variables
!  tin - current time
!  uin - dependent variables
!  cjin- DDASAC variable cj
!  pdin- pdin( i+(j-1)*n ) = cj * delta(i,j) - J(i,j),  
!        J(i,j)=partial(f_i)/partial(u_j).
!  lpd - size of pdin, i.e., n*n

  integer, intent(in)          :: nin, lpd
  real(kind(1.d0)), intent(in) :: tin, uin(nin), cjin, pdin(lpd)

  logical, save    :: initial = .false.
  integer, save    :: lu=4
  integer          :: l, i, j
  real(kind(1.d0)) :: jac(nin,nin)

  if( .not.initial ) then
     open(lu,file='fjac.dat')
     initial = .true.
  endif

  l = 0
  do j = 1, nin
    do i = 1, nin
	  l = l + 1
	  jac(i,j) = -pdin(l)
	end do
	jac(j,j) = jac(j,j) + cjin
  end do

  write(lu,'(1p,10000e25.15)') tin, uin, jac

  return
  end subroutine jac_write

  subroutine jac_dd( t, n, z, p, atol, rtol, dfdz )  !---------------

!  compute the Jacobian using divided differences
!  Based on DDASAC implementation

! Input:
!   t - time
!   n - number of variables
!   z - variables
!   p - pressure
!   atol - absolute error tolerance
!   rtol - relative error tolerance

! Output:
!  dfdz  - Jacobian: dfdz(i,j) = df(i)/dz(j)

! Subroutines called:
!  cires5 to compute f = dzdt

  implicit none
  integer, intent(in)           :: n
  real(kind(1.d0)), intent(in)  :: t, z(n), p, atol, rtol

  real(kind(1.d0)), intent(out) :: dfdz(n,n)

  integer :: j, ieform, ires, ipar(n+1)
  real(kind(1.d0)) :: f(n), zp(n), fp(n), rpar(n+1), zj, del
  real(kind(1.d0)), parameter :: del1=1.d-8

  rpar    = p
  ipar    = 1
  ieform  = 0

  rpar(1) = p
  zp      = z


  call cires5( t, n, z, f, rpar, ipar, ieform, ires )

  do j = 1, n
     del   = del1 * max( abs(z(j)), abs(t*f(j)), atol+rtol*abs(z(j)) )
	 zj    = z(j)
	 zp(j) = z(j) + del
	 del   = zp(j) - zj
	 call cires5( t, n, zp, fp, rpar, ipar, ieform, ires )

	 dfdz(:,j) = (fp-f)/del
	 zp(j) = zj

  end do

  return 
  end subroutine jac_dd

subroutine expm( n, a, expa, maxel )  !--------------------------------

!  compute the matrix exponential of the n x n matrix a.

! input:
!  n    - size of a
!  a    - n x n matrix

! output
!  expa - matrix exponential, expm(a)
!  maxel < 0.  for failure (element of exp(a) too large)
!        >=0   maximum absolute element

  implicit none
  integer, intent(in)           :: n
  real(kind(1.d0)), intent(in)  :: a(n,n)
  real(kind(1.d0)), intent(out) :: expa(n,n), maxel
  
  real(kind(1.d0)) :: aa(n,n), scale(n)
  
  call routine_start(i_expm)

  aa = a
  call expm_scale( 1, n, aa, scale )
  
  if( exp_m == 1 ) then  
     call expm_pade( n, aa, expa, maxel )
  elseif( exp_m == 2 ) then
     call expm_nh( n, aa, expa, maxel )
  else
     call isat_abort('expm',1, mess='invalid exp_m = ', isv=exp_m)
  endif
  
  if( maxel < 0.d0 ) return
  
  call expm_scale( 2, n, expa, scale )
  
  maxel = max( maxval(expa), -minval(expa) ) 
  
  call routine_stop(i_expm)

  return

end subroutine expm

subroutine expm_scale( kmode, n, a, scale )  !-------------------------------

! scale variables to reduce the norm of the matrix a

! input:
!  kmode = 1 - determine the scale factors and transform a 
!        = 2 - invert the transformation 
!  n     - order of the matrix a
!  a     - matrix
!  scale - scale factors

implicit none
integer, intent(in)             :: kmode, n
real(kind(1.d0)), intent(inout) :: a(n,n), scale(n)

integer,          parameter :: pass_max = 5    !  maximum number of passes
real(kind(1.d0)), parameter :: scl_rep = 10.d0  !  threshold scale 

integer :: j, pass
real(kind(1.d0)) :: diag(n), rowmax, colmax, scl, scmax, sclim

if( kmode == 1 ) then

  do j = 1, n  !  store diagonal and remove from a
     diag(j) = a(j,j)
	 a(j,j)  = 0.d0
  end do

  scale = 1.d0
  sclim = 1.d0/float(n)

  pass_loop: do pass = 1, pass_max

  scmax = 0.d0
  do j = 1, n
     rowmax = max( maxval(a(j,:)), -minval(a(j,:)) )
     colmax = max( maxval(a(:,j)), -minval(a(:,j)) )

	 if( max( rowmax, colmax ) < sclim ) then
	    scl = 1.d0
	 elseif( rowmax > colmax ) then
	    scl = max( sclim, sqrt(rowmax*colmax) ) / rowmax
	 else
	    scl = 1.d0 / ( max( sclim, sqrt(rowmax*colmax) ) / colmax )
	 endif
    
     scmax = max( scmax, scl, 1.d0/scl )
	 scale(j) = scale(j)*scl

     a(j,:) = a(j,:) * scl
	 a(:,j) = a(:,j) / scl
  end do

  if( scmax < scl_rep ) exit

  end do pass_loop

  do j = 1, n  !  restore diagonal 
     a(j,j) = diag(j)
  end do

elseif( kmode ==2 ) then

  do j = 1, n
     a(j,:) = a(j,:) / scale(j)
	 a(:,j) = a(:,j) * scale(j)
  end do

else
   call isat_abort('expm_scale',1,mess=' bad kmode = ', isv= kmode )
endif

return

end subroutine expm_scale

!=================================================================================================

subroutine expm_nh( n, a, expa, maxel )  !--------------------------------

!  compute the matrix exponential of the n x n matrix a.

! input:
!  n    - size of a
!  a    - n x n matrix

! output
!  expa - matrix exponential, expm(a)
!  maxel < 0.  for failure (element of exp(a) too large)
!        >=0   maximum absolute element

!  S.B. Pope 6/11/03; based on Algorithm of Najfeld & Havel
!                     Adv. Appl. Math. 16, 321-375 (1995)

implicit none
integer, intent(in)           :: n
real(kind(1.d0)), intent(in)  :: a(n,n)
real(kind(1.d0)), intent(out) :: expa(n,n), maxel

logical, save      :: initialized = .false.

real(kind(1.d0)), save :: f_gamma4, t_log2, num2, num4, den2, den4

integer            :: d, j, info, ipiv(n)
real(kind(1.d0))   :: x(n,n), num(n,n), den(n,n)
real(kind(1.d0))   :: asq_norm, a_to_b, elim

! treat n <=1
if( n <=0 ) then
   call isat_abort('expm_nh',1, mess= 'non-positive n= ', isv=n)
elseif( n == 1 ) then
   expa(1,1) = exp( a(1,1) )
   return
endif

!  on first call, set constants
if( .not.initialized ) then
   f_gamma4 = 4.d0 * 0.9825211d0
   t_log2   = 2.d0 * log(2.d0)
   num2     = 4.d0 / 9.d0
   num4     = 1.d0 / 63.d0
   den2     = 1.d0 / 9.d0
   den4     = 1.d0 / 945.d0
   initialized = .true.
endif

x = matmul( a, a )  ! store A^2 in x

call infnorm( n, n, x, asq_norm )  !  asq_norm = |A^2|

!  treat  exp(0) = I

if( asq_norm == 0.d0 ) then
   expa = 0.d0
   do j = 1, n
      expa(j,j) = 1.d0
   end do
   maxel = 1.d0
   return
endif

!  determine number of squarings, d

if( asq_norm <= f_gamma4 ) then
   d = 0
else
   d = ceiling( log( asq_norm/f_gamma4 ) / t_log2 )
endif

if( diagnostics ) write(lu_err,*)'d = ', d, '  asq_norm = ', asq_norm

a_to_b = 1.d0 / 2.d0**(1+d)  ! factor to convert A to B 
x      = x * a_to_b**2       ! x contains B^2

!  form numerator (num=M) and denominator (den=N) of the Pade approximation

num = num2 * x  ! B^2 terms
den = den2 * x

x  = matmul( x, x )  !  x contains B^4

num = num + num4 * x  ! add B^4 terms
den = den + den4 * x

do j = 1, n
   num(j,j) = num(j,j) + 1.d0  !  add identity
   den(j,j) = den(j,j) + 1.d0
end do
   
x   = matmul( a, den ) * a_to_b  ! x contains B N
den = num - x  ! = M - BN
num = num + x  ! = M + BN

!  solve  exp(2B) = inv(den) * num
call dgesv( n, n, den, n, ipiv, num, n, info )
if( info /= 0 ) then
   call isat_abort('expm_nh',2, mess= 'failed to solve linear system')
endif

!  square exp(2B) d times
elim = .5*sqrt(huge(1.d0))/float(n)  ! upper limit on element
expa = num
do j = 1, d
   expa = matmul( expa, expa )

   maxel = max( maxval(expa), -minval(expa) ) ! test for large element
   if( maxel > elim ) then
	  maxel = -maxel
	  return
   endif

end do

return

end subroutine expm_nh

!=====================================================================================

subroutine expm_pade( n, a, expa, maxel )  !--------------------------------

!  compute the matrix exponential of the n x n matrix a.

! input:
!  n    - size of a
!  a    - n x n matrix

! output
!  expa - matrix exponential, expm(a)
!  maxel < 0.  for failure (element of exp(a) too large)
!        >=0   maximum absolute element

!  S.B. Pope 3/11/03; based on Algorithm 11.3.1 of Golub & Van Loan

implicit none
integer, intent(in)           :: n
real(kind(1.d0)), intent(in)  :: a(n,n)
real(kind(1.d0)), intent(out) :: expa(n,n), maxel

integer            :: q, jsq, dsign, j, k, info, ipiv(n)
real(kind(1.d0))   :: anorm, as(n,n), num(n,n), den(n,n), x(n,n), c, cx(n,n), elim

! treat n <=1
if( n <=0 ) then
   call isat_abort('expm_pade',1, mess= 'non-positive n= ', isv=n)
elseif( n == 1 ) then
   expa(1,1) = exp( a(1,1) )
   return
endif

call infnorm( n, n, a, anorm )

!  treat  exp(0) = I
if( anorm == 0.d0 ) then
   expa = 0.d0
   do j = 1, n
      expa(j,j) = 1.d0
   end do
   maxel = 1.d0
   return
endif

!  scale a: divide by 2**jsq so that infnorm(a)/2**jsq < 1/2

jsq = max(0,1+floor(log(anorm)/log(2.d0))) + 1
if( diagnostics ) write(lu_err,*)' jsq = ', jsq 
as  = a / 2.d0**jsq

!  form Pade approximation:  exp(as) = inv(den) * num
q   = q_pade
x   = as
c   = 0.5d0
den = 0.d0
dsign = 1

do j = 1, n
   den(j,j) = 1.d0
end do

num = den + c*as
den = den - c*as

do k = 2, q
   c = c * float(q-k+1)/float(k*(2*q-k+1))
   x = matmul(as,x)
   cx = c*x
   num = num + cx
   if( dsign > 0 ) then
      den = den + cx
   else
      den = den - cx
   endif
   dsign = -dsign
end do

!  solve  exp(as) = inv(den) * num
call dgesv( n, n, den, n, ipiv, num, n, info )

if( info /= 0 ) call isat_abort('expm_pade',2, mess= 'failed to solve linear system')

!  square exp(as) jsq times
elim = .5*sqrt(huge(1.d0))/float(n)  ! upper limit on element
expa = num
do j = 1, jsq
   expa = matmul( expa, expa )

   maxel = max( maxval(expa), -minval(expa) ) ! test for large element
   if( maxel > elim ) then
	  maxel = -maxel
	  return
   endif

end do

return

end subroutine expm_pade

subroutine infnorm( n, m, a, anorm )
! return the infinity norm of the n x m matrix a

  implicit none
  integer, intent(in)           :: n, m
  real(kind(1.d0)), intent(in)  :: a(n,m)
  real(kind(1.d0)), intent(out) :: anorm

  integer :: j, k
  real(kind(1.d0))::  rowsum

  anorm = 0.d0
  do j = 1, n    ! loop over rows
     rowsum = 0.d0
     do k = 1, m
        rowsum = rowsum + abs( a(j,k) )
     end do
     anorm = max( anorm, rowsum )
  end do

  return
end subroutine infnorm


end module ci_sens
