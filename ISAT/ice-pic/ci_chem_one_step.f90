!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ci_chem_one_step( map, z, T, p, dt, zR, TR, A, iflag )

!  Obtain the reaction mapping {zR,TR} and its gradient A, from {z,T}
!  for time dt, for the case in which the change in composition is very
!  small (e.g., because dt is small), and hence this can be achieved 
!  using a simple, explicit, one-step method.  If this treatment is
!  inaccurate, then iflag=1 is returned, and the other output
!  is not accurate (or set).

! Input:
!   map   = {0,1} for sensitivity required (1) or not (0)
!   z     - species specific moles
!   T     - temperature
!   p     - pressure
!   dt    - mapping time 

! Output:
!   zR    - mapped species specific moles
!   TR    - mapped temperature
!   A     - mapping gradient (not set for map=0)
!   iflag = 0 for success,  = 1 for failure

! Method: 
!  1/ treat dt < 0 (abort) and dt=0 (mapping is the identity)
!  2/ take explicit Euler step for x = {z,T}
!  3/ quit if non-realizable, or if normalized change in x exceeds dxn_max 
!     (which is a specified parameter)
!  4/ take 2nd-order predictor-corrector step (x1)
!  5/ take 2nd-order mid-point step (x2)
!  6/ quit if difference x1-x2 exceeds error tolerance
!  7/ return x1 and A = identity  XXX could use exp(dt*J) instead

use ci_utils
use ci_cksubs
implicit none

integer, intent(in)     :: map
real(k_dp), intent(in)  :: z(ns), T, p, dt

integer, intent(out)    :: iflag
real(k_dp), intent(out) :: zR(ns), TR, A(ns+1,ns+1)

real(k_dp) :: z_ref, T_ref, x(ns+1), dxdt(ns+1), dxn, x1(ns+1), T1, &
              dxdt1(ns+1), dx12(ns+1), errT, errz(ns), dens, x2(ns+1)

real(k_dp), parameter   :: dxn_max = 1.d-6  ! set negative to return with iflag=1 for dt > 0

!  treat dt <= 0

if( dt < 0.d0 ) then
   call isat_abort('ci_chem_one_step', 1, mess= 'negative time step = ', rsv = dt )
   
elseif( dt == 0.d0 ) then
   zR = z
   TR = T
   if( map == 1 ) A = eye(ns+1)
   iflag = 0
   return
   
elseif( dxn_max <= 0.d0 ) then
   !  make no attempt if tolerance make success impossible
   iflag = 1
   return
endif

iflag = 1            ! anticipate failure
z_ref = maxval( z )  ! reference values
T_ref = T

x(1:ns) = z  !  x = {z,T}
x(ns+1) = T

call cidzdt( x, p, dxdt, dens )

!  quit if scaled change too large
dxn = dt * ( maxval( abs( dxdt(1:ns) ) )/z_ref + dxdt(ns+1)/T_ref )
if( dxn > dxn_max ) return

x1 = x + dxdt * dt  ! explicit Euler step

if( minval(x(1:ns)) < 0.d0 ) return  ! quit if negative species

T1 = x1(ns+1)  !  quit if temperature out of range
if ( T1 <= tbadlo  .or.  T1 >= tbadhi ) return

! predictor-corrector
call cidzdt( x1, p, dxdt1, dens )
x1 = 0.5d0 * dt * ( dxdt + dxdt1 ) !  P-C increment

! mid-point
x2 = x + 0.5d0 * dxdt * dt 
call cidzdt( x2, p, dxdt, dens )
x2 = dt * dxdt !  mid-point increment

dx12 = x1 - x2 !  difference in solutions

errT = abs(dx12(ns+1))/(T*(atolc+rtolc))
if( errT > 1.d0 ) return  ! inaccurste

errz = abs( dx12(1:ns) )/( atolc*z_ref + rtolc*z )
if( maxval( errz ) > 1.d0 ) return ! inaccurate

!XXX  write(3,'(1p,10e13.4)') T, dt, dxn, errT, maxval(errz) 

x = x + x1  ! 2nd order predictor corrector
if( minval(x) < 0.d0 ) return  !  quit if negative species

iflag = 0
zR    = x(1:ns)
TR    = x(ns+1)

if( map == 0 ) return

A = eye(ns+1)

return

end subroutine ci_chem_one_step
