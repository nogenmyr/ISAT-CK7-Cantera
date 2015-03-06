!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

! isat-ck additional and extended routines are added in this file.
!=====================================================================

function cistrmIdx(cstrm) result(Idx)
  ! 
  ! Given stream name cstrm, return the stream index
  !    
  use ci_prec
  use streams_mod
  use isat_abort_m
  implicit none
  character*(*),intent(in):: cstrm
  integer                 :: Idx
  logical                 :: err
  character(len=100)      :: mess

  Idx = streams_index(cstrm,err)

  if(err) then
     call streams_ErrMsg(error=mess)
     call isat_abort('streams_index',1,mess=trim(mess))
  endif

end function cistrmIdx

!=====================================================================

subroutine ciz2ze( ncv, z, nev, ze )

  ! rotine to determine specific moles of the elements from specific moles of species

  ! input:
  !  ncv    - number of reduced composition variables (integer)
  !  z      - composition vector in specific moles (length ncv, double)
  !  nev    - number of elements

  ! output:
  !  ze     - specific moles of the elements (length nev)

  use ci_dat
  use ci_dat6

  implicit none
  external ciz2ze_dr

  integer, intent(in)        :: ncv, nev
  real(k_dp), intent(in)     :: z(ncv)
  real(k_dp), intent(out)    :: ze(nev)

  if (modeci > 7) then
     call ciz2ze_dr( ncv, z, nev, ze )
     return
  else
     ze = matmul(devt, z(1:ncv-1))
  endif
end subroutine ciz2ze

!=====================================================================

!===================================================================================
subroutine ci_dpt1( ncv, c, dpt )
  ! routine to return density and temperature for mode 1

  ! input:
  ! ncv    = number of composition variables (integer)
  ! c      = composition vector (length ncv, double)

  ! output:
  ! dpt(1) = density
  ! dpt(2) = pressure
  ! dpt(3) = temperature

  use ci_dat
  implicit none

  integer, intent(in)     :: ncv
  real(k_dp), intent(in)  :: c(ncv)
  real(k_dp), intent(out) :: dpt(3)

  dpt(1) = dptstr(1,1)
  dpt(2) = dptstr(1,2)
  dpt(3) = dptstr(1,3)

end subroutine ci_dpt1

!===================================================================================
subroutine ci_dpt23( ncv, c, dpt )
  ! routine to return density and temperature for modes 2 and 3

  ! input:
  ! ncv    = number of composition variables (integer)
  ! c      = composition vector (length ncv, double)

  ! output:
  ! dpt(1) = density
  ! dpt(2) = pressure
  ! dpt(3) = temperature

  use ci_2
  implicit none

  integer, intent(in)     :: ncv
  real(k_dp), intent(in)  :: c(ncv)
  real(k_dp), intent(out) :: dpt(3)

  call ci_dpt2( c(ncv), dpt )

end subroutine ci_dpt23

!===================================================================================
subroutine ci_dpt67( ncv, c, dpt )
  ! routine to return density and temperature for modes 6 and 7

  ! input:
  ! ncv    = number of composition variables (integer)
  ! c      = composition vector (length ncv, double)
  ! dpt(2) = pressure

  ! output:
  ! dpt(1) = density
  ! dpt(2) = pressure
  ! dpt(3) = temperature

  use ci_dat
  use ci_dat6
  use ci_cksubs
  implicit none

  integer, intent(in)     :: ncv
  real(k_dp), intent(in)  :: c(ncv)
  real(k_dp), intent(inout) :: dpt(3)

  ! local variables
  real(k_dp) :: press, temp, dens, hs, h, y(ns)
  
  hs = c(ncv)
  press = dpt(2)

  call hs2h( hs, href, c(1:ns), ns, h )
  call phi2y( c(1:ns), amolwt, ns, y )
  call temphy( h, y, temp )
  call ctrhoy( press, temp, y, dens, gas )

  dpt(1) = dens
  dpt(3) = temp
  
end subroutine ci_dpt67
!===================================================================================

subroutine ice_temphy( h, y, t, info )

  !  routine to determine temperature, given enthalpy and mass fractions

  !  input:
  !	h  - specific enthalpy of gas mixture
  !	y  - species mass fractions
  !  output:
  !	temperature (K)
  !     info =  0  - temperature found within range (tbadlo < T <= tbadhi)
  !          =  1  - T > tbadhi,   t set to tbadhi
  !          = -1  - T <= tbadlo,  tset to tbadlo
  
  !  routines called:
  !	from chemkin library

  use ci_dat6  
  implicit none
  
  real(k_dp), intent(in)  :: h, y(ns)
  integer,    intent(out) :: info
  real(k_dp), intent(out) :: t
  
  integer, save    :: nitmax = 20
  
  real(k_dp) :: ht, cp, dtemp
  integer    :: niter, last
  
  !  initial guess
  
  t = 0.5 * ( tlow + thigh )
  niter = 0
  info  = 0
  
  !  start newton iterations
  
  last  = 0
  do niter = 1, nitmax
     !  determine h(t) and cp(t)
     call cthbms( t, y, ht, gas )
     call ctcpbs( t, y, cp, gas )
     
     !  temperature increment and new temperature
     dtemp = ( h - ht ) / cp
     t     = t + dtemp
     
     !  test for convergence
     if( last .eq. 1 ) then
        if( t > tbadlo  .and.  t <= tbadhi ) return
        call ice_temphy2( h, y, t, info )
        return
     endif
     
     if( abs(dtemp) .lt. temtol ) last = 1
  enddo
  
  !  Newton's method has failed to converge -- use regula falsi
  call ice_temphy2( h, y, t, info )
  return

end subroutine ice_temphy

!==========================================================================

subroutine ice_temphy2( h, y, t, info )

  !  Routine to determine temperature, given enthalpy and mass fractions.
  !  Robust version used when  ice_temphy  fails.  Uses regula falsi.
  
  !  input:
  !	h  - specific enthalpy of gas mixture
  !	y  - species mass fractions
  !  output:
  !	temperature (K)
  !     info =  0  - temperature found within range (tbadlo < T <= tbadhi)
  !          =  1  - T > tbadhi,   t set to tbadhi
  !          = -1  - T <= tbadlo,  tset to tbadlo

  !  routines called:
  !	from chemkin library

  use ci_dat6
  use isat_abort_m
  implicit none
  
  real(k_dp), intent(in)  :: h, y(ns)
  integer,    intent(out) :: info
  real(k_dp), intent(out) :: t
  
  integer, save :: nitmax=100
  integer       :: iter
  real(k_dp)    :: tlo, thi, hlo, hhi, cpinv, hm, terr, h_tmp

  info = 0

  !  evaluate enthalpy at tbadlo and tbadhi
  tlo = tbadlo
  thi = tbadhi
  call cthbms( tlo, y, hlo, gas )
  call cthbms( thi, y, hhi, gas )
  
  h_tmp = h
  if( h < hlo) then  ! bad low temperature
     h_tmp = hlo
     t     = tlo ! set T = tbadlo
     info  = -1
     return
  elseif (h > hhi) then  ! bad high temperature
     h_tmp = hhi
     t     = thi  !  set T = tbadhi
     info  = 1
     return
  endif

  !  find temperature using regula falsi
  do iter = 1, nitmax
     cpinv = (thi-tlo) / (hhi-hlo)
     t     = tlo + (h_tmp - hlo) * cpinv
     call cthbms( t, y, hm, gas )

     !  test for convergence
     terr = abs(hm-h_tmp) * cpinv
     if( terr < temtol ) return
     
     !  replace hi or low
     if( hm > h_tmp ) then
        hhi = hm
        thi = t
     else
        hlo = hm
        tlo = t
     endif
  enddo
  
  call isat_abort( 'ice_temphy2', 2, mess = 'too many iterations')
  
end subroutine ice_temphy2

!=====================================================================
