subroutine temphz( h, zi, T, check, info )

!  routine to determine temperature, given enthalpy and specific moles

!  input:
!	h  - specific enthalpy of gas mixture (Chemkin/CI units)
!	zi  - species specific moles
!   check = 1   - In case T is out of range, return with T force-set to tbadlo or tbadhi
!         = 0   - In case T is out of range, quit ISAT
!  output:
!	T  - temperature (K)
!   info  =  0  - temperature found within range (tbadlo < T <= tbadhi)
!         =  1  - T > tbadhi,   t set to tbadhi
!         = -1  - T <= tbadlo,  t set to tbadlo

!  Full version of September 2010 includes pre-processing to accelerate calculations of h and cp.
! This "stub" version (for public release of ISAT) merely calls ice_temphy with equivalent parameters.
    use ci_prec
    use ci_dat6
    use isat_abort_m
    use ci_stats
    use ci_cksubs
    implicit none

    real(k_dp) :: y(ns)
    real(k_dp), intent(in)  :: h, zi(ns)
    real(k_dp), intent(out) :: T
    integer, intent(in)     :: check
    integer, intent(out)    :: info

    !============  Use old routine temphy() via ICE-PIC instead of temphz() ==
    call phi2y( zi, amolwt, ns, y )
    call ice_temphy(h, y, T, info )
    if( check == 0 .and. info /= 0 ) then
       if( info == 1 ) then
          call isat_abort( 'temphz', 0, mess = 'Input Enthalpy out of range (h > hhigh)')
       else
          call isat_abort( 'temphz', 0, mess = 'Input Enthalpy out of range (h < hhigh)')
       endif
    endif
    return

end subroutine temphz

