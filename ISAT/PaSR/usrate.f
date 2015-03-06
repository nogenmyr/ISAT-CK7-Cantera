!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

!BEGEXTRACT

       subroutine usrate( ns, p, t, y, wdot, gas )

!  User routine to specify reaction rates.  This version calls (and is
!  functionally equivalent to) the Chemkin routine ckwyp.

!  input:
!      ns      - number of species
!      p       - pressure (Chemkin units)
!      t       - temperature (K)
!      y       - species mass fractions
!  output:
!      wdot    - reaction rates in molar units

      use cantera
	integer ns
	double precision p, t, y(ns), wdot(ns)
      type(phase_t) :: gas

	call ctwyp( p, t, y, wdot, gas )

	return
	end
!ENDEXTRACT
