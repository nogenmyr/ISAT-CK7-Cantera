!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the x2f_mpi      !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

!BEGEXTRACT

       subroutine usrate( ns, p, t, y, liwk, lrwk,
     1			   ickwrk, rckwrk, wdot )

!  User routine to specify reaction rates.  This version calls (and is
!  functionally equivalent to) the Chemkin routine ckwyp.

!  input:
!      ns      - number of species
!      p       - pressure (Chemkin units)
!      t       - temperature (K)
!      y       - species mass fractions
!      liwk    - dimension for Chemkin integer array ickwrk
!      lrwk    - dimension for Chemkin double-precision array rckwrk
!      ickwrk  - Chemkin integer array
!      rckwrk  - Chemkin double-precision array
!  output:
!      wdot    - reaction rates in molar units

	integer ns, liwk, lrwk, ickwrk(liwk)
	double precision p, t, y(ns), rckwrk(lrwk), wdot(ns)

	call ckwyp( p, t, y, ickwrk, rckwrk, wdot )

	return
	end
!ENDEXTRACT
