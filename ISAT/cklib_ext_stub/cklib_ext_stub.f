! If you have your own CKWYP file, you can use it here. If the code is relying on
! the Chemkin arrays ICKWRK or RCKWRK, it cannot be used with this Cantera 
! based version of ISAT-CK7. If the code need subroutines from Chemkin, you can
! call the corresponding subroutine in canteralib.f90 instead (if the corresponding
! subroutine exists)
      SUBROUTINE CKWYP_EXT (P, T, Y, WDOT, GAS)
      use cantera
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION Y(*), WDOT(*)
      type(phase_t) :: GAS
      CALL CTWYP(P, T, Y, WDOT, GAS)
      RETURN
      END
