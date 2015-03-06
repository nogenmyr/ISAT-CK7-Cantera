      SUBROUTINE CTAWT (AWT, gas)
 ! Atomic weight of all species
      use cantera
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION  AWT(*)
      type(phase_t) ::  gas
      call getAtomicWeights(gas, AWT)
      RETURN
      END

      SUBROUTINE CTCOMP (IST, I, gas)
 ! Index of specie named IST (returned is I)
      use cantera
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
      CHARACTER*(*) IST
      type(phase_t) :: gas
      I = speciesIndex(gas, IST)
      RETURN
      END

      SUBROUTINE CTCPBS (T, Y, CPBMS, gas)
 ! Mass based mean specific heat of mixture (constant pressure)
      use cantera
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION  Y(*)
      type(phase_t) :: gas
      call setTemperature(gas, T)
      call setMassfractions(gas, Y)
      CPBMS = cp_mass(gas)*1.d4 ! J/kg*K to ergs/gm*K
      RETURN
      END

      SUBROUTINE CTCPML (T, CPML, gas)
 ! Mole based specific heat of all species (constant pressure)
      use cantera
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION  CPML(*)
      integer :: NSP
      double precision :: p=101325.
      double precision, allocatable :: X(:)
      type(phase_t) :: gas
      NSP = nSpecies(gas)
      allocate(X(NSP))
      X = 0.
      DO 100 K = 1, NSP
        X(K)=1.
        call setState_TPX(gas, T, p, X)
        CPML(K) = cp_mole(gas)*1.d4 !J/kmol*K to erg/mol : 1e4
        X(K)=0.
  100 CONTINUE
      RETURN
      END

      SUBROUTINE CTHBMS (T, Y, HBMS, gas)
 ! Mass based enthalpy of mixture at temperature T
      use cantera
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION  Y(*)
      type(phase_t) :: gas
      DOUBLE PRECISION :: dens = 1.0

      call setTemperature(gas, T)
      call setMassfractions(gas, Y)
!      call setState_TRY(gas, T, dens, Y) ! ideal gas, rho does not matter for H
      HBMS = enthalpy_mass(gas)*1.d4 ! J/kg to erg/gm
      RETURN
      END

      SUBROUTINE CTHML (T, HML, gas)
 ! Mole based enthalpy of all species at temperature T
      use cantera
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION  HML(*)
      integer :: NSP
      double precision :: p=101325.
      double precision, allocatable :: X(:)
      type(phase_t) :: gas
      NSP = nSpecies(gas)
      allocate(X(NSP))
      X = 0.
      DO 100 K = 1, NSP
        X(K)=1.
        call setState_TPX(gas, T, p, X)
        HML(K) = enthalpy_mole(gas)*1.d4 !J/kmol to erg/mol : 1e4
        X(K)=0.
  100 CONTINUE
      RETURN
      END

      SUBROUTINE CTHMS (T, HMS, gas)
 ! Mass based enthalpy of all species at temperature T
      use cantera
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION  HMS(*)
      integer :: NSP
      double precision :: p=101325.
      double precision, allocatable :: X(:)
      type(phase_t) :: gas
      NSP = nSpecies(gas)
      allocate(X(NSP))
      X = 0.
      DO 100 K = 1, NSP
        X(K)=1.
        call setState_TPX(gas, T, p, X)
        HMS(K) = enthalpy_mass(gas)*1.d4 !J/kg to erg/gm : 1e4
        X(K)=0.
  100 CONTINUE
      RETURN
      END

      SUBROUTINE CTINDX (MM, KK, II, NFIT, gas)
 ! MM - nElements, KK - nSpecies, II - nReactions, NFIT==7
      use cantera
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      type(phase_t) ::  gas
        MM = nElements(gas)
        KK = nSpecies(gas)
        II = nReactions(gas)
        NFIT = 7 ! not used?
      RETURN
      END

      SUBROUTINE CTINIT (IFLAG, gas)
 ! Initiation cantera thermophase object gas from file mech.xml 
      USE cantera
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      integer nsp, nrxns
      double precision :: t, p
      type(phase_t) ::  gas
      gas = importPhase('mech.xml','gas')
      KK = nSpecies(gas)
      IFLAG = 0
      if (KK .le. 0) then 
        WRITE(*,*) "Cantera failed to open mechanism: ", &
        & "Check existence of mech.xml file!"
        IFLAG=1
      endif
      call setState_TPX(gas, t, p, 'N2:1')
      RETURN
      END
 
      SUBROUTINE CTNCF (MDIM, NCF, gas)
 ! Matrix of species' elemental composition 
      use cantera
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      type(phase_t) ::  gas
      DIMENSION  NCF(MDIM,*)
      integer :: NELEM, NSP

      NELEM = nElements(gas)
      NSP = nSpecies(gas)

      DO 100 K = 1, NSP
        DO 200 M = 1, NELEM
         NCF(M,K) = nAtoms(gas, K, M)
  200   CONTINUE
  100 CONTINUE
      RETURN
      END
 
      SUBROUTINE CTNU (KDIM, NUKI, gas)
 ! Matrix of reactions' stoichiometric coefficients
      use cantera
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION  NUKI(KDIM,*)
      type(phase_t) :: gas
      integer :: NELEM, NSP

      NRXN = nReactions(gas)
      NSP = nSpecies(gas)

      DO 100 K = 1, NSP
        DO 200 I = 1, NRXN
         NUKI(K,I)=reactantStoichCoeff(gas, K, I)-productStoichCoeff(gas, K, I)
  200   CONTINUE
  100 CONTINUE
      RETURN
      END
 
      SUBROUTINE CTRHOY (P, T, Y, RHO, gas)
 ! Mass density of mixture defined by P, T, Y
      use cantera
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION  Y(*)
      type(phase_t) :: gas
      call setTemperature(gas, T)
      call setMassfractions(gas, Y)
      call setPressure(gas, p/10.) ! from dyn/cm2 to Pa
      RHO = density(gas)/1000. ! kg/m3 / 1000 

      RETURN
      END
 
      SUBROUTINE CTRP (RU, RUC, PA)
 ! Some hardcoded constants (CGS units)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      RU = 83145100.
      RUC = 19872155.
      PA = 1013250.
      
      RETURN
      END
 
      SUBROUTINE CTSYME (LOUT, ENAME, KERR, gas)
 ! Element names
      use cantera
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      CHARACTER*(*) ENAME(*)
      CHARACTER(LEN=10) :: elname
      LOGICAL KERR
      type(phase_t) :: gas
      NE = nElements(gas)
      DO 100 M = 1, NE
        call getElementName(gas, M, elname)
        ENAME(M) = elname
  100 CONTINUE
      RETURN
      END
 
      SUBROUTINE CTSYMS (LOUT, KNAME, KERR, gas)
 ! Species names
      use cantera
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      CHARACTER*(*) KNAME(*)
      CHARACTER(LEN=16) :: spname
      LOGICAL KERR
      type(phase_t) :: gas
      NSP = nSpecies(gas)
      DO 100 K = 1, NSP
        call getSpeciesName(gas, K, spname)
        KNAME(K) = spname
  100 CONTINUE
      RETURN
      END
 
      SUBROUTINE CTWT (WT, gas)
 ! Mole weights of species 
      use cantera
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION  WT(*)
      type(phase_t) ::  gas
      call getMolecularWeights(gas, WT)
      RETURN
      END
 
      SUBROUTINE CTWYP (P, T, Y, WDOT, gas)
 ! Net production rates of species in mixture (P,T,Y)
      use cantera
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION  Y(*), WDOT(*)
      type(phase_t) ::  gas
      
      call setTemperature(gas, T)
      call setMassfractions(gas, Y)
      call setPressure(gas, p/10.)

      NSP = nSpecies(gas)
      call getNetProductionRates(gas, WDOT) ! kmol/m3*s
      DO 100 K = 1, NSP
        WDOT(K) = WDOT(K)/1000.
  100 CONTINUE
      RETURN
      END

