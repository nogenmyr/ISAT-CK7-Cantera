      INTEGER I
      DOUBLE PRECISION SOLVQDRT
      DOUBLE PRECISION K1
      DOUBLE PRECISION K3
      DOUBLE PRECISION K4
      DOUBLE PRECISION P
      DOUBLE PRECISION OOLD
      DOUBLE PRECISION OHOLD
      DOUBLE PRECISION HO2OLD
      DOUBLE PRECISION ACH3
      DOUBLE PRECISION BCH3
      DOUBLE PRECISION CCH3
      DOUBLE PRECISION DN
      DOUBLE PRECISION FN
      DOUBLE PRECISION DHNO
      DOUBLE PRECISION FHNO
      DOUBLE PRECISION GHNO
      DOUBLE PRECISION DNH
      DOUBLE PRECISION GNH
      DOUBLE PRECISION DNH2
      DOUBLE PRECISION DNH3
      DOUBLE PRECISION    EPS
      LOGICAL CONVERGED
C
      EPS = 1.0E-03
      CONVERGED = .FALSE.

      K1 = K(R1F) / K(R1B)
      K3 = K(R3F) / K(R3B)
      K4 = K(R4F) / K(R4B)

C    CALCULATE CONCENTRATIONS OF OH, O, AND HO2 FROM H
 
      OOLD = 0.0
      OHOLD = K(R3B) * C(SH) * C(SH2O) / 
     &     CTCHZERO( K(R3F) * C(SH2) )
      HO2OLD = 0.0
      C(SO) = OOLD 
      C(SOH) = OHOLD
      C(SHO2) = HO2OLD
 
      DO 100 I = 0, 100, 1
 
      C(SO) = SOLVQDRT(
     &	K(R6F) * 2.0 * M(MM1),
     &	K(R1B) * C(SOH)
     &	+ K(R2F) * C(SH2)
     &	+ K(R4B) * C(SH2O)
     &	+ K(R12) * C(SHO2)
     &	+ K(R22) * C(SCO) * M(MM1)
     &	+ K(R111) * C(SC2H2)
     &	+ K(R112) * C(SC2H2),
     &	K(R1F) * C(SO2) * C(SH)
     &	+ K(R2B) * C(SH) * C(SOH)
     &	+ K(R4F) * C(SOH) * C(SOH)
     &	+ K(R6B) * 2.0 * C(SO2) * M(MM1)
     &	+ K(R11) * C(SHO2) * C(SH)
     &	)
 
      C(SOH) = SOLVQDRT(
     &	K(R4F) * 2.0,
     &	K(R1B) * C(SO)
     &	+ K(R2B) * C(SH)
     &	+ K(R3F) * C(SH2)
     &	+ K(R7F) * C(SH) * M(MM1)
     &	+ K(R13F) * C(SHO2)
     &	+ K(R20F) * C(SCO)
     &	+ K(R134F) * C(SC2H4),
     &	K(R1F) * C(SO2) * C(SH)
     &	+ K(R2F) * C(SH2) * C(SO)
     &	+ K(R3B) * C(SH) * C(SH2O)
     &	+ K(R4B) * 2.0 * C(SO) * C(SH2O)
     &	+ K(R7B) * C(SH2O) * M(MM1)
     &	+ K(R12) * C(SHO2) * C(SO)
     &	+ K(R13B) * C(SO2) * C(SH2O)
     &	+ K(R20B) * C(SH) * C(SCO2)
     &	+ K(R21) * C(SCO) * C(SHO2)
     &	)
      C(SHO2) = (
     &	K(R8F) * C(SH) * C(SO2) * M(MM1)
     &	+ K(R13B) * C(SO2) * C(SH2O)
     &	) / ( CTCHZERO( 
     &	K(R8B) * M(MM1)
     &	+ K(R9) * C(SH)
     &	+ K(R10) * C(SH)
     &	+ K(R11) * C(SH)
     &	+ K(R12) * C(SO)
     &	+ K(R13F) * C(SOH)
     &	+ K(R21) * C(SCO) 
     &	+ 1.0E6 
     &	) )
	   
      IF ( DABS( ( C(SO) - OOLD ) / CTCHZERO( C(SO) ) ) .LT. EPS .AND.
     & DABS( ( C(SOH) - OHOLD ) / CTCHZERO( C(SOH) ) ) .LT. EPS .AND.
     & DABS( ( C(SHO2) - HO2OLD ) / CTCHZERO( C(SHO2) ) ) .LT. EPS ) 
     & THEN
      CONVERGED = .TRUE.
      GOTO 101
      ENDIF

      OOLD =  C(SO)
      OHOLD = C(SOH)
      HO2OLD = C(SHO2) 
 100  CONTINUE

 101  IF ( CONVERGED .EQV. .FALSE. ) THEN
         WRITE(900,*)  
         WRITE(901,*)  C(SO),C(SOH),C(SHO2)
 900     FORMAT('#warning: iteration O,OH,HO2 did not converge')
 901     FORMAT('#c[O] = ',f10.8,', c[OH] = ',f10.8,', 
     &       c[HO2] = ',f10.8)
C	 If not converged use simpler expresiions for O,OH,HO2
         C(SOH) = C(SH2O) * C(SH) / CTCHZERO( K3 * C(SH2) ) 
         P = K4 / CTCHZERO( C(SH2O) ) 
         C(SO) = P * C(SOH) * C(SOH) 
      	
         C(SHO2) = ( K(R8F) * C(SH) * C(SO2) * M(MM1) ) / 
     &        ( CTCHZERO( K(R8B) * M(MM1)
     &        + K(R9) * C(SH)
     &        + K(R10) * C(SH)
     &        + K(R11) * C(SH)
     &        + K(R12) * C(SO)
     &        + K(R13F) * C(SOH) ) )
      ENDIF
	   

C    CALCULATE THE CH3 CONCENTRATION
      ACH3 = 2.0 * K(R60) + 
     &    2.0 * K(R36F) * ( K(R164) * C(SH) + K(R166) * C(SOH)  ) / 
     &    CTCHZERO( K(R36B) + K(R164) * C(SH) + K(R166) * C(SOH)  ) 
      BCH3 = K(R53) * C(SO) + K(R84B) * C(SH2) + 
     &    K(R44B) * C(SH) + K(R34F) * C(SH) + K(R55) * C(SOH) +
     &    K(R86B) * C(SH2O) + K(R90B) * C(SOH) *
     &    ( K(R91) * C(SH) + K(R93) * C(SOH)  ) / 
     &    CTCHZERO( K(R90B) + K(R91) * C(SH) + K(R93) * C(SOH)  ) 
      CCH3 = ( K(R34B) + K(R84F) * C(SH) + K(R86F) * C(SOH) ) 
     &     * C(SCH4) 
      C(SCH3) = CCH3 / CTCHZERO( BCH3 / 2.0 + SQRT( BCH3 * BCH3 
     &     / 4.0 + ACH3 * CCH3 ) )
      C(S3XCH2)    = ( ( K(R109) * C(SO2) +  
     &    ( K(R111) + K(R112) ) * C(SO) ) * C(SC2H2)
     &    + K(R44B) * C(SCH3) * C(SH) ) /
     &    CTCHZERO( K(R39) * C(SCH3) + ( K(R40) + K(R41) ) * C(SO2) + 
     &    K(R35F) * C(SH) 
     &     + K(R44F) * ( K(R42B) / K(R42F) ) * C(SH2) ) 
       
      C(SCH)    =     K(R35F) * C(S3XCH2) * C(SH)    / 
     &    CTCHZERO( K(R26) * C(SCO2) + K(R27F) * C(SH2O) 
     &     + K(R25) * C(SO2) ) 
          
      C(SCH3OH) = (
     &    K(R90B) * C(SOH) * C(SCH3) ) / ( CTCHZERO( 
     &    K(R90F)
     &    + K(R91) * C(SH)
     &    + K(R93) * C(SOH) ) )
      C(SCH3O) = (
     &    K(R55) * C(SCH3) * C(SOH) ) / ( CTCHZERO( 
     &    K(R62) * M(MM1)
     &    + K(R63) * C(SH)
     &    + K(R97) * M(MM0) ) )
      C(SCH2OH) = (
     &    K(R91) * C(SCH3OH) * C(SH)
     &    + K(R93) * C(SCH3OH) * C(SOH)
     &    + K(R97) * C(SCH3O) * M(MM0) ) / ( CTCHZERO( 
     &    K(R69) * M(MM1)
     &    + K(R70) * C(SH)
     &    + K(R71) * C(SO2) ) )

      C(S1XCH2) = (
     &    K(R42B) * C(S3XCH2) * M(MM1)
     &    + K(R44B) * C(SH) * C(SCH3) ) / ( CTCHZERO( 
     &    K(R42F) * M(MM1)
     &    + K(R44F) * C(SH2) ) )
C     NITROGEN CHEMISTRY
      C(SN2H) = ( K(RN38B) * C(SH) * C(SN2) ) / ( CTCHZERO(
     &    + K(RN38F)
     &    + K(RN39F) * C(SH)
     &    + K(RN40F) * C(SO)
     &    + K(RN41F) * C(SOH) ) )
      C(SCN) = ( K(RN130B) * C(SH) * C(SHCN)
     &    + K(RN131B) * C(SOH) * C(SHCN) ) / ( CTCHZERO(
     &    K(RN130F) * C(SH2)
     &    + K(RN131F) * C(SH2O)
     &    + K(RN134F) * C(SOH)
     &    + K(RN136F) * C(SO2) ) )
      C(SNCO) = ( K(RN62F) * C(SHNCO) * C(SO)
     &    + K(RN65F) * C(SHNCO) * C(SOH)
     &    + K(RN74B) * C(SH) * C(SHNCO)
     &    + K(RN120F) * C(SHCN) * C(SO)
     &    + K(RN134F) * C(SCN) * C(SOH)
     &    + K(RN136F) * C(SCN) * C(SO2) ) / ( CTCHZERO(
     &    K(RN62B) * C(SOH)
     &    + K(RN65B) * C(SH2O)
     &    + K(RN70F) * M(MM7)
     &    + K(RN71F) * C(SH)
     &    + K(RN72F) * C(SO)
     &    + K(RN74F) * C(SH2)
     &    + K(RN120B) * C(SH)
     &    + K(RN134B) * C(SH)
     &    + K(RN136B) * C(SO) ) )
      DN  =  CTCHZERO(  K(RN16B) * C(SH2)
     &    + K(RN19B) * C(SH2O)
     &    + K(RN35F) * C(SO2)
     &    + K(RN36F) * C(SOH)
     &    + K(RN37F) * C(SNO)
     &    + K(RN70B) * C(SCO) * M(MM7)
     &    + K(RN110B) * C(SHCN) )
      FN  = ( K(RN16B) * C(SH2) + 
     &    K(RN19B) * C(SH2O) ) / DN 
      DHNO  =  CTCHZERO(  K(RN7B) * C(SH)
     &    + K(RN18B) * C(SH)
     &    + K(RN20B) * C(SO)
     &    + K(RN27F) * M(MM3)
     &    + K(RN28F) * C(SH)
     &    + K(RN30F) * C(SOH) ) 
      FHNO  = ( K(RN18B) * C(SH) + 
     &    K(RN20B) * C(SO) ) / DHNO 
      GHNO  = ( K(RN7B) * C(SH) ) / DHNO 
       
      DNH =  CTCHZERO( 
     &    K(RN6B) * C(SH2)
     &    + K(RN9B) * C(SH2O)
     &    + K(RN16F) * C(SH) * ( 1 - FN )
     &    + K(RN17F) * C(SO) 
     &    + K(RN18F) * C(SOH) * ( 1 - FHNO )
     &    + K(RN19F) * C(SOH) * ( 1 - FN )
     &    + K(RN20F) * C(SO2) * ( 1 - FHNO )
     &    + K(RN24XF) * C(SNO)
     &    + K(RN24Y) * C(SNO)
     &    + K(RN25F) * C(SNO)
     &    + K(RN60B) * C(SCO) * M(MM7)
     &    + K(RN63B) * C(SCO2)
     &    + K(RN71B) * C(SCO) )
      GNH =  ( K(RN6B) * C(SH2)
     &    + K(RN9B) * C(SH2O)
     &    + K(RN18F) * C(SOH) * GHNO
     &    + K(RN20F) * C(SO2) * GHNO ) / DNH 
      DNH2 =  CTCHZERO( 
     &    K(RN1B) * C(SH) * M(MM6)
     &    + K(RN2B) * C(SH2)
     &    + K(RN3B) * C(SOH)
     &    + K(RN4B) * C(SH2O)
     &    + K(RN6F) * C(SH) * ( 1 - GNH )
     &    + K(RN7F) * C(SO) * ( 1 - GHNO - GNH * FHNO )
     &    + K(RN9F) * C(SOH) * ( 1 - GNH )
     &    + K(RN13) * C(SNO)
     &    + K(RN14F) * C(SNO)
     &    + K(RN61B) * C(SCO) )
      DNH3 = CTCHZERO( K(RN1F) * M(MM6)
     &    + K(RN2F) * C(SH)
     &    + K(RN3F) * C(SO)
     &    + K(RN4F) * C(SOH) ) 
      C(SNH2) = (
     &    K(RN1F) * C(SNH3) * M(MM6) + K(RN2F) * C(SNH3) * C(SH)
     &    + K(RN3F) * C(SNH3) * C(SO) + K(RN4F) * C(SNH3) * C(SOH)
     &    + K(RN14B) * C(SOH) * C(SN2H) + K(RN61F) * C(SHNCO) * C(SH)
     &    + K(RN60F) * C(SHNCO) * M(MM7) * GNH
     &    + K(RN63F) * C(SHNCO) * C(SO) * GNH
     &    + K(RN17B) * C(SH) * C(SNO) * GNH
     &    + K(RN24XB) * C(SH) * C(SN2O) * GNH
     &    + K(RN25B) * C(SOH) * C(SN2) * GNH
     &    + K(RN71F) * C(SNCO) * C(SH) * GNH 
     &    +  K(RN35B) * C(SO) * C(SNO) * GNH * FN
     &    + K(RN36B) * C(SH) * C(SNO) * GNH * FN
     &    + K(RN37B) * C(SO) * C(SN2) * GNH * FN
     &    + K(RN70F) * C(SNCO) * M(MM7) * GNH * FN
     &    + K(RN110F) * C(SN2) * C(SCH) * GNH * FN
     &    +  K(RN27B) * C(SNO) * C(SH) * M(MM3) * ( GHNO + GNH * FHNO )
     &    + K(RN28B) * C(SH2) * C(SNO) * ( GHNO + GNH * FHNO )
     &    + K(RN30B) * C(SH2O) * C(SNO) * ( GHNO + GNH * FHNO )
     &    ) / DNH2 
       
      C(SNH) = (
     &    K(RN6F) * C(SNH2) * C(SH)
     &    + K(RN9F) * C(SNH2) * C(SOH)
     &    + K(RN17B) * C(SH) * C(SNO)
     &    + K(RN24XB) * C(SH) * C(SN2O)
     &    + K(RN25B) * C(SOH) * C(SN2)
     &    + K(RN60F) * C(SHNCO) * M(MM7)
     &    + K(RN63F) * C(SHNCO) * C(SO)
     &    + K(RN71F) * C(SNCO) * C(SH)  
     &    + K(RN35B) * C(SO) * C(SNO) * FN
     &    + K(RN36B) * C(SH) * C(SNO) * FN
     &    + K(RN37B) * C(SO) * C(SN2) * FN
     &    + K(RN70F) * C(SNCO) * M(MM7) * FN
     &    + K(RN110F) * C(SN2) * C(SCH) * FN
     &    +  K(RN27B) * C(SNO) * C(SH) * M(MM3) * FHNO
     &    + K(RN28B) * C(SH2) * C(SNO) * FHNO
     &    + K(RN30B) * C(SH2O) * C(SNO) * FHNO
     &    ) / DNH 
      C(SN) = ( K(RN16F) * C(SNH) * C(SH)
     &    + K(RN19F) * C(SNH) * C(SOH)
     &    + K(RN35B) * C(SO) * C(SNO)
     &    + K(RN36B) * C(SH) * C(SNO)
     &    + K(RN37B) * C(SO) * C(SN2)
     &    + K(RN70F) * C(SNCO) * M(MM7)
     &    + K(RN110F) * C(SN2) * C(SCH) ) / DN 
      C(SHNO) = ( K(RN7F) * C(SNH2) * C(SO)
     &    + K(RN18F) * C(SNH) * C(SOH)
     &    + K(RN20F) * C(SNH) * C(SO2)
     &    + K(RN27B) * C(SNO) * C(SH) * M(MM3)
     &    + K(RN28B) * C(SH2) * C(SNO)
     &    + K(RN30B) * C(SH2O) * C(SNO) ) / DHNO 


    
