C     Eleven-step mechanism submitted to the 26th Symposium
C     version from March 1996
C     Citation:
C     J. Hewson and M. Bollig, 
C     Reduced Mechanisms for NOx Emissions from Hydrocarbon Diffusion Flames, 
C     Twenty-Sisxth Symposium (International) on Combustion,
C     1996
C
C     Global Reactions
C     I: 	CH4 + 2 H + H2O -> 4H2 + CO		
C     II:	H2O + CO -> H2 + CO2			
C     III: 	2 H -> H2				
C     IV:  	O2 + 3 H2 -> 2 H2O + 2 H 		
C     V:	H2 + 2 CO -> C2H2 + O2 			
C     VI: 	N2 + O2 -> 2 NO				
C     VII:	3 H2 + CO + NO -> HCN + 2 H2O + H	
C     VIII:	H + H2O + CO + NO -> HNCO + O2 + H2 	
C     IX:	NH3 + H + H2O -> NO + 3 H2		
C     X:	N2 + H2O -> N2O + H2			
C     XI:	NO + 2 H2 + O2 -> NO2 + H2O + 2 H	
C
C-----------------------------------------------------------------------
CELEMENTS  H   O C  N END
CSPECIES CH4 H H2O H2 CO CO2 O2 C2H2 N2 NO HCN HNCO NH3 N2O NO2 END
CREACTIONS END
C ======= 11stepF.ckwyp =======
C-----------------------------------------------------------------------
      SUBROUTINE CKWYP_EXT (P, T, Y, WDOT, GAS) 
C  INPUT
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C                   Data type - real scalar
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     Y      - Mass fractions of the species.
C                   cgs units - none
C                   Data type - real array
C                   Dimension Y(*) at least KK, the total number of
C                   species.
C
C  OUTPUT
C     WDOT   - Chemical molar production rates of the species.
C                   cgs units - moles/(cm**3*sec)
C                   Data type - real array
C                   Dimension WDOT(*) at least KK, the total number of
C                   species.
C
C  END PROLOGUE
C-----------------------------------------------------------------------
C ======= 11stepF.f =======
C-----------------------------------------------------------------------
C      SUBROUTINE PRODRATES( CDOT, W, K, C, M, TEMP,
C     &    PRESSURE )
C-----------------------------------------------------------------------
C     THIS SUBROUTINE COMPUTES RATES OF PRODUCTION CDOT
C     IN [KMOLE/(M^3S)]. THE PARAMETERS W ( REACTION RATE ),
C     K ( RATE COEFFICIENT ) AND M ( THIRD BODY CONCENTRATIONS ) ARE
C     JUST WORK SPACE FOR THIS FUNCTION.
C     C CONTAINS THE CONCENTRATIONS OF NON STEADY STATE SPECIES IN
C     [KMOLE/M^3] AND IS WORKSPACE FOR THE STEADY STATE CONCENTRATIONS,
C     WHICH ARE COMPUTED IN THIS FUNCTION.
C     TEMP IS THE TEMPERATURE IN [K] AND
C     PRESSURE_LOC IS THE PRESSURE_LOC IN [PA].
C     CALLED FUNCTIONS ARE 'GETLINDRATECOEFF', 'COMPSTEADYSTATES',
C     'CTCHZERO'
C-----------------------------------------------------------------------
      USE cantera
      IMPLICIT NONE
      type(phase_t) gas
      DOUBLE PRECISION CDOT(15), W(209), K(209), C(40), M(6)
      DOUBLE PRECISION TEMP
      DOUBLE PRECISION PRESSURE_LOC
      INTEGER I
      DOUBLE PRECISION GETLINDRATECOEFF, LT
      DOUBLE PRECISION CTCHZERO
      DOUBLE PRECISION RGAS, RT
      PARAMETER ( RGAS = 8314.34 )
      DOUBLE PRECISION K34F0, K34FINF
      DOUBLE PRECISION FC34
      DOUBLE PRECISION K34B0, K34BINF
      DOUBLE PRECISION K36F0, K36FINF
      DOUBLE PRECISION FC36
      DOUBLE PRECISION K36B0, K36BINF
      DOUBLE PRECISION K51F0, K51FINF
      DOUBLE PRECISION FC51
      DOUBLE PRECISION K51B0, K51BINF
      DOUBLE PRECISION K58F0, K58FINF
      DOUBLE PRECISION FC58
      DOUBLE PRECISION K58B0, K58BINF
      DOUBLE PRECISION KN45F0, KN45FINF
      DOUBLE PRECISION FCN45
      DOUBLE PRECISION KN45B0, KN45BINF
      DOUBLE PRECISION REV36
C- chemkin stuff
      DOUBLE PRECISION WT(15),Y(15)
      DOUBLE PRECISION RU,SUMYOW
      DOUBLE PRECISION WDOT(15)
      DOUBLE PRECISION T, P
C-end of chemkin stuff
      INCLUDE '11stepF.h'
C
C -conversion of stuff from chemkin 
      DATA RU/8.314510D7/
      DATA WT/ 16.043 , 1.008 ,18.015 ,2.016 ,28.011 ,44.010,
     & 31.999, 26.038, 28.013, 30.006, 27.026 ,43.025, 17.0305, 
     & 44.0125, 46.005/
C SPECIES CH4 H H2O H2 CO CO2 O2 C2H2 N2 NO HCN HNCO NH3 N2O NO2 END
      WRITE(*,*) "calling ext!!!"
      TEMP= T
      PRESSURE_LOC=P* 1.E-1
      SUMYOW = 0.0
      DO 10 I = 1,  15
        SUMYOW = SUMYOW + Y(I)/WT(I)
 10   CONTINUE 
      SUMYOW = SUMYOW*T*RU
      DO 20  I = 1,  15
        C(I) = P*Y(I) / ( SUMYOW*WT(I) ) *1.E+3  
 20   CONTINUE
C- end of chemin conversion stuff
      LT = DLOG( TEMP )
      RT = RGAS * TEMP 
C-----------------------------------------------------------------------
C     COMPUTE THE REACTION RATE CONSTANTS
C-----------------------------------------------------------------------
      DO 100 I=1,209
        K(I)=A(I)*DEXP(N(I)*LT-E(I)/RT)
 100  CONTINUE
C
C     PRESSURE DEPENDENT REACTIONS
C
      K34F0 = 6.2570000000D+17 * DEXP(-1.8 * LT)
      K34FINF = 2.1080000000D+11
      FC34 = 0.577000 * DEXP ( -TEMP / 2370.00 )
      K(R34F) = GETLINDRATECOEFF( TEMP, PRESSURE_LOC, K34F0, K34FINF
     &    , FC34 )
      K34B0 = 6.5851871268D+22 * DEXP(-1.8 * LT - 439072242.3 / RT)
      K34BINF = 2.2185671189D+16 * DEXP(-439072242.3 / RT)
      K(R34B) = GETLINDRATECOEFF( TEMP, PRESSURE_LOC, K34B0, K34BINF
     &    , FC34 )
      K36F0 = 1.2720000000D+35 * DEXP(-7 * LT - 11560000 / RT)
      K36FINF = 1.8130000000D+10
      FC36 = 0.380000 * DEXP ( -TEMP / 73.0000 ) 
     &     + 0.620000 * DEXP ( -TEMP / 1180.00 )
      K(R36F) = GETLINDRATECOEFF( TEMP, PRESSURE_LOC, K36F0, K36FINF
     &    , FC36 )
      K36B0 = 5.0084110374D+41 * DEXP(-7 * LT - 384780198.3 / RT)
      K36BINF = 7.1385607004D+16 * DEXP(-373220198.3 / RT)
      K(R36B) = GETLINDRATECOEFF( TEMP, PRESSURE_LOC, K36B0, K36BINF
     &    , FC36 )
      K51F0 = 1.1870000000D+39 * DEXP(-7.5 * LT - 190400000 / RT)
      K51FINF = 2.0000000000D+14 * DEXP(-166290000 / RT)
      FC51 = 0.35
      K(R51F) = GETLINDRATECOEFF( TEMP, PRESSURE_LOC, K51F0, K51FINF
     &    , FC51 )
      K51B0 = 6.2450000000D+35 * DEXP(-7.5 * LT - 27500000 / RT)
      K51BINF = 1.0530000000D+11 * DEXP(-3390000 / RT)
      K(R51B) = GETLINDRATECOEFF( TEMP, PRESSURE_LOC, K51B0, K51BINF
     &    , FC51 )
      K58F0 = 1.0000000000D+13 * DEXP(-126000000 / RT)
      K58FINF = 1.3000000000D+13 * DEXP(-167000000 / RT)
      FC58 = 0.411 * DEXP( -73.4 / TEMP ) + DEXP( -TEMP / 422.8 )
      K(R58F) = GETLINDRATECOEFF( TEMP, PRESSURE_LOC, K58F0, K58FINF
     &    , FC58 )
      K58B0 = 1.5950000000D+10 * DEXP(27390000 / RT)
      K58BINF = 2.0730000000D+10 * DEXP(-13610000 / RT)
      K(R58B) = GETLINDRATECOEFF( TEMP, PRESSURE_LOC, K58B0, K58BINF
     &    , FC58 )
      KN45F0 = 3.9800000000D+11 * DEXP(-237000000 / RT)
      KN45FINF = 1.6000000000D+12 * DEXP(-262000000 / RT)
      FCN45 =  1.0
      K(RN45F) = GETLINDRATECOEFF( TEMP, PRESSURE_LOC, KN45F0, KN45FINF
     &    , FCN45 )
      KN45B0 = 4.7194802843D+06 * DEXP(-73764195.3 / RT)
      KN45BINF = 1.8972785063D+07 * DEXP(-98764195.3 / RT)
      K(RN45B) = GETLINDRATECOEFF( TEMP, PRESSURE_LOC, KN45B0, KN45BINF
     &    , FCN45 )
C
C-----------------------------------------------------------------------
C     COMPUTE THE THIRD-BODY CONCENTRATIONS
C-----------------------------------------------------------------------

      M(MM1) = 3 * C(SCH4) + 6.5 * C(SH2O)
     &     + C(SH2) + 0.75 * C(SCO)
     &     + 0.4 * C(SO2) + 0.4 * C(SN2)
      M(MM6) = C(SCH4) + C(SH)
     &     + C(SH2O) + C(SH2)
     &     + C(SCO) + C(SCO2)
     &     + C(SO2) + C(SC2H2)
     &     + C(SN2) + C(SNO)
     &     + C(SHCN) + C(SHNCO)
     &     + C(SNH3) + C(SN2O)
      M(MM3) = C(SCH4) + C(SH)
     &     + 10 * C(SH2O) + 2 * C(SH2)
     &     + C(SCO) + C(SCO2)
     &     + 2 * C(SO2) + C(SC2H2)
     &     + 2 * C(SN2) + C(SNO)
     &     + C(SHCN) + C(SHNCO)
     &     + C(SNH3) + C(SN2O)
      M(MM8) = C(SCH4) + C(SH)
     &     + C(SH2O) + C(SH2)
     &     + C(SCO) + C(SCO2)
     &     + C(SO2) + C(SC2H2)
     &     + C(SN2) + C(SNO)
     &     + C(SHCN) + C(SHNCO)
     &     + C(SNH3) + C(SN2O)
      M(MM7) = C(SCH4) + C(SH)
     &     + 18.6 * C(SH2O) + C(SH2)
     &     + C(SCO) + C(SCO2)
     &     + 1.5 * C(SO2) + C(SC2H2)
     &     + 1.5 * C(SN2) + C(SNO)
     &     + C(SHCN) + C(SHNCO)
     &     + C(SNH3) + C(SN2O)
C
C-----------------------------------------------------------------------
C     COMPUTE THE CONCENTRATIONS OF STEADY-STATE SPECIES
C-----------------------------------------------------------------------
      CALL COMPSTEADYSTATES( K, C, M )
C
C-----------------------------------------------------------------------
C     COMPUTE THE REACTION RATES
C-----------------------------------------------------------------------
      W(R1F) = K(R1F) * C(SO2) * C(SH)
      W(R1B) = K(R1B) * C(SO) * C(SOH)
      W(R2F) = K(R2F) * C(SH2) * C(SO)
      W(R2B) = K(R2B) * C(SH) * C(SOH)
      W(R3F) = K(R3F) * C(SH2) * C(SOH)
      W(R3B) = K(R3B) * C(SH) * C(SH2O)
      W(R4F) = K(R4F) * C(SOH) * C(SOH)
      W(R4B) = K(R4B) * C(SO) * C(SH2O)
      W(R5F) = K(R5F) * C(SH) * C(SH) * M(MM1)
      W(R5B) = K(R5B) * C(SH2) * M(MM1)
      W(R6F) = K(R6F) * C(SO) * C(SO) * M(MM1)
      W(R6B) = K(R6B) * C(SO2) * M(MM1)
      W(R7F) = K(R7F) * C(SH) * C(SOH) * M(MM1)
      W(R7B) = K(R7B) * C(SH2O) * M(MM1)
      W(R8F) = K(R8F) * C(SH) * C(SO2) * M(MM1)
      W(R8B) = K(R8B) * C(SHO2) * M(MM1)
      W(R9) = K(R9) * C(SHO2) * C(SH)
      W(R10) = K(R10) * C(SHO2) * C(SH)
      W(R11) = K(R11) * C(SHO2) * C(SH)
      W(R12) = K(R12) * C(SHO2) * C(SO)
      W(R13F) = K(R13F) * C(SHO2) * C(SOH)
      W(R13B) = K(R13B) * C(SO2) * C(SH2O)
      W(R20F) = K(R20F) * C(SCO) * C(SOH)
      W(R20B) = K(R20B) * C(SH) * C(SCO2)
      W(R21) = K(R21) * C(SCO) * C(SHO2)
      W(R22) = K(R22) * C(SCO) * C(SO) * M(MM1)
      W(R25) = K(R25) * C(SCH) * C(SO2)
      W(R26) = K(R26) * C(SCH) * C(SCO2)
      W(R27F) = K(R27F) * C(SCH) * C(SH2O)
      W(R27B) = K(R27B) * C(SCH2OH)
      W(R28F) = K(R28F) * C(SCHO) * M(MM1)
      W(R28B) = K(R28B) * C(SH) * C(SCO) * M(MM1)
      W(R35F) = K(R35F) * C(S3XCH2) * C(SH)
      W(R35B) = K(R35B) * C(SH2) * C(SCH)
      W(R38) = K(R38) * C(S3XCH2) * C(S3XCH2)
      W(R39) = K(R39) * C(S3XCH2) * C(SCH3)
      W(R40) = K(R40) * C(S3XCH2) * C(SO2)
      W(R41) = K(R41) * C(S3XCH2) * C(SO2)
      W(R42F) = K(R42F) * C(S1XCH2) * M(MM1)
      W(R42B) = K(R42B) * C(S3XCH2) * M(MM1)
      W(R43) = K(R43) * C(S1XCH2) * C(SO2)
      W(R44F) = K(R44F) * C(S1XCH2) * C(SH2)
      W(R44B) = K(R44B) * C(SH) * C(SCH3)
      W(R45) = K(R45) * C(SCH2O) * M(MM1)
      W(R46) = K(R46) * C(SCH2O) * C(SH)
      W(R48) = K(R48) * C(SCH2O) * C(SOH)
      W(R53) = K(R53) * C(SCH3) * C(SO)
      W(R34F) = K(R34F) * C(SCH3) * C(SH)
      W(R34B) = K(R34B) * C(SCH4)
      W(R55) = K(R55) * C(SCH3) * C(SOH)
      W(R60) = K(R60) * C(SCH3) * C(SCH3)
      W(R36F) = K(R36F) * C(SCH3) * C(SCH3)
      W(R36B) = K(R36B) * C(SC2H6)
      W(R62) = K(R62) * C(SCH3O) * M(MM1)
      W(R63) = K(R63) * C(SCH3O) * C(SH)
      W(R69) = K(R69) * C(SCH2OH) * M(MM1)
      W(R70) = K(R70) * C(SCH2OH) * C(SH)
      W(R71) = K(R71) * C(SCH2OH) * C(SO2)
      W(R84F) = K(R84F) * C(SCH4) * C(SH)
      W(R84B) = K(R84B) * C(SCH3) * C(SH2)
      W(R86F) = K(R86F) * C(SCH4) * C(SOH)
      W(R86B) = K(R86B) * C(SCH3) * C(SH2O)
      W(R90F) = K(R90F) * C(SCH3OH)
      W(R90B) = K(R90B) * C(SOH) * C(SCH3)
      W(R91) = K(R91) * C(SCH3OH) * C(SH)
      W(R93) = K(R93) * C(SCH3OH) * C(SOH)
      W(R97) = K(R97) * C(SCH3O) * M(MM0)
      W(R105F) = K(R105F) * C(SHCCO) * C(SH)
      W(R105B) = K(R105B) * C(SCO) * C(S3XCH2)
      W(R106) = K(R106) * C(SHCCO) * C(SO)
      W(R109) = K(R109) * C(SC2H2) * C(SO2)
      W(R111) = K(R111) * C(SC2H2) * C(SO)
      W(R112) = K(R112) * C(SC2H2) * C(SO)
      W(R51F) = K(R51F) * C(SC2H3)
      W(R51B) = K(R51B) * C(SC2H2) * C(SH)
      W(RA125) = K(RA125) * C(SC2H3) * C(SO2)
      W(R129F) = K(R129F) * C(SC2H4) * M(MM1)
      W(R129B) = K(R129B) * C(SH2) * C(SC2H2) * M(MM1)
      W(R131F) = K(R131F) * C(SC2H4) * C(SH)
      W(R131B) = K(R131B) * C(SH2) * C(SC2H3)
      W(R134F) = K(R134F) * C(SC2H4) * C(SOH)
      W(R134B) = K(R134B) * C(SH2O) * C(SC2H3)
      W(R58F) = K(R58F) * C(SC2H5)
      W(R58B) = K(R58B) * C(SC2H4) * C(SH)
      W(R146) = K(R146) * C(SC2H5) * C(SH)
      W(R164) = K(R164) * C(SC2H6) * C(SH)
      W(R166) = K(R166) * C(SC2H6) * C(SOH)
      W(R170F) = K(R170F) * C(SC2H6) * C(SCH3)
      W(R170B) = K(R170B) * C(SCH4) * C(SC2H5)
      W(RN1F) = K(RN1F) * C(SNH3) * M(MM6)
      W(RN1B) = K(RN1B) * C(SH) * C(SNH2) * M(MM6)
      W(RN2F) = K(RN2F) * C(SNH3) * C(SH)
      W(RN2B) = K(RN2B) * C(SH2) * C(SNH2)
      W(RN3F) = K(RN3F) * C(SNH3) * C(SO)
      W(RN3B) = K(RN3B) * C(SOH) * C(SNH2)
      W(RN4F) = K(RN4F) * C(SNH3) * C(SOH)
      W(RN4B) = K(RN4B) * C(SH2O) * C(SNH2)
      W(RN6F) = K(RN6F) * C(SNH2) * C(SH)
      W(RN6B) = K(RN6B) * C(SH2) * C(SNH)
      W(RN7F) = K(RN7F) * C(SNH2) * C(SO)
      W(RN7B) = K(RN7B) * C(SH) * C(SHNO)
      W(RN9F) = K(RN9F) * C(SNH2) * C(SOH)
      W(RN9B) = K(RN9B) * C(SH2O) * C(SNH)
      W(RN101F) = K(RN101F) * C(SH2NO) * C(SO)
      W(RN101B) = K(RN101B) * C(SO2) * C(SNH2)
      W(RN13) = K(RN13) * C(SNH2) * C(SNO)
      W(RN14F) = K(RN14F) * C(SNH2) * C(SNO)
      W(RN14B) = K(RN14B) * C(SOH) * C(SN2H)
      W(RN16F) = K(RN16F) * C(SNH) * C(SH)
      W(RN16B) = K(RN16B) * C(SH2) * C(SN)
      W(RN17F) = K(RN17F) * C(SNH) * C(SO)
      W(RN17B) = K(RN17B) * C(SH) * C(SNO)
      W(RN18F) = K(RN18F) * C(SNH) * C(SOH)
      W(RN18B) = K(RN18B) * C(SH) * C(SHNO)
      W(RN19F) = K(RN19F) * C(SNH) * C(SOH)
      W(RN19B) = K(RN19B) * C(SH2O) * C(SN)
      W(RN20F) = K(RN20F) * C(SNH) * C(SO2)
      W(RN20B) = K(RN20B) * C(SO) * C(SHNO)
      W(RN24XF) = K(RN24XF) * C(SNH) * C(SNO)
      W(RN24XB) = K(RN24XB) * C(SH) * C(SN2O)
      W(RN24Y) = K(RN24Y) * C(SNH) * C(SNO)
      W(RN25F) = K(RN25F) * C(SNH) * C(SNO)
      W(RN25B) = K(RN25B) * C(SOH) * C(SN2)
      W(RN27F) = K(RN27F) * C(SHNO) * M(MM3)
      W(RN27B) = K(RN27B) * C(SNO) * C(SH) * M(MM3)
      W(RN28F) = K(RN28F) * C(SHNO) * C(SH)
      W(RN28B) = K(RN28B) * C(SH2) * C(SNO)
      W(RN30F) = K(RN30F) * C(SHNO) * C(SOH)
      W(RN30B) = K(RN30B) * C(SH2O) * C(SNO)
      W(RN102F) = K(RN102F) * C(SH2NO) * M(MM6)
      W(RN102B) = K(RN102B) * C(SH) * C(SHNO) * M(MM6)
      W(RN104F) = K(RN104F) * C(SH2NO) * C(SH)
      W(RN104B) = K(RN104B) * C(SOH) * C(SNH2)
      W(RN35F) = K(RN35F) * C(SN) * C(SO2)
      W(RN35B) = K(RN35B) * C(SO) * C(SNO)
      W(RN36F) = K(RN36F) * C(SN) * C(SOH)
      W(RN36B) = K(RN36B) * C(SH) * C(SNO)
      W(RN37F) = K(RN37F) * C(SN) * C(SNO)
      W(RN37B) = K(RN37B) * C(SO) * C(SN2)
      W(RN38F) = K(RN38F) * C(SN2H)
      W(RN38B) = K(RN38B) * C(SH) * C(SN2)
      W(RN39F) = K(RN39F) * C(SN2H) * C(SH)
      W(RN39B) = K(RN39B) * C(SH2) * C(SN2)
      W(RN40F) = K(RN40F) * C(SN2H) * C(SO)
      W(RN40B) = K(RN40B) * C(SH) * C(SN2O)
      W(RN41F) = K(RN41F) * C(SN2H) * C(SOH)
      W(RN41B) = K(RN41B) * C(SH2O) * C(SN2)
      W(RN45F) = K(RN45F) * C(SN2O)
      W(RN45B) = K(RN45B) * C(SO) * C(SN2)
      W(RN46YF) = K(RN46YF) * C(SN2O) * C(SH)
      W(RN46YB) = K(RN46YB) * C(SOH) * C(SN2)
      W(RN47F) = K(RN47F) * C(SN2O) * C(SO)
      W(RN47B) = K(RN47B) * C(SNO) * C(SNO)
      W(RN84F) = K(RN84F) * C(SNO2) * M(MM8)
      W(RN84B) = K(RN84B) * C(SO) * C(SNO) * M(MM8)
      W(RN85F) = K(RN85F) * C(SNO) * C(SHO2)
      W(RN85B) = K(RN85B) * C(SOH) * C(SNO2)
      W(RN87F) = K(RN87F) * C(SNO2) * C(SH)
      W(RN87B) = K(RN87B) * C(SOH) * C(SNO)
      W(RN88F) = K(RN88F) * C(SNO2) * C(SO)
      W(RN88B) = K(RN88B) * C(SO2) * C(SNO)
      W(RN60F) = K(RN60F) * C(SHNCO) * M(MM7)
      W(RN60B) = K(RN60B) * C(SCO) * C(SNH) * M(MM7)
      W(RN61F) = K(RN61F) * C(SHNCO) * C(SH)
      W(RN61B) = K(RN61B) * C(SCO) * C(SNH2)
      W(RN62F) = K(RN62F) * C(SHNCO) * C(SO)
      W(RN62B) = K(RN62B) * C(SOH) * C(SNCO)
      W(RN63F) = K(RN63F) * C(SHNCO) * C(SO)
      W(RN63B) = K(RN63B) * C(SCO2) * C(SNH)
      W(RN65F) = K(RN65F) * C(SHNCO) * C(SOH)
      W(RN65B) = K(RN65B) * C(SH2O) * C(SNCO)
      W(RN70F) = K(RN70F) * C(SNCO) * M(MM7)
      W(RN70B) = K(RN70B) * C(SCO) * C(SN) * M(MM7)
      W(RN71F) = K(RN71F) * C(SNCO) * C(SH)
      W(RN71B) = K(RN71B) * C(SNH) * C(SCO)
      W(RN72F) = K(RN72F) * C(SNCO) * C(SO)
      W(RN72B) = K(RN72B) * C(SCO) * C(SNO)
      W(RN74F) = K(RN74F) * C(SNCO) * C(SH2)
      W(RN74B) = K(RN74B) * C(SH) * C(SHNCO)
      W(RN110F) = K(RN110F) * C(SN2) * C(SCH)
      W(RN110B) = K(RN110B) * C(SN) * C(SHCN)
      W(RN111) = K(RN111) * C(SNO) * C(SCH3)
      W(RN112F) = K(RN112F) * C(SNO) * C(S3XCH2)
      W(RN112B) = K(RN112B) * C(SH) * C(SHNCO)
      W(RN113F) = K(RN113F) * C(SNO) * C(SCH)
      W(RN113B) = K(RN113B) * C(SO) * C(SHCN)
      W(RN120F) = K(RN120F) * C(SHCN) * C(SO)
      W(RN120B) = K(RN120B) * C(SH) * C(SNCO)
      W(RN130F) = K(RN130F) * C(SCN) * C(SH2)
      W(RN130B) = K(RN130B) * C(SH) * C(SHCN)
      W(RN131F) = K(RN131F) * C(SCN) * C(SH2O)
      W(RN131B) = K(RN131B) * C(SOH) * C(SHCN)
      W(RN134F) = K(RN134F) * C(SCN) * C(SOH)
      W(RN134B) = K(RN134B) * C(SH) * C(SNCO)
      W(RN136F) = K(RN136F) * C(SCN) * C(SO2)
      W(RN136B) = K(RN136B) * C(SO) * C(SNCO)
C
C     GLOBAL REACTION RATES
C
      REV36 = ( K(R164) * C(SH) + K(R166) * C(SOH)  ) / 
     &    CTCHZERO( K(R36B) + K(R164) * C(SH) + K(R166) * C(SOH) ) 
      K(RI) =     W(R84F) + W(R86F) + W(R34B) - W(R84B) 
     &          - W(R86B) - W(R34F) 
      K(RII) =     W(R20F) - W(R20B) 
      K(RIII) =     W(R8F) - W(R8B) + W(R41) + W(R34F) - W(R34B) 
     &     - W(R55) + W(R63) + W(R70) + W(R71) 
      K(RIIIN) =     - W(RN1F) + W(RN1B) + W(RN13) + W(RN16F) 
     &     + W(RN17F) - W(RN17B) + W(RN19F) - W(RN19B) + W(RN24XF) 
     &     - W(RN24XB) + W(RN24Y) + W(RN25F) - W(RN25B) + W(RN28F) 
     &     - W(RN28B) + W(RN30F) - W(RN30B)  +W(RN39F) - W(RN39B) 
     &     + W(RN40F) -W(RN40B) + W(RN41F) - W(RN41B) - W(RN45F) 
     &     + W(RN45B) - W(RN60F) + W(RN60B) + W(RN72F) - W(RN72B)  
      K(RIV) =     W(R1F) - W(R1B) + W(R36F)*REV36 + W(R60)  
     &     + W(R25) - W(R35F) + W(R35B) - W(R44F) + W(R44B) 
      K(RIVN) =     + W(RN13) + W(RN14F) - W(RN14B) + W(RN20F) 
     &     - W(RN20B) + W(RN24XF) - W(RN24XB) +  W(RN24Y) + W(RN25F) 
     &     - W(RN25B) + W(RN35F) - W(RN35B) + W(RN37F) - W(RN37B) 
     &     - W(RN47F) - W(RN60F) + W(RN60B) - W(RN61F) + W(RN61B) 
     &     - W(RN62F) + W(RN62B) - W(RN63F) + W(RN63B) - W(RN65F) 
     &     + W(RN65B) + W(RN74F) - W(RN74B) - W(RN110F)  + W(RN136F) 
     &     - W(RN136B) 
      K(RV) =      W(R36F)*REV36 + W(R60) - W(R109) 
     &     - W(R111) - W(R112) 
      K(RVI) =      W(RN37B) + W(RN110F) + W(RN47F)  + W(RN24XB) 
     &      - W(RN13) - W(RN14F) + W(RN14B) - W(RN24XF) 
     &     - W(RN24Y) - W(RN25F) - W(RN37F)  
      K(RVII) =     W(RN110F) + W(RN111) + W(RN113F)  
     &     + W(RN120B)  + W(RN134B) + W(RN136B) 
     &     - W(RN113B) - W(RN120F)  - W(RN134F)  - W(RN136F) 
      K(RVIII) =     + W(RN60B)  + W(RN61B) + W(RN62B)  
     &     + W(RN63B)  + W(RN65B) + W(RN74F) + W(RN112F)   
     &     - W(RN60F)  - W(RN61F) - W(RN62F) 
     &     - W(RN63F) - W(RN65F) - W(RN74B) - W(RN112B)   
      K(RIX) = W(RN1F) + W(RN2F) + W(RN3F) + W(RN4F) 
     &     - W(RN1B) - W(RN2B) - W(RN3B) - W(RN4B)  
      K(RX) = + W(RN24XF)  + W(RN24Y) + W(RN40F) + W(RN45B) 
     &     + W(RN46YB) + W(RN47B) - W(RN24XB) - W(RN40B) 
     &     - W(RN45F) - W(RN46YF) - W(RN47F) 
      K(RXI) =       - W(RN84F) + W(RN84B) + W(RN85F) - W(RN85B) 
     &     - W(RN87F) + W(RN87B) - W(RN88F) + W(RN88B)       
      W(RI) = K(RI) 
      W(RII) = K(RII) 
      W(RIII) = K(RIII) 
      W(RIIIN) = K(RIIIN) 
      W(RIV) = K(RIV) 
      W(RIVN) = K(RIVN) 
      W(RV) = K(RV) 
      W(RVI) = K(RVI) 
      W(RVII) = K(RVII) 
      W(RVIII) = K(RVIII) 
      W(RIX) = K(RIX) 
      W(RX) = K(RX) 
      W(RXI) = K(RXI) 

C-----------------------------------------------------------------------
C     COMPUTE THE SPECIES SOURCE TERMS 
C-----------------------------------------------------------------------

      CDOT(SCH4) = - W(RI) 
      CDOT(SH) = - 2 * W(RI) - 2 * W(RIII) - 2 * W(RIIIN)
     &     + 2 * W(RIV) + 2 * W(RIVN) + W(RVII)
     &     - W(RVIII) - W(RIX) + 2 * W(RXI)
      CDOT(SH2O) = - W(RI) - W(RII) + 2 * W(RIV)
     &     + 2 * W(RIVN) + 2 * W(RVII) - W(RVIII)
     &     - W(RIX) - W(RX) + W(RXI)
      CDOT(SH2) = 4 * W(RI) + W(RII) + W(RIII)
     &     + W(RIIIN) - 3 * W(RIV) - 3 * W(RIVN)
     &     - W(RV) - 3 * W(RVII) + W(RVIII)
     &     + 3 * W(RIX) + W(RX) - 2 * W(RXI)
      CDOT(SCO) = W(RI) - W(RII) - 2 * W(RV)
     &     - W(RVII) - W(RVIII) 
      CDOT(SCO2) = W(RII) 
      CDOT(SO2) = - W(RIV) - W(RIVN) + W(RV)
     &     - W(RVI) + W(RVIII) - W(RXI)
      CDOT(SC2H2) = W(RV) 
      CDOT(SN2) = - W(RVI) - W(RX) 
      CDOT(SNO) = 2 * W(RVI) - W(RVII) - W(RVIII)
     &     + W(RIX) - W(RXI) 
      CDOT(SHCN) = W(RVII) 
      CDOT(SHNCO) = W(RVIII) 
      CDOT(SNH3) = - W(RIX) 
      CDOT(SN2O) = W(RX) 
      CDOT(SNO2) = W(RXI) 
CSPECIES CH4 H H2O H2 CO CO2 O2 C2H2 N2 NO HCN HNCO NH3 N2O NO2 END
C-- transform into chemkin format
      DO 200 I=1,15
       WDOT(I)=1.E-3*CDOT(I)
 200  CONTINUE 
      END
C
C-----------------------------------------------------------------------
C     SOME SUBROUTINES
C-----------------------------------------------------------------------
C     THIS FUNCTION GENERATES TROE FORM PRESSURE DEPENDENT RATE CONSTANTS
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION GETLINDRATECOEFF( TEMP, PRESSURE,
     &    K0, KINF, FC )
C
      IMPLICIT NONE
      DOUBLE PRECISION TEMP, PRESSURE, K0, KINF, FC
      DOUBLE PRECISION R
      PARAMETER (R = 8314.34)
C     [J / kmole K]
      DOUBLE PRECISION NTMP
      DOUBLE PRECISION KL
      DOUBLE PRECISION F
      DOUBLE PRECISION CONC
C
      CONC = PRESSURE / ( R * TEMP )
      NTMP = 0.75 - 1.27 * DLOG10( FC )
      K0 = K0 * CONC / KINF
      KL = K0 / ( 1.0 + K0 )
      F = DLOG10( K0 ) / NTMP
      F = FC ** ( 1.0 / ( F * F + 1.0 ) )
      GETLINDRATECOEFF = KINF * F * KL
      END
C
C-----------------------------------------------------------------------
C     USED TO AVOID POSSIBLE "DIVISION BY ZERO" PROBLEMS IN COMPUTING
C     STEADY-STATE CONCENTRATIONS
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION CTCHZERO( X )
C
      IMPLICIT NONE
      DOUBLE PRECISION X
      IF ( X.EQ.0.0D0 ) THEN
         X = 1.0D-20
      ENDIF
      CTCHZERO = X
      END
C
C
      DOUBLE PRECISION FUNCTION SOLVQDRT( A, B, C )
C-----------------------------------------------------------------------
C	SOLVES THE QUADRATIC EQUATION A x^2 + B x - C = 0
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
C
      DOUBLE PRECISION RAD, A, B, C
C
      IF ( A.EQ.0.0D0 ) THEN
        IF ( B.EQ.0.0D0 ) THEN
          WRITE(*,*) '#WARNING: INVALID ARGUMENTS IN FUNCTION SOLVQDRT '
          SOLVQDRT = 0.0
        ELSE
          SOLVQDRT = C / B
        ENDIF
      ELSE
C
        B = B / A
        C = C / A
C
        RAD = 0.25D0 * B * B + C
        IF ( RAD.GE.0.0D0 ) THEN
          RAD = DSQRT( RAD ) + 0.5 * B
          IF ( RAD.NE.0.0D0 ) THEN
            SOLVQDRT = C / RAD
          ELSE
            SOLVQDRT = -0.5D0 * B + DSQRT( 0.25 * B * B + C )
          ENDIF
        ELSE
          SOLVQDRT = 0.0D0
        ENDIF
      ENDIF
      END
C
C
      SUBROUTINE GETMOLARMASS( MM )
C-----------------------------------------------------------------------
C	FILLS 'MM' WITH SPECIES MOLAR MASS IN KG/KMOLE
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
C
      DOUBLE PRECISION MM(40)
      INCLUDE '11stepF.h'
C
      MM(SCH4) =  1.60430000e+01
      MM(SH) =  1.00800000e+00
      MM(SH2O) =  1.80160000e+01
      MM(SH2) =  2.01600000e+00
      MM(SCO) =  2.80110000e+01
      MM(SCO2) =  4.40110000e+01
      MM(SO2) =  3.20000000e+01
      MM(SC2H2) =  2.60380000e+01
      MM(SN2) =  2.80140000e+01
      MM(SNO) =  3.00100000e+01
      MM(SHCN) =  2.70160000e+01
      MM(SHNCO) =  4.30160000e+01
      MM(SNH3) =  1.70300000e+01
      MM(SN2O) =  4.40160000e+01
      MM(SNO2) =  4.60100000e+01
      MM(SOH) =  1.70080000e+01
      MM(SO) =  1.60000000e+01
      MM(SHO2) =  3.30080000e+01
      MM(SCH) =  1.30190000e+01
      MM(SCHO) =  2.90190000e+01
      MM(SCH2OH) =  3.10240000e+01
      MM(S3XCH2) =  1.40270000e+01
      MM(SCH3) =  1.50350000e+01
      MM(SC2H4) =  2.80540000e+01
      MM(S1XCH2) =  1.40270000e+01
      MM(SCH2O) =  3.00270000e+01
      MM(SCH3O) =  3.10300000e+01
      MM(SC2H6) =  3.00700000e+01
      MM(SCH3OH) =  3.20320000e+01
      MM(SHCCO) =  4.10080000e+01
      MM(SC2H3) =  2.70460000e+01
      MM(SC2H5) =  2.90620000e+01
      MM(SNH2) =  1.60200000e+01
      MM(SNH) =  1.50200000e+01
      MM(SHNO) =  3.10200000e+01
      MM(SH2NO) =  3.20260000e+01
      MM(SN2H) =  2.90280000e+01
      MM(SN) =  1.40100000e+01
      MM(SNCO) =  4.20160000e+01
      MM(SCN) =  2.60160000e+01
C
      END
C
C
      SUBROUTINE GETNAMES( NAMES )
C-----------------------------------------------------------------------
C	FILLS 'NAMES' WITH SPECIES IDENTIFIER/KG
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
C
      CHARACTER *20 NAMES(40)
      INCLUDE '11stepF.h'
C
      NAMES(SCH4)='CH4                 '
      NAMES(SH)='H                   '
      NAMES(SH2O)='H2O                 '
      NAMES(SH2)='H2                  '
      NAMES(SCO)='CO                  '
      NAMES(SCO2)='CO2                 '
      NAMES(SO2)='O2                  '
      NAMES(SC2H2)='C2H2                '
      NAMES(SN2)='N2                  '
      NAMES(SNO)='NO                  '
      NAMES(SHCN)='HCN                 '
      NAMES(SHNCO)='HNCO                '
      NAMES(SNH3)='NH3                 '
      NAMES(SN2O)='N2O                 '
      NAMES(SNO2)='NO2                 '
      NAMES(SOH)='OH                  '
      NAMES(SO)='O                   '
      NAMES(SHO2)='HO2                 '
      NAMES(SCH)='CH                  '
      NAMES(SCHO)='CHO                 '
      NAMES(SCH2OH)='CH2OH               '
      NAMES(S3XCH2)='3-CH2               '
      NAMES(SCH3)='CH3                 '
      NAMES(SC2H4)='C2H4                '
      NAMES(S1XCH2)='1-CH2               '
      NAMES(SCH2O)='CH2O                '
      NAMES(SCH3O)='CH3O                '
      NAMES(SC2H6)='C2H6                '
      NAMES(SCH3OH)='CH3OH               '
      NAMES(SHCCO)='HCCO                '
      NAMES(SC2H3)='C2H3                '
      NAMES(SC2H5)='C2H5                '
      NAMES(SNH2)='NH2                 '
      NAMES(SNH)='NH                  '
      NAMES(SHNO)='HNO                 '
      NAMES(SH2NO)='H2NO                '
      NAMES(SN2H)='N2H                 '
      NAMES(SN)='N                   '
      NAMES(SNCO)='NCO                 '
      NAMES(SCN)='CN                  '
C
      END
C
C
      SUBROUTINE COMPSTEADYSTATES( K, C, M )
C-----------------------------------------------------------------------
C     THIS SUBROUTINE COMPUTES THE STEADY STATE CONCENTRATIONS FROM
C     THE CONCENTRATIONS OF COMPUTED SPECIES AND RATE COEFFICIENTS.
C     CONCENTRATIONS OF COMPUTED SPECIES MAY NOT BE ALTERED.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION CTCHZERO
      DOUBLE PRECISION K(209), C(40), M(6)
      INCLUDE '11stepF.h'
      INCLUDE '11stepF.ss'
      END
C
C
      SUBROUTINE COMPTHERMODATA( H, CP, T )
C-----------------------------------------------------------------------
C     THIS FUNCTION COMPUTES ENTHALPY 'H' AND HEAT CAPACITY 'CP' AS
C     FUNCTION OF TEMPERATURE T FOR ALL NON STEADY STATE SPECIES
C     IN UNITS [J/KG] and [J/KG K], RESPECTIVELY.
C     THE PARAMETER H AND CP SHOULD PROVIDE WORKSPACE OF LENGTH 14
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION H(14), CP(14), T
      INCLUDE '11stepF.h'
C
      IF (T.GT.1.0D3) THEN
      H(SCH4) =  5.18253444e+02 * (
     &      T * (  1.68347880e+00 + T * (  5.11861800e-03
     &      + T * ( -1.29170953e-06 + T * (  1.69639622e-10
     &      + T * ( -9.00684620e-15 ) ) ) ) ) -1.00807870e+04 )
      CP(SCH4) =  5.18253444e+02 * (
     &       1.68347880e+00 + T * (  1.02372360e-02 
     &      + T * ( -3.87512860e-06 + T * (  6.78558490e-10
     &      + T * ( -4.50342310e-14 ) ) ) ) )
      H(SH) =  8.24835317e+03 * (
     &      T * (  2.50000000e+00 + T * (  0.00000000e+00
     &      + T * (  0.00000000e+00 + T * (  0.00000000e+00
     &      + T * (  0.00000000e+00 ) ) ) ) ) +  2.54716270e+04 )
      CP(SH) =  8.24835317e+03 * (
     &       2.50000000e+00 + T * (  0.00000000e+00 
     &      + T * (  0.00000000e+00 + T * (  0.00000000e+00
     &      + T * (  0.00000000e+00 ) ) ) ) )
      H(SH2O) =  4.61497558e+02 * (
     &      T * (  2.61104720e+00 + T * (  1.57815650e-03
     &      + T * ( -3.09951460e-07 + T * (  3.33288450e-11
     &      + T * ( -1.49378702e-15 ) ) ) ) ) -2.98681670e+04 )
      CP(SH2O) =  4.61497558e+02 * (
     &       2.61104720e+00 + T * (  3.15631300e-03 
     &      + T * ( -9.29854380e-07 + T * (  1.33315380e-10
     &      + T * ( -7.46893510e-15 ) ) ) ) )
      H(SH2) =  4.12417659e+03 * (
     &      T * (  3.06670950e+00 + T * (  2.87368775e-04
     &      + T * (  4.64610633e-09 + T * ( -6.37087950e-12
     &      + T * (  5.81971480e-16 ) ) ) ) ) -8.65474120e+02 )
      CP(SH2) =  4.12417659e+03 * (
     &       3.06670950e+00 + T * (  5.74737550e-04 
     &      + T * (  1.39383190e-08 + T * ( -2.54835180e-11
     &      + T * (  2.90985740e-15 ) ) ) ) )
      H(SCO) =  2.96824105e+02 * (
     &      T * (  3.02507810e+00 + T * (  7.21344250e-04
     &      + T * ( -1.87694260e-07 + T * (  2.54645325e-11
     &      + T * ( -1.38219032e-15 ) ) ) ) ) -1.42683500e+04 )
      CP(SCO) =  2.96824105e+02 * (
     &       3.02507810e+00 + T * (  1.44268850e-03 
     &      + T * ( -5.63082780e-07 + T * (  1.01858130e-10
     &      + T * ( -6.91095160e-15 ) ) ) ) )
      H(SCO2) =  1.88915044e+02 * (
     &      T * (  4.45362280e+00 + T * (  1.57008435e-03
     &      + T * ( -4.26136833e-07 + T * (  5.98499175e-11
     &      + T * ( -3.33806640e-15 ) ) ) ) ) -4.89669610e+04 )
      CP(SCO2) =  1.88915044e+02 * (
     &       4.45362280e+00 + T * (  3.14016870e-03 
     &      + T * ( -1.27841050e-06 + T * (  2.39399670e-10
     &      + T * ( -1.66903320e-14 ) ) ) ) )
      H(SO2) =  2.59823125e+02 * (
     &      T * (  3.61221390e+00 + T * (  3.74265830e-04
     &      + T * ( -6.60688233e-08 + T * (  8.43725200e-12
     &      + T * ( -4.78147480e-16 ) ) ) ) ) -1.19781510e+03 )
      CP(SO2) =  2.59823125e+02 * (
     &       3.61221390e+00 + T * (  7.48531660e-04 
     &      + T * ( -1.98206470e-07 + T * (  3.37490080e-11
     &      + T * ( -2.39073740e-15 ) ) ) ) )
      H(SC2H2) =  3.19315616e+02 * (
     &      T * (  4.43677040e+00 + T * (  2.68801955e-03
     &      + T * ( -6.37605567e-07 + T * (  8.21594725e-11
     &      + T * ( -4.31341900e-15 ) ) ) ) ) +  2.56676640e+04 )
      CP(SC2H2) =  3.19315616e+02 * (
     &       4.43677040e+00 + T * (  5.37603910e-03 
     &      + T * ( -1.91281670e-06 + T * (  3.28637890e-10
     &      + T * ( -2.15670950e-14 ) ) ) ) )
      H(SN2) =  2.96792318e+02 * (
     &      T * (  2.85328990e+00 + T * (  8.01106400e-04
     &      + T * ( -2.09789643e-07 + T * (  2.86025550e-11
     &      + T * ( -1.56114930e-15 ) ) ) ) ) -8.90080930e+02 )
      CP(SN2) =  2.96792318e+02 * (
     &       2.85328990e+00 + T * (  1.60221280e-03 
     &      + T * ( -6.29368930e-07 + T * (  1.14410220e-10
     &      + T * ( -7.80574650e-15 ) ) ) ) )
      H(SNO) =  2.77052316e+02 * (
     &      T * (  3.18900000e+00 + T * (  6.69114050e-04
     &      + T * ( -1.76331060e-07 + T * (  2.39798330e-11
     &      + T * ( -1.29695864e-15 ) ) ) ) ) +  9.82832900e+03 )
      CP(SNO) =  2.77052316e+02 * (
     &       3.18900000e+00 + T * (  1.33822810e-03 
     &      + T * ( -5.28993180e-07 + T * (  9.59193320e-11
     &      + T * ( -6.48479320e-15 ) ) ) ) )
      H(SHCN) =  3.07756145e+02 * (
     &      T * (  3.70681210e+00 + T * (  1.66914015e-03
     &      + T * ( -3.97110667e-07 + T * (  4.99822925e-11
     &      + T * ( -2.56529040e-15 ) ) ) ) ) +  1.49626360e+04 )
      CP(SHCN) =  3.07756145e+02 * (
     &       3.70681210e+00 + T * (  3.33828030e-03 
     &      + T * ( -1.19133200e-06 + T * (  1.99929170e-10
     &      + T * ( -1.28264520e-14 ) ) ) ) )
      H(SHNCO) =  1.93284824e+02 * (
     &      T * (  5.13003900e+00 + T * (  2.17756855e-03
     &      + T * ( -5.42300733e-07 + T * (  7.00890125e-11
     &      + T * ( -3.65520740e-15 ) ) ) ) ) -1.41017870e+04 )
      CP(SHNCO) =  1.93284824e+02 * (
     &       5.13003900e+00 + T * (  4.35513710e-03 
     &      + T * ( -1.62690220e-06 + T * (  2.80356050e-10
     &      + T * ( -1.82760370e-14 ) ) ) ) )
      H(SNH3) =  4.88217264e+02 * (
     &      T * (  2.31685770e+00 + T * (  3.14207300e-03
     &      + T * ( -7.08372100e-07 + T * (  8.50467250e-11
     &      + T * ( -4.29400520e-15 ) ) ) ) ) -6.42654870e+03 )
      CP(SNH3) =  4.88217264e+02 * (
     &       2.31685770e+00 + T * (  6.28414600e-03 
     &      + T * ( -2.12511630e-06 + T * (  3.40186900e-10
     &      + T * ( -2.14700260e-14 ) ) ) ) )
      H(SN2O) =  1.88893584e+02 * (
     &      T * (  4.73066790e+00 + T * (  1.41291335e-03
     &      + T * ( -3.85270500e-07 + T * (  5.31592075e-11
     &      + T * ( -2.91281740e-15 ) ) ) ) ) +  8.16176820e+03 )
      CP(SN2O) =  1.88893584e+02 * (
     &       4.73066790e+00 + T * (  2.82582670e-03 
     &      + T * ( -1.15581150e-06 + T * (  2.12636830e-10
     &      + T * ( -1.45640870e-14 ) ) ) ) )
      ELSE
      H(SCH4) =  5.18253444e+02 * (
     &      T * (  7.78741480e-01 + T * (  8.73834200e-03
     &      + T * ( -9.27803000e-06 + T * (  7.62427000e-09
     &      + T * ( -2.44786140e-12 ) ) ) ) ) -9.82522850e+03 )
      CP(SCH4) =  5.18253444e+02 * (
     &       7.78741480e-01 + T * (  1.74766840e-02 
     &      + T * ( -2.78340900e-05 + T * (  3.04970800e-08
     &      + T * ( -1.22393070e-11 ) ) ) ) )
      H(SH) =  8.24835317e+03 * (
     &      T * (  2.50000000e+00 + T * (  0.00000000e+00
     &      + T * (  0.00000000e+00 + T * (  0.00000000e+00
     &      + T * (  0.00000000e+00 ) ) ) ) ) +  2.54716270e+04 )
      CP(SH) =  8.24835317e+03 * (
     &       2.50000000e+00 + T * (  0.00000000e+00 
     &      + T * (  0.00000000e+00 + T * (  0.00000000e+00
     &      + T * (  0.00000000e+00 ) ) ) ) )
      H(SH2O) =  4.61497558e+02 * (
     &      T * (  4.16772340e+00 + T * ( -9.05748500e-04
     &      + T * (  1.98237627e-06 + T * ( -1.21730053e-09
     &      + T * (  3.05839820e-13 ) ) ) ) ) -3.02899690e+04 )
      CP(SH2O) =  4.61497558e+02 * (
     &       4.16772340e+00 + T * ( -1.81149700e-03 
     &      + T * (  5.94712880e-06 + T * ( -4.86920210e-09
     &      + T * (  1.52919910e-12 ) ) ) ) )
      H(SH2) =  4.12417659e+03 * (
     &      T * (  3.35535140e+00 + T * (  2.50680720e-04
     &      + T * ( -7.66896933e-08 + T * ( -1.19763310e-10
     &      + T * (  9.70451700e-14 ) ) ) ) ) -1.01916260e+03 )
      CP(SH2) =  4.12417659e+03 * (
     &       3.35535140e+00 + T * (  5.01361440e-04 
     &      + T * ( -2.30069080e-07 + T * ( -4.79053240e-10
     &      + T * (  4.85225850e-13 ) ) ) ) )
      H(SCO) =  2.96824105e+02 * (
     &      T * (  3.26245170e+00 + T * (  7.55970450e-04
     &      + T * ( -1.29391840e-06 + T * (  1.39548605e-09
     &      + T * ( -4.94990240e-13 ) ) ) ) ) -1.43105390e+04 )
      CP(SCO) =  2.96824105e+02 * (
     &       3.26245170e+00 + T * (  1.51194090e-03 
     &      + T * ( -3.88175520e-06 + T * (  5.58194420e-09
     &      + T * ( -2.47495120e-12 ) ) ) ) )
      H(SCO2) =  1.88915044e+02 * (
     &      T * (  2.27572470e+00 + T * (  4.96103615e-03
     &      + T * ( -3.46970433e-06 + T * (  1.71667170e-09
     &      + T * ( -4.23456020e-13 ) ) ) ) ) -4.83731410e+04 )
      CP(SCO2) =  1.88915044e+02 * (
     &       2.27572470e+00 + T * (  9.92207230e-03 
     &      + T * ( -1.04091130e-05 + T * (  6.86668680e-09
     &      + T * ( -2.11728010e-12 ) ) ) ) )
      H(SO2) =  2.59823125e+02 * (
     &      T * (  3.78371350e+00 + T * ( -1.51168170e-03
     &      + T * (  3.31642503e-06 + T * ( -2.45472752e-09
     &      + T * (  6.60636500e-13 ) ) ) ) ) -1.06381070e+03 )
      CP(SO2) =  2.59823125e+02 * (
     &       3.78371350e+00 + T * ( -3.02336340e-03 
     &      + T * (  9.94927510e-06 + T * ( -9.81891010e-09
     &      + T * (  3.30318250e-12 ) ) ) ) )
      H(SC2H2) =  3.19315616e+02 * (
     &      T * (  2.01356220e+00 + T * (  7.59522300e-03
     &      + T * ( -5.38772967e-06 + T * (  2.26974795e-09
     &      + T * ( -3.82549200e-13 ) ) ) ) ) +  2.61244430e+04 )
      CP(SC2H2) =  3.19315616e+02 * (
     &       2.01356220e+00 + T * (  1.51904460e-02 
     &      + T * ( -1.61631890e-05 + T * (  9.07899180e-09
     &      + T * ( -1.91274600e-12 ) ) ) ) )
      H(SN2) =  2.96792318e+02 * (
     &      T * (  3.70441770e+00 + T * ( -7.10937650e-04
     &      + T * (  9.55679733e-07 + T * ( -3.00722125e-10
     &      + T * ( -2.79093540e-15 ) ) ) ) ) -1.06407950e+03 )
      CP(SN2) =  2.96792318e+02 * (
     &       3.70441770e+00 + T * ( -1.42187530e-03 
     &      + T * (  2.86703920e-06 + T * ( -1.20288850e-09
     &      + T * ( -1.39546770e-14 ) ) ) ) )
      H(SNO) =  2.77052316e+02 * (
     &      T * (  4.04595210e+00 + T * ( -1.70908915e-03
     &      + T * (  2.66063967e-06 + T * ( -1.52848290e-09
     &      + T * (  3.18381520e-13 ) ) ) ) ) +  9.74539340e+03 )
      CP(SNO) =  2.77052316e+02 * (
     &       4.04595210e+00 + T * ( -3.41817830e-03 
     &      + T * (  7.98191900e-06 + T * ( -6.11393160e-09
     &      + T * (  1.59190760e-12 ) ) ) ) )
      H(SHCN) =  3.07756145e+02 * (
     &      T * (  2.45135560e+00 + T * (  4.36041855e-03
     &      + T * ( -3.36473433e-06 + T * (  1.68139245e-09
     &      + T * ( -3.52539180e-13 ) ) ) ) ) +  1.52130020e+04 )
      CP(SHCN) =  3.07756145e+02 * (
     &       2.45135560e+00 + T * (  8.72083710e-03 
     &      + T * ( -1.00942030e-05 + T * (  6.72556980e-09
     &      + T * ( -1.76269590e-12 ) ) ) ) )
      H(SHNCO) =  1.93284824e+02 * (
     &      T * (  2.37221640e+00 + T * (  6.83202000e-03
     &      + T * ( -4.44105267e-06 + T * (  1.61188643e-09
     &      + T * ( -2.08057880e-13 ) ) ) ) ) -1.34370590e+04 )
      CP(SHNCO) =  1.93284824e+02 * (
     &       2.37221640e+00 + T * (  1.36640400e-02 
     &      + T * ( -1.33231580e-05 + T * (  6.44754570e-09
     &      + T * ( -1.04028940e-12 ) ) ) ) )
      H(SNH3) =  4.88217264e+02 * (
     &      T * (  3.77297470e+00 + T * ( -4.14878580e-04
     &      + T * (  3.93396067e-06 + T * ( -3.03171850e-09
     &      + T * (  8.35275800e-13 ) ) ) ) ) -6.69085140e+03 )
      CP(SNH3) =  4.88217264e+02 * (
     &       3.77297470e+00 + T * ( -8.29757160e-04 
     &      + T * (  1.18018820e-05 + T * ( -1.21268740e-08
     &      + T * (  4.17637900e-12 ) ) ) ) )
      H(SN2O) =  1.88893584e+02 * (
     &      T * (  2.61891960e+00 + T * (  4.32198080e-03
     &      + T * ( -2.27035413e-06 + T * (  5.56896925e-10
     &      + T * ( -1.61300660e-14 ) ) ) ) ) +  8.75901230e+03 )
      CP(SN2O) =  1.88893584e+02 * (
     &       2.61891960e+00 + T * (  8.64396160e-03 
     &      + T * ( -6.81106240e-06 + T * (  2.22758770e-09
     &      + T * ( -8.06503300e-14 ) ) ) ) )
      ENDIF
C
      END


