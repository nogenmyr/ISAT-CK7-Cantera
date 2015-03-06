C----------------------------------------------------------------------C
C
C
      SUBROUTINE CKWYP_EXT (P, T, Y, WDOT, GAS)
C
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
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT
C     WDOT   - Chemical molar production rates of the species.
C                   cgs units - moles/(cm**3*sec)
C                   Data type - real array
C                   Dimension WDOT(*) at least KK, the total number of
C                   species.
C
C  END PROLOGUE
C
C*****double precision
      use cantera
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision
C
      PARAMETER (IGLO=5,IREAC=279,KK=9)                                         
      DIMENSION Y(*), WDOT(*)
      DIMENSION XM(IREAC), RF(IREAC), RB(IREAC), W(IREAC)
      DIMENSION RKF(IGLO),XCON(KK)
      type(phase_t) :: GAS
      DATA SMALL/1.D-50/
	if (t .lt. 1200.) then
	 do i=1,9
	  wdot(i)=0.
	 end do
	 return
	end if
      call setTemperature(gas, T)
      call setMassfractions(gas, Y)
      call setPressure(gas, p/10.) ! from dyn/cm2 to Pa
      molDen = molarDensity(gas)/1000. ! kmol/m3 to mol/cm3
      call getMoleFractions(gas, XCON)
      DO i = 1, 9
        XCON(i) = XCON(i)*molDen
      ENDDO
C
C   SET LOCAL VALUES FOR THE CONCENTRATIONS
C
      XH2 = MAX( XCON(1), SMALL )                                               
      XO2 = MAX( XCON(2), SMALL )                                               
      XOH = MAX( XCON(3), SMALL )                                               
      XH2O = MAX( XCON(4), SMALL )                                              
      XCH4 = MAX( XCON(5), SMALL )                                              
      XCO = MAX( XCON(6), SMALL )                                               
      XCO2 = MAX( XCON(7), SMALL )                                              
      XNO = MAX( XCON(8), SMALL )                                               
      XN2 = MAX( XCON(9), SMALL )                                               

      BIG = 0.0
      DO 20 K = 1, KK
        IF( XCON(K) .GT. 0.0 ) THEN
          BIG = MAX( BIG, XCON(K) )
        ENDIF
        WDOT(K) = 0.0
 20   CONTINUE
C
C   EFFECTIVE THIRD BODY FOR ALL REACTIONS
C
      CALL THIRDB( XM, XH2, XO2, XOH, XH2O, XCH4, XCO, XCO2, XNO, XN2 )         
C
C   SET THE ELEMENTARY RATES
C
      CALL ELEMRATE( RF, RB, T, XM )
C
C   EXPRESSIONS FOR STEADY-STATE SPECIES
C
      XC2H = 0.D0                                                               
      XH2CN = 0.D0                                                              
      XHCNN = 0.D0                                                              
      XC2H3 = 0.D0                                                              
      XC2H5 = 0.D0                                                              
      XC = 0.D0                                                                 
      XCH3O = 0.D0                                                              
      XCH = 0.D0                                                                
      XCN = 0.D0                                                                
      XN = 0.D0                                                                 
      XC2H6 = 0.D0                                                              
      XNH = 0.D0                                                                
      XHCCO = 0.D0                                                              
      XNNH = 0.D0                                                               
      XCH2OH = 0.D0                                                             
      XNH2 = 0.D0                                                               
      XHCCOH = 0.D0                                                             
      XCH2S = 0.D0                                                              
      XNCO = 0.D0                                                               
      XC2H2 = 0.D0                                                              
      XHOCN = 0.D0                                                              
      XC2H4 = 0.D0                                                              
      XHNO = 0.D0                                                               
      XHCO = 0.D0                                                               
      XCH2CO = 0.D0                                                             
      XNO2 = 0.D0                                                               
      XCH2 = 0.D0                                                               
      XCH3OH = 0.D0                                                             
      XH2O2 = 0.D0                                                              
      XHCNO = 0.D0                                                              
      XHO2 = 0.D0                                                               
      XHNCO = 0.D0                                                              
      XHCN = 0.D0                                                               
      XCH2O = 0.D0                                                              
      XCH3 = 0.D0                                                               
      XN2O = 0.D0                                                               
      XH = 0.D0                                                                 
      XO = 0.D0                                                                 
      ADJ = 1.D0/BIG

      NITER = 50
      DO 30 N = 1, NITER

      CONMAX = 0.0

C     STEADY-STATE EXPRESSION FOR C2H       

      VOLD = XC2H       
      ABV = RF(22)*XO*XC2H2 +RB(106)*XH*XHCCO +RF(109)*XOH*XC2H2 +              
     &      RB(171)*XCO*XHCO                                                    
      DEN = RB(22)*XOH +RF(106)*XOH +RB(109)*XH2O +RF(171)*XO2                  
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XC2H = ABV/DEN                                                            
      DIFF = ABS( (XC2H-VOLD)/MAX(XC2H,VOLD,SMALL) )                            
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR H2CN      

      VOLD = XH2CN      
      ABV = RF(275)*XCH3*XN +RF(256)*XCH3*XNO +RF(237)*XH*XHCN*XM(237)          
      DEN = RB(275)*XH +RB(256)*XOH +RB(237)*XM(237)                            
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XH2CN = ABV/DEN                                                           
      DIFF = ABS( (XH2CN-VOLD)/MAX(XH2CN,VOLD,SMALL) )                          
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR HCNN      

      VOLD = XHCNN      
      ABV = RB(257)*XH*XCO*XN2 +RB(261)*XCH2*XN2 +RB(260)*XH*XHCO*XN2 +         
     &      RB(259)*XO*XHCO*XN2 +RF(241)*XCH*XN2                                
      DEN = RF(257)*XO +RF(261)*XH +RF(260)*XOH +RF(259)*XO2 +RB(241)           
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XHCNN = ABV/DEN                                                           
      DIFF = ABS( (XHCNN-VOLD)/MAX(XHCNN,VOLD,SMALL) )                          
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR C2H3      

      VOLD = XC2H3      
      ABV = RB(73)*XH2*XC2H2 +RB(111)*XH2O*XC2H2 +RB(24)*XH*XCH2CO +            
     &      RF(75)*XH*XC2H4 +RB(173)*XHCO*XCH2O +RF(112)*XOH*XC2H4 +            
     &      RF(71)*XH*XC2H2                                                     
      DEN = RF(73)*XH +RF(111)*XOH +RF(24)*XO +RB(75)*XH2 +RF(173)*XO2 +        
     &      RB(112)*XH2O +RB(71)                                                
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XC2H3 = ABV/DEN                                                           
      DIFF = ABS( (XC2H3-VOLD)/MAX(XC2H3,VOLD,SMALL) )                          
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR C2H5      

      VOLD = XC2H5      
      ABV = RF(78)*XH*XC2H6 +RF(27)*XO*XC2H6 +RF(113)*XOH*XC2H6 +RB(26)         
     &      *XCH3*XCH2O +RF(159)*XCH3*XCH3 +RF(74)*XH*XC2H4                     
      DEN = RB(78)*XH2 +RB(27)*XOH +RB(113)*XH2O +RF(26)*XO +RB(159)*XH         
     &       +RB(74)                                                            
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XC2H5 = ABV/DEN                                                           
      DIFF = ABS( (XC2H5-VOLD)/MAX(XC2H5,VOLD,SMALL) )                          
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR C         

      VOLD = XC         
      ABV = RB(90)*XH*XCO +RB(122)*XO*XCO +RF(49)*XH*XCH                        
      DEN = RF(90)*XOH +RF(122)*XO2 +RB(49)*XH2                                 
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XC = ABV/DEN                                                              
      DIFF = ABS( (XC-VOLD)/MAX(XC,VOLD,SMALL) )                                
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR CH3O      

      VOLD = XCH3O      
      ABV = RF(69)*XH*XCH3OH +RF(19)*XO*XCH3OH +RB(66)*XOH*XCH3 +RB(170)        
     &      *XHO2*XCH2O +RF(155)*XO2*XCH3 +RF(105)*XOH*XCH3OH +RF(57)*XH        
     &      *XCH2O                                                              
      DEN = RB(69)*XH2 +RB(19)*XOH +RF(66)*XH +RF(170)*XO2 +RB(155)*XO +        
     &      RB(105)*XH2O +RB(57)                                                
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XCH3O = ABV/DEN                                                           
      DIFF = ABS( (XCH3O-VOLD)/MAX(XCH3O,VOLD,SMALL) )                          
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR CH        

      VOLD = XCH        
      ABV = RB(241)*XHCNN +RB(133)*XH*XCH2CO +RB(130)*XH*XC2H4 +RB(247)         
     &      *XH*XNCO +RB(248)*XHCO*XN +RB(246)*XO*XHCN +RB(240)*XN*XHCN         
     &       +RB(49)*XH2*XC +RB(132)*XCO*XHCO +RB(125)*XO*XHCO +RF(93)          
     &      *XOH*XCH2 +RB(126)*XH*XCH2 +RB(127)*XH*XCH2O                        
      DEN = RF(241)*XN2 +RF(133)*XCH2O +RF(130)*XCH4 +RF(247)*XNO +             
     &      RF(248)*XNO +RF(246)*XNO +RF(240)*XN2 +RF(49)*XH +RF(132)           
     &      *XCO2 +RF(125)*XO2 +RB(93)*XH2O +RF(126)*XH2 +RF(127)*XH2O          
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XCH = ABV/DEN                                                             
      DIFF = ABS( (XCH-VOLD)/MAX(XCH,VOLD,SMALL) )                              
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR CN        

      VOLD = XCN        
      ABV = RB(221)*XH*XHCN +RF(233)*XO*XHCN +RB(217)*XCO*XN +RB(218)*XH        
     &      *XNCO +RB(219)*XOH*XHCN +RB(220)*XO*XNCO                            
      DEN = RF(221)*XH2 +RB(233)*XOH +RF(217)*XO +RF(218)*XOH +RF(219)          
     &      *XH2O +RF(220)*XO2                                                  
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XCN = ABV/DEN                                                             
      DIFF = ABS( (XCN-VOLD)/MAX(XCN,VOLD,SMALL) )                              
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR N         

      VOLD = XN         
      ABV = RB(275)*XH*XH2CN +RF(217)*XO*XCN +RF(191)*XH*XNH +RF(227)           
     &      *XNCO*XM(227) +RF(248)*XCH*XNO +RB(179)*XO*XNO +RF(193)*XOH         
     &      *XNH +RF(240)*XCH*XN2 +RB(178)*XO*XN2 +RB(180)*XH*XNO               
      DEN = RF(275)*XCH3 +RB(217)*XCO +RB(191)*XH2 +RB(227)*XCO                 
     &      *XM(227) +RB(248)*XHCO +RF(179)*XO2 +RB(193)*XH2O +RB(240)          
     &      *XHCN +RF(178)*XNO +RF(180)*XOH                                     
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XN = ABV/DEN                                                              
      DIFF = ABS( (XN-VOLD)/MAX(XN,VOLD,SMALL) )                                
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR C2H6      

      VOLD = XC2H6      
      ABV = RB(78)*XH2*XC2H5 +RB(27)*XOH*XC2H5 +RB(113)*XH2O*XC2H5 +            
     &      RF(158)*XCH3*XCH3                                                   
      DEN = RF(78)*XH +RF(27)*XO +RF(113)*XOH +RB(158)                          
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XC2H6 = ABV/DEN                                                           
      DIFF = ABS( (XC2H6-VOLD)/MAX(XC2H6,VOLD,SMALL) )                          
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR NH        

      VOLD = XNH        
      ABV = RF(200)*XO*XNH2 +RF(202)*XH*XNH2 +RF(232)*XO*XHCN +RB(194)          
     &      *XO*XHNO +RF(262)*XO*XHNCO +RB(191)*XH2*XN +RF(203)*XOH*XNH2        
     &       +RF(223)*XH*XNCO +RB(190)*XH*XNO +RF(208)*XO*XNNH +RB(197)         
     &      *XH2*XHNO +RB(192)*XH*XHNO +RB(193)*XH2O*XN +RB(199)*XH*XN2O        
     &                                                                          
      DEN = RB(200)*XOH +RB(202)*XH2 +RB(232)*XCO +RF(194)*XO2 +RB(262)         
     &      *XCO2 +RF(191)*XH +RB(203)*XH2O +RB(223)*XCO +RF(190)*XO +          
     &      RB(208)*XNO +RF(197)*XH2O +RF(192)*XOH +RF(193)*XOH +RF(199)        
     &      *XNO                                                                
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XNH = ABV/DEN                                                             
      DIFF = ABS( (XNH-VOLD)/MAX(XNH,VOLD,SMALL) )                              
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR HCCO      

      VOLD = XHCCO      
      ABV = RF(106)*XOH*XC2H +RF(29)*XO*XCH2CO +RB(274)*XCO*XHCNO +             
     &      RF(80)*XH*XCH2CO +RB(79)*XCH2S*XCO +RF(114)*XOH*XCH2CO +            
     &      RF(21)*XO*XC2H2 +RB(176)*XOH*XCO*XCO +RB(28)*XH*XCO*XCO             
      DEN = RB(106)*XH +RB(29)*XOH +RF(274)*XNO +RB(80)*XH2 +RF(79)*XH +        
     &      RB(114)*XH2O +RB(21)*XH +RF(176)*XO2 +RF(28)*XO                     
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XHCCO = ABV/DEN                                                           
      DIFF = ABS( (XHCCO-VOLD)/MAX(XHCCO,VOLD,SMALL) )                          
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR NNH       

      VOLD = XNNH       
      ABV = RB(207)*XOH*XN2 +RB(209)*XH2*XN2 +RB(210)*XH2O*XN2 +RB(208)         
     &      *XNH*XNO +RB(205)*XH*XN2*XM(205) +RB(206)*XHO2*XN2 +RB(204)         
     &      *XH*XN2                                                             
      DEN = RF(207)*XO +RF(209)*XH +RF(210)*XOH +RF(208)*XO +RF(205)            
     &      *XM(205) +RF(206)*XO2 +RF(204)                                      
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XNNH = ABV/DEN                                                            
      DIFF = ABS( (XNNH-VOLD)/MAX(XNNH,VOLD,SMALL) )                            
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR CH2OH     

      VOLD = XCH2OH     
      ABV = RF(68)*XH*XCH3OH +RF(104)*XOH*XCH3OH +RF(18)*XO*XCH3OH +            
     &      RB(62)*XH2O*XCH2S +RB(169)*XHO2*XCH2O +RB(61)*XOH*XCH3 +            
     &      RF(56)*XH*XCH2O                                                     
      DEN = RB(68)*XH2 +RB(104)*XH2O +RB(18)*XOH +RF(62)*XH +RF(169)*XO2        
     &       +RF(61)*XH +RB(56)                                                 
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XCH2OH = ABV/DEN                                                          
      DIFF = ABS( (XCH2OH-VOLD)/MAX(XCH2OH,VOLD,SMALL) )                        
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR NH2       

      VOLD = XNH2       
      ABV = RB(200)*XOH*XNH +RB(202)*XH2*XNH +RF(268)*XOH*XHNCO +RF(265)        
     &      *XH*XHNCO +RB(203)*XH2O*XNH +RB(201)*XH*XHNO                        
      DEN = RF(200)*XO +RF(202)*XH +RB(268)*XCO2 +RB(265)*XCO +RF(203)          
     &      *XOH +RF(201)*XO                                                    
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XNH2 = ABV/DEN                                                            
      DIFF = ABS( (XNH2-VOLD)/MAX(XNH2,VOLD,SMALL) )                            
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR HCCOH     

      VOLD = XHCCOH     
      ABV = RF(108)*XOH*XC2H2 +RB(82)*XCH2CO*XH                                 
      DEN = RB(108)*XH +RF(82)*XH                                               
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XHCCOH = ABV/DEN                                                          
      DIFF = ABS( (XHCCOH-VOLD)/MAX(XHCCOH,VOLD,SMALL) )                        
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR CH2S      

      VOLD = XCH2S      
      ABV = RB(254)*XH*XHCNO +RF(79)*XH*XHCCO +RB(252)*XH*XHNCO +RB(147)        
     &      *XCH3OH +RB(146)*XH*XCH3 +RF(62)*XH*XCH2OH +RB(145)*XH2O*XCO        
     &       +RB(153)*XCO*XCH2O +RB(148)*XCH2*XH2O +RB(144)*XH*XOH*XCO +        
     &      RB(142)*XCH2*XN2 +RF(97)*XOH*XCH3                                   
      DEN = RF(254)*XNO +RB(79)*XCO +RF(252)*XNO +RF(147)*XH2O +RF(146)         
     &      *XH2 +RB(62)*XH2O +RF(145)*XO2 +RF(153)*XCO2 +RF(148)*XH2O +        
     &      RF(144)*XO2 +RF(142)*XN2 +RB(97)*XH2O                               
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XCH2S = ABV/DEN                                                           
      DIFF = ABS( (XCH2S-VOLD)/MAX(XCH2S,VOLD,SMALL) )                          
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR NCO       

      VOLD = XNCO       
      ABV = RF(218)*XOH*XCN +RF(220)*XO2*XCN +RF(264)*XO*XHNCO +RF(231)         
     &      *XO*XHCN +RF(247)*XCH*XNO +RF(267)*XOH*XHNCO +RB(224)*XH*XCO        
     &      *XNO +RB(227)*XCO*XN*XM(227) +RB(222)*XCO*XNO +RB(223)*XCO          
     &      *XNH                                                                
      DEN = RB(218)*XH +RB(220)*XO +RB(264)*XOH +RB(231)*XH +RB(247)*XH         
     &       +RB(267)*XH2O +RF(224)*XOH +RF(227)*XM(227) +RF(222)*XO +          
     &      RF(223)*XH                                                          
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XNCO = ABV/DEN                                                            
      DIFF = ABS( (XNCO-VOLD)/MAX(XNCO,VOLD,SMALL) )                            
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR C2H2      

      VOLD = XC2H2      
      ABV = RB(108)*XH*XHCCOH +RB(22)*XOH*XC2H +RF(73)*XH*XC2H3 +RF(111)        
     &      *XOH*XC2H3 +RB(107)*XH*XCH2CO +RB(109)*XH2O*XC2H +RF(137)           
     &      *XCH2*XCH2 +RB(71)*XC2H3 +RB(21)*XH*XHCCO +RB(23)*XCH2*XCO          
      DEN = RF(108)*XOH +RF(22)*XO +RB(73)*XH2 +RB(111)*XH2O +RF(107)           
     &      *XOH +RF(109)*XOH +RB(137)*XH2 +RF(71)*XH +RF(21)*XO +RF(23)        
     &      *XO                                                                 
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XC2H2 = ABV/DEN                                                           
      DIFF = ABS( (XC2H2-VOLD)/MAX(XC2H2,VOLD,SMALL) )                          
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR HOCN      

      VOLD = XHOCN      
      ABV = RF(234)*XOH*XHCN +RB(273)*XHNCO*XH                                  
      DEN = RB(234)*XH +RF(273)*XH                                              
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XHOCN = ABV/DEN                                                           
      DIFF = ABS( (XHOCN-VOLD)/MAX(XHOCN,VOLD,SMALL) )                          
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR C2H4      

      VOLD = XC2H4      
      ABV = RB(74)*XC2H5 +RB(75)*XH2*XC2H3 +RF(130)*XCH*XCH4 +RB(25)            
     &      *XCH3*XHCO +RB(112)*XH2O*XC2H3 +RF(138)*XCH2*XCH3                   
      DEN = RF(74)*XH +RF(75)*XH +RB(130)*XH +RF(25)*XO +RF(112)*XOH +          
     &      RB(138)*XH                                                          
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XC2H4 = ABV/DEN                                                           
      DIFF = ABS( (XC2H4-VOLD)/MAX(XC2H4,VOLD,SMALL) )                          
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR HNO       

      VOLD = XHNO       
      ABV = RF(194)*XO2*XNH +RF(201)*XO*XNH2 +RF(197)*XH2O*XNH +RF(192)         
     &      *XOH*XNH +RB(216)*XHO2*XNO +RB(213)*XOH*XNO +RB(214)*XH2*XNO        
     &       +RB(215)*XH2O*XNO +RF(212)*XH*XNO*XM(212)                          
      DEN = RB(194)*XO +RB(201)*XH +RB(197)*XH2 +RB(192)*XH +RF(216)*XO2        
     &       +RF(213)*XO +RF(214)*XH +RF(215)*XOH +RB(212)*XM(212)              
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XHNO = ABV/DEN                                                            
      DIFF = ABS( (XHNO-VOLD)/MAX(XHNO,VOLD,SMALL) )                            
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR HCO       

      VOLD = XHCO       
      ABV = RF(260)*XOH*XHCNN +RF(161)*XCH3*XCH2O +RF(259)*XO2*XHCNN +          
     &      RF(248)*XCH*XNO +RF(171)*XO2*XC2H +RF(15)*XO*XCH2O +RF(173)         
     &      *XO2*XC2H3 +RB(160)*XCH4*XCO +RF(25)*XO*XC2H4 +RF(58)*XH            
     &      *XCH2O +RF(132)*XCH*XCO2 +RF(101)*XOH*XCH2O +RF(125)*XO2*XCH        
     &       +RF(7)*XO*XCH2 +RB(55)*XH2*XCO +RF(135)*XO2*XCH2 +RB(100)          
     &      *XH2O*XCO +RB(168)*XHO2*XCO +RB(167)*XH*XCO*XM(167) +RB(166)        
     &      *XH*XCO*XH2O                                                        
      DEN = RB(260)*XH*XN2 +RB(161)*XCH4 +RB(259)*XO*XN2 +RB(248)*XN +          
     &      RB(171)*XCO +RB(15)*XOH +RB(173)*XCH2O +RF(160)*XCH3 +RB(25)        
     &      *XCH3 +RB(58)*XH2 +RB(132)*XCO +RB(101)*XH2O +RB(125)*XO +          
     &      RB(7)*XH +RF(55)*XH +RB(135)*XOH +RF(100)*XOH +RF(168)*XO2 +        
     &      RF(167)*XM(167) +RF(166)*XH2O                                       
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XHCO = ABV/DEN                                                            
      DIFF = ABS( (XHCO-VOLD)/MAX(XHCO,VOLD,SMALL) )                            
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR CH2CO     

      VOLD = XCH2CO     
      ABV = RF(82)*XHCCOH*XH +RF(107)*XOH*XC2H2 +RB(29)*XOH*XHCCO +             
     &      RB(30)*XCH2*XCO2 +RF(24)*XO*XC2H3 +RB(81)*XCH3*XCO +RB(80)          
     &      *XH2*XHCCO +RF(133)*XCH*XCH2O +RB(114)*XH2O*XHCCO +RF(140)          
     &      *XCH2*XCO                                                           
      DEN = RB(82)*XH +RB(107)*XH +RF(29)*XO +RF(30)*XO +RB(24)*XH +            
     &      RF(81)*XH +RF(80)*XH +RB(133)*XH +RF(114)*XOH +RB(140)              
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XCH2CO = ABV/DEN                                                          
      DIFF = ABS( (XCH2CO-VOLD)/MAX(XCH2CO,VOLD,SMALL) )                        
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR NO2       

      VOLD = XNO2       
      ABV = RB(188)*XO2*XNO +RF(186)*XHO2*XNO +RF(187)*XO*XNO*XM(187) +         
     &      RB(189)*XOH*XNO                                                     
      DEN = RF(188)*XO +RB(186)*XOH +RB(187)*XM(187) +RF(189)*XH                
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XNO2 = ABV/DEN                                                            
      DIFF = ABS( (XNO2-VOLD)/MAX(XNO2,VOLD,SMALL) )                            
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR CH2       

      VOLD = XCH2       
      ABV = RF(261)*XH*XHCNN +RF(30)*XO*XCH2CO +RB(137)*XH2*XC2H2 +             
     &      RF(23)*XO*XC2H2 +RB(140)*XCH2CO +RB(250)*XOH*XHCN +RB(251)          
     &      *XH*XHCNO +RB(138)*XH*XC2H4 +RB(249)*XH*XHNCO +RF(148)*XCH2S        
     &      *XH2O +RF(96)*XOH*XCH3 +RB(93)*XH2O*XCH +RB(7)*XH*XHCO +            
     &      RB(92)*XH*XCH2O +RF(126)*XH2*XCH +RF(142)*XCH2S*XN2 +RB(135)        
     &      *XOH*XHCO                                                           
      DEN = RB(261)*XN2 +RB(30)*XCO2 +RF(137)*XCH2 +RB(23)*XCO +RF(140)         
     &      *XCO +RF(250)*XNO +RF(251)*XNO +RF(138)*XCH3 +RF(249)*XNO +         
     &      RB(148)*XH2O +RB(96)*XH2O +RF(93)*XOH +RF(7)*XO +RF(92)*XOH         
     &       +RB(126)*XH +RB(142)*XN2 +RF(135)*XO2                              
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XCH2 = ABV/DEN                                                            
      DIFF = ABS( (XCH2-VOLD)/MAX(XCH2,VOLD,SMALL) )                            
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR CH3OH     

      VOLD = XCH3OH     
      ABV = RB(69)*XH2*XCH3O +RB(19)*XOH*XCH3O +RF(147)*XH2O*XCH2S +            
     &      RB(68)*XH2*XCH2OH +RB(104)*XH2O*XCH2OH +RB(18)*XOH*XCH2OH +         
     &      RB(105)*XH2O*XCH3O +RF(95)*XOH*XCH3                                 
      DEN = RF(69)*XH +RF(19)*XO +RB(147) +RF(68)*XH +RF(104)*XOH +             
     &      RF(18)*XO +RF(105)*XOH +RB(95)                                      
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XCH3OH = ABV/DEN                                                          
      DIFF = ABS( (XCH3OH-VOLD)/MAX(XCH3OH,VOLD,SMALL) )                        
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR H2O2      

      VOLD = XH2O2      
      ABV = RB(47)*XH2*XHO2 +RB(88)*XH2O*XHO2 +RB(5)*XOH*XHO2 +RF(85)           
     &      *XOH*XOH +RB(89)*XH2O*XHO2                                          
      DEN = RF(47)*XH +RF(88)*XOH +RF(5)*XO +RB(85) +RF(89)*XOH                 
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XH2O2 = ABV/DEN                                                           
      DIFF = ABS( (XH2O2-VOLD)/MAX(XH2O2,VOLD,SMALL) )                          
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR HCNO      

      VOLD = XHCNO      
      ABV = RF(274)*XHCCO*XNO +RF(254)*XCH2S*XNO +RB(271)*XOH*XHCN +            
     &      RB(270)*XHNCO*XH +RF(251)*XCH2*XNO                                  
      DEN = RB(274)*XCO +RB(254)*XH +RF(271)*XH +RF(270)*XH +RB(251)*XH         
     &                                                                          
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XHCNO = ABV/DEN                                                           
      DIFF = ABS( (XHCNO-VOLD)/MAX(XHCNO,VOLD,SMALL) )                          
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR HO2       

      VOLD = XHO2       
      ABV = RF(216)*XO2*XHNO +RF(206)*XO2*XNNH +RB(186)*XOH*XNO2 +              
     &      RF(170)*XO2*XCH3O +RF(47)*XH*XH2O2 +RF(88)*XOH*XH2O2 +RF(5)         
     &      *XO*XH2O2 +RF(169)*XO2*XCH2OH +RF(33)*XH*XO2*XM(33) +RB(45)         
     &      *XH2*XO2 +RF(89)*XOH*XH2O2 +RF(36)*XH*XO2*XN2 +RF(168)*XO2          
     &      *XHCO +RB(4)*XO2*XOH +RB(46)*XOH*XOH +RF(35)*XH*XO2*XH2O +          
     &      RB(87)*XO2*XH2O                                                     
      DEN = RB(216)*XNO +RB(206)*XN2 +RF(186)*XNO +RB(170)*XCH2O +RB(47)        
     &      *XH2 +RB(88)*XH2O +RB(5)*XOH +RB(169)*XCH2O +RB(33)*XM(33) +        
     &      RF(45)*XH +RB(89)*XH2O +RB(36)*XN2 +RB(168)*XCO +RF(4)*XO +         
     &      RF(46)*XH +RB(35)*XH2O +RF(87)*XOH                                  
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XHO2 = ABV/DEN                                                            
      DIFF = ABS( (XHO2-VOLD)/MAX(XHO2,VOLD,SMALL) )                            
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR HNCO      

      VOLD = XHNCO      
      ABV = RF(235)*XOH*XHCN +RF(273)*XHOCN*XH +RF(270)*XHCNO*XH +              
     &      RB(264)*XOH*XNCO +RF(252)*XCH2S*XNO +RB(262)*XCO2*XNH +             
     &      RB(267)*XH2O*XNCO +RB(268)*XCO2*XNH2 +RB(265)*XCO*XNH2 +            
     &      RF(249)*XCH2*XNO                                                    
      DEN = RB(235)*XH +RB(273)*XH +RB(270)*XH +RF(264)*XO +RB(252)*XH +        
     &      RF(262)*XO +RF(267)*XOH +RF(268)*XOH +RF(265)*XH +RB(249)*XH        
     &                                                                          
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XHNCO = ABV/DEN                                                           
      DIFF = ABS( (XHNCO-VOLD)/MAX(XHNCO,VOLD,SMALL) )                          
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR HCN       

      VOLD = XHCN       
      ABV = RB(237)*XH2CN*XM(237) +RF(221)*XH2*XCN +RB(235)*XH*XHNCO +          
     &      RF(271)*XH*XHCNO +RB(233)*XOH*XCN +RB(234)*XH*XHOCN +RB(232)        
     &      *XCO*XNH +RF(255)*XCH3*XNO +RF(219)*XH2O*XCN +RF(250)*XCH2          
     &      *XNO +RF(246)*XCH*XNO +RB(231)*XH*XNCO +RF(240)*XCH*XN2             
      DEN = RF(237)*XH*XM(237) +RB(221)*XH +RF(235)*XOH +RB(271)*XOH +          
     &      RF(233)*XO +RF(234)*XOH +RF(232)*XO +RB(255)*XH2O +RB(219)          
     &      *XOH +RB(250)*XOH +RB(246)*XO +RF(231)*XO +RB(240)*XN               
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XHCN = ABV/DEN                                                            
      DIFF = ABS( (XHCN-VOLD)/MAX(XHCN,VOLD,SMALL) )                            
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR CH2O      

      VOLD = XCH2O      
      ABV = RF(26)*XO*XC2H5 +RB(133)*XH*XCH2CO +RF(173)*XO2*XC2H3 +             
     &      RB(161)*XCH4*XHCO +RB(57)*XCH3O +RF(170)*XO2*XCH3O +RF(83)          
     &      *XH2*XCO +RF(169)*XO2*XCH2OH +RB(56)*XCH2OH +RF(127)*XH2O           
     &      *XCH +RF(153)*XCH2S*XCO2 +RF(92)*XOH*XCH2 +RB(58)*XH2*XHCO +        
     &      RB(15)*XOH*XHCO +RF(10)*XO*XCH3 +RB(101)*XH2O*XHCO                  
      DEN = RB(26)*XCH3 +RF(133)*XCH +RB(173)*XHCO +RF(161)*XCH3 +RF(57)        
     &      *XH +RB(170)*XHO2 +RB(83) +RB(169)*XHO2 +RF(56)*XH +RB(127)         
     &      *XH +RB(153)*XCO +RB(92)*XH +RF(58)*XH +RF(15)*XO +RB(10)*XH        
     &       +RF(101)*XOH                                                       
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XCH2O = ABV/DEN                                                           
      DIFF = ABS( (XCH2O-VOLD)/MAX(XCH2O,VOLD,SMALL) )                          
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR CH3       

      VOLD = XCH3       
      ABV = RB(161)*XCH4*XHCO +RB(275)*XH*XH2CN +RB(256)*XOH*XH2CN +            
     &      RF(26)*XO*XC2H5 +RF(81)*XH*XCH2CO +RB(255)*XH2O*XHCN +              
     &      RB(160)*XCH4*XCO +RB(158)*XC2H6 +RF(25)*XO*XC2H4 +RB(159)*XH        
     &      *XC2H5 +RF(53)*XH*XCH4 +RB(138)*XH*XC2H4 +RF(98)*XOH*XCH4 +         
     &      RF(11)*XO*XCH4 +RF(66)*XH*XCH3O +RB(52)*XCH4 +RB(155)*XO            
     &      *XCH3O +RF(146)*XH2*XCH2S +RB(95)*XCH3OH +RF(61)*XH*XCH2OH +        
     &      RB(96)*XH2O*XCH2 +RB(97)*XH2O*XCH2S +RB(10)*XH*XCH2O                
      DEN = RF(161)*XCH2O +RF(275)*XN +RF(256)*XNO +RB(26)*XCH2O +RB(81)        
     &      *XCO +RF(255)*XNO +RF(160)*XHCO +RF(158)*XCH3 +RB(25)*XHCO +        
     &      RF(159)*XCH3 +RB(53)*XH2 +RF(138)*XCH2 +RB(98)*XH2O +RB(11)         
     &      *XOH +RB(66)*XOH +RF(52)*XH +RF(155)*XO2 +RB(146)*XH +RF(95)        
     &      *XOH +RB(61)*XOH +RF(96)*XOH +RF(97)*XOH +RF(10)*XO                 
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XCH3 = ABV/DEN                                                            
      DIFF = ABS( (XCH3-VOLD)/MAX(XCH3,VOLD,SMALL) )                            
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR N2O       

      VOLD = XN2O       
      ABV = RB(182)*XNO*XNO +RB(181)*XO2*XN2 +RF(199)*XNH*XNO +RB(183)          
     &      *XOH*XN2 +RB(185)*XO*XN2                                            
      DEN = RF(182)*XO +RF(181)*XO +RB(199)*XH +RF(183)*XH +RF(185)             
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XN2O = ABV/DEN                                                            
      DIFF = ABS( (XN2O-VOLD)/MAX(XN2O,VOLD,SMALL) )                            
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR H         

      VOLD = XH         
      ABV = RB(237)*XH2CN*XM(237) +RB(78)*XH2*XC2H5 +RF(221)*XH2*XCN +          
     &      RF(108)*XOH*XC2H2 +RF(260)*XOH*XHCNN +RF(106)*XOH*XC2H +            
     &      RF(235)*XOH*XHCN +RF(107)*XOH*XC2H2 +RF(257)*XO*XHCNN +             
     &      RF(90)*XOH*XC +RF(234)*XOH*XHCN +RB(271)*XOH*XHCN +RF(275)          
     &      *XCH3*XN +RB(261)*XCH2*XN2 +RB(75)*XH2*XC2H3 +RB(80)*XH2            
     &      *XHCCO +RF(205)*XNNH*XM(205) +RB(79)*XCH2S*XCO +RB(81)*XCH3         
     &      *XCO +RB(202)*XH2*XNH +RF(218)*XOH*XCN +RF(204)*XNNH +              
     &      RF(130)*XCH*XCH4 +RB(73)*XH2*XC2H2 +RF(28)*XO*XHCCO +RF(254)        
     &      *XCH2S*XNO +RF(21)*XO*XC2H2 +RF(24)*XO*XC2H3 +RF(231)*XO            
     &      *XHCN +RB(191)*XH2*XN +RB(265)*XCO*XNH2 +RF(133)*XCH*XCH2O +        
     &      RB(74)*XC2H5 +RF(201)*XO*XNH2 +RF(138)*XCH2*XCH3 +RF(252)           
     &      *XCH2S*XNO +RF(247)*XCH*XNO +RB(71)*XC2H3 +RF(3)*XH2*XO +           
     &      RB(53)*XH2*XCH3 +RF(192)*XOH*XNH +RB(223)*XCO*XNH +RF(190)          
     &      *XO*XNH +RB(209)*XH2*XN2 +RF(159)*XCH3*XCH3 +RF(224)*XOH            
     &      *XNCO +RF(180)*XOH*XN +RF(251)*XCH2*XNO +RB(58)*XH2*XHCO +          
     &      RB(69)*XH2*XCH3O +RB(189)*XOH*XNO +RB(214)*XH2*XNO +RF(249)         
     &      *XCH2*XNO +RF(146)*XH2*XCH2S +RF(84)*XH2*XOH +RB(183)*XOH           
     &      *XN2 +RB(68)*XH2*XCH2OH +RF(199)*XNH*XNO +RB(212)*XHNO              
     &      *XM(212) +RB(66)*XOH*XCH3 +RF(7)*XO*XCH2 +RF(92)*XOH*XCH2 +         
     &      RF(126)*XH2*XCH +RB(49)*XH2*XC +RF(127)*XH2O*XCH +RB(47)*XH2        
     &      *XHO2 +RB(57)*XCH3O +RB(62)*XH2O*XCH2S +RF(10)*XO*XCH3 +            
     &      RB(45)*XH2*XO2 +RB(38)*XO*XOH +RB(52)*XCH4 +RF(144)*XO2             
     &      *XCH2S +RB(56)*XCH2OH +RB(55)*XH2*XCO +RB(61)*XOH*XCH3 +            
     &      RF(99)*XOH*XCO +RB(46)*XOH*XOH +RF(167)*XHCO*XM(167) +RB(33)        
     &      *XHO2*XM(33) +RB(36)*XHO2*XN2 +RF(166)*XHCO*XH2O +RB(43)            
     &      *XH2O*XM(43) +RB(35)*XHO2*XH2O                                      
      DEN = RF(237)*XHCN*XM(237) +RF(78)*XC2H6 +RB(221)*XHCN +RB(108)           
     &      *XHCCOH +RB(260)*XHCO*XN2 +RB(106)*XHCCO +RB(235)*XHNCO +           
     &      RB(107)*XCH2CO +RB(257)*XCO*XN2 +RB(90)*XCO +RB(234)*XHOCN +        
     &      RF(271)*XHCNO +RB(275)*XH2CN +RF(261)*XHCNN +RF(75)*XC2H4 +         
     &      RF(80)*XCH2CO +RB(205)*XN2*XM(205) +RF(79)*XHCCO +RF(81)            
     &      *XCH2CO +RF(202)*XNH2 +RB(218)*XNCO +RB(204)*XN2 +RB(130)           
     &      *XC2H4 +RF(73)*XC2H3 +RB(28)*XCO*XCO +RB(254)*XHCNO +RB(21)         
     &      *XHCCO +RB(24)*XCH2CO +RB(231)*XNCO +RF(191)*XNH +RF(265)           
     &      *XHNCO +RB(133)*XCH2CO +RF(74)*XC2H4 +RB(201)*XHNO +RB(138)         
     &      *XC2H4 +RB(252)*XHNCO +RB(247)*XNCO +RF(71)*XC2H2 +RB(3)*XOH        
     &       +RF(53)*XCH4 +RB(192)*XHNO +RF(223)*XNCO +RB(190)*XNO +            
     &      RF(209)*XNNH +RB(159)*XC2H5 +RB(224)*XCO*XNO +RB(180)*XNO +         
     &      RB(251)*XHCNO +RF(58)*XCH2O +RF(69)*XCH3OH +RF(189)*XNO2 +          
     &      RF(214)*XHNO +RB(249)*XHNCO +RB(146)*XCH3 +RB(84)*XH2O +            
     &      RF(183)*XN2O +RF(68)*XCH3OH +RB(199)*XN2O +RF(212)*XNO              
     &      *XM(212) +RF(66)*XCH3O +RB(7)*XHCO +RB(92)*XCH2O +RB(126)           
     &      *XCH2 +RF(49)*XCH +RB(127)*XCH2O +RF(47)*XH2O2 +RF(57)*XCH2O        
     &       +RF(62)*XCH2OH +RB(10)*XCH2O +RF(45)*XHO2 +RF(38)*XO2 +            
     &      RF(52)*XCH3 +RB(144)*XOH*XCO +RF(56)*XCH2O +RF(55)*XHCO +           
     &      RF(61)*XCH2OH +RB(99)*XCO2 +RF(46)*XHO2 +RB(167)*XCO                
     &      *XM(167) +RF(33)*XO2*XM(33) +RF(36)*XO2*XN2 +RB(166)*XCO            
     &      *XH2O +RF(43)*XOH*XM(43) +RF(35)*XO2*XH2O                           
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XH = ABV/DEN                                                              
      DIFF = ABS( (XH-VOLD)/MAX(XH,VOLD,SMALL) )                                
      CONMAX = MAX( CONMAX, DIFF )

C     STEADY-STATE EXPRESSION FOR O         

      VOLD = XO         
      ABV = RB(22)*XOH*XC2H +RB(27)*XOH*XC2H5 +RB(233)*XOH*XCN +RB(257)         
     &      *XH*XCO*XN2 +RB(29)*XOH*XHCCO +RF(259)*XO2*XHCNN +RB(264)           
     &      *XOH*XNCO +RB(30)*XCH2*XCO2 +RB(200)*XOH*XNH +RB(28)*XH*XCO         
     &      *XCO +RF(179)*XO2*XN +RF(122)*XO2*XC +RF(220)*XO2*XCN +             
     &      RF(194)*XO2*XNH +RB(231)*XH*XNCO +RB(232)*XCO*XNH +RB(23)           
     &      *XCH2*XCO +RB(21)*XH*XHCCO +RB(262)*XCO2*XNH +RB(24)*XH             
     &      *XCH2CO +RB(217)*XCO*XN +RB(25)*XCH3*XHCO +RB(26)*XCH3*XCH2O        
     &       +RB(201)*XH*XHNO +RB(188)*XO2*XNO +RB(181)*XO2*XN2 +RB(207)        
     &      *XOH*XN2 +RB(190)*XH*XNO +RB(7)*XH*XHCO +RF(246)*XCH*XNO +          
     &      RB(222)*XCO*XNO +RB(213)*XOH*XNO +RB(19)*XOH*XCH3O +RB(182)         
     &      *XNO*XNO +RB(18)*XOH*XCH2OH +RB(208)*XNH*XNO +RB(15)*XOH            
     &      *XHCO +RB(3)*XH*XOH +RF(185)*XN2O +RF(125)*XO2*XCH +RB(10)          
     &      *XH*XCH2O +RB(11)*XOH*XCH3 +RB(187)*XNO2*XM(187) +RF(178)*XN        
     &      *XNO +RF(155)*XO2*XCH3 +RB(5)*XOH*XHO2 +RF(86)*XOH*XOH +            
     &      RB(4)*XO2*XOH +RF(38)*XH*XO2                                        
      DEN = RF(22)*XC2H2 +RF(27)*XC2H6 +RF(233)*XHCN +RF(257)*XHCNN +           
     &      RF(29)*XCH2CO +RB(259)*XHCO*XN2 +RF(264)*XHNCO +RF(30)              
     &      *XCH2CO +RF(200)*XNH2 +RF(28)*XHCCO +RB(179)*XNO +RB(122)           
     &      *XCO +RB(220)*XNCO +RB(194)*XHNO +RF(231)*XHCN +RF(232)*XHCN        
     &       +RF(23)*XC2H2 +RF(21)*XC2H2 +RF(262)*XHNCO +RF(24)*XC2H3 +         
     &      RF(217)*XCN +RF(25)*XC2H4 +RF(26)*XC2H5 +RF(201)*XNH2 +             
     &      RF(188)*XNO2 +RF(181)*XN2O +RF(207)*XNNH +RF(190)*XNH +RF(7)        
     &      *XCH2 +RB(246)*XHCN +RF(222)*XNCO +RF(213)*XHNO +RF(19)             
     &      *XCH3OH +RF(182)*XN2O +RF(18)*XCH3OH +RF(208)*XNNH +RF(15)          
     &      *XCH2O +RF(3)*XH2 +RB(185)*XN2 +RB(125)*XHCO +RF(10)*XCH3 +         
     &      RF(11)*XCH4 +RF(187)*XNO*XM(187) +RB(178)*XN2 +RB(155)*XCH3O        
     &       +RF(5)*XH2O2 +RB(86)*XH2O +RF(4)*XHO2 +RB(38)*XOH                  
      IF( DEN .LT. 1.0 ) DEN = MAX( ADJ*ABV, DEN, SMALL )
      XO = ABV/DEN                                                              
      DIFF = ABS( (XO-VOLD)/MAX(XO,VOLD,SMALL) )                                
      CONMAX = MAX( CONMAX, DIFF )

      IF( CONMAX .LT. 1.D-5 ) GO TO 35

 30   CONTINUE

 35   CONTINUE
C
C   NET PRODUCTION RATES FOR SKELETAL MECHANISM
C
      CALL NETRATE( W, RF, RB, XM, XH2, XH, XO, XO2, XOH, XH2O, XHO2,           
     &              XH2O2, XC, XCH, XCH2, XCH2S, XCH3, XCH4, XCO, XCO2,         
     &              XHCO, XCH2O, XCH2OH, XCH3O, XCH3OH, XC2H, XC2H2,            
     &              XC2H3, XC2H4, XC2H5, XC2H6, XHCCO, XCH2CO, XHCCOH,          
     &              XN, XNH, XNH2, XNNH, XNO, XNO2, XN2O, XHNO, XCN,            
     &              XHCN, XH2CN, XHCNN, XHCNO, XHOCN, XHNCO, XNCO, XN2 )        
C
C   GLOBAL RATES FOR REDUCED MECHANISM
C
      RKF(1) = +W(4) -W(33) -W(35) -W(36) +W(38) +W(45) + 2.00*W(46)            
     &         - 2.00*W(85) +W(87) +W(122) +W(125) +W(135) +W(144)              
     &         +W(145) +W(155) -W(168) -W(169) -W(170) +W(171) +W(173)          
     &         +W(176) +W(179) -W(181) +W(186) -W(187) +W(189) -W(192)          
     &         -W(197) -W(201) -W(206) -W(212) +W(213) +W(214) +W(215)          
     &         +W(220) +W(259)                                                  
      RKF(2) = +W(7) +W(11) +W(15) +W(33) +W(35) +W(36) -W(38) +W(43)           
     &         -W(46) +W(53) +W(58) +W(61) +W(62) +W(66) +W(71)                 
     &         + 2.00*W(85) +W(93) +W(98) +W(101) -W(122) -W(125)               
     &         -W(126) -W(127) +W(130) +W(137) +W(138) +W(140) -W(144)          
     &         - 2.00*W(155) +W(158) -W(160) -W(166) -W(167) +W(169)            
     &         +W(170) -W(173) -W(176) +W(180) -W(183) - 2.00*W(185)            
     &         + 2.00*W(187) -W(189) +W(190) +W(192) +W(197) +W(199)            
     &         +W(201) -W(204) -W(205) +W(212) -W(220) +W(222) -W(227)          
     &         -W(240) -W(246) -W(247) -W(248) -W(256) -W(259)                  
      RKF(3) = +W(10) -W(11) -W(15) -W(25) +W(52) -W(53) -W(56) -W(57)          
     &         -W(58) -W(75) +W(83) +W(92) -W(98) -W(101) -W(112)               
     &         +W(127) -W(133) +W(138) +W(153) +W(158) +W(159) +W(160)          
     &         +W(169) +W(170) +W(173)                                          
      RKF(4) = +W(23) +W(28) -W(38) -W(46) -W(56) -W(57) +W(61) +W(66)          
     &         -W(75) +W(81) +W(85) -W(95) -W(97) -W(99) -W(112)                
     &         -W(122) -W(125) +W(132) -W(133) -W(135) -W(137) -W(140)          
     &         +W(142) +W(146) +W(148) + 2.00*W(153) - 2.00*W(155)              
     &         +W(169) +W(170) +W(180) -W(183) -W(185) +W(187) -W(189)          
     &         +W(190) +W(192) +W(197) +W(199) +W(201) -W(220) +W(222)          
     &         +W(224) -W(240) -W(246) -W(247) -W(248) -W(249) -W(250)          
     &         -W(251) -W(255) -W(256) -W(259) -W(262) -W(268)                  
      RKF(5) = +W(178) +W(181) +W(183) +W(185) -W(208) -W(240)                  
C
C   REDUCED MECHANISM
C
C   3H2 + O2 + CO2 =  3H2O + CO                                                 
C  H2 +  2OH =  2H2O                                                            
C   3H2 + CO = H2O + CH4                                                        
C  H2 + CO2 = H2O + CO                                                          
C   3H2 + CO2 +  2NO =  3H2O + CO + N2                                          
C
C   SPECIES PRODUCTION RATES
C
C H2                                                                            
      WDOT(1) = - 3.00*RKF(1) - RKF(2) - 3.00*RKF(3) - RKF(4)                   
     &          - 3.00*RKF(5)                                                   
C O2                                                                            
      WDOT(2) = - RKF(1)                                                        
C OH                                                                            
      WDOT(3) = - 2.00*RKF(2)                                                   
C H2O                                                                           
      WDOT(4) = + 3.00*RKF(1) + 2.00*RKF(2) + RKF(3) + RKF(4)                   
     &          + 3.00*RKF(5)                                                   
C CH4                                                                           
      WDOT(5) = + RKF(3)                                                        
C CO                                                                            
      WDOT(6) = + RKF(1) - RKF(3) + RKF(4) + RKF(5)                             
C CO2                                                                           
      WDOT(7) = - RKF(1) - RKF(4) - RKF(5)                                      
C NO                                                                            
      WDOT(8) = - 2.00*RKF(5)                                                   
C N2                                                                            
      WDOT(9) = + RKF(5)                                                        

      RETURN
      END

      SUBROUTINE THIRDB( XM, XH2, XO2, XOH, XH2O, XCH4, XCO, XCO2, XNO,         
     &                   XN2 )                                                  

C*****double precision
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
C       IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision

      DIMENSION XM(*)
C
C   EFFECTIVE THIRD BODY FOR ALL REACTIONS
C
      XM(33) = +XH2+XOH+XCH4+ 0.750*XCO+ 1.500*XCO2+XNO                         
      XM(43) = + 0.730*XH2+XO2+XOH+ 3.650*XH2O+ 2.000*XCH4+XCO+XCO2+XNO         
     &         +XN2                                                             
      XM(52) = + 2.000*XH2+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4+ 1.500*XCO           
     &         + 2.000*XCO2+XNO+XN2                                             
      XM(56) = + 2.000*XH2+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4+ 1.500*XCO           
     &         + 2.000*XCO2+XNO+XN2                                             
      XM(57) = + 2.000*XH2+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4+ 1.500*XCO           
     &         + 2.000*XCO2+XNO+XN2                                             
      XM(71) = + 2.000*XH2+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4+ 1.500*XCO           
     &         + 2.000*XCO2+XNO+XN2                                             
      XM(74) = + 2.000*XH2+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4+ 1.500*XCO           
     &         + 2.000*XCO2+XNO+XN2                                             
      XM(83) = + 2.000*XH2+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4+ 1.500*XCO           
     &         + 2.000*XCO2+XNO+XN2                                             
      XM(85) = + 2.000*XH2+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4+ 1.500*XCO           
     &         + 2.000*XCO2+XNO+XN2                                             
      XM(95) = + 2.000*XH2+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4+ 1.500*XCO           
     &         + 2.000*XCO2+XNO+XN2                                             
      XM(140) = + 2.000*XH2+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4+ 1.500*XCO          
     &          + 2.000*XCO2+XNO+XN2                                            
      XM(147) = + 2.000*XH2+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4+ 1.500*XCO          
     &          + 2.000*XCO2+XNO+XN2                                            
      XM(158) = + 2.000*XH2+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4+ 1.500*XCO          
     &          + 2.000*XCO2+XNO+XN2                                            
      XM(167) = + 2.000*XH2+XO2+XOH+ 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2          
     &          +XNO+XN2                                                        
      XM(185) = + 2.000*XH2+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4+ 1.500*XCO          
     &          + 2.000*XCO2+XNO+XN2                                            
      XM(187) = + 2.000*XH2+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4+ 1.500*XCO          
     &          + 2.000*XCO2+XNO+XN2                                            
      XM(205) = + 2.000*XH2+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4+ 1.500*XCO          
     &          + 2.000*XCO2+XNO+XN2                                            
      XM(212) = + 2.000*XH2+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4+ 1.500*XCO          
     &          + 2.000*XCO2+XNO+XN2                                            
      XM(227) = + 2.000*XH2+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4+ 1.500*XCO          
     &          + 2.000*XCO2+XNO+XN2                                            
      XM(237) = + 2.000*XH2+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4+ 1.500*XCO          
     &          + 2.000*XCO2+XNO+XN2                                            
      XM(241) = + 2.000*XH2+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4+ 1.500*XCO          
     &          + 2.000*XCO2+XNO+XN2                                            

      RETURN
      END

      SUBROUTINE ELEMRATE( RF, RB, T, XM )

C*****double precision
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
C       IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision

      DIMENSION RF(*), RB(*), XM(*)

C   RATE CONSTANTS FOR SKELETAL MECHANSIM

      RUC = 8.31451D0/4.184D0
      ALOGT = DLOG(T)
      RTR = 1.0D3/(RUC*T)

      RF(3) = 5.000D+04 * EXP(  2.670*ALOGT -  6.29000*RTR )                    
      RB(3) = EXP(+(-2.23489D+07)/T/T/T+( 2.80492D+05)/T/T                      
     &        +(-3.94618D+03)/T+( 2.82285D+01)+( 2.13355D-03)*T                 
     &        +(-4.33965D-07)*T*T+( 4.54247D-11)*T*T*T)                         
      RF(4) = 2.000D+13                                                         
      RB(4) = EXP(+( 5.92451D+05)/T/T/T+(-3.41134D+04)/T/T                      
     &        +(-2.65087D+04)/T+( 3.02310D+01)+( 4.08842D-04)*T                 
     &        +(-8.62910D-08)*T*T+( 9.70706D-12)*T*T*T)                         
      RF(5) = 9.630D+06 * EXP(  2.000*ALOGT -  4.00000*RTR )                    
      RB(5) = EXP(+(-1.30325D+07)/T/T/T+( 1.53786D+05)/T/T                      
     &        +(-1.04088D+04)/T+( 2.69031D+01)+( 2.32010D-03)*T                 
     &        +(-4.64925D-07)*T*T+( 4.48793D-11)*T*T*T)                         
      RF(7) = 8.000D+13                                                         
      RB(7) = EXP(+(-1.03571D+07)/T/T/T+( 1.21403D+05)/T/T                      
     &        +(-4.64114D+04)/T+( 3.49933D+01)+(-7.79879D-04)*T                 
     &        +( 2.10360D-07)*T*T+(-2.33865D-11)*T*T*T)                         
      RF(10) = 8.430D+13                                                        
      RB(10) = EXP(+(-1.45276D+07)/T/T/T+( 1.77812D+05)/T/T                     
     &         +(-3.52415D+04)/T+( 3.59934D+01)+(-7.92692D-04)*T                
     &         +( 2.06574D-07)*T*T+(-2.21450D-11)*T*T*T)                        
      RF(11) = 1.020D+09 * EXP(  1.500*ALOGT -  8.60000*RTR )                   
      RB(11) = EXP(+( 1.95738D+07)/T/T/T+(-2.19284D+05)/T/T                     
     &         +(-2.40285D+03)/T+( 2.50874D+01)+( 1.87962D-03)*T                
     &         +(-2.71623D-07)*T*T+( 2.25030D-11)*T*T*T)                        
      RF(15) = 3.900D+13 * EXP(              -  3.54000*RTR )                   
      RB(15) = EXP(+( 1.73749D+07)/T/T/T+(-2.19790D+05)/T/T                     
     &         +(-7.96059D+03)/T+( 2.63453D+01)+( 6.78941D-04)*T                
     &         +(-1.02404D-07)*T*T+( 6.90075D-12)*T*T*T)                        
      RF(18) = 3.880D+05 * EXP(  2.500*ALOGT -  3.10000*RTR )                   
      RB(18) = EXP(+( 8.23274D+06)/T/T/T+(-9.96787D+04)/T/T                     
     &         +(-4.53584D+03)/T+( 2.44522D+01)+( 2.63451D-03)*T                
     &         +(-4.77102D-07)*T*T+( 4.40625D-11)*T*T*T)                        
      RF(19) = 1.300D+05 * EXP(  2.500*ALOGT -  5.00000*RTR )                   
      RB(19) = EXP(+(-1.89709D+07)/T/T/T+( 2.19199D+05)/T/T                     
     &         +(-3.04453D+03)/T+( 2.70742D+01)+( 2.08815D-03)*T                
     &         +(-3.44251D-07)*T*T+( 2.74880D-11)*T*T*T)                        
      RF(21) = 1.020D+07 * EXP(  2.000*ALOGT -  1.90000*RTR )                   
      RB(21) = EXP(+(-1.46891D+07)/T/T/T+( 1.84717D+05)/T/T                     
     &         +(-1.20731D+04)/T+( 3.00062D+01)+( 1.00615D-03)*T                
     &         +(-8.90544D-08)*T*T+( 3.91193D-12)*T*T*T)                        
      RF(22) = 4.600D+19 * EXP( -1.410*ALOGT - 28.95000*RTR )                   
      RB(22) = EXP(+( 1.01289D+07)/T/T/T+(-1.92905D+05)/T/T                     
     &         +( 2.19999D+03)/T+( 3.03345D+01)+(-1.97314D-04)*T                
     &         +(-2.89392D-08)*T*T+( 7.84947D-12)*T*T*T)                        
      RF(23) = 1.020D+07 * EXP(  2.000*ALOGT -  1.90000*RTR )                   
      RB(23) = EXP(+(-2.88069D+07)/T/T/T+( 2.98392D+05)/T/T                     
     &         +(-2.59020D+04)/T+( 2.56523D+01)+( 2.72769D-03)*T                
     &         +(-5.32393D-07)*T*T+( 5.50080D-11)*T*T*T)                        
      RF(24) = 3.000D+13                                                        
      RB(24) = EXP(+( 1.06166D+07)/T/T/T+(-1.35563D+05)/T/T                     
     &         +(-4.50278D+04)/T+( 3.52174D+01)+(-5.71421D-04)*T                
     &         +( 1.66026D-07)*T*T+(-1.94609D-11)*T*T*T)                        
      RF(25) = 1.920D+07 * EXP(  1.830*ALOGT -  0.22000*RTR )                   
      RB(25) = EXP(+(-1.15317D+07)/T/T/T+( 6.96782D+04)/T/T                     
     &         +(-1.41140D+04)/T+( 2.30488D+01)+( 2.51186D-03)*T                
     &         +(-4.79218D-07)*T*T+( 4.75268D-11)*T*T*T)                        
      RF(26) = 1.320D+14                                                        
      RB(26) = EXP(+(-5.58411D+06)/T/T/T+( 1.20929D+04)/T/T                     
     &         +(-3.94571D+04)/T+( 3.11836D+01)+( 9.58585D-04)*T                
     &         +(-2.04178D-07)*T*T+( 2.15354D-11)*T*T*T)                        
      RF(27) = 8.980D+07 * EXP(  1.920*ALOGT -  5.69000*RTR )                   
      RB(27) = EXP(+( 1.19559D+06)/T/T/T+(-1.20562D+04)/T/T                     
     &         +(-4.09241D+03)/T+( 2.49785D+01)+( 2.32012D-03)*T                
     &         +(-4.30151D-07)*T*T+( 4.01341D-11)*T*T*T)                        
      RF(28) = 1.000D+14                                                        
      RB(28) = EXP(+(-1.11426D+07)/T/T/T+( 6.32901D+04)/T/T                     
     &         +(-5.19314D+04)/T+( 3.02585D+01)+( 1.25135D-03)*T                
     &         +(-2.39710D-07)*T*T+( 2.50928D-11)*T*T*T)                        
      RF(29) = 1.000D+13 * EXP(              -  8.00000*RTR )                   
      RB(29) = EXP(+( 1.13827D+07)/T/T/T+(-1.48793D+05)/T/T                     
     &         +(-1.51461D+03)/T+( 2.54177D+01)+( 6.88118D-04)*T                
     &         +(-9.85504D-08)*T*T+( 5.79617D-12)*T*T*T)                        
      RF(30) = 1.750D+12 * EXP(              -  1.35000*RTR )                   
      RB(30) = EXP(+(-8.73406D+05)/T/T/T+(-4.42754D+03)/T/T                     
     &         +(-2.49305D+04)/T+( 2.70007D+01)+( 8.25676D-04)*T                
     &         +(-1.73313D-07)*T*T+( 1.79626D-11)*T*T*T)                        
      RF(33) = 2.800D+18 * EXP( -0.860*ALOGT                )                   
      RB(33) = EXP(+( 1.95829D+06)/T/T/T+(-8.96287D+02)/T/T                     
     &         +(-2.42563D+04)/T+( 3.69564D+01)+(-5.18190D-04)*T                
     &         +( 2.28756D-08)*T*T+( 4.24094D-14)*T*T*T)                        
      RF(35) = 9.380D+18 * EXP( -0.760*ALOGT                )                   
      RB(35) = EXP(+( 1.17500D+06)/T/T/T+( 9.02041D+03)/T/T                     
     &         +(-2.43202D+04)/T+( 3.88448D+01)+(-4.37890D-04)*T                
     &         +( 7.12568D-09)*T*T+( 1.62640D-12)*T*T*T)                        
      RF(36) = 3.750D+20 * EXP( -1.720*ALOGT                )                   
      RB(36) = EXP(+( 8.69457D+06)/T/T/T+(-8.61799D+04)/T/T                     
     &         +(-2.37070D+04)/T+( 3.60112D+01)+(-1.20877D-03)*T                
     &         +( 1.58325D-07)*T*T+(-1.35799D-11)*T*T*T)                        
      RF(38) = 8.300D+13 * EXP(              - 14.41300*RTR )                   
      RB(38) = EXP(+( 6.16394D+06)/T/T/T+(-8.58077D+04)/T/T                     
     &         +( 1.65979D+03)/T+( 2.81102D+01)+( 8.34323D-04)*T                
     &         +(-2.27071D-07)*T*T+( 2.60752D-11)*T*T*T)                        
      RF(43) = 2.200D+22 * EXP( -2.000*ALOGT                )                   
      RB(43) = EXP(+( 2.90282D+06)/T/T/T+(-4.02691D+04)/T/T                     
     &         +(-5.90676D+04)/T+( 4.05688D+01)+(-1.19039D-03)*T                
     &         +( 1.28742D-07)*T*T+(-8.88854D-12)*T*T*T)                        
      RF(45) = 2.800D+13 * EXP(              -  1.06800*RTR )                   
      RB(45) = EXP(+( 2.02756D+06)/T/T/T+(-4.98290D+04)/T/T                     
     &         +(-2.79704D+04)/T+( 3.12979D+01)+( 4.19299D-04)*T                
     &         +(-7.28496D-08)*T*T+( 6.57497D-12)*T*T*T)                        
      RF(46) = 1.340D+14 * EXP(              -  0.63500*RTR )                   
      RB(46) = EXP(+( 6.75639D+06)/T/T/T+(-1.19921D+05)/T/T                     
     &         +(-1.79156D+04)/T+( 2.81934D+01)+( 1.24317D-03)*T                
     &         +(-3.13362D-07)*T*T+( 3.57823D-11)*T*T*T)                        
      RF(47) = 1.210D+07 * EXP(  2.000*ALOGT -  5.20000*RTR )                   
      RB(47) = EXP(+(-1.15974D+07)/T/T/T+( 1.38071D+05)/T/T                     
     &         +(-1.19370D+04)/T+( 2.78618D+01)+( 2.33056D-03)*T                
     &         +(-4.51484D-07)*T*T+( 4.17473D-11)*T*T*T)                        
      RF(49) = 1.100D+14                                                        
      RB(49) = EXP(+( 4.09971D+06)/T/T/T+(-3.92623D+04)/T/T                     
     &         +(-1.17301D+04)/T+( 3.32169D+01)+( 7.79596D-05)*T                
     &         +( 3.94975D-08)*T*T+(-7.78065D-12)*T*T*T)                        
      RFH52 = 1.270D+16 * EXP( -0.630*ALOGT -  0.38300*RTR )                    
      RFL52 = 2.477D+33 * EXP( -4.760*ALOGT -  2.44000*RTR )                    
      PR = RFL52*XM(52)/RFH52                                                   
      PRLOG = DLOG10(PR)
      F4 = 0.7830
      F5 =    74.00
      F6 =  2941.00
      F7 =  6964.00
      FC = (1.0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(FC)
      CPRLOG = PRLOG - (0.4 + 0.67*FCLOG)
      X = CPRLOG/(0.75-1.27*FCLOG-0.14*CPRLOG)
      F = 10.0**( FCLOG/(1.0+X*X) )
      PCOR = PR*F/(1.0+PR)
      RF(52) = RFH52*PCOR                                                       
      RB(52) = EXP(+(-3.05740D+07)/T/T/T+( 3.55834D+05)/T/T                     
     &         +(-5.39874D+04)/T+( 3.85737D+01)+(-5.99776D-04)*T                
     &         +(-6.42667D-08)*T*T+( 1.46495D-11)*T*T*T)                        
      RB(52) = RB(52)*PCOR                                                      
      RF(53) = 6.600D+08 * EXP(  1.620*ALOGT - 10.84000*RTR )                   
      RB(53) = EXP(+( 2.00690D+07)/T/T/T+(-2.23100D+05)/T/T                     
     &         +(-4.53102D+03)/T+( 2.61978D+01)+( 1.98643D-03)*T                
     &         +(-2.77081D-07)*T*T+( 2.12717D-11)*T*T*T)                        
      RF(55) = 7.340D+13                                                        
      RB(55) = EXP(+( 1.05819D+07)/T/T/T+(-1.37230D+05)/T/T                     
     &         +(-4.39297D+04)/T+( 3.19804D+01)+( 9.01373D-04)*T                
     &         +(-1.92155D-07)*T*T+( 1.76229D-11)*T*T*T)                        
      RFH56 = 5.400D+11 * EXP(  0.454*ALOGT -  3.60000*RTR )                    
      RFL56 = 1.270D+32 * EXP( -4.820*ALOGT -  6.53000*RTR )                    
      PR = RFL56*XM(56)/RFH56                                                   
      PRLOG = DLOG10(PR)
      F4 = 0.7187
      F5 =   103.00
      F6 =  1291.00
      F7 =  4160.00
      FC = (1.0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(FC)
      CPRLOG = PRLOG - (0.4 + 0.67*FCLOG)
      X = CPRLOG/(0.75-1.27*FCLOG-0.14*CPRLOG)
      F = 10.0**( FCLOG/(1.0+X*X) )
      PCOR = PR*F/(1.0+PR)
      RF(56) = RFH56*PCOR                                                       
      RB(56) = EXP(+( 7.46324D+06)/T/T/T+(-5.82414D+04)/T/T                     
     &         +(-1.64005D+04)/T+( 2.94363D+01)+( 3.63048D-04)*T                
     &         +(-1.20138D-07)*T*T+( 1.32405D-11)*T*T*T)                        
      RB(56) = RB(56)*PCOR                                                      
      RFH57 = 5.400D+11 * EXP(  0.454*ALOGT -  2.60000*RTR )                    
      RFL57 = 2.200D+30 * EXP( -4.800*ALOGT -  5.56000*RTR )                    
      PR = RFL57*XM(57)/RFH57                                                   
      PRLOG = DLOG10(PR)
      F4 = 0.7580
      F5 =    94.00
      F6 =  1555.00
      F7 =  4200.00
      FC = (1.0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(FC)
      CPRLOG = PRLOG - (0.4 + 0.67*FCLOG)
      X = CPRLOG/(0.75-1.27*FCLOG-0.14*CPRLOG)
      F = 10.0**( FCLOG/(1.0+X*X) )
      PCOR = PR*F/(1.0+PR)
      RF(57) = RFH57*PCOR                                                       
      RB(57) = EXP(+(-1.97404D+07)/T/T/T+( 2.60636D+05)/T/T                     
     &         +(-1.34499D+04)/T+( 3.31518D+01)+(-1.83310D-04)*T                
     &         +( 1.27130D-08)*T*T+(-3.33405D-12)*T*T*T)                        
      RB(57) = RB(57)*PCOR                                                      
      RF(58) = 2.300D+10 * EXP(  1.050*ALOGT -  3.27500*RTR )                   
      RB(58) = EXP(+( 1.05855D+07)/T/T/T+(-1.31380D+05)/T/T                     
     &         +(-9.42217D+03)/T+( 2.67734D+01)+( 1.53255D-03)*T                
     &         +(-2.54337D-07)*T*T+( 2.04006D-11)*T*T*T)                        
      RF(61) = 1.200D+13                                                        
      RB(61) = EXP(+(-6.77346D+05)/T/T/T+(-2.42745D+04)/T/T                     
     &         +(-1.77373D+03)/T+( 2.67779D+01)+( 1.37544D-03)*T                
     &         +(-3.56806D-07)*T*T+( 3.94676D-11)*T*T*T)                        
      RF(62) = 6.000D+12                                                        
      RB(62) = EXP(+(-4.43259D+06)/T/T/T+( 3.43571D+03)/T/T                     
     &         +(-1.58107D+03)/T+( 2.58398D+01)+( 1.82492D-03)*T                
     &         +(-4.46096D-07)*T*T+( 4.74746D-11)*T*T*T)                        
      RF(66) = 3.200D+13                                                        
      RB(66) = EXP(+( 2.65263D+07)/T/T/T+(-3.43152D+05)/T/T                     
     &         +(-4.22115D+03)/T+( 2.40432D+01)+( 1.92180D-03)*T                
     &         +(-4.89656D-07)*T*T+( 5.60422D-11)*T*T*T)                        
      RF(68) = 1.700D+07 * EXP(  2.100*ALOGT -  4.87000*RTR )                   
      RB(68) = EXP(+( 1.28010D+07)/T/T/T+(-1.55061D+05)/T/T                     
     &         +(-6.09539D+03)/T+( 2.62451D+01)+( 2.32377D-03)*T                
     &         +(-4.00661D-07)*T*T+( 3.45945D-11)*T*T*T)                        
      RF(69) = 4.200D+06 * EXP(  2.100*ALOGT -  4.87000*RTR )                   
      RB(69) = EXP(+(-1.44027D+07)/T/T/T+( 1.63817D+05)/T/T                     
     &         +(-3.64797D+03)/T+( 2.85624D+01)+( 1.77741D-03)*T                
     &         +(-2.67810D-07)*T*T+( 1.80199D-11)*T*T*T)                        
      RFH71 = 5.600D+12 * EXP(              -  2.40000*RTR )                    
      RFL71 = 3.800D+40 * EXP( -7.270*ALOGT -  7.22000*RTR )                    
      PR = RFL71*XM(71)/RFH71                                                   
      PRLOG = DLOG10(PR)
      F4 = 0.7507
      F5 =    98.50
      F6 =  1302.00
      F7 =  4167.00
      FC = (1.0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(FC)
      CPRLOG = PRLOG - (0.4 + 0.67*FCLOG)
      X = CPRLOG/(0.75-1.27*FCLOG-0.14*CPRLOG)
      F = 10.0**( FCLOG/(1.0+X*X) )
      PCOR = PR*F/(1.0+PR)
      RF(71) = RFH71*PCOR                                                       
      RB(71) = EXP(+(-2.52082D+07)/T/T/T+( 3.21012D+05)/T/T                     
     &         +(-1.98449D+04)/T+( 2.98925D+01)+(-1.35317D-04)*T                
     &         +(-4.03962D-08)*T*T+( 9.26855D-12)*T*T*T)                        
      RB(71) = RB(71)*PCOR                                                      
      RF(73) = 3.000D+13                                                        
      RB(73) = EXP(+( 2.24578D+07)/T/T/T+(-2.86454D+05)/T/T                     
     &         +(-3.36014D+04)/T+( 3.11510D+01)+( 7.27005D-04)*T                
     &         +(-1.45027D-07)*T*T+( 1.09712D-11)*T*T*T)                        
      RFH74 = 1.080D+12 * EXP(  0.454*ALOGT -  1.82000*RTR )                    
      RFL74 = 1.200D+42 * EXP( -7.620*ALOGT -  6.97000*RTR )                    
      PR = RFL74*XM(74)/RFH74                                                   
      PRLOG = DLOG10(PR)
      F4 = 0.9753
      F5 =   210.00
      F6 =   984.00
      F7 =  4374.00
      FC = (1.0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(FC)
      CPRLOG = PRLOG - (0.4 + 0.67*FCLOG)
      X = CPRLOG/(0.75-1.27*FCLOG-0.14*CPRLOG)
      F = 10.0**( FCLOG/(1.0+X*X) )
      PCOR = PR*F/(1.0+PR)
      RF(74) = RFH74*PCOR                                                       
      RB(74) = EXP(+(-1.67300D+07)/T/T/T+( 1.91195D+05)/T/T                     
     &         +(-1.97183D+04)/T+( 3.08448D+01)+( 3.50640D-04)*T                
     &         +(-1.54781D-07)*T*T+( 2.06667D-11)*T*T*T)                        
      RB(74) = RB(74)*PCOR                                                      
      RF(75) = 1.325D+06 * EXP(  2.530*ALOGT - 12.24000*RTR )                   
      RB(75) = EXP(+(-1.10989D+06)/T/T/T+( 9.16469D+03)/T/T                     
     &         +(-3.20326D+03)/T+( 2.58675D+01)+( 2.80322D-03)*T                
     &         +(-5.11949D-07)*T*T+( 4.69485D-11)*T*T*T)                        
      RF(78) = 1.150D+08 * EXP(  1.900*ALOGT -  7.53000*RTR )                   
      RB(78) = EXP(+( 2.78736D+06)/T/T/T+(-2.97552D+04)/T/T                     
     &         +(-5.92988D+03)/T+( 2.58204D+01)+( 2.31452D-03)*T                
     &         +(-4.13559D-07)*T*T+( 3.66852D-11)*T*T*T)                        
      RF(79) = 1.000D+14                                                        
      RB(79) = EXP(+(-2.25001D+07)/T/T/T+( 1.96976D+05)/T/T                     
     &         +(-9.59617D+03)/T+( 2.88550D+01)+( 1.67052D-03)*T                
     &         +(-4.39047D-07)*T*T+( 5.15413D-11)*T*T*T)                        
      RF(80) = 5.000D+13 * EXP(              -  8.00000*RTR )                   
      RB(80) = EXP(+( 1.28178D+07)/T/T/T+(-1.64508D+05)/T/T                     
     &         +(-2.43893D+03)/T+( 2.77576D+01)+( 6.98575D-04)*T                
     &         +(-8.51090D-08)*T*T+( 2.66409D-12)*T*T*T)                        
      RF(81) = 1.130D+13 * EXP(              -  3.42800*RTR )                   
      RB(81) = EXP(+(-1.59395D+07)/T/T/T+( 1.28264D+05)/T/T                     
     &         +(-1.80334D+04)/T+( 2.51875D+01)+( 1.74353D-03)*T                
     &         +(-4.35699D-07)*T*T+( 4.87500D-11)*T*T*T)                        
      RF(82) = 1.000D+13                                                        
      RB(82) = EXP(+(-1.11016D+07)/T/T/T+( 1.42257D+05)/T/T                     
     &         +(-1.57232D+04)/T+( 2.98741D+01)+( 6.61450D-05)*T                
     &         +( 5.98073D-09)*T*T+(-6.75568D-13)*T*T*T)                        
      RFH83 = 4.300D+07 * EXP(  1.500*ALOGT - 79.60000*RTR )                    
      RFL83 = 5.070D+27 * EXP( -3.420*ALOGT - 84.35000*RTR )                    
      PR = RFL83*XM(83)/RFH83                                                   
      PRLOG = DLOG10(PR)
      F4 = 0.9320
      F5 =   197.00
      F6 =  1540.00
      F7 = 10300.00
      FC = (1.0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(FC)
      CPRLOG = PRLOG - (0.4 + 0.67*FCLOG)
      X = CPRLOG/(0.75-1.27*FCLOG-0.14*CPRLOG)
      F = 10.0**( FCLOG/(1.0+X*X) )
      PCOR = PR*F/(1.0+PR)
      RF(83) = RFH83*PCOR                                                       
      RB(83) = EXP(+(-4.38917D+07)/T/T/T+( 5.56044D+05)/T/T                     
     &         +(-4.22194D+04)/T+( 3.25901D+01)+( 2.05416D-04)*T                
     &         +(-1.40555D-07)*T*T+( 2.26081D-11)*T*T*T)                        
      RB(83) = RB(83)*PCOR                                                      
      RF(84) = 2.160D+08 * EXP(  1.510*ALOGT -  3.43000*RTR )                   
      RB(84) = EXP(+(-2.18402D+07)/T/T/T+( 2.73249D+05)/T/T                     
     &         +(-1.07968D+04)/T+( 3.15026D+01)+( 1.03645D-03)*T                
     &         +(-2.38657D-07)*T*T+( 2.64699D-11)*T*T*T)                        
      RFH85 = 7.400D+13 * EXP( -0.370*ALOGT                )                    
      RFL85 = 2.300D+18 * EXP( -0.900*ALOGT +  1.70000*RTR )                    
      PR = RFL85*XM(85)/RFH85                                                   
      PRLOG = DLOG10(PR)
      F4 = 0.7346
      F5 =    94.00
      F6 =  1756.00
      F7 =  5182.00
      FC = (1.0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(FC)
      CPRLOG = PRLOG - (0.4 + 0.67*FCLOG)
      X = CPRLOG/(0.75-1.27*FCLOG-0.14*CPRLOG)
      F = 10.0**( FCLOG/(1.0+X*X) )
      PCOR = PR*F/(1.0+PR)
      RF(85) = RFH85*PCOR                                                       
      RB(85) = EXP(+(-1.06770D+07)/T/T/T+( 1.78051D+05)/T/T                     
     &         +(-2.63633D+04)/T+( 3.64487D+01)+(-1.67315D-03)*T                
     &         +( 3.22699D-07)*T*T+(-3.14707D-11)*T*T*T)                        
      RB(85) = RB(85)*PCOR                                                      
      RF(86) = 3.570D+04 * EXP(  2.400*ALOGT +  2.11000*RTR )                   
      RB(86) = EXP(+(-2.73763D+07)/T/T/T+( 3.45792D+05)/T/T                     
     &         +(-9.50175D+03)/T+( 2.95716D+01)+( 1.76158D-03)*T                
     &         +(-3.65390D-07)*T*T+( 3.74354D-11)*T*T*T)                        
      RF(87) = 2.900D+13 * EXP(              +  0.50000*RTR )                   
      RB(87) = EXP(+(-7.98495D+06)/T/T/T+( 7.36775D+04)/T/T                     
     &         +(-3.52878D+04)/T+( 3.33863D+01)+( 2.43223D-04)*T                
     &         +(-7.36830D-08)*T*T+( 9.12658D-12)*T*T*T)                        
      RF(88) = 1.750D+12 * EXP(              -  0.32000*RTR )                   
      RB(88) = EXP(+(-5.94416D+06)/T/T/T+( 6.32430D+04)/T/T                     
     &         +(-1.63103D+04)/T+( 2.82097D+01)+( 5.48486D-04)*T                
     &         +(-1.37318D-07)*T*T+( 1.26190D-11)*T*T*T)                        
      RF(89) = 5.800D+14 * EXP(              -  9.56000*RTR )                   
      RB(89) = EXP(+(-5.94416D+06)/T/T/T+( 6.32430D+04)/T/T                     
     &         +(-2.09600D+04)/T+( 3.40131D+01)+( 5.48486D-04)*T                
     &         +(-1.37318D-07)*T*T+( 1.26190D-11)*T*T*T)                        
      RF(90) = 5.000D+13                                                        
      RB(90) = EXP(+(-7.90154D+06)/T/T/T+( 8.65219D+04)/T/T                     
     &         +(-7.83597D+04)/T+( 3.57483D+01)+(-5.09576D-04)*T                
     &         +( 1.23328D-07)*T*T+(-1.16033D-11)*T*T*T)                        
      RF(92) = 2.000D+13                                                        
      RB(92) = EXP(+(-2.77320D+07)/T/T/T+( 3.41193D+05)/T/T                     
     &         +(-4.02322D+04)/T+( 3.85562D+01)+(-1.45882D-03)*T                
     &         +( 3.12764D-07)*T*T+(-3.02872D-11)*T*T*T)                        
      RF(93) = 1.130D+07 * EXP(  2.000*ALOGT -  3.00000*RTR )                   
      RB(93) = EXP(+(-2.02166D+07)/T/T/T+( 2.43039D+05)/T/T                     
     &         +(-1.20690D+04)/T+( 3.05545D+01)+( 1.99349D-03)*T                
     &         +(-4.47011D-07)*T*T+( 4.47198D-11)*T*T*T)                        
      RFH95 = 6.300D+13                                                         
      RFL95 = 2.700D+38 * EXP( -6.300*ALOGT -  3.10000*RTR )                    
      PR = RFL95*XM(95)/RFH95                                                   
      PRLOG = DLOG10(PR)
      F4 = 0.2105
      F5 =    83.50
      F6 =  5398.00
      F7 =  8370.00
      FC = (1.0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(FC)
      CPRLOG = PRLOG - (0.4 + 0.67*FCLOG)
      X = CPRLOG/(0.75-1.27*FCLOG-0.14*CPRLOG)
      F = 10.0**( FCLOG/(1.0+X*X) )
      PCOR = PR*F/(1.0+PR)
      RF(95) = RFH95*PCOR                                                       
      RB(95) = EXP(+(-3.13232D+07)/T/T/T+( 4.22145D+05)/T/T                     
     &         +(-4.81613D+04)/T+( 4.04402D+01)+(-1.42122D-03)*T                
     &         +( 2.41294D-07)*T*T+(-2.05585D-11)*T*T*T)                        
      RB(95) = RB(95)*PCOR                                                      
      RF(96) = 5.600D+07 * EXP(  1.600*ALOGT -  5.42000*RTR )                   
      RB(96) = EXP(+(-7.90555D+06)/T/T/T+( 1.03077D+05)/T/T                     
     &         +(-7.78939D+03)/T+( 2.74931D+01)+( 1.78531D-03)*T                
     &         +(-3.45581D-07)*T*T+( 3.29056D-11)*T*T*T)                        
      RF(97) = 2.501D+13                                                        
      RB(97) = EXP(+(-3.75524D+06)/T/T/T+( 2.77102D+04)/T/T                     
     &         +( 1.92658D+02)/T+( 3.06053D+01)+( 4.49484D-04)*T                
     &         +(-8.92907D-08)*T*T+( 8.00698D-12)*T*T*T)                        
      RF(98) = 1.000D+08 * EXP(  1.600*ALOGT -  3.12000*RTR )                   
      RB(98) = EXP(+( 1.02131D+07)/T/T/T+(-1.01577D+05)/T/T                     
     &         +(-8.73980D+03)/T+( 2.62281D+01)+( 1.79430D-03)*T                
     &         +(-2.74764D-07)*T*T+( 2.35065D-11)*T*T*T)                        
      RF(99) = 4.760D+07 * EXP(  1.228*ALOGT -  0.07000*RTR )                   
      RB(99) = EXP(+(-7.75714D+06)/T/T/T+( 1.52467D+05)/T/T                     
     &         +(-1.37529D+04)/T+( 3.37009D+01)+(-5.97901D-04)*T                
     &         +( 1.75167D-07)*T*T+(-1.94782D-11)*T*T*T)                        
      RF(100) = 5.000D+13                                                       
      RB(100) = EXP(+( 5.69407D+05)/T/T/T+(-1.37231D+04)/T/T                    
     &          +(-5.20361D+04)/T+( 3.36498D+01)+( 7.25298D-04)*T               
     &          +(-1.92988D-07)*T*T+( 2.01745D-11)*T*T*T)                       
      RF(101) = 3.430D+09 * EXP(  1.180*ALOGT +  0.44700*RTR )                  
      RB(101) = EXP(+(-4.45301D+05)/T/T/T+( 5.01807D+03)/T/T                    
     &          +(-1.57386D+04)/T+( 2.78069D+01)+( 1.46086D-03)*T               
     &          +(-2.75646D-07)*T*T+( 2.50114D-11)*T*T*T)                       
      RF(104) = 1.440D+06 * EXP(  2.000*ALOGT +  0.84000*RTR )                  
      RB(104) = EXP(+( 3.57178D+06)/T/T/T+(-4.14713D+04)/T/T                    
     &          +(-1.12645D+04)/T+( 2.51505D+01)+( 2.06739D-03)*T               
     &          +(-3.85744D-07)*T*T+( 3.55621D-11)*T*T*T)                       
      RF(105) = 6.300D+06 * EXP(  2.000*ALOGT -  1.50000*RTR )                  
      RB(105) = EXP(+(-2.36319D+07)/T/T/T+( 2.77406D+05)/T/T                    
     &          +(-9.99464D+03)/T+( 3.03418D+01)+( 1.52103D-03)*T               
     &          +(-2.52894D-07)*T*T+( 1.89875D-11)*T*T*T)                       
      RF(106) = 2.000D+13                                                       
      RB(106) = EXP(+( 1.89217D+06)/T/T/T+( 3.94622D+04)/T/T                    
     &          +(-2.57072D+04)/T+( 3.62692D+01)+(-1.53477D-03)*T               
     &          +( 4.76958D-07)*T*T+(-5.79517D-11)*T*T*T)                       
      RF(107) = 2.180D-04 * EXP(  4.500*ALOGT +  1.00000*RTR )                  
      RB(107) = EXP(+(-4.56540D+07)/T/T/T+( 5.81427D+05)/T/T                    
     &          +(-1.47216D+04)/T+( 2.69374D+01)+( 2.32553D-03)*T               
     &          +(-3.84252D-07)*T*T+( 3.77156D-11)*T*T*T)                       
      RF(108) = 5.040D+05 * EXP(  2.300*ALOGT - 13.50000*RTR )                  
      RB(108) = EXP(+(-1.73200D+07)/T/T/T+( 2.21003D+05)/T/T                    
     &          +(-4.88993D+03)/T+( 3.36121D+01)+( 4.92784D-04)*T               
     &          +(-4.37346D-08)*T*T+( 3.54330D-12)*T*T*T)                       
      RF(109) = 3.370D+07 * EXP(  2.000*ALOGT - 14.00000*RTR )                  
      RB(109) = EXP(+(-2.51587D+07)/T/T/T+( 2.53046D+05)/T/T                    
     &          +(-1.48552D+03)/T+( 2.83426D+01)+( 2.37529D-03)*T               
     &          +(-5.53404D-07)*T*T+( 6.12832D-11)*T*T*T)                       
      RF(111) = 5.000D+12                                                       
      RB(111) = EXP(+( 1.24452D+07)/T/T/T+(-1.62947D+05)/T/T                    
     &          +(-4.17078D+04)/T+( 3.14126D+01)+( 5.50930D-04)*T               
     &          +(-1.45861D-07)*T*T+( 1.35228D-11)*T*T*T)                       
      RF(112) = 3.600D+06 * EXP(  2.000*ALOGT -  2.50000*RTR )                  
      RB(112) = EXP(+(-6.97098D+06)/T/T/T+( 8.01127D+04)/T/T                    
     &          +(-6.06982D+03)/T+( 2.53197D+01)+( 2.20156D-03)*T               
     &          +(-4.29308D-07)*T*T+( 4.11049D-11)*T*T*T)                       
      RF(113) = 3.540D+06 * EXP(  2.120*ALOGT -  0.87000*RTR )                  
      RB(113) = EXP(+(-8.94839D+06)/T/T/T+( 1.15568D+05)/T/T                    
     &          +(-1.08253D+04)/T+( 2.58875D+01)+( 2.31510D-03)*T               
     &          +(-4.49042D-07)*T*T+( 4.27216D-11)*T*T*T)                       
      RF(114) = 7.500D+12 * EXP(              -  2.00000*RTR )                  
      RB(114) = EXP(+( 2.80532D+06)/T/T/T+(-4.10017D+04)/T/T                    
     &          +(-7.52602D+03)/T+( 2.79138D+01)+( 5.22500D-04)*T               
     &          +(-8.59425D-08)*T*T+( 5.21569D-12)*T*T*T)                       
      RF(122) = 5.800D+13 * EXP(              -  0.57600*RTR )                  
      RB(122) = EXP(+(-1.73760D+06)/T/T/T+( 7.14167D+02)/T/T                    
     &          +(-6.97369D+04)/T+( 3.19571D+01)+( 3.24747D-04)*T               
     &          +(-1.03744D-07)*T*T+( 1.44719D-11)*T*T*T)                       
      RF(125) = 3.300D+13                                                       
      RB(125) = EXP(+(-8.21982D+06)/T/T/T+( 9.86815D+04)/T/T                    
     &          +(-3.72474D+04)/T+( 3.22251D+01)+(-4.98666D-04)*T               
     &          +( 1.27909D-07)*T*T+(-1.09317D-11)*T*T*T)                       
      RF(126) = 1.107D+08 * EXP(  1.790*ALOGT -  1.67000*RTR )                  
      RB(126) = EXP(+(-1.94826D+07)/T/T/T+( 2.56311D+05)/T/T                    
     &          +(-8.07975D+02)/T+( 3.20096D+01)+( 8.73802D-04)*T               
     &          +(-1.50745D-07)*T*T+( 1.78651D-11)*T*T*T)                       
      RF(127) = 1.713D+13 * EXP(              +  0.75500*RTR )                  
      RB(127) = EXP(+(-2.31813D+07)/T/T/T+( 2.96488D+05)/T/T                    
     &          +(-3.05702D+04)/T+( 3.76745D+01)+(-1.84631D-03)*T               
     &          +( 4.44776D-07)*T*T+(-4.33272D-11)*T*T*T)                       
      RF(130) = 6.000D+13                                                       
      RB(130) = EXP(+( 1.41369D+07)/T/T/T+(-7.17483D+04)/T/T                    
     &          +(-3.04427D+04)/T+( 3.70705D+01)+(-1.70024D-03)*T               
     &          +( 5.10600D-07)*T*T+(-5.68035D-11)*T*T*T)                       
      RF(132) = 3.400D+12 * EXP(              -  0.69000*RTR )                  
      RB(132) = EXP(+(-1.62454D+07)/T/T/T+( 1.53799D+05)/T/T                    
     &          +(-3.35739D+04)/T+( 2.62121D+01)+( 2.50995D-04)*T               
     &          +(-1.35963D-08)*T*T+( 1.92271D-12)*T*T*T)                       
      RF(133) = 9.460D+13 * EXP(              +  0.51500*RTR )                  
      RB(133) = EXP(+( 2.52301D+07)/T/T/T+(-2.43100D+05)/T/T                    
     &          +(-3.73564D+04)/T+( 3.74812D+01)+(-1.39291D-03)*T               
     &          +( 3.78509D-07)*T*T+(-4.28569D-11)*T*T*T)                       
      RF(135) = 1.320D+13 * EXP(              -  1.50000*RTR )                  
      RB(135) = EXP(+(-4.19319D+06)/T/T/T+( 3.55954D+04)/T/T                    
     &          +(-3.82535D+04)/T+( 2.92518D+01)+( 5.44439D-05)*T               
     &          +(-1.67115D-08)*T*T+( 2.68878D-12)*T*T*T)                       
      RF(137) = 3.200D+13                                                       
      RB(137) = EXP(+( 1.33659D+07)/T/T/T+(-1.15884D+05)/T/T                    
     &          +(-6.66726D+04)/T+( 3.82034D+01)+(-1.00020D-03)*T               
     &          +( 2.35600D-07)*T*T+(-2.90916D-11)*T*T*T)                       
      RF(138) = 4.000D+13                                                       
      RB(138) = EXP(+(-1.31596D+07)/T/T/T+( 2.33201D+05)/T/T                    
     &          +(-3.35769D+04)/T+( 4.04543D+01)+(-1.82225D-03)*T               
     &          +( 4.01354D-07)*T*T+(-4.19261D-11)*T*T*T)                       
      RFH140 = 8.100D+11 * EXP(  0.500*ALOGT -  4.51000*RTR )                   
      RFL140 = 2.690D+33 * EXP( -5.110*ALOGT -  7.09500*RTR )                   
      PR = RFL140*XM(140)/RFH140                                                
      PRLOG = DLOG10(PR)
      F4 = 0.5907
      F5 =   275.00
      F6 =  1226.00
      F7 =  5185.00
      FC = (1.0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(FC)
      CPRLOG = PRLOG - (0.4 + 0.67*FCLOG)
      X = CPRLOG/(0.75-1.27*FCLOG-0.14*CPRLOG)
      F = 10.0**( FCLOG/(1.0+X*X) )
      PCOR = PR*F/(1.0+PR)
      RF(140) = RFH140*PCOR                                                     
      RB(140) = EXP(+(-5.36694D+06)/T/T/T+( 1.34975D+05)/T/T                    
     &          +(-4.25853D+04)/T+( 3.96140D+01)+(-1.42693D-03)*T               
     &          +( 2.64275D-07)*T*T+(-2.56005D-11)*T*T*T)                       
      RB(140) = RB(140)*PCOR                                                    
      RF(142) = 1.500D+13 * EXP(              -  0.60000*RTR )                  
      RB(142) = EXP(+( 8.38231D+06)/T/T/T+(-8.33007D+04)/T/T                    
     &          +(-4.53466D+03)/T+( 2.93664D+01)+( 5.10252D-05)*T               
     &          +(-4.29182D-09)*T*T+(-4.45233D-13)*T*T*T)                       
      RF(144) = 2.800D+13                                                       
      RB(144) = EXP(+( 1.75215D+07)/T/T/T+(-2.19493D+05)/T/T                    
     &          +(-3.34226D+04)/T+( 2.84270D+01)+( 4.15154D-04)*T               
     &          +(-2.77347D-08)*T*T+(-3.73257D-13)*T*T*T)                       
      RF(145) = 1.200D+13                                                       
      RB(145) = EXP(+( 4.75853D+06)/T/T/T+(-6.14284D+04)/T/T                    
     &          +(-9.37675D+04)/T+( 3.02906D+01)+( 8.30767D-04)*T               
     &          +(-2.13992D-07)*T*T+( 2.24181D-11)*T*T*T)                       
      RF(146) = 7.000D+13                                                       
      RB(146) = EXP(+(-6.25728D+06)/T/T/T+( 9.57963D+04)/T/T                    
     &          +(-8.29904D+03)/T+( 3.41778D+01)+(-6.25560D-04)*T               
     &          +( 8.84572D-08)*T*T+(-5.45537D-12)*T*T*T)                       
      RFH147 = 2.000D+13                                                        
      RFL147 = 2.700D+38 * EXP( -6.300*ALOGT -  3.10000*RTR )                   
      PR = RFL147*XM(147)/RFH147                                                
      PRLOG = DLOG10(PR)
      F4 = 0.1507
      F5 =   134.00
      F6 =  2383.00
      F7 =  7265.00
      FC = (1.0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(FC)
      CPRLOG = PRLOG - (0.4 + 0.67*FCLOG)
      X = CPRLOG/(0.75-1.27*FCLOG-0.14*CPRLOG)
      F = 10.0**( FCLOG/(1.0+X*X) )
      PCOR = PR*F/(1.0+PR)
      RF(147) = RFH147*PCOR                                                     
      RB(147) = EXP(+(-2.75679D+07)/T/T/T+( 3.94434D+05)/T/T                    
     &          +(-4.83540D+04)/T+( 3.95377D+01)+(-1.87070D-03)*T               
     &          +( 3.30585D-07)*T*T+(-2.85654D-11)*T*T*T)                       
      RB(147) = RB(147)*PCOR                                                    
      RF(148) = 3.000D+13                                                       
      RB(148) = EXP(+( 8.38231D+06)/T/T/T+(-8.33007D+04)/T/T                    
     &          +(-4.23273D+03)/T+( 3.00595D+01)+( 5.10252D-05)*T               
     &          +(-4.29183D-09)*T*T+(-4.45232D-13)*T*T*T)                       
      RF(153) = 1.400D+13                                                       
      RB(153) = EXP(+(-2.12114D+07)/T/T/T+( 2.27202D+05)/T/T                    
     &          +(-3.15315D+04)/T+( 2.95469D+01)+( 1.76189D-04)*T               
     &          +(-6.01038D-08)*T*T+( 8.19717D-12)*T*T*T)                       
      RF(155) = 2.675D+13 * EXP(              - 28.80000*RTR )                  
      RB(155) = EXP(+(-2.03624D+07)/T/T/T+( 2.57344D+05)/T/T                    
     &          +(-1.35883D+03)/T+( 3.40314D+01)+(-1.08747D-03)*T               
     &          +( 2.62585D-07)*T*T+(-2.99669D-11)*T*T*T)                       
      RFH158 = 2.120D+16 * EXP( -0.970*ALOGT -  0.62000*RTR )                   
      RFL158 = 1.770D+50 * EXP( -9.670*ALOGT -  6.22000*RTR )                   
      PR = RFL158*XM(158)/RFH158                                                
      PRLOG = DLOG10(PR)
      F4 = 0.5325
      F5 =   151.00
      F6 =  1038.00
      F7 =  4970.00
      FC = (1.0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(FC)
      CPRLOG = PRLOG - (0.4 + 0.67*FCLOG)
      X = CPRLOG/(0.75-1.27*FCLOG-0.14*CPRLOG)
      F = 10.0**( FCLOG/(1.0+X*X) )
      PCOR = PR*F/(1.0+PR)
      RF(158) = RFH158*PCOR                                                     
      RB(158) = EXP(+(-2.17658D+07)/T/T/T+( 3.22258D+05)/T/T                    
     &          +(-4.67883D+04)/T+( 4.25668D+01)+(-2.72732D-03)*T               
     &          +( 4.92413D-07)*T*T+(-4.53947D-11)*T*T*T)                       
      RB(158) = RB(158)*PCOR                                                    
      RF(159) = 4.990D+12 * EXP(  0.100*ALOGT - 10.60000*RTR )                  
      RB(159) = EXP(+(-9.72674D+06)/T/T/T+( 1.75635D+05)/T/T                    
     &          +(-1.18240D+03)/T+( 3.51761D+01)+(-1.67098D-03)*T               
     &          +( 3.95001D-07)*T*T+(-4.20964D-11)*T*T*T)                       
      RF(160) = 2.648D+13                                                       
      RB(160) = EXP(+(-2.21763D+07)/T/T/T+( 2.46521D+05)/T/T                    
     &          +(-4.58882D+04)/T+( 3.60767D+01)+( 2.15798D-04)*T               
     &          +(-1.70223D-07)*T*T+( 2.20119D-11)*T*T*T)                       
      RF(161) = 3.320D+03 * EXP(  2.810*ALOGT -  5.86000*RTR )                  
      RB(161) = EXP(+(-3.59586D+07)/T/T/T+( 4.26904D+05)/T/T                    
     &          +(-1.38056D+04)/T+( 2.80950D+01)+( 2.26025D-03)*T               
     &          +(-5.09604D-07)*T*T+( 5.26679D-11)*T*T*T)                       
      RF(166) = 2.244D+18 * EXP( -1.000*ALOGT - 17.00000*RTR )                  
      RB(166) = EXP(+( 2.11652D+07)/T/T/T+(-2.70955D+05)/T/T                    
     &          +( 3.92816D+02)/T+( 3.48570D+01)+(-4.93315D-04)*T               
     &          +( 1.50768D-07)*T*T+(-1.84567D-11)*T*T*T)                       
      RF(167) = 1.870D+17 * EXP( -1.000*ALOGT - 17.00000*RTR )                  
      RB(167) = EXP(+( 2.11652D+07)/T/T/T+(-2.70955D+05)/T/T                    
     &          +( 3.92816D+02)/T+( 3.23721D+01)+(-4.93315D-04)*T               
     &          +( 1.50768D-07)*T*T+(-1.84567D-11)*T*T*T)                       
      RF(168) = 7.600D+12 * EXP(              -  0.40000*RTR )                  
      RB(168) = EXP(+( 8.55436D+06)/T/T/T+(-8.74006D+04)/T/T                    
     &          +(-1.66980D+04)/T+( 2.93780D+01)+( 4.82074D-04)*T               
     &          +(-1.19305D-07)*T*T+( 1.10480D-11)*T*T*T)                       
      RF(169) = 1.800D+13 * EXP(              -  0.90000*RTR )                  
      RB(169) = EXP(+(-1.57974D+07)/T/T/T+( 1.87651D+05)/T/T                    
     &          +(-1.09595D+04)/T+( 3.15071D+01)+( 1.73904D-04)*T               
     &          +(-6.39407D-08)*T*T+( 7.61556D-12)*T*T*T)                       
      RF(170) = 4.280D-13 * EXP(  7.600*ALOGT +  3.53000*RTR )                  
      RB(170) = EXP(+(-4.81236D+07)/T/T/T+( 6.22442D+05)/T/T                    
     &          +(-1.60316D+04)/T+( 2.04228D+01)+( 6.82306D-03)*T               
     &          +(-1.39379D-06)*T*T+( 1.44574D-10)*T*T*T)                       
      RF(171) = 5.000D+13 * EXP(              -  1.50000*RTR )                  
      RB(171) = EXP(+(-1.64188D+07)/T/T/T+( 1.88732D+05)/T/T                    
     &          +(-7.77896D+04)/T+( 3.18721D+01)+( 2.41220D-04)*T               
     &          +( 1.69071D-08)*T*T+(-4.16689D-12)*T*T*T)                       
      RF(173) = 3.980D+12 * EXP(              +  0.24000*RTR )                  
      RB(173) = EXP(+(-2.28333D+07)/T/T/T+( 2.06219D+05)/T/T                    
     &          +(-4.45389D+04)/T+( 2.89944D+01)+( 3.22826D-04)*T               
     &          +(-8.45743D-08)*T*T+( 1.24644D-11)*T*T*T)                       
      RF(176) = 1.600D+12 * EXP(              -  0.85400*RTR )                  
      RB(176) = EXP(+(-4.97861D+06)/T/T/T+(-2.25177D+04)/T/T                    
     &          +(-4.34485D+04)/T+( 2.21836D+01)+( 2.08567D-03)*T               
     &          +(-4.66782D-07)*T*T+( 5.11680D-11)*T*T*T)                       
      RF(178) = 3.500D+13 * EXP(              -  0.33000*RTR )                  
      RB(178) = EXP(+( 3.97977D+05)/T/T/T+(-1.42620D+04)/T/T                    
     &          +(-3.79249D+04)/T+( 3.23101D+01)+( 2.55881D-04)*T               
     &          +(-7.71311D-08)*T*T+( 9.53934D-12)*T*T*T)                       
      RF(179) = 2.650D+12 * EXP(              -  6.40000*RTR )                  
      RB(179) = EXP(+(-1.63369D+06)/T/T/T+( 1.30574D+03)/T/T                    
     &          +(-1.90701D+04)/T+( 2.68215D+01)+( 1.44637D-04)*T               
     &          +(-4.26817D-08)*T*T+( 5.92650D-12)*T*T*T)                       
      RF(180) = 7.333D+13 * EXP(              -  1.12000*RTR )                  
      RB(180) = EXP(+(-7.79764D+06)/T/T/T+( 8.71135D+04)/T/T                    
     &          +(-2.53257D+04)/T+( 3.40816D+01)+(-6.89686D-04)*T               
     &          +( 1.84390D-07)*T*T+(-2.01487D-11)*T*T*T)                       
      RF(181) = 1.400D+12 * EXP(              - 10.81000*RTR )                  
      RB(181) = EXP(+(-8.69622D+06)/T/T/T+( 6.43383D+04)/T/T                    
     &          +(-4.52908D+04)/T+( 2.56512D+01)+( 9.81265D-04)*T               
     &          +(-2.14167D-07)*T*T+( 2.27213D-11)*T*T*T)                       
      RF(182) = 2.900D+13 * EXP(              - 23.15000*RTR )                  
      RB(182) = EXP(+(-1.07279D+07)/T/T/T+( 7.99061D+04)/T/T                    
     &          +(-2.95912D+04)/T+( 2.57743D+01)+( 8.70021D-04)*T               
     &          +(-1.79718D-07)*T*T+( 1.91085D-11)*T*T*T)                       
      RF(183) = 4.400D+14 * EXP(              - 18.88000*RTR )                  
      RB(183) = EXP(+(-2.53227D+06)/T/T/T+(-2.14694D+04)/T/T                    
     &          +(-4.04391D+04)/T+( 2.74619D+01)+( 1.81559D-03)*T               
     &          +(-4.41239D-07)*T*T+( 4.87966D-11)*T*T*T)                       
      RFH185 = 1.300D+11 * EXP(              - 59.62000*RTR )                   
      RFL185 = 6.200D+14 * EXP(              - 56.10000*RTR )                   
      PR = RFL185*XM(185)/RFH185                                                
      PCOR = PR/(1.0+PR)
      RF(185) = RFH185*PCOR                                                     
      RB(185) = EXP(+( 1.65327D+06)/T/T/T+(-7.17434D+04)/T/T                    
     &          +(-9.62595D+03)/T+( 1.94078D+01)+( 1.23436D-03)*T               
     &          +(-2.42374D-07)*T*T+( 2.54247D-11)*T*T*T)                       
      RB(185) = RB(185)*PCOR                                                    
      RF(186) = 2.110D+12 * EXP(              +  0.48000*RTR )                  
      RB(186) = EXP(+( 5.75654D+06)/T/T/T+(-6.01200D+04)/T/T                    
     &          +(-3.20822D+03)/T+( 3.01833D+01)+(-2.64088D-04)*T               
     &          +( 7.82229D-08)*T*T+(-8.97129D-12)*T*T*T)                       
      RF(187) = 1.060D+20 * EXP( -1.410*ALOGT                )                  
      RB(187) = EXP(+( 5.85897D+06)/T/T/T+(-2.97504D+04)/T/T                    
     &          +(-3.62674D+04)/T+( 4.25990D+01)+(-2.05825D-03)*T               
     &          +( 4.14795D-07)*T*T+(-4.37161D-11)*T*T*T)                       
      RF(188) = 3.900D+12 * EXP(              +  0.24000*RTR )                  
      RB(188) = EXP(+(-5.16409D+06)/T/T/T+( 2.60066D+04)/T/T                    
     &          +(-2.29381D+04)/T+( 2.67906D+01)+( 6.72930D-04)*T               
     &          +(-1.64514D-07)*T*T+( 1.86784D-11)*T*T*T)                       
      RF(189) = 1.320D+14 * EXP(              -  0.36000*RTR )                  
      RB(189) = EXP(+( 9.99851D+05)/T/T/T+(-5.98011D+04)/T/T                    
     &          +(-1.43274D+04)/T+( 2.63728D+01)+( 1.50725D-03)*T               
     &          +(-3.91585D-07)*T*T+( 4.47536D-11)*T*T*T)                       
      RF(190) = 5.000D+13                                                       
      RB(190) = EXP(+(-5.46586D+06)/T/T/T+( 6.69650D+04)/T/T                    
     &          +(-3.60094D+04)/T+( 3.42172D+01)+(-5.87394D-04)*T               
     &          +( 1.64822D-07)*T*T+(-1.83099D-11)*T*T*T)                       
      RF(191) = 3.200D+13 * EXP(              -  0.33000*RTR )                  
      RB(191) = EXP(+( 3.76689D+06)/T/T/T+(-3.58641D+04)/T/T                    
     &          +(-1.23377D+04)/T+( 3.23457D+01)+( 1.12749D-04)*T               
     &          +(-6.12622D-09)*T*T+(-1.29322D-12)*T*T*T)                       
      RF(192) = 2.000D+13                                                       
      RB(192) = EXP(+(-1.18344D+07)/T/T/T+( 1.59635D+05)/T/T                    
     &          +(-9.40222D+03)/T+( 3.55399D+01)+(-8.71194D-04)*T               
     &          +( 1.44534D-07)*T*T+(-1.42753D-11)*T*T*T)                       
      RF(193) = 2.000D+09 * EXP(  1.200*ALOGT                )                  
      RB(193) = EXP(+(-1.56451D+07)/T/T/T+( 2.06643D+05)/T/T                    
     &          +(-2.10444D+04)/T+( 3.28711D+01)+( 9.00272D-04)*T               
     &          +(-1.95959D-07)*T*T+( 2.02663D-11)*T*T*T)                       
      RF(194) = 4.610D+05 * EXP(  2.000*ALOGT -  6.50000*RTR )                  
      RB(194) = EXP(+(-2.13362D+07)/T/T/T+( 2.72161D+05)/T/T                    
     &          +(-5.03784D+03)/T+( 2.76020D+01)+( 1.56913D-03)*T               
     &          +(-3.97536D-07)*T*T+( 4.34798D-11)*T*T*T)                       
      RF(197) = 2.000D+13 * EXP(              - 13.85000*RTR )                  
      RB(197) = EXP(+(-1.82189D+06)/T/T/T+( 3.61284D+04)/T/T                    
     &          +(-8.26539D+03)/T+( 3.34866D+01)+(-6.95119D-04)*T               
     &          +( 1.45368D-07)*T*T+(-1.68269D-11)*T*T*T)                       
      RF(199) = 4.160D+14 * EXP( -0.450*ALOGT                )                  
      RB(199) = EXP(+( 8.78682D+06)/T/T/T+(-5.75663D+04)/T/T                    
     &          +(-1.77803D+04)/T+( 3.85027D+01)+(-1.81876D-03)*T               
     &          +( 4.15415D-07)*T*T+(-4.45463D-11)*T*T*T)                       
      RF(200) = 7.000D+12                                                       
      RB(200) = EXP(+( 1.06750D+07)/T/T/T+(-1.27101D+05)/T/T                    
     &          +(-4.87880D+03)/T+( 2.76427D+01)+( 4.13112D-04)*T               
     &          +(-7.03160D-08)*T*T+( 6.37542D-12)*T*T*T)                       
      RF(201) = 4.600D+13                                                       
      RB(201) = EXP(+(-1.15945D+06)/T/T/T+( 3.25343D+04)/T/T                    
     &          +(-1.42810D+04)/T+( 3.44386D+01)+(-4.58082D-04)*T               
     &          +( 7.42184D-08)*T*T+(-7.89992D-12)*T*T*T)                       
      RF(202) = 4.000D+13 * EXP(              -  3.65000*RTR )                  
      RB(202) = EXP(+( 1.21101D+07)/T/T/T+(-1.42816D+05)/T/T                    
     &          +(-7.63986D+03)/T+( 3.01161D+01)+( 4.23569D-04)*T               
     &          +(-5.68746D-08)*T*T+( 3.24333D-12)*T*T*T)                       
      RF(203) = 9.000D+07 * EXP(  1.500*ALOGT +  0.46000*RTR )                  
      RB(203) = EXP(+(-9.65178D+06)/T/T/T+( 1.29441D+05)/T/T                    
     &          +(-1.46360D+04)/T+( 2.93554D+01)+( 1.45199D-03)*T               
     &          +(-2.93957D-07)*T*T+( 2.95548D-11)*T*T*T)                       
      RF(204) = 3.300D+08                                                       
      RB(204) = EXP(+( 8.71180D+06)/T/T/T+(-1.35372D+05)/T/T                    
     &          +(-3.44330D+03)/T+( 1.98045D+01)+( 4.00280D-04)*T               
     &          +(-5.03555D-08)*T*T+( 4.26252D-12)*T*T*T)                       
      RF(205) = 1.300D+14 * EXP( -0.110*ALOGT -  4.98000*RTR )                  
      RB(205) = EXP(+( 9.57341D+06)/T/T/T+(-1.46281D+05)/T/T                    
     &          +(-5.87907D+03)/T+( 3.19411D+01)+( 3.11951D-04)*T               
     &          +(-3.30306D-08)*T*T+( 2.52013D-12)*T*T*T)                       
      RF(206) = 5.000D+12                                                       
      RB(206) = EXP(+( 3.93380D+06)/T/T/T+(-5.09850D+04)/T/T                    
     &          +(-2.82489D+04)/T+( 2.97532D+01)+( 5.72670D-04)*T               
     &          +(-1.62929D-07)*T*T+( 1.79273D-11)*T*T*T)                       
      RF(207) = 2.500D+13                                                       
      RB(207) = EXP(+( 4.52625D+06)/T/T/T+(-8.50984D+04)/T/T                    
     &          +(-5.47575D+04)/T+( 3.09668D+01)+( 9.81512D-04)*T               
     &          +(-2.49220D-07)*T*T+( 2.76343D-11)*T*T*T)                       
      RF(208) = 7.000D+13                                                       
      RB(208) = EXP(+( 1.79650D+06)/T/T/T+(-5.06879D+04)/T/T                    
     &          +(-5.75143D+03)/T+( 3.03542D+01)+( 6.23339D-04)*T               
     &          +(-1.52522D-07)*T*T+( 1.62561D-11)*T*T*T)                       
      RF(209) = 5.000D+13                                                       
      RB(209) = EXP(+( 5.96136D+06)/T/T/T+(-1.00814D+05)/T/T                    
     &          +(-5.56818D+04)/T+( 3.23904D+01)+( 9.91969D-04)*T               
     &          +(-2.35779D-07)*T*T+( 2.45022D-11)*T*T*T)                       
      RF(210) = 2.000D+13                                                       
      RB(210) = EXP(+(-4.05116D+06)/T/T/T+( 2.26925D+04)/T/T                    
     &          +(-6.37882D+04)/T+( 3.35275D+01)+( 8.15893D-04)*T               
     &          +(-2.36612D-07)*T*T+( 2.70539D-11)*T*T*T)                       
      RF(212) = 8.950D+19 * EXP( -1.320*ALOGT -  0.74000*RTR )                  
      RB(212) = EXP(+(-2.14683D+05)/T/T/T+( 1.20434D+04)/T/T                    
     &          +(-2.42363D+04)/T+( 3.91392D+01)+(-7.62528D-04)*T               
     &          +(-1.12534D-08)*T*T+( 6.49763D-12)*T*T*T)                       
      RF(213) = 2.500D+13                                                       
      RB(213) = EXP(+( 6.36855D+06)/T/T/T+(-9.26699D+04)/T/T                    
     &          +(-2.66072D+04)/T+( 2.86109D+01)+( 2.83800D-04)*T               
     &          +( 2.02877D-08)*T*T+(-4.03453D-12)*T*T*T)                       
      RF(214) = 4.500D+11 * EXP(  0.720*ALOGT -  0.66000*RTR )                  
      RB(214) = EXP(+( 2.16398D+06)/T/T/T+(-3.69853D+04)/T/T                    
     &          +(-2.83235D+04)/T+( 3.02154D+01)+( 8.72417D-04)*T               
     &          +(-7.96704D-08)*T*T+( 4.23814D-12)*T*T*T)                       
      RF(215) = 1.300D+07 * EXP(  1.900*ALOGT +  0.95000*RTR )                  
      RB(215) = EXP(+(-1.70913D+07)/T/T/T+( 2.03538D+05)/T/T                    
     &          +(-3.63733D+04)/T+( 2.98332D+01)+( 1.64388D-03)*T               
     &          +(-2.66353D-07)*T*T+( 2.54809D-11)*T*T*T)                       
      RF(216) = 1.000D+13 * EXP(              - 13.00000*RTR )                  
      RB(216) = EXP(+( 5.77610D+06)/T/T/T+(-5.85565D+04)/T/T                    
     &          +(-6.64036D+03)/T+( 2.80903D+01)+(-1.25042D-04)*T               
     &          +( 1.06579D-07)*T*T+(-1.37416D-11)*T*T*T)                       
      RF(217) = 7.700D+13                                                       
      RB(217) = EXP(+(-1.15195D+05)/T/T/T+( 3.70437D+03)/T/T                    
     &          +(-3.91870D+04)/T+( 3.34825D+01)+( 1.28335D-04)*T               
     &          +(-6.09734D-08)*T*T+( 1.17563D-11)*T*T*T)                       
      RF(218) = 4.000D+13                                                       
      RB(218) = EXP(+( 3.58094D+06)/T/T/T+( 1.61790D+03)/T/T                    
     &          +(-1.56856D+04)/T+( 3.72638D+01)+(-1.77215D-03)*T               
     &          +( 4.19799D-07)*T*T+(-4.20533D-11)*T*T*T)                       
      RF(219) = 8.000D+12 * EXP(              -  7.46000*RTR )                  
      RB(219) = EXP(+( 1.12112D+07)/T/T/T+(-1.04447D+05)/T/T                    
     &          +(-6.18117D+03)/T+( 3.04007D+01)+(-5.92674D-04)*T               
     &          +( 1.09753D-07)*T*T+(-6.78097D-12)*T*T*T)                       
      RF(220) = 6.140D+12 * EXP(              +  0.44000*RTR )                  
      RB(220) = EXP(+( 9.74489D+06)/T/T/T+(-8.41898D+04)/T/T                    
     &          +(-6.55149D+03)/T+( 3.14501D+01)+(-9.37827D-04)*T               
     &          +( 1.92728D-07)*T*T+(-1.59781D-11)*T*T*T)                       
      RF(221) = 2.100D+13 * EXP(              -  4.71000*RTR )                  
      RB(221) = EXP(+( 1.19870D+06)/T/T/T+( 1.90598D+04)/T/T                    
     &          +(-1.29037D+04)/T+( 3.34191D+01)+(-7.68749D-04)*T               
     &          +( 1.08919D-07)*T*T+(-4.22936D-12)*T*T*T)                       
      RF(222) = 2.350D+13                                                       
      RB(222) = EXP(+(-1.14938D+07)/T/T/T+( 8.92000D+04)/T/T                    
     &          +(-4.82636D+04)/T+( 2.85074D+01)+( 1.21080D-03)*T               
     &          +(-2.96383D-07)*T*T+( 3.36609D-11)*T*T*T)                       
      RF(223) = 5.400D+13                                                       
      RB(223) = EXP(+(-6.02791D+06)/T/T/T+( 2.22350D+04)/T/T                    
     &          +(-1.22542D+04)/T+( 2.66653D+01)+( 1.79819D-03)*T               
     &          +(-4.61205D-07)*T*T+( 5.19708D-11)*T*T*T)                       
      RF(224) = 2.500D+12                                                       
      RB(224) = EXP(+(-7.30823D+06)/T/T/T+( 3.89260D+04)/T/T                    
     &          +( 3.05061D+03)/T+( 2.63396D+01)+( 6.29567D-04)*T               
     &          +(-9.75183D-08)*T*T+( 1.02891D-11)*T*T*T)                       
      RF(227) = 8.800D+16 * EXP( -0.500*ALOGT - 48.00000*RTR )                  
      RB(227) = EXP(+( 4.40585D+06)/T/T/T+(-9.77710D+04)/T/T                    
     &          +( 3.97767D+03)/T+( 3.12560D+01)+( 9.17754D-04)*T               
     &          +(-2.03158D-07)*T*T+( 2.25178D-11)*T*T*T)                       
      RF(231) = 1.107D+04 * EXP(  2.640*ALOGT -  4.98000*RTR )                  
      RB(231) = EXP(+(-1.97317D+07)/T/T/T+( 2.60075D+05)/T/T                    
     &          +(-8.41981D+03)/T+( 2.97173D+01)+( 1.10606D-03)*T               
     &          +(-1.18359D-07)*T*T+( 7.12559D-12)*T*T*T)                       
      RF(232) = 2.767D+03 * EXP(  2.640*ALOGT -  4.98000*RTR )                  
      RB(232) = EXP(+(-2.57596D+07)/T/T/T+( 2.82310D+05)/T/T                    
     &          +(-2.06740D+04)/T+( 2.33761D+01)+( 2.90425D-03)*T               
     &          +(-5.79565D-07)*T*T+( 5.90964D-11)*T*T*T)                       
      RF(233) = 2.134D+09 * EXP(  1.580*ALOGT - 26.60000*RTR )                  
      RB(233) = EXP(+(-1.50098D+07)/T/T/T+( 1.53340D+05)/T/T                    
     &          +(-2.93679D+03)/T+( 2.87413D+01)+( 2.02703D-03)*T               
     &          +(-3.71210D-07)*T*T+( 3.23885D-11)*T*T*T)                       
      RF(234) = 1.100D+06 * EXP(  2.030*ALOGT - 13.37000*RTR )                  
      RB(234) = EXP(+(-2.34782D+07)/T/T/T+( 2.90353D+05)/T/T                    
     &          +(-4.62820D+03)/T+( 3.22297D+01)+( 4.27065D-04)*T               
     &          +(-4.66065D-10)*T*T+(-1.12570D-12)*T*T*T)                       
      RF(235) = 4.400D+03 * EXP(  2.260*ALOGT -  6.40000*RTR )                  
      RB(235) = EXP(+(-2.88055D+07)/T/T/T+( 3.55222D+05)/T/T                    
     &          +(-1.42509D+04)/T+( 2.89048D+01)+( 1.30570D-04)*T               
     &          +( 7.93026D-08)*T*T+(-1.07104D-11)*T*T*T)                       
      RF(237) = 1.400D+26 * EXP( -3.400*ALOGT -  1.90000*RTR )                  
      RB(237) = EXP(+( 8.34664D+04)/T/T/T+(-1.75641D+04)/T/T                    
     &          +(-1.25066D+04)/T+( 3.91120D+01)+(-3.47261D-03)*T               
     &          +( 7.10790D-07)*T*T+(-7.50371D-11)*T*T*T)                       
      RF(240) = 2.857D+08 * EXP(  1.100*ALOGT - 20.40000*RTR )                  
      RB(240) = EXP(+(-3.70446D+06)/T/T/T+( 9.88473D+04)/T/T                    
     &          +(-9.88353D+03)/T+( 2.99907D+01)+(-1.15964D-05)*T               
     &          +( 5.22099D-08)*T*T+(-7.33633D-12)*T*T*T)                       
      RFH241 = 3.100D+12 * EXP(  0.150*ALOGT                )                   
      RFL241 = 1.300D+25 * EXP( -3.160*ALOGT -  0.74000*RTR )                   
      PR = RFL241*XM(241)/RFH241                                                
      PRLOG = DLOG10(PR)
      F4 = 0.6670
      F5 =   235.00
      F6 =  2117.00
      F7 =  4536.00
      FC = (1.0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(FC)
      CPRLOG = PRLOG - (0.4 + 0.67*FCLOG)
      X = CPRLOG/(0.75-1.27*FCLOG-0.14*CPRLOG)
      F = 10.0**( FCLOG/(1.0+X*X) )
      PCOR = PR*F/(1.0+PR)
      RF(241) = RFH241*PCOR                                                     
      RB(241) = EXP(+( 4.31061D+06)/T/T/T+( 2.47258D+04)/T/T                    
     &          +(-1.64590D+04)/T+( 3.52285D+01)+(-1.91889D-03)*T               
     &          +( 4.61926D-07)*T*T+(-5.00167D-11)*T*T*T)                       
      RB(241) = RB(241)*PCOR                                                    
      RF(246) = 5.000D+13                                                       
      RB(246) = EXP(+( 5.30969D+06)/T/T/T+(-2.44984D+04)/T/T                    
     &          +(-3.66742D+04)/T+( 3.57140D+01)+(-6.39014D-04)*T               
     &          +( 1.48328D-07)*T*T+(-1.52209D-11)*T*T*T)                       
      RF(247) = 2.000D+13                                                       
      RB(247) = EXP(+( 6.25683D+06)/T/T/T+(-2.62247D+04)/T/T                    
     &          +(-4.09018D+04)/T+( 3.72676D+01)+(-1.65287D-03)*T               
     &          +( 4.45767D-07)*T*T+(-4.99128D-11)*T*T*T)                       
      RF(248) = 3.000D+13                                                       
      RB(248) = EXP(+(-6.58612D+06)/T/T/T+( 9.73758D+04)/T/T                    
     &          +(-2.13979D+04)/T+( 3.39138D+01)+(-6.43304D-04)*T               
     &          +( 1.70590D-07)*T*T+(-1.68582D-11)*T*T*T)                       
      RF(249) = 3.100D+17 * EXP( -1.380*ALOGT -  1.27000*RTR )                  
      RB(249) = EXP(+( 9.04250D+06)/T/T/T+(-9.33300D+04)/T/T                    
     &          +(-4.62700D+04)/T+( 3.81757D+01)+(-2.87825D-03)*T               
     &          +( 6.56308D-07)*T*T+(-6.99683D-11)*T*T*T)                       
      RF(250) = 2.900D+14 * EXP( -0.690*ALOGT -  0.76000*RTR )                  
      RB(250) = EXP(+( 1.47410D+07)/T/T/T+(-1.56010D+05)/T/T                    
     &          +(-3.68672D+04)/T+( 3.07272D+01)+(-6.39974D-04)*T               
     &          +( 1.12383D-07)*T*T+(-1.25300D-11)*T*T*T)                       
      RF(251) = 3.800D+13 * EXP( -0.360*ALOGT -  0.58000*RTR )                  
      RB(251) = EXP(+( 9.12667D+06)/T/T/T+(-6.04661D+04)/T/T                    
     &          +(-1.16237D+04)/T+( 3.57747D+01)+(-2.22801D-03)*T               
     &          +( 5.37507D-07)*T*T+(-5.84946D-11)*T*T*T)                       
      RF(252) = 3.100D+17 * EXP( -1.380*ALOGT -  1.27000*RTR )                  
      RB(252) = EXP(+( 1.74248D+07)/T/T/T+(-1.76631D+05)/T/T                    
     &          +(-5.05028D+04)/T+( 3.72030D+01)+(-2.82723D-03)*T               
     &          +( 6.52016D-07)*T*T+(-7.04135D-11)*T*T*T)                       
      RF(254) = 3.800D+13 * EXP( -0.360*ALOGT -  0.58000*RTR )                  
      RB(254) = EXP(+( 1.75090D+07)/T/T/T+(-1.43767D+05)/T/T                    
     &          +(-1.58564D+04)/T+( 3.48020D+01)+(-2.17698D-03)*T               
     &          +( 5.33215D-07)*T*T+(-5.89399D-11)*T*T*T)                       
      RF(255) = 9.600D+13 * EXP(              - 28.80000*RTR )                  
      RB(255) = EXP(+( 1.39634D+07)/T/T/T+(-1.43175D+05)/T/T                    
     &          +(-5.54582D+04)/T+( 3.30916D+01)+( 4.14605D-04)*T               
     &          +(-8.98745D-08)*T*T+( 5.96125D-12)*T*T*T)                       
      RF(256) = 1.000D+12 * EXP(              - 21.75000*RTR )                  
      RB(256) = EXP(+( 1.77989D+05)/T/T/T+( 1.83639D+04)/T/T                    
     &          +(-5.28756D+03)/T+( 2.78234D+01)+(-7.43415D-04)*T               
     &          +( 2.71674D-07)*T*T+(-3.80114D-11)*T*T*T)                       
      RF(257) = 2.200D+13                                                       
      RB(257) = EXP(+(-1.07225D+07)/T/T/T+( 5.31245D+04)/T/T                    
     &          +(-7.28023D+04)/T+( 2.96353D+01)+( 1.59726D-03)*T               
     &          +(-3.36167D-07)*T*T+( 3.61408D-11)*T*T*T)                       
      RF(259) = 1.200D+13                                                       
      RB(259) = EXP(+(-1.37054D+07)/T/T/T+( 8.88308D+04)/T/T                    
     &          +(-2.08842D+04)/T+( 2.57664D+01)+( 1.54067D-03)*T               
     &          +(-3.57642D-07)*T*T+( 4.14610D-11)*T*T*T)                       
      RF(260) = 1.200D+13                                                       
      RB(260) = EXP(+(-1.98693D+07)/T/T/T+( 1.74639D+05)/T/T                    
     &          +(-2.97969D+04)/T+( 2.97061D+01)+( 7.06348D-04)*T               
     &          +(-1.30571D-07)*T*T+( 1.53858D-11)*T*T*T)                       
      RF(261) = 1.000D+14                                                       
      RB(261) = EXP(+(-1.36977D+07)/T/T/T+( 1.03509D+05)/T/T                    
     &          +(-3.46997D+04)/T+( 2.87733D+01)+( 2.06746D-03)*T               
     &          +(-5.39796D-07)*T*T+( 6.21440D-11)*T*T*T)                       
      RF(262) = 9.800D+07 * EXP(  1.410*ALOGT -  8.50000*RTR )                  
      RB(262) = EXP(+(-3.16029D+06)/T/T/T+( 5.99194D+04)/T/T                    
     &          +(-2.50062D+04)/T+( 2.80130D+01)+( 2.01679D-03)*T               
     &          +(-4.52515D-07)*T*T+( 4.71923D-11)*T*T*T)                       
      RF(264) = 2.200D+06 * EXP(  2.110*ALOGT - 11.40000*RTR )                  
      RB(264) = EXP(+(-4.47705D+06)/T/T/T+( 7.64111D+04)/T/T                    
     &          +(-1.72510D+03)/T+( 2.62469D+01)+( 2.36468D-03)*T               
     &          +(-4.70136D-07)*T*T+( 4.52391D-11)*T*T*T)                       
      RF(265) = 2.250D+07 * EXP(  1.700*ALOGT -  3.80000*RTR )                  
      RB(265) = EXP(+(-1.79684D+07)/T/T/T+( 1.85088D+05)/T/T                    
     &          +(-5.01419D+03)/T+( 2.27660D+01)+( 3.42053D-03)*T               
     &          +(-7.96450D-07)*T*T+( 8.43401D-11)*T*T*T)                       
      RF(267) = 4.650D+12 * EXP(              -  6.85000*RTR )                  
      RB(267) = EXP(+( 3.47294D+06)/T/T/T+(-2.50404D+04)/T/T                    
     &          +(-7.11856D+03)/T+( 2.92598D+01)+( 5.04732D-04)*T               
     &          +(-1.25204D-07)*T*T+( 1.12364D-11)*T*T*T)                       
      RF(268) = 1.550D+12 * EXP(              -  6.85000*RTR )                  
      RB(268) = EXP(+(-2.79088D+06)/T/T/T+( 4.71944D+04)/T/T                    
     &          +(-1.83966D+04)/T+( 3.00369D+01)+( 4.71447D-04)*T               
     &          +(-1.60125D-07)*T*T+( 1.84826D-11)*T*T*T)                       
      RF(270) = 2.100D+15 * EXP( -0.690*ALOGT -  2.85000*RTR )                  
      RB(270) = EXP(+(-2.66903D+06)/T/T/T+(-1.38787D+02)/T/T                    
     &          +(-3.59440D+04)/T+( 3.09169D+01)+(-3.85255D-04)*T               
     &          +( 6.68268D-08)*T*T+(-6.24650D-12)*T*T*T)                       
      RF(271) = 2.700D+11 * EXP(  0.180*ALOGT -  2.12000*RTR )                  
      RB(271) = EXP(+( 1.61957D+06)/T/T/T+(-4.49685D+04)/T/T                    
     &          +(-2.65455D+04)/T+( 2.27067D+01)+( 1.99756D-03)*T               
     &          +(-5.05449D-07)*T*T+( 5.40429D-11)*T*T*T)                       
      RF(273) = 2.000D+07 * EXP(  2.000*ALOGT -  2.00000*RTR )                  
      RB(273) = EXP(+(-1.91915D+07)/T/T/T+( 2.40395D+05)/T/T                    
     &          +(-1.52670D+04)/T+( 3.10327D+01)+( 1.12481D-03)*T               
     &          +(-1.99005D-07)*T*T+( 1.84519D-11)*T*T*T)                       
      RF(274) = 2.350D+13                                                       
      RB(274) = EXP(+(-7.81095D+06)/T/T/T+( 8.89089D+04)/T/T                    
     &          +(-2.53906D+04)/T+( 3.33860D+01)+(-2.17385D-04)*T               
     &          +( 3.74679D-08)*T*T+(-1.69619D-12)*T*T*T)                       
      RF(275) = 6.100D+14 * EXP( -0.310*ALOGT -  0.29000*RTR )                  
      RB(275) = EXP(+(-5.19145D+06)/T/T/T+( 7.47356D+04)/T/T                    
     &          +(-1.90527D+04)/T+( 3.42864D+01)+(-1.68203D-03)*T               
     &          +( 5.04889D-07)*T*T+(-6.30705D-11)*T*T*T)                       

      RETURN
      END

      SUBROUTINE NETRATE( W, RF, RB, XM, XH2, XH, XO, XO2, XOH, XH2O,           
     &                    XHO2, XH2O2, XC, XCH, XCH2, XCH2S, XCH3, XCH4,        
     &                    XCO, XCO2, XHCO, XCH2O, XCH2OH, XCH3O, XCH3OH,        
     &                    XC2H, XC2H2, XC2H3, XC2H4, XC2H5, XC2H6,              
     &                    XHCCO, XCH2CO, XHCCOH, XN, XNH, XNH2, XNNH,           
     &                    XNO, XNO2, XN2O, XHNO, XCN, XHCN, XH2CN,              
     &                    XHCNN, XHCNO, XHOCN, XHNCO, XNCO, XN2 )               

C*****double precision
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
C       IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision

      DIMENSION W(*), RF(*), RB(*), XM(*)

C   NET PRODUCTION RATES FOR SKELETAL MECHANSIM

      W(3) = RF(3)*XH2*XO - RB(3)*XH*XOH                                        
      W(4) = RF(4)*XO*XHO2 - RB(4)*XO2*XOH                                      
      W(5) = RF(5)*XO*XH2O2 - RB(5)*XOH*XHO2                                    
      W(7) = RF(7)*XO*XCH2 - RB(7)*XH*XHCO                                      
      W(10) = RF(10)*XO*XCH3 - RB(10)*XH*XCH2O                                  
      W(11) = RF(11)*XO*XCH4 - RB(11)*XOH*XCH3                                  
      W(15) = RF(15)*XO*XCH2O - RB(15)*XOH*XHCO                                 
      W(18) = RF(18)*XO*XCH3OH - RB(18)*XOH*XCH2OH                              
      W(19) = RF(19)*XO*XCH3OH - RB(19)*XOH*XCH3O                               
      W(21) = RF(21)*XO*XC2H2 - RB(21)*XH*XHCCO                                 
      W(22) = RF(22)*XO*XC2H2 - RB(22)*XOH*XC2H                                 
      W(23) = RF(23)*XO*XC2H2 - RB(23)*XCH2*XCO                                 
      W(24) = RF(24)*XO*XC2H3 - RB(24)*XH*XCH2CO                                
      W(25) = RF(25)*XO*XC2H4 - RB(25)*XCH3*XHCO                                
      W(26) = RF(26)*XO*XC2H5 - RB(26)*XCH3*XCH2O                               
      W(27) = RF(27)*XO*XC2H6 - RB(27)*XOH*XC2H5                                
      W(28) = RF(28)*XO*XHCCO - RB(28)*XH*XCO*XCO                               
      W(29) = RF(29)*XO*XCH2CO - RB(29)*XOH*XHCCO                               
      W(30) = RF(30)*XO*XCH2CO - RB(30)*XCH2*XCO2                               
      W(33) = RF(33)*XH*XO2*XM(33) - RB(33)*XHO2*XM(33)                         
      W(35) = RF(35)*XH*XO2*XH2O - RB(35)*XHO2*XH2O                             
      W(36) = RF(36)*XH*XO2*XN2 - RB(36)*XHO2*XN2                               
      W(38) = RF(38)*XH*XO2 - RB(38)*XO*XOH                                     
      W(43) = RF(43)*XH*XOH*XM(43) - RB(43)*XH2O*XM(43)                         
      W(45) = RF(45)*XH*XHO2 - RB(45)*XH2*XO2                                   
      W(46) = RF(46)*XH*XHO2 - RB(46)*XOH*XOH                                   
      W(47) = RF(47)*XH*XH2O2 - RB(47)*XH2*XHO2                                 
      W(49) = RF(49)*XH*XCH - RB(49)*XH2*XC                                     
      W(52) = RF(52)*XH*XCH3 - RB(52)*XCH4                                      
      W(53) = RF(53)*XH*XCH4 - RB(53)*XH2*XCH3                                  
      W(55) = RF(55)*XH*XHCO - RB(55)*XH2*XCO                                   
      W(56) = RF(56)*XH*XCH2O - RB(56)*XCH2OH                                   
      W(57) = RF(57)*XH*XCH2O - RB(57)*XCH3O                                    
      W(58) = RF(58)*XH*XCH2O - RB(58)*XH2*XHCO                                 
      W(61) = RF(61)*XH*XCH2OH - RB(61)*XOH*XCH3                                
      W(62) = RF(62)*XH*XCH2OH - RB(62)*XH2O*XCH2S                              
      W(66) = RF(66)*XH*XCH3O - RB(66)*XOH*XCH3                                 
      W(68) = RF(68)*XH*XCH3OH - RB(68)*XH2*XCH2OH                              
      W(69) = RF(69)*XH*XCH3OH - RB(69)*XH2*XCH3O                               
      W(71) = RF(71)*XH*XC2H2 - RB(71)*XC2H3                                    
      W(73) = RF(73)*XH*XC2H3 - RB(73)*XH2*XC2H2                                
      W(74) = RF(74)*XH*XC2H4 - RB(74)*XC2H5                                    
      W(75) = RF(75)*XH*XC2H4 - RB(75)*XH2*XC2H3                                
      W(78) = RF(78)*XH*XC2H6 - RB(78)*XH2*XC2H5                                
      W(79) = RF(79)*XH*XHCCO - RB(79)*XCH2S*XCO                                
      W(80) = RF(80)*XH*XCH2CO - RB(80)*XH2*XHCCO                               
      W(81) = RF(81)*XH*XCH2CO - RB(81)*XCH3*XCO                                
      W(82) = RF(82)*XHCCOH*XH - RB(82)*XCH2CO*XH                               
      W(83) = RF(83)*XH2*XCO - RB(83)*XCH2O                                     
      W(84) = RF(84)*XH2*XOH - RB(84)*XH*XH2O                                   
      W(85) = RF(85)*XOH*XOH - RB(85)*XH2O2                                     
      W(86) = RF(86)*XOH*XOH - RB(86)*XO*XH2O                                   
      W(87) = RF(87)*XOH*XHO2 - RB(87)*XO2*XH2O                                 
      W(88) = RF(88)*XOH*XH2O2 - RB(88)*XH2O*XHO2                               
      W(89) = RF(89)*XOH*XH2O2 - RB(89)*XH2O*XHO2                               
      W(90) = RF(90)*XOH*XC - RB(90)*XH*XCO                                     
      W(92) = RF(92)*XOH*XCH2 - RB(92)*XH*XCH2O                                 
      W(93) = RF(93)*XOH*XCH2 - RB(93)*XH2O*XCH                                 
      W(95) = RF(95)*XOH*XCH3 - RB(95)*XCH3OH                                   
      W(96) = RF(96)*XOH*XCH3 - RB(96)*XH2O*XCH2                                
      W(97) = RF(97)*XOH*XCH3 - RB(97)*XH2O*XCH2S                               
      W(98) = RF(98)*XOH*XCH4 - RB(98)*XH2O*XCH3                                
      W(99) = RF(99)*XOH*XCO - RB(99)*XH*XCO2                                   
      W(100) = RF(100)*XOH*XHCO - RB(100)*XH2O*XCO                              
      W(101) = RF(101)*XOH*XCH2O - RB(101)*XH2O*XHCO                            
      W(104) = RF(104)*XOH*XCH3OH - RB(104)*XH2O*XCH2OH                         
      W(105) = RF(105)*XOH*XCH3OH - RB(105)*XH2O*XCH3O                          
      W(106) = RF(106)*XOH*XC2H - RB(106)*XH*XHCCO                              
      W(107) = RF(107)*XOH*XC2H2 - RB(107)*XH*XCH2CO                            
      W(108) = RF(108)*XOH*XC2H2 - RB(108)*XH*XHCCOH                            
      W(109) = RF(109)*XOH*XC2H2 - RB(109)*XH2O*XC2H                            
      W(111) = RF(111)*XOH*XC2H3 - RB(111)*XH2O*XC2H2                           
      W(112) = RF(112)*XOH*XC2H4 - RB(112)*XH2O*XC2H3                           
      W(113) = RF(113)*XOH*XC2H6 - RB(113)*XH2O*XC2H5                           
      W(114) = RF(114)*XOH*XCH2CO - RB(114)*XH2O*XHCCO                          
      W(122) = RF(122)*XO2*XC - RB(122)*XO*XCO                                  
      W(125) = RF(125)*XO2*XCH - RB(125)*XO*XHCO                                
      W(126) = RF(126)*XH2*XCH - RB(126)*XH*XCH2                                
      W(127) = RF(127)*XH2O*XCH - RB(127)*XH*XCH2O                              
      W(130) = RF(130)*XCH*XCH4 - RB(130)*XH*XC2H4                              
      W(132) = RF(132)*XCH*XCO2 - RB(132)*XCO*XHCO                              
      W(133) = RF(133)*XCH*XCH2O - RB(133)*XH*XCH2CO                            
      W(135) = RF(135)*XO2*XCH2 - RB(135)*XOH*XHCO                              
      W(137) = RF(137)*XCH2*XCH2 - RB(137)*XH2*XC2H2                            
      W(138) = RF(138)*XCH2*XCH3 - RB(138)*XH*XC2H4                             
      W(140) = RF(140)*XCH2*XCO - RB(140)*XCH2CO                                
      W(142) = RF(142)*XCH2S*XN2 - RB(142)*XCH2*XN2                             
      W(144) = RF(144)*XO2*XCH2S - RB(144)*XH*XOH*XCO                           
      W(145) = RF(145)*XO2*XCH2S - RB(145)*XH2O*XCO                             
      W(146) = RF(146)*XH2*XCH2S - RB(146)*XH*XCH3                              
      W(147) = RF(147)*XH2O*XCH2S - RB(147)*XCH3OH                              
      W(148) = RF(148)*XCH2S*XH2O - RB(148)*XCH2*XH2O                           
      W(153) = RF(153)*XCH2S*XCO2 - RB(153)*XCO*XCH2O                           
      W(155) = RF(155)*XO2*XCH3 - RB(155)*XO*XCH3O                              
      W(158) = RF(158)*XCH3*XCH3 - RB(158)*XC2H6                                
      W(159) = RF(159)*XCH3*XCH3 - RB(159)*XH*XC2H5                             
      W(160) = RF(160)*XCH3*XHCO - RB(160)*XCH4*XCO                             
      W(161) = RF(161)*XCH3*XCH2O - RB(161)*XCH4*XHCO                           
      W(166) = RF(166)*XHCO*XH2O - RB(166)*XH*XCO*XH2O                          
      W(167) = RF(167)*XHCO*XM(167) - RB(167)*XH*XCO*XM(167)                    
      W(168) = RF(168)*XO2*XHCO - RB(168)*XHO2*XCO                              
      W(169) = RF(169)*XO2*XCH2OH - RB(169)*XHO2*XCH2O                          
      W(170) = RF(170)*XO2*XCH3O - RB(170)*XHO2*XCH2O                           
      W(171) = RF(171)*XO2*XC2H - RB(171)*XCO*XHCO                              
      W(173) = RF(173)*XO2*XC2H3 - RB(173)*XHCO*XCH2O                           
      W(176) = RF(176)*XO2*XHCCO - RB(176)*XOH*XCO*XCO                          
      W(178) = RF(178)*XN*XNO - RB(178)*XO*XN2                                  
      W(179) = RF(179)*XO2*XN - RB(179)*XO*XNO                                  
      W(180) = RF(180)*XOH*XN - RB(180)*XH*XNO                                  
      W(181) = RF(181)*XO*XN2O - RB(181)*XO2*XN2                                
      W(182) = RF(182)*XO*XN2O - RB(182)*XNO*XNO                                
      W(183) = RF(183)*XH*XN2O - RB(183)*XOH*XN2                                
      W(185) = RF(185)*XN2O - RB(185)*XO*XN2                                    
      W(186) = RF(186)*XHO2*XNO - RB(186)*XOH*XNO2                              
      W(187) = RF(187)*XO*XNO*XM(187) - RB(187)*XNO2*XM(187)                    
      W(188) = RF(188)*XO*XNO2 - RB(188)*XO2*XNO                                
      W(189) = RF(189)*XH*XNO2 - RB(189)*XOH*XNO                                
      W(190) = RF(190)*XO*XNH - RB(190)*XH*XNO                                  
      W(191) = RF(191)*XH*XNH - RB(191)*XH2*XN                                  
      W(192) = RF(192)*XOH*XNH - RB(192)*XH*XHNO                                
      W(193) = RF(193)*XOH*XNH - RB(193)*XH2O*XN                                
      W(194) = RF(194)*XO2*XNH - RB(194)*XO*XHNO                                
      W(197) = RF(197)*XH2O*XNH - RB(197)*XH2*XHNO                              
      W(199) = RF(199)*XNH*XNO - RB(199)*XH*XN2O                                
      W(200) = RF(200)*XO*XNH2 - RB(200)*XOH*XNH                                
      W(201) = RF(201)*XO*XNH2 - RB(201)*XH*XHNO                                
      W(202) = RF(202)*XH*XNH2 - RB(202)*XH2*XNH                                
      W(203) = RF(203)*XOH*XNH2 - RB(203)*XH2O*XNH                              
      W(204) = RF(204)*XNNH - RB(204)*XH*XN2                                    
      W(205) = RF(205)*XNNH*XM(205) - RB(205)*XH*XN2*XM(205)                    
      W(206) = RF(206)*XO2*XNNH - RB(206)*XHO2*XN2                              
      W(207) = RF(207)*XO*XNNH - RB(207)*XOH*XN2                                
      W(208) = RF(208)*XO*XNNH - RB(208)*XNH*XNO                                
      W(209) = RF(209)*XH*XNNH - RB(209)*XH2*XN2                                
      W(210) = RF(210)*XOH*XNNH - RB(210)*XH2O*XN2                              
      W(212) = RF(212)*XH*XNO*XM(212) - RB(212)*XHNO*XM(212)                    
      W(213) = RF(213)*XO*XHNO - RB(213)*XOH*XNO                                
      W(214) = RF(214)*XH*XHNO - RB(214)*XH2*XNO                                
      W(215) = RF(215)*XOH*XHNO - RB(215)*XH2O*XNO                              
      W(216) = RF(216)*XO2*XHNO - RB(216)*XHO2*XNO                              
      W(217) = RF(217)*XO*XCN - RB(217)*XCO*XN                                  
      W(218) = RF(218)*XOH*XCN - RB(218)*XH*XNCO                                
      W(219) = RF(219)*XH2O*XCN - RB(219)*XOH*XHCN                              
      W(220) = RF(220)*XO2*XCN - RB(220)*XO*XNCO                                
      W(221) = RF(221)*XH2*XCN - RB(221)*XH*XHCN                                
      W(222) = RF(222)*XO*XNCO - RB(222)*XCO*XNO                                
      W(223) = RF(223)*XH*XNCO - RB(223)*XCO*XNH                                
      W(224) = RF(224)*XOH*XNCO - RB(224)*XH*XCO*XNO                            
      W(227) = RF(227)*XNCO*XM(227) - RB(227)*XCO*XN*XM(227)                    
      W(231) = RF(231)*XO*XHCN - RB(231)*XH*XNCO                                
      W(232) = RF(232)*XO*XHCN - RB(232)*XCO*XNH                                
      W(233) = RF(233)*XO*XHCN - RB(233)*XOH*XCN                                
      W(234) = RF(234)*XOH*XHCN - RB(234)*XH*XHOCN                              
      W(235) = RF(235)*XOH*XHCN - RB(235)*XH*XHNCO                              
      W(237) = RF(237)*XH*XHCN*XM(237) - RB(237)*XH2CN*XM(237)                  
      W(240) = RF(240)*XCH*XN2 - RB(240)*XN*XHCN                                
      W(241) = RF(241)*XCH*XN2 - RB(241)*XHCNN                                  
      W(246) = RF(246)*XCH*XNO - RB(246)*XO*XHCN                                
      W(247) = RF(247)*XCH*XNO - RB(247)*XH*XNCO                                
      W(248) = RF(248)*XCH*XNO - RB(248)*XHCO*XN                                
      W(249) = RF(249)*XCH2*XNO - RB(249)*XH*XHNCO                              
      W(250) = RF(250)*XCH2*XNO - RB(250)*XOH*XHCN                              
      W(251) = RF(251)*XCH2*XNO - RB(251)*XH*XHCNO                              
      W(252) = RF(252)*XCH2S*XNO - RB(252)*XH*XHNCO                             
      W(254) = RF(254)*XCH2S*XNO - RB(254)*XH*XHCNO                             
      W(255) = RF(255)*XCH3*XNO - RB(255)*XH2O*XHCN                             
      W(256) = RF(256)*XCH3*XNO - RB(256)*XOH*XH2CN                             
      W(257) = RF(257)*XO*XHCNN - RB(257)*XH*XCO*XN2                            
      W(259) = RF(259)*XO2*XHCNN - RB(259)*XO*XHCO*XN2                          
      W(260) = RF(260)*XOH*XHCNN - RB(260)*XH*XHCO*XN2                          
      W(261) = RF(261)*XH*XHCNN - RB(261)*XCH2*XN2                              
      W(262) = RF(262)*XO*XHNCO - RB(262)*XCO2*XNH                              
      W(264) = RF(264)*XO*XHNCO - RB(264)*XOH*XNCO                              
      W(265) = RF(265)*XH*XHNCO - RB(265)*XCO*XNH2                              
      W(267) = RF(267)*XOH*XHNCO - RB(267)*XH2O*XNCO                            
      W(268) = RF(268)*XOH*XHNCO - RB(268)*XCO2*XNH2                            
      W(270) = RF(270)*XHCNO*XH - RB(270)*XHNCO*XH                              
      W(271) = RF(271)*XH*XHCNO - RB(271)*XOH*XHCN                              
      W(273) = RF(273)*XHOCN*XH - RB(273)*XHNCO*XH                              
      W(274) = RF(274)*XHCCO*XNO - RB(274)*XCO*XHCNO                            
      W(275) = RF(275)*XCH3*XN - RB(275)*XH*XH2CN                               

      RETURN
      END

