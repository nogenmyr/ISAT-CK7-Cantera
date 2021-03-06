C-----------------------------------------------------------------------
C ======= 11stepF.h =======
C-----------------------------------------------------------------------

      INTEGER  SCH4, SH
     &    , SH2O, SH2
     &    , SCO, SCO2
     &    , SO2, SC2H2
     &    , SN2, SNO
     &    , SHCN, SHNCO
     &    , SNH3, SN2O
     &    , SNO2, SOH
     &    , SO, SHO2
     &    , SCH, SCHO
     &    , SCH2OH, S3XCH2
     &    , SCH3, SC2H4
     &    , S1XCH2, SCH2O
     &    , SCH3O, SC2H6
     &    , SCH3OH, SHCCO
      INTEGER  SC2H3, SC2H5
     &    , SNH2, SNH
     &    , SHNO, SH2NO
     &    , SN2H, SN
     &    , SNCO, SCN
     &    , SEND
      INTEGER  RI, RII
     &    , RIII, RIIIN
     &    , RIV, RIVN
     &    , RV, RVI
     &    , RVII, RVIII
     &    , RIX, RX
     &    , RXI, R1F
     &    , R1B, R2F
     &    , R2B, R3F
     &    , R3B, R4F
     &    , R4B, R5F
     &    , R5B, R6F
     &    , R6B, R7F
     &    , R7B, R8F
     &    , R8B, R9
      INTEGER  R10, R11
     &    , R12, R13F
     &    , R13B, R20F
     &    , R20B, R21
     &    , R22, R25
     &    , R26, R27F
     &    , R27B, R28F
     &    , R28B, R35F
     &    , R35B, R38
     &    , R39, R40
     &    , R41, R42F
     &    , R42B, R43
     &    , R44F, R44B
     &    , R45, R46
     &    , R48, R53
      INTEGER  R34F, R34B
     &    , R55, R60
     &    , R36F, R36B
     &    , R62, R63
     &    , R69, R70
     &    , R71, R84F
     &    , R84B, R86F
     &    , R86B, R90F
     &    , R90B, R91
     &    , R93, R97
     &    , R105F, R105B
     &    , R106, R109
     &    , R111, R112
     &    , R51F, R51B
     &    , RA125, R129F
      INTEGER  R129B, R131F
     &    , R131B, R134F
     &    , R134B, R58F
     &    , R58B, R146
     &    , R164, R166
     &    , R170F, R170B
     &    , RN1F, RN1B
     &    , RN2F, RN2B
     &    , RN3F, RN3B
     &    , RN4F, RN4B
     &    , RN6F, RN6B
     &    , RN7F, RN7B
     &    , RN9F, RN9B
     &    , RN101F, RN101B
     &    , RN13, RN14F
      INTEGER  RN14B, RN16F
     &    , RN16B, RN17F
     &    , RN17B, RN18F
     &    , RN18B, RN19F
     &    , RN19B, RN20F
     &    , RN20B, RN24XF
     &    , RN24XB, RN24Y
     &    , RN25F, RN25B
     &    , RN27F, RN27B
     &    , RN28F, RN28B
     &    , RN30F, RN30B
     &    , RN102F, RN102B
     &    , RN104F, RN104B
     &    , RN35F, RN35B
     &    , RN36F, RN36B
      INTEGER  RN37F, RN37B
     &    , RN38F, RN38B
     &    , RN39F, RN39B
     &    , RN40F, RN40B
     &    , RN41F, RN41B
     &    , RN45F, RN45B
     &    , RN46YF, RN46YB
     &    , RN47F, RN47B
     &    , RN84F, RN84B
     &    , RN85F, RN85B
     &    , RN87F, RN87B
     &    , RN88F, RN88B
     &    , RN60F, RN60B
     &    , RN61F, RN61B
     &    , RN62F, RN62B
      INTEGER  RN63F, RN63B
     &    , RN65F, RN65B
     &    , RN70F, RN70B
     &    , RN71F, RN71B
     &    , RN72F, RN72B
     &    , RN74F, RN74B
     &    , RN110F, RN110B
     &    , RN111, RN112F
     &    , RN112B, RN113F
     &    , RN113B, RN120F
     &    , RN120B, RN130F
     &    , RN130B, RN131F
     &    , RN131B, RN134F
     &    , RN134B, RN136F
     &    , RN136B
     &    , REND
      INTEGER  MM1, MM0
     &    , MM6, MM3
     &    , MM8, MM7
     &    , MEND
      PARAMETER (  SCH4 = 1, SH = 2
     &    , SH2O = 3, SH2 = 4
     &    , SCO = 5, SCO2 = 6
     &    , SO2 = 7, SC2H2 = 8
     &    , SN2 = 9, SNO = 10
     &    , SHCN = 11, SHNCO = 12
     &    , SNH3 = 13, SN2O = 14
     &    , SNO2 = 15, SOH = 16
     &    , SO = 17, SHO2 = 18
     &    , SCH = 19, SCHO = 20
     &    , SCH2OH = 21, S3XCH2 = 22
     &    , SCH3 = 23, SC2H4 = 24
     &    , S1XCH2 = 25, SCH2O = 26
     &    , SCH3O = 27, SC2H6 = 28
     &    , SCH3OH = 29, SHCCO = 30 )
      PARAMETER (  SC2H3 = 31, SC2H5 = 32
     &    , SNH2 = 33, SNH = 34
     &    , SHNO = 35, SH2NO = 36
     &    , SN2H = 37, SN = 38
     &    , SNCO = 39, SCN = 40
     &    , SEND = 41 )
      PARAMETER (  RI = 1, RII = 2
     &    , RIII = 3, RIIIN = 4
     &    , RIV = 5, RIVN = 6
     &    , RV = 7, RVI = 8
     &    , RVII = 9, RVIII = 10
     &    , RIX = 11, RX = 12
     &    , RXI = 13, R1F = 14
     &    , R1B = 15, R2F = 16
     &    , R2B = 17, R3F = 18
     &    , R3B = 19, R4F = 20
     &    , R4B = 21, R5F = 22
     &    , R5B = 23, R6F = 24
     &    , R6B = 25, R7F = 26
     &    , R7B = 27, R8F = 28
     &    , R8B = 29, R9 = 30 )
      PARAMETER (  R10 = 31, R11 = 32
     &    , R12 = 33, R13F = 34
     &    , R13B = 35, R20F = 36
     &    , R20B = 37, R21 = 38
     &    , R22 = 39, R25 = 40
     &    , R26 = 41, R27F = 42
     &    , R27B = 43, R28F = 44
     &    , R28B = 45, R35F = 46
     &    , R35B = 47, R38 = 48
     &    , R39 = 49, R40 = 50
     &    , R41 = 51, R42F = 52
     &    , R42B = 53, R43 = 54
     &    , R44F = 55, R44B = 56
     &    , R45 = 57, R46 = 58
     &    , R48 = 59, R53 = 60 )
      PARAMETER (  R34F = 61, R34B = 62
     &    , R55 = 63, R60 = 64
     &    , R36F = 65, R36B = 66
     &    , R62 = 67, R63 = 68
     &    , R69 = 69, R70 = 70
     &    , R71 = 71, R84F = 72
     &    , R84B = 73, R86F = 74
     &    , R86B = 75, R90F = 76
     &    , R90B = 77, R91 = 78
     &    , R93 = 79, R97 = 80
     &    , R105F = 81, R105B = 82
     &    , R106 = 83, R109 = 84
     &    , R111 = 85, R112 = 86
     &    , R51F = 87, R51B = 88
     &    , RA125 = 89, R129F = 90 )
      PARAMETER (  R129B = 91, R131F = 92
     &    , R131B = 93, R134F = 94
     &    , R134B = 95, R58F = 96
     &    , R58B = 97, R146 = 98
     &    , R164 = 99, R166 = 100
     &    , R170F = 101, R170B = 102
     &    , RN1F = 103, RN1B = 104
     &    , RN2F = 105, RN2B = 106
     &    , RN3F = 107, RN3B = 108
     &    , RN4F = 109, RN4B = 110
     &    , RN6F = 111, RN6B = 112
     &    , RN7F = 113, RN7B = 114
     &    , RN9F = 115, RN9B = 116
     &    , RN101F = 117, RN101B = 118
     &    , RN13 = 119, RN14F = 120 )
      PARAMETER (  RN14B = 121, RN16F = 122
     &    , RN16B = 123, RN17F = 124
     &    , RN17B = 125, RN18F = 126
     &    , RN18B = 127, RN19F = 128
     &    , RN19B = 129, RN20F = 130
     &    , RN20B = 131, RN24XF = 132
     &    , RN24XB = 133, RN24Y = 134
     &    , RN25F = 135, RN25B = 136
     &    , RN27F = 137, RN27B = 138
     &    , RN28F = 139, RN28B = 140
     &    , RN30F = 141, RN30B = 142
     &    , RN102F = 143, RN102B = 144
     &    , RN104F = 145, RN104B = 146
     &    , RN35F = 147, RN35B = 148
     &    , RN36F = 149, RN36B = 150 )
      PARAMETER (  RN37F = 151, RN37B = 152
     &    , RN38F = 153, RN38B = 154
     &    , RN39F = 155, RN39B = 156
     &    , RN40F = 157, RN40B = 158
     &    , RN41F = 159, RN41B = 160
     &    , RN45F = 161, RN45B = 162
     &    , RN46YF = 163, RN46YB = 164
     &    , RN47F = 165, RN47B = 166
     &    , RN84F = 167, RN84B = 168
     &    , RN85F = 169, RN85B = 170
     &    , RN87F = 171, RN87B = 172
     &    , RN88F = 173, RN88B = 174
     &    , RN60F = 175, RN60B = 176
     &    , RN61F = 177, RN61B = 178
     &    , RN62F = 179, RN62B = 180 )
      PARAMETER (  RN63F = 181, RN63B = 182
     &    , RN65F = 183, RN65B = 184
     &    , RN70F = 185, RN70B = 186
     &    , RN71F = 187, RN71B = 188
     &    , RN72F = 189, RN72B = 190
     &    , RN74F = 191, RN74B = 192
     &    , RN110F = 193, RN110B = 194
     &    , RN111 = 195, RN112F = 196
     &    , RN112B = 197, RN113F = 198
     &    , RN113B = 199, RN120F = 200
     &    , RN120B = 201, RN130F = 202
     &    , RN130B = 203, RN131F = 204
     &    , RN131B = 205, RN134F = 206
     &    , RN134B = 207, RN136F = 208
     &    , RN136B = 209
     &    , REND = 210 )
      PARAMETER (  MM1 = 1, MM0 = 2
     &    , MM6 = 3, MM3 = 4
     &    , MM8 = 5, MM7 = 6
     &    , MEND = 7 )

      DOUBLE PRECISION A(REND),N(REND),E(REND)


      INTEGER IH

      DATA (A(IH),IH=1,30)  /  1.0000000000D-09, 1.0000000000D-03
     &    , 1.0000000000D-03, 1.0000000000D-03
     &    , 1.0000000000D-09, 1.0000000000D-09
     &    , 1.0000000000D-06, 1.0000000000D-03
     &    , 1.0000000000D-12, 1.0000000000D-09
     &    , 1.0000000000D-06, 1.0000000000D-03
     &    , 1.0000000000D-09, 2.0000000000D+11
     &    , 1.1581967524D+10, 5.0600000000D+01
     &    , 2.2748368939D+01, 1.0000000000D+05
     &    , 4.6548161312D+05, 1.5000000000D+06
     &    , 1.5530807739D+07, 1.8000000000D+12
     &    , 5.8737227298D+15, 2.9000000000D+11
     &    , 7.3465989712D+15, 2.2000000000D+16
     &    , 3.3416899160D+20, 2.3000000000D+12
     &    , 3.2908221995D+15, 1.5000000000D+11 /
      DATA (A(IH),IH=31,60)  /  2.5000000000D+10, 3.0000000000D+10
     &    , 1.8000000000D+10, 6.0000000000D+10
     &    , 6.3696876351D+11, 6.0000000000D+03
     &    , 1.7317292123D+06, 1.5000000000D+11
     &    , 7.1000000000D+07, 6.0000000000D+10
     &    , 3.4000000000D+09, 5.7000000000D+09
     &    , 4.4395258170D+14, 7.1000000000D+11
     &    , 9.2879665816D+08, 6.0000000000D+09
     &    , 3.7632719367D+09, 1.1000000000D+11
     &    , 4.2000000000D+10, 1.3000000000D+10
     &    , 1.2000000000D+10, 1.2000000000D+10
     &    , 4.8901418456D+09, 3.1000000000D+10
     &    , 7.2000000000D+10, 2.2424953360D+11
     &    , 5.0000000000D+13, 2.3000000000D+07
     &    , 3.4000000000D+06, 8.4300000000D+10 /
      DATA (A(IH),IH=61,90)  /  6.2570000000D+20, 6.5851871268D+25
     &    , 2.2600000000D+11, 1.0000000000D+13
     &    , 1.2720000000D+38, 5.0084110374D+44
     &    , 5.0000000000D+10, 1.8000000000D+10
     &    , 5.0000000000D+10, 3.0000000000D+10
     &    , 1.0000000000D+10, 1.3000000000D+01
     &    , 4.0307171516D-01, 1.6000000000D+04
     &    , 2.3091996575D+03, 4.7900000000D+24
     &    , 1.1516823069D+19, 4.0000000000D+10
     &    , 1.0000000000D+10, 2.0000000000D+08
     &    , 1.5000000000D+11, 3.5087770750D+10
     &    , 9.6000000000D+10, 2.0000000000D+05
     &    , 1.7200000000D+01, 1.7200000000D+01
     &    , 1.1870000000D+42, 6.2450000000D+38
     &    , 5.4200000000D+09, 2.5000000000D+14 /
      DATA (A(IH),IH=91,120)  /  6.9024230683D+09, 1.7000000000D+12
     &    , 6.2938833387D+10, 6.5000000000D+10
     &    , 1.1201744294D+10, 1.0000000000D+16
     &    , 1.5950000000D+13, 3.0000000000D+10
     &    , 1.4000000000D+06, 7.2000000000D+03
     &    , 1.5000000000D-10, 6.2741122782D-11
     &    , 2.2000000000D+13, 6.9636088998D+08
     &    , 6.4000000000D+02, 6.6104740002D+01
     &    , 9.4000000000D+03, 4.3649595891D+02
     &    , 2.0400000000D+03, 9.8081099473D+02
     &    , 4.0000000000D+10, 2.5151931669D+10
     &    , 9.9000000000D+11, 8.1645968235D+12
     &    , 4.0000000000D+03, 1.1707761726D+04
     &    , 7.5000000000D+10, 4.8785029128D+10
     &    , 2.0000000000D+17, 9.3000000000D+08 /
      DATA (A(IH),IH=121,150)  /  5.2424793425D+08, 1.0000000000D+10
     &    , 4.3883408182D+10, 9.2000000000D+10
     &    , 6.2536947345D+11, 4.0000000000D+10
     &    , 1.1669397696D+12, 5.0000000000D+08
     &    , 1.0213459815D+10, 4.6000000000D+02
     &    , 7.7713886451D+02, 2.9400000000D+11
     &    , 1.3787467866D+14, -2.2000000000D+10
     &    , 2.2000000000D+10, 1.7947817742D+11
     &    , 1.5000000000D+13, 2.3823796850D+09
     &    , 4.4000000000D+08, 2.2804120708D+08
     &    , 3.6000000000D+10, 8.6849172759D+10
     &    , 5.0000000000D+13, 1.0587819354D+10
     &    , 5.0000000000D+10, 1.8834220767D+09
     &    , 6.4000000000D+06, 1.2769711218D+06
     &    , 3.8000000000D+10, 1.3092794501D+11 /
      DATA (A(IH),IH=151,180)  /  3.3000000000D+09, 1.3645919783D+10
     &    , 1.0000000000D+08, 2.7887244677D+05
     &    , 1.0000000000D+11, 9.1001079407D+11
     &    , 1.0000000000D+11, 2.3517681424D+13
     &    , 5.0000000000D+10, 2.1179664619D+12
     &    , 3.9800000000D+14, 4.7194802843D+09
     &    , 2.2300000000D+11, 3.8793293308D+09
     &    , 2.9000000000D+10, 4.2034865839D+08
     &    , 1.0000000000D+13, 8.9133740829D+07
     &    , 2.1000000000D+09, 9.5357076078D+09
     &    , 3.5000000000D+11, 4.5766833006D+09
     &    , 1.0000000000D+10, 2.2580339609D+09
     &    , 1.1000000000D+13, 5.5002486353D+07
     &    , 2.2000000000D+04, 5.7087615208D+02
     &    , 2.2000000000D+03, 2.9032897137D+02 /
      DATA (A(IH),IH=181,209)  /  9.6000000000D+04, 2.0324987418D+05
     &    , 6.4000000000D+02, 8.7448115110D+02
     &    , 3.1000000000D+13, 2.3173078445D+09
     &    , 5.0000000000D+10, 2.7792820942D+09
     &    , 4.7000000000D+10, 1.7758624881D+10
     &    , 7.6000000000D-01, 2.5890810378D+00
     &    , 4.4000000000D+09, 3.2717637336D+10
     &    , 8.3000000000D+08, 2.9000000000D+09
     &    , 7.6702477804D+11, 1.1000000000D+11
     &    , 3.3822898074D+12, 1.4000000000D+03
     &    , 5.6360507690D+03, 3.6000000000D+05
     &    , 1.8309122037D+06, 7.8000000000D+09
     &    , 8.5223053492D+09, 4.2000000000D+10
     &    , 1.9127641509D+12, 7.2000000000D+09
     &    , 1.8988776237D+10
     &    /
      DATA (N(IH),IH=1,30)  /  0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 2.67
     &    , 2.67, 1.6
     &    , 1.6, 1.14
     &    , 1.14, -1
     &    , -1, -1
     &    , -1, -2
     &    , -2, -0.8
     &    , -0.8, 0 /
      DATA (N(IH),IH=31,60)  /  0, 0
     &    , 0, 0
     &    , 0, 1.5
     &    , 1.5, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 1.05
     &    , 1.2, 0 /
      DATA (N(IH),IH=61,90)  /  -1.8, -1.8
     &    , 0, 0
     &    , -7, -7
     &    , 0, 0
     &    , 0, 0
     &    , 0, 3
     &    , 3, 1.83
     &    , 1.83, -2.5
     &    , -2.5, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 1.5
     &    , 2.8, 2.8
     &    , -7.5, -7.5
     &    , 0, 0 /
      DATA (N(IH),IH=91,120)  /  0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 1.5, 2
     &    , 6, 6
     &    , 0, 0
     &    , 2.39, 2.39
     &    , 1.94, 1.94
     &    , 2.04, 2.04
     &    , 0, 0
     &    , -0.5, -0.5
     &    , 2, 2
     &    , 0, 0
     &    , -2.6, 0 /
      DATA (N(IH),IH=121,150)  /  0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0.5
     &    , 0.5, 2
     &    , 2, -0.4
     &    , -0.4, -0.23
     &    , -0.23, -0.23
     &    , 0, 0
     &    , 0.72, 0.72
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 1, 1
     &    , 0, 0 /
      DATA (N(IH),IH=151,180)  /  0.3, 0.3
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 1.7, 1.7
     &    , 2.11, 2.11 /
      DATA (N(IH),IH=181,209)  /  1.41, 1.41
     &    , 2, 2
     &    , -0.5, -0.5
     &    , 0, 0
     &    , 0, 0
     &    , 3, 3
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 2.1
     &    , 2.1, 1.55
     &    , 1.55, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0
     &    /
      DATA (E(IH),IH=1,30)  /  0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 0
     &    , 0, 70300000
     &    , 629666.5776, 26300000
     &    , 18522968.99, 13800000
     &    , 77295528.59, 420000
     &    , 71692559.6, 0
     &    , 434773499.6, 0
     &    , 496666802, 0
     &    , 498269028.2, 0
     &    , 195694881.3, 4200000 /
      DATA (E(IH),IH=31,60)  /  2900000, 7200000
     &    , -1700000, 0
     &    , 302574146.9, -3100000
     &    , 98652595.4, 98700000
     &    , -19000000, 0
     &    , 2900000, -3200000
     &    , 362863324.2, 70300000
     &    , 8029703.452, -7500000
     &    , 3163324.36, 3400000
     &    , 0, 6200000
     &    , 6200000, 0
     &    , 37483317.28, 0
     &    , 0, 60517145.94
     &    , 320000000, 13700000
     &    , -1900000, 0 /
      DATA (E(IH),IH=61,90)  /  0, 439072242.3
     &    , 64800000, 134000000
     &    , 11560000, 384780198.3
     &    , 105000000, 0
     &    , 105000000, 0
     &    , 30000000, 33600000
     &    , 29301257.27, 11600000
     &    , 70796785.86, 400100000
     &    , 15387108.02, 25500000
     &    , 7100000, 29300000
     &    , 0, 123728804.6
     &    , 0, 126000000
     &    , 2100000, 2100000
     &    , 190400000, 27500000
     &    , 0, 319800000 /
      DATA (E(IH),IH=91,120)  /  144563031.6, 62900000
     &    , 47548385.41, 24900000
     &    , 73043914, 126000000
     &    , -27390000, 0
     &    , 31100000, 3600000
     &    , 25400000, 47625898.51
     &    , 391000000, -62319272.88
     &    , 42600000, 24054226.73
     &    , 27050000, 727195.7218
     &    , 2370000, 47319755.32
     &    , 15280000, 46623576.77
     &    , 0, 120895988.7
     &    , 4190000, 99029105.36
     &    , 0, 125859513.2
     &    , 3870000, 0 /
      DATA (E(IH),IH=121,150)  /  2448122.663, 0
     &    , 121969739.6, 0
     &    , 317102676.6, 0
     &    , 97329442.98, 8370000
     &    , 193835268.2, 27200000
     &    , 54859109.55, 0
     &    , 164404735.9, 0
     &    , 0, 428165399.8
     &    , 204000000, -3223234.985
     &    , 2720000, 230270264.6
     &    , 0, 291045793.2
     &    , 209000000, -40911300.13
     &    , 0, 56189179.73
     &    , 26300000, 159539634.6
     &    , 0, 202909968 /
      DATA (E(IH),IH=151,180)  /  0, 313972691.2
     &    , 0, 22287354.32
     &    , 0, 457060853.9
     &    , 0, 185523159
     &    , 0, 520556382.5
     &    , 237000000, 73764195.3
     &    , 70100000, 333860663.9
     &    , 96900000, 249597940.7
     &    , 276000000, -27723594.19
     &    , -2010000, 36348379.51
     &    , 6280000, 129552874.4
     &    , 2510000, 195453207.8
     &    , 360000000, -3474999.338
     &    , 15900000, 55854923.5
     &    , 47900000, -3014915.138 /
      DATA (E(IH),IH=181,209)  /  35700000, 200974064.7
     &    , 10700000, 31057644.46
     &    , 201000000, 2632624.376
     &    , 0, 114436384.4
     &    , 0, 431539061
     &    , 16700000, 59837884.13
     &    , 92000000, 77189279.49
     &    , 67300000, -2500000
     &    , 355755324.2, 0
     &    , 299161970.7, 25600000
     &    , 30892144.94, 12600000
     &    , 93293824.28, 31200000
     &    , 48398295.69, 0
     &    , 93763000.23, -1750000
     &    , 22342666.81
     &    /
