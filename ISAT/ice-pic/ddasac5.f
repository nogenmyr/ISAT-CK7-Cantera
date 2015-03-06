c
C***********************************************************************
C
C     VERSION OF FEBRUARY 1995. NEW FEATURES SINCE MARCH 1993 INCLUDE:
C
C***********************************************************************
C
C       1. New calling sequences for DDASAC and user subroutines.
C       2. Initialization at t=t0 for u, u', W, and W' when needed, with
C          improved initial conditions and a robust damped Newton method.
C       3. Optional use of tout=t0 on the initial call of DDASAC, 
C          to permit a display of the (possibly adjusted) Y and 
C          initialized YPRIME before the first integration step.
C       4. Inclusion of diagonal, dense, and band forms for user-defined
C          E-matrices. 
C       5. Extension of the sensitivity analysis to parameter-dependent
C          u0_ and E.
C       6. Optimization of the initial stepsize on request.
C       7. Modification of DDASTP to use an Euler predictor and corrector
C          for the first step, with the initial Y and YPRIME saved in PHI.
C          PHI is converted to a backward difference table at the end of 
C          the first integration step. Phase 0 (the simultaneous increase 
C          of order and stepsize on the first few steps) is discontinued 
C          in view of item 6. 
C       8. Scaling of each iteration matrix to row maxnorms of 1.D0, to
C          reduce rounding error in the LU factorizations and solutions.
C       9. Use of LAPACK for the linear algebra computations.
C      10. Enforcement of nonnegativity and YSTOP criteria (when requested)
C          in each prediction and each corrector iteration, as well as in the 
C          initialization and initial step selection. Linear interpolation of 
C          Y(IYstop) - Cstop is used to detect crossing of Cstop; 
C          this will suffice as long as the stepsize |H| is not too large.
C      11. Computation of space requirements directly from the pointers.
C      12. Protection against overflow in DDANRM.
C      13. Stricter convergence tests in DDASTP, using the undamped
C          Newton step and starting the corrector loop with S = 100.D0 .
C          Weights for the tests of iteration convergence and truncation 
C          error are computed at the predicted states.
C      14. Reporting whenever any variable has RTOL and ATOL both zero.
C      15. Adoption of the max-norm as the default for error testing.
C      16. Use of magnifying perturbations DEL=SIGN(DEL,Y(N)) in DDAJAC 
C          to avoid inadvertent crossings through Y(N)=0.
C      17. Provision of an error message whenever the rank of the initial
C          matrix E is less than the number of nonzero rows.
C      18. New argument lists and mnemonic names.
C      19. Cosmetic changes to clarify the code.
C      20. As of May 29, 1994, nonlinear algebraic systems can be solved
C          by DDASAC with Info(11)=1, Info(13)=-1 and a dummy Esub. As
C          usual with Newton schemes, a good starting guess is needed.
C      21. Partial factorizations, guided by rounding error estimates
C          for finite-difference Jacobians, are now provided (5-31-1994) 
C          to assist the early stages of initialization.
C      22. V. Hiremath, 07/22/11, Add variable LOG to disable logging
C
C***********************************************************************
C      
C     COPYRIGHT 1995 BY MIKE CARACOTSIOS AND W. E. STEWART
C
C***********************************************************************
C
      SUBROUTINE DDASAC5( t, tout, NSTVAR, Y, YPRIME, RTOL, ATOL, INFO,
     1                    RWORK, LRW, IWORK, LIW, RPAR, IPAR, IDID, LUN,
     2                    Ieform, fsub, Esub, Jac, Bsub )
C
C***********************************************************************
C                                                                      *
C  This code solves mixed differential-algebraic equation systems      *
C  of the form                                                         *
C                                                                      *
C             E(t,U;p)*dU/dt = f(t,U;p)      ( t in [t0,tfinal] )  (1) *
C                                                                      *
C  with a coefficient matrix function E of rank .ge. 1                 *
C  and initial conditions on U(t;p) given by                           *
C                                                                      *
C             E(t,U;p)[U(t;p) - U(t0-;p)] = 0      at t=t0         (2) *
C                                                                      *
C  together with any equations in (1) that have only zeros in the      *
C  initial E. The nonzero rows of the initial E are tested for linear  *
C  independence, and an error stop is made if E fails this test. Here  *
C                                                                      *
C  t           is the timelike coordinate of the problem               *
C  t0_         is a t-value just before t0                             *
C  U(t;p)      is the state vector, approximated by the first          *
C              NSTVAR elements of the array Y(i,1+NSPAR)               *
C  E(t,U;p)    is a user-defined variable matrix of order NSTVAR       *
C  f(t,U;p)    is a user-defined vector function                       *
C  p           is a vector of parameters independent of t .            *
C                                                                      *
C  If requested, the code will solve concurrently for the parametric   *
C  sensitivity functions W(t):=partial U(t;p)/partial p .              *
C                                                                      *
C  This is an extension of Linda Petzold's implicit integrator DDASSL. *
C  The extension to parametric sensitivity analysis was done by        *
C  Mike Caracotsios.   Warren Stewart implemented Items 1-21.          *
C                                                                      *
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL DONE, LYSTOP
      DIMENSION Y(*),YPRIME(*),INFO(18)
      DIMENSION RWORK(*),IWORK(*)
      DIMENSION RTOL(*),ATOL(*)
      DIMENSION RPAR(*),IPAR(*)
      PARAMETER (ZERO=0.0D0, TENM3=1.0D-3, HALF=0.5D0, ONE=1.0D0)
      PARAMETER (FOUR=4.0D0, HUNDRD=1.0D2)
C     Fixed pointers into integer array IWORK
      PARAMETER (LIWM=1, LML=1, LMU=2, LMXORD=3, LJCALC=5,
     1           LPHASE=6, LK=7, LKOLD=8, LNS=9, LNSTL=10, LNST=11,
     2           LNRE=12, LNJE=13, LETF=14, LCTF=15, LITSTP=16, 
     3           LITINI=17, LJCINI=18, LIPVT=24)
C     Fixed pointers into double precision array RWORK
      PARAMETER (Ltstop=1, LHMAX=2, LH=3, Ltn=4, LCJ=5, LCJOLD=6,
     1           LHOLD=7, LS=8, LROUND=9, LALPHA=11, LBETA=17,
     2           LGAMMA=23, LPSI=29, LSIGMA=35,
     3           LUDELY=41, LCstop=42, LYstol=43, Ldt= 44, LDELTA=45)
      EXTERNAL  fsub, Esub, Jac, Bsub
      COMMON/DDA002/IBAND,ML,MU,MBAND,LDBAND,LDBM1,LENEWK,NALG,LENPD,
     1           NRESEV,IJAC,NONNEG,NSPAR,IDSYS,IDFDP,IRELAX,LUN0,LOG
      SAVE
      IF(INFO(1).NE.0) GO TO 100
      INFO(1)=1
      tn=t
      NRESEV=0

C     set LOG=1 to enable logging
      LOG=0
C

C
C     The following block is executed for the initial call only.
C     It checks the inputs and performs initializations.
C
C     Check the inputs INFO(2),...,INFO(11).
C
      DO 10 I=2,11
        IF(INFO(I).NE.0.AND.INFO(I).NE.1) GO TO 610
10    CONTINUE
      IJAC=INFO(5)
      NONNEG=INFO(10)
C
C     See if parametric sensitivity calculations have been requested
C     by the user. Read INFO(12) = NSPAR, the number of sensitivity
C     parameters; make sure that the parameter pointers supplied in IPAR
C     are positive and unique. Determine NEQ, the total number of
C     discretized state and sensitivity equations. Set common block
C     DDA002.
C
      LUN0=LUN
      IF(INFO(12).LT.0) GO TO 790
      NSPAR=INFO(12)
      IF(NSPAR.GT.0) THEN
        DO 30 I=1,NSPAR
          LI=I
          IF(IPAR(I).LE.0) GO TO 830
          IF(I.GT.1) THEN
            IM1=I-1
            DO 20 J=1,IM1
              LJ=J
              IF(IPAR(J).EQ.IPAR(I)) GO TO 840
20          CONTINUE
          END IF
30      CONTINUE
      END IF
      NEQ = NSTVAR*(NSPAR+1)
      IF(NEQ.LE.0) GO TO 620
      IF(INFO(13).LT.-2.OR.INFO(13).GT.2) GO TO 800
      IDSYS=INFO(13)
      IF(INFO(14).LT.0.OR.INFO(14).GT.1) GO TO 810
      IDFDP=INFO(14)
      IF(INFO(15).LT.0.OR.INFO(15).GT.1) GO TO 810
C
C     Set NWTS, the dimension of the error weight vector WT
C
      IF(INFO(15).EQ.0) THEN
        NWTS=NSTVAR
      ELSE
        NWTS=NEQ
      END IF
C
C     Set the index inorm (0 for max-norm, 1 for Euclidean)
C
      IF(INFO(16).LT.0.OR.INFO(16).GT.1) GO TO 810
      INORM=INFO(16)
C
C     See that the input value IYstop is meaningful
C
      IYstop=INFO(17)
      IF(IYstop.LT.0.OR.IYstop.GT.NEQ) GO TO 820
C
C     Set the range of the factor TEMP1 for corrector Newton steps.
C
      IF(INFO(18).LT.0.OR.INFO(18).GT.1) GO TO 810
      IRELAX=INFO(18)
C
C     Finish the pointers for IWORK. Check the dimension LIW. 
C     The array IEROW is not used if IDSYS=0. 
C
      LIEROW=LIPVT+NSTVAR
      IF(IDSYS.EQ.0) LENIW=LIEROW      
      IF(IDSYS.NE.0) LENIW=LIEROW+NSTVAR-1
      IF(LIW.LT.LENIW) GO TO 650
C
C     Initialize counters.
C
      IWORK(LNST)=0
      IWORK(LNRE)=0
      IWORK(LNJE)=0
      IWORK(LNSTL)=0
C
C     Check and store maximum integration order.
C
      MXORD=5
      IF(INFO(9).EQ.1) THEN
        MXORD=IWORK(LMXORD)
        IF(MXORD.LT.1.OR.MXORD.GT.5) GO TO 630
      END IF
      IWORK(LMXORD)=MXORD
      IF(INFO(6).EQ.0) THEN
        Iband=0
        LENPD=NSTVAR*NSTVAR
        MSAVE=0
      ELSE
C
C       Set constants for band storage of the iteration matrix.
C       Check bandwidths ML and MU.
C
        IBAND=1
        ML=IWORK(LML)
        MU=IWORK(LMU)
        IF(ML.LT.0.OR.ML.GE.NSTVAR) GO TO 760
        IF(MU.LT.0.OR.MU.GE.NSTVAR) GO TO 770
        MBAND=ML+MU+1
        LDBAND=2*ML+MU+1                
        LDBM1=LDBAND-1
        LENPD=LDBAND*NSTVAR
        MSAVE=NSTVAR/MBAND + 1
      END IF
C
C     Compute the length LENEWK of EWK, the work array for matrix E.
C
      IF(IDSYS.EQ.0) THEN
        LENEWK=1
      ELSE IF(ABS(IDSYS).EQ.1) THEN
        LENEWK=NSTVAR
      ELSE IF(INFO(6).EQ.0) THEN
        LENEWK=NSTVAR*NSTVAR
      ELSE
        LENEWK=LENPD
      END IF
C
C     Set pointers for RWORK storage. 
C
      LERLOC=LDELTA+NEQ
      LWT=LERLOC+NEQ
      LPHI=LWT+NWTS
      LPD2=LPHI+2*NEQ
      LENPHI=(IWORK(LMXORD)+1)*NEQ
      IF(INFO(8).EQ.0.AND.NSPAR.GT.0.AND.INFO(15).EQ.1) THEN
C       Here the truncation error of W affects the calculation of
C       initial stepsize. The space after column 2 of PHI will be 
C       overlaid with a second Jacobian PD2 during calculations 
C       of dW/dt for optimization of the initial H.
        LROMAX=LPHI+MAX(LENPHI,(2*NEQ+LENPD+2*MSAVE))
      ELSE
        LROMAX=LPHI+LENPHI
      END IF
      LDTEM=LROMAX+NSTVAR
      LEWK=LDTEM+NSTVAR
      LPD=LEWK+LENEWK
      LENRW=LPD+LENPD+2*MSAVE
C
C     Check the RWORK dimension LRW. Initialize DELTA,... to zero.
C
      IF(LRW.LT.LENRW) GO TO 640      
      DO 60 I=LDELTA,LENRW
        RWORK(I)=ZERO
60    CONTINUE
C
C     Compute machine constants and hmin.
C
      UROUND=DLAMCH('E')
      RMAX=DLAMCH('O')      
      RWORK(LROUND) = UROUND
      HMIN = FOUR*UROUND*MAX(ABS(t),ABS(tout))
C
C     Check hmax.
C
      IF(INFO(7).EQ.1) THEN
        HMAX=ABS(RWORK(LHMAX))
        IF(HMAX.LE.ZERO) GO TO 700
      END IF
      IF(INFO(8).EQ.1.AND.ABS(RWORK(LH)).LE.ZERO) GO TO 720
C
C     Set stopping state and tolerance when requested.
C
      IF(IYstop.GT.0) THEN
        Cstop = RWORK(LCstop)
        Ystol = ABS(RWORK(LYstol))
        IF(ABS(Ystol).LE.ZERO) Ystol=TENM3
        IF(ABS(Cstop).GT.ZERO) Ystol=MAX(Ystol,HUNDRD*UROUND)
      ELSE
        Cstop = ZERO
        Ystol = ZERO
      END IF
      LYSTOP = .FALSE.
C
C     Check entry state values for nonnegativity when requested:
C
      IDID=0
      IF(NONNEG.EQ.1) CALL DDFEAS(NSTVAR,NSTVAR,t,ZERO,DUM,DUM,Y,Y,
     1  IYstop,Cstop,.FALSE.,DUM,-1,IDUM,IWORK(LIWM),IDID)
      IF(IDID.LT.0) RETURN
C
C     Check the error tolerance settings:
C
      RTOLI=RTOL(1)
      ATOLI=ATOL(1)
      DO 80 I=1,NEQ
        IF(INFO(2).EQ.1)RTOLI=RTOL(I)
        IF(INFO(2).EQ.1)ATOLI=ATOL(I)
        IF(RTOLI.LE.ZERO.AND.ATOLI.LE.ZERO) GO TO 680
        IF(RTOLI.LT.ZERO) GO TO 660
        IF(ATOLI.LT.ZERO) GO TO 670
80    CONTINUE
C
C     Initialize the error weight vector.
C
      CALL DDAWTS(NWTS,INFO(2),RTOL,ATOL,Y,RWORK(LWT),IDID,LUN0)
      IF(IDID.LT.0) RETURN
C
C     See if the precision requested is obtainable near the initial state.
C
      UNORM=DDANRM(NSTVAR,Y,RWORK(lWT),UROUND,RMAX,INORM)
      R=UNORM*UROUND*HUNDRD
      IF(R.GT.ONE) GO TO 345
C
C     Check the direction and distance to tstop, and direction to tout.
C
      IF(INFO(4).GT.0) THEN
        tstop=RWORK(1)
        IF(ABS(tstop-t).LE.ZERO) GO TO 740
        IF(RWORK(LH)*(tstop-t).LT.ZERO) GO TO 750
      END IF
      IF(RWORK(LH)*(tout-t).LT.ZERO) GO TO 710
C
C     If INFO(11)=1, the following dt value is used in computing YPRIME.
C     t-differencing of algebraic equations is suppressed if dt=ZERO.
C
      dt=RWORK(Ldt)
      IF(ABS(dt).GT.ZERO) THEN
        DENOM=MAX(ABS(dt),SQRT(UROUND)*ABS(t))
        dt=SIGN(DENOM,dt)
      END IF
C
C     Compute consistent initial Y and YPRIME if requested,
C     and compute or adjust the initial H for DDASTP in RWORK(LH).
C
      CALL DDAINI(t,tstop,Y,YPRIME,NSTVAR,NEQ,RWORK(LH),HMAX,HMIN,
     1  dt,RTOL,ATOL,RMAX,INORM,NWTS,INFO(2),INFO(4),INFO(7),
     2  INFO(8),INFO(11),IDID,RPAR,IPAR,RWORK(LROUND),RWORK(LUDELY),
     3  RWORK(LDELTA),RWORK(LERLOC),RWORK(LWT),RWORK(LPHI),
     4  RWORK(LPD2),RWORK(LROMAX),RWORK(LDTEM),RWORK(LEWK),
     5  RWORK(LPD),IWORK(LIWM),
     6  IWORK(LJCALC),IWORK(LITINI),IWORK(LJCINI),IWORK(LIEROW),
     7  IYstop,Cstop,Ystol,LYSTOP,Ieform,fsub,Esub,Jac,Bsub)
C
C     Terminate run if IDID is negative:
C
      IF (IDID .LE. -33) RETURN
      IF (IDID .LT. 0) GO TO 490
      H=RWORK(LH)
C
C     If initial point is within the stopping tolerance, defer
C     Y-stopping until YPRIME(IYstop) changes sign.
C
      IF(IYstop.GT.0) THEN
        IF(ABS(Y(IYstop)-Cstop).LE.Ystol*(ONE+ABS(Cstop))) THEN
          LYSTOP = .FALSE.
          STSIGN = SIGN(ONE,YPRIME(IYstop))
        ELSE
          LYSTOP = .TRUE.
        END IF
      END IF
C
C     Test for return to report initialization at t0:
C
      IF(ABS(tout-t).LE.ZERO.OR.NALG.EQ.NSTVAR) THEN
        IDID=0
        RETURN
      END IF
      GO TO 310
C
C     If no step has been taken yet, proceed to call of DDASTP.
C
100   CONTINUE
      H = RWORK(LH)
      IF(IWORK(LNST).EQ.0) GO TO 340
C
C     The following block is for continuation calls
C     only. Here we check INFO(1), and if the
C     last step was interrupted we check whether
C     appropriate action was taken.
C
      IF(INFO(1).EQ.1) GO TO 110
      IF(INFO(1).NE.-1) GO TO 610
C
C     If control reaches here, the last step was interrupted
C     by an error condition from DDASTP and
C     appropriate action was not taken. This
C     is a fatal error.
C
      CALL XERRDD(
     1     49HDDASAC-- THE LAST STEP TERMINATED WITH A NEGATIVE,
     2     49,0,0,0,0,0,ZERO,ZERO,LUN)
      CALL XERRDD(
     1     47HDDASAC-- VALUE (=I1) OF IDID AND NO APPROPRIATE,
     2     47,0,1,IDID,0,0,ZERO,ZERO,LUN)
      CALL XERRDD(
     1     41HDDASAC-- ACTION WAS TAKEN. RUN TERMINATED,
     2     41,1,0,0,0,0,ZERO,ZERO,LUN)
      RETURN
110   CONTINUE
      IWORK(LNSTL)=IWORK(LNST)
C
C     This block is for continuation calls only. Its purpose
C     is to check stop conditions before taking a step.
C     Adjust h if necessary to meet hmax bound.
C
210   CONTINUE
      DONE = .FALSE.
      tn=RWORK(Ltn)
      H=RWORK(LH)
      IF(INFO(7) .EQ. 1) THEN
        IF(ABS(H) .GT. HMAX) H = SIGN(HMAX,H)
      END IF
      IF(ABS(t - tout) .LE. ZERO) GO TO 780
      IF((t - tout)*H .GT. ZERO) GO TO 710
      IF(INFO(4) .EQ. 1) GO TO 250
      IF(INFO(3) .EQ. 1) GO TO 230
      IF((tn-tout)*H.LT.ZERO) GO TO 300
      CALL DDATRP (tn,tout,Y,YPRIME,NEQ,IWORK(LKOLD),
     1             RWORK(LPHI),RWORK(LPSI))
      t=tout
      IDID = 3
      DONE = .TRUE.
      GO TO 300
230   IF((tn-t)*H .LE. ZERO) GO TO 300
      IF((tn - tout)*H .GT. ZERO) GO TO 240
      CALL DDATRP (tn,tn,Y,YPRIME,NEQ,IWORK(LKOLD),
     1             RWORK(LPHI),RWORK(LPSI))
      t = tn
      IDID = 1
      DONE = .TRUE.
      GO TO 300
240   CONTINUE
      CALL DDATRP (tn,tout,Y,YPRIME,NEQ,IWORK(LKOLD),
     1             RWORK(LPHI),RWORK(LPSI))
      t = tout
      IDID = 3
      DONE = .TRUE.
      GO TO 300
250   IF(INFO(3) .EQ. 1) GO TO 260
      tstop=RWORK(Ltstop)
      IF(INFO(4).GT.0.AND.(tn-tstop)*H.GT.ZERO) GO TO 750
      IF(INFO(4).GT.0.AND.(tstop-tout)*H.LT.ZERO) GO TO 690
      IF((tn-tout)*H.LT.ZERO) GO TO 280
      CALL DDATRP (tn,tout,Y,YPRIME,NEQ,IWORK(LKOLD),
     1             RWORK(LPHI),RWORK(LPSI))
      t=tout
      IDID = 3
      DONE = .TRUE.
      GO TO 300
260   tstop = RWORK(Ltstop)
      IF(INFO(4).GT.0.AND.(tn-tstop)*H .GT. ZERO) GO TO 750
      IF(INFO(4).GT.0.AND.(tstop-tout)*H .LT. ZERO) GO TO 690
      IF((tn-t)*H .LE. ZERO) GO TO 280
      IF((tn - tout)*H .GT. ZERO) GO TO 270
      CALL DDATRP (tn,tn,Y,YPRIME,NEQ,IWORK(LKOLD),
     1             RWORK(LPHI),RWORK(LPSI))
      t = tn
      IDID = 1
      DONE = .TRUE.
      GO TO 300
270   CONTINUE
      CALL DDATRP (tn,tout,Y,YPRIME,NEQ,IWORK(LKOLD),
     1             RWORK(LPHI),RWORK(LPSI))
      t = tout
      IDID = 3
      DONE = .TRUE.
      GO TO 300
280   CONTINUE
C
C     Check whether we are within roundoff of tstop.
C
      IF(INFO(4).GT.0) THEN
        IF(ABS(tn-tstop).LE.HUNDRD*UROUND*(ABS(tn)+ABS(H))) THEN
          CALL DDATRP(tn,tstop,Y,YPRIME,NEQ,IWORK(LKOLD),
     1    RWORK(LPHI),RWORK(LPSI))
          IDID=2
          t=tstop
          DONE = .TRUE.
        END IF
      ELSE
        tnext=tn+H
        IF((tnext-tstop)*H.LE.ZERO) GO TO 300
        H=(tstop-tn)
        RWORK(LH)=H
      END IF
300   IF (DONE) GO TO 480
C
C     The next block contains the call to the
C     one-step integrator DDASTP.
C     This is a looping point for the integration
C     steps.
C
C     Check for too many steps.
C
310   IF(IWORK(LNST)-IWORK(LNSTL).GT.500) THEN
        IDID=-1
        GO TO 370
      END IF
C
C     Test for too much accuracy requested.
C
340   R=DDANRM(NWTS,RWORK(LPHI),RWORK(LWT),UROUND,RMAX,INORM)
     1  *HUNDRD*UROUND
345   IF(R.GT.ONE) THEN
C
C     Multiply RTOL and ATOL by r and return.
C
        IF(INFO(2).EQ.0) THEN
          RTOL(1)=R*RTOL(1)
          ATOL(1)=R*ATOL(1)
          IDID=-2
          GO TO 370
        ELSE
          CALL DSCAL(NWTS,R,RTOL,1)
          CALL DSCAL(NWTS,R,ATOL,1)
        END IF
        IDID=-2
        GO TO 370
      END IF
C
C     Compute minimum stepsize.
C
      HMIN=FOUR*UROUND*MAX(ABS(tn),ABS(tout))
      IF(IYstop.GT.0) Ystold=Y(IYstop)
      CALL DDASTP (tn,H,NSTVAR,NEQ,Y,YPRIME,RTOL,ATOL,RMAX,INORM,NWTS,
     1  INFO(2),IDID,RPAR,IPAR,RWORK(LDELTA),RWORK(LERLOC),RWORK(LWT),
     2  RWORK(LPHI),RWORK(LROMAX),RWORK(LDTEM),RWORK(LEWK),RWORK(LPD),
     3  IWORK(LIWM),RWORK(LALPHA),
     4  RWORK(LBETA),RWORK(LGAMMA),RWORK(LPSI),RWORK(LSIGMA),
     5  RWORK(LCJ),RWORK(LCJOLD),RWORK(LHOLD),RWORK(LS),HMIN,
     6  RWORK(LROUND),RWORK(LUDELY),IWORK(LITSTP),IWORK(LPHASE),
     7  IWORK(LJCALC),IWORK(LK),IWORK(LKOLD),IWORK(LNS),INFO(7),HMAX,
     8  IYstop,Cstop,LYSTOP,IWORK(LIEROW),Ieform,fsub,Esub,Jac,Bsub)
      IF(IDID.LE.-33) THEN
        t = tn
        RETURN
      END IF
370   IF(IDID.LT.0) GO TO 490
C
C     This block handles successful returns from DDASTP (IDID.ge.0). 
C     Test for stop conditions.
C
      IF(IYstop.GT.0) THEN
C
C       Analyze the algebraic constraint imposed on Y(IYstop).
C
        IF(IYstop.GT.0.AND.STSIGN*YPRIME(IYstop).LT.ZERO) 
     1                                    LYSTOP = .TRUE.
        IF(LYSTOP) THEN
          IF(ABS(Y(IYstop)-Cstop).LE.Ystol*(ONE+ABS(Cstop)).OR.
     1    (Y(IYstop)-Cstop)*SIGN(ONE,Ystold-Cstop).LE.ZERO) THEN
            IDID = 4
            t = tn
            GO TO 480
          END IF
        END IF
      END IF
      IF(INFO(4).NE.0) GO TO 410
      IF(INFO(3).NE.0) GO TO 390
      IF((tn-tout)*H.LT.ZERO) GO TO 310
      CALL DDATRP(tn,tout,Y,YPRIME,NEQ,IWORK(LKOLD),
     1            RWORK(LPHI),RWORK(LPSI))
      IDID=3
      t=tout
      GO TO 470
390   IF((tn-tout)*H.GE.ZERO) GO TO 400
      t=tn
      IDID=1
      GO TO 470
400   CALL DDATRP(tn,tout,Y,YPRIME,NEQ,IWORK(LKOLD),
     1            RWORK(LPHI),RWORK(LPSI))
      IDID=3
      t=tout
      GO TO 470
410   IF(INFO(3).NE.0) GO TO 440
      IF((tn-tout)*H.LT.ZERO) GO TO 420
      CALL DDATRP(tn,tout,Y,YPRIME,NEQ,IWORK(LKOLD),
     1            RWORK(LPHI),RWORK(LPSI))
      t=tout
      IDID=3
      GO TO 470
420   IF(ABS(tn-tstop).LE.HUNDRD*UROUND*(ABS(tn)+ABS(H))) GO TO 430
      tnext=tn+H*(ONE+FOUR*UROUND)
      IF((tnext-tstop)*H.LE.ZERO) GO TO 310
      H=(tstop-tn)*(ONE-FOUR*UROUND)
      GO TO 310
430   IDID=2
      t=tstop
      GO TO 470
440   IF((tn-tout)*H.GE.ZERO) GO TO 460
      IF(ABS(tn-tstop).LE.HUNDRD*UROUND*(ABS(tn)+ABS(H))) GO TO 450
      t=tn
      IDID=1
      GO TO 470
450   IDID=2
      t=tstop
      GO TO 470
460   CALL DDATRP(tn,tout,Y,YPRIME,NEQ,IWORK(LKOLD),
     1            RWORK(LPHI),RWORK(LPSI))
      t=tout
      IDID=3
470   CONTINUE
C
C     All successful returns from DDASAC are made from
C     this block. Return to the calling program.
C
480   CONTINUE
      RWORK(Ltn)=tn
      RWORK(LH)=H
      RETURN
C
C     This block handles all unsuccessful
C     returns other than for illegal input.
C
490   CONTINUE
      ITEMP=-IDID
      GO TO (500,510,520,600,600,530,540,550,560,570,580,590) ITEMP
C
C     The maximum number of steps was taken before
C     reaching tout.
C
500   CALL XERRDD(
     1     38HDDASAC-- AT CURRENT t (=R1)  500 STEPS,
     2     38,0,0,0,0,1,tn,ZERO,LUN)
      CALL XERRDD(
     1     37HDDASAC-- HAVE BEEN TAKEN ON THIS CALL,
     2     37,0,0,0,0,0,ZERO,ZERO,LUN)
      GO TO 600
C
C     Tolerances too small for machine precision.
C
510   CALL XERRDD(
     1     47HDDASAC-- AT t (=R1) TOO MUCH ACCURACY REQUESTED,
     2     47,0,0,0,0,1,tn,ZERO,LUN)
      CALL XERRDD(
     1     48HDDASAC-- FOR PRECISION OF MACHINE. RTOL AND ATOL,
     2     48,0,0,0,0,0,ZERO,ZERO,LUN)
      CALL XERRDD(
     1     45HDDASAC-- WERE INCREASED TO APPROPRIATE VALUES,
     2     45,0,0,0,0,0,ZERO,ZERO,LUN)
      GO TO 600
C
C     WT(i) .le. ZERO for some i (not at start of problem).
C
520   CALL XERRDD(
     1     38HDDASAC-- AT t (=R1) SOME ELEMENT OF WT,
     2     38,0,0,0,0,1,tn,ZERO,LUN)
      CALL XERRDD(
     1     28HDDASAC-- HAS BECOME .LE. 0.0,
     2     28,0,0,0,0,0,ZERO,ZERO,LUN)
      GO TO 600
C
C     Error test failed repeatedly or with h=hmin.
C
530   CALL XERRDD(
     1     44HDDASAC-- AT t (=R1) AND STEPSIZE H (=R2) THE,
     2     44,0,0,0,0,2,tn,H,LUN)
      CALL XERRDD(
     1     57HDDASAC-- ERROR TEST FAILED REPEATEDLY OR WITH ABS(H)=HMIN,
     2     57,0,0,0,0,0,ZERO,ZERO,LUN)
      GO TO 600
C
C     Corrector convergence failed repeatedly or with h=hmin.
C
540   CALL XERRDD(
     1     44HDDASAC-- AT t (=R1) AND STEPSIZE H (=R2) THE,
     2     44,0,0,0,0,2,tn,H,LUN)
      CALL XERRDD(
     1     48HDDASAC-- CORRECTOR FAILED TO CONVERGE REPEATEDLY,
     2     48,0,0,0,0,0,ZERO,ZERO,LUN)
      CALL XERRDD(
     1     28HDDASAC-- OR WITH ABS(H)=HMIN,
     2     28,0,0,0,0,0,ZERO,ZERO,LUN)
      GO TO 600
C
C     The iteration matrix is singular.
C
550   CALL XERRDD(
     1     44HDDASAC-- AT t (=R1) AND STEPSIZE H (=R2) THE,
     2     44,0,0,0,0,2,tn,H,LUN)
      CALL XERRDD(
     1     37HDDASAC-- ITERATION MATRIX IS SINGULAR,
     2     37,0,0,0,0,0,ZERO,ZERO,LUN)
      GO TO 600
C
C     Corrector failure preceded by error test failures.
C
560   CALL XERRDD(
     1     44HDDASAC-- AT t (=R1) AND STEPSIZE H (=R2) THE,
     2     44,0,0,0,0,2,tn,H,LUN)
      CALL XERRDD(
     1     49HDDASAC-- CORRECTOR COULD NOT CONVERGE.  ALSO, THE,
     2     49,0,0,0,0,0,ZERO,ZERO,LUN)
      CALL XERRDD(
     1     38HDDASAC-- ERROR TEST FAILED REPEATEDLY.,
     2     38,0,0,0,0,0,ZERO,ZERO,LUN)
      GO TO 600
C
C     Corrector failure because IRES = -1.
C
570   CALL XERRDD(
     1     44HDDASAC-- AT t (=R1) AND STEPSIZE H (=R2) THE,
     2     44,0,0,0,0,2,tn,H,LUN)
      CALL XERRDD(
     1     45HDDASAC-- CORRECTOR COULD NOT CONVERGE BECAUSE,
     2     45,0,0,0,0,0,ZERO,ZERO,LUN)
      CALL XERRDD(
     1     36HDDASAC-- IRES WAS EQUAL TO MINUS ONE,
     2     36,0,0,0,0,0,ZERO,ZERO,LUN)
      GO TO 600
C
C     Failure because IRES = -2.
C
580   CALL XERRDD(
     1     40HDDASAC-- AT t (=R1) AND STEPSIZE H (=R2),
     2     40,0,0,0,0,2,tn,H,LUN)
      CALL XERRDD(
     1     36HDDASAC-- IRES WAS EQUAL TO MINUS TWO,
     2     36,0,0,0,0,0,ZERO,ZERO,LUN)
      GOTO 600
C
C     DDAINI subroutine failed to complete the initialization.
C
590   CALL XERRDD(
     1     48HDDASAC-- AT t (=R1) THE INITIAL Y OR YPRIME OR h,
     2     48,0,0,0,0,1,t,ZERO,LUN)
      CALL XERRDD(
     1     31HDDASAC-- COULD NOT BE COMPLETED,
     2     31,0,0,0,0,0,ZERO,ZERO,LUN)
C
C     Restore values and return
C
600   CONTINUE
      INFO(1)=-1
      t=tn
      RWORK(Ltn)=tn
      RWORK(LH)=H
      RETURN
C
C     This block handles all error returns due
C     to illegal input, as detected before calling
C     DDASTP. First the error message routine is
C     called. If this happens twice in
C     succession, execution is terminated.
C
610   CALL XERRDD(
     1     51HDDASAC-- AN ELEMENT OF INFO(1,...,11) IS NOT 0 OR 1,
     2     51,0,0,0,0,0,ZERO,ZERO,LUN)
      GO TO 850
620   CALL XERRDD(
     1     25HDDASAC-- NEQ (=I1) .LE. 0,
     2     25,0,1,NEQ,0,0,ZERO,ZERO,LUN)
      GO TO 850
630   CALL XERRDD(
     1     34HDDASAC-- MAXORD (=I1) NOT IN RANGE,
     2     34,0,1,MXORD,0,0,ZERO,ZERO,LUN)
      GO TO 850
640   CALL XERRDD(
     160HDDASAC-- RWORK LENGTH NEEDED, LENRW (=I1), EXCEEDS LRW (=I2),
     260,0,2,LENRW,LRW,0,ZERO,ZERO,LUN)
      GO TO 850
650   CALL XERRDD(
     160HDDASAC-- IWORK LENGTH NEEDED, LENIW (=I1), EXCEEDS LIW (=I2),
     260,0,2,LENIW,LIW,0,ZERO,ZERO,LUN)
      GO TO 850
660   CALL XERRDD(
     1     39HDDASAC-- SOME ELEMENT OF RTOL IS .LT. 0,
     2     39,0,0,0,0,0,ZERO,ZERO,LUN)
      GO TO 850
670   CALL XERRDD(
     1     39HDDASAC-- SOME ELEMENT OF ATOL IS .LT. 0,
     2     39,0,0,0,0,0,ZERO,ZERO,LUN)
      GO TO 850
680   CALL XERRDD(
     1     55HDDASAC-- SOME VARIABLE HAS RTOL AND ATOL BOTH .LE. ZERO,
     2     47,0,0,0,0,0,ZERO,ZERO,LUN)
      GO TO 850
690   CALL XERRDD(
     1     54HDDASAC-- INFO(4) = 1 AND tstop (=R1) BEFORE tout (=R2),
     2     54,0,0,0,0,2,tstop,tout,LUN)
      GO TO 850
700   CALL XERRDD(
     1     28HDDASAC-- HMAX (=R1) .LE. 0.0,
     2     28,0,0,0,0,1,HMAX,ZERO,LUN)
      GO TO 850
710   CALL XERRDD(
     1     34HDDASAC-- tout (=R1) BEFORE t (=R2),
     2     34,0,0,0,0,2,tout,t,LUN)
      GO TO 850
720   CALL XERRDD(
     1     29HDDASAC-- INFO(8)=1 AND H0=0.0,
     2     29,0,0,0,0,0,ZERO,ZERO,LUN)
      GO TO 850
740   CALL XERRDD(
     1 58HDDASAC-- t (=R1) TOO NEAR tstop (=R2) TO START INTEGRATION,
     2 58,0,0,0,0,2,t,tstop,LUN)
      GO TO 850
750   CALL XERRDD(
     1     49HDDASAC-- INFO(4)=1 AND tstop (=R1) BEFORE t (=R2),
     2     49,0,0,0,0,2,tstop,t,LUN)
      GO TO 850
760   CALL XERRDD(
     1     52HDDASAC-- ML (=I1) ILLEGAL. EITHER .LT. 0 OR .GT. NEQ,
     2     52,0,1,IWORK(LML),0,0,ZERO,ZERO,LUN)
      GO TO 850
770   CALL XERRDD(
     1     52HDDASAC-- MU (=I1) ILLEGAL. EITHER .LT. 0 OR .GT. NEQ,
     2     52,0,1,IWORK(LMU),0,0,ZERO,ZERO,LUN)
      GO TO 850
780   CALL XERRDD(
     1     49HDDASAC-- STEP FROM t (=R1) TO tout (=R2) IS ZERO.,
     2     49,0,0,0,0,2,t,tout,LUN)
      GO TO 850
790   CALL XERRDD(
     1     38HDDASAC-- INFO(12) (=NSPAR) IS NEGATIVE,
     2     38,0,0,0,0,0,ZERO,ZERO,LUN)
      GO TO 850
800   CALL XERRDD(
     1     44HDDASAC-- INVALID INFO(13): .LT. -2 OR .GT. 2,
     2     44,0,1,NEQ,0,0,ZERO,ZERO,LUN)
      GO TO 850
810   CALL XERRDD(
     1     49HDDASAC-- INVALID INFO(14,15,16 OR 18): NOT 0 OR 1,
     2     49,0,0,0,0,0,ZERO,ZERO,LUN)
      GO TO 850
820   CALL XERRDD(
     1     51HDDASAC-- INVALID INFO(17): .LT. 1 OR .GT. NEQ (=I1),
     2     51,0,1,NEQ,0,0,ZERO,ZERO,LUN)
      GO TO 850
830   CALL XERRDD(
     1     46HDDASAC-- NONPOSITIVE POINTER VALUE AT IPAR(I1),
     2     46,0,1,LI,0,0,ZERO,ZERO,LUN)
      GO TO 850
840   CALL XERRDD(
     1     52HDDASAC-- DUPLICATE POINTERS AT IPAR(I1) AND IPAR(I2),
     2     52,0,2,LI,LJ,0,ZERO,ZERO,LUN)
      GO TO 850
850   IF(INFO(1).EQ.-1) GO TO 860
      INFO(1)=-1
      IDID=-33
      RETURN
860   CALL XERRDD(
     1     46HDDASAC-- REPEATED OCCURRENCES OF ILLEGAL INPUT,
     2     46,0,0,0,0,0,ZERO,ZERO,LUN)
870   CALL XERRDD(
     1     47HDDASAC-- RUN TERMINATED. APPARENT INFINITE LOOP,
     2     47,1,0,0,0,0,ZERO,ZERO,LUN)
      RETURN
C
C     End of Subroutine DDASAC
C
      END
      SUBROUTINE DDASTP (X,H,NSTVAR,NEQ,Y,YPRIME,RTOL,ATOL,RMAX,INORM,
     1           NWTS,INFO2,IDID,RPAR,IPAR,DELTA,ERLOC,WT,
     2           PHI,ROMAX,DTEM,EWK,PD,IWM,ALPHA,BETA,GAMMA,PSI,SIGMA,
     3           CJ,CJOLD,HOLD,S,HMIN,UROUND,UDEL,ITMAX,IPHASE,
     4           JCALC,K,KOLD,NS,INFO7,HMAX,IYstop,Cstop,LYSTOP,
     5           IEROW,Ieform,fsub,Esub,Jac,Bsub)
C
C     DDASTP solves a system of differential and algebraic equations
C     for one t-step. The initial Y and YPRIME are saved in PHI
C     for predictions during the first t-step. Backward differences are
C     formed in PHI at the successful completion of the first t-step, 
C     and are used for all interpolations and predictions from then on.
C     
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL CONVGD, LYSTOP
      DIMENSION Y(*),YPRIME(*),RTOL(*),ATOL(*)
      DIMENSION RPAR(*),IPAR(*),DELTA(*),ERLOC(*),WT(*),PHI(NEQ,*)
      DIMENSION ROMAX(*),DTEM(*),EWK(*),PD(*),IWM(*),DUM1(1)
      DIMENSION ALPHA(*),BETA(*),GAMMA(*),PSI(*),SIGMA(*),IEROW(*)
      PARAMETER (LMXORD=3, LNST=11, LNRE=12, LNJE=13, LETF=14, LCTF=15)
      PARAMETER (ZERO=0.0D0, TENM4=1.0D-4, TENM3=1.0D-3, DOT33=0.33D0)
      PARAMETER (XRATE=0.25D0, QUARTR=0.25D0, HALF=0.5D0, DOT9=0.9D0)
      PARAMETER (ONE=1.0D0, TWO=2.0D0, HUNDRD=1.0D2)
      EXTERNAL  fsub, Esub, Jac, Bsub
      COMMON/DDA002/IBAND,ML,MU,MBAND,LDBAND,LDBM1,LENEWK,NALG,LENPD,
     1              NRESEV,IJAC,NONNEG,NSPAR,IDSYS,IDFDP,IRELAX,LUN0
C
C     Block 1.
C     Initializations for all calls. 
C
      IDID=1
      MAXIT=MAX(4,ITMAX)
      XOLD=X
      NCF=0
      NSF=0
      NEF=0
      IF(IWM(LNST).EQ.0) THEN
C
C       Initialize constants for the Euler corrector. Set iphase=1
C       since we are starting with a near-maximal step. 
C
        IWM(LETF) = 0
        IWM(LCTF) = 0
        K=1
        KOLD=0
        HOLD=ZERO
        JCALC=-1
        DELNRM=HUNDRD
        KP1=K+1
        KP2=K+2
        KM1=K-1
        NS=1
        ALPHA(1)=ONE
        ALPHA(2)=HALF
        BETA(1)=ONE
        BETA(2)=ONE
        SIGMA(1)=ONE
        SIGMA(2)=HALF
        CK=HALF  
        JSTART=0
      ELSE
        JSTART=1
      END IF
C
C     Block 2.
C     Compute constants for this step.
C
20    CONTINUE
      IF(JSTART.GT.0) THEN
        KP1=K+1
        KP2=K+2
        KM1=K-1
C        XOLD=X
        IF(H.NE.HOLD.OR.K .NE. KOLD) NS = 0
        NS=MIN(NS+1,KOLD+2)
        NSP1=NS+1
        IF(KP1 .GE. NS) THEN
          BETA(1)=ONE
          ALPHA(1)=ONE
          TEMP1=H
          GAMMA(1)=ZERO
          SIGMA(1)=ONE
          DO 30 I=2,KP1
            TEMP2=PSI(I-1)
            PSI(I-1)=TEMP1
            BETA(I)=BETA(I-1)*PSI(I-1)/TEMP2
            TEMP1=TEMP2+H
            ALPHA(I)=H/TEMP1
            SIGMA(I)=DBLE(I-1)*SIGMA(I-1)*ALPHA(I)
            GAMMA(I)=GAMMA(I-1)+ALPHA(I-1)/H
30        CONTINUE
          PSI(KP1)=TEMP1
        END IF
        ALPHAS = ZERO
        ALPHA0 = ZERO
        DO 50 I = 1,K
          ALPHAS = ALPHAS - ONE/DBLE(I)
          ALPHA0 = ALPHA0 - ALPHA(I)
50      CONTINUE
C
C       Compute leading coefficient cj.
C
        CJLAST = CJ
        CJ = -ALPHAS/H
C
C       Compute variable stepsize error coefficient ck.
C
        CK = ABS(ALPHA(KP1) + ALPHAS - ALPHA0)
        CK = MAX(CK,ALPHA(KP1))
C
C       Change phi to phi star.
C
        IF(JSTART.GT.0) THEN
          IF(KP1 .GE. NSP1) THEN
            DO 70 J=NSP1,KP1
              DO 60 I=1,NEQ
                PHI(I,J)=BETA(J)*PHI(I,J)
60            CONTINUE
70          CONTINUE
          END IF
        END IF
      END IF
C
C     Block 3.
C     Predict the solution and derivative. 
C
120   CONTINUE
      IF(JSTART.EQ.0) THEN
C       Prediction from initial Y and YPRIME, saved in PHI:
        CALL DCOPY(NEQ,PHI,1,Y,1)
        CALL DCOPY(NEQ,PHI(1,2),1,YPRIME,1)
        CALL DAXPY(NEQ,H,YPRIME,1,Y,1)
C
C       Test feasibility of predicted step, and shorten if necessary.
C      
        IF(LYSTOP.OR.NONNEG.GT.0) THEN
          CALL DDFEAS(NSTVAR,NEQ,XOLD,H,HFEAS,HMIN,PHI,Y,
     1    IYstop,Cstop,LYSTOP,DOT9,1,IDUM,IWM,IDID)
          IF(IDID.LT.0) GO TO 610
          H=HFEAS
        END IF
        PSI(1)=H    
        PSI(2)=TWO*H
        CJOLD=ONE/H 
        CJ=CJOLD    
      ELSE
C       Prediction from backward differences, saved in PHI:
        DO 130 I=1,NEQ
          Y(I)=PHI(I,1)
          YPRIME(I)=ZERO
130     CONTINUE
        DO 150 J=2,KP1
          DO 140 I=1,NEQ
            Y(I)=Y(I)+PHI(I,J)
            YPRIME(I)=YPRIME(I)+GAMMA(J)*PHI(I,J)
140       CONTINUE
150     CONTINUE
C
C       Test feasibility of predicted step, and shorten if necessary.
C
        IF(LYSTOP.OR.NONNEG.GT.0) THEN
          CALL DDFEAS(NSTVAR,NEQ,XOLD,H,HFEAS,HMIN,PHI,Y,
     1    IYstop,Cstop,LYSTOP,DOT9,0,IDUM,IWM,IDID)
          IF(ABS(HFEAS).LT.ABS(H)) GO TO 515
        END IF
      END IF
C
C     Update time. 
C
      X = XOLD + H
C
C     Compute the weights and PNORM at the predicted point.
C
      CALL DDAWTS(NWTS,INFO2,RTOL,ATOL,Y,WT,IDID,LUN0)
      IF(IDID.LT.0) RETURN
      PNORM = DDANRM(NSTVAR,Y,WT,UROUND,RMAX,INORM)
C
C     Decide whether new jacobian is needed.
C
      TEMP1 = (ONE - XRATE)/(ONE + XRATE)
      TEMP2 = ONE/TEMP1
      IF (CJ/CJOLD .LT. TEMP1 .OR. CJ/CJOLD .GT. TEMP2) JCALC = -1
C
C     Block 4. Solve the corrector equation by a
C     modified Newton scheme.
C
      CONVGD= .TRUE.
      M=0
      IRES = 0
      IWM(LNRE)=IWM(LNRE)+1
      CALL DDARES(3,X,ZERO,NSTVAR,DUM1,Y,YPRIME,DELTA,DUM1,IRES,RPAR,
     1  IPAR,EWK,IEROW,PD,IWM,Ieform,fsub,Esub)
      IF (IRES .LT. 0) GO TO 290
C
C     If indicated, reevaluate the
C     iteration matrix PD(i,j) = dF(i)/dY(j) + CJ*dF(i)/dY'(j)
C     Set jcalc to 0 as an indicator that this has been done.
C
      IF(JCALC .EQ. -1) THEN
        JCALC=0
        IWM(LNJE)=IWM(LNJE)+1
        CALL DDAJAC(3,X,ZERO,NSTVAR,DUM1,Y,YPRIME,DELTA,CJ,H,
     1       DTEM,DUM1,IER,WT,ROMAX,PD,IWM,IRES,UROUND,UDEL,LRANK,
     2       RPAR,IPAR,EWK,IEROW,Ieform,fsub,Esub,Jac)
        CJOLD=CJ
        IF (IRES .LT. 0) GO TO 290
        IF(IER .NE. 0) GO TO 290
        NSF=0
      END IF
C
C     Initialize the error accumulation vector ERLOC.
C
      DO 170 I=1,NEQ
        ERLOC(I)=ZERO
170   CONTINUE
C
C     Begin corrector loop.
C
180   CONTINUE
C
C     Compute the Newton correction vector for the current iteration.
C     Store it in DELTA, and compute its norm for the convergence test.
C
      CALL DDASLV(NSTVAR,NSTVAR,DELTA,DELTA,ROMAX,PD,IWM,IER)
      DELNRM=DDANRM(NSTVAR,DELTA,WT,UROUND,RMAX,INORM)
C
C     Compute multiplier temp1 to accelerate convergence.
C
      TEMP1 = TWO/(ONE + CJ/CJOLD)
      IF(IRELAX.NE.0)TEMP1=MIN(TEMP1,ONE)
C
C     Reduce TEMP1 if necessary for feasibility
C     and multiply the vector DELTA by TEMP1.
C
      IF(LYSTOP.OR.NONNEG.GT.0) THEN
        ZTEMP=TEMP1
        CALL DCOPY(NSTVAR,Y,1,DTEM,1)
        CALL DAXPY(NSTVAR,-TEMP1,DELTA,1,DTEM,1)
        CALL DDFEAS(NSTVAR,NSTVAR,X,ZTEMP,TEMP1,TENM3,Y,DTEM,
     1    IYstop,Cstop,LYSTOP,DOT9,2,M,IWM,IDID)
        IF(IDID.LT.0) GO TO 610
      END IF
      CALL DSCAL(NSTVAR,TEMP1,DELTA,1)
C
C     Update Y,ERLOC,and YPRIME.
C
      DO 200 I=1,NSTVAR
        Y(I)=Y(I)-DELTA(I)
        ERLOC(I)=ERLOC(I)-DELTA(I)
        YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
200   CONTINUE
      IF(NONNEG.GT.0) THEN
        DO 210 I=1,NSTVAR
          Y(I)=MAX(Y(I),ZERO)
210     CONTINUE
      END IF
C
C     Test convergence, using the norm DELNRM of the full Newton step.
C
      IF (DELNRM .LE. HUNDRD*UROUND*PNORM) GO TO 240
      IF (M .EQ. 0) THEN
        OLDNRM = DELNRM
        S = HUNDRD
      ELSE
        RATE = (DELNRM/OLDNRM)**(ONE/DBLE(M))
        IF (RATE .GT. DOT9) GO TO 230
        S = RATE/(ONE - RATE)
      END IF
      IF (S*DELNRM .LE. DOT33) GO TO 240
C
C     The corrector has not yet converged.
C     Update m and test whether the
C     maximum number of iterations have
C     been tried.
C
      M=M+1
      IF(M.GE.MAXIT) GO TO 230
C
C     Evaluate the residual
C     and go back to do another iteration.
C
      IRES = 0
      IWM(LNRE)=IWM(LNRE)+1
      CALL DDARES (3,X,ZERO,NSTVAR,DUM1,Y,YPRIME,DELTA,DUM1,IRES,RPAR,
     1     IPAR,EWK,IEROW,PD,IWM,Ieform,fsub,Esub)
      IF (IRES .LT. 0) GO TO 290
      GO TO 180
C
C     The corrector failed to converge in maxit
C     iterations. If the iteration matrix
C     is not current,re-do the step with
C     a new iteration matrix.
C
230   CONTINUE
      IF(JCALC.EQ.0) GO TO 290
      JCALC=-1
      GO TO 120
240   CONTINUE
C
C     Block 5. Parametric sensitivity computation, done on request
C     after convergence of the modified Newton method. 
C
270   IF(NSPAR.LE.0) GO TO 300
C
C     Update the right hand side of the state equations
C     and the Jacobian matrix only if necessary.
C     Store the update of the right hand side in DELTA(i,1).
C
      IRES=0
      IF(IJAC.EQ.0) THEN
        IWM(LNRE)=IWM(LNRE)+1
        CALL DDARES (3,X,ZERO,NSTVAR,DUM1,Y,YPRIME,DELTA,DUM1,IRES,
     1     RPAR,IPAR,EWK,IEROW,PD,IWM,Ieform,fsub,Esub)
        IF(IRES.LT.0) GO TO 560
      END IF
      IF(JCALC.EQ.0) GO TO 280
      IWM(LNJE)=IWM(LNJE)+1
      CALL DDAJAC (3,X,ZERO,NSTVAR,DUM1,Y,YPRIME,DELTA,CJ,H,
     1     DTEM,DUM1,IER,WT,ROMAX,PD,IWM,IRES,UROUND,UDEL,LRANK,
     2     RPAR,IPAR,EWK,IEROW,Ieform,fsub,Esub,Jac)
      IF(IRES.LT.0) GO TO 560
      IF(IER.NE.0) GO TO 290
280   CONTINUE
      CALL DDSENS (3,X,ZERO,NSTVAR,DUM1,Y,YPRIME,ERLOC,DELTA,DTEM,
     1     DUM1,RPAR,IPAR,EWK,IEROW,ROMAX,PD,IWM,CJ,UROUND,UDEL,
     2     IRES,Ieform,fsub,Esub,Bsub)
      IF(IRES.LT.0) GO TO 560
      GO TO 300
C
C     Singular iteration matrix or no convergence
C     with current iteration matrix.
C
290   CONVGD= .FALSE.
      IF(IRES.EQ.-2) GO TO 560
300   JCALC = 1
      IF(.NOT.CONVGD) GO TO 510
C
C     Block 6.
C     Estimate the errors at orders k,k-1,k-2
C     as if stepsize were constant. Estimate
C     the local error at order k and test
C     whether the current step is successful.
C
      ENORM = DDANRM(NWTS,ERLOC,WT,UROUND,RMAX,INORM)
      ERK = SIGMA(K+1)*ENORM
      TERK = DBLE(K+1)*ERK
      EST = ERK
      KNEW=K
      IF(K .EQ. 1) GO TO 350
      DO 310 I = 1,NEQ
        DELTA(I) = PHI(I,KP1) + ERLOC(I)
310   CONTINUE
      ERKM1=SIGMA(K)*DDANRM(NWTS,DELTA,WT,UROUND,RMAX,INORM)
      TERKM1 = DBLE(K)*ERKM1
      IF(K .GT. 2) GO TO 320
      IF(TERKM1 .LE. HALF*TERK) GO TO 340
      GO TO 350
320   CONTINUE
      DO 330 I = 1,NEQ
        DELTA(I) = PHI(I,K) + DELTA(I)
330   CONTINUE
      ERKM2=SIGMA(K-1)*DDANRM(NWTS,DELTA,WT,UROUND,RMAX,INORM)
      TERKM2 = DBLE(K-1)*ERKM2
      IF(MAX(TERKM1,TERKM2).GT.TERK) GO TO 350
C
C     Lower the order.
C
340   CONTINUE
      KNEW=K-1
      EST = ERKM1
C
C     Calculate the local error for the current step
C     to see if the step was successful.
C
350   CONTINUE
      ERR = CK * ENORM
      IF(ERR .GT. ONE) GO TO 510
C
C     Block 7.
C     The step is successful. 
C     Determine the best order and stepsize for
C     the next step. Update the differences
C     for the next step.
C
      IDID=1
      IWM(LNST)=IWM(LNST)+1
      KDIFF=K-KOLD
      KOLD=K
      HOLD=H
C
C     Estimate the error at order k+1 unless
C     already decided to lower order, or
C     already using maximum order, or
C     stepsize not constant, or
C     order raised in previous step. 
C
      IF(KNEW.EQ.KM1) GO TO 390
      IF(K.EQ.IWM(LMXORD)) GO TO 420
      IF(KP1.GE.NS.OR.KDIFF.EQ.1) GO TO 420
      DO 360 I=1,NEQ
        DELTA(I)=ERLOC(I)-PHI(I,KP2)
360   CONTINUE
      ERKP1=(ONE/DBLE(K+2))*DDANRM(NWTS,DELTA,WT,UROUND,RMAX,INORM)
      TERKP1 = DBLE(K+2)*ERKP1
      IF(K.GT.1) GO TO 370
      IF(TERKP1.GE.HALF*TERK) GO TO 420
      GO TO 380
370   IF(TERKM1.LE.MIN(TERK,TERKP1)) GO TO 390
      IF(TERKP1.GE.TERK.OR.K.EQ.IWM(LMXORD)) GO TO 420
C
C     Raise order.
C
380   K=KP1
      EST = ERKP1
      GO TO 420
C
C     Lower order.
C
390   K=KM1
      EST = ERKM1
C
C     Determine the appropriate stepsize for the next step.
C
420   HNEW=H
      R=(TWO*EST+TENM4)**(-ONE/DBLE(K+1))
      IF(R .LT. TWO) GO TO 430
      HNEW = TWO*H
      GO TO 440
430   IF(R .GT. ONE) GO TO 440
      R = MAX(HALF,MIN(DOT9,R))
      HNEW = H*R
440   H=HNEW
C
C     Adjust h if necessary to meet hmax bound
C
450   IF(INFO7.EQ.1.AND.ABS(H).GT.HMAX) H = SIGN(HMAX,H)
C
C     Update differences for next step.
C
460   CONTINUE
      IF(JSTART.GT.0) THEN
        IF(KOLD.NE.IWM(LMXORD)) THEN
          DO 470 I=1,NEQ
            PHI(I,KP2)=ERLOC(I)
470       CONTINUE
480     END IF
        DO 490 I=1,NEQ
          PHI(I,KP1)=PHI(I,KP1)+ERLOC(I)
490     CONTINUE
        DO 500 J1=2,KP1
          J=KP1-J1+1
          DO 500 I=1,NEQ
            PHI(I,J)=PHI(I,J)+PHI(I,J+1)
500     CONTINUE
      ELSE
        DO 505 I=1,NEQ
          PHI(I,2)=Y(I)-PHI(I,1)
          PHI(I,1)=Y(I)
505     CONTINUE
      END IF
      RETURN
C
C     Block 8.
C     The step is unsuccessful. Restore x,psi,phi.
C     Determine appropriate stepsize for another
C     attempt, or exit with an error flag 
C     if there have been many failures.
C
C     Restore x,phi,psi.
C
510   IPHASE = 1
515   X=XOLD
      IF(JSTART.GT.0) THEN
        IF(KP1.GE.NSP1) THEN
          DO 530 J=NSP1,KP1
            TEMP1=ONE/BETA(J)
            DO 520 I=1,NEQ
              PHI(I,J)=TEMP1*PHI(I,J)
520         CONTINUE
530       CONTINUE
        END IF
        DO 550 I=2,KP1
          PSI(I-1)=PSI(I)-H
550     CONTINUE
      END IF
      IF(IDID.EQ.-33) GO TO 610
      IF(LYSTOP.OR.NONNEG.GT.0) THEN
        IF(ABS(HFEAS).LT.ABS(H)) THEN
          H=HFEAS
          GO TO 20
        END IF
      END IF
C
C     Test whether failure is due to corrector iteration
C     or error test.
C
      IF(CONVGD) GO TO 580
      IWM(LCTF)=IWM(LCTF)+1
C
C     The newton iteration failed to converge with
C     a current iteration matrix.  Determine the cause
C     of the failure and take appropriate action.
C
      IF(IER.EQ.0) GO TO 560
C
C     The iteration matrix is singular. Reduce
C     the stepsize by a factor of 4. If
C     this happens three times in a row on
C     the same step, return with an error flag.
C
      NSF=NSF+1
      R = QUARTR
      H=H*R
      IF (NSF .LT. 3 .AND. ABS(H) .GE. HMIN) GO TO 620
      IDID=-8
      GO TO 610
C
C     The newton iteration failed to converge for a reason
C     other than a singular iteration matrix.  If IRES = -2, then
C     return.  Otherwise, reduce the stepsize and try again, unless
C     too many failures have occured.
C
560   CONTINUE
      IF (IRES .GT. -2) GO TO 570
      IDID = -11
      GO TO 610
570   NCF = NCF + 1
      R = QUARTR
      H = H*R
      IF (NCF .LT. 10 .AND. ABS(H) .GE. HMIN) GO TO 620
      IDID = -7
      IF (IRES .LT. 0) IDID = -10
      IF (NEF .GE. 3) IDID = -9
      GO TO 610
C
C     The newton scheme converged, and the cause
C     of the failure was the error estimate
C     exceeding the tolerance.
C
580   NEF=NEF+1
      IWM(LETF)=IWM(LETF)+1
      IF (NEF .EQ. 1) THEN
C
C     On first error test failure, keep current order or lower
C     order by one.  Compute new stepsize based on differences
C     of the solution.
C
        K = KNEW
        TEMP2 = K + 1
        R = DOT9*(TWO*EST+TENM4)**(-ONE/TEMP2)
        R = MAX(QUARTR,MIN(DOT9,R))
        H = H*R
        IF (ABS(H) .GE. HMIN) GO TO 620
        IDID = -6
590   ELSE IF (NEF .EQ. 2) THEN
C
C     On second error test failure, use the current order or
C     decrease order by one.  Reduce the stepsize by a factor of
C     one quarter.
C
        K = KNEW
        H = QUARTR*H
        IF (ABS(H) .GE. HMIN) GO TO 620
        IDID = -6
600   ELSE IF (NEF .GE. 3) THEN
C
C     On third and subsequent error test failures, set the order to
C     one and reduce the stepsize by a factor of one quarter.
C
        K = 1
        H = QUARTR*H
        IF (ABS(H) .GE. HMIN) GO TO 620
        IDID = -6
      END IF
C
C     For all crashes, restore Y to its value at X=XOLD,
C     interpolate to find YPRIME there, and return.
C
610   CONTINUE
      IF(JSTART.GT.0) THEN
        CALL DDATRP(X,X,Y,YPRIME,NEQ,K,PHI,PSI)
      ELSE
        CALL DCOPY(NEQ,PHI,1,Y,1)
        CALL DCOPY(NEQ,PHI(1,2),1,YPRIME,1)
      END IF
      RETURN
C
C     Go back and try this step again.
C
620   GO TO 20
C
C     End of Subroutine DDASTP
C
      END
      SUBROUTINE DDAINI(t0,tstop,Y,YPRIME,NSTVAR,NEQ,H,HMAX,HMIN,
     1  dt,RTOL,ATOL,RMAX,INORM,NWTS,INFO2,INFO4,INFO7,INFO8,INFO11,
     2  IDID,RPAR,IPAR,UROUND,UDEL,DELTA,ERLOC,WT,PHI,PD2,ROMAX,DTEM,
     3  EWK,PD,IWM,JCALC,ITINI,JCINI,IEROW,
     4  IYstop,Cstop,Ystol,LYSTOP,Ieform,fsub,Esub,Jac,Bsub)
C
C     DDAINI computes, on request, initial values U and dU/dt consistent 
C     with the DAE system. The parametric sensitivity matrix W is initialized
C     when present. The initial dW/dt is calculated if the error control 
C     includes W, that is, if DDASAC was called with INFO(12)>0 and
C     INFO(15) = 1. On request, an optimal initial step for the integration 
C     is computed by a call to DFINDH. CJ is set to zero throughout DDAINI 
C     since no backward differences occur in the initialization equations.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL LYSTOP
      PARAMETER (MGSAVE=10, NRMIJO=2, LNST=11, LNRE=12, LNJE=13)
      PARAMETER (ZERO=0.0D0, TENM3=1.0D-3, FIFTH=0.2D0, HALF=0.5D0)
      PARAMETER (DOT9=0.9D0, ONE=1.0D0, CHSIGN=-1.0D0, TWO=2.0D0)
      DIMENSION Y(*),YPRIME(*),RTOL(*),ATOL(*),RPAR(*),IPAR(*)
      DIMENSION DELTA(*),ERLOC(*),WT(*),PHI(NEQ,*),PD2(*),ROMAX(*)
      DIMENSION DTEM(*),EWK(*),PD(*),IWM(*),IEROW(*)
      DIMENSION CORNRM(MGSAVE),RESNRM(MGSAVE),DUM1(1)
      EXTERNAL  fsub, Esub, Jac, Bsub, IDAMAX
      COMMON/DDA002/IBAND,ML,MU,MBAND,LDBAND,LDBM1,LENEWK,NALG,LENPD,
     1              NRESEV,IJAC,NONNEG,NSPAR,IDSYS,IDFDP,IRELAX,LUN0
      IDID=1
      IF(IDSYS.EQ.0) LRANK=NSTVAR
C
C     Save input state vector u0 in phi(i,1). 
C
      CALL DCOPY(NSTVAR,Y,1,PHI,1)
      IF(INFO11.EQ.0) THEN
        IF(INFO8.EQ.1) GO TO 100
        IF(IDSYS.NE.0) JCALC=1
        GO TO 60
      END IF
C
C     Block 1a. Complete the state vector u (first NSTVAR elements of Y) 
C     at t=t0. Start the computation by calling DDARES with JOB=0.
C              
      MAXIT=ITINI
      IF(ITINI.EQ.0) MAXIT=50
      MJAC=JCINI
      IF(JCINI.EQ.0) MJAC=5
      DO 30 J=1,MAXIT
        IRES=0
        IWM(LNRE)=IWM(LNRE)+1
C       DDARES alters only the arguments DELTA, EWK and IEROW when JOB=0.
        CALL DDARES (0,t0,ZERO,NSTVAR,PHI,Y,DUM1,DELTA,DUM1,IRES,
     1      RPAR,IPAR,EWK,IEROW,PD,IWM,Ieform,fsub,Esub)
        IF (IRES.LT.0) THEN
          IDID=-12
          RETURN
        END IF
C       Evaluate current weights for error tests and perturbation step
C       selection.
        IDID=0
        CALL DDAWTS(NWTS,INFO2,RTOL,ATOL,Y,WT,IDID,LUN0)
        IF(IDID.LT.0) RETURN
C
C       Evaluate and factorize the initial iteration matrix, if necessary.
C       CJ and H are set to zero so that YPRIME (which does not appear
C       in the initial conditions) will be ignored in this call of DDAJAC.
C
        IF(IDSYS.NE.0.AND.(MOD(J,MJAC).EQ.1.OR.JCALC.NE.0)) THEN
          JCALC=0
          IWM(LNJE)=IWM(LNJE)+1
          CALL DDAJAC(0,t0,ZERO,NSTVAR,PHI,Y,YPRIME,DELTA,ZERO,ZERO,
     1         DTEM,DUM1,IER,WT,ROMAX,PD,IWM,IRES,UROUND,UDEL,LRANK,
     2         RPAR,IPAR,EWK,IEROW,Ieform,fsub,Esub,Jac)
          IF(LRANK.EQ.0.OR.IRES.LT.0.OR.IER.LT.0) THEN
            IDID=-12
            RETURN
          END IF
          IF(LRANK.LT.NSTVAR) JCALC=1
        END IF
        JSAVE=MOD(J,MGSAVE)
        RESNRM(JSAVE)=DDANRM(NSTVAR,DELTA,ROMAX,UROUND,RMAX,INORM)
        CALL DCOPY(NSTVAR,DELTA,1,ERLOC,1)
        IF(IDSYS.NE.0) THEN
          CALL DDASLV(LRANK,NSTVAR,DELTA,ERLOC,ROMAX,PD,IWM,IER)
          IF(IER.NE.0) THEN
            IDID=-12
            RETURN
          END IF
        END IF
        CORNRM(JSAVE)=DDANRM(NSTVAR,DELTA,WT,UROUND,RMAX,INORM)
        IF(LRANK.EQ.NSTVAR.AND.CORNRM(JSAVE).LT.ONE) THEN
C         Convergence achieved. Latest residual vector is saved in ERLOC.
          GO TO 40 
        END IF
C
C     If the Newton step exceeds the convergence specification, seek
C     a descent of either the residual norm or the correction norm.
C     For J.LE.NRMIJO(:=2 here), the comparison norms are taken from 
C     the start of the current Newton step; for larger J, they are
C     the largest of the values saved in RESNRM AND CORNRM.
C
        LSAVE=MIN(J,MGSAVE)
        IF(J.LE.NRMIJO) THEN
          RSNRM0=RESNRM(JSAVE)
          DLNRM0=CORNRM(JSAVE)
        ELSE
          RSNRM0=RESNRM(IDAMAX(LSAVE,RESNRM,1))
          DLNRM0=CORNRM(IDAMAX(LSAVE,CORNRM,1))
        END IF
C       Initialize damping factor
        DAMP=CHSIGN
C       Save the origin of the current step in DTEM:
        CALL DCOPY(NSTVAR,Y,1,DTEM,1) 
        IF(NONNEG.GT.0) THEN
C         Adjust damping factor if necessary:
          IDID=0
          CALL DAXPY(NSTVAR,DAMP,DELTA,1,Y,1)
          CALL DDFEAS(NSTVAR,NSTVAR,t0,CHSIGN,DAMP,TENM3,DTEM,Y,
     1         IDUM,DUM1(1),.FALSE.,DOT9,2,J,IWM,IDID)
          IF(IDID.LT.0) RETURN
        END IF
        DO 20 IDAMP=1,10
C         Set Y to origin of iteration:
          CALL DCOPY(NSTVAR,DTEM,1,Y,1) 
C         Compute Y for damped step:
          CALL DAXPY(NSTVAR,DAMP,DELTA,1,Y,1) 
          IRES=0
          IWM(LNRE)=IWM(LNRE)+1
          CALL DDARES (0,t0,ZERO,NSTVAR,PHI,Y,DUM1,ERLOC,DUM1,IRES,
     1        RPAR,IPAR,EWK,IEROW,PD,IWM,Ieform,fsub,Esub)
          IF(IRES.EQ.-1) GO TO 10
          IF(IRES.LT.-1) THEN
            IDID=-12
            RETURN
          END IF
C         ERLOC now contains the residual vector at the latest Y.
          CALL DCOPY(NSTVAR,ERLOC,1,PHI(1,2),1)
          TRES=DDANRM(NSTVAR,ERLOC,ROMAX,UROUND,RMAX,INORM)
          IER=0
C         Compute a trial correction to the state, in PHI(*,2):
          IF(IDSYS.NE.0) THEN
            CALL DDASLV(LRANK,NSTVAR,PHI(1,2),ERLOC,ROMAX,PD,IWM,IER)
            IF(IER.NE.0) THEN
              IDID=-12
              RETURN
            END IF
          END IF
          TDY=DDANRM(NSTVAR,PHI(1,2),WT,UROUND,RMAX,INORM)
C         Full-rank convergence of Y within specified tolerances:
          IF(TDY.LT.ONE.AND.LRANK.EQ.NSTVAR) GO TO 40 
C         Descent achieved. Do another iteration if permitted:
          IF((TRES.LT.RSNRM0).OR.(TDY.LT.DLNRM0)) GO TO 30 
C         Reduce the damping factor and try again:
10        DAMP=DAMP*HALF
20      CONTINUE
C       Request new Jacobian if line search failed:
        JCALC=-1 
30    CONTINUE
C
C     Exit for failure to converge on initial state.
C
      IF(LOG.EQ.1) WRITE(LUN0,35) MAXIT, CORNRM(JSAVE)
35    FORMAT(' INITIALIZATION DID NOT CONVERGE IN',I6,' ITERATIONS.'
     1  ,/,' NORM OF LAST NEWTON STEP IS ',1P,D12.2,' VS. 1.0D0',
     2  ' PERMITTED.',/,' USE OF A LARGER ITERATION LIMIT',
     3  ' IN IWORK(17) MAY HELP.',/,1X)
      IDID=-12
      RETURN
40    CONTINUE
C     
C     Test max-norm of U-U0_ for control of Jacobian updating
C
      DO 50 I=1,NSTVAR
        DTEM(I)=Y(I)-PHI(I,1)
50    CONTINUE
      TOTDY=DDANRM(NSTVAR,DTEM,WT,UROUND,RMAX,0)
C
C     Block 1b. Initialize W if parametric sensitivities were requested.
C
      IF(NSPAR.GT.0) THEN
C       Update the Jacobian if necessary for the computation of W.
        IF(TOTDY.GE.ONE.AND.(IDSYS.GT.0.OR.NALG.GT.0)) THEN
C         Load the latest residual vector into the front of DELTA 
C         as a base value for subroutine DDAJAC:
          CALL DCOPY(NSTVAR,ERLOC,1,DELTA,1)
          IWM(LNJE)=IWM(LNJE)+1
          CALL DDAJAC(0,t0,ZERO,NSTVAR,PHI,Y,YPRIME,DELTA,ZERO,ZERO,
     1         DTEM,DUM1,IER,WT,ROMAX,PD,IWM,IRES,UROUND,UDEL,LRANK,
     2         RPAR,IPAR,EWK,IEROW,Ieform,fsub,Esub,Jac)
          IF(IRES.LT.0.OR.IER.NE.0) THEN
            IDID=-12
            RETURN
          END IF
        END IF
C       Load the same residual vector into the front of DELTA 
C       as a base value for sensitivity computations:
        CALL DCOPY(NSTVAR,ERLOC,1,DELTA,1)
C       This call of DDSENS computes W in elements NSTVAR+1,... of DELTA.
        CALL DDSENS(0,t0,ZERO,NSTVAR,PHI,Y,DUM1,DUM1,DELTA,DTEM,
     1       DUM1,RPAR,IPAR,EWK,IEROW,ROMAX,PD,IWM,ZERO,UROUND,UDEL,
     2       IRES,Ieform,fsub,Esub,Bsub)
        IF(IRES.LT.0) THEN
          IDID=-12
          RETURN
        END IF
C       Insert the sensitivity elements into Y.
        LCOPY=NSPAR*NSTVAR
        CALL DCOPY(LCOPY,DELTA(NSTVAR+1),1,Y(NSTVAR+1),1)
      END IF 
C     If E is zero at the input state, return to report the
C     solution of the algebraic problem.
      IF(NALG.EQ.NSTVAR) RETURN
C
C     Block 2. Analysis of YPRIME at and/or near t=t0.
C     Start by computing a new Jacobian if necessary, with JOB=-1:
      IF(TOTDY.GE.ONE.AND.(IDSYS.GT.0.OR.NALG.GT.0)) JCALC=1
60    IF(JCALC.EQ.1) THEN
        IF(INFO11.EQ.1) THEN
C         Load the latest residual vector into the front of DELTA 
C         as a base value for subroutine DDAJAC:
          CALL DCOPY(NSTVAR,ERLOC,1,DELTA,1)
        ELSE
C         Compute a base residual vector DELTA at the given initial state:
          IRES=0
          IWM(LNRE) = IWM(LNRE)+1
          CALL DDARES (0,t0,ZERO,NSTVAR,PHI,Y,DUM1,DELTA,DUM1,IRES,
     1        RPAR,IPAR,EWK,IEROW,PD,IWM,Ieform,fsub,Esub)
          IF(IRES.LT.0) THEN
            IDID=-12
            RETURN
          END IF
        END IF
C       Use JOB=-1 in this call of DDAJAC to suppress derivatives of EWK:
        IWM(LNJE)=IWM(LNJE)+1
        CALL DDAJAC(-1,t0,ZERO,NSTVAR,PHI,Y,YPRIME,DELTA,ZERO,ZERO,
     1       DTEM,DUM1,IER,WT,ROMAX,PD,IWM,IRES,UROUND,UDEL,LRANK,
     2       RPAR,IPAR,EWK,IEROW,Ieform,fsub,Esub,Jac)
        IF(IRES.LT.0.OR.IER.NE.0) THEN
          IDID=-12
          RETURN
        END IF
      END IF
      IF(INFO11.EQ.1.OR.INFO8.EQ.0) THEN
        ITIME=0
        CALL DDAYPR(ITIME,t0,Y,YPRIME,NSTVAR,NEQ,H,dt,WT,NWTS,IDID,
     1  RPAR,IPAR,PHI,PD2,DELTA,ERLOC,ROMAX,PD,IWM,UROUND,UDEL,EWK,
     2  IEROW,DTEM,Ieform,fsub,Esub,Jac,Bsub) 
        IF(IDID.LT.0) RETURN
      END IF
C
C     Save Y(t0+) and YPRIME(t0+) in PHI.
C
100   CALL DCOPY(NEQ,Y,1,PHI,1)
      CALL DCOPY(NEQ,YPRIME,1,PHI(1,2),1)
C
C     Block 3. Determine a feasible initial step.
C     First find an upper bound HUP on the stepsize.
C
      HUP=SQRT(RMAX)
      IF(INFO4.GT.0) HUP=MIN(HUP,ABS(tstop-t0))
      IF(INFO7.GT.0) HUP=MIN(HUP,HMAX)
      IF(IYstop.GT.0) THEN
        IF(ABS(Y(IYstop)-Cstop).GT.Ystol*(ONE+ABS(Cstop))) THEN
          LYSTOP = .TRUE.
          CRIT=SIGN(ONE,H)*YPRIME(IYstop)*SIGN(ONE,Cstop-Y(IYstop))
          IF(CRIT.GT.ZERO) THEN
            HSTOP=(Cstop-Y(IYstop))/YPRIME(IYstop)
            HUP=MIN(HUP,ABS(HSTOP))
          END IF
        END IF
      END IF
      IF(INFO8.EQ.0) THEN
C       Compute a trial stepsize H from Eq. (6) of Shampine, Applied
C       Numerical Mathematics 3, North-Holland (1987), 331-337.
        YPNRM0=DDANRM(NSTVAR,YPRIME,WT,UROUND,RMAX,INORM)
        IF(YPNRM0.GT.FIFTH/HUP) THEN
          H=SIGN(FIFTH/YPNRM0,H)
        ELSE
          H=SIGN(HUP,H)
        END IF
C
C       Use Shampine's method to approximate the largest feasible H 
C       such that all shorter steps satisfy the local error tolerances.
C    
        CALL DFINDH(t0,Y,YPRIME,NSTVAR,NEQ,dt,H,HUP,HMIN,WT,
     1    RTOL,ATOL,RMAX,INORM,NWTS,INFO2,IDID,
     2    RPAR,IPAR,PHI,PD2,DELTA,ERLOC,ROMAX,PD,IWM,UROUND,UDEL,EWK,
     3    IEROW,DTEM,IYstop,Cstop,LYSTOP,Ieform,fsub,Esub,Jac,Bsub)
        CALL DCOPY(NEQ,PHI,1,Y,1)
        CALL DCOPY(NEQ,PHI(1,2),1,YPRIME,1)
        IF(IDID.LT.0) RETURN
      ELSE
        HMAG=MIN(HUP,ABS(H))
        H=SIGN(HMAG,H)
      END IF
      RETURN
C
C     End of Subroutine DDAINI
C
      END
      SUBROUTINE DDAYPR(ITIME,t,Y,YPRIME,NSTVAR,NEQ,H,dt,WT,NWTS,
     1  IDID,RPAR,IPAR,PHI,PD2,DELTA,ERLOC,ROMAX,PD,IWM,UROUND,UDEL,
     2  EWK,IEROW,DTEM,Ieform,fsub,Esub,Jac,Bsub) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(*),YPRIME(*),WT(*),RPAR(*),IPAR(*)
      DIMENSION PHI(NEQ,*),PD2(*),DELTA(*),ERLOC(*),ROMAX(*),PD(*)
      DIMENSION IWM(*),EWK(*),IEROW(*),DTEM(*),DUM1(1)
      PARAMETER (LNRE=12, LNJE=13)
      PARAMETER (ZERO=0.0D0, HALF=0.5D0, ONE=1.0D0, CHSIGN=-1.0D0)
      PARAMETER (TWO=2.0D0, HUNDRD=100.0D0)
      EXTERNAL  fsub, Esub, Jac, Bsub, IDAMAX
      COMMON/DDA002/IBAND,ML,MU,MBAND,LDBAND,LDBM1,LENEWK,NALG,LENPD,
     1              NRESEV,IJAC,NONNEG,NSPAR,IDSYS,IDFDP,IRELAX,LUN0
C
C     Block 1. Compute u'= du/dt at time t, using as reference the vector 
C     DELTA(1,...NSTVAR) returned by DDARES with JOB=1. The input dt, if 
C     nonzero, must be small enough to give a good secant approximation to 
C     partial FA/partial t. If dt is zero, partial FA/partial t is taken
C     to be zero, shortening the execution of DDARES.
C
      IWM(LNRE)=IWM(LNRE)+1
      CALL DDARES (1,t,dt,NSTVAR,DUM1,Y,DUM1,DELTA,DTEM,IRES,RPAR,IPAR,
     1      EWK,IEROW,PD,IWM,Ieform,fsub,Esub)
      IF (IRES.LT.0) THEN
        IDID=-12
        RETURN
      END IF
C
C     The vector [-FD -dFA/dt] has been returned in DELTA. Now 
C     copy DELTA to YPRIME and compute du/dt there: 
C
      CALL DCOPY(NSTVAR,DELTA,1,YPRIME,1)
      CALL DSCAL(NSTVAR,CHSIGN,YPRIME,1)
      IF(IDSYS.NE.0) THEN
        CALL DDASLV(NSTVAR,NSTVAR,YPRIME,YPRIME,ROMAX,PD,IWM,IER) 
        IF(IER.NE.0) THEN
          IDID=-12
          RETURN
        END IF
      END IF
C     Return if parametric sensitivity analysis is not requested.
      IF(NSPAR.EQ.0) RETURN
C
C     Block 2. Initialization of WPRIME for sensitivity analysis.
C
C     If the sensitivity matrix is not to be used in the error control,
C     initialize WPRIME to zero on the first call.
C
      IF(NWTS.EQ.NSTVAR) THEN
        IF(ITIME.EQ.0) THEN
          DO 20 I=NSTVAR+1,NEQ
            YPRIME(I)=ZERO
20        CONTINUE
        END IF
      ELSE
C
C       The initial WPRIME is needed for error control.
C       Compute WPRIME(t0) when ITIME=0, and
C       approximate WPRIME(t0+H) to order H when ITIME=1.
C
C       Compute a base residual vector DELTA(1,...,NSTVAR) with JOB=2;
C       save this vector for use by DDSENS (and by DDAJAC if needed).
        IWM(LNRE)=IWM(LNRE)+1
        CALL DDARES (2,t,dt,NSTVAR,DUM1,Y,YPRIME,DELTA,DTEM,IRES,
     1            RPAR,IPAR,EWK,IEROW,PD,IWM,Ieform,fsub,Esub)
        IF (IRES.LT.0) THEN
          IDID=-12
          RETURN
        END IF
C       Compute parametric derivatives of -[Eu'-FD-partial FA/partial t]
C       in DELTA(NSTVAR+1,...), using ERLOC to store RHS2 during DDARES.
        CALL DDSENS (2,t,dt,NSTVAR,PHI,Y,YPRIME,DUM1,DELTA,DTEM,
     1       ERLOC,RPAR,IPAR,EWK,IEROW,ROMAX,PD,IWM,ZERO,UROUND,UDEL,
     2       IRES,Ieform,fsub,Esub,Bsub)
        IF(IRES.LT.0) THEN
          IDID=-12
          RETURN
        END IF
        LW=NSTVAR*NSPAR
        CALL DCOPY(LW,DELTA(NSTVAR+1),1,YPRIME(NSTVAR+1),1)        
        IF(ITIME.EQ.0) THEN
C         Determine the biggest initial sensitivity coefficient.
          WBIG=ABS(Y(NSTVAR+IDAMAX(LW,Y(NSTVAR+1),1)))    
        END IF
        IF(WBIG.GT.ZERO.OR.ITIME.EQ.1) THEN
C         Compute the Jacobian of [ Eu' - FD - partial FA/partial t ];
C         store it in PD2.
          IWM(LNJE)=IWM(LNJE)+1
          CALL DDAJAC(2,t,dt,NSTVAR,DUM1,Y,YPRIME,DELTA,ZERO,ZERO,
     1         DTEM,ERLOC,IER,WT,ROMAX,PD2,IWM,IRES,UROUND,UDEL,LRANK,
     2         RPAR,IPAR,EWK,IEROW,Ieform,fsub,Esub,Jac)
          IF(IRES.LT.0.OR.IER.NE.0) THEN
            IDID=-12
            RETURN
          END IF
C         Subtract PD2*W(t0) from YPRIME(NSTVAR+1,...). 
          L=1
          DO 40 J=2,NSPAR+1
            L=L+NSTVAR
            IF(IBAND.EQ.0) THEN
              CALL DGEMV('N',NSTVAR,NSTVAR,
     1                      CHSIGN,PD2,NSTVAR,Y(L),1,ONE,YPRIME(L),1)
            ELSE
              CALL DGBMV('N',NSTVAR,NSTVAR,ML,MU,
     1                CHSIGN,PD2(1+ML),LDBAND,Y(L),1,ONE,YPRIME(L),1)
            END IF
40        CONTINUE
C         Solve for the first approximation to 
C         [WPRIME:=YPRIME(NSTVAR+1,...)] at time t.
C         If ITIME=0, this result is YPRIME(t0).
          IF(IDSYS.NE.0) THEN
            L=1
            DO 60 J=2,NSPAR+1
              L=L+NSTVAR
              CALL DDASLV
     1          (NSTVAR,NSTVAR,YPRIME(L),YPRIME(L),ROMAX,PD,IWM,IER) 
              IF(IER.NE.0) THEN
                IDID=-12
                RETURN
              END IF
  60        CONTINUE
          END IF
        END IF
        IF(ITIME.EQ.0) THEN
          RETURN  
        ELSE
C         Compute the second approximation to dW/dt.
          CALL DCOPY(LW,PHI(NSTVAR+1,1),1,ERLOC(NSTVAR+1),1)
          CALL DAXPY(LW,H,PHI(NSTVAR+1,2),1,ERLOC(NSTVAR+1),1)
          L=1
          DO 80 J=2,NSPAR+1
            L=L+NSTVAR
            IF(IBAND.EQ.0) THEN
              CALL DGEMV('N',NSTVAR,NSTVAR,
     1                      CHSIGN,PD2,NSTVAR,ERLOC(L),1,ZERO,DTEM,1)
            ELSE
              CALL DGBMV('N',NSTVAR,NSTVAR,ML,MU,
     1                CHSIGN,PD2(1+ML),LDBAND,ERLOC(L),1,ZERO,DTEM,1)
            END IF
            IF(IDSYS.NE.0) THEN
              CALL DDASLV(NSTVAR,NSTVAR,DTEM,DTEM,ROMAX,PD,IWM,IER)
              IF(IER.NE.0) THEN
                IDID=-12
                RETURN
              END IF
            END IF
            CALL DAXPY(NSTVAR,ONE,DTEM,1,YPRIME(L),1)
80        CONTINUE
        END IF
      END IF
      RETURN
C
C     End of Subroutine DDAYPR
C
      END
      SUBROUTINE DFINDH(t0,Y,YPRIME,NSTVAR,NEQ,dt,H,HUP,HMIN,WT,
     1  RTOL,ATOL,RMAX,INORM,NWTS,INFO2,IDID,
     2  RPAR,IPAR,PHI,PD2,DELTA,ERLOC,ROMAX,PD,IWM,UROUND,UDEL,EWK,
     3  IEROW,DTEM,IYstop,Cstop,LYSTOP,Ieform,fsub,Esub,Jac,Bsub)
C
C     This subroutine is called when INFO(8)=0 to find a good
C     initial step H for the predictor-corrector computation in DDASTP.
C     The approach of Shampine, Applied Numerical Mathematics 3 (1987)
C     133-137 is used, with du/dt or dY/dt replacing his right-hand 
C     function f(x,y). The tolerance vector WT is evaluated at the 
C     current forward-Euler state, Y1 = Y0 + h*Y'0.
C

C     Modified by Qing Tang, 4/99
C     ---------------------------
C     Shampine iteration may fail if the right-hand 
C     function f(x,y) is discontinuous/non-monotonic.
C     This modification makes the algorithm robust.
C     If the Shampine iteration can not find an optimal
C     initial time step, a largest value will be chosed 
C     as the return value from the failed trial time steps
C     which satisfied the condition 4 and 5.
C     
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

	integer, save :: n_warn=0, lim_warn=10

      LOGICAL LYSTOP, PASSD5
      DIMENSION Y(*),YPRIME(*),WT(*),RTOL(*),ATOL(*),RPAR(*),IPAR(*)
      DIMENSION PHI(NEQ,*),DELTA(*),ERLOC(*),ROMAX(*),PD(*),IWM(*)
      DIMENSION EWK(*),IEROW(*),DTEM(*), PD2(*)
      PARAMETER (LNRE=12, LNJE=13)
      PARAMETER (ZERO=0.0D0, TENM3=1.0D-3, TENM2=1.0D-2, TENTH=0.1D0)
      PARAMETER (HALF=0.5D0, DOT8=0.8D0, DOT9=0.9D0)
      PARAMETER (ONE=1.0D0, ONEP2=1.2D0, TWO=2.0D0, HUNDRD=100.D0)
      EXTERNAL  fsub, Esub, Jac, Bsub
      COMMON/DDA002/IBAND,ML,MU,MBAND,LDBAND,LDBM1,LENEWK,NALG,LENPD,
     1              NRESEV,IJAC,NONNEG,NSPAR,IDSYS,IDFDP,IRELAX,LUN0

      PASSD5=.FALSE.
C
C     The vectors y and dy/dt at t=t0 are ready in PHI on entry.
C     Begin by predicting y1 and y'1 with the current h; put these values
C     in Y and YPRIME.
C

C     Modified by Qing, Tang. 4/99
C     ----------------------------
      Hbest=-1.d0
C     ----------------------------
	     
      DO 50 K=1,100
        CALL DCOPY(NEQ,PHI,1,Y,1)
C       Compute y1 in Y:
        CALL DAXPY(NEQ,H,PHI(1,2),1,Y,1)    
C       Check y1 for feasibility; adjust to shorter H if necessary.
        t1=t0 + H
        IF(LYSTOP.OR.NONNEG.GT.0) THEN
          CALL DDFEAS(NSTVAR,NEQ,t1,H,HFEAS,HMIN,PHI,Y,
     1    IYstop,Cstop,LYSTOP,DOT9,1,IDUM,IWM,IDID)
          IF(IDID.LT.0) RETURN
          IF(ABS(HFEAS).LT.ABS(H)) HUP=MIN(HUP,ABS(H))
          H=HFEAS
          t1=t0 + HFEAS    
        END IF
C       Compute y'1 in YPRIME:
        IF(ABS(dt).GT.ZERO) dt=MAX(HMIN,TENM3*ABS(H))
        ITIME=1
        CALL DDAYPR(ITIME,t1,Y,YPRIME,NSTVAR,NEQ,H,dt,WT,NWTS,IDID,
     1  RPAR,IPAR,PHI,PD2,DELTA,ERLOC,ROMAX,PD,IWM,UROUND,UDEL,EWK,
     2  IEROW,DTEM,Ieform,fsub,Esub,Jac,Bsub)
        IF(IDID.LT.0) RETURN
C
C       Compute weight vector WT at current y1.
C
        CALL DDAWTS(NWTS,INFO2,RTOL,ATOL,Y,WT,IDID,LUN0)
        IF(IDID.LT.0) RETURN
        DO 10 I=1,NWTS
          ERLOC(I)=YPRIME(I)-PHI(I,2)
10      CONTINUE
        DIFNRM=DDANRM(NWTS,ERLOC,WT,UROUND,RMAX,INORM)
        IF(.NOT.PASSD5) THEN
          IF(ABS(H)*DIFNRM.GT.TWO) HUP=MIN(HUP,ABS(H))
          YPN0=DDANRM(NWTS,PHI(1,2),WT,UROUND,RMAX,INORM)
          YPN1=DDANRM(NWTS,YPRIME,WT,UROUND,RMAX,INORM)
          IF(ABS(H)*(YPN0+YPN1).GT.TWO) THEN
            H=TENTH*H
            IF(ABS(H).LE.HMIN) THEN
              IF(LOG.EQ.1) WRITE(LUN0,20) HMIN
20            FORMAT(' THE TOLERANCES REQUIRE AN INITIAL STEP',
     1             ' SHORTER THAN HMIN=',1P,D12.3)
              IDID=-33
              RETURN
            END IF
          ELSE
            PASSD5=.TRUE.
          END IF
        END IF

        IF(PASSD5) THEN
	
C       Modified by Qing, Tang. 4/99
C       -- Store the larger trial H in Hbest  
 	  Hbest=MAX(Hbest,H)
C       ------------------------------------

          RALPHA=SQRT(DIFNRM*ABS(H))
          IF(RALPHA .GT. TENM2) THEN
            ALPHA=ONE/RALPHA
          ELSE
            ALPHA=HUNDRD
          END IF
          IF(ALPHA.GT.ONEP2) THEN
            HSIZE=MIN(ALPHA*ABS(H),HALF*(ABS(HUP)+ABS(H)))
            H=SIGN(HSIZE,H)
          ELSE IF(ALPHA.LT.ONE) THEN
            HUP=ABS(H)
            H=MAX(TENTH,ALPHA)*H
            IF(ABS(H).LE.HMIN) THEN
              IF(LOG.EQ.1) WRITE(LUN0,30) K,H,HMIN
30            FORMAT(1X,' STEP INITIALIZATION FAILED IN CYCLE',I4,/,
     1          ' WITH H=',1P,D12.4,' VS. HMIN=',1PD12.4) 
              IDID=-33
              RETURN
            END IF
          ELSE
            H=DOT8*H
c            WRITE(LUN0,40) H,K
40          FORMAT(1X,/,' Initial h =',1P,D9.2,' found in',I3,
     1             ' modified Shampine iterations.',/,1X)
            RETURN
          END IF
        END IF
        N=K
50    CONTINUE
      
C     Modified by Qing, Tang. 4/99
C     ----------------------------
      IF(Hbest.GT.Hmin)THEN
	H=Hbest
	  if( n_warn < lim_warn ) then
	     n_warn = n_warn+1
           IF(LOG.EQ.1) WRITE(LUN0,70)H
	  endif
70      FORMAT(' SHAMPINE FAILED AND CHOOSE THE BEST H ', D9.2,'
     + AS THE INITIAL H')
	RETURN
      ELSE
	  if( n_warn < lim_warn ) then
	     n_warn = n_warn+1
           IF(LOG.EQ.1) WRITE(LUN0,*)'SHAMPINE FAILED AND Hbest < HMIN'
	  endif
        IDID=-33	
        RETURN
      END IF

C      WRITE(LUN0,60) N
C60    FORMAT(' NO INITIAL  H  FOUND IN',I4,' SHAMPINE ITERATIONS.')
C      IDID=-33
C      RETURN
C     ----------------------------


C
C     End of Subroutine DFINDH.
C
      END
      SUBROUTINE DDFEAS(NSTVAR,N2,t,H,HFEAS,HMIN,YB,YH,
     1  IYstop,Cstop,LYSTOP,FACTOR,IYFIX,M,IWM,IDID)
C
C     This subroutine tests YB and YH, the Y-vectors before and after a
C     proposed solution increment, whenever nonnegative state elements are
C     requested or a stopping value Cstop of a state or sensitivity variable
C     Y(IYstop) needs watching. A move, HFEAS, is computed for which Cstop
C     is not crossed with LYSTOP=.TRUE., and/or for which the nonnegativity 
C     is maintained when requested. The argument FACTOR is the permitted 
C     fractional move toward the governing limit. 
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LYSTOP
      DIMENSION YB(*), YH(*), IWM(*)
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, LNST=11)
      COMMON/DDA002/IBAND,ML,MU,MBAND,LDBAND,LDBM1,LENEWK,NALG,LENPD,
     1              NRESEV,IJAC,NONNEG,NSPAR,IDSYS,IDFDP,IRELAX,LUN0
C
      HFEAS=H
      IF(NONNEG.NE.0) THEN
        DO 130 I=1,NSTVAR
          IF(YB(I).LT.ZERO) THEN
            IF(LOG.EQ.1) WRITE(LUN0,120) I,YB(I),t
120   FORMAT(' NEGATIVE Y(',I6,')=',1P,D12.3,' AT t=',
     1  1P,D12.5,/,' VIOLATES THE SPECIFICATION INFO(10) = 1 .')
            IDID=-33
          END IF
130     CONTINUE
        IF(IYFIX.LT.0.OR.IDID.EQ.-33) RETURN
C       Shorten the step when necessary to maintain nonnegativity.
        DO 140 I=1,NSTVAR
          IF(YH(I).LT.(ONE-FACTOR)*YB(I)) THEN
            TEST= FACTOR*H*(YB(I)/(YB(I)-YH(I)))
            IF(ABS(TEST).LT.ABS(HFEAS)) THEN
              ILIM=I
              HFEAS=TEST
            END IF
          END IF
140     CONTINUE
        IF(ABS(HFEAS).LE.HMIN) THEN
          HLIM=HFEAS/FACTOR
          HMF=ABS(HLIM)/HMIN
          IF(IYFIX.LE.1) THEN
             IF(LOG.EQ.1) WRITE(LUN0,150) ILIM, HMF, ABS(HLIM)
          END IF
150   FORMAT(1X,/,' Y(',I6,') IS NEGATIVE FOR ANY STEPSIZE ABS(H) ',
     1 '.GE.',1P,D9.2,'*HMIN =',1P,D9.2)
          IF(IYFIX.EQ.2) THEN
             IF(LOG.EQ.1) WRITE(LUN0,160) ILIM, HMF, M, IWM(LNST)
          END IF
160   FORMAT(1X,/,' Y(',I6,') IS NEGATIVE FOR ANY CORRECTION .GE.',
     1  1P,D9.2,' OF NEWTON STEP',I4,/,' IN INTEGRATION STEP',I5,'.')
          IDID=-33
          RETURN
        END IF
      END IF
      IF(LYSTOP) THEN
        TEST=(YH(IYstop)-Cstop)*SIGN(ONE,YB(IYstop)-Cstop)
        IF(TEST.LT.ZERO) THEN
          TEMP=(Cstop-YB(IYstop))/(YH(IYstop)-YB(IYstop))
          TEMP=H*TEMP
          IF(ABS(TEMP).LT.ABS(HFEAS)) HFEAS = TEMP
        END IF
      END IF
      IF(IYFIX.GE.1.AND.ABS(HFEAS).LT.ABS(H)) THEN
C
C       Compute new YH for feasible step:
C      
        W2=HFEAS/H
        W1=ONE-W2
        DO 170 I=1,NSTVAR
          YH(I)=MAX(ZERO,W1*YB(I)+W2*YH(I))
170     CONTINUE
        IF(N2.GT.NSTVAR) THEN
          DO 180 I = NSTVAR+1,N2
            YH(I) = W1*YB(I)+W2*YH(I)
180       CONTINUE
        END IF
      END IF
      RETURN
C     
C     End of Subroutine DDFEAS
C
      END
      SUBROUTINE DDAJAC (JOB,t,dt,NSTVAR,U0,U,UPRIME,DELTA,CJ0,H,
     1           DTEM,RHS2,IER,WT,ROMAX,PD,IWM,IRES,UROUND,UDEL,LRANK,
     2           RPAR,IPAR,EWK,IEROW,Ieform,fsub,Esub,Jac)
C
C     This routine computes one of the following matrices, according to 
C     the JOB specification:
C
C     JOB=0:  PD = partial { E[u-u0_] - fA } / partial u^T
C     JOB=-1: PD =                   E        - partial f_A / partial u^T
C     Job=2:  PD = partial { Eu' - fD - partial fA / partial t } / partial u^T
C     Job=3:  PD = CJ*E + partial { Eu' - f } / over partial u^T
C
C     Row scaling and LU factorization of PD are then performed, except 
C     for JOB=2. In that case, PD denotes the array PD2 
C     which is used without factorization in Subroutine DDAYPR. 
C
C     Finite-difference evaluation of PD is available for each JOB type.
C     User-provided derivatives can be used for JOB=0, -1, and 3  if
C     partial E / partial u^T = 0, and for JOB=2 if additionally 
C     partial**2 fA/(partial theta)(partial t) = 0 at the initial state
C     (as will always be true if E is nonsingular there).
C     Checking these validity conditions is the user's responsibility.
C     PD is computed as a dense matrix of order NSTVAR if IBAND=0, 
C     or a band matrix with leading dimension LDBAND if IBAND=1.
C
	use ci_sens, only: sens_store

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION U0(*),U(*),UPRIME(*),DELTA(*),DTEM(*),RHS2(*),WT(*)
      DIMENSION ROMAX(*),PD(*),IWM(*),RPAR(*),IPAR(*),EWK(*),IEROW(*)
      PARAMETER (LNRE=12, LIPVT=24)
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, CHSIGN=-1.0D0, TWO=2.0D0)
      EXTERNAL  fsub, Esub, Jac
      COMMON/DDA002/IBAND,ML,MU,MBAND,LDBAND,LDBM1,LENEWK,NALG,LENPD,
     1              NRESEV,IJAC,NONNEG,NSPAR,IDSYS,IDFDP,IRELAX,LUN0
      LRANK=NSTVAR
      IER = 0
      SQUR = SQRT(UROUND)
      SQUR = MAX(SQUR,UDEL)
      IF(IBAND.EQ.1) THEN
        MBA=MIN(MBAND,NSTVAR)
        MSAVE=(NSTVAR/MBAND)+1
        ISAVE=LENPD 
        IPSAVE=ISAVE+MSAVE
      END IF
      DO 10 I=1,LENPD
        PD(I)=ZERO
10    CONTINUE
      CJ=CJ0
      IRES=0
      IF(IBAND.EQ.0) THEN
        IF(IJAC.EQ.1) THEN
C         Dense matrix, including user-provided derivatives:
          IF(JOB.GT.0.OR.NALG.GT.0) THEN
            CALL Jac(t,NSTVAR,U,PD,RPAR,IPAR,Ieform,IRES)
            IF(IRES.NE.1) THEN
              IF(LOG.EQ.1) WRITE(LUN0,15) IRES
15            FORMAT(1X/' IRES=',I3,' RETURNED BY USER SUBROUTINE JAC',
     1        ' INDICATES AN UNSUCCESSFUL CALCULATION.')
              IRES=-2
              RETURN
            END IF
            CALL DSCAL(LENPD,CHSIGN,PD,1)
            IF(JOB.LE.0.AND.NALG.LT.NSTVAR) THEN
              DO 17 I=1,NSTVAR
                IF(IEROW(I).NE.0) CALL DSCAL(NSTVAR,ZERO,PD(I),NSTVAR)
17            CONTINUE
            ELSE IF(JOB.EQ.2.AND.NALG.GT.0) THEN
              DO 18 I=1,NSTVAR
                IF(IEROW(I).EQ.0) CALL DSCAL(NSTVAR,ZERO,PD(I),NSTVAR)
18            CONTINUE
            END IF
          END IF
          IF(JOB.NE.3) CJ=ONE
          IF(IDSYS.EQ.0) THEN
            DO 20 I=1,LENPD,NSTVAR+1
              PD(I)=PD(I)+CJ
20          CONTINUE
          ELSE IF(ABS(IDSYS).EQ.1) THEN
            INCPD=NSTVAR+1
            CALL DAXPY(NSTVAR,CJ,EWK,1,PD,INCPD)
          ELSE IF(ABS(IDSYS).EQ.2) THEN
            CALL DAXPY(LENPD,CJ,EWK,1,PD,1)
          END IF
        ELSE
C         Dense finite-difference-generated matrix.
          IRES=0
          NROW=0
          DO 60 I=1,NSTVAR
            DEL=SQUR*MAX(ABS(U(I)),ABS(H*UPRIME(I)),ABS(WT(I)))
            DEL=SIGN(DEL,U(I))
            DEL=(U(I)+DEL)-U(I)
            USAVE=U(I)
            UPSAVE=UPRIME(I)
            U(I)=U(I)+DEL
            UPRIME(I)=UPRIME(I)+CJ*DEL
            CALL DDARES (JOB,t,dt,NSTVAR,U0,U,UPRIME,DTEM,RHS2,IRES,
     1                RPAR,IPAR,EWK,IEROW,PD,IWM,Ieform,fsub,Esub)
            IF (IRES .LT. 0) RETURN
            U(I)=USAVE
            UPRIME(I)=UPSAVE
            DELINV=ONE/DEL
            DO 50 L=1,NSTVAR
              PD(NROW+L)=(DTEM(L)-DELTA(L))*DELINV
50          CONTINUE
            NROW=NROW+NSTVAR
60        CONTINUE
        END IF
        IF(JOB.EQ.2) RETURN  
	  
	  if( job == 3 ) call sens_store( nstvar, t, u(1:nstvar), cj,
     1							 pd(1:lenpd), lenpd )    
C
C       Row scaling of dense iteration matrix:
C
        DO 70 I=1,NSTVAR
          INDMAX=(IDAMAX(NSTVAR,PD(I),NSTVAR)-1)*NSTVAR + I
          DENOM=ABS(PD(INDMAX))
          IF(DENOM.LE.ZERO) THEN
            IF(LOG.EQ.1) WRITE(LUN0,65) I, t
65          FORMAT(1X,/,' FATAL ERROR: ROW ',I7,' OF ITERATION MATRIX',
     1      '  PD IS ZERO AT t =',1P,D12.4,/,1X)
            IRES=-2
            RETURN  
          ELSE
            ROMAX(I)=DENOM
            ALPHA=ONE/DENOM
            CALL DSCAL(NSTVAR,ALPHA,PD(I),NSTVAR)
          END IF
70      CONTINUE
C       LU decomposition of scaled dense iteration matrix.
        CALL DGETRF(NSTVAR,NSTVAR,PD,NSTVAR,IWM(LIPVT),IER)
        IF(JOB.EQ.0.AND.IER.GE.0.AND.IJAC.EQ.0) THEN
C         Compute effective rank of finite-difference-generated PD, 
C         based on rounding precision of pivotal difference quotients:
          LRANK=0
          IF(IER.EQ.0) LIM=NSTVAR
          IF(IER.GT.0) LIM=IER-1
          DO 80 J=1,LIM
            LOC=IWM(LIPVT+J-1)
            DEL=SQUR*MAX(ABS(U(J)),ABS(WT(J)))         
            DEL=(U(J)+DEL)-U(J)
            ENUM=UROUND*ABS(DELTA(LOC))*TWO
            IF(ENUM.GT.DEL*ABS(PD((J-1)*NSTVAR+J))) RETURN
            LRANK=LRANK+1
80        CONTINUE
        ELSE
          IF(IER.GT.0) LRANK=IER-1
        END IF
        RETURN
      ELSE
        IF(IJAC.EQ.1) THEN
C         Band matrix, including user-provided derivatives:
          IF(JOB.GT.0.OR.NALG.GT.0) THEN
            CALL Jac(t,NSTVAR,U,PD,RPAR,IPAR,Ieform,IRES)
            IF(IRES.NE.1) THEN
              IF(LOG.EQ.1) WRITE(LUN0,15) IRES
              IRES=-2
              RETURN
            END IF
            CALL DSCAL(LENPD,CHSIGN,PD,1)
            IF(JOB.LE.0.AND.NALG.LT.NSTVAR) THEN
              DO 90 I=1,NSTVAR
                IF(IEROW(I).NE.0) THEN
                  J1=MAX(1,I-ML)
                  J2=MIN(NSTVAR,I+MU)
                  N=J2-J1+1
                  INDX1=I-J1+MBAND+(J1-1)*LDBAND
                  CALL DSCAL(N,ZERO,PD(INDX1),LDBM1)
                END IF
90            CONTINUE
            ELSE IF(JOB.EQ.2.AND.NALG.GT.ZERO) THEN
              DO 100 I=1,NSTVAR
                IF(IEROW(I).EQ.0) THEN
                  J1=MAX(1,I-ML)
                  J2=MIN(NSTVAR,I+MU)
                  N=J2-J1+1
                  INDX1=I-J1+MBAND+(J1-1)*LDBAND
                  CALL DSCAL(N,ZERO,PD(INDX1),LDBM1)
                END IF
100           CONTINUE
            END IF                
          END IF
          IF(JOB.NE.3) CJ=ONE
          IF(IDSYS.EQ.0) THEN
            INDX=MBAND
            DO 120 I=1,NSTVAR
              PD(INDX)=PD(INDX)+CJ
              INDX=INDX+LDBAND
120         CONTINUE
          ELSE IF(ABS(IDSYS).EQ.1) THEN
            CALL DAXPY(NSTVAR,CJ,EWK,1,PD(MBAND),LDBAND)
          ELSE IF(ABS(IDSYS).EQ.2) THEN
            CALL DAXPY(LENPD,CJ,EWK,1,PD,1)
          END IF
        ELSE
C         Band finite-difference-generated matrix.
          DO 250 J=1,MBA
            DO 220 N=J,NSTVAR,MBAND
              K= (N-J)/MBAND + 1
              PD(ISAVE+K)=U(N)
              PD(IPSAVE+K)=UPRIME(N)
              DEL=SQUR*MAX(ABS(U(N)),ABS(H*UPRIME(N)),ABS(WT(N)))
              DEL=SIGN(DEL,U(N))
              DEL=(U(N)+DEL)-U(N)
              U(N)=U(N)+DEL
              UPRIME(N)=UPRIME(N)+CJ*DEL
220         CONTINUE
            CALL DDARES (JOB,t,dt,NSTVAR,U0,U,UPRIME,DTEM,RHS2,IRES,
     1               RPAR,IPAR,EWK,IEROW,PD,IWM,Ieform,fsub,Esub)
            IF(IRES.LT.0) RETURN
            DO 240 N=J,NSTVAR,MBAND
              K= (N-J)/MBAND + 1
              U(N)=PD(ISAVE+K)
              UPRIME(N)=PD(IPSAVE+K)
              DEL=SQUR*MAX(ABS(U(N)),ABS(H*UPRIME(N)),ABS(WT(N)))
              DEL=SIGN(DEL,U(N))
              DEL=(U(N)+DEL)-U(N)
              DELINV=ONE/DEL
              I1=MAX(1,(N-MU))
              I2=MIN(NSTVAR,(N+ML))
              II=N*LDBM1-ML
              DO 230 I=I1,I2
                PD(II+I)=(DTEM(I)-DELTA(I))*DELINV
230           CONTINUE
240         CONTINUE
250       CONTINUE
        END IF
      END IF
      IF(JOB.EQ.2) RETURN
C
C     Row scaling of band iteration matrix:
C
      DO 300 I=1,NSTVAR
        J1=MAX(1,I-ML)
        J2=MIN(NSTVAR,I+MU)
        N=J2-J1+1
        INDX1=I-J1+MBAND + (J1-1)*LDBAND
        INDMAX=(IDAMAX(N,PD(INDX1),LDBM1)-1)*LDBM1+INDX1
        DENOM=ABS(PD(INDMAX))
        IF(DENOM.LE.ZERO) THEN
          IF(LOG.EQ.1) WRITE(LUN0,65) I, t
          IRES=-2
          RETURN  
        ELSE
          ROMAX(I)=DENOM
          ALPHA=ONE/DENOM
          CALL DSCAL(N,ALPHA,PD(INDX1),LDBM1)
        END IF
300   CONTINUE        
C
C     LU decomposition of scaled band iteration matrix:
C
      CALL DGBTRF(NSTVAR,NSTVAR,ML,MU,PD,LDBAND,IWM(LIPVT),IER)
      IF(JOB.EQ.0.AND.IER.GE.0.AND.IJAC.EQ.0) THEN
C       Compute effective rank of finite-difference-generated PD, 
C       based on rounding precision of pivotal difference quotients:
        KPIV=MBAND
        LRANK=0
        IF(IER.EQ.0) LIM=NSTVAR
        IF(IER.GT.0) LIM=IER-1
        DO 320 J=1,LIM
          LOC=IWM(LIPVT+J-1)
          DEL=SQUR*MAX(ABS(U(J)),ABS(WT(J)))         
          DEL=(U(J)+DEL)-U(J)
          ENUM=UROUND*ABS(DELTA(LOC))*TWO
          IF(ENUM.GT.DEL*ABS(PD(KPIV))) RETURN
          LRANK=LRANK+1
          KPIV=KPIV+LDBAND
320     CONTINUE
      ELSE
        IF(IER.GT.0) LRANK=IER-1
      END IF
      RETURN
C
C     End of Subroutine DDAJAC
C
      END
      SUBROUTINE DDSENS (JOB,t,dt,NSTVAR,Y0,Y,YPRIME,ERLOC,DELTA,DTEM,
     1           RHS2,RPAR,IPAR,EWK,IEROW,ROMAX,PD,IWM,CJ,UROUND,UDEL,
     2           IRES,Ieform,fsub,Esub,Bsub)
C
C     This routine controls the solution of the sensitivity equations.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y0(NSTVAR,*),Y(NSTVAR,*),YPRIME(NSTVAR,*)
      DIMENSION ERLOC(NSTVAR,*),DELTA(NSTVAR,*),DTEM(*),RHS2(*)
      DIMENSION ROMAX(*),PD(*),RPAR(*),IPAR(*),EWK(*),IEROW(*),IWM(*)
      EXTERNAL  fsub, Esub, Bsub
      COMMON/DDA002/IBAND,ML,MU,MBAND,LDBAND,LDBM1,LENEWK,NALG,LENPD,
     1              NRESEV,IJAC,NONNEG,NSPAR,IDSYS,IDFDP,IRELAX,LUN0
      PARAMETER (ZERO=0.0D0, CHSIGN=-1.0D0)
C
C     Clear the sensitivity elements of DELTA. Note that DDSENS is 
C     not called unless NSPAR.GE.1 .
C
      DO 10 J=2,NSPAR+1
        DO 10 I=1,NSTVAR
          DELTA(I,J)=ZERO
10    CONTINUE
C
C     Indicate if DELTA(1,...NSTVAR) is an updated residual vector.
C
      IF(IJAC.EQ.0.OR.JOB.EQ.2) THEN
        IDEL0=1
      ELSE
        IDEL0=0
      END IF
C
C     Compute the partial derivatives of the right-hand functions 
C     with respect to the parameter vector p. Store these values
C     in DELTA(*,IPARM+1). 
C
      DO 20 IPARM=1,NSPAR
        CALL DDFDPJ (JOB,t,dt,NSTVAR,Y0,Y,YPRIME,IPARM,DELTA(1,1),
     1  DELTA(1,IPARM+1),DTEM,RHS2,RPAR,IPAR,EWK,IEROW,PD,IWM,UROUND,
     2  UDEL,IRES,Ieform,fsub,Esub,Bsub,IDEL0)
        IF(IRES.LT.0) RETURN
        IF(JOB.EQ.3) CALL DDSRHS (Y(1,IPARM+1),YPRIME(1,IPARM+1),
     1               DELTA(1,IPARM+1),EWK,NSTVAR,CJ)
20    CONTINUE
      IF(JOB.EQ.2) RETURN
C
C     Now use the LU decomposition of the iteration
C     matrix to solve the sensitivity equations.
C     Store the solution elements in DELTA(i,IPARM+1).
C
      IF(IDSYS.NE.0.OR.JOB.GT.2) THEN
        DO 70 IPARM=1,NSPAR       
          CALL DDASLV (NSTVAR,NSTVAR,DELTA(1,IPARM+1),
     1                             DELTA(1,IPARM+1),ROMAX,PD,IWM,IER)
          IF(IER.NE.0)RETURN
70      CONTINUE
      END IF
      IF(JOB.LT.3) RETURN
C
C     The remaining calculations are for corrector formulas in DDASTP.
C     Update the sensitivity solution and error vector and return.
C
      DO 80 I=1,NSTVAR
        DO 80 J=2,NSPAR+1
          ERLOC(I,J)=-Y(I,J)+DELTA(I,J)
          Y(I,J)=DELTA(I,J)
          YPRIME(I,J)=YPRIME(I,J)+CJ*ERLOC(I,J)
80    CONTINUE
      RETURN
C
C     End of Subroutine DDSENS
C
      END
      SUBROUTINE DDFDPJ (JOB,t,dt,NSTVAR,Y0,Y,YPRIME,IPARM,DEL0,
     1           DELTA,DTEM,RHS2,RPAR,IPAR,EWK,IEROW,PD,IWM,UROUND,
     2           UDEL,IRES,Ieform,fsub,Esub,Bsub,IDEL0)
C
C     This routine computes partial derivatives of the right-hand 
C     functions with respect to parameter number IPARM, all except
C     the predictor terms which are done for JOB=3 in Subroutine DDSRHS.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y0(*),Y(*),YPRIME(*),DTEM(*),RHS2(*)
      DIMENSION DEL0(*),DELTA(*),IWM(*),RPAR(*),IPAR(*)
      DIMENSION EWK(*),IEROW(*),PD(*)
      EXTERNAL  fsub, Esub, Bsub
      COMMON/DDA002/IBAND,ML,MU,MBAND,LDBAND,LDBM1,LENEWK,NALG,LENPD,
     1              NRESEV,IJAC,NONNEG,NSPAR,IDSYS,IDFDP,IRELAX,LUN0
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, CHSIGN=-1.0D0, LNRE=12)
      IW=NSTVAR*IPARM      
      IF(IDFDP.EQ.1.OR.(JOB.EQ.0.AND.NALG.EQ.0)) THEN
C       Here the function sensitivities are available directly.
        IF(JOB.EQ.0) THEN
          IF(NALG.GT.0) THEN
C           Compute partial (f vector)/partial (Sensitivity parameter IPARM) 
C           using subroutine Bsub. Accept if IRES comes back reset to 1; 
C           otherwise use a divided difference of (-Residual).
            IRES=0
            CALL Bsub(t,NSTVAR,Y,DELTA,IPARM,RPAR,IPAR,Ieform,IRES)
            IF(IRES.NE.1) GO TO 50
            DO 10 I=1,NSTVAR
              IF(IEROW(I).NE.0) DELTA(I)=ZERO
10          CONTINUE
          END IF
C         Add the contributions of the sensitivity values at t0_:
          IF(IDSYS.EQ.0) THEN
C           E is a unit matrix.
            DO 20 I=1,NSTVAR
              DELTA(I)=DELTA(I) + Y(IW+I)
20          CONTINUE
          ELSE IF(ABS(IDSYS).EQ.1) THEN
C           E is a diagonal matrix.
            DO 30 I=1,NSTVAR
              IF(IEROW(I).NE.0) DELTA(I)=DELTA(I) + EWK(I)*Y(IW+I)
30          CONTINUE
          ELSE IF(ABS(IDSYS).EQ.2) THEN
            IF(IBAND.EQ.0) THEN
C             E is a dense matrix.
              CALL DGEMV('N',NSTVAR,NSTVAR,
     1                       ONE,EWK,NSTVAR,Y(IW+1),1,ONE,DELTA,1)
            ELSE
C             E is a band matrix.
              CALL DGBMV('N',NSTVAR,NSTVAR,ML,MU,
     1                 ONE,EWK(1+ML),LDBAND,Y(IW+1),1,ONE,DELTA,1)
            END IF
          END IF
        ELSE
          IRES=0
          CALL Bsub(t,NSTVAR,Y,DELTA,IPARM,RPAR,IPAR,Ieform,IRES)
          IF(IRES.NE.1) GO TO 50
          IF(JOB.EQ.2.AND.NALG.GT.0) THEN
C           The dFA/dtheta elements are inappropriate here, and are 
C           replaced by zeros, thus neglecting d**2 FA/(dtheta)(dt) .
            DO 40 I=1,NSTVAR
              IF(IEROW(I).EQ.0) DELTA(I)=ZERO
40          CONTINUE
          END IF
        END IF
C       Direct sensitivity calculations completed.
        RETURN          
      END IF
50    CONTINUE
C     Set reference vector DEL0 when needed for differencing:
      IF(IDEL0.EQ.0) THEN
        IRES=0
        IWM(LNRE)=IWM(LNRE)+1
        CALL DDARES (JOB,t,dt,NSTVAR,Y0,Y,YPRIME,DEL0,RHS2,IRES,
     1    RPAR,IPAR,EWK,IEROW,PD,IWM,Ieform,fsub,Esub)
        IF(IRES.LT.0) RETURN
        IDEL0=1
      END IF
C     Compute a perturbation step and perturb the parameter
C     RPAR(jpar) to compute the derivatives.
      SQUR=SQRT(UROUND)
      SQUR=MAX(SQUR,UDEL)
      JPAR=IPAR(IPARM)
      DEL=SQUR*(ONE+ABS(RPAR(JPAR)))
      PSAVE=RPAR(JPAR)
      RPAR(JPAR)=RPAR(JPAR)+DEL
      DEL=RPAR(JPAR)-PSAVE
      IRES=0
      IWM(LNRE)=IWM(LNRE)+1
      CALL DDARES (JOB,t,dt,NSTVAR,Y0,Y,YPRIME,DTEM,RHS2,IRES,
     1             RPAR,IPAR,EWK,IEROW,PD,IWM,Ieform,fsub,Esub)
      IF(IRES.LT.0)RETURN
      DELINV=-ONE/DEL
      DO 60 I=1,NSTVAR
        DELTA(I)=(DTEM(I)-DEL0(I))*DELINV
60    CONTINUE
      RPAR(JPAR)=PSAVE
      RETURN
C
C     End of Subroutine DDFDPJ
C
      END
      SUBROUTINE DDSRHS(W,WPRIME,DELTA,EWK,NSTVAR,CJ)
C
C     This subroutine completes the right-hand sides of the 
C     sensitivity BDF equations by adding the predictor terms.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION W(*),WPRIME(*),DELTA(*),EWK(*)
      PARAMETER (CHSIGN=-1.0D0, ONE=1.0D0)
      COMMON/DDA002/IBAND,ML,MU,MBAND,LDBAND,LDBM1,LENEWK,NALG,LENPD,
     1              NRESEV,IJAC,NONNEG,NSPAR,IDSYS,IDFDP,IRELAX,LUN0
      IF(IDSYS.EQ.0) THEN
        DO 10 K=1,NSTVAR
          DELTA(K)=DELTA(K)-(WPRIME(K)-CJ*W(K))
10      CONTINUE
      ELSE IF(ABS(IDSYS).EQ.1) THEN
        DO 20 K=1,NSTVAR
          DELTA(K)=DELTA(K) - EWK(K)*(WPRIME(K)-CJ*W(K))
20      CONTINUE
      ELSE IF(ABS(IDSYS).EQ.2) THEN
        IF(IBAND.EQ.0) THEN
          CALL DGEMV('N',NSTVAR,NSTVAR,
     1                     CHSIGN,EWK,NSTVAR,WPRIME,1,ONE,DELTA,1)
          CALL DGEMV('N',NSTVAR,NSTVAR,
     1                     CJ,EWK,NSTVAR,W,1,ONE,DELTA,1)
        ELSE
          CALL DGBMV('N',NSTVAR,NSTVAR,ML,MU,
     1               CHSIGN,EWK(1+ML),LDBAND,WPRIME,1,ONE,DELTA,1)
          CALL DGBMV('N',NSTVAR,NSTVAR,ML,MU,
     1                   CJ,EWK(1+ML),LDBAND,W,1,ONE,DELTA,1)
        END IF            
      END IF
      RETURN
C
C     End of Subroutine DDSRHS
C
      END
      SUBROUTINE DDARES (JOB,t,dt,NSTVAR,U0,U,UPRIME,DELTA,RHS2,IRES,
     1      RPAR,IPAR,EWK,IEROW,PD,IWM,Ieform,fsub,Esub)
C
C     DDARES computes residuals or right-hand terms, according to JOB.
C     JOB=0 calls for residuals of the initial equations for u:
C           DELTA = E[u-u0] - fA
C     Here fA is a vector of right-hand function values for the algebraic
C     equations, and zeros for the differential equations. fD, used below,
C     is the right-hand vector  f  minus fA. If E is nonsingular, fA is
C     empty and is not computed. An error return is given if the nonzero
C     rows of the initial E are linearly dependent, since this indicates
C     an ill-posed problem.
C     JOB=-1 is like JOB=0 except that E is left at its previous value;
C            Subroutine Esub is not called.
C     JOB=1 calls for the vector
C           DELTA = -fD - partial fA/partial t
C     which appears in the equation for the initial du/dt, and 
C     in u-differenced form in the equation for the initial dW/dt.
C     JOB=2 calls for the vector
C           DELTA = Edu/dt - fD - partial fA/partial t
C     whose theta-derivatives appear in the equation for the initial dW/dt.
C     JOB=3 calls for the corrector residuals 
C           DELTA = Edu/dt - f
C     for the integrator DDASTP. 
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION U0(*),U(*),UPRIME(*),DELTA(*)
      DIMENSION RHS2(*),RPAR(*),IPAR(*),EWK(*),IEROW(*),PD(*),IWM(*)
      COMMON/DDA002/IBAND,ML,MU,MBAND,LDBAND,LDBM1,LENEWK,NALG,LENPD,
     1              NRESEV,IJAC,NONNEG,NSPAR,IDSYS,IDFDP,IRELAX,LUN0
      PARAMETER (ZERO=0.0D0, HALF=0.5D0, ONE=1.0D0, CHSIGN=-1.0D0)
      PARAMETER (LNRE=12, LIPVT=24)
      EXTERNAL  fsub, Esub
C     NRESEV is the number of times this subroutine has been called.
      NRESEV=NRESEV+1
C
C     Default setting of NALG (number of zero rows in the global matrix E).
C
      IF(IDSYS.EQ.0) NALG=0
C
C     Initialization or update of work array EWK for global user-defined 
C     matrix E. Test the rank of the initial E and set the matrix 
C     IEROW (needed only when IDSYS.NE.0), whose nonzero elements denote 
C     nonzero rows in the initial E. Also initialize NALG.
C
      IF((JOB.GE.0.AND.IDSYS.GT.0).OR.(NRESEV.EQ.1.AND.IDSYS.NE.0))
     1  THEN
        DO 20 I=1,LENEWK
          EWK(I)=ZERO
20      CONTINUE
C       Update EWK with the user's subroutine Esub.
        CALL Esub(t,NSTVAR,U,EWK,RPAR,IPAR,Ieform,IRES)
        IF(IRES.LT.0) RETURN
        IF(NRESEV.EQ.1) THEN
          NDROWS=0
C         On the first call, set the Ierow vector.
          IF(ABS(IDSYS).EQ.1) THEN
C           E is a general diagonal matrix.
            DO 30 I=1,NSTVAR
              IF(ABS(EWK(I)).GT.ZERO) THEN
                IEROW(I)=1
                NDROWS=NDROWS+1
              END IF
30          CONTINUE
          ELSE IF(ABS(IDSYS).EQ.2) THEN
            CALL DCOPY(LENEWK,EWK,1,PD,1)
            NRANK=0
            IF(IBAND.EQ.0) THEN
C             Label and count nonzero rows of matrix E.
              DO 40 I=1,NSTVAR
                INDX=I+NSTVAR*(IDAMAX(NSTVAR,EWK(I),NSTVAR)-1)
                IF(EWK(INDX).GT.ZERO) THEN
                  IEROW(I)=1
                  NDROWS=NDROWS+1
                END IF
40            CONTINUE
C             Compute rank of dense matrix PD.
              I1=1
              J1=1
              L1=1
              M=NSTVAR
              N=NSTVAR
50            CALL DGETRF(M,N,PD(L1),NSTVAR,IWM(LIPVT),INFO)
              IF(INFO.GT.0) THEN
                INC=INFO-1
                NRANK = NRANK + INC
                I1=I1+INC
                J1=J1+INFO
                L1=I1+(J1-1)*NSTVAR
                M=NSTVAR+1-I1
                N=NSTVAR+1-J1
                IF(J1.LE.NSTVAR) GO TO 50
              ELSE
                NRANK=NRANK+N
              END IF
            ELSE
C             Label and count nonzero rows of E.
              DO 60 I=1,NSTVAR
                J1=MAX(1,I-ML)
                J2=MIN(NSTVAR,I+MU)
                N=J2-J1+1
                INDX1=I-J1+MBAND + (J1-1)*LDBAND
                KBIG=INDX1 + (IDAMAX(N,EWK(INDX1),LDBM1) - 1)*LDBM1
                IF(ABS(EWK(KBIG)).GT.ZERO) THEN
                  IEROW(I)=1
                  NDROWS=NDROWS+1
                END IF
60            CONTINUE            
C           Compute rank of band matrix PD:
              L1=MBAND
              J1=1
              N=NSTVAR
70            CALL DGBTRF
     1             (N,N,ML,MU,PD(L1+1-MBAND),LDBAND,IWM(LIPVT),INFO)
              IF(INFO.GT.0) THEN
                INC=INFO-1
                NRANK = NRANK + INC
                J1=J1+INFO
                L1=MBAND+(J1-1)*LDBAND
                IF(J1.LE.NSTVAR) THEN
                  N=NSTVAR+1-J1
C                 Normalize and combine rows J1-1 and J1:
                  LNGTHO=MIN(MBAND-1,N)
                  LNGTHN=MIN(MBAND,N)
                  IOLD=IDAMAX(LNGTHO,PD(L1-1),LDBM1)
                  DIVOLD=PD(L1-1+(IOLD-1)*LDBM1)
                  INEW=IDAMAX(LNGTHN,PD(L1),LDBM1)
                  DIVNEW=PD(L1+(INEW-1)*LDBM1)
                  ALPHA=ONE
                  BETA=ONE
                  DO 80 I=1,LNGTHO
                    LNEW=L1+(I-1)*LDBM1
                    TOLD=PD(LNEW-1)
                    TNEW=PD(LNEW)
                    IF(ABS(TOLD).GT.ZERO.AND.ABS(TNEW).GT.ZERO) THEN
                      ALPHA=SIGN(ONE/ABS(DIVOLD),TOLD)
                      BETA =SIGN(ONE/ABS(DIVNEW),TNEW)
                      GO TO 90
                    END IF
80                CONTINUE
90                CALL DSCAL(LNGTHN,BETA,PD(L1),LDBM1)
                  CALL DAXPY(LNGTHO,ALPHA,PD(L1-1),LDBM1,PD(L1),LDBM1)
                  GO TO 70
                END IF 
              ELSE
                NRANK=NRANK+N
              END IF
            END IF
C           Report if the rank of E is less than NDROWS:
            IF(NRANK.LT.NDROWS) THEN
              IF(LOG.EQ.1) WRITE(LUN0,100) NDROWS,NRANK
100           FORMAT(1X/' INITIAL MATRIX E HAS',I6,' NONZERO ROWS, ',
     1        'BUT ITS RANK IS ONLY',I6,'.',/,' THE INITIAL Y, ',
     2        'YPRIME or h CANNOT BE COMPLETED. REVIEW THE PROBLEM ',
     3        'FORMULATION.')
              IRES=-2
              RETURN
            END IF
          END IF
          NALG=NSTVAR-NDROWS
        END IF
      END IF
C
C     See how many evaluations of the right-hand vector f are needed.
C
      IF(JOB.EQ.3.OR.JOB.LE.0.OR.IDSYS.EQ.0.OR.ABS(dt).LE.ZERO) THEN
        IREP = 0
      ELSE
        IREP = 1
      END IF
C     Start analysis of right-hand function.
      taug=t
110   CONTINUE
      IF(IDSYS.NE.0.OR.JOB.GT.0) THEN
        CALL fsub(taug,NSTVAR,U,DELTA,RPAR,IPAR,Ieform,IRES)
        IF(IRES.LT.0) RETURN
      ELSE
        GO TO 130
      END IF
      IF(IREP.EQ.1) THEN
        CALL DCOPY(NSTVAR,DELTA,1,RHS2,1)
        taug = t + dt
        IREP = 2
        GO TO 110
      ELSE IF(IREP.EQ.2) THEN
        DENOM=taug-t
        DO 120 K=1,NSTVAR
          IF(IEROW(K).EQ.0) THEN
            DELTA(K)=(DELTA(K)-RHS2(K))/DENOM
          ELSE
            DELTA(K)=RHS2(K)
          END IF
120     CONTINUE
      END IF
C     Assemble the residuals according to the JOB specification.
      CALL DSCAL(NSTVAR,CHSIGN,DELTA,1)
130   IF(JOB.LE.0) THEN 
C       DELTA = E[u-u0_] - fA
        IF(IDSYS.NE.0) THEN
          DO 140 I=1,NSTVAR
            IF(IEROW(I).NE.0) DELTA(I)=ZERO
140       CONTINUE
        END IF
        IF(IDSYS.EQ.0) THEN
C         Unit E-matrix
          CALL DCOPY(NSTVAR,U,1,DELTA,1)
          CALL DAXPY(NSTVAR,CHSIGN,U0,1,DELTA,1)
        ELSE IF(ABS(IDSYS).EQ.1) THEN
C         General diagonal E-matrix
          DO 150 K=1,NSTVAR
            DELTA(K)=EWK(K)*(U(K)-U0(K))+DELTA(K)
150       CONTINUE
        ELSE IF(ABS(IDSYS).EQ.2) THEN          
          IF(IBAND.EQ.0) THEN
C           Dense E-matrix
            CALL DGEMV('N',NSTVAR,NSTVAR,ONE,EWK,NSTVAR,U,1,ONE,
     1                                                   DELTA,1)
            CALL DGEMV('N',NSTVAR,NSTVAR,CHSIGN,EWK,NSTVAR,U0,1,ONE,
     1                                                   DELTA,1)
          ELSE
C           Band E-matrix.
            CALL DGBMV('N',NSTVAR,NSTVAR,ML,MU,
     1                        ONE,EWK(1+ML),LDBAND,U,1,ONE,DELTA,1)
            CALL DGBMV('N',NSTVAR,NSTVAR,ML,MU,
     1                     CHSIGN,EWK(1+ML),LDBAND,U0,1,ONE,DELTA,1)
          END IF
        END IF
      ELSE IF (JOB.EQ.1) THEN
        CONTINUE
      ELSE IF (JOB.GE.2) THEN
        IF(IDSYS.EQ.0) THEN
C         Unit E-matrix.
          CALL DAXPY(NSTVAR,ONE,UPRIME,1,DELTA,1)
        ELSE IF(ABS(IDSYS).EQ.1) THEN
C         Diagonal E-matrix.
          DO 160 K=1,NSTVAR
            DELTA(K)=DELTA(K)+EWK(K)*UPRIME(K)
160       CONTINUE
        ELSE IF(ABS(IDSYS).EQ.2) THEN
          IF(IBAND.EQ.0) THEN
C           Dense E-matrix.
            CALL DGEMV('N',NSTVAR,NSTVAR,ONE,EWK,NSTVAR,
     1                                         UPRIME,1,ONE,DELTA,1)
          ELSE
C           Band E-matrix.
            CALL DGBMV('N',NSTVAR,NSTVAR,ML,MU,ONE,EWK(1+ML),LDBAND,
     1                                         UPRIME,1,ONE,DELTA,1)
          END IF
        END IF
      END IF
      RETURN
C
C     End of Subroutine DDARES
C
      END
      SUBROUTINE DDAWTS(NWTS,IWT,RTOL,ATOL,Y,WT,IDID,LUN0)
C
C     This subroutine sets the error weight vector
C     WT according to WT(i)=RTOL(i)*abs(Y(i))+ATOL(i),
C     i=1,...,n.
C     RTOL and ATOL are scalars if IWT = 0,
C     and vectors if IWT = 1.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION RTOL(*),ATOL(*),Y(*),WT(*)
      PARAMETER (ZERO=0.0D0)
      RTOLI=RTOL(1)
      ATOLI=ATOL(1)
      DO 30 I=1,NWTS
        IF (IWT.EQ.0) GO TO 10
        RTOLI=RTOL(I)
        ATOLI=ATOL(I)
10      WT(I)=RTOLI*ABS(Y(I))+ATOLI
        IF(WT(I).LE.ZERO) THEN
          IF(LOG.EQ.1) WRITE(LUN0,20) I, WT(I)
20        FORMAT(1X,/,' DDASAC STOPPED WITH NONPOSITIVE WEIGHT WT(',
     1    I6,')=',1P,D12.3,' IN SUBROUTINE DDAWTS.')
          IDID=-33
          RETURN
        END IF
30    CONTINUE
      RETURN
C
C     End of Subroutine DDAWTS
C
      END
      SUBROUTINE DDATRP(X,XOUT,YOUT,YPOUT,NEQ,KOLD,PHI,PSI)
C
C     The methods in subroutine DDASTP use polynomials
C     to approximate the solution. DDATRP approximates the
C     solution and its derivative at time XOUT by evaluating
C     one of these polynomials and its derivative there.
C     Information defining this polynomial is passed from
C     DDASTP, so DDATRP cannot be used alone.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION YOUT(*),YPOUT(*)
      DIMENSION PHI(NEQ,*),PSI(*)
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
      KOLDP1=KOLD+1
      TEMP1=XOUT-X
      DO 10 I=1,NEQ
        YOUT(I)=PHI(I,1)
        YPOUT(I)=ZERO
10    CONTINUE
      C=ONE
      D=ZERO
      GAMMA=TEMP1/PSI(1)
      DO 30 J=2,KOLDP1
        D=D*GAMMA+C/PSI(J-1)
        C=C*GAMMA
        GAMMA=(TEMP1+PSI(J-1))/PSI(J)
        DO 20 I=1,NEQ
          YOUT(I)=YOUT(I)+C*PHI(I,J)
          YPOUT(I)=YPOUT(I)+D*PHI(I,J)
20      CONTINUE
30    CONTINUE
      RETURN
C
C     End of Subroutine DDATRP
C
      END
      DOUBLE PRECISION FUNCTION DDANRM(N,V,WT,UROUND,RMAX,INORM)
C
C     This function routine computes a weighted norm of a vector v, 
C     with weights from the vector WT. The max-norm
C        DDANRM=max{V(i)/WT(i), i=1,...,n}
C     is returned if inorm.eq.0, and the Euclidean norm
C        DDANRM=sqrt((1/n)*sum(V(i)/WT(i))**2)
C     is returned if inorm.eq.1 .  If any V(i)/WT(i) exceeds the
C     machine limit RMAX, then RMAX is returned as DDANRM.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, INORM
      DIMENSION V(*), WT(*)
      PARAMETER (ZERO=0.0D0, ONE=1.0D0) 
      QMAX = ZERO
      DO 10 I = 1,N
        IF(WT(I).LE.ONE) THEN
          IF(ABS(V(I)).GE.WT(I)*RMAX) THEN
            DDANRM = RMAX
            RETURN
          ELSE
            QMAX = MAX(ABS(V(I))/WT(I),QMAX)
          END IF
        ELSE
          QMAX = MAX(ABS(V(I))/WT(I),QMAX)
        END IF       
10    CONTINUE
      IF(INORM.EQ.0.OR.QMAX.LE.ZERO) THEN
        DDANRM = QMAX
        RETURN
      END IF
      SMALL = QMAX*SQRT(UROUND/DBLE(N))
      DENOM = SQRT(QMAX*DBLE(N)*SQRT(UROUND))
      SUM = ZERO
      DO 20 I = 1,N
       IF(ABS(V(I)/WT(I)).GE.SMALL)
     1       SUM = SUM + ((V(I)/WT(I))/DENOM)**2
20    CONTINUE
      DDANRM = DENOM*SQRT(SUM/DBLE(N))
      RETURN
C
C     End of Function DDANRM
C
      END
      SUBROUTINE DDASLV(LRANK,NSTVAR,DELTA,DSAVE,ROMAX,PD,IWM,IER)
C
C     This routine computes solution vectors, or Newton correction vectors,
C     for the various equation systems encountered in DDASAC. The LU-factored
C     coefficient matrix is provided in the work array PD, and the list of
C     row interchanges performed in the factorization is provided in the
C     integer work array IWM. Reduced-rank Newton corrections are permitted 
C     in the initial state computation, but full rank (LRANK=NSTVAR) 
C     is demanded ultimately for convergence. 
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DELTA(*),DSAVE(*),ROMAX(*),PD(*),IWM(*)
      PARAMETER (LIPVT=24, ZERO=0D0)
      COMMON/DDA002/IBAND,ML,MU,MBAND,LDBAND,LDBM1,LENEWK,NALG,LENPD,
     1              NRESEV,IJAC,NONNEG,NSPAR,IDSYS,IDFDP,IRELAX,LUN0
C
C     Apply the row scaling to the right-hand vector. Delete the
C     right-hand elements of any equations in which pivots were not found.
C
      IF(LRANK.EQ.NSTVAR) THEN
        DO 10 I=1,NSTVAR
          DELTA(I)=DELTA(I)/ROMAX(I)
10      CONTINUE
      ELSE
        DO 20 I=1,NSTVAR
          DELTA(I)=ZERO
20      CONTINUE
        DO 30 I=1,LRANK
          LOC=IWM(LIPVT+I-1)
          DELTA(LOC)=DSAVE(LOC)/ROMAX(LOC)
30      CONTINUE
      END IF
      IF(IBAND.EQ.0) THEN
C       Solve via LU-factored dense matrix. 
        CALL DGETRS('N',LRANK,1,PD,NSTVAR,IWM(LIPVT),DELTA,NSTVAR,IER)
      ELSE
C       Solve via LU-factored band matrix. 
        CALL DGBTRS('N',LRANK,ML,MU,
     1                         1,PD,LDBAND,IWM(LIPVT),DELTA,NSTVAR,IER)
      END IF
      RETURN
C
C     End of Subroutine DDASLV
C
      END
C
C***********************************************************************
C                   E R R O R  Handling Subroutines                    *
C***********************************************************************
C
      SUBROUTINE XERRDD (MSG, NMES, IERT, NI, I1, I2, NR, R1, R2, LUN)
      INTEGER MSG, NMES, IERT, NI, I1, I2, NR
      INTEGER I, LUN, NWDS
      DOUBLE PRECISION R1, R2
      DIMENSION MSG(NMES)
C
C Subroutine XERRDD, as given here, constitutes
C a simplified version of the slatec error handling package
C written by A. C. Hindmarsh at LLL.  version of January 23, 1980,
C modified by L. R. Petzold, April 1982 and by W. E. Stewart, January 1994.
C This version is in double precision.
C
C All arguments are input arguments.
C
C msg    = the message (Hollerith literal or integer array).
C nmes   = the length of msg (number of characters).
C iert   = the error type..
C          1 means recoverable (control returns to caller).
C          2 means fatal (run is aborted--see note below).
C ni     = number of integers (0, 1, or 2) to be printed with message.
C i1,i2  = integers to be printed, depending on ni.
C nr     = number of reals (0, 1, or 2) to be printed with message.
C r1,r2  = reals to be printed, depending on ni.
C
C note..  this routine is machine-dependent and specialized for use
C in limited context, in the following ways..
C 1. the number of hollerith characters stored per word, denoted
C    by ncpw below, is set in a data statement below.
C 2. the value of nmes is assumed to be at most 60.
C    (multi-line messages are generated by repeated calls.)
C 3. if iert = 2, control passes to the statement   stop
C    to abort the run.  this statement may be machine-dependent.
C 4. r1 and r2 are assumed to be in real and are printed
C    in d21.13 format.
C***************************************************************
C 5. THIS VERSION TAKES THE FOLLOWING INTEGERS FROM THE ARGUMENT LUN
C    OF DDASAC, RATHER THAN FROM A DATA STATEMENT:
C       mesflg = print control flag..  replaced here by LUN.
C                print all messages IF LUN.GT.0
C                no printing        IF LUN.LE.0
C       LUN     = logical unit number for messages.
C
C The following are instructions for installing this routine
C in different machine environments.
C
C For a different number of characters per word, change the
C data statement setting ncpw below.
C Alternatives for various computers are shown in comment
C statements.
C
C For a different run-abort command, change the statement following
C statement 100 at the end.
C
C The following value of ncpw is valid for the cdc-6600 and
C cdc-7600 computers.
C     data ncpw/10/
C The following is valid for the cray-1 computer.
C     data ncpw/8/
C The following is valid for the burroughs 6700 and 7800 computers.
C     data ncpw/6/
C The following is valid for the pdp-10 computer.
C     data ncpw/5/
C The following is valid for the vax computer with 4 bytes per integer,
C and for the ibm-360, ibm-303x, and ibm-43xx computers.
C     data ncpw/4/
C The following is valid for the pdp-11, or vax with 2-byte integers.
C     data ncpw/2/
C
      DIMENSION NFORM(13)
      DATA NFORM(1)/1H(/,NFORM(2)/1H1/,NFORM(3)/1HX/,NFORM(4)/1H,/
      DATA NFORM(7)/1HA/,NFORM(10)/1H,/,NFORM(11)/1HA/,NFORM(13)/1H)/
      DATA NCPW/4/
      IF(LUN.GT.0) THEN
        NCH = MIN(NMES,60)
        NWDS = NCH/NCPW
        CALL D88FMT(2,NWDS,NFORM(5))
        CALL D88FMT(2,NCPW,NFORM(8))
        NREM = NCH - NWDS*NCPW
        IF (NREM .GT. 0) NWDS = NWDS + 1
        IF (NREM .LT. 1) NREM = 1
        CALL D88FMT(1,NREM,NFORM(12))
        IF(LOG.EQ.1) WRITE(LUN,NFORM) (MSG(I),I=1,NWDS)
        IF (NI.EQ.1 .AND. LOG.EQ.1) WRITE (LUN, 10) I1
10      FORMAT(6X,23HIN ABOVE MESSAGE,  I1 =,I10)
        IF (NI.EQ.2 .AND. LOG.EQ.1) WRITE (LUN, 20) I1,I2
20      FORMAT(6X,23HIN ABOVE MESSAGE,  I1 =,I10,3X,4HI2 =,I10)
        IF (NR.EQ.1 .AND. LOG.EQ.1) WRITE (LUN, 30) R1
30      FORMAT(6X,23HIN ABOVE MESSAGE,  R1 =,D21.13)
        IF (NR.EQ.2 .AND. LOG.EQ.1) WRITE (LUN, 40) R1,R2
40      FORMAT(6X,15HIN ABOVE,  R1 =,D21.13,3X,4HR2 =,D21.13)
      END IF
50    IF (IERT .NE. 2) RETURN
      STOP
C
C     End of Subroutine XERRDD
C
      END
      SUBROUTINE D88FMT(N,IVALUE,IFMT)
C
C     D88FMT replaces ifmt(1), ... ,ifmt(n) with the
C     characters corresponding to the n least significant
C     digits of ivalue.
C     Taken from the Bell laboratories port library error handler
C     latest revision ---  7 June 1978.
C
C     Jones R.E., *SLATEC common mathematical library error handling
C     package*, SAND78-1189, Sandia Laboratories, 1978.
C
      DIMENSION IFMT(N),IDIGIT(10)
      DATA IDIGIT(1),IDIGIT(2),IDIGIT(3),IDIGIT(4),IDIGIT(5),
     1     IDIGIT(6),IDIGIT(7),IDIGIT(8),IDIGIT(9),IDIGIT(10)
     2     /1H0,1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9/
      NT = N
      IT = IVALUE
10    IF (NT .EQ. 0) RETURN
      INDEX = MOD(IT,10)
      IFMT(NT) = IDIGIT(INDEX+1)
      IT = IT/10
      NT = NT - 1
      GO TO 10
C
C     End of Subroutine D88FMT
C
      END
