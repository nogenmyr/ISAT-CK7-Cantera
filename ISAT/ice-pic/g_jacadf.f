!  This file contains all of the routines needed for calculating the analytical
! Jacobian matrix, namely: g_cires5, g_cidzdt, and g_cklib

c=======================subroutine  g_cires5============================
C                           DISCLAIMER
C
C   This file was generated on 12/07/04 by the version of
C   ADIFOR compiled on June, 1998.
C
C   ADIFOR was prepared as an account of work sponsored by an
C   agency of the United States Government, Rice University, and
C   the University of Chicago.  NEITHER THE AUTHOR(S), THE UNITED
C   STATES GOVERNMENT NOR ANY AGENCY THEREOF, NOR RICE UNIVERSITY,
C   NOR THE UNIVERSITY OF CHICAGO, INCLUDING ANY OF THEIR EMPLOYEES
C   OR OFFICERS, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES
C   ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETE-
C   NESS, OR USEFULNESS OF ANY INFORMATION OR PROCESS DISCLOSED, OR
C   REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
C
C
      subroutine g_cires5(g_p_, time, n, z, g_z, ldg_z, dzdt, g_dzdt, ld
     *g_dzdt, rpar, ipar, ieform, ires)
C
C  subroutine to evaluate residual for 1995 ddasac
C     delta = zp - dz/dt
C
C  input:
C     time   - not used
C     n - number of equations (must be ns+1)
C     z      - {phi,T}
C     rpar   - not used  except rpar(n+1) = press
C     ipar   - not used   except ipar(n+1) = n+1
C       ieform - not used
C  output:
C     dzdt   - dz(t)/dt
C     ires   = 0 for successful evaluation
C
C  note:
C     ipar(n+1) must be set to n+1; and rpar(n+1) must be the pressure
C
C
C
        implicit double precision (a-h, o-z), integer (i-n)
        dimension z(n), dzdt(n), rpar(n + 1), ipar(n + 1)
C
        integer g_p_, ldg_z, ldg_dzdt
        double precision g_z(ldg_z, n), g_dzdt(ldg_dzdt, n), g_dens(g_p
     *_)
        integer g_ehfid
!Laniu        save g_dens
        external g_cidzdt
        data g_ehfid /0/
C
claniu        call ehsfid(g_ehfid, 'cires5','g_cires5.f')
C

        press = rpar(ipar(n + 1))
C
C
        call g_cidzdt(g_p_, z, g_z, ldg_z, press, dzdt, g_dzdt, ldg_dzdt
     *, dens, g_dens, g_p_)
        ires = 0
C
        return
      end subroutine g_cires5
C
C
C=========================end of subroutine g_cires5=====================================================


c=============== subroutine g_cidzdt=====================================     
C                     DISCLAIMER
C
C   This file was generated on 12/07/04 by the version of
C   ADIFOR compiled on June, 1998.
C
C   ADIFOR was prepared as an account of work sponsored by an
C   agency of the United States Government, Rice University, and
C   the University of Chicago.  NEITHER THE AUTHOR(S), THE UNITED
C   STATES GOVERNMENT NOR ANY AGENCY THEREOF, NOR RICE UNIVERSITY,
C   NOR THE UNIVERSITY OF CHICAGO, INCLUDING ANY OF THEIR EMPLOYEES
C   OR OFFICERS, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES
C   ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETE-
C   NESS, OR USEFULNESS OF ANY INFORMATION OR PROCESS DISCLOSED, OR
C   REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
C
      subroutine g_cidzdt(g_p_, z, g_z, ldg_z, press, dzdt, g_dzdt, ldg_
     *dzdt, dens, g_dens, ldg_dens)
C
C  evaluate time-rate-of-change of  z(t)
CRen in present implementation, usrate is not implemented.
C  input:
C     z    - {phi,T}
C  output:
C     dzdt - dz(t) / dt
C     dens - density
C SBP modified 2/5/98 to ensure realizability
C
C          
!         use chem_init
         use ci_dat6
        implicit double precision (a-h, o-z), integer (i-n)
!ren        parameter (ns = 1000, liwk = 6000, lrwk = 6000)
C
        dimension z(ns + 1), dzdt(ns + 1), y(ns), hh(ns)
!ren        dimension amolwt(ns), rwk(lrwk)
 !ren       integer iwk(liwk)
C
C
C  form mass fractions (y)
C

        integer g_i_, g_p_, ldg_z, ldg_dzdt, ldg_dens
        double precision d6_b, d5_b, d4_b, d3_v, d3_b, d4_v, d2_b, g_y(g
     *_p_, ns), g_z(ldg_z, ns + 1), g_t(g_p_)
        double precision g_tsum(g_p_), g_dzdt(ldg_dzdt, ns + 1), g_de
     *ns(ldg_dens), g_hh(g_p_, ns), g_cp(g_p_), g_rwk(g_p_, lrw
     *k)
        integer g_ehfid
!laniu        save g_y, g_t, g_tsum, g_hh, g_cp, g_rwk
!Laniu        save  g_t, g_tsum,  g_cp
        external g_ckhms
        external g_ckwyp
        external g_ckcpbs
        external g_ckrhoy
        data g_ehfid /0/
C
claniu        call ehsfid(g_ehfid, 'cidzdt','g_cidzdt.f')
        g_rwk=0.d0
C

        do ii = 1, ns
          do g_i_ = 1, g_p_
            g_y(g_i_, ii) = amolwt(ii) * g_z(g_i_, ii)
          enddo
          y(ii) = z(ii) * amolwt(ii)
C--------
        enddo
C  evaluate: dens, cp, dxdt, hh
C
        do g_i_ = 1, g_p_
          g_t(g_i_) = g_z(g_i_, ns + 1)
        enddo
        t = z(ns + 1)
C--------
        call g_ckrhoy(g_p_, press, t, g_t, g_p_, y, g_y, g_p_, iwk
     *, rwk, g_rwk, g_p_, dens, g_dens, ldg_dens)
        call g_ckcpbs(g_p_, t, g_t, g_p_, y, g_y, g_p_, iwk, rwk, 
     *g_rwk, g_p_, cp, g_cp, g_p_)
C
        call g_ckwyp(g_p_, press, t, g_t, g_p_, y, g_y, g_p_, iwk,
     * rwk, g_rwk, g_p_, dzdt, g_dzdt, ldg_dzdt)
C
        call g_ckhms(g_p_, t, g_t, g_p_, iwk, rwk, g_rwk, g_p_, hh
     *, g_hh, g_p_)
C
C  form rates for species and sums
C
        do g_i_ = 1, g_p_
          g_tsum(g_i_) = 0.0d0
        enddo
        tsum = 0.d0
C--------
C  ensure realizable
        do ii = 1, ns
C
Cren  if( z(ii) .le. 0. ) dzdt(ii) = max( dzdt(ii), 0.d0 )
C
          d3_v = dzdt(ii) / dens
          d2_b = 1.0d0 / dens
          d3_b = (-d3_v) / dens
          do g_i_ = 1, g_p_
            g_dzdt(g_i_, ii) = d3_b * g_dens(g_i_) + d2_b * g_dzdt(g_i_,
     * ii)
          enddo
          dzdt(ii) = d3_v
C--------
          d5_b = amolwt(ii) * dzdt(ii)
          d6_b = amolwt(ii) * hh(ii)
          do g_i_ = 1, g_p_
            g_tsum(g_i_) = d6_b * g_dzdt(g_i_, ii) + d5_b * g_hh(g_i_, i
     *i) + g_tsum(g_i_)
          enddo
          tsum = tsum + hh(ii) * dzdt(ii) * amolwt(ii)
C--------
        enddo
C      
        d4_v = (-tsum) / cp
        d3_b = (-d4_v) / cp
        d4_b = -(1.0d0 / cp)
        do g_i_ = 1, g_p_
          g_dzdt(g_i_, ns + 1) = d3_b * g_cp(g_i_) + d4_b * g_tsum(g_i_)
        enddo
        dzdt(ns + 1) = d4_v
C--------
C
        return
      end subroutine g_cidzdt

c================================end of subroutine g_cidzdt============


c======================subroutine g_cklib==================== 

C                           DISCLAIMER
C
C   This file was generated on 12/07/04 by the version of
C   ADIFOR compiled on June, 1998.
C
C   ADIFOR was prepared as an account of work sponsored by an
C   agency of the United States Government, Rice University, and
C   the University of Chicago.  NEITHER THE AUTHOR(S), THE UNITED
C   STATES GOVERNMENT NOR ANY AGENCY THEREOF, NOR RICE UNIVERSITY,
C   NOR THE UNIVERSITY OF CHICAGO, INCLUDING ANY OF THEIR EMPLOYEES
C   OR OFFICERS, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES
C   ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETE-
C   NESS, OR USEFULNESS OF ANY INFORMATION OR PROCESS DISCLOSED, OR
C   REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
C
C
C       SCCS Version  = 3.18
C       Date of delta = 08/10/94
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckabs' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckabe' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckabml' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckabms' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckaml' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckams' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckathm' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckawt' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckcdc' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckcdxp' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckcdxr' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckcdyp' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckcdyr' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckchrg' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckcomp' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckcont' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckcpbl' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      subroutine g_ckcpbs(g_p_, t, g_t, ldg_t, y, g_y, ldg_y, ickwrk, rc
     *kwrk, g_rckwrk, ldg_rckwrk, cpbms, g_cpbms, ldg_cpbms)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCPBS (T, Y, ICKWRK, RCKWRK, CPBMS)
C     Returns the mean specific heat at constant pressure;
C     see Eq. (34).
C
C  INPUT
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
C     CPBMS  - Mean specific heat at constant pressure in mass units.
C                   cgs units - ergs/(gm*K)
C                   Data type - real scalar
C
C  END PROLOGUE
C
C*****precision > double
        implicit double precision (a-h, o-z), integer (i-n)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
        dimension ickwrk(*), rckwrk(*), y(*)
C       Start of include file /home6/laniu/adifor/adiforMatab-Pre/ckstrt.h
C
C
C     include file for cklib.f
C
        common /ckstrt/ nmm, nkk, nii, mxsp, mxtb, mxtp, ncp, ncp1, ncp2
     *, ncp2t, npar, nlar, nfar, nlan, nfal, nrev, nthb, nrlt, nwl, neim
     *, njan, njar, nft1, nf1r, nexc, nrnu, nord, mxord, icmm, ickk, icn
     *c, icph, icch, icnt, icnu, icnk, icns, icnr, iclt, icrl, icrv, icw
     *l, icfl, icfo, ickf, ictb, ickn, ickt, icei, ictd, icjn, icf1, ice
     *x, icrnu, icord, ickor, ncaw, ncwt, nctt, ncaa, ncco, ncrv, nclt, 
     *ncrl, ncfl, nckt, ncwl, ncjn, ncf1, ncex, ncru, ncrc, ncpa, nckf, 
     *nckr, ncrnu, nckor, nck1, nck2, nck3, nck4, nci1, nci2, nci3, nci4
C
C     END include file for cklib.f
C
C

        integer g_i_, g_p_, ldg_cpbms, ldg_y, ldg_rckwrk, ldg_t
        double precision g_cpbms(ldg_cpbms), g_y(ldg_y, *), g_rckwrk(ldg
     *_rckwrk, *), g_t(ldg_t)
        integer g_ehfid
        external g_ckcpms
        intrinsic dble
        data g_ehfid /0/
C
claniu        call ehsfid(g_ehfid, 'ckcpbs','g_cklib.f')
C
 
        call g_ckcpms(g_p_, t, g_t, ldg_t, ickwrk, rckwrk, g_rckwrk, ldg
     *_rckwrk, rckwrk(nck1), g_rckwrk(1, nck1), ldg_rckwrk)
C
        do g_i_ = 1, g_p_
          g_cpbms(g_i_) = 0.0d0
        enddo
        cpbms = 0.0d0
C--------
        do 99956 k = 1, nkk
          do g_i_ = 1, g_p_
            g_cpbms(g_i_) = y(k) * g_rckwrk(g_i_, nck1 + k - 1) + rckwrk
     *(nck1 + k - 1) * g_y(g_i_, k) + g_cpbms(g_i_)
          enddo
          cpbms = cpbms + y(k) * rckwrk(nck1 + k - 1)
C--------
100       continue
99956   continue
C
        return
      end subroutine g_ckcpbs
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckcpml' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      subroutine g_ckcpms(g_p_, t, g_t, ldg_t, ickwrk, rckwrk, g_rckwrk,
     * ldg_rckwrk, cpms, g_cpms, ldg_cpms)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCPMS (T, ICKWRK, RCKWRK, CPMS)
C     Returns the specific heats at constant pressure in mass units;
C     see Eq. (26).
C
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT
C     CPMS   - Specific heats at constant pressure in mass units
C              for the species.
C                   cgs units - ergs/(gm*K)
C                   Data type - real array
C                   Dimension CPMS(*) at least KK, the total number of
C                   species.
C
C  END PROLOGUE
C
C*****precision > double
        implicit double precision (a-h, o-z), integer (i-n)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
        dimension ickwrk(*), rckwrk(*), cpms(*), tn(10)
C       Start of include file /home6/laniu/adifor/adiforMatab-Pre/ckstrt.h
C
C
C     include file for cklib.f
C
        common /ckstrt/ nmm, nkk, nii, mxsp, mxtb, mxtp, ncp, ncp1, ncp2
     *, ncp2t, npar, nlar, nfar, nlan, nfal, nrev, nthb, nrlt, nwl, neim
     *, njan, njar, nft1, nf1r, nexc, nrnu, nord, mxord, icmm, ickk, icn
     *c, icph, icch, icnt, icnu, icnk, icns, icnr, iclt, icrl, icrv, icw
     *l, icfl, icfo, ickf, ictb, ickn, ickt, icei, ictd, icjn, icf1, ice
     *x, icrnu, icord, ickor, ncaw, ncwt, nctt, ncaa, ncco, ncrv, nclt, 
     *ncrl, ncfl, nckt, ncwl, ncjn, ncf1, ncex, ncru, ncrc, ncpa, nckf, 
     *nckr, ncrnu, nckor, nck1, nck2, nck3, nck4, nci1, nci2, nci3, nci4
C
C     END include file for cklib.f
C
C

        integer g_i_, g_p_, i1_b, ldg_t, ldg_rckwrk, ldg_cpms
        double precision d5_b, d4_b, d3_b, d5_v, d2_v, d1_p, d2_b, g_tn(
     *g_p_, 10), g_t(ldg_t), g_sum(g_p_)
        double precision g_rckwrk(ldg_rckwrk, *), g_cpms(ldg_cpms, *)
        integer g_ehfid
!Laniu        save g_tn, g_sum
        intrinsic dble
        data g_ehfid /0/
C
claniu        call ehsfid(g_ehfid, 'ckcpms','g_cklib.f')
C

        do g_i_ = 1, g_p_
          g_tn(g_i_, 1) = 0.0d0
        enddo
        tn(1) = 1.0d0
C--------
        do 99951 n = 2, ncp

               if ( ( n - 1) .ne. 0 ) then
                  d2_v = t ** (( n - 1)-1)
                  d1_p = ( n - 1) *  d2_v
                  d2_v =  d2_v * t
               else
C            Maybe this should be  d2_v = t ** ( n - 1)
             d2_v = 1.0d0
                  d1_p = 0.0d0
               endif
     
     
     
          i1_b = 0
          do g_i_ = 1, g_p_
            g_tn(g_i_, n) = d1_p * g_t(g_i_)
          enddo
          tn(n) = d2_v
C--------
150       continue
99951   continue
C
        do 99948 k = 1, nkk
          l = 1
          do 99950 n = 2, ickwrk(icnt + k - 1) - 1
            temp = rckwrk(nctt + (k - 1) * mxtp + n - 1)
            if (t .gt. temp) then
              l = l + 1
            endif
220         continue
99950     continue
C
          na1 = ncaa + (l - 1) * ncp2 + (k - 1) * ncp2t
          do g_i_ = 1, g_p_
            g_sum(g_i_) = 0.0d0
          enddo
          sum = 0.0d0
C--------
          do 99949 n = 1, ncp
            do g_i_ = 1, g_p_
              g_sum(g_i_) = tn(n) * g_rckwrk(g_i_, na1 + n - 1) + rckwrk
     *(na1 + n - 1) * g_tn(g_i_, n) + g_sum(g_i_)
            enddo
            sum = sum + tn(n) * rckwrk(na1 + n - 1)
C--------
240         continue
99949     continue
          d5_v = rckwrk(ncru) * sum / rckwrk(ncwt + k - 1)
          d2_b = 1.0d0 / rckwrk(ncwt + k - 1)
          d3_b = (-d5_v) / rckwrk(ncwt + k - 1)
          d4_b = d2_b * sum
          d5_b = d2_b * rckwrk(ncru)
          do g_i_ = 1, g_p_
            g_cpms(g_i_, k) = d3_b * g_rckwrk(g_i_, ncwt + k - 1) + d5_b
     * * g_sum(g_i_) + d4_b * g_rckwrk(g_i_, ncru)
          enddo
          cpms(k) = d5_v
C--------
C
250       continue
99948   continue
        return
      end subroutine g_ckcpms
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckcpor' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckcray' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckctc' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckctx' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckctxp' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckctxr' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckcty' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckctyp' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckctyr' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckcvbl' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckcvbs' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckcvml' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckcvms' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckeqc' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckeqxp' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckeqxr' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckeqyp' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckeqyr' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckfal' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckgbml' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckgbms' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckgml' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckgms' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckhbml' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckhbms' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckhml' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      subroutine g_ckhms(g_p_, t, g_t, ldg_t, ickwrk, rckwrk, g_rckwrk, 
     *ldg_rckwrk, hms, g_hms, ldg_hms)
C
C  START PROLOGUE
C
C  SUBROUTINE CKHMS  (T, ICKWRK, RCKWRK, HMS)
C     Returns the enthalpies in mass units;  see Eq. (27).
C
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C  OUTPUT
C     HMS    - Enthalpies in mass units for the species.
C                   cgs units - ergs/gm
C                   Data type - real array
C                   Dimension HMS(*) at least KK, the total number of
C                   species.
C
C  END PROLOGUE
C
C*****precision > double
        implicit double precision (a-h, o-z), integer (i-n)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
        dimension ickwrk(*), rckwrk(*), hms(*), tn(10)
C       Start of include file /home6/laniu/adifor/adiforMatab-Pre/ckstrt.h
C
C
C     include file for cklib.f
C
        common /ckstrt/ nmm, nkk, nii, mxsp, mxtb, mxtp, ncp, ncp1, ncp2
     *, ncp2t, npar, nlar, nfar, nlan, nfal, nrev, nthb, nrlt, nwl, neim
     *, njan, njar, nft1, nf1r, nexc, nrnu, nord, mxord, icmm, ickk, icn
     *c, icph, icch, icnt, icnu, icnk, icns, icnr, iclt, icrl, icrv, icw
     *l, icfl, icfo, ickf, ictb, ickn, ickt, icei, ictd, icjn, icf1, ice
     *x, icrnu, icord, ickor, ncaw, ncwt, nctt, ncaa, ncco, ncrv, nclt, 
     *ncrl, ncfl, nckt, ncwl, ncjn, ncf1, ncex, ncru, ncrc, ncpa, nckf, 
     *nckr, ncrnu, nckor, nck1, nck2, nck3, nck4, nci1, nci2, nci3, nci4
C
C     END include file for cklib.f
C
C

        integer g_i_, g_p_, i1_b, ldg_t, ldg_rckwrk, ldg_hms
        double precision d9_b, d8_b, d3_b, d5_v, d9_v, d1_p, d6_v, d2_v,
     * d5_b, d4_b
        double precision d2_b, g_rut(g_p_), g_t(ldg_t), g_rckwrk(ldg_
     *rckwrk, *), g_tn(g_p_, 10), g_sum(g_p_), g_hms(ldg_hms, *)
        integer g_ehfid
!Laniu        save g_rut, g_tn, g_sum
        intrinsic dble
        data g_ehfid /0/
C
claniu        call ehsfid(g_ehfid, 'ckhms','g_cklib.f')
C

        do g_i_ = 1, g_p_
          g_rut(g_i_) = t * g_rckwrk(g_i_, ncru) + rckwrk(ncru) * g_t(g_
     *i_)
        enddo
        rut = t * rckwrk(ncru)
C--------
        do g_i_ = 1, g_p_
          g_tn(g_i_, 1) = 0.0d0
        enddo
        tn(1) = 1.0d0
C--------
        do 99910 n = 2, ncp

               if ( ( n - 1) .ne. 0 ) then
                  d2_v = t ** (( n - 1)-1)
                  d1_p = ( n - 1) *  d2_v
                  d2_v =  d2_v * t
               else
C            Maybe this should be  d2_v = t ** ( n - 1)
             d2_v = 1.0d0
                  d1_p = 0.0d0
               endif
     
     
     
          i1_b = 0
          d3_b = 1.0d0 / dble(n) * d1_p
          do g_i_ = 1, g_p_
            g_tn(g_i_, n) = d3_b * g_t(g_i_)
          enddo
          tn(n) = d2_v / dble(n)
C--------
150       continue
99910   continue
C
        do 99907 k = 1, nkk
          l = 1
          do 99909 n = 2, ickwrk(icnt + k - 1) - 1
            temp = rckwrk(nctt + (k - 1) * mxtp + n - 1)
            if (t .gt. temp) then
              l = l + 1
            endif
220         continue
99909     continue
C
          na1 = ncaa + (l - 1) * ncp2 + (k - 1) * ncp2t
          do g_i_ = 1, g_p_
            g_sum(g_i_) = 0.0d0
          enddo
          sum = 0.0d0
C--------
          do 99908 n = 1, ncp
            do g_i_ = 1, g_p_
              g_sum(g_i_) = tn(n) * g_rckwrk(g_i_, na1 + n - 1) + rckwrk
     *(na1 + n - 1) * g_tn(g_i_, n) + g_sum(g_i_)
            enddo
            sum = sum + tn(n) * rckwrk(na1 + n - 1)
C--------
225         continue
99908     continue
          d5_v = rckwrk(na1 + ncp1 - 1) / t
          d6_v = sum + d5_v
          d9_v = rut * d6_v / rckwrk(ncwt + k - 1)
          d2_b = 1.0d0 / rckwrk(ncwt + k - 1)
          d3_b = (-d9_v) / rckwrk(ncwt + k - 1)
          d4_b = d2_b * d6_v
          d5_b = d2_b * rut
          d8_b = d5_b * (1.0d0 / t)
          d9_b = d5_b * ((-d5_v) / t)
          do g_i_ = 1, g_p_
            g_hms(g_i_, k) = d3_b * g_rckwrk(g_i_, ncwt + k - 1) + d9_b 
     ** g_t(g_i_) + d8_b * g_rckwrk(g_i_, na1 + ncp1 - 1) + d5_b * g_sum
     *(g_i_) + d4_b * g_rut(g_i_)
          enddo
          hms(k) = d9_v
C--------
250       continue
99907   continue
        return
      end subroutine g_ckhms
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckhort' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C
CInactive procedure or block data 'ckhrx' excluded from output
C
C----------------------------------------------------------------------C
C
CInactive procedure or block data 'ckieim' excluded from output
C
C----------------------------------------------------------------------C
C
CInactive procedure or block data 'ckiexc' excluded from output
C
C----------------------------------------------------------------------C
C
CInactive procedure or block data 'cki2ch' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckindx' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckinit' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckinu' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckiord' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckirnu' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckitr' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckkfkr' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckkfrt' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'cklen' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckmmwc' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckmmwx' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckmmwy' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckmxtp' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckncf' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'cknpar' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'cknu' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'cknuf' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckpc' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckphaz' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckpnt' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckpx' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckpy' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckqc' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckqxp' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckqxr' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckqyp' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckqyr' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckr2ch' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckraex' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckrat' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      subroutine g_ckratt(g_p_, rckwrk, g_rckwrk, ldg_rckwrk, ickwrk, ii
     *, maxsp, ru, g_ru, ldg_ru, patm, g_patm, ldg_patm, t, g_t, ldg_t, 
     *nspec, nu, nunk, npar, par, g_par, ldg_par, nrev, irev, rpar, g_rp
     *ar, ldg_rpar, nlan, nlar, ilan, plt, g_plt, ldg_plt, nrlt, irlt, r
     *plt, g_rplt, ldg_rplt, smh, g_smh, ldg_smh, nrnu, irnu, rnu, g_rnu
     *, ldg_rnu, neim, ieim, itdep, njan, njar, ijan, pjan, g_pjan, ldg_
     *pjan, nft1, nf1r, ift1, pf1, g_pf1, ldg_pf1, rkft, g_rkft, ldg_rkf
     *t, rkrt, g_rkrt, ldg_rkrt, eqk, g_eqk, ldg_eqk)
C
C  START PROLOGUE
C
C  SUBROUTINE CKRATT (RCKWRK, ICKWRK, II, MAXSP, RU, PATM, T, NSPEC,
C 1                   NU, NUNK, NPAR, PAR, NREV, IREV, RPAR, NLAN,
C 2                   NLAR, ILAN, PLT, NRLT, IRLT, RPLT, SMH, NRNU,
C 3                   IRNU, RNU, NEIM, IEIM, ITDEP, NJAN, NJAR, IJAN,
C 4                   PJAN, NFT1, NF1R, IFT1, PF1, RKFT, RKRT, EQK)
C
C  END PROLOGUE
C
C*****precision > double
        implicit double precision (a-h, o-z), integer (i-n)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
        dimension rckwrk(*), ickwrk(*), nspec(*), nu(maxsp, *), nunk(max
     *sp, *), par(npar, *), irev(*), rpar(npar, *), ilan(*), irlt(*), pl
     *t(nlar, *), rplt(nlar, *), smh(*), rkft(*), rkrt(*), eqk(*), irnu(
     **), rnu(maxsp, *), t(*), ieim(*), itdep(*), ijan(*), pjan(njar, *)
     *, ift1(*), pf1(nf1r, *)
C
        common /mach/ small, big, exparg
C

        integer g_i_, g_p_, i1_b, ldg_t, ldg_par, ldg_rkft, ldg_plt, ldg
     *_smh, ldg_eqk, ldg_rnu
        integer ldg_patm, ldg_ru, ldg_rkrt, ldg_rpar, ldg_rplt, ldg_pjan
     *, ldg_pf1, ldg_rckwrk
        double precision d6_b, d7_b, d8_v, d2_p, d10_b, d9_v, d2_v, d8_b
     *, d2_b, d1_p
        double precision d1_w, d3_v, d4_v, d5_v, d6_v, d7_v, d3_b, d4_b,
     * d5_b, g_alogt(g_p_)
        double precision g_t(ldg_t, *), g_d1_w(g_p_), g_par(ldg_par, 
     *npar, *), g_rkft(ldg_rkft, *), g_tfac(g_p_), g_plt(ldg_plt, nla
     *r, *), g_sumsmh(g_p_), g_smh(ldg_smh, *), g_eqk(ldg_eqk, *), g_
     *rnu(ldg_rnu, maxsp, *)
        double precision g_pfac(g_p_), g_patm(ldg_patm), g_ru(ldg_ru)
     *, g_rnusum(g_p_), g_pfr(g_p_), g_rkrt(ldg_rkrt, *), g_rpar(l
     *dg_rpar, npar, *), g_rplt(ldg_rplt, nlar, *), g_temp(g_p_), g_p
     *fac2(g_p_)
        double precision g_tev(g_p_), g_sumj(g_p_), g_pjan(ldg_pja
     *n, njar, *), g_pf1(ldg_pf1, nf1r, *), g_rckwrk(ldg_rckwrk, *)
        integer g_ehfid
!Laniu        save g_sumj
!Laniu        save g_alogt, g_d1_w, g_tfac, g_sumsmh, g_pfac, g_rnusum, g_pfr,
!Laniu     * g_temp, g_pfac2, g_tev
        external g_cksmh
        data g_ehfid /0/
C
claniu        call ehsfid(g_ehfid, 'ckratt','g_cklib.f')
C

        d2_v = log(t(1))
             d1_p = 1.0d0 / t(1)
     
        do g_i_ = 1, g_p_
          g_alogt(g_i_) = d1_p * g_t(g_i_, 1)
        enddo
        alogt = d2_v
C--------
C
        do 99832 i = 1, ii
          d6_v = par(3, i) / t(1)
          d4_b = -(1.0d0 / t(1))
          d5_b = -((-d6_v) / t(1))
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d5_b * g_t(g_i_, 1) + d4_b * g_par(g_i_, 3, i
     *) + par(2, i) * g_alogt(g_i_) + alogt * g_par(g_i_, 2, i)
          enddo
          d1_w = par(2, i) * alogt - d6_v
          d3_v = exp(d1_w)
               d1_p =  d3_v
     
          d4_b = par(1, i) * d1_p
          do g_i_ = 1, g_p_
            g_rkft(g_i_, i) = d4_b * g_d1_w(g_i_) + d3_v * g_par(g_i_, 1
     *, i)
          enddo
          rkft(i) = par(1, i) * d3_v
C--------
20        continue
99832   continue
C
C        Landau-Teller reactions
C
        do 99831 n = 1, nlan
          i = ilan(n)
          d3_v = dble(1.0 / 3.0)

               if ( t(1) .ne. 0.0d0 ) then
                  d4_v = t(1) ** ( d3_v - 2.0d0)
                  d4_v =  d4_v * t(1)
                  d2_p =  d3_v *  d4_v
                  d4_v =  d4_v * t(1)
               else
C            (t(1) = 0)
             d4_v = t(1) **  d3_v
     
                  if (  d3_v .lt. 1.0d0 ) then
                     call ehbfDO (10,t(1), d3_v, d4_v, d2_p, 0.0d0,
     +g_ehfid,847)
                 else if (  d3_v .lt. 2.0d0 ) then
                     d2_p = 0.0d0
                     call ehbfDO (10,t(1), d3_v, d4_v, d2_p, 0.0d0,
     +g_ehfid,852)
                  else
                     d2_p = 0.0d0
                  endif
               endif
     
          d5_v = plt(1, n) / d4_v
          d7_v = dble(2.0 / 3.0)

               if ( t(1) .ne. 0.0d0 ) then
                  d8_v = t(1) ** ( d7_v - 2.0d0)
                  d8_v =  d8_v * t(1)
                  d1_p =  d7_v *  d8_v
                  d8_v =  d8_v * t(1)
               else
C            (t(1) = 0)
             d8_v = t(1) **  d7_v
     
                  if (  d7_v .lt. 1.0d0 ) then
                     call ehbfDO (10,t(1), d7_v, d8_v, d1_p, 0.0d0,
     +g_ehfid,873)
                  else if (  d7_v .lt. 2.0d0 ) then
                     d1_p = 0.0d0
                     call ehbfDO (10,t(1), d7_v, d8_v, d1_p, 0.0d0,
     +g_ehfid,878)
                  else
                     d1_p = 0.0d0
                  endif
               endif
     
          d9_v = plt(2, n) / d8_v
          d7_b = 0.0d0
          d10_b = 0.0d0
          d4_b = 1.0d0 / d8_v
          d8_b = 1.0d0 / d4_v
          d6_b = (-d9_v) / d8_v * d1_p + (-d5_v) / d4_v * d2_p
          do g_i_ = 1, g_p_
            g_tfac(g_i_) = d4_b * g_plt(g_i_, 2, n) + d6_b * g_t(g_i_, 1
     *) + d8_b * g_plt(g_i_, 1, n)
          enddo
          tfac = d5_v + d9_v
C--------
          d3_v = exp(tfac)
               d1_p =  d3_v
     
          d4_b = rkft(i) * d1_p
          do g_i_ = 1, g_p_
            g_rkft(g_i_, i) = d4_b * g_tfac(g_i_) + d3_v * g_rkft(g_i_, 
     *i)
          enddo
          rkft(i) = rkft(i) * d3_v
C--------
25        continue
99831   continue
C
        call g_cksmh(g_p_, t(1), g_t(1, 1), ldg_t, ickwrk, rckwrk, g_rck
     *wrk, ldg_rckwrk, smh, g_smh, ldg_smh)
        do 99829 i = 1, ii
          do g_i_ = 1, g_p_
            g_sumsmh(g_i_) = 0.0d0
          enddo
          sumsmh = 0.0d0
C--------
          do 99830 n = 1, maxsp
            if (nunk(n, i) .ne. 0) then
              d4_b = dble(nu(n, i))
              do g_i_ = 1, g_p_
                g_sumsmh(g_i_) = d4_b * g_smh(g_i_, nunk(n, i)) + g_sums
     *mh(g_i_)
              enddo
              sumsmh = sumsmh + dble(nu(n, i)) * smh(nunk(n, i))
C--------
            endif
C*****precision > double
C*****END precision > double
C*****precision > single
C     1           SUMSMH = SUMSMH + REAL(NU(N,I))*SMH(NUNK(N,I))
C*****END precision > single
40          continue
99830     continue
          if (sumsmh .ne. 0.0) then
            d2_v = min (sumsmh, exparg)
     
                 if (sumsmh .lt.  exparg) then
                    d1_p = 1.0d0
                    d2_p = 0.0d0
                 else if (sumsmh .gt.  exparg) then
                    d1_p = 0.0d0
                    d2_p = 1.0d0
                 else
                    call ehbfDO (8,sumsmh, exparg, d2_v, d1_p, d2_p,
     +g_ehfid,946)
                    d2_p = 1.0d0 -  d1_p
                 endif
     
     
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = d1_p * g_sumsmh(g_i_)
            enddo
            d1_w = d2_v
            d2_v = exp(d1_w)
                 d1_p =  d2_v
     
            do g_i_ = 1, g_p_
              g_eqk(g_i_, i) = d1_p * g_d1_w(g_i_)
            enddo
            eqk(i) = d2_v
C--------
          endif
50        continue
99829   continue
C
        do 99827 n = 1, nrnu
          do g_i_ = 1, g_p_
            g_sumsmh(g_i_) = 0.0d0
          enddo
          sumsmh = 0.0d0
C--------
          i = irnu(n)
          do 99828 l = 1, maxsp
            if (nunk(l, i) .ne. 0) then
              do g_i_ = 1, g_p_
                g_sumsmh(g_i_) = rnu(l, n) * g_smh(g_i_, nunk(l, i)) + s
     *mh(nunk(l, i)) * g_rnu(g_i_, l, n) + g_sumsmh(g_i_)
              enddo
              sumsmh = sumsmh + rnu(l, n) * smh(nunk(l, i))
C--------
            endif
45          continue
99828     continue
          if (sumsmh .ne. 0.0) then
            d2_v = min (sumsmh, exparg)
     
                 if (sumsmh .lt.  exparg) then
                    d1_p = 1.0d0
                    d2_p = 0.0d0
                 else if (sumsmh .gt.  exparg) then
                    d1_p = 0.0d0
                    d2_p = 1.0d0
                 else
                    call ehbfDO (8,sumsmh, exparg, d2_v, d1_p, d2_p,
     +g_ehfid,997)
                    d2_p = 1.0d0 -  d1_p
                 endif
    
     
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = d1_p * g_sumsmh(g_i_)
            enddo
            d1_w = d2_v
            d2_v = exp(d1_w)
                 d1_p =  d2_v
     
            do g_i_ = 1, g_p_
              g_eqk(g_i_, i) = d1_p * g_d1_w(g_i_)
            enddo
            eqk(i) = d2_v
C--------
          endif
55        continue
99827   continue
C
        d4_v = ru * t(1)
        d5_v = patm / d4_v
        d2_b = 1.0d0 / d4_v
        d3_b = (-d5_v) / d4_v
        d4_b = d3_b * t(1)
        d5_b = d3_b * ru
        do g_i_ = 1, g_p_
          g_pfac(g_i_) = d5_b * g_t(g_i_, 1) + d4_b * g_ru(g_i_) + d2_b 
     ** g_patm(g_i_)
        enddo
        pfac = d5_v
C--------
        do 99826 i = 1, ii
          nusumk = nu(1, i) + nu(2, i) + nu(3, i) + nu(4, i) + nu(5, i) 
     *+ nu(6, i)

               if (  nusumk .ne. 0 ) then
                  d3_v = pfac ** ( nusumk-1)
                  d1_p =  nusumk *  d3_v
                  d3_v =  d3_v * pfac
               else
C            Maybe this should be  d3_v = pfac **  nusumk
             d3_v = 1.0d0
                  d1_p = 0.0d0
               endif
     
     
     
          d4_b = eqk(i) * d1_p
          do g_i_ = 1, g_p_
            g_eqk(g_i_, i) = d4_b * g_pfac(g_i_) + d3_v * g_eqk(g_i_, i)
          enddo
          eqk(i) = eqk(i) * d3_v
C--------
60        continue
99826   continue
        do 99825 n = 1, nrnu
          do g_i_ = 1, g_p_
            g_rnusum(g_i_) = g_rnu(g_i_, 6, n) + g_rnu(g_i_, 5, n) + g_r
     *nu(g_i_, 4, n) + g_rnu(g_i_, 3, n) + g_rnu(g_i_, 2, n) + g_rnu(g_i
     *_, 1, n)
          enddo
          rnusum = rnu(1, n) + rnu(2, n) + rnu(3, n) + rnu(4, n) + rnu(5
     *, n) + rnu(6, n)
C--------
          i = irnu(n)

               if ( pfac .ne. 0.0d0 ) then
                  d3_v = pfac ** ( rnusum-2.0d0)
     
c            d3_v = pfac ** ( rnusum-1)
             d3_v =  d3_v * pfac
     
C            d1_p =  rnusum * (pfac ** ( rnusum-1) )
             d1_p =  rnusum *  d3_v
     
                  if ( pfac .lt. 0.0d0 ) then
c               d3_v = pfac **  rnusum
                d3_v =  d3_v * pfac
     
                     call ehbfDO (12,pfac, rnusum, d3_v, d1_p, d2_p,
     +g_ehfid,1080)
     
                  else
c               (pfac > 0 here)
C               Use  d2_p for scratch.
                d2_p = log (pfac)
     
     
C               This makes  d3_v = pfac **  rnusum
                d3_v =  d3_v * pfac
     
C               This makes  d2_p = log(pfac) * pfac **  rnusum
                d2_p =  d2_p *  d3_v
                  endif
               else
C            (pfac == 0)

                  d3_v = pfac **  rnusum
     
                  d1_p = 0.0d0
                  call ehbfDO (12,pfac, rnusum, d3_v, d1_p, d2_p,
     +g_ehfid,1102)
     
               endif
     
          do g_i_ = 1, g_p_
            g_pfr(g_i_) = d2_p * g_rnusum(g_i_) + d1_p * g_pfac(g_i_)
          enddo
          pfr = d3_v
C--------
          do g_i_ = 1, g_p_
            g_eqk(g_i_, i) = eqk(i) * g_pfr(g_i_) + pfr * g_eqk(g_i_, i)
          enddo
          eqk(i) = eqk(i) * pfr
C--------
65        continue
99825   continue
C
        do 99824 i = 1, ii
C
C     RKR=0.0 for irreversible reactions, else RKR=RKF/MAX(EQK,SMALL)
C
          do g_i_ = 1, g_p_
            g_rkrt(g_i_, i) = 0.0d0
          enddo
          rkrt(i) = 0.0d0
C--------
          if (nspec(i) .gt. 0) then
            d3_v = max (eqk(i), small)
     
                 if (eqk(i) .gt.  small) then
                    d1_p = 1.0d0
                    d2_p = 0.0d0
                 else if (eqk(i) .lt.  small) then
                    d1_p = 0.0d0
                    d2_p = 1.0d0
                 else
                    call ehbfDO (7,eqk(i), small, d3_v, d1_p, d2_p,
     +g_ehfid,1140)
                    d2_p = 1.0d0 -  d1_p
                 endif
     
     
            d4_v = rkft(i) / d3_v
            d2_b = 1.0d0 / d3_v
            d4_b = (-d4_v) / d3_v * d1_p
            do g_i_ = 1, g_p_
              g_rkrt(g_i_, i) = d4_b * g_eqk(g_i_, i) + d2_b * g_rkft(g_
     *i_, i)
            enddo
            rkrt(i) = d4_v
C--------
          endif
68        continue
99824   continue
C
C     if reverse parameters have been given:
C
        do 99823 n = 1, nrev
          i = irev(n)
          d6_v = rpar(3, n) / t(1)
          d4_b = -(1.0d0 / t(1))
          d5_b = -((-d6_v) / t(1))
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d5_b * g_t(g_i_, 1) + d4_b * g_rpar(g_i_, 3, 
     *n) + rpar(2, n) * g_alogt(g_i_) + alogt * g_rpar(g_i_, 2, n)
          enddo
          d1_w = rpar(2, n) * alogt - d6_v
          d3_v = exp(d1_w)
               d1_p =  d3_v
     
          d4_b = rpar(1, n) * d1_p
          do g_i_ = 1, g_p_
            g_rkrt(g_i_, i) = d4_b * g_d1_w(g_i_) + d3_v * g_rpar(g_i_, 
     *1, n)
          enddo
          rkrt(i) = rpar(1, n) * d3_v
C--------
          d3_v = rkft(i) / rkrt(i)
          d2_b = 1.0d0 / rkrt(i)
          d3_b = (-d3_v) / rkrt(i)
          do g_i_ = 1, g_p_
            g_eqk(g_i_, i) = d3_b * g_rkrt(g_i_, i) + d2_b * g_rkft(g_i_
     *, i)
          enddo
          eqk(i) = d3_v
C--------
70        continue
99823   continue
C
C     if reverse Landau-Teller parameters have been given:
C
        do 99822 n = 1, nrlt
          i = irlt(n)
          d3_v = dble(1.0 / 3.0)

               if ( t(1) .ne. 0.0d0 ) then
                  d4_v = t(1) ** ( d3_v - 2.0d0)
                  d4_v =  d4_v * t(1)
                  d2_p =  d3_v *  d4_v
                  d4_v =  d4_v * t(1)
               else
C            (t(1) = 0)
             d4_v = t(1) **  d3_v
     
                  if (  d3_v .lt. 1.0d0 ) then
                     call ehbfDO (10,t(1), d3_v, d4_v, d2_p, 0.0d0,
     +g_ehfid,1210)
                  else if (  d3_v .lt. 2.0d0 ) then
                     d2_p = 0.0d0
                     call ehbfDO (10,t(1), d3_v, d4_v, d2_p, 0.0d0,
     +g_ehfid,1215)
                  else
                     d2_p = 0.0d0
                  endif
               endif
     
          d5_v = rplt(1, n) / d4_v
          d7_v = dble(2.0 / 3.0)

               if ( t(1) .ne. 0.0d0 ) then
                  d8_v = t(1) ** ( d7_v - 2.0d0)
                  d8_v =  d8_v * t(1)
                  d1_p =  d7_v *  d8_v
                  d8_v =  d8_v * t(1)
               else
C            (t(1) = 0)
             d8_v = t(1) **  d7_v
     
                  if (  d7_v .lt. 1.0d0 ) then
                     call ehbfDO (10,t(1), d7_v, d8_v, d1_p, 0.0d0,
     +g_ehfid,1236)
                  else if (  d7_v .lt. 2.0d0 ) then
                     d1_p = 0.0d0
                     call ehbfDO (10,t(1), d7_v, d8_v, d1_p, 0.0d0,
     +g_ehfid,1241)
                  else
                     d1_p = 0.0d0
                  endif
               endif
     
          d9_v = rplt(2, n) / d8_v
          d7_b = 0.0d0
          d10_b = 0.0d0
          d4_b = 1.0d0 / d8_v
          d8_b = 1.0d0 / d4_v
          d6_b = (-d9_v) / d8_v * d1_p + (-d5_v) / d4_v * d2_p
          do g_i_ = 1, g_p_
            g_tfac(g_i_) = d4_b * g_rplt(g_i_, 2, n) + d6_b * g_t(g_i_, 
     *1) + d8_b * g_rplt(g_i_, 1, n)
          enddo
          tfac = d5_v + d9_v
C--------
          d3_v = exp(tfac)
               d1_p =  d3_v
     
          d4_b = rkrt(i) * d1_p
          do g_i_ = 1, g_p_
            g_rkrt(g_i_, i) = d4_b * g_tfac(g_i_) + d3_v * g_rkrt(g_i_, 
     *i)
          enddo
          rkrt(i) = rkrt(i) * d3_v
C--------
          d3_v = rkft(i) / rkrt(i)
          d2_b = 1.0d0 / rkrt(i)
          d3_b = (-d3_v) / rkrt(i)
          do g_i_ = 1, g_p_
            g_eqk(g_i_, i) = d3_b * g_rkrt(g_i_, i) + d2_b * g_rkft(g_i_
     *, i)
          enddo
          eqk(i) = d3_v
C--------
75        continue
99822   continue
C
C     electron-impact reactions
C
        do 99819 n = 1, neim
          i = ieim(n)
          do g_i_ = 1, g_p_
            g_rkrt(g_i_, i) = 0.0d0
          enddo
          rkrt(i) = 0.0d0
C--------
          do g_i_ = 1, g_p_
            g_temp(g_i_) = g_t(g_i_, itdep(n))
          enddo
          temp = t(itdep(n))
C--------
          d3_v = log(temp)
               d1_p = 1.0d0 / temp
     
          d6_v = par(3, i) / temp
          d4_b = -(1.0d0 / temp)
          d5_b = -((-d6_v) / temp) + par(2, i) * d1_p
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d4_b * g_par(g_i_, 3, i) + d5_b * g_temp(g_i_
     *) + d3_v * g_par(g_i_, 2, i)
          enddo
          d1_w = par(2, i) * d3_v - d6_v
          d3_v = exp(d1_w)
               d1_p =  d3_v
     
          d4_b = par(1, i) * d1_p
          do g_i_ = 1, g_p_
            g_rkft(g_i_, i) = d4_b * g_d1_w(g_i_) + d3_v * g_par(g_i_, 1
     *, i)
          enddo
          rkft(i) = par(1, i) * d3_v
C--------
          d4_v = ru * temp
          d5_v = patm / d4_v
          d2_b = 1.0d0 / d4_v
          d3_b = (-d5_v) / d4_v
          d4_b = d3_b * temp
          d5_b = d3_b * ru
          do g_i_ = 1, g_p_
            g_pfac2(g_i_) = d5_b * g_temp(g_i_) + d4_b * g_ru(g_i_) + d2
     *_b * g_patm(g_i_)
          enddo
          pfac2 = d5_v
C--------
          nusumk = nu(1, i) + nu(2, i) + nu(3, i) + nu(4, i) + nu(5, i) 
     *+ nu(6, i)
          d4_v = pfac2 / pfac

               if (  nusumk .ne. 0 ) then
                  d5_v = d4_v ** ( nusumk-1)
                  d1_p =  nusumk *  d5_v
                  d5_v =  d5_v * d4_v
               else
C            Maybe this should be  d5_v = d4_v **  nusumk
             d5_v = 1.0d0
                  d1_p = 0.0d0
               endif
     
     
     
          d4_b = eqk(i) * d1_p
          d5_b = d4_b * (1.0d0 / pfac)
          d6_b = d4_b * ((-d4_v) / pfac)
          do g_i_ = 1, g_p_
            g_eqk(g_i_, i) = d6_b * g_pfac(g_i_) + d5_b * g_pfac2(g_i_) 
     *+ d5_v * g_eqk(g_i_, i)
          enddo
          eqk(i) = eqk(i) * d5_v
C--------
C
          do 99821 l = 1, nrnu
            if (irnu(l) .eq. i) then
              do g_i_ = 1, g_p_
                g_rnusum(g_i_) = g_rnu(g_i_, 6, l) + g_rnu(g_i_, 5, l) +
     * g_rnu(g_i_, 4, l) + g_rnu(g_i_, 3, l) + g_rnu(g_i_, 2, l) + g_rnu
     *(g_i_, 1, l)
              enddo
              rnusum = rnu(1, l) + rnu(2, l) + rnu(3, l) + rnu(4, l) + r
     *nu(5, l) + rnu(6, l)
C--------
              d3_v = pfac2 / pfac

                   if ( d3_v .ne. 0.0d0 ) then
                      d5_v = d3_v ** ( rnusum-2.0d0)
     
c                d5_v = d3_v ** ( rnusum-1)
                 d5_v =  d5_v * d3_v
     
C                d1_p =  rnusum * (d3_v ** ( rnusum-1) )
                 d1_p =  rnusum *  d5_v
     
                      if ( d3_v .lt. 0.0d0 ) then
c                   d5_v = d3_v **  rnusum
                    d5_v =  d5_v * d3_v
     
                         call ehbfDO (12,d3_v, rnusum, d5_v, d1_p, d2_p,
     +g_ehfid,1381)
     
                      else
c                   (d3_v > 0 here)
C                   Use  d2_p for scratch.
                    d2_p = log (d3_v)
     
     
C                   This makes  d5_v = d3_v **  rnusum
                    d5_v =  d5_v * d3_v
     
C                   This makes  d2_p = log(d3_v) * d3_v **  rnusum
                    d2_p =  d2_p *  d5_v
                      endif
                   else
C                (d3_v == 0)

                      d5_v = d3_v **  rnusum
     
                      d1_p = 0.0d0
                      call ehbfDO (12,d3_v, rnusum, d5_v, d1_p, d2_p,
     +g_ehfid,1403)
     
                   endif
     
              d4_b = d1_p * (1.0d0 / pfac)
              d5_b = d1_p * ((-d3_v) / pfac)
              do g_i_ = 1, g_p_
                g_pfr(g_i_) = d2_p * g_rnusum(g_i_) + d5_b * g_pfac(g_i_
     *) + d4_b * g_pfac2(g_i_)
              enddo
              pfr = d5_v
C--------
              do g_i_ = 1, g_p_
                g_eqk(g_i_, i) = eqk(i) * g_pfr(g_i_) + pfr * g_eqk(g_i_
     *, i)
              enddo
              eqk(i) = eqk(i) * pfr
C--------
            endif
80          continue
99821     continue
C
          if (nspec(i) .gt. 0) then
            d3_v = max (eqk(i), small)
     
                 if (eqk(i) .gt.  small) then
                    d1_p = 1.0d0
                    d2_p = 0.0d0
                 else if (eqk(i) .lt.  small) then
                    d1_p = 0.0d0
                    d2_p = 1.0d0
                 else
                    call ehbfDO (7,eqk(i), small, d3_v, d1_p, d2_p,
     +g_ehfid,1437)
                    d2_p = 1.0d0 -  d1_p
                 endif
     
     
            d4_v = rkft(i) / d3_v
            d2_b = 1.0d0 / d3_v
            d4_b = (-d4_v) / d3_v * d1_p
            do g_i_ = 1, g_p_
              g_rkrt(g_i_, i) = d4_b * g_eqk(g_i_, i) + d2_b * g_rkft(g_
     *i_, i)
            enddo
            rkrt(i) = d4_v
C--------
          endif
          do 99820 n2 = 1, nrev
            i2 = irev(n2)
            if (i2 .eq. i) then
              d3_v = log(temp)
                   d1_p = 1.0d0 / temp
     
              d6_v = rpar(3, n2) / temp
              d4_b = -(1.0d0 / temp)
              d5_b = -((-d6_v) / temp) + rpar(2, n2) * d1_p
              do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d4_b * g_rpar(g_i_, 3, n2) + d5_b * g_tem
     *p(g_i_) + d3_v * g_rpar(g_i_, 2, n2)
              enddo
              d1_w = rpar(2, n2) * d3_v - d6_v
              d3_v = exp(d1_w)
                   d1_p =  d3_v
     
              d4_b = rpar(1, n2) * d1_p
              do g_i_ = 1, g_p_
                g_rkrt(g_i_, i) = d4_b * g_d1_w(g_i_) + d3_v * g_rpar(g_
     *i_, 1, n2)
              enddo
              rkrt(i) = rpar(1, n2) * d3_v
C--------
              d3_v = max (rkrt(i), small)
     
                   if (rkrt(i) .gt.  small) then
                      d1_p = 1.0d0
                      d2_p = 0.0d0
                   else if (rkrt(i) .lt.  small) then
                      d1_p = 0.0d0
                      d2_p = 1.0d0
                   else
                      call ehbfDO (7,rkrt(i), small, d3_v, d1_p, d2_p,
     +g_ehfid,1487)
                      d2_p = 1.0d0 -  d1_p
                   endif
     
     
              d4_v = rkft(i) / d3_v
              d2_b = 1.0d0 / d3_v
              d4_b = (-d4_v) / d3_v * d1_p
              do g_i_ = 1, g_p_
                g_eqk(g_i_, i) = d4_b * g_rkrt(g_i_, i) + d2_b * g_rkft(
     *g_i_, i)
              enddo
              eqk(i) = d4_v
C--------
            endif
83          continue
99820     continue
85        continue
99819   continue
C
C      jannev, langer, evans & post - type reactions
C
        do 99817 n = 1, njan
          i = ijan(n)
          do g_i_ = 1, g_p_
            g_rkrt(g_i_, i) = 0.0d0
          enddo
          rkrt(i) = 0.0d0
C--------
C
C       CONVERT E- TEMPERATURE TO eV's
C
          do g_i_ = 1, g_p_
            g_temp(g_i_) = g_t(g_i_, 2)
          enddo
          temp = t(2)
C--------
          d2_b = 1.0d0 / dble(11600.)
          do g_i_ = 1, g_p_
            g_tev(g_i_) = d2_b * g_temp(g_i_)
          enddo
          tev = temp / dble(11600.)
C--------
          d3_v = log(temp)
               d1_p = 1.0d0 / temp
     
          d6_v = par(3, i) / temp
          d4_b = -(1.0d0 / temp)
          d5_b = -((-d6_v) / temp) + par(2, i) * d1_p
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d4_b * g_par(g_i_, 3, i) + d5_b * g_temp(g_i_
     *) + d3_v * g_par(g_i_, 2, i)
          enddo
          d1_w = par(2, i) * d3_v - d6_v
          d3_v = exp(d1_w)
               d1_p =  d3_v
     
          d4_b = par(1, i) * d1_p
          do g_i_ = 1, g_p_
            g_rkft(g_i_, i) = d4_b * g_d1_w(g_i_) + d3_v * g_par(g_i_, 1
     *, i)
          enddo
          rkft(i) = par(1, i) * d3_v
C--------
          do g_i_ = 1, g_p_
            g_sumj(g_i_) = 0.0d0
          enddo
          sumj = 0.0d0
C--------
          do 99818 j = 1, njar
            d4_v = log(tev)
                 d2_p = 1.0d0 / tev
     

                 if ( ( j - 1) .ne. 0 ) then
                    d5_v = d4_v ** (( j - 1)-1)
                    d1_p = ( j - 1) *  d5_v
                    d5_v =  d5_v * d4_v
                 else
C              Maybe this should be  d5_v = d4_v ** ( j - 1)
               d5_v = 1.0d0
                    d1_p = 0.0d0
                 endif
     
     
     
            i1_b = 0
            d7_b = pjan(j, n) * d1_p * d2_p
            do g_i_ = 1, g_p_
              g_sumj(g_i_) = d7_b * g_tev(g_i_) + d5_v * g_pjan(g_i_, j,
     * n) + g_sumj(g_i_)
            enddo
            sumj = sumj + pjan(j, n) * d5_v
C--------
90          continue
99818     continue
          d3_v = exp(sumj)
               d1_p =  d3_v
     
          d4_b = rkft(i) * d1_p
          do g_i_ = 1, g_p_
            g_rkft(g_i_, i) = d4_b * g_sumj(g_i_) + d3_v * g_rkft(g_i_, 
     *i)
          enddo
          rkft(i) = rkft(i) * d3_v
C--------
95        continue
99817   continue
C
C      reactions using fit#1:  k = A * T^B * exp(v1/T+v2/T^2+v3/T^3...)
C
        do 99814 n = 1, nft1
          i = ift1(n)
          do g_i_ = 1, g_p_
            g_rkrt(g_i_, i) = 0.0d0
          enddo
          rkrt(i) = 0.0d0
C--------
C
C         CHECK IF REACTION IS ALSO AN ELECTRON-IMPACT REAX
C
          do g_i_ = 1, g_p_
            g_temp(g_i_) = g_t(g_i_, 1)
          enddo
          temp = t(1)
C--------
          do 99816 n2 = 1, neim
            if (ieim(n2) .eq. i) then
              do g_i_ = 1, g_p_
                g_temp(g_i_) = g_t(g_i_, itdep(n2))
              enddo
              temp = t(itdep(n2))
C--------
            endif
100         continue
99816     continue
C

               if ( temp .ne. 0.0d0 ) then
                  d4_v = temp ** ( par(2, i)-2.0d0)
     
c            d4_v = temp ** ( par(2, i)-1)
             d4_v =  d4_v * temp
     
C            d1_p =  par(2, i) * (temp ** ( par(2, i)-1) )
             d1_p =  par(2, i) *  d4_v
     
                  if ( temp .lt. 0.0d0 ) then
c               d4_v = temp **  par(2, i)
                d4_v =  d4_v * temp
     
                     call ehbfDO (12,temp, par(2, i), d4_v, d1_p, d2_p,
     +g_ehfid,1640)
     
                  else
c               (temp > 0 here)
C               Use  d2_p for scratch.
                d2_p = log (temp)
     
     
C               This makes  d4_v = temp **  par(2, i)
                d4_v =  d4_v * temp
     
C               This makes  d2_p = log(temp) * temp **  par(2, i)
                d2_p =  d2_p *  d4_v
                  endif
               else
C            (temp == 0)

                  d4_v = temp **  par(2, i)
     
                  d1_p = 0.0d0
                  call ehbfDO (12,temp, par(2, i), d4_v, d1_p, d2_p,
     +g_ehfid,1662)
     
               endif
     
          d4_b = par(1, i) * d1_p
          d5_b = par(1, i) * d2_p
          do g_i_ = 1, g_p_
            g_rkft(g_i_, i) = d5_b * g_par(g_i_, 2, i) + d4_b * g_temp(g
     *_i_) + d4_v * g_par(g_i_, 1, i)
          enddo
          rkft(i) = par(1, i) * d4_v
C--------
          do g_i_ = 1, g_p_
            g_sumj(g_i_) = 0.0d0
          enddo
          sumj = 0.0d0
C--------
          do 99815 j = 1, nf1r
            rootj = 1. / float(j)
            if (temp .ge. big ** rootj .or. sumj .ge. log(big)) then
              do g_i_ = 1, g_p_
                g_sumj(g_i_) = 0.0d0
              enddo
              sumj = log(big)
C--------
            else
              d4_v = dble(float(j))

                   if ( temp .ne. 0.0d0 ) then
                      d5_v = temp ** ( d4_v - 2.0d0)
                      d5_v =  d5_v * temp
                      d1_p =  d4_v *  d5_v
                      d5_v =  d5_v * temp
                   else
C                (temp = 0)
                 d5_v = temp **  d4_v
     
                      if (  d4_v .lt. 1.0d0 ) then
                         call ehbfDO (10,temp, d4_v, d5_v, d1_p, 0.0d0,
     +g_ehfid,1702)
                      else if (  d4_v .lt. 2.0d0 ) then
                         d1_p = 0.0d0
                         call ehbfDO (10,temp, d4_v, d5_v, d1_p, 0.0d0,
     +g_ehfid,1707)
                      else
                         d1_p = 0.0d0
                      endif
                   endif
     
              d6_v = pf1(j, n) / d5_v
              d7_b = 0.0d0
              d4_b = 1.0d0 / d5_v
              d6_b = (-d6_v) / d5_v * d1_p
              do g_i_ = 1, g_p_
                g_sumj(g_i_) = d6_b * g_temp(g_i_) + d4_b * g_pf1(g_i_, 
     *j, n) + g_sumj(g_i_)
              enddo
              sumj = sumj + d6_v
C--------
            endif
103         continue
99815     continue
          if (sumj .ge. log(big)) then
            do g_i_ = 1, g_p_
              g_sumj(g_i_) = 0.0d0
            enddo
            sumj = log(big)
C--------
          endif
          d3_v = exp(sumj)
               d1_p =  d3_v
     
          d4_b = rkft(i) * d1_p
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d4_b * g_sumj(g_i_) + d3_v * g_rkft(g_i_, i)
          enddo
          d1_w = rkft(i) * d3_v
          d2_v = min (big, d1_w)
     
               if (big .lt.  d1_w) then
                  d1_p = 1.0d0
                  d2_p = 0.0d0
               else if (big .gt.  d1_w) then
                  d1_p = 0.0d0
                  d2_p = 1.0d0
               else
                  call ehbfDO (8,big, d1_w, d2_v, d1_p, d2_p,
     +g_ehfid,1752)
                  d2_p = 1.0d0 -  d1_p
               endif
     
     
          do g_i_ = 1, g_p_
            g_rkft(g_i_, i) = d2_p * g_d1_w(g_i_)
          enddo
          rkft(i) = d2_v
C--------
105       continue
99814   continue
C
        return
      end subroutine g_ckratt
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      subroutine g_ckratx(g_p_, ii, kk, maxsp, maxtb, t, g_t, ldg_t, c, 
     *g_c, ldg_c, nu, nunk, npar, par, g_par, ldg_par, nfal, ifal, ifop,
     * kfal, nfar, fpar, g_fpar, ldg_fpar, nthb, ithb, ntbs, aik, g_aik,
     * ldg_aik, nktb, rkft, g_rkft, ldg_rkft, rkrt, g_rkrt, ldg_rkrt, rk
     *f, g_rkf, ldg_rkf, rkr, g_rkr, ldg_rkr, ctb, g_ctb, ldg_ctb, nrnu,
     * irnu, rnu, g_rnu, ldg_rnu, nord, iord, mxord, kord, rord, g_rord,
     * ldg_rord)
C
C  START PROLOGUE
C
C  SUBROUTINE CKRATX (II, KK, MAXSP, MAXTB, T, C, NU, NUNK, NPAR,
C 1                   PAR, NFAL, IFAL, IFOP, KFAL, NFAR, FPAR, NTHB,
C 2                   ITHB, NTBS, AIK, NKTB, RKFT, RKRT, RKF, RKR,
C 3                   CTB, NRNU, IRNU, RNU, NORD, IORD, MXORD, KORD,
C 4                   RORD)
C
C  END PROLOGUE
C
C*****precision > double
        implicit double precision (a-h, o-z), integer (i-n)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
        dimension c(*), nu(maxsp, *), nunk(maxsp, *), par(npar, *), ifal
     *(*), ifop(*), kfal(*), fpar(nfar, *), ithb(*), ntbs(*), aik(maxtb,
     * *), nktb(maxtb, *), rkft(*), rkrt(*), rkf(*), rkr(*), ctb(*), irn
     *u(*), rnu(maxsp, *), iord(*), kord(mxord, *), rord(mxord, *)
C
        common /mach/ small, big, exparg
C

        integer g_i_, g_p_, i1_b, ldg_ctb, ldg_rkf, ldg_rkr, ldg_c, ldg_
     *aik, ldg_t, ldg_fpar
        integer ldg_rkft, ldg_rkrt, ldg_par, ldg_rnu, ldg_rord
        double precision d8_b, d6_p, d5_p, d4_p, d3_p, d15_b, d14_v, d2_
     *v, d3_v, d13_b
        double precision d2_b, d3_b, d4_v, d5_v, d6_v, d4_b, d5_b, d6_b,
     * d1_p, d1_w
        double precision d7_v, d7_b, d2_p, d2_w, d8_v, d9_v, d12_b, d11_
     *v, d9_b, g_ctb(ldg_ctb, *)
        double precision g_rkf(ldg_rkf, *), g_rkr(ldg_rkr, *), g_ctot(g_
     *p_), g_c(ldg_c, *), g_aik(ldg_aik, maxtb, *), g_alogt(g_p_),
     * g_t(ldg_t), g_d1_w(g_p_), g_fpar(ldg_fpar, nfar, *), g_rklow(g
     *_p_)
        double precision g_pr(g_p_), g_rkft(ldg_rkft, *), g_pcor(g_p
     *_), g_prlog(g_p_), g_xp(g_p_), g_d2_w(g_p_), g_fc(g_p
     *_), g_fcent(g_p_), g_fclog(g_p_), g_xn(g_p_)
        double precision g_cprlog(g_p_), g_flog(g_p_), g_rkrt(ldg_
     *rkrt, *), g_par(ldg_par, npar, *), g_c1(g_p_), g_rnu(ldg_rnu, m
     *axsp, *), g_c4(g_p_), g_c2(g_p_), g_c3(g_p_), g_c5(g_p
     *_)
        double precision g_c6(g_p_), g_cnk(g_p_), g_rord(ldg_rord,
     * mxord, *)
        integer g_ehfid
!Laniu        save g_c6, g_cnk
!Laniu        save g_fcent, g_fclog, g_xn, g_cprlog, g_flog, g_c1, g_c4, g_c2,
!Laniu     * g_c3, g_c5
!Laniu        save g_ctot, g_alogt, g_d1_w, g_rklow, g_pr, g_pcor, g_prlog, g_
!Laniu     *xp, g_d2_w, g_fc
        intrinsic dble
        data g_ehfid /0/
C
claniu        call ehsfid(g_ehfid, 'ckratx','g_cklib.f')
C

        do 99813 i = 1, ii
          do g_i_ = 1, g_p_
            g_ctb(g_i_, i) = 0.0d0
          enddo
          ctb(i) = 1.0d0
C--------
          do g_i_ = 1, g_p_
            g_rkf(g_i_, i) = 0.0d0
          enddo
          rkf(i) = 0.0d0
C--------
          do g_i_ = 1, g_p_
            g_rkr(g_i_, i) = 0.0d0
          enddo
          rkr(i) = 0.0d0
C--------
20        continue
99813   continue
C
C     third-body reactions
C
        if (nthb .gt. 0) then
          do g_i_ = 1, g_p_
            g_ctot(g_i_) = 0.0d0
          enddo
          ctot = 0.0d0
C--------
          do 99812 k = 1, kk
            do g_i_ = 1, g_p_
              g_ctot(g_i_) = g_c(g_i_, k) + g_ctot(g_i_)
            enddo
            ctot = ctot + c(k)
C--------
10          continue
99812     continue
          do 99810 n = 1, nthb
            do g_i_ = 1, g_p_
              g_ctb(g_i_, ithb(n)) = g_ctot(g_i_)
            enddo
            ctb(ithb(n)) = ctot
C--------
            do 99811 l = 1, ntbs(n)
              d3_v = aik(l, n) - 1.0d0
              do g_i_ = 1, g_p_
                g_ctb(g_i_, ithb(n)) = d3_v * g_c(g_i_, nktb(l, n)) + c(
     *nktb(l, n)) * g_aik(g_i_, l, n) + g_ctb(g_i_, ithb(n))
              enddo
              ctb(ithb(n)) = ctb(ithb(n)) + d3_v * c(nktb(l, n))
C--------
80            continue
99811       continue
99810     continue
        endif
C
C     If fall-off (pressure correction):
C
        if (nfal .gt. 0) then
          d2_v = log(t)
               d1_p = 1.0d0 / t
     
          do g_i_ = 1, g_p_
            g_alogt(g_i_) = d1_p * g_t(g_i_)
          enddo
          alogt = d2_v
C--------
C
          do 99809 n = 1, nfal
C
            d6_v = fpar(3, n) / t
            d4_b = -(1.0d0 / t)
            d5_b = -((-d6_v) / t)
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = d5_b * g_t(g_i_) + d4_b * g_fpar(g_i_, 3, n
     *) + fpar(2, n) * g_alogt(g_i_) + alogt * g_fpar(g_i_, 2, n)
            enddo
            d1_w = fpar(2, n) * alogt - d6_v
            d3_v = exp(d1_w)
                 d1_p =  d3_v
     
            d4_b = fpar(1, n) * d1_p
            do g_i_ = 1, g_p_
              g_rklow(g_i_) = d4_b * g_d1_w(g_i_) + d3_v * g_fpar(g_i_, 
     *1, n)
            enddo
            rklow = fpar(1, n) * d3_v
C--------
C
C        CONCENTRATION OF THIRD BODY
C
            if (kfal(n) .eq. 0) then
              d5_v = rklow * ctb(ifal(n)) / rkft(ifal(n))
              d2_b = 1.0d0 / rkft(ifal(n))
              d3_b = (-d5_v) / rkft(ifal(n))
              d4_b = d2_b * ctb(ifal(n))
              d5_b = d2_b * rklow
              do g_i_ = 1, g_p_
                g_pr(g_i_) = d3_b * g_rkft(g_i_, ifal(n)) + d5_b * g_ctb
     *(g_i_, ifal(n)) + d4_b * g_rklow(g_i_)
              enddo
              pr = d5_v
C--------
              do g_i_ = 1, g_p_
                g_ctb(g_i_, ifal(n)) = 0.0d0
              enddo
              ctb(ifal(n)) = 1.0d0
C--------
            else
              d5_v = rklow * c(kfal(n)) / rkft(ifal(n))
              d2_b = 1.0d0 / rkft(ifal(n))
              d3_b = (-d5_v) / rkft(ifal(n))
              d4_b = d2_b * c(kfal(n))
              d5_b = d2_b * rklow
              do g_i_ = 1, g_p_
                g_pr(g_i_) = d3_b * g_rkft(g_i_, ifal(n)) + d5_b * g_c(g
     *_i_, kfal(n)) + d4_b * g_rklow(g_i_)
              enddo
              pr = d5_v
C--------
            endif
C
            d2_v = 1.0d0 + pr
            d3_v = pr / d2_v
            d2_b = 1.0d0 / d2_v + (-d3_v) / d2_v
            do g_i_ = 1, g_p_
              g_pcor(g_i_) = d2_b * g_pr(g_i_)
            enddo
            pcor = d3_v
C--------
C
            if (ifop(n) .gt. 1) then
              d2_v = max (pr, small)
     
                   if (pr .gt.  small) then
                      d1_p = 1.0d0
                      d2_p = 0.0d0
                   else if (pr .lt.  small) then
                      d1_p = 0.0d0
                      d2_p = 1.0d0
                   else
                      call ehbfDO (7,pr, small, d2_v, d1_p, d2_p,
     +g_ehfid,1982)
                      d2_p = 1.0d0 -  d1_p
                   endif
     
     
              do g_i_ = 1, g_p_
                g_d1_w(g_i_) = d1_p * g_pr(g_i_)
              enddo
              d1_w = d2_v
              d2_v = log10 (d1_w)
                   d1_p = 1.0d0 / ( d1_w * log (10.0d0) )
     
              do g_i_ = 1, g_p_
                g_prlog(g_i_) = d1_p * g_d1_w(g_i_)
              enddo
              prlog = d2_v
C--------
C
              if (ifop(n) .eq. 2) then
C
C              8-PARAMETER SRI FORM
C
                d2_v = prlog * prlog
                     d1_p = 2.0d0 * prlog
     
                d3_v = 1.0d0 + d2_v
                d4_v = 1.0d0 / d3_v
                d4_b = (-d4_v) / d3_v * d1_p
                do g_i_ = 1, g_p_
                  g_xp(g_i_) = d4_b * g_prlog(g_i_)
                enddo
                xp = d4_v
C--------
                d4_v = (-fpar(5, n)) / t
                d3_b = (-d4_v) / t
                d4_b = -(1.0d0 / t)
                do g_i_ = 1, g_p_
                  g_d1_w(g_i_) = d3_b * g_t(g_i_) + d4_b * g_fpar(g_i_, 
     *5, n)
                enddo
                d1_w = d4_v
                d4_v = (-t) / fpar(6, n)
                d3_b = (-d4_v) / fpar(6, n)
                d4_b = -(1.0d0 / fpar(6, n))
                do g_i_ = 1, g_p_
                  g_d2_w(g_i_) = d3_b * g_fpar(g_i_, 6, n) + d4_b * g_t(
     *g_i_)
                enddo
                d2_w = d4_v
                d3_v = exp(d1_w)
                     d6_p =  d3_v
     
               d6_v = exp(d2_w)
                     d5_p =  d6_v
     

                     if ( (fpar(4, n) * d3_v + d6_v) .ne. 0.0d0 ) then
                        d9_v = (fpar(4, n) * d3_v + d6_v) ** ( xp-2.0d0)
     
c                  d9_v = (fpar(4, n) * d3_v + d6_v) ** ( xp-1)
                   d9_v =  d9_v * (fpar(4, n) * d3_v + d6_v)
     
C                  d3_p =  xp * ((fpar(4, n) * d3_v + d6_v) ** ( xp-1) )
                   d3_p =  xp *  d9_v
     
                      if ( (fpar(4, n) * d3_v + d6_v) .lt. 0.0d0 ) then
c                     d9_v = (fpar(4, n) * d3_v + d6_v) **  xp
                      d9_v =  d9_v * (fpar(4, n) * d3_v + d6_v)
     
                     call ehbfDO (12,fpar(4, n) * d3_v + d6_v, xp, d9_v
     +, d3_p, d4_p,g_ehfid,2054)
     
                        else
c                     ((fpar(4, n) * d3_v + d6_v) > 0 here)
C                     Use  d4_p for scratch.
                      d4_p = log ((fpar(4, n) * d3_v + d6_v))
     
     
C                     This makes  d9_v = (fpar(4, n) * d3_v + d6_v) **  xp
                      d9_v =  d9_v * (fpar(4, n) * d3_v + d6_v)
     
C                     This makes  d4_p = log((fpar(4, n) * d3_v + d6_v)) * (fpar(4, n) * d3_v + d6_v) **  xp
                      d4_p =  d4_p *  d9_v
                        endif
                     else
C                  ((fpar(4, n) * d3_v + d6_v) == 0)

                        d9_v = (fpar(4, n) * d3_v + d6_v) **  xp
     
                        d3_p = 0.0d0
                call ehbfDO (12,fpar(4, n) * d3_v + d6_v, xp, d9_v, d
     +3_p, d4_p,g_ehfid,2077)
     
                     endif
     
                d11_v = d9_v * fpar(7, n)

                     if ( t .ne. 0.0d0 ) then
                        d14_v = t ** ( fpar(8, n)-2.0d0)
     
c                  d14_v = t ** ( fpar(8, n)-1)
                   d14_v =  d14_v * t
     
C                  d1_p =  fpar(8, n) * (t ** ( fpar(8, n)-1) )
                   d1_p =  fpar(8, n) *  d14_v
     
                        if ( t .lt. 0.0d0 ) then
c                     d14_v = t **  fpar(8, n)
                      d14_v =  d14_v * t
     
               call ehbfDO (12,t, fpar(8, n), d14_v, d1_p, d2_p,
     +g_ehfid,2098)
     
                        else
c                     (t > 0 here)
C                     Use  d2_p for scratch.
                      d2_p = log (t)
     
     
C                     This makes  d14_v = t **  fpar(8, n)
                      d14_v =  d14_v * t
     
C                     This makes  d2_p = log(t) * t **  fpar(8, n)
                      d2_p =  d2_p *  d14_v
                        endif
                     else
C                  (t == 0)

                        d14_v = t **  fpar(8, n)
     
                        d1_p = 0.0d0
                call ehbfDO (12,t, fpar(8, n), d14_v, d1_p, d2_p,
     +g_ehfid,2120)
     
                     endif
     
                d4_b = d11_v * d1_p
                d5_b = d11_v * d2_p
                d6_b = d14_v * fpar(7, n)
                d7_b = d14_v * d9_v
                d8_b = d6_b * d3_p
                d9_b = d6_b * d4_p
                d12_b = d8_b * d5_p
                d13_b = d8_b * d3_v
                d15_b = d8_b * fpar(4, n) * d6_p
                do g_i_ = 1, g_p_
                  g_fc(g_i_) = d5_b * g_fpar(g_i_, 8, n) + d4_b * g_t(g_
     *i_) + d7_b * g_fpar(g_i_, 7, n) + d9_b * g_xp(g_i_) + d12_b * g_d2
     *_w(g_i_) + d15_b * g_d1_w(g_i_) + d13_b * g_fpar(g_i_, 4, n)
                enddo
                fc = d11_v * d14_v
C--------
C
              else
C
C              6-PARAMETER TROE FORM
C
                d4_v = (-t) / fpar(5, n)
                d3_b = (-d4_v) / fpar(5, n)
                d4_b = -(1.0d0 / fpar(5, n))
                do g_i_ = 1, g_p_
                  g_d1_w(g_i_) = d3_b * g_fpar(g_i_, 5, n) + d4_b * g_t(
     *g_i_)
                enddo
                d1_w = d4_v
                d4_v = (-t) / fpar(6, n)
                d3_b = (-d4_v) / fpar(6, n)
                d4_b = -(1.0d0 / fpar(6, n))
                do g_i_ = 1, g_p_
                  g_d2_w(g_i_) = d3_b * g_fpar(g_i_, 6, n) + d4_b * g_t(
     *g_i_)
                enddo
                d2_w = d4_v
                d2_v = 1.0d0 - fpar(4, n)
                d4_v = exp(d1_w)
                     d2_p =  d4_v
     
                d7_v = exp(d2_w)
                     d1_p =  d7_v
     
                d6_b = fpar(4, n) * d1_p
                d9_b = d2_v * d2_p
                d4_b = d7_v + (-d4_v)
                do g_i_ = 1, g_p_
                  g_fcent(g_i_) = d6_b * g_d2_w(g_i_) + d9_b * g_d1_w(g_
     *i_) + d4_b * g_fpar(g_i_, 4, n)
                enddo
                fcent = d2_v * d4_v + fpar(4, n) * d7_v
C--------
C
C              7-PARAMETER TROE FORM
C
                if (ifop(n) .eq. 4) then
                  d4_v = (-fpar(7, n)) / t
                  d3_b = (-d4_v) / t
                  d4_b = -(1.0d0 / t)
                  do g_i_ = 1, g_p_
                    g_d1_w(g_i_) = d3_b * g_t(g_i_) + d4_b * g_fpar(g_i_
     *, 7, n)
                  enddo
                  d1_w = d4_v
                  d3_v = exp(d1_w)
                       d1_p =  d3_v
     
                  do g_i_ = 1, g_p_
                    g_fcent(g_i_) = d1_p * g_d1_w(g_i_) + g_fcent(g_i_)
                  enddo
                  fcent = fcent + d3_v
C--------
                endif
C
                d2_v = max (fcent, small)
     
                     if (fcent .gt.  small) then
                        d1_p = 1.0d0
                        d2_p = 0.0d0
                     else if (fcent .lt.  small) then
                        d1_p = 0.0d0
                        d2_p = 1.0d0
                     else
                        call ehbfDO (7,fcent, small, d2_v, d1_p, d2_p,
     +g_ehfid,2210)
                        d2_p = 1.0d0 -  d1_p
                     endif
     
     
                do g_i_ = 1, g_p_
                  g_d1_w(g_i_) = d1_p * g_fcent(g_i_)
                enddo
                d1_w = d2_v
                d2_v = log10 (d1_w)
                     d1_p = 1.0d0 / ( d1_w * log (10.0d0) )
     
                do g_i_ = 1, g_p_
                  g_fclog(g_i_) = d1_p * g_d1_w(g_i_)
                enddo
                fclog = d2_v
C--------
                d3_b = -dble(1.27)
                do g_i_ = 1, g_p_
                  g_xn(g_i_) = d3_b * g_fclog(g_i_)
                enddo
                xn = dble(0.75) - dble(1.27) * fclog
C--------
                d5_b = -dble(0.67)
                do g_i_ = 1, g_p_
                  g_cprlog(g_i_) = d5_b * g_fclog(g_i_) + g_prlog(g_i_)
                enddo
                cprlog = prlog - (dble(0.4) + dble(0.67) * fclog)
C--------
                d5_v = xn - dble(0.14) * cprlog
                d6_v = cprlog / d5_v
                d7_v = d6_v * d6_v
                     d1_p = 2.0d0 * d6_v
     
                d8_v = 1.0d0 + d7_v
                d9_v = fclog / d8_v
                d2_b = 1.0d0 / d8_v
                d5_b = (-d9_v) / d8_v * d1_p
                d7_b = d5_b * ((-d6_v) / d5_v)
                d6_b = d5_b * (1.0d0 / d5_v) + (-d7_b) * dble(0.14)
                do g_i_ = 1, g_p_
                  g_flog(g_i_) = d7_b * g_xn(g_i_) + d6_b * g_cprlog(g_i
     *_) + d2_b * g_fclog(g_i_)
                enddo
                flog = d9_v
C--------
                d2_v = dble(10.0) **  flog
                     if ( dble(10.0) .gt. 0.0d0 ) then
                        
C                  Use  d1_p to store `log(dble(10.0))'
                   d1_p = log (dble(10.0))
     
                        d1_p =  d1_p *  d2_v
                     else
C                  ( dble(10.0) <= 0 )

                        
            if ( (dble(10.0) .eq. 0.0d0) .and. ( flog .ne. 0.0d0)) then
                           d1_p = 0.0d0
                        else
                    call ehbfDO (11,dble(10.0), flog, d2_v, 0.0d0, d1_
     +p,g_ehfid,2274)
                        endif
     
                     endif
     
                do g_i_ = 1, g_p_
                  g_fc(g_i_) = d1_p * g_flog(g_i_)
                enddo
                fc = d2_v
C--------
              endif
              do g_i_ = 1, g_p_
                g_pcor(g_i_) = fc * g_pcor(g_i_) + pcor * g_fc(g_i_)
              enddo
              pcor = fc * pcor
C--------
            endif
C
            do g_i_ = 1, g_p_
              g_rkft(g_i_, ifal(n)) = rkft(ifal(n)) * g_pcor(g_i_) + pco
     *r * g_rkft(g_i_, ifal(n))
            enddo
            rkft(ifal(n)) = rkft(ifal(n)) * pcor
C--------
            do g_i_ = 1, g_p_
              g_rkrt(g_i_, ifal(n)) = rkrt(ifal(n)) * g_pcor(g_i_) + pco
     *r * g_rkrt(g_i_, ifal(n))
            enddo
            rkrt(ifal(n)) = rkrt(ifal(n)) * pcor
C--------
90          continue
99809     continue
        endif
C
C     Multiply by the product of reactants and product of products
C     PAR(4,I) is a perturbation factor
C
        do 99808 i = 1, ii
          d3_v = rkft(i) * ctb(i)
          d4_b = par(4, i) * ctb(i)
          d5_b = par(4, i) * rkft(i)
          do g_i_ = 1, g_p_
            g_rkft(g_i_, i) = d3_v * g_par(g_i_, 4, i) + d5_b * g_ctb(g_
     *i_, i) + d4_b * g_rkft(g_i_, i)
          enddo
          rkft(i) = d3_v * par(4, i)
C--------
          d3_v = rkrt(i) * ctb(i)
          d4_b = par(4, i) * ctb(i)
          d5_b = par(4, i) * rkrt(i)
          do g_i_ = 1, g_p_
            g_rkrt(g_i_, i) = d3_v * g_par(g_i_, 4, i) + d5_b * g_ctb(g_
     *i_, i) + d4_b * g_rkrt(g_i_, i)
          enddo
          rkrt(i) = d3_v * par(4, i)
C--------
C
          if (nu(1, i) .ne. 0) then

                 if (  iabs(nu(1, i)) .ne. 0 ) then
                    d3_v = c(nunk(1, i)) ** ( iabs(nu(1, i))-1)
                    d1_p =  iabs(nu(1, i)) *  d3_v
                    d3_v =  d3_v * c(nunk(1, i))
                 else
C              Maybe this should be  d3_v = c(nunk(1, i)) **  iabs(nu(1, i))
               d3_v = 1.0d0
                    d1_p = 0.0d0
                 endif
     
     
     
            i1_b = 0
            d4_b = rkft(i) * d1_p
            do g_i_ = 1, g_p_
              g_rkf(g_i_, i) = d4_b * g_c(g_i_, nunk(1, i)) + d3_v * g_r
     *kft(g_i_, i)
            enddo
            rkf(i) = rkft(i) * d3_v
C--------

                 if (  nu(4, i) .ne. 0 ) then
                    d3_v = c(nunk(4, i)) ** ( nu(4, i)-1)
                    d1_p =  nu(4, i) *  d3_v
                    d3_v =  d3_v * c(nunk(4, i))
                 else
C              Maybe this should be  d3_v = c(nunk(4, i)) **  nu(4, i)
               d3_v = 1.0d0
                    d1_p = 0.0d0
                 endif
     
     
     
            d4_b = rkrt(i) * d1_p
            do g_i_ = 1, g_p_
              g_rkr(g_i_, i) = d4_b * g_c(g_i_, nunk(4, i)) + d3_v * g_r
     *krt(g_i_, i)
            enddo
            rkr(i) = rkrt(i) * d3_v
C--------
            if (nunk(2, i) .ne. 0) then

                   if (  iabs(nu(2, i)) .ne. 0 ) then
                      d3_v = c(nunk(2, i)) ** ( iabs(nu(2, i))-1)
                      d1_p =  iabs(nu(2, i)) *  d3_v
                      d3_v =  d3_v * c(nunk(2, i))
                   else
C                Maybe this should be  d3_v = c(nunk(2, i)) **  iabs(nu(2, i))
                 d3_v = 1.0d0
                      d1_p = 0.0d0
                   endif
    
     
     
              i1_b = 0
              d4_b = rkf(i) * d1_p
              do g_i_ = 1, g_p_
                g_rkf(g_i_, i) = d4_b * g_c(g_i_, nunk(2, i)) + d3_v * g
     *_rkf(g_i_, i)
              enddo
              rkf(i) = rkf(i) * d3_v
C--------
              if (nunk(3, i) .ne. 0) then

                     if (  iabs(nu(3, i)) .ne. 0 ) then
                        d3_v = c(nunk(3, i)) ** ( iabs(nu(3, i))-1)
                        d1_p =  iabs(nu(3, i)) *  d3_v
                        d3_v =  d3_v * c(nunk(3, i))
                     else
C                  Maybe this should be  d3_v = c(nunk(3, i)) **  iabs(nu(3, i))
                   d3_v = 1.0d0
                        d1_p = 0.0d0
                     endif
     
     
     
                i1_b = 0
                d4_b = rkf(i) * d1_p
                do g_i_ = 1, g_p_
                  g_rkf(g_i_, i) = d4_b * g_c(g_i_, nunk(3, i)) + d3_v *
     * g_rkf(g_i_, i)
                enddo
                rkf(i) = rkf(i) * d3_v
C--------
              endif
            endif
            if (nunk(5, i) .ne. 0) then

                   if (  nu(5, i) .ne. 0 ) then
                      d3_v = c(nunk(5, i)) ** ( nu(5, i)-1)
                      d1_p =  nu(5, i) *  d3_v
                      d3_v =  d3_v * c(nunk(5, i))
                   else
C                Maybe this should be  d3_v = c(nunk(5, i)) **  nu(5, i)
                 d3_v = 1.0d0
                      d1_p = 0.0d0
                   endif
     
     
     
              d4_b = rkr(i) * d1_p
              do g_i_ = 1, g_p_
                g_rkr(g_i_, i) = d4_b * g_c(g_i_, nunk(5, i)) + d3_v * g
     *_rkr(g_i_, i)
              enddo
              rkr(i) = rkr(i) * d3_v
C--------
              if (nunk(6, i) .ne. 0) then

                     if (  nu(6, i) .ne. 0 ) then
                        d3_v = c(nunk(6, i)) ** ( nu(6, i)-1)
                        d1_p =  nu(6, i) *  d3_v
                        d3_v =  d3_v * c(nunk(6, i))
                     else
C                  Maybe this should be  d3_v = c(nunk(6, i)) **  nu(6, i)
                   d3_v = 1.0d0
                        d1_p = 0.0d0
                     endif
     
     
     
                d4_b = rkr(i) * d1_p
                do g_i_ = 1, g_p_
                  g_rkr(g_i_, i) = d4_b * g_c(g_i_, nunk(6, i)) + d3_v *
     * g_rkr(g_i_, i)
                enddo
                rkr(i) = rkr(i) * d3_v
C--------
              endif
            endif
          endif
150       continue
99808   continue
C
        do 99807 n = 1, nrnu
          i = irnu(n)
          d3_v = abs(rnu(1, n))
     
               if (rnu(1, n) .gt. 0.0d0) then
                  d3_p =  1.0d0
               else if (rnu(1, n) .lt. 0.0d0) then
                  d3_p = -1.0d0
               else
                  call ehufDO (3,rnu(1, n), d3_v, d3_p,
     +g_ehfid,2478)
               endif
     

               if ( c(nunk(1, i)) .ne. 0.0d0 ) then
                  d4_v = c(nunk(1, i)) ** ( d3_v-2.0d0)
     
c            d4_v = c(nunk(1, i)) ** ( d3_v-1)
             d4_v =  d4_v * c(nunk(1, i))
     
C            d1_p =  d3_v * (c(nunk(1, i)) ** ( d3_v-1) )
             d1_p =  d3_v *  d4_v
     
                  if ( c(nunk(1, i)) .lt. 0.0d0 ) then
c               d4_v = c(nunk(1, i)) **  d3_v
                d4_v =  d4_v * c(nunk(1, i))
     
              call ehbfDO (12,c(nunk(1, i)), d3_v, d4_v, d1_p, d2_p,
     +g_ehfid,2497)
     
                  else
c               (c(nunk(1, i)) > 0 here)
C               Use  d2_p for scratch.
                d2_p = log (c(nunk(1, i)))
     
     
C               This makes  d4_v = c(nunk(1, i)) **  d3_v
                d4_v =  d4_v * c(nunk(1, i))
     
C               This makes  d2_p = log(c(nunk(1, i))) * c(nunk(1, i)) **  d3_v
                d2_p =  d2_p *  d4_v
                  endif
               else
C            (c(nunk(1, i)) == 0)

                  d4_v = c(nunk(1, i)) **  d3_v
     
                  d1_p = 0.0d0
                  call ehbfDO (12,c(nunk(1, i)), d3_v, d4_v, d1_p, d2_p,
     +g_ehfid,2519)
     
               endif
     
          d4_b = d2_p * d3_p
          do g_i_ = 1, g_p_
            g_c1(g_i_) = d4_b * g_rnu(g_i_, 1, n) + d1_p * g_c(g_i_, nun
     *k(1, i))
          enddo
          c1 = d4_v
C--------

               if ( c(nunk(4, i)) .ne. 0.0d0 ) then
                  d3_v = c(nunk(4, i)) ** ( rnu(4, n)-2.0d0)
     
c            d3_v = c(nunk(4, i)) ** ( rnu(4, n)-1)
             d3_v =  d3_v * c(nunk(4, i))
     
C            d1_p =  rnu(4, n) * (c(nunk(4, i)) ** ( rnu(4, n)-1) )
             d1_p =  rnu(4, n) *  d3_v
     
                  if ( c(nunk(4, i)) .lt. 0.0d0 ) then
c               d3_v = c(nunk(4, i)) **  rnu(4, n)
                d3_v =  d3_v * c(nunk(4, i))
     
              call ehbfDO (12,c(nunk(4, i)), rnu(4, n), d3_v, d1_p, d2
     +_p,g_ehfid,2547)
     
                  else
c               (c(nunk(4, i)) > 0 here)
C               Use  d2_p for scratch.
                d2_p = log (c(nunk(4, i)))
     
     
C               This makes  d3_v = c(nunk(4, i)) **  rnu(4, n)
                d3_v =  d3_v * c(nunk(4, i))
     
C               This makes  d2_p = log(c(nunk(4, i))) * c(nunk(4, i)) **  rnu(4, n)
                d2_p =  d2_p *  d3_v
                  endif
               else
C            (c(nunk(4, i)) == 0)

                  d3_v = c(nunk(4, i)) **  rnu(4, n)
     
                  d1_p = 0.0d0
            call ehbfDO (12,c(nunk(4, i)), rnu(4, n), d3_v, d1_p, d2_p,
     +g_ehfid,2570)
     
               endif
     
          do g_i_ = 1, g_p_
            g_c4(g_i_) = d2_p * g_rnu(g_i_, 4, n) + d1_p * g_c(g_i_, nun
     *k(4, i))
          enddo
          c4 = d3_v
C--------
          do g_i_ = 1, g_p_
            g_rkf(g_i_, i) = rkft(i) * g_c1(g_i_) + c1 * g_rkft(g_i_, i)
          enddo
          rkf(i) = rkft(i) * c1
C--------
          do g_i_ = 1, g_p_
            g_rkr(g_i_, i) = rkrt(i) * g_c4(g_i_) + c4 * g_rkrt(g_i_, i)
          enddo
          rkr(i) = rkrt(i) * c4
C--------
          if (nunk(2, i) .ne. 0) then
            d3_v = abs(rnu(2, n))
     
                 if (rnu(2, n) .gt. 0.0d0) then
                    d3_p =  1.0d0
                 else if (rnu(2, n) .lt. 0.0d0) then
                    d3_p = -1.0d0
                 else
                    call ehufDO (3,rnu(2, n), d3_v, d3_p,
     +g_ehfid,2600)
                 endif
     

                 if ( c(nunk(2, i)) .ne. 0.0d0 ) then
                    d4_v = c(nunk(2, i)) ** ( d3_v-2.0d0)
     
c              d4_v = c(nunk(2, i)) ** ( d3_v-1)
               d4_v =  d4_v * c(nunk(2, i))
     
C              d1_p =  d3_v * (c(nunk(2, i)) ** ( d3_v-1) )
               d1_p =  d3_v *  d4_v
     
                    if ( c(nunk(2, i)) .lt. 0.0d0 ) then
c                 d4_v = c(nunk(2, i)) **  d3_v
                  d4_v =  d4_v * c(nunk(2, i))
     
              call ehbfDO (12,c(nunk(2, i)), d3_v, d4_v, d1_p, d2_p,
     +g_ehfid,2620)
     
                    else
c                 (c(nunk(2, i)) > 0 here)
C                 Use  d2_p for scratch.
                  d2_p = log (c(nunk(2, i)))
    
     
C                 This makes  d4_v = c(nunk(2, i)) **  d3_v
                  d4_v =  d4_v * c(nunk(2, i))
     
C                 This makes  d2_p = log(c(nunk(2, i))) * c(nunk(2, i)) **  d3_v
                  d2_p =  d2_p *  d4_v
                    endif
                 else
C              (c(nunk(2, i)) == 0)

                    d4_v = c(nunk(2, i)) **  d3_v
     
                    d1_p = 0.0d0
                call ehbfDO (12,c(nunk(2, i)), d3_v, d4_v, d1_p, d2_p,
     +g_ehfid,2642)
     
                 endif
     
            d4_b = d2_p * d3_p
            do g_i_ = 1, g_p_
              g_c2(g_i_) = d4_b * g_rnu(g_i_, 2, n) + d1_p * g_c(g_i_, n
     *unk(2, i))
            enddo
            c2 = d4_v
C--------
            do g_i_ = 1, g_p_
              g_rkf(g_i_, i) = rkf(i) * g_c2(g_i_) + c2 * g_rkf(g_i_, i)
            enddo
            rkf(i) = rkf(i) * c2
C--------
            if (nunk(3, i) .ne. 0) then
              d3_v = abs(rnu(3, n))
     
                   if (rnu(3, n) .gt. 0.0d0) then
                      d3_p =  1.0d0
                   else if (rnu(3, n) .lt. 0.0d0) then
                      d3_p = -1.0d0
                  else
                      call ehufDO (3,rnu(3, n), d3_v, d3_p,
     +g_ehfid,2668)
                   endif
     

                   if ( c(nunk(3, i)) .ne. 0.0d0 ) then
                      d4_v = c(nunk(3, i)) ** ( d3_v-2.0d0)
     
c                d4_v = c(nunk(3, i)) ** ( d3_v-1)
                 d4_v =  d4_v * c(nunk(3, i))
     
C                d1_p =  d3_v * (c(nunk(3, i)) ** ( d3_v-1) )
                 d1_p =  d3_v *  d4_v
     
                      if ( c(nunk(3, i)) .lt. 0.0d0 ) then
c                   d4_v = c(nunk(3, i)) **  d3_v
                    d4_v =  d4_v * c(nunk(3, i))
     
              call ehbfDO (12,c(nunk(3, i)), d3_v, d4_v, d1_p, d2_
     +p,g_ehfid,2688)
     
                      else
c                   (c(nunk(3, i)) > 0 here)
C                   Use  d2_p for scratch.
                    d2_p = log (c(nunk(3, i)))
     
     
C                   This makes  d4_v = c(nunk(3, i)) **  d3_v
                    d4_v =  d4_v * c(nunk(3, i))
     
C                   This makes  d2_p = log(c(nunk(3, i))) * c(nunk(3, i)) **  d3_v
                    d2_p =  d2_p *  d4_v
                      endif
                   else
C                (c(nunk(3, i)) == 0)

                      d4_v = c(nunk(3, i)) **  d3_v
     
                      d1_p = 0.0d0
             call ehbfDO (12,c(nunk(3, i)), d3_v, d4_v, d1_p, d2_p,
     +g_ehfid,2710)
     
                   endif
     
              d4_b = d2_p * d3_p
              do g_i_ = 1, g_p_
                g_c3(g_i_) = d4_b * g_rnu(g_i_, 3, n) + d1_p * g_c(g_i_,
     * nunk(3, i))
              enddo
              c3 = d4_v
C--------
              do g_i_ = 1, g_p_
                g_rkf(g_i_, i) = rkf(i) * g_c3(g_i_) + c3 * g_rkf(g_i_, 
     *i)
              enddo
              rkf(i) = rkf(i) * c3
C--------
            endif
          endif
          if (nunk(5, i) .ne. 0) then

                 if ( c(nunk(5, i)) .ne. 0.0d0 ) then
                    d3_v = c(nunk(5, i)) ** ( rnu(5, n)-2.0d0)
     
c              d3_v = c(nunk(5, i)) ** ( rnu(5, n)-1)
               d3_v =  d3_v * c(nunk(5, i))
     
C              d1_p =  rnu(5, n) * (c(nunk(5, i)) ** ( rnu(5, n)-1) )
               d1_p =  rnu(5, n) *  d3_v
     
                    if ( c(nunk(5, i)) .lt. 0.0d0 ) then
c                 d3_v = c(nunk(5, i)) **  rnu(5, n)
                  d3_v =  d3_v * c(nunk(5, i))
     
               call ehbfDO (12,c(nunk(5, i)), rnu(5, n), d3_v, d1_p, 
     +d2_p,g_ehfid,2747)
     
                    else
c                 (c(nunk(5, i)) > 0 here)
C                 Use  d2_p for scratch.
                  d2_p = log (c(nunk(5, i)))
     
     
C                 This makes  d3_v = c(nunk(5, i)) **  rnu(5, n)
                  d3_v =  d3_v * c(nunk(5, i))
     
C                 This makes  d2_p = log(c(nunk(5, i))) * c(nunk(5, i)) **  rnu(5, n)
                  d2_p =  d2_p *  d3_v
                    endif
                 else
C              (c(nunk(5, i)) == 0)

                    d3_v = c(nunk(5, i)) **  rnu(5, n)
     
                    d1_p = 0.0d0
             call ehbfDO (12,c(nunk(5, i)), rnu(5, n), d3_v, d1_p, d2_
     +p,g_ehfid,2770)
     
                 endif
     
            do g_i_ = 1, g_p_
              g_c5(g_i_) = d2_p * g_rnu(g_i_, 5, n) + d1_p * g_c(g_i_, n
     *unk(5, i))
            enddo
            c5 = d3_v
C--------
            do g_i_ = 1, g_p_
              g_rkr(g_i_, i) = rkr(i) * g_c5(g_i_) + c5 * g_rkr(g_i_, i)
            enddo
            rkr(i) = rkr(i) * c5
C--------
            if (nunk(6, i) .ne. 0) then

                   if ( c(nunk(6, i)) .ne. 0.0d0 ) then
                      d3_v = c(nunk(6, i)) ** ( rnu(6, n)-2.0d0)
     
c                d3_v = c(nunk(6, i)) ** ( rnu(6, n)-1)
                 d3_v =  d3_v * c(nunk(6, i))
     
C                d1_p =  rnu(6, n) * (c(nunk(6, i)) ** ( rnu(6, n)-1) )
                 d1_p =  rnu(6, n) *  d3_v
     
                      if ( c(nunk(6, i)) .lt. 0.0d0 ) then
c                   d3_v = c(nunk(6, i)) **  rnu(6, n)
                    d3_v =  d3_v * c(nunk(6, i))
     
                call ehbfDO (12,c(nunk(6, i)), rnu(6, n), d3_v, d1_p
     +, d2_p,g_ehfid,2803)
     
                      else
c                   (c(nunk(6, i)) > 0 here)
C                   Use  d2_p for scratch.
                    d2_p = log (c(nunk(6, i)))
     
     
C                   This makes  d3_v = c(nunk(6, i)) **  rnu(6, n)
                    d3_v =  d3_v * c(nunk(6, i))
     
C                   This makes  d2_p = log(c(nunk(6, i))) * c(nunk(6, i)) **  rnu(6, n)
                    d2_p =  d2_p *  d3_v
                      endif
                   else
C                (c(nunk(6, i)) == 0)

                      d3_v = c(nunk(6, i)) **  rnu(6, n)
     
                      d1_p = 0.0d0
              call ehbfDO (12,c(nunk(6, i)), rnu(6, n), d3_v, d1_p, d
     +2_p,g_ehfid,2826)
     
                   endif
     
              do g_i_ = 1, g_p_
                g_c6(g_i_) = d2_p * g_rnu(g_i_, 6, n) + d1_p * g_c(g_i_,
     * nunk(6, i))
              enddo
              c6 = d3_v
C--------
              do g_i_ = 1, g_p_
                g_rkr(g_i_, i) = rkr(i) * g_c6(g_i_) + c6 * g_rkr(g_i_, 
     *i)
              enddo
              rkr(i) = rkr(i) * c6
C--------
            endif
          endif
160       continue
99807   continue
C
        do 99805 n = 1, nord
          i = iord(n)
          do g_i_ = 1, g_p_
            g_rkf(g_i_, i) = g_rkft(g_i_, i)
          enddo
          rkf(i) = rkft(i)
C--------
          do g_i_ = 1, g_p_
            g_rkr(g_i_, i) = g_rkrt(g_i_, i)
          enddo
          rkr(i) = rkrt(i)
C--------
C
          do 99806 l = 1, mxord
            nk = kord(l, n)
            if (nk .lt. 0) then
              nk = iabs(nk)

                   if ( c(nk) .ne. 0.0d0 ) then
                      d3_v = c(nk) ** ( rord(l, n)-2.0d0)
     
c                d3_v = c(nk) ** ( rord(l, n)-1)
                 d3_v =  d3_v * c(nk)
     
C                d1_p =  rord(l, n) * (c(nk) ** ( rord(l, n)-1) )
                 d1_p =  rord(l, n) *  d3_v
     
                      if ( c(nk) .lt. 0.0d0 ) then
c                   d3_v = c(nk) **  rord(l, n)
                    d3_v =  d3_v * c(nk)
     
              call ehbfDO (12,c(nk), rord(l, n), d3_v, d1_p, d2_p,
     +g_ehfid,2881)
     
                      else
c                   (c(nk) > 0 here)
C                   Use  d2_p for scratch.
                    d2_p = log (c(nk))
     
     
C                   This makes  d3_v = c(nk) **  rord(l, n)
                    d3_v =  d3_v * c(nk)
     
C                   This makes  d2_p = log(c(nk)) * c(nk) **  rord(l, n)
                    d2_p =  d2_p *  d3_v
                      endif
                   else
C                (c(nk) == 0)

                      d3_v = c(nk) **  rord(l, n)
     
                      d1_p = 0.0d0
              call ehbfDO (12,c(nk), rord(l, n), d3_v, d1_p, d2_p,
     +g_ehfid,2903)
     
                   endif
     
              do g_i_ = 1, g_p_
                g_cnk(g_i_) = d2_p * g_rord(g_i_, l, n) + d1_p * g_c(g_i
     *_, nk)
              enddo
              cnk = d3_v
C--------
              do g_i_ = 1, g_p_
                g_rkf(g_i_, i) = rkf(i) * g_cnk(g_i_) + cnk * g_rkf(g_i_
     *, i)
              enddo
              rkf(i) = rkf(i) * cnk
C--------
            else
              if (nk .gt. 0) then

                     if ( c(nk) .ne. 0.0d0 ) then
                        d3_v = c(nk) ** ( rord(l, n)-2.0d0)
     
c                  d3_v = c(nk) ** ( rord(l, n)-1)
                   d3_v =  d3_v * c(nk)
     
C                  d1_p =  rord(l, n) * (c(nk) ** ( rord(l, n)-1) )
                   d1_p =  rord(l, n) *  d3_v
     
                        if ( c(nk) .lt. 0.0d0 ) then
c                     d3_v = c(nk) **  rord(l, n)
                      d3_v =  d3_v * c(nk)
     
                   call ehbfDO (12,c(nk), rord(l, n), d3_v, d1_p, d2_
     +p,g_ehfid,2938)
     
                        else
c                     (c(nk) > 0 here)
C                     Use  d2_p for scratch.
                      d2_p = log (c(nk))
     
     
C                     This makes  d3_v = c(nk) **  rord(l, n)
                      d3_v =  d3_v * c(nk)
     
C                     This makes  d2_p = log(c(nk)) * c(nk) **  rord(l, n)
                      d2_p =  d2_p *  d3_v
                        endif
                     else
C                  (c(nk) == 0)

                        d3_v = c(nk) **  rord(l, n)
     
                        d1_p = 0.0d0
                call ehbfDO (12,c(nk), rord(l, n), d3_v, d1_p, d2_p,
     +g_ehfid,2960)
     
                     endif
     
                do g_i_ = 1, g_p_
                  g_cnk(g_i_) = d2_p * g_rord(g_i_, l, n) + d1_p * g_c(g
     *_i_, nk)
                enddo
                cnk = d3_v
C--------
                do g_i_ = 1, g_p_
                  g_rkr(g_i_, i) = rkr(i) * g_cnk(g_i_) + cnk * g_rkr(g_
     *i_, i)
                enddo
                rkr(i) = rkr(i) * cnk
C--------
              endif
            endif
190         continue
99806     continue
200       continue
99805   continue
C
        return
      end  subroutine g_ckratx
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckrdex' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckrhex' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckrhoc' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckrhox' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      subroutine g_ckrhoy(g_p_, p, t, g_t, ldg_t, y, g_y, ldg_y, ickwrk,
     * rckwrk, g_rckwrk, ldg_rckwrk, rho, g_rho, ldg_rho)
C
C  START PROLOGUE
C
C  SUBROUTINE CKRHOY (P, T, Y, ICKWRK, RCKWRK, RHO)
C     Returns the mass density of the gas mixture given the pressure,
C     temperature and mass fractions;  see Eq. (2).
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
C     RHO    - Mass density.
C                   cgs units - gm/cm**3
C                   Data type - real scalar
C
C  END PROLOGUE
C
C*****precision > double
        implicit double precision (a-h, o-z), integer (i-n)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
        dimension ickwrk(*), rckwrk(*), y(*)
C       Start of include file /home6/laniu/adifor/adiforMatab-Pre/ckstrt.h
C
C
C     include file for cklib.f
C
        common /ckstrt/ nmm, nkk, nii, mxsp, mxtb, mxtp, ncp, ncp1, ncp2
     *, ncp2t, npar, nlar, nfar, nlan, nfal, nrev, nthb, nrlt, nwl, neim
     *, njan, njar, nft1, nf1r, nexc, nrnu, nord, mxord, icmm, ickk, icn
     *c, icph, icch, icnt, icnu, icnk, icns, icnr, iclt, icrl, icrv, icw
     *l, icfl, icfo, ickf, ictb, ickn, ickt, icei, ictd, icjn, icf1, ice
     *x, icrnu, icord, ickor, ncaw, ncwt, nctt, ncaa, ncco, ncrv, nclt, 
     *ncrl, ncfl, nckt, ncwl, ncjn, ncf1, ncex, ncru, ncrc, ncpa, nckf, 
     *nckr, ncrnu, nckor, nck1, nck2, nck3, nck4, nci1, nci2, nci3, nci4
C
C     END include file for cklib.f
C
C

        integer g_i_, g_p_, ldg_y, ldg_rckwrk, ldg_rho, ldg_t
        double precision d6_b, d6_v, d5_b, d4_b, d3_b, d2_b, d3_v, d5_v,
     * d4_v, g_sumyow(g_p_)
        double precision g_y(ldg_y, *), g_rckwrk(ldg_rckwrk, *), g_rho(l
     *dg_rho), g_t(ldg_t)
        integer g_ehfid
!Laniu        save g_sumyow
        intrinsic dble
        data g_ehfid /0/
C
claniu        call ehsfid(g_ehfid, 'ckrhoy','g_cklib.f')
C

        do g_i_ = 1, g_p_
          g_sumyow(g_i_) = 0.0d0
        enddo
        sumyow = 0.0d0
C--------
        do 99801 k = 1, nkk
          d4_v = y(k) / rckwrk(ncwt + k - 1)
          d4_b = 1.0d0 / rckwrk(ncwt + k - 1)
          d5_b = (-d4_v) / rckwrk(ncwt + k - 1)
          do g_i_ = 1, g_p_
            g_sumyow(g_i_) = d5_b * g_rckwrk(g_i_, ncwt + k - 1) + d4_b 
     ** g_y(g_i_, k) + g_sumyow(g_i_)
          enddo
          sumyow = sumyow + d4_v
C--------
150       continue
99801   continue
        d3_v = sumyow * t
        d5_v = d3_v * rckwrk(ncru)
        d6_v = p / d5_v
        d2_b = (-d6_v) / d5_v
        d3_b = d2_b * rckwrk(ncru)
        d4_b = d2_b * d3_v
        d5_b = d3_b * t
        d6_b = d3_b * sumyow
        do g_i_ = 1, g_p_
          g_rho(g_i_) = d4_b * g_rckwrk(g_i_, ncru) + d6_b * g_t(g_i_) +
     * d5_b * g_sumyow(g_i_)
        enddo
        rho = d6_v
C--------
        return
      end subroutine g_ckrhoy
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckrp' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'cksave' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'cksbml' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'cksbms' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      subroutine g_cksmh(g_p_, t, g_t, ldg_t, ickwrk, rckwrk, g_rckwrk, 
     *ldg_rckwrk, smh, g_smh, ldg_smh)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSMH  (T, ICKWRK, RCKWRK, SMH)*
C     Returns the array of entropies minus enthalpies for the species.
C     It is normally not called directly by the user.
C
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT
C     SMH    - Entropy minus enthalpy for the species,
C              SMH(K) = S(K)/R - H(K)/RT.
C                   cgs units - none
C                   Data type - real array
C                   Dimension SMH(*) at least KK, the total number of
C                   species.
C
C  END PROLOGUE
C
C*****precision > double
        implicit double precision (a-h, o-z), integer (i-n)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
        dimension ickwrk(*), rckwrk(*), smh(*), tn(10)
C       Start of include file /home6/laniu/adifor/adiforMatab-Pre/ckstrt.h
C
C
C     include file for cklib.f
C
        common /ckstrt/ nmm, nkk, nii, mxsp, mxtb, mxtp, ncp, ncp1, ncp2
     *, ncp2t, npar, nlar, nfar, nlan, nfal, nrev, nthb, nrlt, nwl, neim
     *, njan, njar, nft1, nf1r, nexc, nrnu, nord, mxord, icmm, ickk, icn
     *c, icph, icch, icnt, icnu, icnk, icns, icnr, iclt, icrl, icrv, icw
     *l, icfl, icfo, ickf, ictb, ickn, ickt, icei, ictd, icjn, icf1, ice
     *x, icrnu, icord, ickor, ncaw, ncwt, nctt, ncaa, ncco, ncrv, nclt, 
     *ncrl, ncfl, nckt, ncwl, ncjn, ncf1, ncex, ncru, ncrc, ncpa, nckf, 
     *nckr, ncrnu, nckor, nck1, nck2, nck3, nck4, nci1, nci2, nci3, nci4
C
C     END include file for cklib.f
C
C

        integer g_i_, g_p_, i1_b, ldg_t, ldg_rckwrk, ldg_smh
        double precision d1_p, d2_v, d6_v, d5_b, d4_b, d3_b, g_tn(g_p
     *_, 10), g_t(ldg_t), g_sum(g_p_), g_rckwrk(ldg_rckwrk, *)
        double precision g_smh(ldg_smh, *)
        integer g_ehfid
 !Laniu       save g_tn, g_sum
        intrinsic dble
        data g_ehfid /0/
C
claniu        call ehsfid(g_ehfid, 'cksmh','g_cklib.f')
C

        d2_v = log(t)
             d1_p = 1.0d0 / t
     
        do g_i_ = 1, g_p_
          g_tn(g_i_, 1) = d1_p * g_t(g_i_)
        enddo
        tn(1) = d2_v - 1.0d0
C--------
        do 99798 n = 2, ncp

               if ( ( n - 1) .ne. 0 ) then
                  d2_v = t ** (( n - 1)-1)
                  d1_p = ( n - 1) *  d2_v
                  d2_v =  d2_v * t
               else
C            Maybe this should be  d2_v = t ** ( n - 1)
             d2_v = 1.0d0
                  d1_p = 0.0d0
               endif
     
     
     
          i1_b = 0
          d3_b = 1.0d0 / dble((n - 1) * n) * d1_p
          do g_i_ = 1, g_p_
            g_tn(g_i_, n) = d3_b * g_t(g_i_)
          enddo
          tn(n) = d2_v / dble((n - 1) * n)
C--------
150       continue
99798   continue
C
        do 99795 k = 1, nkk
          l = 1
          do 99797 n = 2, ickwrk(icnt + k - 1) - 1
            temp = rckwrk(nctt + (k - 1) * mxtp + n - 1)
            if (t .gt. temp) then
              l = l + 1
            endif
220         continue
99797     continue
C
          na1 = ncaa + (l - 1) * ncp2 + (k - 1) * ncp2t
          do g_i_ = 1, g_p_
            g_sum(g_i_) = 0.0d0
          enddo
          sum = 0.0d0
C--------
          do 99796 n = 1, ncp
            do g_i_ = 1, g_p_
              g_sum(g_i_) = tn(n) * g_rckwrk(g_i_, na1 + n - 1) + rckwrk
     *(na1 + n - 1) * g_tn(g_i_, n) + g_sum(g_i_)
            enddo
            sum = sum + tn(n) * rckwrk(na1 + n - 1)
C--------
225         continue
99796     continue
          d6_v = rckwrk(na1 + ncp1 - 1) / t
          d4_b = -(1.0d0 / t)
          d5_b = -((-d6_v) / t)
          do g_i_ = 1, g_p_
            g_smh(g_i_, k) = d5_b * g_t(g_i_) + d4_b * g_rckwrk(g_i_, na
     *1 + ncp1 - 1) + g_rckwrk(g_i_, na1 + ncp2 - 1) + g_sum(g_i_)
          enddo
          smh(k) = sum + rckwrk(na1 + ncp2 - 1) - d6_v
C--------
C
250       continue
99795   continue
        return
      end subroutine g_cksmh
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'cksml' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'cksms' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'cksnum' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'cksor' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'cksubs' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'cksyme' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'cksymr' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'cksyms' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckthb' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckubml' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckubms' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckuml' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckums' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckwc' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckwl' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckwt' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckwxp' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckwxr' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      subroutine g_ckwyp(g_p_, p, t, g_t, ldg_t, y, g_y, ldg_y, ickwrk, 
     *rckwrk, g_rckwrk, ldg_rckwrk, wdot, g_wdot, ldg_wdot)
C
C  START PROLOGUE
C
C  SUBROUTINE CKWYP  (P, T, Y, ICKWRK, RCKWRK, WDOT)
C     Returns the molar production rates of the species given the
C     pressure, temperature and mass fractions;  see Eq. (49).
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
C*****precision > double
        implicit double precision (a-h, o-z), integer (i-n)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
        dimension t_a(1)
        dimension ickwrk(*), rckwrk(*), y(*), wdot(*)
C       Start of include file /home6/laniu/adifor/adiforMatab-Pre/ckstrt.h
C
C
C     include file for cklib.f
C
        common /ckstrt/ nmm, nkk, nii, mxsp, mxtb, mxtp, ncp, ncp1, ncp2
     *, ncp2t, npar, nlar, nfar, nlan, nfal, nrev, nthb, nrlt, nwl, neim
     *, njan, njar, nft1, nf1r, nexc, nrnu, nord, mxord, icmm, ickk, icn
     *c, icph, icch, icnt, icnu, icnk, icns, icnr, iclt, icrl, icrv, icw
     *l, icfl, icfo, ickf, ictb, ickn, ickt, icei, ictd, icjn, icf1, ice
     *x, icrnu, icord, ickor, ncaw, ncwt, nctt, ncaa, ncco, ncrv, nclt, 
     *ncrl, ncfl, nckt, ncwl, ncjn, ncf1, ncex, ncru, ncrc, ncpa, nckf, 
     *nckr, ncrnu, nckor, nck1, nck2, nck3, nck4, nci1, nci2, nci3, nci4
C
C     END include file for cklib.f
C
C

        integer g_i_, g_p_, ldg_t, ldg_wdot, ldg_rckwrk, ldg_y
        double precision d4_b, g_t_a(g_p_, 1), g_t(ldg_t), g_wdot(ldg
     *_wdot, *), g_rop(g_p_), g_rckwrk(ldg_rckwrk, *), g_y(ldg_y, *)
        integer g_ehfid
!Laniu        save g_t_a, g_rop
        external g_ckratx
        external g_ckytcp
        external g_ckratt
        data g_ehfid /0/
C
claniu        call ehsfid(g_ehfid, 'ckwyp','g_cklib.f')
C

        do g_i_ = 1, g_p_
          g_t_a(g_i_, 1) = g_t(g_i_)
        enddo
        t_a(1) = t
C--------
        call g_ckratt(g_p_, rckwrk, g_rckwrk, ldg_rckwrk, ickwrk, nii, m
     *xsp, rckwrk(ncru), g_rckwrk(1, ncru), ldg_rckwrk, rckwrk(ncpa), g_
     *rckwrk(1, ncpa), ldg_rckwrk, t_a, g_t_a, g_p_, ickwrk(icns), ic
     *kwrk(icnu), ickwrk(icnk), npar + 1, rckwrk(ncco), g_rckwrk(1, ncco
     *), ldg_rckwrk, nrev, ickwrk(icrv), rckwrk(ncrv), g_rckwrk(1, ncrv)
     *, ldg_rckwrk, nlan, nlar, ickwrk(iclt), rckwrk(nclt), g_rckwrk(1, 
     *nclt), ldg_rckwrk, nrlt, ickwrk(icrl), rckwrk(ncrl), g_rckwrk(1, n
     *crl), ldg_rckwrk, rckwrk(nck1), g_rckwrk(1, nck1), ldg_rckwrk, nrn
     *u, ickwrk(icrnu), rckwrk(ncrnu), g_rckwrk(1, ncrnu), ldg_rckwrk, n
     *eim, ickwrk(icei), ickwrk(ictd), njan, njar, ickwrk(icjn), rckwrk(
     *ncjn), g_rckwrk(1, ncjn), ldg_rckwrk, nft1, nf1r, ickwrk(icf1), rc
     *kwrk(ncf1), g_rckwrk(1, ncf1), ldg_rckwrk, rckwrk(nckf), g_rckwrk(
     *1, nckf), ldg_rckwrk, rckwrk(nckr), g_rckwrk(1, nckr), ldg_rckwrk,
     * rckwrk(nci1), g_rckwrk(1, nci1), ldg_rckwrk)
C
        call g_ckytcp(g_p_, p, t, g_t, ldg_t, y, g_y, ldg_y, ickwrk, rck
     *wrk, g_rckwrk, ldg_rckwrk, rckwrk(nck1), g_rckwrk(1, nck1), ldg_rc
     *kwrk)
C
        call g_ckratx(g_p_, nii, nkk, mxsp, mxtb, t, g_t, ldg_t, rckwrk(
     *nck1), g_rckwrk(1, nck1), ldg_rckwrk, ickwrk(icnu), ickwrk(icnk), 
     *npar + 1, rckwrk(ncco), g_rckwrk(1, ncco), ldg_rckwrk, nfal, ickwr
     *k(icfl), ickwrk(icfo), ickwrk(ickf), nfar, rckwrk(ncfl), g_rckwrk(
     *1, ncfl), ldg_rckwrk, nthb, ickwrk(ictb), ickwrk(ickn), rckwrk(nck
     *t), g_rckwrk(1, nckt), ldg_rckwrk, ickwrk(ickt), rckwrk(nckf), g_r
     *ckwrk(1, nckf), ldg_rckwrk, rckwrk(nckr), g_rckwrk(1, nckr), ldg_r
     *ckwrk, rckwrk(nci1), g_rckwrk(1, nci1), ldg_rckwrk, rckwrk(nci2), 
     *g_rckwrk(1, nci2), ldg_rckwrk, rckwrk(nci3), g_rckwrk(1, nci3), ld
     *g_rckwrk, nrnu, ickwrk(icrnu), rckwrk(ncrnu), g_rckwrk(1, ncrnu), 
     *ldg_rckwrk, nord, ickwrk(icord), mxord, ickwrk(ickor), rckwrk(ncko
     *r), g_rckwrk(1, nckor), ldg_rckwrk)
C
        do 99747 k = 1, nkk
          do g_i_ = 1, g_p_
            g_wdot(g_i_, k) = 0.0d0
          enddo
          wdot(k) = 0.0d0
C--------
50        continue
99747   continue
        do 99745 n = 1, mxsp
          do 99746 i = 1, nii
            k = ickwrk(icnk + (i - 1) * mxsp + n - 1)
            if (k .ne. 0) then
              do g_i_ = 1, g_p_
                g_rop(g_i_) = -g_rckwrk(g_i_, nci2 + i - 1) + g_rckwrk(g
     *_i_, nci1 + i - 1)
              enddo
              rop = rckwrk(nci1 + i - 1) - rckwrk(nci2 + i - 1)
C--------
              d4_b = dble(ickwrk(icnu + (i - 1) * mxsp + n - 1))
              do g_i_ = 1, g_p_
                g_wdot(g_i_, k) = d4_b * g_rop(g_i_) + g_wdot(g_i_, k)
              enddo
              wdot(k) = wdot(k) + rop * dble(ickwrk(icnu + (i - 1) * mxs
     *p + n - 1))
C--------
C*****precision > double
C*****END precision > double
C*****precision > single
C     1         REAL (ICKWRK(IcNU + (I-1)*MXSP + N - 1))
C*****END precision > single
            endif
100         continue
99746     continue
99745   continue
C
        if (nrnu .le. 0) then
          return
        endif
        do 99743 l = 1, nrnu
          i = ickwrk(icrnu + l - 1)
          do g_i_ = 1, g_p_
            g_rop(g_i_) = -g_rckwrk(g_i_, nci2 + i - 1) + g_rckwrk(g_i_,
     * nci1 + i - 1)
          enddo
          rop = rckwrk(nci1 + i - 1) - rckwrk(nci2 + i - 1)
C--------
          do 99744 n = 1, mxsp
            k = ickwrk(icnk + (i - 1) * mxsp + n - 1)
            if (k .ne. 0) then
              do g_i_ = 1, g_p_
                g_wdot(g_i_, k) = rop * g_rckwrk(g_i_, ncrnu + (l - 1) *
     * mxsp + n - 1) + rckwrk(ncrnu + (l - 1) * mxsp + n - 1) * g_rop(g_
     *i_) + g_wdot(g_i_, k)
              enddo
              wdot(k) = wdot(k) + rop * rckwrk(ncrnu + (l - 1) * mxsp + 
     *n - 1)
C--------
            endif
200         continue
99744     continue
99743   continue
        return
      end subroutine g_ckwyp
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckwypk' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckwyr' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckxnum' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckxtcp' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckxtcr' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckxty' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      subroutine g_ckytcp(g_p_, p, t, g_t, ldg_t, y, g_y, ldg_y, ickwrk,
     * rckwrk, g_rckwrk, ldg_rckwrk, c, g_c, ldg_c)
C
C  START PROLOGUE
C
C  SUBROUTINE CKYTCP (P, T, Y, ICKWRK, RCKWRK, C)
C     Returns the molar concentrations given the pressure,
C     temperature and mass fractions;  see Eq. (7).
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
C     C      - Molar concentrations of the species.
C                   cgs units - mole/cm**3
C                   Data type - real array
C                   Dimension C(*) at least KK, the total number of
C                   species.
C
C  END PROLOGUE
C
C*****precision > double
        implicit double precision (a-h, o-z), integer (i-n)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
        dimension ickwrk(*), rckwrk(*), y(*), c(*)
C       Start of include file /home6/laniu/adifor/adiforMatab-Pre/ckstrt.h
C
C
C     include file for cklib.f
C
        common /ckstrt/ nmm, nkk, nii, mxsp, mxtb, mxtp, ncp, ncp1, ncp2
     *, ncp2t, npar, nlar, nfar, nlan, nfal, nrev, nthb, nrlt, nwl, neim
     *, njan, njar, nft1, nf1r, nexc, nrnu, nord, mxord, icmm, ickk, icn
     *c, icph, icch, icnt, icnu, icnk, icns, icnr, iclt, icrl, icrv, icw
     *l, icfl, icfo, ickf, ictb, ickn, ickt, icei, ictd, icjn, icf1, ice
     *x, icrnu, icord, ickor, ncaw, ncwt, nctt, ncaa, ncco, ncrv, nclt, 
     *ncrl, ncfl, nckt, ncwl, ncjn, ncf1, ncex, ncru, ncrc, ncpa, nckf, 
     *nckr, ncrnu, nckor, nck1, nck2, nck3, nck4, nci1, nci2, nci3, nci4
C
C     END include file for cklib.f
C
C

        integer g_i_, g_p_, ldg_y, ldg_rckwrk, ldg_t, ldg_c
        double precision d6_b, d6_v, d5_b, d4_b, d3_b, d4_v, d3_v, d5_v,
     * g_sumyow(g_p_), g_y(ldg_y, *)
        double precision g_rckwrk(ldg_rckwrk, *), g_t(ldg_t), g_c(ldg_c,
     * *)
        integer g_ehfid
!Laniu        save g_sumyow
        intrinsic dble
        data g_ehfid /0/
C
claniu        call ehsfid(g_ehfid, 'ckytcp','g_cklib.f')
C

        do g_i_ = 1, g_p_
          g_sumyow(g_i_) = 0.0d0
        enddo
        sumyow = 0.0d0
C--------
        do 99725 k = 1, nkk
          d4_v = y(k) / rckwrk(ncwt + k - 1)
          d4_b = 1.0d0 / rckwrk(ncwt + k - 1)
          d5_b = (-d4_v) / rckwrk(ncwt + k - 1)
          do g_i_ = 1, g_p_
            g_sumyow(g_i_) = d5_b * g_rckwrk(g_i_, ncwt + k - 1) + d4_b 
     ** g_y(g_i_, k) + g_sumyow(g_i_)
          enddo
          sumyow = sumyow + d4_v
C--------
150       continue
99725   continue
        d3_v = sumyow * t
        d4_b = rckwrk(ncru) * t
        d5_b = rckwrk(ncru) * sumyow
        do g_i_ = 1, g_p_
          g_sumyow(g_i_) = d3_v * g_rckwrk(g_i_, ncru) + d5_b * g_t(g_i_
     *) + d4_b * g_sumyow(g_i_)
        enddo
        sumyow = d3_v * rckwrk(ncru)
C--------
        do 99724 k = 1, nkk
          d5_v = sumyow * rckwrk(ncwt + k - 1)
          d6_v = p * y(k) / d5_v
          d3_b = (-d6_v) / d5_v
          d4_b = d3_b * rckwrk(ncwt + k - 1)
          d5_b = d3_b * sumyow
          d6_b = 1.0d0 / d5_v * p
          do g_i_ = 1, g_p_
            g_c(g_i_, k) = d5_b * g_rckwrk(g_i_, ncwt + k - 1) + d4_b * 
     *g_sumyow(g_i_) + d6_b * g_y(g_i_, k)
          enddo
          c(k) = d6_v
C--------
200       continue
99724   continue
        return
      end subroutine g_ckytcp
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckytcr' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ckytx' excluded from output
C
CInactive procedure or block data 'ifirch' excluded from output
C
CInactive procedure or block data 'ilasch' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ippari' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'ipparr' excluded from output
CInactive procedure or block data 'ipplen' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'pkcdyp' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'pkcpbs' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'pkcpms' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'pkcpor' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'pkctyp' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'pkcvbs' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'pkcvms' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'pkhml' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'pkhms' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'pkqyp' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'pkrhoc' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'pkrhox' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'pkrhoy' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'pkwyp' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'pkwypk' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'pkxtcp' excluded from output
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
CInactive procedure or block data 'pkytcp' excluded from output

c==================end of subroutine in g_cklib====================
