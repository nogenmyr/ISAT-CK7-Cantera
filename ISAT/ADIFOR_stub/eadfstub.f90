      SUBROUTINE EADFSTUB(RNAME)
      USE isat_abort_m
      CHARACTER*6 RNAME
      CHARACTER*52 S1, S2, S3
      CHARACTER*162 SBIG
      S1 = ' is an ADIFOR stub library routine. You must link to'
      S2 = ' the real ADIFOR software to use this functionality.'
      S3 = ' Please refer to the README files for instructions. '
      WRITE (SBIG,'(A6,3A52)') RNAME, S1, S2, S3
      CALL isat_abort('ADIFOR STUB',0,mess=SBIG)
      RETURN
      END
