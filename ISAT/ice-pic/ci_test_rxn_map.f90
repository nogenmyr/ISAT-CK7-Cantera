!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ci_test_rxn_map( press, dt )
  ! Computes the reaction mapping error under modes 7,8 and 9 using
  ! the exact mappings cc, ct directly evaluated under mode 6

  use ci_dat
  use ci_dat8
  use ci_ice_cksubs
  implicit none

  real(k_sp), intent(in) :: press, dt

  ! variables for modes 6 to 9
  real(k_sp) :: dpt(3)
  real(k_sp) :: cc6(ns+1), ct6(ns+1), cc7(ns+1), ct7(ns+1), cc89(nrc+1), ct89(nrc+1)
  real(k_dp) :: h, hs, hs_n, tmp, dcc6(ns+1), dcc89(nrc+1)

  ! others
  integer    :: lu_map, lu_ers, lu_out, lu_sts, imap, i
  real(k_dp) :: cpu_start, cpu_secs

  ! isat stats
  integer    :: info(100)
  real(k_sp) :: rinfo(50), stats(100)

  call isat_lu( lu_map )
  call isat_lu( lu_ers )
  call isat_lu( lu_out )
  call isat_lu( lu_sts )

  open( lu_map, file='rxn_map.op')

  call cpu_time(cpu_start)

  if( modeci == 7 ) then

     open( lu_ers, file='rxn_errs_7.op' )
     open( lu_out, file='rxn_mapi_7.op' )
     open( lu_sts, file='isat_stats_7.op' )

     do imap = 1, huge(1)
        read(lu_map, *, err=100, end=100) i, cc6, ct6, h, tmp, dpt
        cc7 = cc6

        dpt(2) = press
        call scirxn( dt, ns+1, cc7, ct7, dpt )
        write(lu_out, '(i8,1p1000e20.11)') imap, cc7, ct7

        call errors( imap, ns, cc6(1:ns), ct6(1:ns), cc7(1:ns), ct7(1:ns) )
     enddo

  elseif( modeci == 8 .or. modeci == 9 ) then

     if( modeci == 8) then
        open( lu_ers, file='rxn_errs_8.op' )
        open( lu_out, file='rxn_mapi_8.op' )
        open( lu_sts, file='isat_stats_8.op' )
     else
        open( lu_ers, file='rxn_errs_9.op' )
        open( lu_out, file='rxn_mapi_9.op' )
        open( lu_sts, file='isat_stats_9.op' )
     endif

     do imap = 1, huge(1)
        read(lu_map, *, err=100, end=100) i, cc6, ct6, h, tmp, dpt

        ! get hs
        hs = cc6(ns+1)
        dcc6 = cc6
        call hs2h( hs, href, dcc6(1:ns), ns, h )

        ! compute cc89 = {r, hs_n}
        cc89(1:nrc) = matmul(BBT, cc6(1:ns))
        dcc89 = cc89
        call h2hs_n( h, href_n, dcc89(1:nrc), nrc, hs_n )
        cc89(nrc+1) = hs_n

        dpt(2) = press
        call scirxn( dt, nrc+1, cc89, ct89, dpt )
        write(lu_out, '(i8,1p1000e20.11)') imap, cc89, ct89

        call errors( imap, nrs, cc6(CS), ct6(CS), cc89(1:nrs), ct89(1:nrs) )

     enddo

  endif

100 continue

  call cpu_time(cpu_secs)
  cpu_secs = cpu_secs - cpu_start

  info(81) = 2
  call scisat( 6, info, rinfo, stats )

  write(lu_sts,*)
  write(lu_sts,'(a,1p,e13.4)')'queries  = ', stats(1)
  write(lu_sts,'(a,1p,e13.4)')'p_retvs  = ', stats(2)
  write(lu_sts,'(a,1p,e13.4)')'s_retvs  = ', stats(3)
  write(lu_sts,'(a,1p,e13.4)')'grows    = ', stats(4)
  write(lu_sts,'(a,1p,e13.4)')'adds     = ', stats(5)
  write(lu_sts,'(a,1p,e13.4)')'replaces = ', stats(6)
  write(lu_sts,'(a,1p,e13.4)')'dir eval = ', stats(7)
  write(lu_sts,*)
  write(lu_sts,'(a,1p,e13.4)')'CPU (sec) =', cpu_secs
  write(lu_sts,'(a,1p,e13.4)')'CPU (mu s)/query = ', 1.e6*cpu_secs/stats(1)


contains

  subroutine errors( imap, size, cc, ct, cci, cti )

    use ci_utils
    implicit none
    
    integer, intent(in) :: imap, size
    real(k_sp), intent(in) :: cc(size), ct(size), cci(size), cti(size)

    ! local
    real(k_sp) :: errs(10)

    errs = 0.0
    errs(1) = snorm(ct-cc)                 ! dc =  |ct - cc|
    errs(2) = snorm(cti-cci)               ! dce = |cti - cci|
    errs(3) = dot_product(ct-cc, cti-cci)  ! angle
    errs(4) = snorm(cti - ct)              ! e = |cti - ct|

    write(lu_ers,'(i8,1p,20e15.6)') imap, errs(1:4)

  end subroutine errors

end subroutine ci_test_rxn_map
