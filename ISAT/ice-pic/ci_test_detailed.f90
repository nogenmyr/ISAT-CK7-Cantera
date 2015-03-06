!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ci_test_detailed( ns_skel, press, dt )
  ! Computes the reaction mapping error for a skeletal mechanism
  ! relative to the detailed mechanism

  use ci_dat
  use ci_dat8
  use ci_ice_cksubs
  use ci_utils
  implicit none

  integer, intent(in) :: ns_skel
  real(k_sp), intent(in) :: press, dt

  ! variables for modes 6 to 9
  real(k_sp) :: dpt(3)
  real(k_sp) :: cc_skel(ns_skel+1), ct_skel(ns_skel+1)
  real(k_sp) :: ccr(ns), ctr(ns), ccu(ns-ns_skel), ctu(ns-ns_skel)
  real(k_sp) :: cc67(ns+1), ct67(ns+1)
  real(k_dp) :: tmp, h, hs, hs_n, dcc67(ns+1), y(ns)

  ! others
  integer    :: lu_sps, lu_map, lu_ers, lu_out, lu_sts
  integer    :: i, j, k, l, imap, nusp
  integer    :: skel_sps(ns_skel), sort_sps(ns_skel), unrep_sps(ns-ns_skel)
  real(k_dp) :: cpu_start, cpu_secs


  ! isat stats
  integer    :: info(100)
  real(k_sp) :: rinfo(50), stats(100)

  call isat_lu( lu_sps )
  call isat_lu( lu_map )
  call isat_lu( lu_ers )
  call isat_lu( lu_out )
  call isat_lu( lu_sts )

  open( lu_map, file='rxn_map.op')
  open( lu_sps, file='skel_sps.in')

  ! read indices of skeletal species in the detailed mechanism
  read(lu_sps, *) skel_sps

  sort_sps = skel_sps
  call sort( sort_sps )

  k = 1
  l = 1
  do j = 1, ns
     if( j /= sort_sps(l) ) then
        unrep_sps(k) = j
        k = k + 1
     else
        l = l + 1
     endif
  enddo

  nusp = ns - ns_skel

  call cpu_time(cpu_start)

  if( modeci == 6 .or. modeci == 7 ) then

     if( modeci == 6) then
        open( lu_ers, file='rxn_errs_skel_6.op' )
        open( lu_out, file='rxn_mapi_skel_6.op' )
        open( lu_sts, file='isat_stats_skel_6.op' )
     else
        open( lu_ers, file='rxn_errs_skel_7.op' )
        open( lu_out, file='rxn_mapi_skel_7.op' )
        open( lu_sts, file='isat_stats_skel_7.op' )
     endif

     do imap = 1, huge(1)
        read(lu_map, *, err=100, end=100) i, cc_skel, ct_skel, h, tmp, dpt

        ! set the skeletal species
        ccr = 0.0
        ccr(skel_sps) = cc_skel(1:ns_skel)

        ctr = 0.0
        ctr(skel_sps) = ct_skel(1:ns_skel)

        cc67(1:ns) = ccr(1:ns)
        dcc67(1:ns) = cc67(1:ns)
        call phi2y( dcc67(1:ns), amolwt, ns, y )
        call cthbms( tmp, y,   h, gas )
        call h2hs( h, href, dcc67(1:ns), ns, hs )
        cc67(ns+1) = hs

        call scirxn( dt, ns+1, cc67, ct67, dpt )
        write(lu_out, '(i8,1p1000e20.11)') imap, cc67, ct67

        ccu(1:nusp) = cc67(unrep_sps)
        ctu(1:nusp) = ct67(unrep_sps)

        call errors( imap, ns_skel, cc67(skel_sps), ct67(skel_sps), cc_skel(1:ns_skel), ct_skel(1:ns_skel) )
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

    ! -- additional errors pertaining to the skeletal mode
    errs(5) = snorm(cc67(1:ns))
    errs(6) = snorm(ccu)
    errs(7) = snorm(ct67(1:ns))
    errs(8) = snorm(ctu)

    write(lu_ers,'(i8,1p,20e15.6)') imap, errs(1:8)

  end subroutine errors

end subroutine ci_test_detailed
