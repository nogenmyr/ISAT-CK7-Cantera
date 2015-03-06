!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ci_test_interface( dt, ncv, cc, ct, dpt )
  ! Test the interface subroutines

  use ci_dat,  only: modeci
  use ci_prec
  use ci_utils

  implicit none

  integer, intent(in)     :: ncv
  real(k_dp), intent(in)  :: dt, cc(ncv)
  real(k_dp), intent(out) :: ct(ncv)
  real(k_dp), intent(inout) :: dpt(3)

  integer  :: nsv, nrsv, nev, nfullv, krep, nfull_ext
  real(k_dp) :: mf
  real(k_dp), allocatable :: comp(:), ze(:), workspace(:)
  character(16), allocatable :: cname(:), enames(:)
  character(16) :: esnames(3)

  ! ciinit has already been called
  write(0,*) ' This is test run for modeci = ', modeci
  call print_var( cc, 'cc', 0 )

  ! calling cicomp
  krep = 3
  nfullv = ncv + 3
  allocate( cname(nfullv) )
  allocate( comp(nfullv) )

  comp(1) = dpt(2)
  call cicomp( 0, cc, krep, nfullv, comp, cname )
  call print_var( cname, 'cicomp (ncv = 0): cname', 0)

  comp(1) = dpt(2)
  call cicomp( ncv, cc, krep, nfullv, comp, cname )
  call print_var( comp, 'cicomp: comp', 0)

  ! call cicomp_sr for modeci > 7
  if( modeci > 7) then
       call cisize_full_ext( nfull_ext )
       write(0,*) 'cisize_full_ext = ', nfull_ext

       krep = 3
       nfullv = nfull_ext
       deallocate( cname )
       deallocate( comp )

       allocate( cname(nfullv) )
       allocate( comp(nfullv) )

       comp(1) = dpt(2)
       call cicomp_ext( 0, cc, krep, nfullv, comp, cname )
       call print_var( cname, 'cicomp_sr (ncv = 0): cname', 0)

       comp(1) = dpt(2)
       call cicomp_ext( ncv, cc, krep, nfullv, comp, cname )
       call print_var( comp, 'cicomp_ext: comp', 0)
  endif

  ! calling cirxn
  call cirxn( dt, ncv, cc, ct, dpt )
  call print_var( ct, 'cirxn: ct', 0 )
  
  ! call ci_dens_temp
  call ci_dens_temp( ncv, cc, dpt )
  call print_var( dpt, 'ci_dens_temp: dpt', 0 )

  ! call ci_mix_frac
  esnames(1) = 'CHO'
  esnames(2) = 'FUEL'
  esnames(3) = 'PILOT'
  allocate( workspace(2*ncv) )
  workspace(1) = 3
  workspace(2) = 0.4
  workspace(3) = 0.2
  workspace(4) = 0.4

  call ci_mix_frac( -1, esnames, ncv, cc, workspace, mf )
  call print_var( workspace, "workspace", 0 )
  call ci_mix_frac( 1, esnames, ncv, cc, workspace, mf )

  write(0,*) 'ci_mix_frac, mf = ', mf

  return
end subroutine ci_test_interface
