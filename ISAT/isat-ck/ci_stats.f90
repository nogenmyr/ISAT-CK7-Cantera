!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module ci_stats

  ! Notes on usage of ci_stats.  

  ! To calculate the no. of calls and cputime spent on any operation,
  ! simply define a unique key for that operation in this module and
  ! then do the following:

  ! At the beginning of the operation:
  ! call routine_start(<key variable>)
  !    -> This updates the no. of calls
  !    -> Stores the starting time

  ! At the end of the operation:
  ! call routine_stop(<key variable>)
  !    -> This updates the cputime spent

  ! Alternately, you may also use update_calls(<key variable>) and
  ! update_time(<key variable>, cputime) to explicitly do the above
  ! things.

  use ci_prec
  use ci_dat,  only: modeci

  type :: operation
     integer :: calls
     integer :: failures
     real(k_dp) :: cputime   ! stored in micro seconds
     real(k_dp) :: starttime ! used to store starting time
  end type operation

  ! Initialise key variables to be used for storing stats
  integer, parameter :: i_ci_ceq                    = 1
  integer, parameter :: i_lp                        = 2
  integer, parameter :: i0_ci_ICE_chem_map          = 3
  integer, parameter :: i1_ci_ICE_chem_map          = 4
  integer, parameter :: i_ci_cem_tan                = 5
  integer, parameter :: i_ci_ice_tan                = 6
  integer, parameter :: i0_rmap2                    = 7
  integer, parameter :: i1_rmap2                    = 8
  integer, parameter :: i_ceq_nrs                   = 9
  integer, parameter :: i_ci_linrmap                = 12
  integer, parameter :: i_expm                      = 13
  integer, parameter :: i_ci_ice_facet_newton       = 14
  integer, parameter :: i_ci_ice_facet_newton_lmap  = 15
  integer, parameter :: i_ci_ice_pic_newton         = 16
  integer, parameter :: i_ci_ice_pic_newton_lmap    = 17
  integer, parameter :: i_ci_ice_recon_pic          = 18
  integer, parameter :: i_ci_ice_recon_pic_lmap     = 19
  integer, parameter :: i_ci_ice_recon_bc           = 20
  integer, parameter :: i_ci_ice_recon_bc_lmap      = 21
  integer, parameter :: i_ci_ice_recon              = 22
  integer, parameter :: i_ci_ice_pic_bound          = 23
  integer, parameter :: i_ci_dens_temp              = 24
  integer, parameter :: i_cidpt_dr                  = 25

  ! Set n_ops to maximum no. of operations.
  integer, parameter       :: n_ops = 30

  type(operation), save    :: op_stats(n_ops)

contains

  subroutine update_calls(index)
    ! Increments no. of calls for operation(index) by 1

    implicit none
    integer, intent(in) :: index

    op_stats(index)%calls = op_stats(index)%calls + 1

  end subroutine update_calls

  subroutine update_time(index, cputime)
    ! Increments time for operation(index) by time

    implicit none
    integer, intent(in)    :: index
    real(k_dp), intent(in) :: cputime

    op_stats(index)%cputime = op_stats(index)%cputime + 1.e6*cputime

  end subroutine update_time

  subroutine routine_start(index)
    ! Does two things:
    !     - Increments calls by 1
    !     - Stores starting time to be used by routine_stop

    implicit none
    integer, intent(in)    :: index

    call cpu_time(op_stats(index)%starttime)
    call update_calls(index)

  end subroutine routine_start

  subroutine routine_stop(index)
    ! Updates cputime for operation(index) by calling cpu_time to find
    ! end_time and using the stored starting time.

    implicit none
    integer, intent(in)    :: index
    real(k_dp)       :: stop_time

    call cpu_time(stop_time)
    call update_time(index,(stop_time - op_stats(index)%starttime))

  end subroutine routine_stop

  subroutine routine_failed(index)
    ! Updates no. of failures

    implicit none
    integer, intent(in)    :: index

    op_stats(index)%failures = op_stats(index)%failures + 1
  end subroutine routine_failed

  subroutine print_ci_stats(fracof, luout)
    ! Print stats, optionally as a fraction of operation(fracof)

    implicit none
    integer, intent(in), optional:: fracof
    integer, intent(in), optional:: luout
    integer :: lu

    if(present(luout)) lu = luout

    if(present(fracof)) then
       if(op_stats(fracof)%calls>0) then
          write(lu,'(a100)')'--------------------------------------------------------------------------------------------------'
          write(lu,'(a25,2a10,2a20,2a10)') 'Operation','Calls', 'Failures', 'Time (mu s)', 'Time/call (mu s)', 'Percent', '%/call'
          write(lu,'(a100)')'--------------------------------------------------------------------------------------------------'
          
          call frac_istats("ci_ceq", i_ci_ceq, fracof, lu)
          call frac_istats("ci_ice_chem_map (map=0)", i0_ci_ICE_chem_map, fracof, lu)
          call frac_istats("ci_ice_chem_map (map=1)", i1_ci_ICE_chem_map, fracof, lu)
          call frac_istats("ci_cem_tan", i_ci_cem_tan , fracof, lu)
          call frac_istats("ci_ice_tan", i_ci_ice_tan , fracof, lu)
          call frac_istats("ci_ice_pic_bound", i_ci_ice_pic_bound, fracof, lu)
          call frac_istats("lp", i_lp, fracof, lu)
          call frac_istats("rmap2 (map = 0)", i0_rmap2, fracof, lu)
          call frac_istats("rmap2 (map = 1)", i1_rmap2, fracof, lu)
          call frac_istats("ceq_nrs", i_ceq_nrs, fracof, lu)
          call frac_istats("ci_linrmap", i_ci_linrmap, fracof, lu)
          call frac_istats("expm", i_expm, fracof, lu)
          call frac_istats("ci_ice_facet_newton", i_ci_ice_facet_newton, fracof, lu)
          call frac_istats("ci_ice_facet_newton_lmap", i_ci_ice_facet_newton_lmap, fracof, lu)
          call frac_istats("ci_ice_pic_newton", i_ci_ice_pic_newton, fracof, lu)
          call frac_istats("ci_ice_pic_newton_lmap", i_ci_ice_pic_newton_lmap, fracof, lu)
          call frac_istats("ci_ice_recon_pic", i_ci_ice_recon_pic, fracof, lu)
          call frac_istats("ci_ice_recon_pic_lmap", i_ci_ice_recon_pic_lmap, fracof, lu)
          call frac_istats("ci_ice_recon_bc", i_ci_ice_recon_bc, fracof, lu)
          call frac_istats("ci_ice_recon_bc_lmap", i_ci_ice_recon_bc_lmap, fracof, lu)
          call frac_istats("ci_ice_recon", i_ci_ice_recon, fracof, lu)
          call frac_istats("ci_dens_temp", i_ci_dens_temp, fracof, lu)
          call frac_istats("cidpt_dr", i_cidpt_dr, fracof, lu)
          
          write(lu,'(a100)')'---------------------------------------------------------------------------------------------------'
       endif
    else
       write(lu,*)'--------------------------------------------------------------------------------'
       write(lu,'(a25,2a10,2a20)') 'Operation','Calls', 'Failures', 'Time (mu s)', 'Time/call (mu s)'
       write(lu,*)'--------------------------------------------------------------------------------'

       call istats("ci_ceq", i_ci_ceq, lu)
       call istats("ci_ice_chem_map (map=0)", i0_ci_ICE_chem_map, lu)
       call istats("ci_ice_chem_map (map=1)", i1_ci_ICE_chem_map, lu)
       call istats("ci_cem_tan", i_ci_cem_tan , lu)
       call istats("ci_ice_tan", i_ci_ice_tan , lu)
       call istats("ci_ice_pic_bound", i_ci_ice_pic_bound, lu)
       call istats("lp", i_lp, lu)
       call istats("rmap2 (map = 0)", i0_rmap2, lu)
       call istats("rmap2 (map = 1)", i1_rmap2, lu)
       call istats("ceq_nrs", i_ceq_nrs, lu)
       call istats("ci_linrmap", i_ci_linrmap, lu)
       call istats("expm", i_expm, lu)
       call istats("ci_ice_facet_newton", i_ci_ice_facet_newton, lu)
       call istats("ci_ice_facet_newton_lmap", i_ci_ice_facet_newton_lmap, lu)
       call istats("ci_ice_pic_newton", i_ci_ice_pic_newton, lu)
       call istats("ci_ice_pic_newton_lmap", i_ci_ice_pic_newton_lmap, lu)
       call istats("ci_ice_recon_pic", i_ci_ice_recon_pic, lu)
       call istats("ci_ice_recon_pic_lmap", i_ci_ice_recon_pic_lmap, lu)
       call istats("ci_ice_recon_bc", i_ci_ice_recon_bc, lu)
       call istats("ci_ice_recon_bc_lmap", i_ci_ice_recon_bc_lmap, lu)
       call istats("ci_ice_recon", i_ci_ice_recon, lu)
       call istats("ci_dens_temp", i_ci_dens_temp, lu)
       call istats("cidpt_dr", i_cidpt_dr, lu)

       write(lu,*)'------------------------------------------------------------------------------'
    endif
  end subroutine print_ci_stats

  subroutine istats(name, index, luout)
    ! helper function: called from print_ci_stats
    implicit none
    character(*), intent(in) :: name
    integer, intent(in) :: index
    integer, intent(in) :: luout

    if(op_stats(index)%calls>0) then
       write(luout,'(a25,2i10,1p,3e20.2)') name, op_stats(index)%calls, op_stats(index)%failures, &
            op_stats(index)%cputime, op_stats(index)%cputime/op_stats(index)%calls
    else
       write(luout,'(a25,2i10,1p,3e20.2)') name, op_stats(index)%calls, op_stats(index)%failures, &
            op_stats(index)%cputime, 0.d0
    endif

  end subroutine istats

  subroutine frac_istats(name, i, ir, luout)
    ! helper function: called from print_ci_stats
    ! prints stats of op_stats(i) as a fraction of op_stats(ir)
    implicit none
    character(*), intent(in) :: name
    integer, intent(in) :: i, ir
    integer, intent(in) :: luout

    if(op_stats(i)%calls>0) then
          write(luout,'(a25,2i10,1p,2e20.5,0p,2f10.2)') name, op_stats(i)%calls, op_stats(i)%failures, &
               op_stats(i)%cputime, op_stats(i)%cputime/op_stats(i)%calls, &
               op_stats(i)%cputime/op_stats(ir)%cputime*100, &
               op_stats(i)%cputime/op_stats(ir)%cputime*op_stats(ir)%calls/op_stats(i)%calls*100
    else
       write(luout,'(a25,2i10,1p,2e20.2,0p,2f10.2)') name, op_stats(i)%calls, op_stats(i)%failures, &
            op_stats(i)%cputime, 0.d0, 0.d0, 0.d0
    endif

  end subroutine frac_istats
  
end module ci_stats
