!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the x2f_mpi      !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

! Author: Varun Hiremath <vh63@cornell.edu>
! Date: Wed,  9 May 2012 14:36:25 -0500

! This file provides subroutines to read/write ISAT tables

!---------------------------------------------------------------
subroutine x2f_pt_index(index)
  ! returns the partition index for the current rank
  use x2f_mpi_grpdata
  implicit none

  integer, intent(out) :: index

  if( .not. impinit) then
     print *, 'x2f_init_table: ISAT_MP not yet initialized.'
     stop
  endif

  index = pt_index
  return
end subroutine x2f_pt_index
!---------------------------------------------------------------
subroutine x2f_save_table(itab_opp)
  ! save ISAT table 
  ! itab_opp = 0, save all the tables
  ! itab_opp = 1, save only one table (first) per partition

  use x2f_mpi_grpdata
  implicit none

  integer, intent(in) :: itab_opp

  integer :: info(100)
  real(kind(1.d0)) :: rinfo(50), stats(100)

  info     = -12345
  rinfo    = -12345.
  
  if(itab_opp == 0) then
     info(81) = 1 ! force writing the table
     call cisat( 12, info, rinfo, stats )
  elseif( mod(rank, set_pt) == 0 ) then
     ! simple method for now: 
     ! save ISAT table on the first rank of each partition
     info(81) = 1 ! force writing the table
     call cisat( 12, info, rinfo, stats )
     print *, 'x2f_save_table: saving table on rank ', rank
     print *, 'x2f_save_table: using partition index ', pt_index
  endif

  return
end subroutine x2f_save_table
!---------------------------------------------------------------
