!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

!  This file contains dummy MPI routine for use with serial computations

subroutine isat_mpi_rank( myrank, nprocs )

!  return rank of this process and number of processes

   implicit none
   integer, intent(out) :: myrank, nprocs

   myrank = 0
   nprocs = 1

   return
  
end subroutine isat_mpi_rank  !-------------------

subroutine mpi_initialized( flag, ierr )

   implicit none
   logical, intent(out) :: flag
   integer, intent(out) :: ierr
   
   flag = .false.  !  MPI not initialized
   ierr = 0
   
   return
end subroutine mpi_initialized  !-----------------

subroutine mpi_abort( i, j, k )
!  This should not be called 
   implicit none
   integer :: i, j, k
   
   write(0,*)'mpi_abort called in serial run: stopping'
   stop

end subroutine mpi_abort