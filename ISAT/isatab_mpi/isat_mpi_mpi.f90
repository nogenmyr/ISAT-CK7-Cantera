!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine isat_mpi_rank( myrank, nprocs )

!  return rank of this process and number of processes

   use isat_abort_m
   use mpi
   implicit none
   integer, intent(out) :: myrank, nprocs

   logical :: flag
   integer :: ierr

!  check that MPI has been initialized

   call MPI_INITIALIZED( flag, ierr )
   if( .not.flag ) then
      write(lu_err,*)'isat_mpi_rank: MPI not initialized; flag, ierr = ', flag, ierr
	  stop
   endif

!  get myrank and nprocs

   call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr ) 
   call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

   return
end subroutine isat_mpi_rank