!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the x2f_mpi      !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine x2f_mpi_nproc( nproc )

use MPI
implicit none
integer, intent(inout)  :: nproc
integer                 :: ierr 

call MPI_COMM_size( MPI_COMM_WORLD, nproc, ierr )

end subroutine x2f_mpi_nproc