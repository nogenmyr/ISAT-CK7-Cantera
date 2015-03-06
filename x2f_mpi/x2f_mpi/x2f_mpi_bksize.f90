!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the x2f_mpi      !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ISAT_MP_bksize(mpicomm,n_uneval,ldxf,blocksize,nb,nbmax,np)
!
! Decide the number of blocks for n_uneval particles
! Input:
!   mpicomm    - MPI communicator
!   n_uneval   - number of unevaluated particles
!   ldxf       - leading dimension of xf array
!   blocksize  - the size of array specified by user which is in terms of
!                Mega bytes
! Output:
!    nb        - number of blocks
!   nbmax      - maximum number of blocks
!    np        - average number of particles in each block
!                notice: the last block might have different number of particles
!
use MPI
implicit none
integer, intent(in)    :: mpicomm, n_uneval, ldxf
real, intent(in)       :: blocksize
integer, intent(out)   :: nb, nbmax, np

integer                :: ierr

np  = blocksize*1.048576e6/(kind(1.d0)*ldxf)
nb  = ceiling(n_uneval/float(np))
if ( np>n_uneval ) np = n_uneval
call MPI_Allreduce( nb, nbmax, 1, MPI_INTEGER, MPI_MAX, mpicomm, ierr )

end subroutine ISAT_MP_bksize