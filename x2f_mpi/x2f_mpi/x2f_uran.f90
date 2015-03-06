!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the x2f_mpi      !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module x2f_uran

! Module for x2f_mpi to distribute the vectors uniformly and randomly.
! This guarantees load balancing but has a high degree of message passing.
! The current code was revised by SRL, 5/2011, as indicated.

use MPI
implicit none

! The following are now defined as module variables - SRL, 5/2011
integer, allocatable, save :: seednum(:), seed_on_entry(:)
integer, save :: seedsize = 0

contains  !-------------------------------------------------------------

subroutine x2f_uran_init( myrank )

! x2f_uran_init takes myrank as an argument to ensure that each process
! will generate a different seednum for subsequent calls to random_number

integer, intent(in)  :: myrank
!integer, allocatable :: seednum(:)  - see above for seednum and seedsize
integer :: i

if ( seedsize == 0 ) call random_seed( size = seedsize )
if ( .not.allocated(seednum) ) allocate( seednum(seedsize) )
seednum = (/ ( 123 * myrank + 7 * i, i = 1, seedsize ) /)
!seednum = (/ ( 123 * myrank + 6 * i, i = 1, seedsize ) /)
!call random_seed( put = seednum(1:seedsize) ) - moved by SRL, 5/2011
if ( .not.allocated(seed_on_entry) ) allocate( seed_on_entry(seedsize) )

return
end subroutine x2f_uran_init

!-----------------------------------------------------------------------

subroutine x2f_uran_p( nproc, mpicomm, myrank, ntogo, nv, p, &
               n_incoming, n_outgoing )

! (Comments modified throughout by SRL, 5/2001)
! The goal: reassign p(i) with a random number 1:nproc wherever it is initially -1.
! Only reassign p(i) if its initial value is -1; elements equal to zero are unchanged.
! Note, ntogo is the number of elements of p(1:nv) initially equal to -1.
! Elements of p equal 0 where evaluations have already been done via quick try.

integer, intent(in)    :: nproc, mpicomm, myrank, ntogo, nv
integer, intent(inout) :: p(nv)
integer, intent(out)   :: n_incoming(:), n_outgoing(:)

integer :: n_incoming1(nproc,nproc),ntogo_array(nproc)
integer :: msgsize, total_active, active_pool(nproc), pool_count(nproc)
integer :: j, i, j_pick, p_pick, ierr, n_left, n_offset, n_left_prfsum
real    :: roffset, x1
integer :: iip

! Note - we now make multiple calls to random_seed each time this subroutine is called.
! Why?  To keep other routines from altering the seed for this routine, and vice versa.
! Previously there was just a single initialization call in x2f_uran_init - SRL, 5/2011

call random_seed( get = seed_on_entry(1:seedsize) )
call random_seed( put = seednum(1:seedsize) )

! Start by figuring out how many particles will go to each of the other processes.
! n_outgoing(j) is the number of particles that are destined for process j.
! It is logical to initialize it ntogo/nproc (which could be zero), for all j.
! But there could be some particles left over...

msgsize = ntogo / nproc
n_outgoing(1:nproc) = msgsize
n_left = ntogo - msgsize * nproc

! Now figure out which processes get the leftovers.  Global communication is needed.
! To spread the leftovers evenly, single particles are assigned to sequential ranks.
! On rank 0, the starting point of the sequence is 1 plus a random offset < nproc.
! For higher ranks, the starting point is found through the "prefix sum" of n_left.
! This is the sum of leftovers on ranks < myrank.  We use MPI_Scan to compute this.
! It gives the rank (modulo nproc) of the last process to be assigned a particle.
! The random offset is necessary to ensure leftover particles do not keep returning
! to the same range of ranks each time this routine is called (added by SRL 5/2011).

if ( myrank == 0 ) then
   call random_number( roffset )
   n_offset = roffset * nproc  ! want a random offset between 0 and nproc-1
   n_left = n_left + n_offset
end if
call MPI_Scan( n_left, n_left_prfsum, 1, MPI_INTEGER, MPI_SUM, mpicomm, ierr)
if ( myrank == 0 ) then
   n_left = n_left - n_offset  ! remove offset after calculating prefix sum
end if
! We want the prefix sum to include only the lower ranks, excluding the current rank;
! on all ranks, the prefix sum will retain the random offset.
n_left_prfsum = n_left_prfsum - n_left
call iuranwor2( msgsize, nproc, n_left, n_left_prfsum, n_outgoing )

! At this point, all processes know n_outgoing, so an all-to-all determines n_incoming
call MPI_Alltoall( n_outgoing, 1, MPI_INTEGER, n_incoming, 1, MPI_INTEGER, mpicomm, ierr )

! All that's left to do is to assign the exact particles that are going to each process.
! Imagine drawing processor assignments one at a time out of a pool of slips of paper.
! The slips are numbered 1 to nproc.  A slip is replaced in the pool after it is drawn.
! But: after some number j has been drawn n_outgoing(j) times, we must remove its slip.
! pool_count(j) tells us how many more times a given number j is allowed to be drawn...

pool_count = n_outgoing

! pool_count(j) goes down by one each time j is drawn.  As drawing continues, it will
! eventually reach zero for some j.  When that occurs, total_active is reduced by 1,
! and active_pool is rearranged so that active_pool(k) /= j for k in 1:total_active.
! In the usual case, nv > ntogo > nproc.  We therefore start with total_active = nproc
! and active_pool(j) = j.  But we also want to handle situations where ntogo < nproc.
! That means some slips will not be present in the pool from the start:

total_active = 0
do j = 1, nproc
   if( n_outgoing(j) /= 0 ) then
      total_active = total_active + 1
      active_pool(total_active) = j
   end if
end do

do i = 1, nv
   if( p(i) == -1 ) then
      call random_number( x1 )
      ! since 0 <= x1 < 1, j_pick is a bounded integer, 1 <= j_pick <= total_active
      j_pick = int( x1 * total_active ) + 1
      p_pick = active_pool(j_pick)
      pool_count(p_pick) = pool_count(p_pick) - 1
      if( pool_count(p_pick) == 0 ) then
         if( j_pick /= total_active ) then
            ! move an active pool number from the last slot to the newly empty slot
            active_pool(j_pick) = active_pool(total_active)
         end if
         total_active = total_active - 1
      end if
      p(i) = p_pick
   end if
end do

call random_seed( get = seednum(1:seedsize) )
call random_seed( put = seed_on_entry(1:seedsize) )

return
end subroutine x2f_uran_p

subroutine iuranwor2( msgsize, nproc, nleft, nleftpsum, n_outgoing )

integer, intent(in)      :: msgsize, nproc, nleft, nleftpsum
integer, intent(inout)   :: n_outgoing(nproc)

integer   :: nleft_sub, nwraps

nwraps = nleftpsum/nproc
nleft_sub = nleftpsum - nwraps*nproc
if (nleft_sub+nleft <= nproc) then
   n_outgoing(nleft_sub+1:nleft_sub+nleft) = msgsize + 1
else
   n_outgoing(nleft_sub+1:nproc) = msgsize + 1
   n_outgoing(1:nleft-(nproc-nleft_sub)) = msgsize + 1
end if

return
end subroutine iuranwor2

end module x2f_uran
