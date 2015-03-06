!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the x2f_mpi      !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module x2f_loclev

! This load-balancing strategy attempts to minimize message passing by
! pushing vectors to nearby processors in the sequence 1, 2, 3,... nproc
! until all have close to the same number of function evaluations to perform

use MPI
implicit none

contains  !-------------------------------------------------------------

subroutine x2f_loclev_init
return
end subroutine x2f_loclev_init


!-----------------------------------------------------------------------

subroutine x2f_loclev_p( nproc, mpicomm, myrank, ntogo, nv, p, &
               n_incoming, n_outgoing )

integer, intent(in)    :: nproc, mpicomm, myrank, ntogo, nv
integer, intent(inout) :: p(nv)
integer, intent(out)   :: n_incoming(:), n_outgoing(:)

integer :: k_me, ntogo_k, ntogoarray(nproc), ntogo_ave, nstaying
integer :: ntagged, n_pass(nproc), n_vac(nproc)
integer :: k, k_to, k_from, n_to_pass, i, ierr

n_incoming = 0
n_outgoing = 0

k_me = myrank + 1

! message pass to obtain ntogoarray(1:nproc) = [ ntogo on processes 1:nproc ]
! note, ntogo <= nv because some evaluations may have been done via quick try

call MPI_Allgather( ntogo, 1, MPI_INTEGER, ntogoarray, 1, MPI_INTEGER, &
                    mpicomm, ierr )

ntogo_ave = 1 + sum( ntogoarray(1:nproc) ) / nproc  ! note upward rounding

nstaying = min( ntogo, ntogo_ave )  ! process up to ntogo_ave vectors locally
n_outgoing(k_me) = nstaying
n_incoming(k_me) = nstaying

if( ntogo == ntogo_ave ) then  ! easy case: no vacancies, nothing to pass...
   where( p == -1 ) p = k_me   ! ...done in one step, nothing more to calculate!
   return
end if

! ...otherwise need to calculate arrays n_incoming and n_outgoing...
! first, set n_vac and n_pass for everyone (identical calculation on all procs)

do k = 1, nproc
   ntogo_k = ntogoarray(k)
   if( ntogo_k < ntogo_ave ) then
      n_vac(k) = ntogo_ave - ntogo_k   ! nvac is the number of "vacancies"
      n_pass(k) = 0                    ! no vectors need to be passed from here
   else 
      n_vac(k)  = 0                    ! no vacancies here
      n_pass(k) = ntogo_k - ntogo_ave  ! n_pass vectors need to be passed
   endif
end do

k_to = 1  ! processor 1 is first destination (maybe no vacancies there?)...
          ! for each k_from, the following loop identifies one or more k_to

do k_from = 1, nproc

   do while ( n_pass(k_from) > 0 ) 
      ! find one or more receivers for  n_pass(k_from)  vectors

      do while( n_vac(k_to) == 0 )  ! find next processor with vacancies
         k_to = k_to + 1
      end do

      n_to_pass = min( n_pass(k_from), n_vac(k_to)  )
      n_pass(k_from) = n_pass(k_from) - n_to_pass
      n_vac(k_to)    = n_vac(k_to)    - n_to_pass

      ! this is the only part of this loop that isn't identical everywhere
      if( k_to == k_me ) n_incoming(k_from) = n_to_pass
      if( k_from == k_me ) n_outgoing(k_to) = n_to_pass

   end do

end do

! now we're ready to make the processor assignments

i = 1
do k_to = 1, nproc
   n_to_pass = n_outgoing(k_to)
   ntagged = 0
   do while ( ntagged < n_to_pass )
      if( i > nv ) then
         write(0,*) 'x2f_loclev_p: could not make all processor assignments'
         return
      end if
      if( p(i) == -1 ) then
         p(i) = k_to
         ntagged = ntagged + 1
      end if
      i = i + 1
   end do
end do

return
end subroutine x2f_loclev_p

end module x2f_loclev
