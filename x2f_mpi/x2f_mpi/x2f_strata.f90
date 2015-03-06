!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the x2f_mpi      !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module x2f_strata

! Module for x2f_mpi based on stratification of the input.
!
! "Stratification" here means:
! - project each vector x onto the direction determined by unit vector u
! - sort resulting scalars s into nbin bins, each having equal extent along u
!   (Note: by definition, bin boundaries are the same on every processor)
! - group bins into nproc batches so each batch has an equal number of x-vectors
!   (Note: assign x-vectors to bins locally; assign bins to batches globally)
! - if desired, assign batches to processors so as to maximize local processing

use MPI  ! MPI coding indicated by XXX in comments below

! Usage:
!   in calling routine, include the statements:
!     use x2f_strata
!     type (strata_type), save :: strat
!
!   call x2f_strata_init(...,strat) to initialize
!   call x2f_strata_p(strat,...) to assign x-vectors to processors

implicit none
type :: strata_type  !--------------------------------------------------

real, pointer    :: u(:)
integer, pointer :: bin(:), bp(:,:), a(:,:), bat2proc(:), bin2proc(:)
integer :: mpicomm, myrank, nproc, nbin, mode
logical :: set_smax, set_bat2proc
real    :: smin, smax

end type strata_type

real, allocatable, save, target    :: u_sv(:)
integer, allocatable, save, target :: bin_sv(:), bp_sv1(:,:), a_sv(:,:)
integer, allocatable, save, target :: bat2proc_sv(:), bin2proc_sv(:)


contains  !-------------------------------------------------------------

subroutine x2f_strata_init( nx, mbin, mpicomm, myrank, nproc, mode, uset, strat )

! Initialize data structure strat

! Input:
!   nx      - number of components of x-vector
!   mbin    - multiplier for nproc, making nbin = mbin * nprocs
!   mpicomm - MPI communicator
!   mode    - determines specification of vector u
!           = 0:  u = uset
!           = k, (1 <= k <= nx):  u(k) = 1, u(j) = 0, j/=k

integer, intent(in)              :: nx, mbin, mpicomm, myrank, nproc, mode
real, intent(in)                 :: uset(nx)
type(strata_type), intent(inout) :: strat

integer :: nbin, ierr, nproc0, nbin0

nbin   = mbin * nproc

strat%mpicomm = mpicomm
strat%myrank  = myrank
strat%nproc   = nproc
strat%nbin    = nbin
strat%mode    = mode

strat%set_smax     = .true.  ! smax and smin are to be set
strat%set_bat2proc = .true.  ! bat2proc is to be set

if( allocated( u_sv ) ) deallocate( u_sv )
if( allocated( bin_sv ) ) deallocate( bin_sv )
if( allocated( bin2proc_sv ) ) deallocate( bin2proc_sv )
allocate( u_sv(nx) )
allocate( bin_sv(nbin) )
allocate( bin2proc_sv(nbin) )
strat%u => u_sv
strat%bin => bin_sv
strat%bin2proc => bin2proc_sv

! (SRL) The arrays below only need their true, full size on rank 0 -
!       yet pointers must associated on all ranks, or an error may occur
if( myrank == 0 ) then
   nproc0 = nproc
   nbin0  = nbin
else
   nproc0 = 1
   nbin0  = 1
endif

if( allocated( bat2proc_sv ) ) deallocate( bat2proc_sv )
if( allocated( bp_sv1 ) ) deallocate( bp_sv1 )
if( allocated( a_sv ) ) deallocate( a_sv )
allocate( bat2proc_sv(nproc0) )
allocate( bp_sv1(nbin0,nproc0) )
allocate( a_sv(nproc0,nproc0) )
strat%bat2proc => bat2proc_sv
strat%bp => bp_sv1
strat%a => a_sv

if( mode == 0 ) then
   strat%u(1:nx) = uset(1:nx)
   if( abs(dot_product(uset(1:nx),uset(1:nx)) - 1.) > 1.e-5 ) then
      write(0,*)'x2f_strata_init: norm of uset is not equal to 1.0'
   end if
elseif( mode >= 1  .and.  mode <= nx ) then
   strat%u(1:nx) = 0.
   strat%u(mode) = 1.
else
   write(0,*)'x2f_strata_init: bad value of mode = ', mode
   stop
endif

return

end subroutine x2f_strata_init


!-------------------------------------------------------------

subroutine x2f_strata_p( strat, nx, ldxf, xf, ntogo, nv, p, &
               n_incoming, n_outgoing )

!  Assign x-vectors (hereafter, "particles") to processors

type (strata_type)           :: strat
integer, intent(in)          :: nx, ldxf, ntogo, nv
real(kind(1.d0)), intent(in) :: xf(ldxf,nv)
integer, intent(inout)       :: p(nv)
integer, intent(out)         :: n_incoming(:), n_outgoing(:)

integer :: i, ieval, islot(ntogo), jj, jp(ntogo), j, k, nbin, nproc, total, bnum
integer :: mpicomm, myrank, ierr, istat(MPI_STATUS_SIZE)
!real    :: s(ntogo), slocalmin, slocalmax, smin, smax, scale, sibounded
real(kind(1.d0)) :: s(ntogo), slocalmin, slocalmax, smin, smax, scale, sibounded
real    :: batch, sump, nmax

integer, dimension(strat%nproc) :: proc2bat

nbin    = strat%nbin
nproc   = strat%nproc
mpicomm = strat%mpicomm
myrank  = strat%myrank

!  evaluate s for all vectors that need it  -----------
i = 1
ieval = 0
do while ( ieval < ntogo )
   if ( i > nv ) then
      write(0,*) 'x2f_strata_p: found less than ntogo vectors to evaluate'
      return
   end if
   if( p(i) == -1 ) then
      ieval = ieval + 1
      s(ieval) = dot_product( strat%u(1:nx), xf(1:nx,i) ) 
      islot(ieval) = i
   end if
   i = i + 1
end do

!  evaluate smax and smin if need be  -----------
strat%set_smax = .true.
if( strat%set_smax ) then
   if ( ntogo>0 ) then
   slocalmin = minval( s(1:ntogo) )
   slocalmax = maxval( s(1:ntogo) )
   else
   slocalmin =  1.0d01
   slocalmax = -1.0d01 
   end if
   
   ! XXX message pass to obtain global max and min
   call MPI_Allreduce( slocalmin, smin, 1, MPI_REAL8, MPI_MIN, mpicomm, ierr ) 
   call MPI_Allreduce( slocalmax, smax, 1, MPI_REAL8, MPI_MAX, mpicomm, ierr )
   
   strat%smin = smin
   strat%smax = smax
   
   strat%set_smax = .false. 
   
! by default, evaluate smin and smax only once
! can be over-ridden by setting:  strat%set_smax = .true.  in calling routine
endif

!  assign particles to bins, and count particles in bins  ---------
strat%bin = 0
smax = strat%smax
smin = strat%smin
if( smin /= smax ) then
   !scale = nbin / ( smax - smin )
   scale = (1-1/nproc)*nbin/( smax - smin )
else
   scale = 1.0  ! number is arbitrary: s(i) = smin = smax for all i
end if

!nmax = ntogo*nproc/nbin
if (ntogo < 2*nproc ) then
   nmax = 1.0
else
   nmax = ntogo/(2*nproc)
end if
do i  = 1, ntogo
   sibounded = max( min( s(i) , 0.999999*smax ), smin )  ! result must be < smax
   jj = int( ( sibounded - smin ) * scale ) + 1
   do while (strat%bin(jj)>=nmax) 
      jj = jj+1
   end do
   if (jj>nbin) jj = nbin
   jp(i) = jj
   strat%bin(jj) = strat%bin(jj) + 1
end do

!  on node 0, assemble the array  bp  !----------------
!  bp(j,k) = number of particles in bin j that are resident on processor k
!  (definition of bp changes later)
!!!!!!!!!!!!!!!!!!!!!!!!!!
!! -- NOTATION ALERT -- !!  must distinguish (rank, 0:nproc-1) from (processor, 1:nproc)
!!!!!!!!!!!!!!!!!!!!!!!!!!
!  XXX rank 0 receives  strat%bin  from other processors
  
call MPI_Gather(strat%bin, nbin, MPI_INTEGER, strat%bp, nbin, MPI_INTEGER, 0, mpicomm, ierr)

if( myrank == 0 ) then
   
   total = sum( strat%bp(1:nbin,1:nproc) )  ! total number of particles
   batch = total / float(nproc)             ! number in each batch
   
!  redefine bp(j,k) to be the batch number of bin j on processor k  ------------
   sump = 0.
   bnum = 1
   strat%a(1:nproc,1:nproc) = 0  ! a(k,j) = number of particles in batch j on processor k

   do j = 1, nbin
      do k = 1, nproc
         strat%a(k,bnum) = strat%a(k,bnum) + strat%bp(j,k)
         sump            = sump + strat%bp(j,k)
         strat%bp(j,k)   = bnum
         if( sump > batch ) then
            sump = sump - batch
            bnum = min( bnum+1, nproc )
         endif
      end do
   end do

   strat%bat2proc = (/ (i, i = 1, nproc) /)

!  set up the inverse to bat2proc: composition with inverse is identity  --------------
!  ...these vectors are just shorthand ways of writing permutation matrices

   proc2bat(strat%bat2proc) = (/ (i, i = 1, nproc) /)

!  set bin2proc for each processor and pass back -----------------
!  set n_incoming so each processor knows how many total particles to expect from others
!  ...accordingly, for processor k, n_incoming = strat%a(:,proc2bat(k))

   do k = nproc, 1, -1
      do j = 1, nbin
         strat%bin2proc(j) = strat%bat2proc( strat%bp(j,k) )
      end do
      if ( k > 1 ) then  ! XXX pass to processor k [rank = k-1]
         call MPI_Send( strat%bin2proc, nbin, MPI_INTEGER, k-1, 101, mpicomm, ierr )
         call MPI_Send( strat%a(1,proc2bat(k)), nproc, MPI_INTEGER, k-1, 102, mpicomm, ierr )
      end if
   end do

!  assign n_incoming on processor 1 [rank = 0]; bin2proc is already set

   n_incoming(1:nproc) = strat%a(1:nproc,proc2bat(1))

else  ! XXX receive on rank > 0

   call MPI_Recv( strat%bin2proc, nbin, MPI_INTEGER, 0, 101, mpicomm, istat, ierr )
   call MPI_Recv( n_incoming, nproc, MPI_INTEGER, 0, 102, mpicomm, istat, ierr )

endif

!  set n_outgoing on all processors based on bin2proc

n_outgoing(1:nproc) = 0
do j = 1, nbin
   k = strat%bin2proc(j)
   n_outgoing(k) = n_outgoing(k) + strat%bin(j)
end do

!  assign particles to processors

do i = 1, ntogo
   p(islot(i)) = strat%bin2proc( jp(i) )
end do

return

end subroutine x2f_strata_p


!-------------------------------------------------------------

subroutine max_trace( n, a, p )

!  Given an n x n matrix a, determine the permutation of the columns
!  which maximizes the trace.  Column k is moved to p(k).

implicit none
integer, intent(in)  :: n, a(n,n)
integer, intent(out) :: p(n)

integer :: max_sweeps=1, j, jj, k, kk, nswap, sweep, gain

do j = 1, n  !  set permutation vector
   p(j)  = j
end do

do sweep = 1, max_sweeps  !  sweep over columns
   nswap = 0

   do jj = 1, n-1
      do kk = jj+1, n
         j = p(jj)
         k = p(kk)
         gain  = a(j,kk) + a(k,jj) - ( a(j,jj) +a (k,kk) )

         if( gain > 0 ) then
            p(jj) = k
            p(kk) = j
            nswap = nswap + 1
         endif
      end do
   end do

   if( nswap == 0 ) then
      return
   endif

end do

return
end subroutine max_trace

end module x2f_strata
