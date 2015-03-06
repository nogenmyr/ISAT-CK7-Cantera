!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the x2f_mpi      !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module x2f_mpi_adapt
! 
! Module for ISAT_MP to do adaptive strategy - plus related utility routines
!
! "adaptation" here means:
!    find out the optimum number of processors and group the processors into 
!    different groups. The adaptation strategy is based on the prediction of 
!    wall clock time for the future stage.
!
use MPI
implicit none
integer              :: a_nolog
real, allocatable    :: t_Ri_ary(:), t_Di_ary(:), P_Di_ary(:), P_ij_mtx(:,:), &
                        a_i_ary(:), tRij_mtx(:,:), P_Di_PR1(:)
integer, allocatable :: N_i_ary(:), match(:)
integer, allocatable, save  :: Nij_ary(:)
!
! Statistics used to estimate the wall clock time in the future stage
!   t_Ri_ary(i)  - average retrieve time for group i 
!   t_Di_ary(i)  - average direct evaluation time for group i
!   P_Di_ary(i)  - percentage of doing direct evaluation for group i
!    P_ij(i,j)   - probability of a particle in group i being retrievable
!                  from group j
!    P_ij_mtx    - matrix of P_ij
!    a_i_ary(i)  - number of table entries for group i
!    N_i_ary(i)  - average number of particles in i-th group
!    match(i)    = j (j>0.and.j/=0) group i should match with group j   
!                = 0 no group pairing performed
!                
contains
   
   subroutine x2f_mpi_ptinit( n_base, pt_index, mpicomm_pt, lu_op )
   !
   ! Divide the processors used in the simulation into several partitions, 
   ! return the partition index for each processor and the communicator used
   ! within each partition
   !
   ! Input:
   !    n_base  - the desired partition size (don't partition if == 0)
   
   ! Output:
   !    pt_index  - indicate which partition the processor belongs to
   !   mpicomm_pt - indicate the communicator used within partition
   !
   use x2f_mpi_grpdata, only : set_pt
   implicit none

   integer, intent(in)   :: n_base, lu_op
   integer, intent(out)  :: pt_index, mpicomm_pt
   integer  :: nproc, iproc, ierr, n_pt, i_pt

   call MPI_COMM_RANK( MPI_COMM_WORLD, iproc, ierr )
   call MPI_COMM_SIZE( MPI_COMM_WORLD, nproc, ierr )

   if ( nproc<n_base .or. n_base==0 ) then
      pt_index   = 1
      set_pt     = nproc ! partition size = nproc
      mpicomm_pt = MPI_COMM_WORLD
   else
      if ( mod(nproc,n_base)==0 ) then
         n_pt = nproc/n_base
         do i_pt = 1, n_pt
            if ( iproc>=(i_pt-1)*n_base.and.iproc<=i_pt*n_base-1 ) then
               pt_index = i_pt
               exit
            end if
         end do
      else
         print *, 'the number of processors can not be partitioned!'
         stop
      end if
   end if
   
   call MPI_COMM_SPLIT( MPI_COMM_WORLD, pt_index, iproc, mpicomm_pt, ierr )
   if ( n_base > 0 .and. n_base < nproc .and. iproc == 0 .and. ierr == MPI_SUCCESS) then
      write(lu_op,'(A,I5,A,I5,A)') 'x2f_mpi_ptinit:', nproc, ' procs partitioned into blocks of', n_base, ' procs.'
   endif

   end subroutine x2f_mpi_ptinit

   subroutine iflocalstg( localstg, mpicomm_pt, tstg1, lu_op )
   !
   ! Test whether to stop using the PLP mode and switch to a different one
   ! Added for use with new mode-switching strategies, SRL 2/2011
   !
   implicit none
   integer, intent(in)    :: mpicomm_pt
   integer, intent(in)    :: tstg1, lu_op
   logical, intent(inout) :: localstg
   double precision, save :: wtime0 = -1
   double precision       :: wtimer
   real    :: a_c(2), a_star, a_cmax(2)
   integer :: ierr, iproc

   localstg = .true.
   
   call ci_stats( 12, a_c(1) )
   call ci_stats( 13, a_star )
   a_c(1) = a_c(1)/a_star

   if ( wtime0 < 0 ) wtime0 = MPI_WTIME()
   wtimer = MPI_WTIME() - wtime0
   ! Convert elapsed wall time from seconds to hours
   a_c(2) = wtimer/3600.
   
   call MPI_ALLREDUCE(a_c, a_cmax, 2, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr)
   
   ! Check if a table is 25% full or wall clock time exceeds specified time tstg1 
   if ((a_cmax(1) > 0.25) .or. (a_cmax(2) > tstg1/60.)) localstg = .false.

   ! Log data
   call MPI_COMM_RANK( MPI_COMM_WORLD, iproc, ierr )
   if ( iproc == 0 .and. ierr == MPI_SUCCESS ) then
      write(lu_op, '(A,10F8.2)') 'iflocalstg: table(%), time', a_cmax(1)*100, a_cmax(2)
   endif

   end subroutine iflocalstg

   subroutine x2f_mpi_grpinit( mpicomm_pt, grp_index, g, ng )
   !
   ! Return the group index(grp_index) for each processor, number of processors
   ! (g) within each group and number of group(ng) initially. Notice some 
   ! processors might have same group index, but they belong to different 
   ! partition. 
   !
   ! Input:
   !   mpicomm_pt  - communicator used within group
   ! Output:
   !   grp_index   - group index for each processor
   !      g        - number of processor within each group, g=1 initially
   !     ng        - number of groups

   implicit none
   integer, intent(in)    :: mpicomm_pt
   integer, intent(inout) :: grp_index, g, ng
   integer                :: ierr

   call MPI_COMM_SIZE( mpicomm_pt, ng, ierr )
   call MPI_COMM_RANK( mpicomm_pt, grp_index, ierr )
   grp_index = grp_index + 1
   g         = 1

   end subroutine x2f_mpi_grpinit 

   subroutine ci_stats( index, stat_index )
   !
   ! Return ISAT statistics: stat_index = stats(index)
   !
   ! Input:
   !    index      - the index of statistics requested
   ! Output:
   !   stat_index  - ISAT statistics, stats(index)
   !
   implicit none
   integer, intent(in)    :: index
   real, intent(inout) :: stat_index

   integer                :: info(100)
   real(kind(1.d0))       :: rinfo(100)
   real(kind(1.d0))       :: stats(100)

   info     = -12345
   rinfo    = -12345.
   info(81) = 0
   call cisat( 6, info, rinfo, stats )
   stat_index = real(stats(index))

   end subroutine ci_stats
   
   subroutine x2f_mpi_setTabEntry( msstg, ips, mps, A_max_entry, a_i_limit, a_0 )
   !
   ! Set the number of table entries per processor on i-th pairing stage limit  
   ! to a_i_limit
   !
   ! Input:
   !     msstg     - if previous stage is measuring stage or not
   !      ips      - the i-th pairing stage
   !      mps      - total number of pairing stage
   !      A_max    - maximum number of table entries
   !      
   ! Output:
   !   a_i_limit   - the number of table entries per processor on i-th pairing
   !                 stage

   implicit none
   logical, intent(in)    :: msstg
   integer, intent(in)    :: ips, mps
   real, intent(in)       :: A_max_entry
   real, intent(inout)    :: a_i_limit, a_0
  
   integer  :: info(100), test
   real     :: rinfo(100), stats(100), q

   test = mps
      if ( ips<test ) then
         q = 1./2.
         a_i_limit = -(A_max_entry*q*(1-q**ips)/(1-q))
         
      elseif ( ips==test ) then
         a_i_limit  = -A_max_entry
      end if

      info     = -12345
      rinfo    = -12345.
      info(80) = 0
      rinfo(8) = real(a_i_limit) 

      call scisat( 1, info, rinfo, stats )
      
      call ci_stats( 12, a_0 )

   end subroutine x2f_mpi_setTabEntry

   subroutine pairstg_init( mpicomm_pt, ips, mps )
   !
   ! Initialize ips, and return the total number of pairing stage needed for 
   ! the simulation
   !
   ! Input:
   !   mpicomm_pt  - communication used within partition
   ! Output:
   !      ips      - the index of current pairing stage
   !      mps      - total number of pairing stages
   !
   implicit none
   integer, intent(in)     :: mpicomm_pt
   integer, intent(inout)  :: ips, mps
   
   integer                 :: nproc, ierr

   ips = 1
   call MPI_COMM_SIZE( mpicomm_pt, nproc, ierr )
   mps = int(log(real(nproc))/log(2.))

   end subroutine pairstg_init

   subroutine ifpairstg( ips, mps, msstg, pairstg, mpicomm_pt, lu_log )
   !
   ! Check to see if this is pairing stage or not
   ! Input:
   !      mps      - total number of pairing stages
   !     msstg     - if last step is measuring stage or not
   !   mpicomm_pt  - communicator within partition
   ! Output:
   !    pairstg    - pairing stage or not
   !      ips      - the i-th pairing stage

   implicit none
   integer, intent(in)    :: mps, mpicomm_pt, lu_log
   integer, intent(inout) :: ips
   logical, intent(in)    :: msstg
   logical, intent(inout) :: pairstg
   
   integer  :: test
   
   test = mps
      
   if ( ips<=test ) then
      pairstg = .true.
      if ( msstg ) then
         ips = ips + 1
         if ( ips>test ) pairstg = .false.
      end if
   else
      pairstg = .false.
   end if

   end subroutine ifpairstg
   
   subroutine x2f_mpi_grpgather( g, intracomm, ng, mpicomm_pt, wall_R, wall_D, a_0, nv_un, nv_un_PR1, nv, lu_log )
   !
   ! Gather some statistics information from ng groups: t_Di_ary(ng), t_Ri_ary(ng)
   ! P_Di_ary(ng), a_i_ary(ng)
   !
   ! Input:
   !wall_R  - wall clock time spent in retrieve attempt
   !wall_D  - wall clock time spent in direct evaluation attempt
   !  nv_un - number of unevaluated particles
   !   nv   - total number of particles
   !    g   - number of processors in each group
   !   ng   - number of groups
   ! intracomm  - intracommunicator 
   ! mpicomm_pt - communicator within partition
   !
   ! Output:
   ! t_Di_ary, t_Ri_ary, P_Di_ary and a_i_ary
   !
   implicit none
   integer, intent(in)    :: g, intracomm, ng, mpicomm_pt, nv_un, nv_un_PR1, nv, lu_log
   real, intent(in)       :: wall_R, wall_D, a_0
   
   real    :: t_D, t_R, P_D, a_i, P_D_PR1
   integer :: ierr, rank_in_grp, N_i, newcomm, rank
   
   if(a_nolog == 0) write(lu_log,'(A,10I10)') 'grpgather:', nv_un, nv_un_PR1, nv
   call x2f_mpi_measure( intracomm, g, wall_R, wall_D, a_0, nv, nv_un, nv_un_PR1, t_R, t_D, P_D, a_i, N_i, P_D_PR1, lu_log )
   call MPI_COMM_RANK( intracomm, rank_in_grp, ierr )
   
   call x2f_mpi_crtintra( rank_in_grp, mpicomm_pt, newcomm )
   if ( rank_in_grp==0 ) then
      call MPI_GATHER( t_D, 1, MPI_REAL, t_Di_ary, 1, MPI_REAL, 0, newcomm, ierr )
	  call MPI_GATHER( t_R, 1, MPI_REAL, t_Ri_ary, 1, MPI_REAL, 0, newcomm, ierr )
	  call MPI_GATHER( P_D, 1, MPI_REAL, P_Di_ary, 1, MPI_REAL, 0, newcomm, ierr )
	  call MPI_GATHER( a_i, 1, MPI_REAL, a_i_ary, 1, MPI_REAL, 0, newcomm, ierr )
	  call MPI_GATHER( N_i, 1, MPI_INTEGER, N_i_ary, 1, MPI_INTEGER, 0, newcomm, ierr )
	  call MPI_GATHER( P_D_PR1, 1, MPI_REAL, P_Di_PR1, 1, MPI_REAL, 0, newcomm, ierr )
	  	                   
	  call MPI_COMM_RANK( newcomm, rank, ierr )
	  if ( rank==0 .and. a_nolog == 0 ) then
	     write(lu_log,'(A,10000e15.6)') 't_Di_ary', t_Di_ary
	     write(lu_log,'(A,10000e15.6)') 't_Ri_ary', t_Ri_ary
	     write(lu_log,'(A,10000e15.6)') 'P_Di_ary', P_Di_ary
	     write(lu_log,'(A,10000e15.6)') 'P_Di_PR1', P_Di_PR1
	     write(lu_log,'(A,10000e15.6)') 'a_i_ary ', a_i_ary
	     write(lu_log,'(A,10000I8)')    'N_i_ary ', N_i_ary
	  end if  
   end if
   
   call x2f_mpi_freecomm( newcomm )
   
   end subroutine x2f_mpi_grpgather
   
   subroutine x2f_mpi_measure( intracomm, g, wall_R, wall_D, a_0_local, nv, nv_un, nv_un_PR1, t_R, t_D, &
                               P_D, a_i, N_i, P_D_PR1, lu_log )
   !
   ! Measure average statistics of t_R, t_D, P_D and a_i(sum of entries in group i) for group i
   !
   implicit none
   integer, intent(in)    :: intracomm, nv, nv_un, nv_un_PR1, lu_log, g
   real, intent(in)       :: wall_R, wall_D, a_0_local
   real, intent(inout)    :: t_R, t_D, P_D, a_i, P_D_PR1
   integer, intent(inout) :: N_i
   
   real     :: wall_R_max, wall_D_max, n_a_sum, n_d_sum, a_i_local, n_q, n_r, n_g, n_a, n_d 
   real     :: n_f, cpu_f, t_D_local, a_0
   integer  :: nv_ave, nv_un_ave, ierr, nv_un_PR1_ave
   
   call ci_stats( 77, n_f )
   call ci_stats( 91, cpu_f )
   t_D_local = cpu_f/n_f
   
   call isat_stats( 2, n_q, n_r, n_g, n_a, n_d, lu_log )
   call ci_stats( 12, a_i_local )
   call MPI_REDUCE( n_a, n_a_sum, 1, MPI_REAL, MPI_SUM, 0, intracomm, ierr )
   call MPI_REDUCE( n_d, n_d_sum, 1, MPI_REAL, MPI_SUM, 0, intracomm, ierr )
   call MPI_REDUCE( a_i_local, a_i, 1, MPI_REAL, MPI_SUM, 0, intracomm, ierr )
   call MPI_REDUCE( a_0_local, a_0, 1, MPI_REAL, MPI_SUM, 0, intracomm, ierr )
   call MPI_REDUCE( wall_R, wall_R_max, 1, MPI_REAL, MPI_MAX, 0, intracomm, ierr )
   call MPI_REDUCE( nv, nv_ave, 1, MPI_INTEGER, MPI_SUM, 0, intracomm, ierr )
   call MPI_REDUCE( nv_un_PR1, nv_un_PR1_ave, 1, MPI_INTEGER, MPI_SUM, 0, intracomm, ierr )
   call MPI_REDUCE( t_D_local, t_D, 1, MPI_REAL, MPI_MAX, 0, intracomm, ierr )
   t_R = wall_R_max*g/nv_ave
   !P_D = (n_a_sum+n_d_sum*n_d_sum/nv_ave)/nv_ave
   P_D = (n_a_sum+n_d_sum/2.0)/nv_ave
   P_D_PR1 = nv_un_PR1_ave/real(nv_ave)
   !write(lu_log,'(A,10E15.6)') 'check1', n_a_sum, n_d_sum, a_0, a_i, p_D
   !if ( a_0_local>1.0 ) p_D = p_D*a_0/a_i
   !write(lu_log,'(A,10E15.6)') 'check2', n_a_sum, n_d_sum, a_0, a_i, p_D
   N_i = int(nv_ave*1.0/g)
   
   end subroutine x2f_mpi_measure
   
   subroutine isat_stats( inout, n_q, n_r, n_g, n_a, n_d, lu_log )
   
   implicit none
   integer, intent(in)  :: inout, lu_log
   real, intent(inout)  :: n_q, n_r, n_g, n_a, n_d
   real, save  :: n_q_pre, n_r_pre, n_g_pre, n_a_pre, n_d_pre

   real        :: n_pr, n_sr
   
   call ci_stats(1, n_q)
   call ci_stats(2, n_pr)
   call ci_stats(3, n_sr)
   call ci_stats(4, n_g)
   call ci_stats(5, n_a)
   call ci_stats(7, n_d)
   n_r = n_pr + n_sr
   
   if ( inout==1 ) then
      n_q_pre = n_q
      n_r_pre = n_r
      n_g_pre = n_g
      n_a_pre = n_a
      n_d_pre = n_d
   elseif ( inout==2 ) then
      n_q = n_q - n_q_pre
      n_r = n_r - n_r_pre
      n_g = n_g - n_g_pre
      n_a = n_a - n_a_pre
      n_d = n_d - n_d_pre
   end if   
   !write(lu_log,'(A,I4,10E15.6)') 'inout', inout, n_q, n_r, n_a, n_g
   
   end subroutine isat_stats

   subroutine x2f_mpi_allocate_adapt( size )
   !
   ! Allocate arrays for later use in adaptive strategy
   !
   implicit none
   integer, intent(in)  :: size
   
   if ( allocated(t_Ri_ary)  )  deallocate( t_Ri_ary )
   if ( allocated(t_Di_ary)  )  deallocate( t_Di_ary )
   if ( allocated(P_Di_ary)  )  deallocate( P_Di_ary )
   if ( allocated(a_i_ary )  )  deallocate( a_i_ary  )
   if ( allocated(P_ij_mtx)  )  deallocate( P_ij_mtx )
   if ( allocated(tRij_mtx)  )  deallocate( tRij_mtx )
   if ( allocated(N_i_ary)   )  deallocate( N_i_ary  )
   if ( allocated( match  )  )  deallocate(  match   )
   if ( allocated(P_Di_PR1)  )  deallocate( P_Di_PR1 )
   allocate(    t_Ri_ary(size)    )
   allocate(    t_Di_ary(size)    )
   allocate(    P_Di_ary(size)    )
   allocate(    a_i_ary(size)     )
   allocate(  P_ij_mtx(size,size) )
   allocate(  tRij_mtx(size,size) )
   allocate(    N_i_ary(size)     )
   allocate(      match(size)     )
   allocate(    P_Di_PR1(size)    )
   t_Ri_ary = 0.0
   t_Di_ary = 0.0
   P_Di_ary = 0.0
   a_i_ary  = 0.0
   P_ij_mtx = 0.0
   tRij_mtx = 0.0
   N_i_ary  = 0
   match    = 0
   P_Di_PR1 = 0.0
   
   end subroutine x2f_mpi_allocate_adapt
   
   subroutine x2f_mpi_Pij( nv_un, nv_nolocal, mpicomm_pt, intracomm, i )
   !
   ! Measure the quantity Pij(grp_index, grp_h), the probability that a particle 
   ! from group grp_index can be retrieved (without DE) from the ISAT table from
   ! group grp_h 
   !
   implicit none
   integer, intent(in)  :: nv_un, nv_nolocal, mpicomm_pt, intracomm, i

   integer  :: nv_un_sum, nv_nolocal_sum, ierr, rank_in_grp, newcomm, rank
   real     :: P_ij

   P_ij = 0
   call MPI_REDUCE( nv_un, nv_un_sum, 1, MPI_INTEGER, MPI_SUM, 0, intracomm, &
                    ierr )
   call MPI_REDUCE( nv_nolocal, nv_nolocal_sum, 1, MPI_INTEGER, MPI_SUM, 0, &
                    intracomm, ierr )
   call MPI_COMM_RANK( intracomm, rank_in_grp, ierr )                                
   if (rank_in_grp==0) then
      P_ij = 1-real(nv_un_sum)/real(nv_nolocal_sum)
   end if
   
   call x2f_mpi_crtintra( rank_in_grp, mpicomm_pt, newcomm )
   
   call MPI_ALLGATHER( P_ij, 1, MPI_REAL, P_ij_mtx(i,:), 1, MPI_REAL, newcomm, &
                       ierr )
   
   call x2f_mpi_freecomm( newcomm )
   
   end subroutine x2f_mpi_Pij

   subroutine x2f_mpi_Pij_mtx( ng, lu_log )
   !
   ! convert Pij matrix to the pattern of how pairing is formed
   ! 
   implicit none
   
   integer, intent(in)  :: ng, lu_log
   integer  :: grp_mtx(ng, ng), i
   real     :: temp(ng)

   call x2f_mpi_crtgrpmtx( ng, grp_mtx )
   call x2f_mpi_cvtPij( ng, grp_mtx, P_ij_mtx )
      
   if(a_nolog == 0) then
      write(lu_log,'(A)') 'P_ij_mtx'
      do i = 1, ng
         temp = P_ij_mtx(i,:)
         write(lu_log,'(10000E15.6)') temp
      end do
   endif

   end subroutine x2f_mpi_Pij_mtx
   
   subroutine x2f_mpi_tRij( t_i, nv_nolocal, mpicomm_pt, intracomm, i, g )
   !
   ! Measure the quantity Pij(grp_index, grp_h), the probability that a particle 
   ! from group grp_index can be retrieved (without DE) from the ISAT table from
   ! group grp_h 
   !
   implicit none
   real, intent(in)     :: t_i
   integer, intent(in)  :: nv_nolocal, mpicomm_pt, intracomm, i, g

   real     :: t_i_max, tRij
   integer  :: ierr, rank_in_grp, newcomm, nv_nolocal_sum

   tRij = 0.0
   call MPI_REDUCE( nv_nolocal, nv_nolocal_sum, 1, MPI_INTEGER, MPI_SUM, 0, &
                    intracomm, ierr )
   call MPI_REDUCE( t_i, t_i_max, 1, MPI_REAL, MPI_MAX, 0, intracomm, ierr )
   call MPI_COMM_RANK( intracomm, rank_in_grp, ierr )                                
   if (rank_in_grp==0) then
      tRij = t_i_max/(nv_nolocal_sum/g)
   end if
   
   call x2f_mpi_crtintra( rank_in_grp, mpicomm_pt, newcomm )
   
   call MPI_ALLGATHER( tRij, 1, MPI_REAL, tRij_mtx(i,:), 1, MPI_REAL, newcomm, &
                       ierr )
   
   call x2f_mpi_freecomm( newcomm )
   
   end subroutine x2f_mpi_tRij

   subroutine x2f_mpi_tRij_mtx( ng, lu_log )
   !
   ! convert tRij matrix to the pattern of how pairing is formed
   ! 
   implicit none
   
   integer, intent(in)  :: ng, lu_log
   integer  :: grp_mtx(ng, ng), i
   real     :: temp(ng)

   call x2f_mpi_crtgrpmtx( ng, grp_mtx )
   call x2f_mpi_cvtPij( ng, grp_mtx, tRij_mtx )
      
   if(a_nolog == 0) then
      write(lu_log,'(A)') 'tRij_mtx'
      do i = 1, ng
         temp = tRij_mtx(i,:)
         write(lu_log,'(10000E15.6)') temp
      end do
   endif

   end subroutine x2f_mpi_tRij_mtx
   
   subroutine x2f_mpi_crtgrpmtx( ng, grp_mtx )
   !
   ! Based on the pattens how each pair form, create the corresponding matrix 
   !
   integer, intent(in)    :: ng
   integer, intent(inout) :: grp_mtx(ng,ng)
   
   integer  :: i, j, k, grp_ary(ng)
   do i = 1, ng
      grp_ary(i) = i
   end do
   grp_mtx(1,1:ng) = grp_ary(1:ng)
   
   do i = 2, ng
      do j = 1, ng
         call x2f_mpi_pair( i, ng, j, k )
         grp_ary(j) = k
      end do
      grp_mtx(i,1:ng) = grp_ary(1:ng)
   end do
      
   end subroutine x2f_mpi_crtgrpmtx
   
   subroutine x2f_mpi_cvtPij( ng, grp_mtx, P_ij_mtx )
   !
   ! Convert collected Pij matrix to its actual form 
   !
   implicit none
   integer, intent(in)    :: ng
   integer, intent(inout) :: grp_mtx(ng,ng) 
   real, intent(inout)    :: P_ij_mtx(ng,ng)
   
   integer :: i, j, rcd1, rcd2
   real    :: rcd3
   do i = 1, ng
      do j = 1, ng
         rcd1 = i
         rcd2 = grp_mtx(i,j)
         if ( rcd2.ne.i ) then
            rcd3 = P_ij_mtx(i,j)
            call swap( rcd1, rcd2, rcd3, j, ng, grp_mtx, P_ij_mtx )
         end if
      end do
   end do   
   
   end subroutine x2f_mpi_cvtPij
   
   recursive subroutine swap(rcd1, rcd2, rcd3, j, N, grp_mtx, rdmtx)

   implicit none
   integer, intent(in)     :: rcd1, j, N
   integer, intent(inout)  :: rcd2, grp_mtx(N,N)
   real, intent(inout)     :: rcd3, rdmtx(N,N)

   integer  :: temp
   real     :: temprd

   temp   = grp_mtx(rcd2,j)
   temprd = rdmtx(rcd2,j)
   grp_mtx(rcd2,j) = rcd2
   rdmtx(rcd2,j)   = rcd3
   rcd2 = temp
   rcd3 = temprd

   if (rcd2==rcd1) then
      grp_mtx(rcd1,j) = rcd1
      rdmtx(rcd1,j)   = rcd3
   else
      call swap(rcd1, rcd2, rcd3, j, N, grp_mtx, rdmtx)
   end if

   end subroutine swap

   subroutine x2f_mpi_crtintra( grp_index, mpicomm_pt, intracomm )
   !
   ! Based on group index (grp_index), create intra-communicator among each group
   ! Input:
   !   grp_index  - group index
   !   mpicomm_pt - communicator within group
   ! Output:
   !   intracomm  - intra-communicator created
   !
   integer, intent(in)    :: grp_index, mpicomm_pt
   integer, intent(inout) :: intracomm
   integer                :: rank_in_mpicomm, ierr

   call MPI_COMM_RANK( mpicomm_pt, rank_in_mpicomm, ierr )
   call MPI_COMM_SPLIT( mpicomm_pt, grp_index, rank_in_mpicomm, intracomm, ierr )

   end subroutine x2f_mpi_crtintra
   
   subroutine x2f_mpi_freecomm( commname )
   !
   ! Any temporary communicators created by x2f_mpi_crtintra or x2f_mpi_crtinter
   ! should be freed after they are no longer needed - it is possible to run out
   !
   integer, intent(inout) :: commname
   integer                :: ierr
   
   if ( commname/=MPI_COMM_NULL ) then
      call MPI_COMM_FREE( commname, ierr )
   end if
   
   end subroutine x2f_mpi_freecomm
   
   subroutine x2f_mpi_setcomms( initial_comm, final_comm, gat, n_spec1, info )
   !
   ! Set a communicator for each of the x2f_mpi attempts, via the info array
   ! (This subroutine was renamed and revamped by SRL 2/2011)
   !
   implicit none
   ! For adaptive mode the following are usually intracomm, mpicomm_uran
   integer, intent(in)    :: initial_comm, final_comm
   integer, intent(in)    :: gat, n_spec1
   integer, intent(inout) :: info(gat, n_spec1)

   info(1:gat-1,5) = initial_comm
   info(gat,5)     = final_comm
   ! Whenever communicator changes, x2f_mpi needs to be reinitialized
   if ( initial_comm /= final_comm )   info(gat,6) = 1
      
   end subroutine x2f_mpi_setcomms
   
   subroutine x2f_mpi_crtinter( grp_index, grp_h, ng, intracomm, mpicomm_pt, &
                                intercomm )
   !
   ! Created inter-communicator between local group and group grp_h
   ! Input:
   !   grp_h  - the group index which local group wants to communicate
   !      ng  - number of group
   !intracomm - intra-communicator within each group
   ! Output:
   !intercomm - inter-communicator between local group and group grp_h
   !
   implicit none
   integer, intent(in)    :: grp_index, grp_h, ng, intracomm, mpicomm_pt
   integer, intent(inout) :: intercomm
   
   integer  :: rank_in_mpicomm, rank_in_grp, rank_in_mpicomm_ary(ng)
   integer  :: mpicomm, newcomm, ierr, leader
   
   call MPI_COMM_RANK( mpicomm_pt, rank_in_mpicomm, ierr )
   call MPI_COMM_RANK( intracomm, rank_in_grp, ierr )
   ! Split the group based on rank_in_grp
   call MPI_COMM_SPLIT( mpicomm_pt, rank_in_grp, rank_in_mpicomm, newcomm, ierr )
   
   if ( rank_in_grp==0 ) then
      ! Gather the leader of each group on rank(rank_in_mpicomm) 0
      call MPI_Gather( rank_in_mpicomm, 1, MPI_INTEGER, rank_in_mpicomm_ary, &
                       1, MPI_INTEGER, 0, newcomm, ierr )
   end if
      
   call x2f_mpi_freecomm( newcomm )
   
   ! Broadcast the leader information to each processor in the simulation
   call MPI_Bcast( rank_in_mpicomm_ary, ng, MPI_INTEGER, 0, mpicomm_pt, ierr )
   ! Create inter-communicator between group grp_index and grp_h
   leader = rank_in_mpicomm_ary(grp_h)
   call MPI_INTERCOMM_CREATE( intracomm,0,mpicomm_pt,leader,1,intercomm,ierr )
                              
   end subroutine x2f_mpi_crtinter
   
   subroutine x2f_mpi_crtinter0( grp_index, grp_h, ng, intracomm, mpicomm_pt, &
                                intercomm )
   !
   ! Created inter-communicator between local group and group grp_h
   ! Input:
   !   grp_h  - the group index which local group wants to communicate
   !      ng  - number of group
   !intracomm - intra-communicator within each group
   ! Output:
   !intercomm - inter-communicator between local group and group grp_h
   !
   implicit none
   integer, intent(in)    :: grp_index, grp_h, ng, intracomm, mpicomm_pt
   integer, intent(inout) :: intercomm
   
   integer  :: rank_in_mpicomm, rank_in_grp, rank_in_mpicomm_ary(ng)
   integer  :: mpicomm, newcomm, ierr, leader
   
   call MPI_COMM_RANK( mpicomm_pt, rank_in_mpicomm, ierr )
   call MPI_COMM_RANK( intracomm, rank_in_grp, ierr )
   ! Split the group based on rank_in_grp
   call MPI_COMM_SPLIT( mpicomm_pt, rank_in_grp, rank_in_mpicomm, newcomm, ierr )
   
   if ( rank_in_grp==0 ) then
      ! Gather the leader of each group on rank(rank_in_mpicomm) 0
      call MPI_Gather( rank_in_mpicomm, 1, MPI_INTEGER, rank_in_mpicomm_ary, &
                       1, MPI_INTEGER, 0, newcomm, ierr )
   end if
   ! Broadcast the leader information to each processor in the simulation
   call MPI_Bcast( rank_in_mpicomm_ary, ng, MPI_INTEGER, 0, mpicomm_pt, ierr )
   ! Create inter-communicator between group grp_index and grp_h
   leader = rank_in_mpicomm_ary(grp_h)
   call MPI_INTERCOMM_CREATE( intracomm,0,mpicomm_pt,leader,1,intercomm,ierr )
                              
   end subroutine x2f_mpi_crtinter0
   
   subroutine x2f_mpi_getxfsize( nv_local, intracomm, intercomm, nv_nolocal )
   !
   ! Exchange the size of xf on the corresponding processor for 2 groups
   !
   implicit none
   integer, intent(in)     :: nv_local, intracomm, intercomm
   integer, intent(inout)  :: nv_nolocal
   
   integer  :: rank_in_grp, ierr, stats(MPI_STATUS_SIZE)
   
   call MPI_COMM_RANK( intracomm, rank_in_grp, ierr )
   call MPI_SEND( nv_local, 1, MPI_INTEGER, rank_in_grp, 1, intercomm, ierr )
   call MPI_RECV( nv_nolocal, 1, MPI_INTEGER, rank_in_grp, 1, intercomm, stats, &
                  ierr )
   
   end subroutine x2f_mpi_getxfsize
   
   subroutine x2f_mpi_xfswap( nv_local, ldxf, xf_local, intracomm, &
                              intercomm, nv_nolocal, xf_nolocal )
   !
   ! Exchange the particle compositions on the corresponding processor 
   ! for 2 groups 
   !
   implicit none
   integer, intent(in)             :: nv_local, ldxf, intracomm, intercomm
   integer, intent(in)             :: nv_nolocal
   real(kind(1.d0)), intent(in)    :: xf_local(ldxf, nv_local)
   real(kind(1.d0)), intent(inout) :: xf_nolocal(ldxf, nv_nolocal)  
      
   integer  :: rank_in_grp, ierr, stats(MPI_STATUS_SIZE), req
   
   call MPI_COMM_RANK( intracomm, rank_in_grp, ierr )
   call MPI_ISEND( xf_local, ldxf*nv_local, MPI_DOUBLE_PRECISION, rank_in_grp, 2, &
                   intercomm, req, ierr )
   call MPI_REQUEST_FREE( req, ierr )
   call MPI_IRECV( xf_nolocal, ldxf*nv_nolocal, MPI_DOUBLE_PRECISION, & 
                   rank_in_grp, 2, intercomm, req, ierr )
   call MPI_WAIT( req, stats, ierr )
   
   end subroutine x2f_mpi_xfswap
   
   subroutine x2f_mpi_pair( i, ng, grp_g, grp_h )
   !
   ! In the i-th pairing, for the given group grp_g, find the group grp_h with
   ! which it should be paired, so that the pairing is unique for 1 <= i <= ng.
   ! As usual, ng is the total number of groups.
   !
   ! Note that this function must be an involution, i.e., it has to be its own
   ! inverse, for all i.  This is because when A is paired with B, B must also
   ! be paired with A.  Thus f_i(f_i(A)) = f_i(B) = A.  (Comments by SRL, 6/26/10)
   !
   implicit none
   integer, intent(in)    :: i, ng, grp_g
   integer, intent(inout) :: grp_h
   
   integer  :: p
   
   if ( i==1 ) then
      ! This the identity pairing
      grp_h = grp_g
      return
   else
      ! For i>1, we don't want grp_g to be paired with itself again.
      ! This turns out to be a problem whenever grp_g==p, as defined here
      if ( mod(i,2)==0 ) then
         p = (ng+i)/2
      else
         p = (i+1)/2
      end if
   
      if (grp_g==i) then
         ! This special case comes first to handle grp_g==ng correctly for i=ng
         grp_h = 1
      elseif (grp_g==p) then
         ! The "problem" group p (which depends on i) will be paired with group ng
         grp_h = ng
      elseif (grp_g==ng) then
         grp_h = p
      elseif (grp_g<i) then
         ! Done with special cases; for most grp_g, one of the next two will apply
         grp_h = i-grp_g+1
      elseif (grp_g>i) then
         grp_h = i-grp_g+ng
      end if
      return
   end if
   
   end subroutine x2f_mpi_pair
   
   subroutine construct_ary( dist_strat, atmps, n, gat, op_frq, n_spec1, &
                             n_spec2, nolocal, info, rinfo, info_xf, rinfo_xf )
   implicit none
   integer, intent(in)    :: n, dist_strat(n), atmps(n), gat, op_frq, n_spec1, &
                             n_spec2, nolocal
   integer, intent(inout) :: info(gat, n_spec1), info_xf(:)
   real(kind(1.d0)), intent(inout) :: rinfo(gat,n_spec2), rinfo_xf(:)
   integer                :: igat, i, j  

   info     = 0
   rinfo    = 0.d0
   info_xf  = 0
   rinfo_xf = 0.d0   

   ! Sep 16, 2007
   info_xf(1) = atmps(2)
   ! Sep 16, 2007
   igat = 0
   do i = 1, n
      do j = 1, atmps(i)
         igat = igat + 1
		 info(igat,1) = i
		 info(igat,2) = dist_strat(i)
		 info(igat,7) = op_frq
         ! Sep 24, 2007
		 if (j==1) info(igat,6) = 1
		 if (info(igat,2)==4) info(igat,3)=nolocal
	  end do
   end do

   end subroutine construct_ary

   subroutine x2f_mpi_tc( nx, nf, ldxf, nv, mpicomm_pt, grp_index, Tc, lu_log )

   implicit none
   integer, intent(in)      :: nx, nf, ldxf, nv, mpicomm_pt, grp_index, lu_log
   real, intent(inout)      :: Tc
   
   external                      :: x2f_cirxn, x2f_setmethod 
   real(kind(1.d0)), allocatable :: xf(:,:) 
   real(kind(1.d0)), allocatable :: rinfo(:,:)
   real(kind(1.d0))              :: rinfo_xf(1), stats(100)
   integer, allocatable          :: info(:,:)
   integer                       :: info_xf(1)
   integer  :: dist_strat(3), atmps(3), k_pos, gat, n, n_spec1, n_spec2, op_frq
   integer  :: nolocal, i
   
   real(kind(1.d0)) :: starttime, endtime, cpu2, cpu1
   real             :: Tc_temp 
   integer          :: ierr, nproc, nv_test

   call MPI_ALLREDUCE( nv, nv_test, 1, MPI_INTEGER, MPI_SUM, mpicomm_pt, ierr )
   call MPI_COMM_SIZE( mpicomm_pt, nproc, ierr )
   nv_test = nv_test/nproc
   allocate(  xf(ldxf,nv_test)  )

   dist_strat(1:3) = (/0, 0, 1/)
   atmps(1:3)      = (/0, 0, 1/)
   k_pos           = 1
   op_frq          = -1
   nolocal         = 0
   gat             = sum( atmps )
   n               = size( atmps, 1 )
   n_spec1         = 10
   n_spec2         = 50
   allocate( info(gat,n_spec1)  )
   allocate( rinfo(gat,n_spec2) )
   call construct_ary( dist_strat, atmps, n, gat, op_frq, n_spec1, n_spec2, &
                       nolocal, info, rinfo, info_xf, rinfo_xf )
   call x2f_mpi_setcomms( mpicomm_pt, mpicomm_pt, gat, n_spec1, info )
   info_xf = -1
   starttime = MPI_Wtime()
   call cpu_time( cpu1 )   
   do i = 1, 10
   xf = 1.d0
   call x2f_mpi( nx, nf, ldxf, nv_test, xf, k_pos, info, rinfo, info_xf, &
                 rinfo_xf, 1, n_spec1, n_spec2, x2f_cirxn, x2f_setmethod, stats )
   end do
   call cpu_time( cpu2 )
   endtime = MPI_Wtime()
   
   Tc_temp = (endtime-starttime-(cpu2-cpu1))/(10*nv_test)
   if(a_nolog == 0) write(lu_log,*) 'Tc_temp', Tc_temp, nv_test
   call MPI_REDUCE( Tc_temp, Tc, 1, MPI_REAL, MPI_MAX, 0, mpicomm_pt, ierr )
   
   deallocate( info  )
   deallocate( rinfo )
   deallocate(  xf   )

   end subroutine x2f_mpi_tc
   
   subroutine x2f_mpi_newmatch( mpicomm_pt, ips, ng, g, r, A_max, t_C, lu_log )
   !
   ! Estimate wall clock time for final stage for pairing and not pairing cases
   ! respectively
   ! Some inputs:
   !      r    - parameter for power law
   !    A_max  - the maximum number of table entries per processor
   !    t_C  - communication time
   !
   implicit none
   integer, intent(in)    :: mpicomm_pt, ips, ng, g, lu_log
   real, intent(in)       :: r, A_max
   real, intent(inout)    :: t_C
   real, save             :: T_dec = 0.0, phi = 0.5, rat = 2.0
   
   integer  :: rank_in_grp, ierr, rank, newcomm
   integer  :: i, j, N_i, N_j, N_i_pri, cnter, edge(ng*(ng-1)/2,2)
   real     :: t_Ri, t_Rj, t_Ri_pri, t_Rj_pri, t_Di_pri, t_Dj_pri, temp
   real     :: P_Di, P_Dj, P_Di_pri, P_Dj_pri, P_ij, P_ji, P_ii, P_jj
   real     :: a_i, a_i_pri, a_j, a_j_pri, a_ij, a_ji, Nij, N_F
   real     :: T_NP(ng), T_P((ng-1)*ng/2), T_NP_max, T_P_max, P_D_PR1_i, P_D_PR1_j
   real     :: tempI, tempJ, tempK, tempL, t_Rij, t_Rji
   
   call MPI_COMM_RANK( mpicomm_pt, rank, ierr )
   if ( rank==0 ) then
   !
   ! Estimate time on the next stage
   ! Variable names having "name_pri" extension means they are estimated values
   !
   ! Not pairing case:
   ! 1/ Retrieval time, direct evaluation time and number of particles are
   !    estimated based on previous stage's information
   ! 2/ a_i_pri: number of table entries estimated for final stage; we assume 
   !    that all processors have the maximum number of table entries:
   !             a_i_pri = g*A_max
   ! 3/ P_Di_pri: fraction of direct evaluation estimated for next stage; we 
   !    assume that P_Di varies inversely with the r-th power of the number of 
   !    table entries:
   !             P_Di_pri = P_Di*(a_i/a_i_pri)^r
   ! 4/ the estimate time is then
   !             T' = T_NP(i) = N_i_pri*t_Ri_pri+N_i_pri*p_Di_pri*t_Di_pri
   !
   do i = 1, ng
      t_Ri_pri = tRij_mtx(i,i)
      t_Di_pri = t_Di_ary(i)
      a_i      = a_i_ary(i)
      a_i_pri  = g*A_max
      P_Di     = P_Di_ary(i)
      P_Di_pri = P_Di*(a_i/a_i_pri)**r
      N_i_pri  = N_i_ary(i)
      T_NP(i)  = N_i_pri*t_Ri_pri + N_i_pri*p_Di_pri*t_Di_pri
      if(a_nolog==0) write(lu_log,'(A,6e15.6,2I8)') 'T_NP', &
           T_NP(i), t_Ri_pri, t_Di_pri, t_Ri_pri*N_i_pri, N_i_pri*p_Di_pri*t_Di_pri, P_Di_pri
   end do
   T_NP_max = maxval(T_NP(1:ng))
   if(a_nolog==0) write(lu_log,'(A,E15.6)') 'T_NP_max', T_NP_max
   
   !
   ! Pairing case:
   ! 1/ Retrieval time, direct evaluation time and number of particles are
   !    estimated based on previous stage's information
   ! 2/ a_i_pri and a_j_pri are the total number of table entries added based
   !    on particles from groups i and j (regardless of whether these entries
   !    are in ISAT tables stored in group i or j).
   !    
   !    We assume:
   !       Define a_ij as the "effective number of table entries" (on the final
   !       stage in the group formed by pairing groups i and j) for retrieving 
   !       particles from the original group i. 
   !       If the particles on group i and j are disjoint, then a_ij = A_max
   !       If they are coincident, then a_ij = 2*A_max
   !       In general, we model a_ij as
   !                 a_ij = (A_max + A_max * (P_ij/max(P_ij,1-P_Di)))*g
   ! 3/ P_Di_pri, the probability of a group i particle requiring DE on the next
   !    stage is estimated as
   !                 P_Di_pri = P_Di*(a_i/a_ij)^r
   ! 4/ t_Ri, t_Rj are the average wall clock time for a retrieve;
   !    t_Rij, t_Rji are
   ! 5/ the detail information about how to estimate total time can be found in Liuyan 
   !    thesis
   !                  
   T_P   = 0   
   cnter = 0
   do i = 1, ng-1
      a_i      = a_i_ary(i)
      P_Di     = P_Di_ary(i)
      P_D_PR1_i= P_Di_PR1(i)
      t_Ri     = tRij_mtx(i,i)
      t_Di_pri = t_Di_ary(i)
      N_i      = N_i_ary(i)
      temp     = 2*g*A_max
      do j = i+1, ng
         cnter    = cnter + 1
         a_j      = a_i_ary(j)         
         P_Dj     = P_Di_ary(j)
         P_D_PR1_j= P_Di_PR1(j)
         t_Rj     = tRij_mtx(j,j)
         t_Rij    = tRij_mtx(i,j)
         t_Rji    = tRij_mtx(j,i)
         t_Dj_pri = t_Di_ary(j)
         N_j      = N_i_ary(j)         
         P_ij     = P_ij_mtx(i,j)
         P_ji     = P_ij_mtx(j,i)
         
         P_ii     = P_ij_mtx(i,i)
         P_jj     = P_ij_mtx(j,j)
                  
         a_i_pri  = temp*a_i/(a_i+a_j)
         a_j_pri  = temp*a_j/(a_i+a_j)
         
         a_ij     = a_i_pri + a_j_pri * P_ij/(max(P_ij,1-P_Dj))
         a_ji     = a_j_pri + a_i_pri * P_ji/(max(P_ji,1-P_Di))
         
         P_Di_pri = P_Di * (a_i/a_ij)**r
         P_Dj_pri = P_Dj * (a_j/a_ji)**r
         
         Nij      = (N_i+N_j)/2.0
         N_F      = (N_i*P_Di_pri+N_j*P_Dj_pri)/2.0
         
         if ( N_i/real(N_j)>rat.or.N_j/real(N_i)>rat ) then
            tempI = N_i*t_Ri/2.0 + N_j*t_Rij/2.0 + N_j*t_C/2.0
            tempJ = N_j*t_Rj/2.0 + N_i*t_Rji/2.0 + N_i*t_C/2.0
            tempK = N_i*(1-P_ij)*t_Ri/2.0 + N_j*(1-P_jj)*t_Rji/2.0
            tempL = N_j*(1-P_ji)*t_Rj/2.0 + N_i*(1-P_ii)*t_Rij/2.0
            if (i==1.and.j==2 .and. a_nolog == 0) write(lu_log,'(A,e15.6)') 'enter1'
         else
            tempI = N_i*t_Ri
            tempJ = N_j*t_Rj
            tempK = N_j*P_D_PR1_j*(t_Rji+t_C)*g/(2*g-1)
            tempL = N_i*P_D_PR1_i*(t_Rij+t_C)*g/(2*g-1)
            if (i==1.and.j==2 .and. a_nolog == 0) write(lu_log,'(A,e15.6)') 'enter2'
         end if   
         
         T_P(cnter) = max(tempI,tempJ) + max(tempK,tempL) + N_F * max(t_Di_pri, t_Dj_pri)
                                    
         edge(cnter,1) = i
         edge(cnter,2) = j
         
         if(a_nolog==0) write(lu_log,'(A,20e15.6)') 'T_P2', T_P(cnter), max(tempI, tempJ), max(tempK, tempL), N_F * max(t_Di_pri, t_Dj_pri)
      end do 
   end do 

   call findmat( ng, cnter, edge, T_P, match, T_P_max )
   if(a_nolog==0) write(lu_log,'(A,e15.6)') 'T_P_max', T_P_max
   
   if ( T_NP_max<T_P_max ) then
      match = 0
      T_dec = T_NP_max
   elseif ( T_NP_max>T_P_max ) then
      if ( ips==1 ) then
         T_dec = T_P_max
      else
         temp = phi*T_NP_max + (1.-phi)*T_dec
         if(a_nolog==0) write(lu_log,'(A,10e15.6)') 'phi', phi, temp, T_NP_max, T_dec
         if ( temp>T_P_max ) then
            T_dec = T_P_max
         else
            match = 0
            T_dec = T_NP_max
         end if
      end if
   end if
      
   if(a_nolog==0) write(lu_log,'(A,10000I4)') 'match', match
   
   end if

   end subroutine x2f_mpi_newmatch 
   
   subroutine findmat( nv, M, edge, wt, mate, wtmin )
   !
   ! Given nv(nv can be divided by 2) projects, partition these into nv/2 pairs 
   ! of objects. There is a known positive cost(i,j) associated with pairing 
   ! objects i and j. The partition selected minimizes the maximum value of 
   ! cost(i,j) (for i and j being a pair).
   !
   ! The problem is close related to the Maximum Cardinality Match Problem (MCM): 
   ! Given an undirected (and unweighted) graph, determine the largest possible 
   ! matching. A "matching" is a set M of edges such that each vertex belongs to 
   ! at most one member of M. A matching M is said to be "perfect" if each vertex 
   ! belongs to exactly one member of M (i.e., all the vertices are used).
   ! 
   ! However our graph is weighted and our goal is to find the perfect matching  
   ! that minimizes the maximum weight among all edges of M. The algorithm for 
   ! doing this that makes use of MCM is the following:
   !
   ! 1/ Let G=(V,E) be the graph. Sort the edges by weight. The idea is to use 
   !    binary search on the sorted list of edges. Let w be the middle weight in 
   !    the sorted list of edges.
   !
   ! 2/ Form G', an unweighted graph by retaining only the edges with weight less 
   !    than w.
   !
   ! 3/ Run the MCM algorithm to see if G' has a perfect match. If it does, make w 
   !    smaller (as in binary search). If it does not, make w larger (as in binary
   !    search)
   !
   ! 4/ If the binary search range is down to a single weight, report that weight; 
   !    otherwise go to Step 2.
   !
   ! The best algorithm for MCM runs in time O(n^3). Thus, the above algorithm would
   ! run in time O(n^3*log n).
   !
   ! The website with links to free software for MCM:
   !     http://www.cs.sunysb.edu/~algorith/files/matching.shtml
   !
   !   nv    - number of vertices
   !    M    - number of edges
   !  edge   - edge(i,1) and edge(i,2) are indices of the i-th edge
   !   wt    - wt(i), weight for the i-th edges
   !  mate   - mate(i), the mate for the i-th vertex
   !   
   !   dgr   - dgr(i), edge degree of the i-th vertex  
   !   nds   - index of the adjacent vertex
   !
   implicit none
   integer, intent(in)    :: nv, M
   real, intent(inout)    :: wt(M), wtmin
   integer, intent(inout) :: edge(M,2), mate(nv)
   
   integer  :: Mh, Ml, Mt, istep, flag, mate0(nv), dgr(nv), nds(nv*(nv-1))
   
   ! sort the edges by weight
   call sortwt( M, wt, edge )
   
   ! use binary search on the sorted list
   Mh = M
   Ml = 1
   Mt = M
   istep = 0
   do while ( istep==0.or.Mh-Ml>0 )
      istep = istep + 1
      !print *, 'begin:', istep, Mh, Ml, Mt
   
      ! change the data form to the one required by c subroutine
      call cvtdata( nv, M, Mt, edge, dgr, nds )

      ! call c subroutine
      call trymatch( nv, dgr, nds, mate0 )
            
      flag = 1
      ! test to see if the matched indices are valid or not
      call testmatch( nv, mate0, flag )
      
      if ( Mh-Ml>1 ) then
         if ( flag==1 ) then
            Mh = Mt
            wtmin = wt(Mt)
            mate  = mate0
         else
            Ml = Mt
         end if
      else
         if ( flag==1 ) then
            Mh = Mt - 1
            wtmin = wt(Mt)
            mate  = mate0
         else
            Ml = Ml + 1
         end if
      end if
      
      Mt = (Mh+Ml)/2
      !print *, 'end', istep
   end do
   
   end subroutine findmat    
   
   subroutine sortwt( M, wt, edge )
   !
   ! Sort the edges by weight
   !
   implicit none
   integer, intent(in)     :: M
   integer, intent(inout)  :: edge(M,2)
   real, intent(inout)     :: wt(M)
   
   integer  :: i, j, index
   real     :: tmp, tmpary(2)
   
   do i = 1, M - 1
      tmp = wt(i)
      index = i
      do j = i + 1, M
         if ( tmp>wt(j) ) then
            index = j
            tmp   = wt(j)
         end if
      end do
      if ( index>i ) then
         wt(index) = wt(i)
         wt(i)     = tmp
         
         tmpary(1:2) = edge(index,1:2)
         edge(index,1:2) = edge(i,1:2)
         edge(i,1:2) = tmpary(1:2)
      end if
   end do
   
   end subroutine sortwt
   
   subroutine cvtdata( nv, M, Mt, edge, dgr, nds )
   !
   ! Change the data to the form required by c subroutine
   !
   !   wt    - wt(i), weight for the i-th edges
   !  mate   - mate(i), the mate for the i-th vertex
   !   
   !   dgr   - dgr(i), edge degree of the i-th vertex  
   !   nds   - index of the adjacent vertex
   !         - nds(sum(dgr(1:i-1))+1:sum(dgr(1:i))) stores all the vertices
   !           that connected with the i-th vertex
   !
   implicit none
   integer, intent(in)    :: nv, M, Mt
   integer, intent(in)    :: edge(M,2)
   integer, intent(inout) :: dgr(nv), nds(nv*(nv-1))
   
   integer  :: i, j, cnt, vertex
   
   nds = 0
   cnt = 0
   do i = 1, nv
      dgr(i) = 0
      do j = 1, Mt
         vertex = 0
         if ( edge(j,1)==i ) then
            vertex = edge(j,2)
         elseif ( edge(j,2)==i ) then
            vertex = edge(j,1)
         end if
         
         if ( vertex/=0 ) then
            dgr(i) = dgr(i) + 1
            cnt = cnt + 1
            nds(cnt) = vertex
         end if
      end do
   end do 
   
   end subroutine cvtdata
   
   subroutine testmatch( nv, mate, flag )
   !
   ! check the matched vertices are valid or not
   !
   implicit none
   integer, intent(in)    :: nv, mate(nv)
   integer, intent(inout) :: flag
   
   integer  :: i
   
   do i = 1, nv
      if ( mate(i)==0 ) then
         flag = 0
         exit
      end if
   end do
   
   end subroutine
   
   subroutine x2f_mpi_setcomm( mpicomm_pt, ng, g, match, grp_index, istep )
   !                      
   ! Based on the match information, reset the group index grp_index
   ! Input:
   !   ng   - number of group
   ! match  - size of ng, match(i) indicates group i is paired with group match(i)
   ! Output:
   ! grp_index - new group index 
   !
   implicit none
   integer, intent(in)    :: mpicomm_pt, match(ng), istep
   integer, intent(inout) :: ng, g, grp_index
   
   integer  :: mpicomm, ierr, ngd2, rank 
   integer, allocatable  :: mate1(:), mate2(:)
   
   call MPI_Bcast( match, ng, MPI_INTEGER, 0, mpicomm_pt, ierr )
   
   if ( match(1)==0 ) then
      ! No pairing is performed
      return
   else
      if ( mod(ng,2)/=0 ) stop 'Err: using unapproiate number of processors'
      ngd2 = ng/2
      if ( allocated(mate1) ) deallocate( mate1 )
      if ( allocated(mate2) ) deallocate( mate2 )
      allocate( mate1(ngd2) )
      allocate( mate2(ngd2) )
      
      call cvt_match( ngd2, match, mate1, mate2 )
      call set_grpindex( ngd2, mate1, mate2, grp_index )
      ng = ngd2
      g  = g*2
      
   end if
   
   end subroutine x2f_mpi_setcomm
   
   subroutine cvt_match( ngd2, mate, mate1, mate2 )
   !
   !  input:
   !    mate  - size of ngd2*2, mate(i) indicates group i is paired with group 
   !            mate(i). e.g., we might have the following data set
   !                  i:  1,  2,  3,  4,  5,  6
   !            mate(i):  3,  6,  1,  5,  4,  2
   !  output: 
   !    mate1 - size of ngd2, stores smaller group index 
   !    mate2 - size of ngd2, stores larger group index
   !            mate1(i) and mate2(i) form one new pair.
   !            e.g., based on the previous example
   !            mate1(i):  1,  2,  4
   !            mate2(i):  3,  6,  5
   !

   implicit none
   integer, intent(in)    :: ngd2, mate(ngd2*2)
   integer, intent(inout) :: mate1(ngd2), mate2(ngd2)
   
   integer  :: index1, index2, i
   logical  :: appear
   
   index1 = 1
   index2 = 1
   mate1(index1) = index2
   mate2(index1) = mate(index2)
   do while (index1<ngd2)
      index1 = index1+1
      appear = .true.
      do while (appear)
         index2 = index2 + 1
         do i = 1, index1-1
            if ( index2==mate1(i).or.index2==mate2(i) ) then
               Exit
            end if
         end do
         if ( i==index1 ) appear=.false.
      end do
      mate1(index1) = index2
      mate2(index1) = mate(index2)
   end do

   end subroutine cvt_match
   
   subroutine set_grpindex(ng,mate1,mate2,color)
   !
   ! Reset the group index(color) for new group
   ! Input:
   !   ng   - number of new group 
   ! mate1  - size of ng, mate1(i) stores one group index belonged to new group i
   ! mate2  - size of ng, mate2(i) stores another group index belonged to new 
   !          group i. Here group mate1(i) and mate2(i) belong to new group i which
   !          is the result of pairing.
   ! 
   ! Output:
   ! color  - reset group index 
   !        

   implicit none
   integer, intent(in)    :: ng, mate1(ng), mate2(ng)
   integer, intent(inout) :: color
   
   integer  :: color_old, i

   color_old = color
   do i = 1, ng
      if ( color_old==mate1(i).or.color_old==mate2(i) ) then
         color = i
         Exit
      end if
   end do

   end subroutine set_grpindex
   
   subroutine fileinit( lu_op, lu_log, lu_time )
   !
   ! Initialize the file which can be written to 
   !
   implicit none
   integer, intent(inout)  :: lu_op, lu_log, lu_time
   
   integer                 :: rank, ierr
   character(30)           :: filenm, head, tail
   
   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

   if( rank == 0 ) then
      call isat_lu( lu_op )
      open( lu_op, file='isatmp.op', status='replace', action = 'write' )
   endif

   if( a_nolog == 0 ) then
      ! create log files only if a_nolog = 0
      head = 'isatmp'
      tail = 'log'
      call x2fmpi_file_name(head,rank,-1,tail,filenm)

      call isat_lu( lu_log )
      open( lu_log, file=filenm, status='replace', action = 'write' )
   endif
   
   if( a_nolog < 2 ) then
      ! output time stats only if a_nolog < 2
      head = 'time'
      tail = 'op'
      call x2fmpi_file_name(head,rank,-1,tail,filenm)
      call isat_lu( lu_time )
      open( lu_time, file=filenm, status='replace', action = 'write' )
   endif
   
   end subroutine fileinit
   
   subroutine ifmsstg( mpicomm_pt, intracomm, g, msstg, lu_log )
   !
   !Check to see if this is measuring stage or not
   ! Input:
   !   mpicomm_pt  - communicator within partition
   !
   ! Output:
   !    msstg      - measuring stage or not
   !
   integer, intent(in)    :: mpicomm_pt, intracomm, g, lu_log
   logical, intent(inout) :: msstg
   
   integer  :: ierr

   call reach_astar( mpicomm_pt, intracomm, g, msstg, lu_log )
   call MPI_BCAST( msstg, 1, MPI_LOGICAL, 0, mpicomm_pt, ierr )
      
   end subroutine ifmsstg
   
   subroutine reach_astar( mpicomm_pt, intracomm, g, msstg, lu_log )
   
   integer, intent(in)    :: mpicomm_pt, intracomm, g, lu_log
   logical, intent(inout) :: msstg
   
   real     :: a_i, a_i_sum, a_star, a_max
   integer  :: rank_in_grp, ierr, newcomm, rank 
   
   msstg = .false.
   
   call ci_stats( 12, a_i )
   call ci_stats( 13, a_star )
   call MPI_COMM_RANK( mpicomm_pt, rank, ierr )

   if ( g>1 ) then
      call MPI_REDUCE( a_i, a_i_sum, 1, MPI_REAL, MPI_SUM, 0, intracomm, ierr )
      a_i = a_i_sum
      
      call MPI_COMM_RANK( intracomm, rank_in_grp, ierr )
      call x2f_mpi_crtintra( rank_in_grp, mpicomm_pt, newcomm )
      
      if ( rank_in_grp==0 ) then
         call MPI_ALLGATHER( a_i, 1, MPI_REAL, a_i_ary, 1, MPI_REAL, newcomm, ierr )     
         a_max = maxval(a_i_ary)
      end if
         
      call x2f_mpi_freecomm( newcomm )
   
   elseif ( g==1 ) then
      !call MPI_REDUCE( a_i, a_max, 1, MPI_REAL, MPI_MAX, 0, intracomm, ierr )
      call MPI_ALLGATHER( a_i, 1, MPI_REAL, a_i_ary, 1, MPI_REAL, mpicomm_pt, ierr )
      a_max = maxval(a_i_ary)      
   end if
   
   if ( rank==0 ) then
      if ( a_max>=g*a_star ) then
         msstg = .true.
      else
         msstg = .false.
      end if
      if(a_nolog==0) write(lu_log,'(A,10000e15.6)') 'a_i_ary ', a_i_ary
   end if 
   
   end subroutine reach_astar
   
   subroutine x2f_mpi_tabfd( mpicomm_pt, lu_log, tabfd )
   !
   ! Check to see if all of the tables are fully developed or not
   !
   integer, intent(in)     :: mpicomm_pt, lu_log
   logical, intent(inout)  :: tabfd
   
   real     :: n_q, n_r, n_g, n_a, n_d
   real     :: temp_fd, thrsh_tabfd
   integer  :: fd, fdmin, ierr  
   
   call isat_stats( 2, n_q, n_r, n_g, n_a, n_d, lu_log )
   thrsh_tabfd = 1.e-4
   temp_fd = (n_g + n_a)/n_q
   if ( temp_fd < thrsh_tabfd ) then
      fd = 1
   else
      fd = 0
   end if
   
   call MPI_ALLREDUCE( fd, fdmin, 1, MPI_INTEGER, MPI_MIN, mpicomm_pt, ierr )
   if ( fdmin==1 ) then
      tabfd = .true.
   else
      tabfd = .false.
   end if
         
   end subroutine x2f_mpi_tabfd
   
end module x2f_mpi_adapt
