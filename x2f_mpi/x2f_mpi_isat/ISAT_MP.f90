!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the x2f_mpi      !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ISAT_MP( nx, nf, ldxf, nv, xf, mode, nSR2, info_imp, istep )
!
! ------------ OVERVIEW
!
! The software ISAT_MP is an extension of the ISAT algorithm (Pope 1997) to the
! multi-processor enviroment for efficient function evaluation. It tries to
! minimize the wall clock time spent in the evaluation of the computationally
! expensive function f(x) for an ensemble of x by using the x2f_mpi library
! and the serial ISAT algorithm. Given an ensemble of x, ISAT_MP returns the
! approximate value of f for each x which (with high probability) is within a
! specified tolerance of the exact value.
!
! ISAT_MP calls x2f_mpi, which contains three different distribution strategies
! and is responsible for distributing the ensemble of x among some or all of the
! processors in the simulation, and then x2f_mpi calls ISAT, which performs the
! function evaluation.
!
! In the multi-processor environment, before performing function evaluations,
! the ensemble of x may be distributed among some or all of the processors used
! in the simulation. Therefore some x originally located on one processor might
! be sent to the other processors and some x might stay locally. Then function
! evaluations are performed for all of the x resident on the processor at that
! time using the local ISAT table. After function evaluations, each f(x) is sent
! back to the processor from which the corresponding x originated.
!
! To evaluate f(x), the ISAT algorithm provides a sequence of methods, namely
! PR (primary retrieve), PR+SR (primary and secondary retrieve), and complete
! ISAT (including primary retrieve, secondary retrieve, grow, add and direct
! evaluation). The first two methods(PR and PR+SR) may or may not be successful
! in function evaluation, but they are cheap compared to the last one
! (complete ISAT) which always succeeds.
!
! In general, the computationally cheap methods are tried first to perform
! function evaluation. In the multi-processor environment, the retrieve
! methods can be tried on ISAT tables on different processors. By setting the
! values of atmps array in ISAT_MP, the user can specify the number of attempts
! for each of the methods, namely PR, PR+SR and complete ISAT. Setting atmps(i)=0
! corresponds to not using the method. Notice that more attempts on retrieval
! methods do not necessarily mean more efficiency.
!
! In x2f_mpi, a number of different distribution strategies are implemented,
! including: purely local processing (PLP), uniform random distribution (URAN),
! and preferential distribution (PREF). Different distribution strategies are
! referenced by the index(e.g, 0 refers to PLP, 1 refers to URAN and 4 refers to
! PREF). In ISAT_MP, the user can specify the distribution strategies by setting
! the values in array dist_strat. By specifying dist_strat(i), the distribution
! strategy for the i-th method is set.
!
! When using ISAT_MP, the key vectors one needs to specify are the number of
! attempts for each method (atmps) and the distribution strategies for each method
! (dist_strat). In the current implementation, the user chooses presets for these
! vectors by setting mode to a value from 1 to 7.  A future implementation should
! allow advanced users to specify values of atmps and dist_strat when mode=0.
!
! ------------ INPUT
!
!  nx        - length of x vectors
!  nf        - length of f vectors
!  ldxf      - the leading dimension of the array xf ( ldxf >=max(nx,nf) )
!  nv        - number of vectors of x
!  xf        - array which contains in xf(1:nx,1:nv) the vectors of x on input,
!            - and the vectors of f on output (see below)
! mode       - integer determing the action to be taken (see below)
! nSR2       - integer determing the maximum number of secondary retrieve attempts
! info_imp   - integer array dimensioned 10. The contents of info_imp(1:2)
!              specify control parameters for using PREF strategy.
! info_imp(1)= nolocal = 1, do not try ISAT table on local processor
!                      = 0, attempt to retrieve from ISAT table on local processor
! info_imp(2)= att_max, maximum number of attempts for retrieving method
! info_imp(3)= adap = 1, use adaptive parallel ISAT strategy
!                   = 0, use pre-set strategies for parallel ISAT (see below)
! info_imp(4)= set_pt > 0, split MPI_COMM_WORLD into partitions of this size
!                     = 0, do not split MPI_COMM_WORLD (added 2/2011)
! info_imp(5)= tstg1  = time spent in stage 1 of PURAN (in minutes)
!
! info_imp(10)= nolog - control output files
!                     = 0, log all data
!                     > 0, disable isatmp_# log files
!                     > 1, disable time_# log files

! info_imp(6:9) ---- reserved for future use
!
! ------------ OUTPUT
!
!  xf        - array which contains in xf(1:nf,1:nv) the vectors of f
!
! ------------ MODE
!
! Various pre-set strategies can be chosen by specifying a "mode" (1<=mode<=7).
! (FUTURE) To create a user-set strategy, set mode=0 and specify the distribution
! strategies (dist_strat) and number of attempts (atmps) for each method; also,
! modify k_pos to satisfy certain requirements (see below).

!    mode   = -1 - only initialize partitions and return
!
!    mode    = 1 - PLP (purely local processing)
!                  with dist_strat(3)=(/0,0,0/) and
!                       atmps(3)     =(/0,0,1/) specified
!            = 2 - URAN (uniform random distribution)
!                  with dist_strat(3)=(/0,0,1/) and
!                       atmps(3)     =(/0,0,1/) specified
!            = 3 - PREF(nproc)/URAN with nolocal=0
!                  with dist_strat(3)=(/0,4,1/) and
!                       atmps(3)     =(/0,nproc,1/) specified
!                  nproc - the number of processor used in the simulation
!                - PREF(nproc-1)/URAN with nolocal=1
!                  with dist_strat(3)=(/0,4,1/) and
!                       atmps(3)     =(/0,nproc-1,1/) specified
!                  nproc - the number of processor used in the simulation
!            = 4 - PREF(ntmp)/URAN with nolocal=0 and attempts limited to att_max
!                  with dist_strat(3)=(/0,4,1/) and
!                       atmps(3)     =(/0,ntmp,1/) specified
!                  ntmp - determined by the number of processor(nproc), maximum
!                         number of limits(att_max) and nolocal
!                - PREF(ntmp)/URAN with nolocal=1 and attempts limited to att_max
!                  with dist_strat(3)=(/0,4,1/) and
!                       atmps(3)     =(/0,ntmp,1/) specified
!                  ntmp - determined by the number of processor(nproc), maximum
!                         number of limits(att_max) and nolocal
!            = 5 - ? (documentation needs to be added here)
!            = 6 - PLP at first, switching over to URAN within partitions, after
!                  some specified fraction of the table is full (added 2/2011)
!            = 7 - PLP at first, switching over to PREF/URAN within partitions,
!                  after some specified fraction of the table is full (added 2/2011)
!
! ------------ PARAMETERS
!
! ------------ (FUTURE) Parameters must be specified by user only when mode=0
! dist_strat - distribution strategies used for each method, dist_strat(1:3)
! dist_strat(i)
!            - indicates the distribution used for i-th method
!            = 0, PLP
!            = 1, URAN
!            = 4, PREF
! atmps      - number of attempts for each method, atmps(1:3)
! atmps(1)   - set the number of attempts for PR method
! atmps(2)   - set the number of attempts for PR+SR method
! atmps(3)   - set the number of attempts for complete ISAT method
! k_pos      - some integer value 1<=k_pos<=min(nx,nf)
!              the value of xf(k_pos,i) is used to check if the particle has been
!              evaluated or not. if xf(k_pos,i)>0, the particle has been evaluated,
!              otherwise it's not! For a general problem, it requires that there
!              exists one direction whose value is strictly positive otherwise the
!              user needs to expand the x array to include such direction. For
!              combustion problems, usually k_pos=nx-2 which corresponds to
!              enthalpy.
! op_frq     > 0, integer, output statistics from x2f_mpi to x2fstats_#.out every
!              op_frq step
!            < 0, no output
!
! ------------ Other parameters
!    gat     - global number of attempts
!     n      - total number of methods
!  nspec1    - used to specify the size of second dimension of info array
!  nspec2    - used to specify the size of second dimension of rinfo array
!att_max_deft- default maximum number of attempts
!

use x2f_mpi_grpdata
use x2f_mpi_adapt
implicit none

integer, parameter               :: att_max_deft = 8
integer, intent(in)              :: nx, nf, ldxf, nv, mode, nSR2, istep
integer, intent(inout)           :: info_imp(10)
real(kind(1.d0)), intent(inout)  :: xf(ldxf,nv)

external                         :: x2f_cirxn, x2f_setmethod
logical, save                    :: init=.false.
integer, save                    :: dist_strat(3), atmps(3), k_pos, gat, &
                                    n, n_spec1, n_spec2, nolocal, ntmp, op_frq, &
                                    adap, tstg1, nolog, lu_op, lu_log, lu_time
real, save                       :: r, A_max, A_max_entry
logical, save                    :: msstg = .false., tabfd = .false., localstg = .true.

integer, allocatable             :: info(:,:)
real(kind(1.d0)), allocatable    :: rinfo(:,:), xf_nolocal(:,:)
real(kind(1.d0))                 :: rinfo_xf(1), stats(100), xf_temp(ldxf,nv)
integer                          :: info_xf(1)
integer                          :: i, ierr, nproc, att_max
integer                          :: grp_h, nv_nolocal, nv_un, nv_un_PR1
integer                          :: wall_e, wall_b, wall_rate, wall_max
real                             :: t_i, a_star, a_0, n_q, n_r, n_g, n_a, n_d, wall_R, wall_D
integer, save                    :: next_op, mpicomm_uran
integer, save                    :: time_steps = 0

r = 1.0

if ( .not.impinit ) then
   nolocal = info_imp(1)
   att_max = info_imp(2)
   adap    = info_imp(3)
   set_pt  = info_imp(4)
   tstg1   = info_imp(5)
   nolog   = info_imp(10)
   a_nolog = nolog ! set adaptive module's no logging variable
   ! compute rank
   call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
   ! initialize the log files
   call fileinit( lu_op, lu_log, lu_time )
   ! partition the communicator as specified by set_pt, added by SRL 2/2011
   call x2f_mpi_ptinit( set_pt, pt_index, mpicomm_pt, lu_op )
   impinit = .true.
end if

if (mode == -1) return

time_steps = time_steps + 1

if ( adap==0 ) then
   ! Initialize pre-set strategies
   if ( .not.init .or. localstg ) then
      if ( mode==1 ) then
         dist_strat(1:3) = (/0, 0, 0/)
         atmps(1:3)      = (/0, 0, 1/)
      elseif ( mode==2 ) then
         dist_strat(1:3) = (/0, 0, 1/)
         atmps(1:3)      = (/0, 0, 1/)
      elseif ( mode==3 ) then
         call x2f_mpi_nproc( nproc )
         ntmp = nproc-nolocal
         dist_strat(1:3) = (/0, 4, 1/)
         atmps(1:3)      = (/0, ntmp, 1/)
      elseif ( mode==4 ) then
         if ( att_max<=0 ) att_max = att_max_deft
         call x2f_mpi_nproc( nproc )
         if ( nproc-nolocal>att_max ) then
            ntmp = att_max
         else
            ntmp = nproc-nolocal
         end if
         dist_strat(1:3) = (/0, 4, 1/)
         atmps(1:3)      = (/0, ntmp, 1/)
      elseif ( mode==5 ) then
         call x2f_mpi_nproc( nproc )
         ntmp = nproc-nolocal
         dist_strat(1:3) = (/4, 4, 1/)
         atmps(1:3)      = (/ntmp, min(nSR2,ntmp), 1/)
      ! new code added here for modes 6 and 7, 2/2011
      elseif ( mode==6 ) then
         call iflocalstg( localstg, mpicomm_pt, tstg1, lu_op )
         if ( localstg ) then
            dist_strat(1:3) = (/0, 0, 0/)
            atmps(1:3)      = (/0, 0, 1/)
         else
            dist_strat(1:3) = (/0, 0, 1/)
            atmps(1:3)      = (/0, 0, 1/)
         end if
      elseif ( mode==7 ) then
         print *, 'mode =', mode, 'has not been implemented yet!'
         stop
      elseif ( mode>7 .or. mode<1 ) then
         print *, 'mode =', mode, 'has not been implemented yet!'
         stop
      end if

      k_pos   = nx - 2
      op_frq  = -1
      gat     = sum( atmps )
      n       = size( atmps, 1 )
      n_spec1 = 10
      n_spec2 = 50

      if ( .not.init .and. rank == 0 ) then
         !write the mode, distribution strategies and number of attempts to file
         write(lu_op,'(A,10I6)') 'mode=', mode
         write(lu_op,'(A,10I6)') 'nolocal', nolocal
         write(lu_op,'(A,10I6)') 'distribution strategies:', dist_strat
         write(lu_op,'(A,10I6)') 'number of attempts     :', atmps
         if(mode == 6) write(lu_op,'(A,10I6)') 'P-URAN: time stage 1 (PLP), min. = ', tstg1
      endif
      init = .true.

      if( mode==6 .and. rank == 0 .and. .not.localstg) then
         write(lu_op,'(A)') 'P-URAN: finished PLP stage, switching strategy to URAN'
         write(lu_op,'(A,10I6)') 'P-URAN: numer of time steps in PLP = ', time_steps
         write(lu_op,'(A,10I6)') 'distribution strategies:', dist_strat
         write(lu_op,'(A,10I6)') 'number of attempts     :', atmps
      endif

      if(mode < 6) localstg = .false.
   end if

   !specify info, rinfo, info_xf and rinfo_xf arrays
   allocate(  info(gat,n_spec1) )
   allocate( rinfo(gat,n_spec2) )

   call construct_ary( dist_strat, atmps, n, gat, op_frq, n_spec1, n_spec2, &
                       nolocal, info, rinfo, info_xf, rinfo_xf )
   !added the following call related to partitioning, SRL 2/2011
   call x2f_mpi_setcomms( mpicomm_pt, mpicomm_pt, gat, n_spec1, info )

   call x2f_mpi( nx, nf, ldxf, nv, xf, k_pos, info, rinfo, info_xf, rinfo_xf, &
                 gat, n_spec1, n_spec2, x2f_cirxn, x2f_setmethod, stats )

   deallocate( info  )
   deallocate( rinfo )

   if( nolog < 2 ) then
      ! VH - 03/01/11 - write stats for each attempt
      if ( istep==1 ) write(lu_time,'(A,I8)') '%istep, wall_R, wall_D, cpu_R, cpu_D, [wall_igat, cpu_igat], gat = ', gat
      write(lu_time,'(I8,100E15.6)') istep, stats(50), stats(51), stats(40), stats(41), stats(52:52+gat-2), stats(42:42+gat-2)
   endif

else
   !initialize the adaptive strategy
   if ( .not.pairinit ) then
      !initialize communicator for each partition
      !this is now done above for all modes, not just adaptive, SRL 2/2011
      !call x2f_mpi_ptinit( pt_index, mpicomm_pt )

      !initialize grp_index for each processor
      call x2f_mpi_grpinit( mpicomm_pt, grp_index, g, ng )

      !determine the number of stage needed
      call pairstg_init( mpicomm_pt, ips, mps )

      !measure communication time: t_C
      call x2f_mpi_tc( nx, nf, ldxf, nv, mpicomm_pt, grp_index, t_C, lu_log )
      if( nolog == 0 ) write(lu_log,'(A,E15.6,I8)') 't_C', t_C

      !return the number of table entries
      call ci_stats( 13, A_max_entry )
      pairinit = .true.
   end if
   !check if the current stage is still pairing stage
   call ifpairstg( ips, mps, msstg, pairstg, mpicomm_pt, lu_log )

  !pairstg stays .TRUE. until the final pairing stage is over, which occurs when
  !the average table size in any one group exceeds a predefined max table size;
  !then pairstg tips over to .FALSE. and stays that way for the rest of the run
   if ( pairstg ) then
      !save xf array in xf_temp for measuring Pij matrix
      xf_temp = xf

      !set method and number of attempts based on g
      if ( g==1 ) then
         dist_strat = (/ 0,   0,  0 /)
         atmps      = (/ 1,   0,  1 /)
	  else
         dist_strat = (/ 4,   0,  1 /)
         atmps      = (/ g,   0,  1 /)
      end if

      k_pos   = nx - 2
      op_frq  = -1
      n_spec1 = 10
      n_spec2 = 50
      gat     = sum( atmps )
      n       = size( atmps,1 )

      !specify info, rinfo, info_xf and rinfo_xf arrays
      allocate(  info(gat,n_spec1)  )
      allocate( rinfo(gat,n_spec2)  )

      call construct_ary( dist_strat, atmps, n, gat, op_frq, n_spec1, n_spec2, &
                          0, info, rinfo, info_xf, rinfo_xf )

      !set table entry limit a_star
      call x2f_mpi_setTabEntry( msstg, ips, mps, A_max_entry, a_star, a_0 )
      if( nolog == 0 ) write(lu_log,'(A,5I8)') 'A_max, a_i_limit', int(A_max_entry), int(a_star), int(a_0)

      !get intra-communicator
      call x2f_mpi_crtintra( grp_index, mpicomm_pt, intracomm )

      !set communicator for each attempt, and store it in info array
      if ( .not.tabfd ) mpicomm_uran = intracomm
      call x2f_mpi_setcomms( intracomm, mpicomm_uran, gat, n_spec1, info )

      call isat_stats( 1, n_q, n_r, n_g, n_a, n_d, lu_log )

      !evaluate particles' composition after reaction, and measure the time spent
      call system_clock( count = wall_b)
      call x2f_mpi( nx, nf, ldxf, nv, xf, k_pos, info, rinfo, info_xf, rinfo_xf, &
                    gat, n_spec1, n_spec2, x2f_cirxn, x2f_setmethod, stats )
      call system_clock( count=wall_e, count_rate = wall_rate, count_max = wall_max )
      if ( wall_e<wall_b ) wall_e = wall_e + wall_max
      t_i = float( wall_e-wall_b )/float( wall_rate )
      if( nolog == 0 ) write(lu_log,'(A,I8,10E15.6)') 't_i', istep, t_i

      if( nolog < 2 ) then
         if ( istep==1 ) write(lu_time,'(A)') '%istep, wall_R, wall_D, cpu_R, cpu_D'
         write(lu_time,'(I8,100E15.6)') istep, stats(50), stats(51), stats(40), stats(41)
      endif

      deallocate( info )
      deallocate( rinfo )

      !initialize arrays storing time statistics
      call x2f_mpi_allocate_adapt( ng )

      call ifmsstg( mpicomm_pt, intracomm, g, msstg, lu_log )
      if(nolog==0) write(lu_log,'(A,2I8,3L4)') 'istep', istep, ips, pairstg, msstg, tabfd

      if ( .not.tabfd ) then
         call x2f_mpi_tabfd( mpicomm_pt, lu_log, tabfd )
         if ( tabfd ) then
            msstg = .true.
            mpicomm_uran = mpicomm_pt
         end if
      end if
      if(nolog==0) write(lu_log,'(A,L4)') 'tabfd', tabfd

      !gather all statistics from each group
      if ( msstg ) then
         nv_un  = int(stats(34))
         nv_un_PR1 = int(stats(62))
         wall_R = stats(50)
         wall_D = stats(21)
         call x2f_mpi_grpgather( g, intracomm, ng, mpicomm_pt, wall_R, wall_D, a_0, nv_un, nv_un_PR1, nv, lu_log )
         !measure Pij matrix: start
         if ( g==1 ) then
            dist_strat = (/ 0,   0,  0 /)
            atmps      = (/ 1,   0,  0 /)
	     else
            dist_strat = (/ 4,   0,  0 /)
            atmps      = (/ g,   0,  0 /)
         end if
         gat        = sum( atmps )

         !specify info, rinfo, info_xf and rinfo_xf arrays
         allocate(  info(gat,n_spec1)  )
         allocate( rinfo(gat,n_spec2)  )

         op_frq = -1
	     call construct_ary( dist_strat, atmps, n, gat, op_frq, n_spec1, n_spec2, &
                             0, info, rinfo, info_xf, rinfo_xf )

         !set communicator for each attempt, and store it in info array
	     call x2f_mpi_setcomms( intracomm, intracomm, gat, n_spec1, info )

         !in the following loop, pair two groups (grp_index and grp_h), exchange
         !the particles on two groups, and measure the probability of a particle
         !in group grp_index being retrievable from other group grp_h
         do i = 1, ng
            !for group grp_index, pick another group grp_h to form a pair
            call x2f_mpi_pair( i, ng, grp_index, grp_h )

            if ( grp_index.ne.grp_h ) then
               !create intercommunicator between group grp_index and group grp_h
               call x2f_mpi_crtinter( grp_index, grp_h, ng, intracomm, mpicomm_pt, &
                                      intercomm )
               !get the number of particles on the other group
               call x2f_mpi_getxfsize( nv, intracomm, intercomm, nv_nolocal )

               allocate( xf_nolocal(ldxf,nv_nolocal) )
               !exchange the particles
               call x2f_mpi_xfswap( nv, ldxf, xf_temp, intracomm, intercomm, &
                                    nv_nolocal, xf_nolocal )
               call x2f_mpi_freecomm( intercomm )
            elseif ( grp_index.eq.grp_h ) then
               nv_nolocal = nv
               allocate( xf_nolocal(ldxf,nv_nolocal) )
               xf_nolocal = xf_temp
            end if

            !try to retrieve particles' compositions after reaction step
            call system_clock( count = wall_b)
            call x2f_mpi( nx, nf, ldxf, nv_nolocal, xf_nolocal, k_pos, info, &
                          rinfo, info_xf, rinfo_xf, gat, n_spec1, n_spec2, &
                          x2f_cirxn, x2f_setmethod, stats )
            call system_clock( count=wall_e, count_rate = wall_rate, count_max = wall_max )
            if ( wall_e<wall_b ) wall_e = wall_e + wall_max
            t_i = float( wall_e-wall_b )/float( wall_rate )

            deallocate( xf_nolocal )

            nv_un = int(stats(33))
            if(nolog==0) write(lu_log,'(A,2I8,10E15.6)')  'nv_un', grp_h, nv_un
            !measure the quantity "Pij(grp_index, grp_h)"
            call x2f_mpi_Pij( nv_un, nv_nolocal, mpicomm_pt, intracomm, i )
            call x2f_mpi_tRij( t_i, nv_nolocal, mpicomm_pt, intracomm, i, g )
         end do

         deallocate( info )
         deallocate( rinfo )

         !gather Pij into matrix P_ij_mtx
         call x2f_mpi_Pij_mtx( ng, lu_log )
         call x2f_mpi_tRij_mtx( ng, lu_log )

         call x2f_mpi_newmatch( mpicomm_pt, ips, ng, g, r, A_max_entry, t_C, lu_log )

         call x2f_mpi_setcomm( mpicomm_pt, ng, g, match, grp_index, istep )
         if(nolog==0) write( lu_log, '(A, 3I6)' ) 'grp_index', grp_index, ng, g

      end if

      call x2f_mpi_freecomm( intracomm )

   else
      if ( .not.init ) then
         if ( g==1 ) then
            dist_strat = (/ 0,   0,  1 /)
            atmps      = (/ 1,   0,  1 /)
	     else
           !dist_strat = (/ 4,   0,  1 /)
           !atmps      = (/ g,   0,  1 /)
           dist_strat = (/ 4,   4,  1 /)
           atmps      = (/ g-nolocal, min(g,nSR2), 1 /)

           call x2f_mpi_crtintra( grp_index, mpicomm_pt, intracomm )
         end if

         gat = sum( atmps )
         init = .true.
         if( rank == 0 ) then
            write(lu_op,'(A,10I4)') 'adaptive init', g, ng, grp_index
            write(lu_op,'(A,10I6)') 'nolocal', nolocal
            write(lu_op,'(A,10I6)') 'distribution strategies:', dist_strat
            write(lu_op,'(A,10I6)') 'number of attempts     :', atmps
         endif
         next_op = istep
      end if
      allocate(  info(gat,n_spec1)  )
      allocate( rinfo(gat,n_spec2)  )
      op_frq = -1
      call construct_ary( dist_strat, atmps, n, gat, op_frq, n_spec1, n_spec2, &
                          nolocal, info, rinfo, info_xf, rinfo_xf )

      call x2f_mpi_setcomms( intracomm, mpicomm_pt, gat, n_spec1, info )

      call x2f_mpi( nx, nf, ldxf, nv, xf, k_pos, info, rinfo, info_xf, rinfo_xf, &
                    gat, n_spec1, n_spec2, x2f_cirxn, x2f_setmethod, stats )

      deallocate(  info  )
      deallocate( rinfo  )

      if(nolog < 2) write(lu_time,'(I8,100E15.6)') istep, stats(50), stats(51), stats(40), stats(41)

   end if
end if

end subroutine ISAT_MP
