!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the x2f_mpi      !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine x2f_mpi( nx, nf, ldxf, nv, xf, k_pos, info, rinfo, info_xf, rinfo_xf, &
                    gat, n_spec1, n_spec2, x2f, x2f_setmethod, stats )
!
! Overview:  there is a user-supplied subroutine x2f which evaluates a function f(x),
!     where x is an nx-vector, and f is an nf-vector. The calling sequence of x2f is:
!        call x2f( qt_mode, nx, nf, x, f, info_xf, rinfo_xf )
!     where info_xf and rinfo_xf can be used to communicate additional (x-independent)
!     data, and qt_mode concerns "quick-try mode" and is explained below. 
!
!     Given nv x-vectors stored in an array xf, x2f_mpi returns the corresponding f
!     for all x-vectors (as evaluated by x2f), replacing the original x-vectors in xf.
!     When used in an MPI environment, various message-passing strategies can be used
!     to balance the load of the evaluation of f for all values of x across processors,
!     because the amount of work required to compute f might vary strongly with x.
!
! Input:
!   nx         - length of x vectors
!   nf         - length of f vectors
!   ldxf       - the leading dimension of the array xf ( ldxf >= max(nx,nf) )
!   nv         - number of vectors of x
!   xf         - array which (on input) contains in xf(1:nx,1:nv) the vectors of x
!   k_pos      - xf(k_pos,i) indicates if the i-th particle has been evaluated or not
!                >  0, evaluated
!                <= 0, not evaluated
!   info       - integer array of parameters controlling x2f_mpi operation
!   rinfo      - real array of parameters controlling x2f_mpi operation 
!   info_xf    - integer array passed to x2f
!   rinfo_xf   - real array passed to x2f
!   x2f        - the name of the user-supplied subroutine
!x2f_setmethod - the name of the user-supplied subroutine; to activate the specified 
!                method so that it can be used
!
! Output:
!   xf        - on output, the values of f are contained in xf(1:nf,1:nv)
!   stats     - real array containing information about the performance of x2f_mpi
!
! info        - integer array of parameters controlling x2f_mpi operation
!   i         - the i-th global attempt 
!   info(i,1) - indicate the method to be used in the i-th global attempt
!   info(i,2) = mode_p:  mode of parallel implementation
!             = 0 - purely local processing (PLP) - coded, no message passing needed
!             = 1 - uniform random distribution (URAN) - coded with MPI
!             = 2 - local (in rank), level distribution (LOCLEV) - coded with MPI
!             = 3 - balancing by input stratification (STRATA) - coded with MPI
!                  * REQUIRES info(i,3), info(i,4) and possibly rinfo(11:nx+10) set *
!             = 4 - preference to some processor(PREF) - coded with MPI
!                  * REQUIRES info(i,3) set *
!   info(i,3) - parameter for the specific distribution strategy
!               'STRATA': ( when info(i,2)==3 )
!                     u(info(i,3))=1, where u(1:nx) is the unit vector along which the
!                     projection is taken
!               'PREF': ( when info(i,2)==4 )
!                     0 - local processor needs to be tried
!                     1 - local processor has been tried
!   info(i,4) - mbin for STRATA: number of bins assigned to each processor
!   info(i,5) - MPI communicator
!   info(i,6) = 1 (re)initialize: call the "init" routine appropriate for mode_p,
!                 will (re)initialize the criteria by which x-vectors are distributed
!   info(i,7) = nn >  0, output from x2f_mpi every nn steps
!               nn <= 0, no output
!   
! rinfo:
! rinfo(i,11:nx+10) = used for STRATA ( applied when info(i,2)==3,info(i,3)==0 )
!
! stats:
!  stats(1:12) : CPU times for each attempt
!  stats(1)  = CPU time for this attempt
!  stats(2)  = CPU time for evaluating non-local particles
!  stats(3)  = CPU time for evaluating local particles
!  stats(4)  = CPU time spend in x2f_***_p
!  stats(5)  = CPU time spend in sorting particles 1st place
!  stats(6)  = CPU time spend in sorting particles 2nd place
!  stats(7)  = CPU time spend in MPI_Alltoallv communication 1st place
!  stats(8)  = CPU time spend in MPI_Alltoallv communication 2nd place
!  stats(9)  = CPU time spend in MPI_Waitany 1st place
!  stats(10) = CPU time spend in MPI_Waitany 2nd place
!  stats(11) = CPU time spend in MPI_Waitany 3rd place  
!  stats(12) = CPU time spend in retrieve attempts
!
!  stats(21:29) : Wall-clock time for each attempt
!  stats(21) = Wall-clock time for this attempt
!  stats(22) = Wall-clock time for evaluating non-local particles
!  stats(23) = Wall-clock time for evaluating local particles
!  stats(24) = Wall-clock time spend in x2f_***_p
!  stats(25) = Wall-clock time spend in sorting particles 1st place
!  stats(26) = Wall-clock time spend in sorting particles 2nd place 
!  stats(27) = Wall-clock time spend in MPI_Alltoallv communication 1st place
!  stats(28) = Wall-clock time spend in MPI_Alltoallv communication 2nd place
!  stats(29) = Wall-clock time spend in retrieve attempts
!
!  stats(31) = number of particle stayed local
!  stats(32) = number of particle from other processor
!  stats(33) = number of particles on the processor still needs to be evaluated
!              after this attempt
!  stats(34) = number of particles on the processor still needs to be evaluated
!              on this attempt
!
!  stats(40)      = CPU time spent in retrieve attempts
!  stats(41)      = CPU time spent in DE attempt
!  stats(41+igat) = CPU time spent in "igat" retrieve attempt
!  stats(50)      = Wall-clock time spent in retrieve attempts
!  stats(51)      = Wall-clock time spent in DE attempt
!  stats(51+igat) = Wall-clock time spent in "igat" retrieve attempt
!  stats(60+igat) = number of particles on the processor still needs to be evaluated
!  stats(70+igat) = number of particles from other processor 
!
! Notes on timing:
! The quantities wall_R, wall_D, cpu_R, cpu_D are written into log files, time_N.log.
! Retrieves include all primary and secondary retrieves that are attempted by rank N.
! The retrieves may have been successful or not; all attempts contribute to the totals.
! Time spent in communication while particles are being shuttled around is included, too.
! Attempted particles may be nonresident, i.e., they may have come from other processes.
!
! The quantity (wall minus cpu) = (system plus idle), rather than the communcation time.
! The above becomes invalid under multithreading, where we can easily have cpu > wall.
! System time is likely I/O, as the fast implementations of MPI operate in user space.
! (The special mpd's or Multi-Purpose Daemons are just used for MPI process initiation.)
! Idle time doesn't result from MPI blocking states: waiting uses cpu time by polling!
!
use MPI
use x2f_uran
use x2f_loclev
use x2f_strata
use x2f_pref
use sort_tools

implicit none
integer, parameter         :: k_dp = kind(1.d0)

external                   :: x2f, x2f_setmethod
integer,    intent(in)     :: nx, nf, nv, ldxf, k_pos, gat, n_spec1, n_spec2
integer,    intent(in)     :: info(gat,n_spec1), info_xf(*)
real(k_dp), intent(in)     :: rinfo(gat,n_spec2), rinfo_xf(*)
real(k_dp), intent(inout)  :: xf(ldxf,nv)
real(k_dp), intent(out)    :: stats(100)

logical, save              :: initialized = .false., pt2pt = .true.
integer, save              :: mode_p, myrank, nproc, mpicomm
type (strata_type), save   :: strat
type (pref_type), save     :: pref

integer, allocatable, save :: n_nonres(:), n_distrib(:), n_nonres_recv(:), &
                              n_distrib_recv(:)
integer, allocatable, save :: start_nr_buf(:), start_d_buf(:)
integer, allocatable, save :: start_nr_buf1(:), start_d_buf1(:)
integer, allocatable, save :: req_d(:), req_nr1(:), req_nr2(:)
real(k_dp), allocatable    :: xf_nr_buf(:,:)
!real(k_dp)                 :: f(nf)
real(k_dp)                 :: f(ldxf)

integer    :: p(nv), sourceof(nv), sortedto(nv), bucket_total(1), sourceof1(nv)
integer    :: start_marker, start_marker1, count_nr, count_d, startat, stopat
integer    :: mode_u, mbin, ict, index, ierr, istatus(MPI_STATUS_SIZE)

integer    :: nv_un, ntogo, nlp, nout, nnr, igat, method, ipar, j, i, pref_at
integer    :: wall0, wall, wall_rate, wall_max, wall_b, wall_e
real       :: cpu0, cpu, cpu_b, cpu_e, cpu_nr, cpu_nr0, uset(nx)

logical, save      :: fileinit0 = .false.
integer, save      :: istep = 1, unit2 = 15
character(30),save :: blank, head, tail, times_op2
integer            :: count

nv_un = nv
do i = 1, nv_un
   sourceof(i)  = i
   sourceof1(i) = i
end do

! Loop over all attempts

do igat = 1, gat
   if (igat==1) stats  = 0.d0
   ! Set method for i-th attempt
   if (igat==1) then
      method = info(igat,1)
      call x2f_setmethod(method,info_xf,rinfo_xf)     
   else
      if ( info(igat,1)/=info(igat-1,1) ) then
         method = info(igat,1)
         call x2f_setmethod(method,info_xf,rinfo_xf)    
      end if
   end if
   
   if ( .not.fileinit0.and.igat==1 ) then
      ! MPI global rank of this processor
      call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr )
            
      if ( info(igat,7)>0 ) then       
         blank = repeat(' ',30)  
         tail  = blank
         tail  = 'out'
      
         head = blank
         head = 'x2fstats'
         call x2fmpi_file_name(head, myrank, -1, tail, times_op2)
         open( unit2,file=times_op2,status='replace',action = 'write' )
         fileinit0 = .true.
      end if
   end if    

   ! Set distribution strategy for i-th attempt and initialize it
   if ( (.not.initialized ).or.(info(igat,2)/= mode_p).or.(info(igat,6)==1) ) then
      mode_p = info(igat,2)
      if ( mode_p < 0 .or. mode_p > 4 ) then
         write(0,*) 'x2f_mpi: not yet coded for mode_p = ', mode_p
         return
      elseif( mode_p == 0 ) then        
         myrank = 0
         nproc  = 1
      else  ! get myrank and nproc from MPI and initialize load balancing method 
         if ( info(igat,5) == 0 ) then
            mpicomm = MPI_COMM_WORLD  ! default value
         else
            mpicomm = info(igat,5)  ! bad value of mpicomm will throw an exception in MPI
         end if
         
         call MPI_Comm_rank( mpicomm, myrank, ierr )  ! MPI rank of this processor
         call MPI_Comm_size( mpicomm, nproc, ierr )   ! number of processors   
         
         if ( mode_p == 1 ) then
            call x2f_uran_init( myrank )
         elseif ( mode_p == 2 ) then
            call x2f_loclev_init
         elseif ( mode_p == 3 ) then
            mode_u = info(igat,3)
            if ( mode_u == 0 ) uset(1:nx) = rinfo(igat,11:nx+10) 
            mbin   = info(igat,4)
            call x2f_strata_init(nx,mbin,mpicomm,myrank,nproc,mode_u,uset,strat)
         elseif( mode_p == 4 ) then
            pref_at = 0
            call x2f_check_atmpts(gat,n_spec1,info,igat,nproc)
            call x2f_pref_init( pref, myrank, nv_un, info(igat,3))
         end if
      end if

      ! allocate several constant-size, persistent arrays for later use
      if( allocated( n_nonres      ) ) deallocate( n_nonres      )
      if( allocated( n_nonres_recv ) ) deallocate( n_nonres_recv )
      if( allocated( n_distrib     ) ) deallocate( n_distrib     )
      if( allocated( n_distrib_recv) ) deallocate( n_distrib_recv)
      if( allocated( start_nr_buf  ) ) deallocate( start_nr_buf  )
      if( allocated( start_d_buf   ) ) deallocate( start_d_buf   )
      if( allocated( start_nr_buf1 ) ) deallocate( start_nr_buf1 )
      if( allocated( start_d_buf1  ) ) deallocate( start_d_buf1  )
      if( allocated( req_d         ) ) deallocate( req_d         )
      if( allocated( req_nr1       ) ) deallocate( req_nr1       )
      if( allocated( req_nr2       ) ) deallocate( req_nr2       )
      
      allocate( n_nonres      (nproc) )
      allocate( n_nonres_recv (nproc) )
      allocate( n_distrib     (nproc) )
      allocate( n_distrib_recv(nproc) )
      allocate( start_nr_buf  (nproc) )
      allocate( start_d_buf   (nproc) )
      allocate( start_nr_buf1 (nproc) )
      allocate( start_d_buf1  (nproc) )
      allocate( req_d         (nproc) )
      allocate( req_nr1       (nproc) )
      allocate( req_nr2       (nproc) )

      initialized = .true.
   end if
   
   pt2pt = .true.
   !if ( mode_p==1 ) pt2pt = .false.

   !Call barrier before get wall clock time
   ! -- VH 03/04/2011 - comment out BARRIER calls, not needed --
   !if (mode_p/=0) then
   !   call MPI_BARRIER(mpicomm,ierr)
   !end if
   call cpu_time( cpu0 )  
   call system_clock( count = wall0 )   

   p(1:nv_un) = -1
   ntogo      = nv_un
   stats(34)  = nv_un
   if (igat<9) stats(60+igat) = stats(34)

   !write(unit2,*) 'here0', igat
   !=========== phase 1: distribute x vectors to other processors ===================
   call system_clock( count = wall_b )
   call cpu_time( cpu_b )
   if ( mode_p == 0 ) then
      nlp        = nv_un       ! number of vectors processed locally (all)
      nout       = 0           ! number of vectors processed elsewhere (none)
      nnr        = 0           ! number of nonresident vectors to process (none)
      p(1:nv_un) = -p(1:nv_un) ! causes all remaining vectors to be processed locally
      n_distrib(1) = nv_un
      n_nonres(1)  = nv_un
         
   elseif ( mode_p == 1 ) then
      call x2f_uran_p( nproc, mpicomm, myrank, ntogo, nv_un, p, n_nonres, n_distrib )

   elseif ( mode_p == 2 ) then
      call x2f_loclev_p( nproc, mpicomm, myrank, ntogo, nv_un, p, n_nonres, n_distrib )

   elseif ( mode_p == 3 ) then
      call x2f_strata_p( strat, nx, ldxf, xf, ntogo, nv_un, p, n_nonres, n_distrib )
   
   elseif ( mode_p == 4 ) then
      pref_at = pref_at + 1
      call x2f_pref_p( nproc, myrank, mpicomm, pref_at, pref, p, nv_un, n_nonres, &
                       n_distrib )
      do i = 1, nv_un
         if ( p(i)<=0.or.p(i)>nproc ) stop 'p(i) error'
      end do
   end if
   !write(unit2,*) 'here01', igat
   call cpu_time( cpu_e )
   call system_clock( count=wall_e, count_rate=wall_rate, count_max=wall_max )
   stats(4)  = cpu_e-cpu_b
   if ( wall_e<wall_b ) wall_e = wall_e + wall_max
   stats(24) = float( wall_e-wall_b )/float( wall_rate )

   ! the displacement used for MPI_Alltoallv
   n_nonres_recv(1:nproc)   = n_nonres(1:nproc)
   n_nonres_recv(myrank+1)  = 0
   n_distrib_recv(1:nproc)  = n_distrib(1:nproc)
   n_distrib_recv(myrank+1) = 0

   ! Create an input buffer big enough to hold all the incoming x vectors
   nnr = sum( n_nonres ) - n_nonres(myrank+1)
   allocate( xf_nr_buf(ldxf,nnr) )

   ! no need to create an outgoing buffer; xf will be sorted for that purpose
   nout = sum( n_distrib ) - n_distrib(myrank+1)
   nlp  = nv_un - nout
   stats(31) = nlp
   stats(32) = nnr
   
   if (igat<9) stats(70+igat) = stats(32)
   if (igat<9) stats(80+igat) = stats(31)
   
!-------------------------------------------------------------------
! ...if storage gets tight, here's what to try in the future:
! - make the incoming buffer smaller and reuse it for each message
!   i.e., try nin = maxval( n_nonres )  ! smaller than nnr
! - post all nonblocking sends and receives related to xf
! - start local processing, periodically interrupting with MPI_Testany
! - if a batch arrives, process it immediately and send it back
!-------------------------------------------------------------------

   ! in the following arrays, non-participating processors are flagged with 0
   ! this array marks the starting point of each message in xf_nr_buf
   start_nr_buf(1:nproc)  = 0
   start_nr_buf1(1:nproc) = 0

   ! this array marks the starting point of each message distributed from xf
   start_d_buf(1:nproc)  = 0
   start_d_buf1(1:nproc) = 0

   ! arrays in the following group hold request handles for nonblocking calls
   !req_d(1:nproc)   = 0
   !req_nr1(1:nproc) = 0
   !req_nr2(1:nproc) = 0
   req_d(1:nproc)   = MPI_REQUEST_NULL
   req_nr1(1:nproc) = MPI_REQUEST_NULL
   req_nr2(1:nproc) = MPI_REQUEST_NULL
 
   ! start communication: post nonblocking receives for all expected messages
   start_marker  = 1
   start_marker1 = 0
   count_nr      = 0
   do j = 1, nproc
      if ( pt2pt ) then
         if ( n_nonres_recv(j)/=0 ) then
            call MPI_Irecv( xf_nr_buf(1,start_marker), ldxf*n_nonres(j), &
                            MPI_DOUBLE_PRECISION, j-1, 200+j, mpicomm, req_nr1(j), ierr )
            start_nr_buf(j) = start_marker
            start_marker = start_marker + n_nonres(j)
            count_nr = count_nr + 1
         end if
      else
         start_nr_buf(j)  = start_marker
         start_marker     = start_marker + n_nonres_recv(j)
         if ( n_nonres_recv(j)/=0 ) count_nr = count_nr + 1
         start_nr_buf1(j) = start_marker1
         start_marker1    = start_marker1 + ldxf*n_nonres_recv(j)
      end if
   end do

   ! paint some boundaries in xf; portions of it will be message buffers
   start_marker  = 1
   start_marker1 = 0
   count_d       = 0
   do j = 1, nproc
      if ( pt2pt ) then
         if ( n_distrib(j) /= 0 ) then
            start_d_buf(j) = start_marker
            start_marker = start_marker + n_distrib(j)
            if ( j /= myrank + 1 ) count_d = count_d + 1
         end if
      else
         start_d_buf(j)   = start_marker
         start_marker     = start_marker + n_distrib(j)
         start_d_buf1(j)  = start_marker1
         start_marker1    = start_marker1 + ldxf*n_distrib(j)
      end if
   end do
    
! rearrange data to make an output buffer in xf, using a bucket sort algorithm
! the number of buckets equals nproc plus 1 for those already done by quick try
! sourceof is an array to link sorted locations in xf back to source locations
! for nout < i <= nv, sourceof(i) identifies a source for local processing
! note, each vector to be fed as input to the x2f function will be of length nx
   
   call system_clock( count=wall_b )
   call cpu_time( cpu_b )
   call bucket_sort( nv_un, ldxf, nx, xf, p, nproc, n_distrib, sourceof, sourceof1 )

   if ( mode_p==4 ) then
      p(1:nv_un) = pref%b(1:nv_un)
      do j = 1, nv_un
         pref%b(j) = p(sourceof1(j))
      end do
   end if
   
   call cpu_time( cpu_e )
   call system_clock( count=wall_e, count_rate=wall_rate, count_max=wall_max )
   stats(5) = cpu_e-cpu_b
   if ( wall_e<wall_b ) wall_e = wall_e + wall_max
   stats(25) = float( wall_e-wall_b )/float( wall_rate )

   ! buffer is set: post nonblocking sends AND receives for all outgoing messages
   do j = 1, nproc
      if ( ( start_d_buf(j) /= 0 ) .and. ( j /= myrank + 1 ) .and. pt2pt) then
         call MPI_Isend( xf(1,start_d_buf(j)), ldxf*n_distrib(j), &
              MPI_DOUBLE_PRECISION, j-1, 201+myrank, mpicomm, req_d(j), ierr )
         call MPI_Request_free( req_d(j), ierr )
         call MPI_Irecv( xf(1,start_d_buf(j)), ldxf*n_distrib(j), &
              MPI_DOUBLE_PRECISION, j-1, 300+j, mpicomm, req_d(j), ierr )
      end if
   end do 
   
   if ( .not.pt2pt ) then
      call system_clock( count = wall_b )
      call cpu_time( cpu_b )
      call MPI_Alltoallv( xf,ldxf*n_distrib_recv,start_d_buf1,MPI_DOUBLE_PRECISION, &
         xf_nr_buf,ldxf*n_nonres_recv,start_nr_buf1,MPI_DOUBLE_PRECISION,mpicomm,ierr )
      call cpu_time( cpu_e )
      call system_clock( count=wall_e, count_rate = wall_rate, count_max = wall_max )
      if ( wall_e<wall_b ) wall_e = wall_e+wall_max
      stats(7)  = cpu_e-cpu_b
      stats(27) = float( wall_e-wall_b )/float( wall_rate )
   end if
 
   !write(unit2,*) 'here1', igat
   !=========== phase 3: local processing of resident vectors ========================
   ! identify chunk of sorted xf array holding vectors for local processing
   startat = start_d_buf(myrank+1)
   stopat  = start_d_buf(myrank+1) + n_distrib(myrank+1) - 1
   ! (SRL) insert missing wall timing code
   call system_clock( count = wall_b )  
   call cpu_time( cpu_b )
   do i = startat, stopat   
      call x2f( ldxf, nx, xf(1:nx,i), f(1:ldxf), k_pos, info_xf, rinfo_xf )
      xf(1:ldxf,i) = f(1:ldxf)
   end do

   call cpu_time( cpu_e )
   call system_clock( count=wall_e, count_rate = wall_rate, count_max = wall_max )
   if ( wall_e<wall_b ) wall_e = wall_e + wall_max
   stats(3)  = cpu_e - cpu_b 
   stats(23) = float( wall_e-wall_b )/float( wall_rate )

   !write(unit2,*) 'here2', igat
   !=========== phase 4: processing of non-resident vectors ==========================
   call system_clock( count = wall_b )
   call cpu_time( cpu_nr0 )    
   index = 1
   do ict = 1, count_nr
      if ( pt2pt ) then
         call cpu_time( cpu_b )
         call MPI_Waitany( nproc, req_nr1, index, istatus, ierr )
         call cpu_time( cpu_e )
         stats(9) = stats(9)+cpu_e-cpu_b 
        
         startat = start_nr_buf(index)
         stopat  = start_nr_buf(index) + n_nonres(index) - 1
      else
         if ( n_nonres_recv(index)==0 ) then
            do while ( n_nonres_recv(index)==0 )
               index = index + 1
            end do 
         end if
         startat = start_nr_buf(index)
         stopat  = start_nr_buf(index) + n_nonres_recv(index) - 1
         index = index + 1
      end if

      do i = startat, stopat    
         call x2f(ldxf,nx,xf_nr_buf(1:nx,i),f(1:ldxf),k_pos,info_xf,rinfo_xf)  
         xf_nr_buf(1:ldxf,i) = f(1:ldxf)
      end do
   
      if ( pt2pt ) then
         call MPI_Isend( xf_nr_buf(1,startat), ldxf*n_nonres(index), &
         MPI_DOUBLE_PRECISION, index-1, 301+myrank, mpicomm, req_nr2(index), ierr )
      end if
   end do

   call cpu_time( cpu_nr )
   call system_clock( count=wall_e, count_rate = wall_rate, count_max = wall_max )
   if ( wall_e<wall_b ) wall_e = wall_e + wall_max
   stats(2)  = cpu_nr - cpu_nr0 
   stats(22) = float( wall_e-wall_b )/float( wall_rate )
 
   !=========== phase 5: receive f vectors from other processors =====================
   !write(unit2,*) 'here3', igat
   if ( pt2pt ) then
      call cpu_time( cpu_b )
      do ict = 1, count_d
         call MPI_Waitany( nproc, req_d, index, istatus, ierr )
      end do
      call cpu_time( cpu_e )
      stats(10) = cpu_e-cpu_b
   else
      call system_clock( count = wall_b )
      call cpu_time( cpu_b )
      call MPI_Alltoallv( xf_nr_buf, ldxf*n_nonres_recv, start_nr_buf1, &
           MPI_DOUBLE_PRECISION, xf,ldxf*n_distrib_recv, start_d_buf1, &
           MPI_DOUBLE_PRECISION,mpicomm,ierr )
      call cpu_time( cpu_e )
      call system_clock( count=wall_e, count_rate = wall_rate, count_max = wall_max )
      if ( wall_e<wall_b ) wall_e = wall_e+wall_max
      stats(8)  = cpu_e-cpu_b
      stats(28) = float( wall_e-wall_b )/float( wall_rate )
   end if
 
   call system_clock( count=wall_b )
   call cpu_time( cpu_b )
   ntogo = 0
   if ( mode_p/=4 ) then
      do i = 1, nv_un
         if ( xf(k_pos,i)>0 ) then
            p(i) = 0
         else
            xf(k_pos,i) = -xf(k_pos,i)
            p(i) = 1
            ntogo  = ntogo + 1
         end if
      end do
      bucket_total(1) = ntogo
      call bucket_sort( nv_un, ldxf, ldxf, xf, p, 1, bucket_total, sourceof, sourceof1 )
   else
      do i = 1, sum(n_distrib)
         if ( xf(k_pos,i)>0 ) then
            pref%b(i) = 0
         else
            xf(k_pos,i) = -xf(k_pos,i)
            ntogo  = ntogo + 1
         end if
      end do
      
      call x2f_pref_update( pref, pref_at, nproc, n_distrib, myrank )
      
      ! Sort xf array based on b(i)
      call bucket_sort( nv_un, ldxf, ldxf, xf, pref%b, pref%nbat, pref%npar, sourceof, sourceof1 )
      startat = 1
      do i = 1, pref%nbat
         stopat = startat+pref%npar(i)-1
         pref%b(startat:stopat) = i 
         startat = stopat+1
      end do     
      if ( sum(pref%npar)<nv_un ) pref%b(sum(pref%npar)+1:nv_un)=0      
   end if 

   call cpu_time( cpu_e )
   call system_clock( count=wall_e, count_rate=wall_rate, count_max=wall_max )
   if ( wall_e<wall_b ) wall_e = wall_e + wall_max
   stats(6) = cpu_e - cpu_b
   stats(26) = float( wall_e-wall_b )/float( wall_rate )
   
   nv_un = ntogo + nv_un - sum(n_distrib)
   stats(33) = nv_un
   
   ! mop up the nonblocking calls that sent back computed, nonresident f-vectors
   if ( pt2pt ) then
      call cpu_time( cpu_b )
      do ict = 1, count_nr
         call MPI_Waitany( nproc, req_nr2, index, istatus, ierr )
      end do
      call cpu_time( cpu_e )
      stats(11) = cpu_e-cpu_b
   end if
   call cpu_time( cpu ) 
   stats(1) = cpu - cpu0
   if (igat<9) stats(41+igat) = stats(1)

   ! -- VH 03/04/2011 - comment out BARRIER calls, not needed --
   !if (mode_p/=0) then
   !   call MPI_BARRIER(mpicomm,ierr)
   !end if

   call system_clock( count = wall, count_rate = wall_rate, count_max = wall_max )
   if ( wall<wall0 ) wall = wall + wall_max
   stats(21) = float( wall-wall0 )/float( wall_rate ) 
   if (igat<9) stats(51+igat) = stats(21)
   
   if ( igat<gat ) then
      stats(40) = stats(40) + stats(1)
      stats(50) = stats(50) + stats(21)
   elseif ( igat==gat ) then
      stats(41) = stats(1)
      stats(51) = stats(21)
   end if
   
   deallocate ( xf_nr_buf )
      
   if ( info(igat,7)>0 ) then
      if ( mod(istep,info(igat,7))==0 ) then
         ! (SRL) include format in write statement
         write(unit2,'(100(E16.6))') float(igat),stats(1:33),float(istep)
      end if
   end if

end do ! end of looping over attempts

istep = istep + 1

call pigeonhole_sort( nv, ldxf, nf, xf, sourceof, sortedto)

return
end subroutine x2f_mpi
