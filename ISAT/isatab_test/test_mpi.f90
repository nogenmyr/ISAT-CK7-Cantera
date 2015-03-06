program test_mpi  

!  ISAT  test program which calls the simple driver isatab.

!  S.B. Pope 10/7/2006

!  Table 1: 3D, linear, fixed domain
!  Table 2: 5D, non-linear, fixed domain (100Mbytes), ifull=1
!  Table 3: 4D, non-linear, shifting domain

   use mpi
   use isat_rnu
   implicit none

   integer, parameter :: ntab = 3, nn = 10
   integer, parameter :: nxt(ntab)  = (/ 3, 5, 4 /)
   integer, parameter :: nft(ntab)  = (/ 7, 6, 2 /)
   integer, parameter :: nht(ntab)  = (/ 9, 6, 0 /)
   integer, parameter :: ifsin(ntab)  = (/ 0, 1, 1 /)

   double precision :: rinfo(50,ntab), stats(100), rusr(1), err(ntab)
   double precision :: xmin(nn,ntab), xmax(nn,ntab), dx(nn,ntab), &
                       x0(nn,ntab)
   double precision :: errtol(ntab) = (/ 1.d-3, 1.d-2, 1.d-2 /)
   double precision, pointer :: x(:), f(:), g(:,:), h(:), &
               ft(:), gt(:,:), ht(:)

   double precision :: errmax(3,ntab), errsum(3,ntab), xtab
   double precision :: tstart, tend
   integer :: info(100,3), iusr(2), need(3), itest(ntab), inctes, iq, id
   integer :: nx, nf, nh, nhd
   integer :: nquery = 10000, ntest = 10000
   integer :: myid, ierr, numprocs, world
   integer :: lu = 6
   logical :: if_mpi = .false.

   external fghex

!---  set ranges randomly (but the same for each processor)

   call rnu( xmin )
   call rnu( xmax )

!----  initialization: set x ranges randomly  ----------------------------

   x0   = xmin
   xmax = xmax + .5
   dx   = 0.

!-- specifications for all tables
   info  = 0
   rinfo = 0
   info(2,:)  = 1	! if_g
   info(12,:) = 2	! isatop
   info(28,:) = 1   ! idites
   info(70,:) = 1   ! cdf_err
   rinfo(1,:) = errtol
   rusr       = 1.d0
!-- special settings: Table 2
   info(21,2) = 1   ! ifull
   rinfo(8,2) = 100.  !  stomby

   dx(:,3)    = 1.d-5

   errsum = 0.
   errmax = 0.
   inctes = max( nquery/ntest, 1 )
   itest  = 0

!----  initialize MPI  --------------------------------------------------

      if_mpi = .true.

      call MPI_INIT( ierr )
      tstart = MPI_WTIME()
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      world = MPI_COMM_WORLD

!---- set different random seed on each processor

   call rnused( myid )
   call rnu( tend )

   if( myid == 0 ) write(0,*)'test_mpi started: numprocs = ',numprocs

!----  loop over queries  ------------------------------------------------

   do iq = 1, nquery
   call rnu( xtab )
   id = 1 + ntab*xtab

   nx = nxt(id)
   nf = nft(id)
   nh = nht(id)
   nhd = max( 1, nh )
   iusr(1) = id
   iusr(2) = ifsin(id)

   allocate( x(nx) )
   allocate( f(nf) )
   allocate( g(nf,nx) )
   allocate( h(nhd) )

!--- generate x randomly

   x0(:,id) = x0(:,id) + dx(:,id)
   call rnu( x )
   x = x0(1:nx,id) + xmax(1:nx,id) * x

   call isatab( id, 0, nx, x, nf, nh, nhd, fghex, iusr, rusr, &
                info(:,id), rinfo(:,id), f, g, h, stats )

!---  test for accuracy

   if( mod(iq-1,inctes) .eq. 0 ) then
      need = 1
      allocate( ft(nf) )
      allocate( gt(nf,nx) )
      allocate( ht(nh) )
      call fghex( need, nx, x, nf, nh, iusr, rusr, ft, gt, ht )
      err(1) = sqrt( dot_product( f-ft, f-ft) )
      gt = abs(gt-g)
      err(2) = sqrt( sum( gt * gt ) )
      if( nh == 0 ) then 
         err(3) = 0.
      else
         err(3) = sqrt( dot_product( h-ht, h-ht) )
      endif
      deallocate( ft, gt, ht)

      itest(id)  = itest(id) + 1
      errsum(:,id) = errsum(:,id) + err
      errmax(1,id) = max( err(1), errmax(1,id) )
      errmax(2,id) = max( err(2), errmax(2,id) )
      errmax(3,id) = max( err(3), errmax(3,id) )
   endif
   
   deallocate( x, f, g, h )
   end do  !  end of loop over queries

   write(lu,*)'testex completed: myid=', myid
   write(lu,600) nquery

   do id = 1, 3
      if( itest(id) == 0 ) cycle
      write(lu,605) id
      write(lu,610) itest(id)
      write(lu,620) errsum(:,id)/itest(id)
      write(lu,630) errmax(:,id)
   end do

600   format('number of queries =', i8)
605   format('table number      =', i8)
610   format('number of tests   =', i8)
620   format('average errors    =', 1p,3e13.4)
630   format('maximum errors    =', 1p,3e13.4)

   if( .not.if_mpi ) stop

   tend = MPI_WTIME() - tstart
!SER   tend = 0.
   if( myid == 0 ) write(0,"(a,1pe13.4)")  &
                   'test_mpi finished: elapsed time = ', tend
  call MPI_FINALIZE(ierr)

   stop
   end program test_mpi
