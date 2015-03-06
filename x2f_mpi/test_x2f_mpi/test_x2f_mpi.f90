!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the x2f_mpi      !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

program test_x2f_mpi

! A battery of tests calling subroutine x2f_mpi just to make sure it works
! ...see the source for x2f_mpi to learn the definitions of its arguments
! likewise, subroutine x2f is given particularly simple forms--see below

use MPI

implicit none
integer, parameter       :: k_dp = kind(1.d0)
integer, parameter       :: k_pos = 1, gat = 1, n_spec1 = 12, n_spec2 = 1

real(k_dp), allocatable  :: x(:,:), xf(:,:), f0(:,:), errorvec(:)
integer, allocatable     :: seednum(:)
integer    :: ierr, myrank, testnum, seedsize, i, nhit, sumnhit, sumnv
real       :: x1, x2
real(k_dp) :: approxpi, avtenth
logical    :: errorfree, parallel_error = .false.

integer    :: j, selectx2f, iat, nv, nx, ldxf, nf, highnf
external   :: car2pol, tenthx1, onemethod
integer    :: info(gat,n_spec1), info_xf(1)
real(k_dp) :: rinfo(gat,n_spec2), rinfo_xf(1)
real(k_dp) :: stats(100)

call MPI_Init( ierr )
call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr)  ! set myrank with this call

nv = 600000
! purposely imbalance numbers of values across processors
nv = nv * ( myrank + 1 ) - int( 0.3 * ( myrank + 1 ) )

! only one global attempt because supplied methods work every time
iat = gat

info(iat,2)  = 0  ! MPI communicator defaults to MPI_COMM_WORLD if info(iat,2) is 0
info(iat,3)  = 2  ! u(2) = 1, establishes unit vector for projection in STRATA
info(iat,4)  = 4  ! number of bins per processor in STRATA

! ---------------------- first test set begins here ----------------------
!
! x2f is car2pol; nx = 2; ldxf = 2; nf = 1, 2, 3
! try to generate error for low ldxf, then reset ldxf = 3

nx     = 2
ldxf   = 2
highnf = 3

allocate( x(nx,nv) )
allocate( xf(ldxf,nv) )
allocate( f0(highnf,nv) )
allocate( errorvec(highnf) )

call random_seed( size = seedsize )
allocate( seednum(seedsize) )
seednum = (/ ( 100 * myrank + i, i = 1, seedsize ) /)
call random_seed( put = seednum(1:seedsize) )
do i = 1, nv
   call random_number( x1 )
   call random_number( x2 )
   x(1,i) = 2.0 * ( x1 - 0.5 )
   x(2,i) = 2.0 * ( x2 - 0.5 )
   call car2pol( ldxf, nx, x(1:nx,i), f0(1:ldxf,i), info_xf, rinfo_xf )
end do

do j = 0, 3  ! -------- first test loop begins here ----------------------

info(iat,2) = j  ! testing mode_p = 0, 1, 2, 3
print "('starting car2pol tests for mode_p = ',i2)", info(iat,2)

!do nf = 1, highnf
do nf = 2, 2

   errorfree = .false.
   !do while ( .not. errorfree )

      xf(1:nx,1:nv) = x  ! prepare for x2f_mpi function evaluation
      call x2f_mpi( nx, nf, ldxf, nv, xf, k_pos, info, rinfo, info_xf, rinfo_xf, &
                    gat, n_spec1, n_spec2, car2pol, onemethod, stats )

      if ( stats(1) /= -4.d0 ) then
         errorfree = .true.
      else
         print "('error -4 detected!  nx =',i2,',  nf =',i2,',  ldxf =',i2)", nx, nf, ldxf
         deallocate( xf )
         ldxf = nf
         allocate( xf(ldxf,nv) )
      end if

   !end do

! begin little calculation to compute an approximation to pi

   nhit = 0
   do i = 1, nv
      if( xf(1,i) .le. 1.0d0 ) then
         nhit = nhit + 1
      end if
   end do

   approxpi = 4.0d0 * dfloat( nhit ) / dfloat( nv )

   errorvec(1:nf) = sum( dabs( xf(1:nf,1:nv) - f0(1:nf,1:nv) ), dim = 2 )

! print values, statistics, and error due to parallelized computation

   print "('rank ',i2,':',f12.8,' is the value of pi based on',i8,' points')", &
           myrank, approxpi, nv
   print "('processed--',i8,' local,',i8,' n.r.,',i8,' tot.')", &
           int( stats(31) ), int( stats(32) ), int( stats(31) + stats(32) )
   print "('CPU,  sec--',f8.3,' local,',f8.3,' n.r.,',f8.3,' tot.')", &
           stats(3), stats(2), stats(1)
   print "('wall, sec--',f8.3,' local,',f8.3,' n.r.,',f8.3,' tot.')", &
           stats(23), stats(22), stats(21)

   if( maxval( errorvec(1:nf) ) > 0.0 ) then
      parallel_error = .true.
      print "('TEST FAILED!')"
      print "('error sum-- ',3d19.12)", errorvec(1:nf)
   end if

! get a better estimate of pi by combining results from all processors

   call MPI_Reduce( nhit, sumnhit, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   call MPI_Reduce( nv, sumnv, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   if ( myrank == 0 ) then
      approxpi = 4.0d0 * dfloat( sumnhit ) / dfloat( sumnv )
      print "('for all processors: ',f12.8,' is the value of pi')", approxpi
   end if

end do

end do  ! ------------ first test loop ends here -------------------------


! ------------------ end first test set; begin second --------------------
!
! x2f is tenthx1; nx = 2; ldxf = 3; nf = 2, 1
! verify there is no error generated for high ldxf, so ldxf is never reset

nx     = 2
ldxf   = 2
highnf = 2

do i = 1, nv
   x(1,i) = dfloat(i)
   x(2,i) = dfloat(mod(i,10))
!   x(2,i) = 1.d0
   call tenthx1( highnf, nx, x(1:nx,i), f0(1:highnf,i), info_xf, rinfo_xf )
end do

do j = 0, 3  ! -------- second test loop begins here ---------------------

info(iat,2) = j  ! testing mode_p = 0, 1, 2, 3
print "('starting tenthx1 tests for mode_p = ',i2)", info(iat,2)

!do nf = highnf, 1, -1
do nf = 2, 2

   errorfree = .false.
   do while ( .not. errorfree )

      xf(1:nx,1:nv) = x  ! prepare for x2f_mpi function evaluation
      call x2f_mpi( nx, nf, ldxf, nv, xf, k_pos, info, rinfo, info_xf, rinfo_xf, &
                    gat, n_spec1, n_spec2, tenthx1, onemethod, stats )
   
      if ( stats(1) /= -4.d0 ) then
         errorfree = .true.
      else
         print "('error -4 detected!  nx =',i2,',  nf =',i2,',  ldxf =',i2)", nx, nf, ldxf
         deallocate( xf )
         ldxf = nf
         allocate( xf(ldxf,nv) )
      end if

   end do

   avtenth = sum( xf(1,:) )/dfloat( nv )

   errorvec(1:nf) = sum( dabs( xf(1:nf,1:nv) - f0(1:nf,1:nv) ), dim = 2 )
   print "('rank ',i2,':',f11.4,' is 1/10 the average value of x1 for nv = ',i8)", &
           myrank, avtenth, nv
   print "('processed--',i8,' local,',i8,' n.r.,',i8,' tot.')", &
           int( stats(31) ), int( stats(32) ), int( stats(31) + stats(32) )
   print "('CPU,  sec--',f8.3,' local,',f8.3,' n.r.,',f8.3,' tot.')", &
           stats(3), stats(2), stats(1)
   print "('wall, sec--',f8.3,' local,',f8.3,' n.r.,',f8.3,' tot.')", &
           stats(23), stats(22), stats(21)

   if( maxval( errorvec(1:nf) ) > 0.0 ) then
      parallel_error = .true.
      print "('TEST FAILED!')"
      print "('error sum-- ',3d19.12)", errorvec(1:nf)
   end if

end do

end do  ! ------------ second test loop ends here ------------------------

if( .not. parallel_error) then
   print "('All parallel tests succeeded!')"
else
   print "('Some parallel tests failed...')"
end if

call MPI_Barrier( MPI_COMM_WORLD, ierr )
call MPI_Finalize( ierr )

stop
end program test_x2f_mpi


!------------------------------------------------------------------------

subroutine car2pol( nf, nx, car, pol, info_xf, rinfo_xf )

! takes car(1) and car(2) as 2d Cartesian coords, converts them to polar
! on return, pol(1) is the radial coord, pol(2) is the azimuthal
! if nx == 1, the result is them same as if car(2) were equal to 0
! if nf == 1, polar angle is not returned; radial coord only is in pol(1)

implicit none
integer, parameter :: k_dp = kind(1.d0)

integer    :: qt_mode = 0, nx, nf
real(k_dp) :: car(nx), pol(nf), ytest, ycoord
integer    :: info_xf(*)
real(k_dp) :: rinfo_xf(*)

if( qt_mode == -1 ) then
   qt_mode = 0
   return
elseif( qt_mode == 1 ) then
   ! want quick try only when ( car(1), car(2) ) is in second quadrant
   ! but prevent quick try for nx == 1; in effect, put point in quadrant 1 or 4
   ytest = car(min(2,nx))
   if( ( car(1) >= 0 ) .or. ( ytest <= 0 ) ) then
      qt_mode = -1
      return
   else
      qt_mode = 0
   end if
end if

pol(1:nf) = 0.0d0

if ( nx > 1 ) then
   ycoord = car(2)
else
   ycoord = 0.d0
end if

pol(1) = dsqrt( car(1)*car(1) + ycoord*ycoord )

if ( nf > 1 ) then
   pol(2) = atan2( ycoord, car(1) )
end if    

return
end subroutine car2pol


!------------------------------------------------------------------------

subroutine tenthx1( ldxf, nx, vec, vec1d10, info_xf, rinfo_xf )

implicit none
integer, parameter :: k_dp = kind(1.d0)

integer    :: qt_mode = 0, nx, ldxf, i
real(k_dp) :: vec(nx), vec1d10(ldxf)
integer    :: info_xf(*)
real(k_dp) :: rinfo_xf(*)

if( qt_mode == -1 ) then
   qt_mode = 0
   return
elseif( qt_mode == 1 ) then
   ! quick try for only about 1/3 of input if broadly distributed
   if( mod(int(vec(1)),3) /= 0 ) then
      qt_mode = -1
      return
   else
      qt_mode = 0
   end if
end if

do i = 1, ldxf
   if ( i <= nx ) then
      vec1d10(i) = vec(i)
   else
      vec1d10(i) = 0.0d0
   end if
end do

vec1d10(1) = 1.0d-1 * vec(1)

return
end subroutine tenthx1

subroutine onemethod(method,info_xf,rinfo_xf)

implicit none
integer :: method, info_xf(*), rinfo_xf(*)

return
end  