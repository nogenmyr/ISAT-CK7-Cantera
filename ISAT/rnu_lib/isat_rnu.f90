!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module isat_rnu

!  contains the following user-callable subroutines to generate 
!  pseudo-random numbers uniformly in the  exclusive interval [ 0 , 1 ].

!     rnu( x )           returns random numbers in x, where x is single
!                        or double precision, and of rank 0, 1 or 2.
!     rnuget( i1, i2 )   get random number generator seeds
!     rnuput( i1, i2 )   set random number generator seeds
!     rnused( k )        set seeds to k-th initialization values

   interface rnu
      module procedure rnus0, rnus1, rnus2, rnud0, rnud1, rnud2, rnui1
   end interface

   integer, save, private :: i1=12345, i2=67890
   integer, save           :: lu_err = 0
   contains

!---------------------------------------------------------------------------

   subroutine rnus1( x )

!  routine to generate pseudo-random numbers uniformly in the
!   exclusive interval [ 0 , 1 ].

!  november 1990;  fortran90 version may 2000.

!   x(i) , i = 1, n  array in which random numbers are returned

!   the entries  rnuget  and  rnuput  can be used to get and set
!   the seeds i1 and i2.

!  method: see f. james, computer physics communications, vol.60
!   p.329, 1990.

   implicit none

   real(kind(1.0e0)), intent(out) :: x(:)

   integer :: n, i, ik, ix

   n = size(x)

   do i = 1, n
      ik = i1 / 53668
      i1 = 40014 * ( i1 - ik * 53668 ) - ik * 12211
      if( i1 < 0 ) i1 = i1 + 2147483563

      ik = i2 / 52774
      i2 = 40692 * ( i2 - ik * 52774 ) - ik * 3791
      if( i2 < 0 ) i2 = i2 + 2147483399

      ix = i1 - i2
      if( ix < 1 ) ix = ix + 2147483562

      x(i) = float( ix ) * 4.656612e-10
   end do

   return
   end subroutine rnus1

!---------------------------------------------------------------------------

   subroutine rnud1( x )

   implicit none

   real(kind(1.0d0)), intent(out) :: x(:)

   integer :: n, i, ik, ix

   n = size(x)

   do i = 1, n
      ik = i1 / 53668
      i1 = 40014 * ( i1 - ik * 53668 ) - ik * 12211
      if( i1 < 0 ) i1 = i1 + 2147483563

      ik = i2 / 52774
      i2 = 40692 * ( i2 - ik * 52774 ) - ik * 3791
      if( i2 < 0 ) i2 = i2 + 2147483399

      ix = i1 - i2
      if( ix < 1 ) ix = ix + 2147483562

      x(i) = float( ix ) * 4.656612e-10
   end do

   return
   end subroutine rnud1

!---------------------------------------------------------------------------

   subroutine rnui1( j )

   implicit none

   integer, intent(out) :: j(:)

   integer :: n, i, ik, ix

   n = size(j)

   do i = 1, n
      ik = i1 / 53668
      i1 = 40014 * ( i1 - ik * 53668 ) - ik * 12211
      if( i1 < 0 ) i1 = i1 + 2147483563

      ik = i2 / 52774
      i2 = 40692 * ( i2 - ik * 52774 ) - ik * 3791
      if( i2 < 0 ) i2 = i2 + 2147483399

      ix = i1 - i2
      if( ix < 1 ) ix = ix + 2147483562

      j(i) = ix
   end do

   return
   end subroutine rnui1

!---------------------------------------------------------------------------

   subroutine rnus2( x )

   implicit none

   real(kind(1.0e0)), intent(out) :: x(:,:)

   integer :: j

   do j = 1, size(x,2)
      call rnus1( x(:,j) )
   end do

   return
   end subroutine rnus2


!---------------------------------------------------------------------------

   subroutine rnud2( x )

   implicit none

   real(kind(1.0d0)), intent(out) :: x(:,:)

   integer :: j

   do j = 1, size(x,2)
      call rnud1( x(:,j) )
   end do

   return
   end subroutine rnud2

!---------------------------------------------------------------------------

   subroutine rnus0( x )

   implicit none
   real(kind(1.0e0)), intent(out) :: x
   real(kind(1.0e0)) :: xx(1)

   call rnus1( xx )
   x = xx(1)

   return
   end subroutine rnus0

!---------------------------------------------------------------------------

   subroutine rnud0( x )

   implicit none
   real(kind(1.0d0)), intent(out) :: x
   real(kind(1.0e0)) :: xx(1)

   call rnus1( xx )
   x = xx(1)

   return
   end subroutine rnud0

!---------------------------------------------------------------------------

   subroutine rnuget( ig1, ig2 )

   implicit none
   integer, intent(out) :: ig1, ig2

   ig1 = i1
   ig2 = i2

   return
   end subroutine rnuget

!---------------------------------------------------------------------------

   subroutine rnuput( ip1, ip2 )

   implicit none
   integer, intent(in) :: ip1, ip2

   i1 = ip1
   i2 = ip2

   return
   end subroutine rnuput

!---------------------------------------------------------------------------

   subroutine rnused( ks )

   implicit none

   integer, intent(in) :: ks
   integer :: i, ik, i1t
   
   i2 = 12345
   i1 = 67890 

   if( ks < 0 ) then
      write(lu_err,*)' rnused: warning, ks < 0 ', ks
      return
   endif

   do i = 1, ks
      ik = i1 / 53668
      i1 = 40014 * ( i1 - ik * 53668 ) - ik * 12211
      if( i1 < 0 ) i1 = i1 + 2147483563

      ik = i2 / 52774
      i2 = 40692 * ( i2 - ik * 52774 ) - ik * 3791
      if( i2 < 0 ) i2 = i2 + 2147483399
   end do

   i1t = i1
   i1  = i2
   i2  = i1t

   return
   end subroutine rnused

!---------------------------------------------------------------------------

end module isat_rnu
