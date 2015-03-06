!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ceq_min_pert( nz, nc, B, c, z0, eps, z, iret )

!  Determine z = z0 + d which satisfies:
!     1) z(i) >= eps
!     2) B' * z = c
!     3) t = max_i( |d(i)| ) is minimized.

! Input:
!  nz  - length of z
!  nc  - length of c
!  B   - nz x nc equality constraint matrix
!  c   - constraint vector
!  z0  - nz-vector, z0
!  eps - positive threshold

! Output:
!  z    - solution
!  iret = 0 for success

   implicit none
   integer, intent(in)           :: nz, nc
   real(kind(1.d0)), intent(in)  :: B(nz,nc), c(nc), z0(nz), eps
   real(kind(1.d0)), intent(out) :: z(nz)
   integer, intent(out)          :: iret

   real(kind(1.d0)) :: A(nc+nz,3*nz), x(3*nz), f(3*nz), r(nc+nz), d(nz), tp, tm, res, &
                       Bsc, csc, zsc, ceq_norm
   integer :: i, nx, nr
   logical :: linear = .false.

! x = [ u v w ] = [ z - eps (t+d)/2   (t-d)/2  ]
! f = [ 0 1 1 1   [ 0 0 0     1 1 1    1 1 1   ] 
! minimize sum( f*x ) = nz * t subject to  A x = r

   nx = 3*nz
   nr = nc+nz

   Bsc = max( maxval(B), -minval(B) )  !  scale factors
   csc = max( maxval(c), -minval(c) )

   zsc = csc / Bsc

   A(1:nz+nc,1:3*nz) = 0.d0

   if( linear ) then  !  min. sum of dz
      do i = 1, nz
         A(i,i)      =  1.d0
	     A(i,nz+i)   = -1.d0
	     A(i,2*nz+i) =  1.d0
      end do
   else  !  min. sum of dz / z0
      do i = 1, nz
         A(i,i)      =  1.d0
	     A(i,nz+i)   = -(max( z0(i), eps )/zsc)**0.0d0
	     A(i,2*nz+i) =  (max( z0(i), eps )/zsc)**0.0d0
      end do
   endif

   A(nz+1:nz+nc,1:nz)   = transpose( B ) / Bsc

   r(1:nz)       = ( z0(1:nz) -eps ) /zsc

   do i = 1, nc
      r(nz+i) = ( c(i) - eps * sum(B(:,i)) ) / csc
   end do

   f(1:nz)      = 0.d0
   do i = 1, nz
      f(nz+i)   = ( eps/(eps+z0(i)) )**0.0d0
	  f(2*nz+i) = f(nz+i)
   end do

   call ceq_linprog( nx, nr, f, A, r, x, iret )

   if( iret == 0 ) then
      do i = 1, nz
	     z(i) = max( x(i), 0.d0 ) * zsc + eps
	  end do
   endif

   if( .false. ) return

!  diagnostics

   write(0,*)' '
   write(0,*)'ceq_min_pert: iret = ', iret
   write(0,*)'ceq_min_pert: z0, z, z0-z'
   do i = 1, nz
      write(0,'(1p,10e11.2)') z0(i), z(i), z0(i)-z(i), x(nz+i), x(2*nz+i), &
	                          x(nz+i)+x(2*nz+i), x(nz+i)-x(2*nz+i)
   end do

   write(0,*) 'ceq_min_pert: min(x) = ', minval(x)
   r(1:nc) = matmul( z(1:nz), B(1:nz,1:nc) ) - c(1:nc)
   write(0,*) 'ceq_min_pert: residual = ', ceq_norm( nc, r(1:nc) )
   write(0,*) 'ceq_min_pert: sum of perts. = ', sum( x(nz+1:3*nz) )

   return
end subroutine ceq_min_pert
