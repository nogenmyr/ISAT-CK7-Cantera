!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module ci_utils
  use ci_prec

  interface diagmul

     ! This module provides 'diagmul' routine to compute product of two
     ! matrices (similar to fortran's inbuilt 'matmul' routine) which is
     ! optimized for diagonal matrices.

     ! One of the arguments (x) to the 'diagmul' routine should be a 1D
     ! array containing the diagonal matrix entries, and the other
     ! argument could be a 1D array (a) or any general 2d matrix (A).
     ! The matrices are multiplied in the same order in which they are
     ! passed to the routine. So, for example:
     !
     ! diagmul(x,A) returns diag(x)*A ;and
     ! diagmul(A,x) returns A*diag(x)
     !
     ! where of course if both the matrices are 1D the order doesn't
     ! matter. So,
     !
     ! diagmul(a,x) = diagmul(x,a) = a*diag(x) = diag(x)*a

     module procedure diag1d, diag_pre2d, diag_post2d
  end interface

  interface print_var
     module procedure print_real1d, print_real2d, print_int1d, print_int2d, print_char
  end interface
contains

  function diag1d( x, a )

    !  Multiply the vector a by X=diag(x), i.e., return diag(x)*a
    implicit none

    real(k_dp), intent(in)             :: x(:), a(:)
    real(k_dp)                         :: diag1d(size(a))

    if (size(a) /= size(x)) then
       write(0,*) 'The shapes of the arguments are inconsistent or nonconformable.   [DIAGMUL]'
       stop
    endif

    diag1d = x*a

    return
  end function diag1d

  function diag_pre2d( x, A )

    !  Pre-multiply the m x n matrix A by X=diag(x), i.e., return diag(x)*A
    implicit none

    real(k_dp), intent(in)             :: x(:), A(:,:)
    real(k_dp)                         :: diag_pre2d(size(A,1),size(A,2))

    integer                            :: j

    if (size(A,1) /= size(x)) then
       write(0,*) 'The shapes of the arguments are inconsistent or nonconformable.   [DIAGMUL]'
       stop
    endif

    do j = 1, size(A,2)
       diag_pre2d(:,j) = x * A(:,j)
    end do

    return
  end function diag_pre2d

  function diag_post2d( B, y)

    !  Post-multiply the m x n matrix B by Y=diag(y), i.e., return B*diag(y)
    implicit none

    real(k_dp), intent(in)             :: y(:), B(:,:)
    real(k_dp)                         :: diag_post2d(size(B,1),size(B,2))

    integer                            :: j

    if (size(B,2) /= size(y)) then
       write(0,*) 'The shapes of the arguments are inconsistent or nonconformable.   [DIAGMUL]'
       stop
    endif

    do j = 1, size(B,2)
       diag_post2d(:,j) = B(:,j) * y(j)
    end do

    return
  end function diag_post2d

  subroutine print_char(cn, name, luout)
    ! print 1d real variable to luout
    implicit none

    character(*), intent(in)           :: cn(:)
    character(*), optional, intent(in) :: name
    integer, optional, intent(in)      :: luout

    ! local variables
    integer :: lu
    if(present(luout)) lu = luout

    write(lu,*) '  '
    write(lu,*) name
    write(lu,'(100A)') '  ',cn

  end subroutine print_char

  subroutine print_real1d(x, name, luout)
    ! print 1d real variable to luout
    implicit none

    real(k_dp), intent(in)             :: x(:)
    character(*), optional, intent(in) :: name
    integer, optional, intent(in)      :: luout

    ! local variables
    integer :: lu
    if(present(luout)) lu = luout

    write(lu,*) '  '
    write(lu,*) name
    write(lu,'(1p,200e10.2)') x

  end subroutine print_real1d

  subroutine print_int1d(x, name, luout)
    ! print 1d integer variable to luout
    implicit none

    integer, intent(in)             :: x(:)
    character(*), optional, intent(in) :: name
    integer, optional, intent(in)      :: luout

    ! local variables
    integer :: lu
    if(present(luout)) lu = luout

    write(lu,*) '  '
    write(lu,*) name
    write(lu,'(1p,200i10)') x

  end subroutine print_int1d

  subroutine print_real2d(A, name, luout)
    ! print 2d real variable to luout
    implicit none

    real(k_dp), intent(in)             :: A(:,:)
    character(*), optional, intent(in) :: name
    integer, optional, intent(in)      :: luout

    ! local variables
    integer :: i, j, lu
    if(present(luout)) lu = luout

    write(lu,*) '  '
    write(lu,*) name
    do i = 1, size(A,1)
       write(lu,'(1p,200e10.2)') (A(i,j),j=1,size(A,2))
    end do
  end subroutine print_real2d

  subroutine print_int2d(A, name, luout)
    ! print 2d integer variable to luout
    implicit none

    integer, intent(in)                :: A(:,:)
    character(*), optional, intent(in) :: name
    integer, optional, intent(in)      :: luout

    ! local variables
    integer :: i, j, lu
    if(present(luout)) lu = luout

    write(lu,*) '  '
    write(lu,*) name
    do i = 1, size(A,1)
       write(lu,'(1p,200i10)') (A(i,j),j=1,size(A,2))
    end do
  end subroutine print_int2d

  function norm(x)
    ! Returns 2-norm of a vector x
    implicit none
    real(k_dp), intent(in)   :: x(:)
    real(k_dp)               :: norm

    norm = sqrt(dot_product(x,x))
    return
  end function norm

  function snorm(sx)
    ! Returns 2-norm of a vector x (single precision)
    implicit none
    real(k_sp), intent(in)   :: sx(:)
    real(k_sp)               :: snorm

    snorm = sqrt(dot_product(sx,sx))
    return
  end function snorm

  subroutine sort(x)
    ! sort input array x in increasing order
    implicit none
    integer :: x(:)
    integer :: i, j, t, done

    do i = 1,size(x)
       done = 1
       do j = 2,size(x)
          if(x(j)<x(j-1)) then
             t = x(j-1)
             x(j-1) = x(j)
             x(j) = t
             done = 0
          endif
       enddo
       if(done==1) return
    enddo
  end subroutine sort


  function norm_dot(x,y)   !------------------------------------------------------
    ! Returns the normalized dot product of two vectors defined by:
    !   norm_dot(x,y) = x.y/max(|x|,|y|)^2
    ! If both x and y are zero, then 0 is returned
    implicit none
    real(k_dp), intent(in)   :: x(:), y(:)
    real(k_dp)               :: norm_dot, longest

    longest = max( norm(x), norm(y) )
    if( longest > 0.d0 ) then
       norm_dot = dot_product(x,y) / longest**2
    else
       norm_dot = 0.d0
    endif
    return
  end function norm_dot

  function eye(n)
    ! Returns n x n identity matrix
    implicit none
    integer, intent(in) :: n
    real(k_dp)          :: eye(n,n)
    integer :: i

    eye = 0.d0
    do i = 1, n
       eye(i,i) = 1.d0
    end do

    return
  end function eye

  subroutine re_norm( x, w, beta )

    !  Re-normalize a rate-of change of composition vector, x, 
    !  so that  w^T x = 0.
    !  x  is replaced by x + beta * abs(x), where beta = -w^T x / w^t abs(x).
    !  Optionally, beta is returned.

    use isat_abort_m

    real(k_dp), intent(inout) :: x(:)
    real(k_dp), intent(in)    :: w(:)
    real(k_dp), optional, intent(out) :: beta

    real(k_dp)       :: beta1, denom

    denom = dot_product( w, abs(x) )

    if( denom <= 0.d0 ) return  !  x=0 satisfies normalization

    beta1 = -dot_product( w, x ) / denom
    x    = x + beta1 * abs(x)

    if( present(beta) ) beta = beta1

    return
  end subroutine re_norm

  subroutine principal_angles( m, na, nb, A, B, rankA, rankB, theta )  !------------------

    !  Determine the principal angles  theta(1:nth)  between the subspaces spanned 
    !  by the columns of the m x na matrix A and the m x nb matrix B, where
    !  nth = min( rankA, rankB ), with rankA=rank(A) and rankB=rank(B). 
    !  This implementation is based on Algorithm 3.1 of Knyazev & Argentati (2002),
    !  SIAM J Sci Comput 23, 2009-2041.

    use isat_abort_m

    integer,    intent(in)  :: m, na, nb
    real(k_dp), intent(in)  :: A(m,na), B(m,nb)
    integer,    intent(out) :: rankA, rankB
    real(k_dp), intent(out) :: theta( min(na,nb) )

    real(k_dp), parameter :: epsilon = 1.d-12

    integer    :: lwork, infoA, infoB, i
    real(k_dp) :: AA(m,max(na,nb)), BB(m,max(na,nb)), SA(min(m,na)), SB(min(m,nb)), &
         UA(m,na), UB(m,nb), VT(1,1), work(2*m*(na+nb)+20*(m+na+nb)), &
         tol, sigma(max(na,nb)), mu(max(na,nb))

    lwork = 2*m*(na+nb)+20*(m+na+nb)

    AA  = A
    BB  = B

    ! form UA and UB: orthogonal bases for span(A) and span(B)
    call dgesvd( 'S', 'N', m, na, AA, m, SA, UA, m, VT, 1, &
         work, lwork, infoA )

    call dgesvd( 'S', 'N', m, nb, BB, m, SB, UB, m, VT, 1, &
         work, lwork, infoB )

    if( infoA /= 0  .or.  infoB /= 0 ) call isat_abort( 'principal_angles', 1, &
         mess=' SVD failed, infoA, infoB = ', &
         ivar = (/ infoA, infoB /) )

    !  determine ranks
    tol = max(m,na) * SA(1) * epsilon
    do i = 1, min(m,na)
       if( SA(i) < tol ) exit
       rankA = i
    end do

    tol = max(m,nb) * SB(1) * epsilon
    do i = 1, min(m,nb)
       if( SB(i) < tol ) exit
       rankB = i
    end do

    theta = 0.d0
    nth   = min( rankA, rankB )
    if( nth == 0 ) return

    AA(1:rankA,1:rankB) = matmul( transpose(UA(:,1:rankA)), UB(:,1:rankB) )
    BB(1:rankA,1:rankB) = AA(1:rankA,1:rankB)

    call dgesvd( 'N', 'N', rankA, rankB, BB, m, sigma, UB, m, VT, 1, &
         work, lwork, infoA )

    if( rankA >= rankB ) then
       BB(:,1:rankB) = UB(:,1:rankB) - matmul( UA(:,1:rankA), AA(1:rankA,1:rankB) )
    else
       BB(:,1:rankA) = UA(:,1:rankA) - matmul( UB(:,1:rankB), transpose(AA(1:rankA,1:rankB)) )
    endif

    call dgesvd( 'N', 'N', m, nth, BB, m, mu, UB, m, VT, 1, &
         work, lwork, infoB )

    if( infoA /= 0  .or.  infoB /= 0 ) call isat_abort( 'principal_angles', 2, &
         mess=' SVD failed, infoA, infoB = ', &
         ivar = (/ infoA, infoB /) )  

    do i = 1, nth
       j = nth+1-i
       if( sigma(j)**2 < 0.5d0 ) then
          theta(i) = acos( min( sigma(j), 1.d0 ) )
       else
          theta(i) = asin( min( mu(i), 1.d0 ) )
       endif
    end do

    return                                  

  end subroutine principal_angles

  subroutine test_principal_angles  !-------------------------------------------------------

    !  test principal_angles

    integer    :: lu, n, na, nb, nth, rankA, rankB
    real(k_dp) :: fn, fna, fnb
    real(k_dp), allocatable :: A(:,:), B(:,:), theta(:)

    call isat_lu(lu)
    open( lu, file = 'nab.txt' )
    read(lu,*) fn, fna, fnb  !  file generated by Matlab test code subspace.m
    close(lu)
    n  = nint(fn)
    na = nint(fna)
    nb = nint(fnb)

    write(0,*)'test_principal_angles: n, na, nb = ', n, na, nb

    allocate( A(n,na) )
    open( lu, file = 'Amat.txt' )
    read(lu,*) A
    close(lu) 

    allocate( B(n,nb) )
    open( lu, file = 'Bmat.txt' )
    read(lu,*) B
    close(lu)

    allocate( theta( min(na,nb) ) )

    call principal_angles( n, nb, na, B, A, rankA, rankB, theta )

    write(0,*)'test_principal_angles: rankA, rankB = ', rankA, rankB
    write(0,'((1p,e13.4))') (theta(1:nth))

    return

  end subroutine test_principal_angles

  subroutine svd( m, n, A, S )  !------------------------------------------

    !  Determine the singular values  S(1:min(m,n))  of the m x n matrix A.

    use isat_abort_m

    integer,    intent(in)  :: m, n
    real(k_dp), intent(in)  :: A(m,n)
    real(k_dp), intent(out) :: S( min(m,n) )


    integer    :: lwork, info
    real(k_dp) :: AA(m,n), U(1,1), VT(1,1), work(2*m*n+20*(m+n))

    lwork = 2*m*n+20*(m+n)

    AA  = A

    call dgesvd( 'N', 'N', m, n, AA, m, S, U, 1, VT, 1, &
         work, lwork, info )

    if( info /= 0  ) call isat_abort( 'svd', 1, mess=' SVD failed, info= ', isv=info )

    return                                  

  end subroutine svd


  subroutine BT_red( nrs, CS, BT )
    ! Computes the transpose of the reduction matrix for a given set
    ! of represented species.

    ! Input
    !   nrs: number of represented species
    !   CS: index of represented species
    
    ! Output
    !   BT: Transpose of Reduction matrix

    use ci_dat6
    implicit none

    integer, intent(in) :: nrs, CS(nrs)
    real(kind(1.d0)), intent(out) :: BT(nrs+ne, ns)

    ! Local variables
    integer :: i, ev(ns, ne)

    ! Element Vector
    ev(1:ns, 1:ne) = dev(1:ns, 1:ne)

    BT = 0.d0
    do i = 1, nrs
       BT(i, CS(i)) = 1
       ev(CS(i), 1:ne) = 0 ! set ev to zero for repr. species
    enddo

    BT(nrs+1:nrs+ne,1:ns) = transpose(ev(1:ns, 1:ne))

  end subroutine BT_red

  subroutine Mappings( nrs, CS, MR, MU )
    ! Rep. and Unrep. Mapping Matricies

    use ci_dat6
    implicit none

    integer, intent(in) :: nrs, CS(nrs)
    real(kind(1.d0)), intent(out) :: MR(nrs, ns), MU(ns-nrs, ns)
    
    ! local
    integer :: i, j, indices(ns)

    ! Mark the represented species
    indices = 0
    do i = 1, ns
       indices(CS(i)) = 1
    enddo

    ! Represented mapping
    MR = 0.d0
    do i = 1, nrs
       MR(i, CS(i)) = 1
    enddo

    ! Unrepresented mapping
    MU = 0.d0
    j = 1
    do i = 1, ns
       if(indices(i) == 0) then
          MU(j, i) = 1
          j = j + 1
       endif
    enddo

  end subroutine Mappings

end module ci_utils
