!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ceq_ind_col(nr,nc,ncs,B,thresh,indcol,info)
!-----------------
  implicit none
  integer, intent(in) :: nr, nc, ncs
  real(kind(1.d0)), intent(in) :: B(nr,nc), thresh
  integer, intent(out) :: indcol(nc), info
!!---------------- 
!  determine independent columns of the matrix B,
!   given that columns 1:ncs are independent.

! Input:
!	nr	- number of rows of B
!	nc	- number of columns of B
!   ncs - index, such that B(:,1:ncs) has full column rank
!   B   - matrix
!   thresh  - threshold for determining rank

! Output:
!   indcol(k)=0 if k-th column is dependent of columns 1:k-1
!   indcol(k)=1 if k-th column is independent of columns 1:k-1
!   info < 0 indicates failure

integer :: lwork, k, jpvt(nc)
real(kind(1.d0)) :: R(nr,nc), tau(nr+nc), work(3*(nr+nc))

indcol=0
lwork=size( work )
R=B
jpvt=0
!  perform QR with column pivoting:  B P = Q R

call dgeqpf(nr,nc,R(1:nr,1:nc),nr,jpvt,tau(1:nr+nc), work(1:lwork), info)

if( info /=0 ) then
   !write(0,*)'ceq_ind_col2: QR failed, info= ', info
   return
end if

do k=1,min(nc,nr)   !  loop over possibly dependent columns
    if( abs(R(k,k)) >= thresh ) then
	   indcol( jpvt(k) ) = 1	   
	else
	   exit
    endif
end do

!  check that the first ncs columns are dependent

do k=1,ncs
   if( indcol(k) /=1 ) then
      info = -2
      !write(0,*)'ceq_ind_col: 1st ncs columns not independent'
      return
   endif
end do

return

end subroutine ceq_ind_col