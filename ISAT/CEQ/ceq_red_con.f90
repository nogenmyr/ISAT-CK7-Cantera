!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ceq_red_con(ns,ne,ncs,ng,nc,Ein,CS,Bg,ifop,lud, &
	B,ned,neu,el_order,nsd,nsu,sp_order,E,nrc,nb,BR,A,iret)

!-----------------
  implicit none
  integer, intent(in) :: ns, ne, ncs, ng, nc, CS(ncs), ifop, lud 
  real(kind(1.d0)), intent(in) :: Ein(ns,ne), Bg(ns,ng)
  integer, intent(out) :: ned, neu, el_order(ne), nsd, nsu, &
     sp_order(ns), nrc, nb, iret
  real(kind(1.d0)), intent(out) :: B(ns,nc), E(ns,ne), BR(ns,nc), A(nc,nc)
!-----------------
! Set the reduced constraint matrix BR, and determine the ordering of 
!   the elements and species.

! Input:
!	ns		- number of species
!	ne		- number of elements
!	ncs		- number of constrained species
!	ng		- number of general linear constraints
!	nc		- number of constraints = ne+ncs+ng
!	Ein     - element matrix (with user-input ordering of elements and 
!               species) (ns x ne)
!	CS		- column vector of indices of constrained species (ncs x 1)
!	Bg		- general constraint matrix (ns x ng)
!	ifop	> 0 for diagnostic output
!	lud		- logical unit for diagnostic output
!
!  Output:
!   B           - basic constraint matrix (ns x nc)
!	ned     	- number of determined elements
!	neu     	- number of undetermined elements
!	el_order	- the j-th ordered element is el_order(j) (ne x 1)
!	nsd     	- number of determined species
!	nsu	        - number of undetermined species
!	sp_order	- the k-th ordered species is sp_order(k) (ns x 1)
!	E			- the element matrix with re-ordered elements and 
!               species (ns x ne)
!	nrc			- number of reduced constraints
!	nb			= number of independent constraints =nsd+nrc
!	BR			- the reduced constraint matrix (nsu x nrc)
!   A           - modified constraint transformation matrix (nb x nc)
!   iret        < 0 for failure
!   iret        = -1, error in one of first 5 arguments
!   iret        = -2, the specified constrained species are not distinct
!   iret        = -3, (SVD failure)
!   iret        = -4, (ceq_ind_col failed)


integer :: k, info, lwork, rank_def, sp_det(ns), kk, ndebg, el_det(ne), &
  ngi, is, js, ie, je, j, indcol(ns+ne+ng)
real(kind(1.d0)) :: Sc(ns,ncs),S(ns+nc),U(ns,ns),VT(nc,nc), &
  work(4*(ns+nc)*(ns+nc)), DEBg(ns,ns+ne+ng), EU(ns,ne), BS(ns,nc), &
  PU(ns,nc), BTTPU(nc,nc), BB(ns,nc), ceq_norm
real(kind(1.d0)), parameter ::  thresh=1.d-10 ! threshold for singular values and vectors

iret = 0  ! anticipate success

! print out dimensions and input
if( ifop>=5 ) then
    write(lud,*)' '
    write(lud,*)'Input to ceq_red_con'
    write(lud,*)'ne = ', ne
    write(lud,*)'ns = ', ns
    write(lud,*)'ncs= ', ncs
    write(lud,*)'ng=  ', ng
    write(lud,*)'nc = ', nc
    write(lud,*)'Constrained species'
    write(lud,'((20i4))')CS
endif
    
!  check input

call check_input
if( iret < 0 ) return  !  bad input

if( .true. ) then !  new test for constrained species being distinct SBP 2/6/09 

do j = 1, ncs-1
   do k = j+1, ncs
      if( CS(j) == CS(k) ) then
         iret=-2
		 write(0,*)'ceq_red_con: CS not distinct'
		 return
      endif
   end do
end do

else  ! original test which inappropriately assumes ordeering

if( ncs>1 ) then
    if( minval(CS(2:ncs)-CS(1:ncs-1)) <1 ) then
        iret=-2
		write(0,*)'ceq_red_con: CS not distinct'
		return
    endif
endif

endif  !  end of distinct test

iret=0

B=0.d0     ! set basic constraint matrix
! set constrained species matrix
if( ncs>0 ) then
  Sc=0.d0
  do k=1,ncs 
    Sc(CS(k),k)=1.d0
  end do 
  B(:,1:ncs)=Sc              ! first ncs columns - constrained species
endif

B(1:ns,ncs+1:ncs+ne)=Ein(1:ns,1:ne)      ! next ne columns - element vectors

if( ng>0 ) &
B(1:ns,ne+ncs+1:nc)=Bg(1:ns,1:ng)        ! last ng columns - general constraints

BB(1:ns,1:nc)=B(1:ns,1:nc)

lwork = size(work)
call dgesvd('A','A',ns,nc,BB(1:ns,1:nc),ns,S(1:ns+nc),U(1:ns,1:ns),ns, &
             VT(1:nc,1:nc),nc,work(1:lwork),lwork,info)

if( info /=0 ) then
  if( ifop >=1 ) write(lud,*)'ceq_red_con: svd failed, info= ', info
  iret = -3
  return
endif

rank_def=0     !  determine rank deficiency of B
do j=1,nc
    if( S(j)<S(1)*thresh ) then
        rank_def=nc-j+1
        exit
    endif
end do

if( rank_def>0 .and. ifop>=3 ) then
    write(lud,*)' '
    write(lud,*)'ceq_red_con: basic constraint matrix is singular'
    write(lud,*)'rank deficiency= ', rank_def
endif

nb=nc-rank_def

!  identify all determined species
!  species k is determined if U(k,nb+1:ns) = 0
sp_det=0
sp_order=0
kk=0
do k=1,ns
    if( ceq_norm( ns-nb, U(k,nb+1:ns) ) <thresh ) then
        sp_det(k)=1	! species k is determined
        kk=kk+1
        sp_order(kk)=k	! determined species are first in ordering
    endif
end do
nsd=sum(sp_det)
nsu=ns-nsd

do k=1,ns
   if( sp_det(k)==0 ) then
      kk=kk+1
      sp_order(kk)=k	! undetermined species are last in ordering
   endif
end do

if( ifop>=3 ) then
    write(lud,*)' '
    write(lud,'(a,i5)')'nsd = ', nsd
    write(lud,'(a,i5)')'nsu = ', nsu
    if( ifop>=4 ) then
        write(lud,*)'sp_order = '
        write(lud,'((20i4))')sp_order
    endif
endif

!  form DEBg = [D E Bg]
ndebg=nsd+ne+ng
DEBg=0.d0
DEBg(1:ns,nsd+1:ndebg)=B(1:ns,ncs+1:nc)  ! set [E Bg]
do kk=1,nsd
    k=sp_order(kk)
    DEBg(k,kk)=1.d0   ! set [D]
    DEBg(k,nsd+1:ndebg)=0.d0  ! set rows of [E Bg] to zero for determined species
end do

! determine independent columns of DEBg
call ceq_ind_col(ns,ndebg,nsd,DEBg(1:ns,1:ndebg),thresh,indcol(1:ndebg),info)

if( info < 0 ) then
   if( ifop >=1 ) write(lud,*) 'ceq_red_con: ceq_ind_col failed; info = ', info
   iret = -4
   return
endif

!  identify all determined elements
el_det=0
el_order=0
neu=0
ned=0

do k=1,ne
    if( indcol(nsd+k)==0 ) then	!	element k is determined
        ned=ned+1
        el_det(k)=1
        el_order(ned)=k	! determined elements are first in ordering
    else
        neu=neu+1			! element is undetermined
    endif
end do

kk=0
EU=0.d0
do k=1,ne
   if( el_det(k)==0 ) then
      kk=kk+1
      el_order(kk+ned)=k		! undetermined elements are last in ordering
      EU(:,kk)=Ein(:,k)
   endif
end do

if( ifop>=3 ) then
    write(lud,*)' '
    write(lud,*)'ned = ', ned
    write(lud,*)'neu = ', neu
    if( ifop >=4 ) then
        write(lud,*)'el_order = '
        write(lud,'((20i4))')el_order
    endif
endif

!  identify any linearly dependent columns of Bg
if( ng > 0 ) then
   ngi = sum(indcol(nsd+ne+1:ndebg))
else
   ngi = 0
endif

!  assemble BS (BR prior to species re-ordering)
nrc=neu+ngi		! number of reduced constraints
BS=0.d0
kk=0
do k=1,ne
   if( el_det(k)==0 ) then
      kk=kk+1
      BS(:,kk)=Ein(:,k)
   endif
end do

do k=1,ng
   if( indcol(nsd+ne+k)==1 ) then
      kk=kk+1
      BS(:,kk)=Bg(:,k)
   endif
end do

!  assemble BR
BR=0.d0
do kk=1,nsu
   k=sp_order(nsd+kk)
   BR(kk,1:nrc)=BS(k,1:nrc)
end do

!  assemble E
E=0.d0
do is=1,ns
   js=sp_order(is)
   do ie=1,ne
      je=el_order(ie)
      E(is,ie)=Ein(js,je)
   end do
end do

!  determine modified constraint transformation matrix

call ceq_reorder(ns,nb,U(1:ns,1:nb),sp_order,PU(1:ns,1:nb))

BTTPU=0.d0
if( nsd > 0 ) BTTPU(1:nsd,1:nb)=PU(1:nsd,1:nb)
if( nsu > 0 ) BTTPU(nsd+1:nb,1:nb)=matmul(transpose(BR(1:nsu,1:nrc)), PU(nsd+1:ns,1:nb) )
do j=1,nb
    BTTPU(1:nb,j)=BTTPU(1:nb,j)/S(j)
end do
A=0.d0
A(1:nb,1:nc)=matmul(BTTPU(1:nb,1:nc),VT)


!  set to zero near-zero components of A
do j=1,nb
    do k=1,nc
        if( abs(A(j,k))<thresh ) A(j,k)=0.d0
    end do
end do

if( ifop<3 ) return

write(lud,*)'number of elements,              ne= ' ,ne
write(lud,*)'number of species ,              ns= ' ,ns
write(lud,*)'number of constrained species , ncs= ' ,ncs
write(lud,*)'number of general constraints , ng=  ' ,ng
write(lud,*)'number of determined elements , ned= ' ,ned
write(lud,*)'number of determined species ,  nsd= ' ,nsd
write(lud,*)'  '
write(lud,*)'number of undetermined species ,  nsu= ' ,nsu
write(lud,*)'number of undetermined elements , neu= ' ,neu
write(lud,*)'number of indep. gen. constr.,    ngi =' ,ngi
write(lud,*)'number of reduced constraints,    nrc= ' ,nrc
write(lud,*)' '

return

contains	!----------------------------------------------------

subroutine check_input	
  integer e(5), i

  e=0

  if( ns < 1 ) e(1)=1
  if( ne < 1 ) e(2)=1
  if( ncs < 0 .or. ncs > ns ) e(3)=1
  if( ng < 0 ) e(4)=1
  if( nc /= ne+ncs+ng ) e(5)=1

  if( maxval(e) == 0 ) return

  if( ifop >=1 ) then
     write(lud,*)'ceq_red_con: error in input'
     do i=1,5
       if( e(i) >0 ) write(lud,'(a,i3)') '   error in argument ', i
     end do
  endif
 
  iret = -1
  return

end subroutine check_input
end subroutine ceq_red_con