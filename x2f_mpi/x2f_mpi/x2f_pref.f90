!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the x2f_mpi      !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module x2f_pref
!
! "pref" here means:
!  Particles have preference to go to the processors those they have not been to
!  before. In another word, the distribution strategy tries to distribute particles
!  based on their previous steps' (history) information. This means if the particle
!  has been to the particular processor but not been evaluated, then it can not be
!  assigned to the processor again. The strategy tries to make the work load balanced,
!  at the same time distriute particles as possible as it can.
!  
use MPI  ! MPI coding indicated by XXX in comments below
!
! Usage:
!   in calling routine, include the statements:
!     use x2f_pref
!     type (pref_type), save :: pref
!
!   call x2f_pref_init(pref,...) to initialize 
!   call x2f_pref_p(pref,...) to assign x-vectors to processors
!   call x2f_updat(pref,... ) to update the arrays containing history information
!
implicit none
type :: pref_type  !--------------------------------------------------
!
! nbat    - number of batches 
! npar    - npar(1:nbat),number of particles in each batch
! hist    - hist(igat,nbat), save the processor indices the batch has been to
! b       - b(1:sum(npar)), b(ip)=jb indicates the particle ip belong to batch jb.
!           Particles being to the same processors are grouped into one batch.
!
    integer          :: nbat
    integer, pointer :: npar(:), hist(:,:), b(:)

end type pref_type

integer, allocatable, save, target :: npar0(:), hist0(:,:), b0(:)

contains  !-------------------------------------------------------------

subroutine x2f_pref_init( pref, myrank, nv, nolocal )
!
! Initialize data structure pref
! 
! Input:
!    myrank   - the rank of the processor
!    nv       - number of vector
!    nolocal  - try local processor or not
! Output:
!    pref     - data structure being initialized
!
integer, intent(in) :: myrank, nv, nolocal
type(pref_type), intent(inout) :: pref

if ( allocated(npar0) )  deallocate( npar0 )
if ( allocated(hist0) )  deallocate( hist0 )
if (  allocated(b0)   )  deallocate(  b0   )
allocate(  npar0(1)  )
allocate( hist0(1,1) )
allocate(   b0(nv)   )

pref%npar => npar0
pref%hist => hist0
pref%b    => b0

! Initially every processor has just one batch and all the particles on the 
! processor belong to the batch
pref%npar(1) = nv
pref%nbat    = 1
pref%b(1:nv) = 1

! Depending on whether to try ISAT table on local processor or not, set the 
! initial pref%hist array. If do not want to try local processor, by setting
! pref%hist(1,1)=myrank+1 indicates the processor has tried local processor; else 
! setting pref%hist(1,1)=0 indicates none of the processor (within given 
! communicator) has tried.
if ( nolocal == 1 ) then
   pref%hist(1,1) = myrank + 1
else
   pref%hist(1,1) = 0 
end if

end subroutine x2f_pref_init

!-------------------------------------------------------------

subroutine x2f_pref_p( nproc, myrank, mpicomm, igat, pref, p, nv, &
                       n_incoming, n_outgoing )
!
! Assign particles to processors
!
! Input:
!    mpicomm  - communicator used  
!    igat     - the index of attempt
!    nv       - number of unevaluated particle
!
! Output:
!    p(i)     - indicate i-th particle goes to p(i) processor 
! n_incoming  - size of nproc, indicate the number of particle expected from every
!               processor
! n_outgoing  - size of nproc, indicate the number of particle left (from local 
!               processor) for every processor
!
type (pref_type)            :: pref
integer, intent(in)         :: nproc, myrank, mpicomm, igat, nv
integer, intent(inout)      :: p(nv), n_incoming(nproc), n_outgoing(nproc)

real, parameter             :: para = 0.08
integer, allocatable :: npar(:), npar_array(:), Mb0(:), Mb(:), hist_array(:,:)
integer, allocatable :: bat_index(:,:), bat_start(:)
integer   :: nbat, nbat_array(nproc), nbat_tot, aveN, Vp(nproc)
integer   :: a, left, ierr, ibstart, ibend, ip, ib, tempn
integer   :: disp(nproc), ipstart, stmker, recv(nproc), maxN_grp, minN_grp

n_incoming(1:nproc) = 0
n_outgoing(1:nproc) = 0
nbat = pref%nbat

if ( igat==1 .and. pref%hist(1,1)==0 ) then

call MPI_Allreduce( nv, maxN_grp, 1, MPI_INTEGER, MPI_MAX, mpicomm, ierr )
call MPI_Allreduce( nv, minN_grp, 1, MPI_INTEGER, MPI_MIN, mpicomm, ierr )
!print *, 'ratio:', maxN_grp/(1.0*minN_grp)
if ( minN_grp>0.and.real(maxN_grp/minN_grp)<2.0 ) then
   n_incoming(myrank+1) = nv
   n_outgoing(myrank+1) = nv
   p(1:nv) = myrank + 1
   return
end if  

end if   

! Gather the batch information from each processor
call MPI_Allgather( nbat, 1, MPI_INTEGER, nbat_array, 1, MPI_INTEGER, mpicomm, ierr )
! Number of total batches over all processors
! nbat_array  - size of nproc, indicate the number of batch on each processor
nbat_tot = sum(nbat_array(1:nproc))

allocate(   npar_array(nbat_tot)    )
allocate(      Mb0(nbat_tot)        )
allocate(      Mb(nbat_tot)         )
allocate( hist_array(igat,nbat_tot) )
allocate(  bat_index(nbat_tot,2)    )
allocate(     bat_start(nbat)       )

! bat_index(j,1)=ip : global j-th batch belongs to ip processor
! bat_index(j,2)=i  : global j-th batch corresponds to local ith batch 
ibstart = 1
do ip = 1, nproc
   ibend = sum(nbat_array(1:ip))
   do ib = ibstart, ibend
      bat_index(ib,1) = ip
      bat_index(ib,2) = ib-ibstart+1
   end do
   ibstart = 1 + ibend
end do

! bat_start(ib) = jp  : local ib-th batch starts from jp-th particle
if ( nbat>0 ) then
   bat_start(1) = 1
   if ( nbat>1 ) then
      do ib = 2, nbat
         bat_start(ib) = bat_start(ib-1)+pref%npar(ib-1) 
      end do
   end if
end if

! Gather the number of particles in each batch: npar_array(1:nbat_tot)
! npar_array(nbat_tot)  - indicate number of particles in each batch
disp(1:nproc) = 0
stmker = 0
do ip = 1, nproc
   disp(ip) = stmker
   stmker   = stmker + nbat_array(ip)
   recv(ip) = nbat_array(ip)
end do

call MPI_Allgatherv( pref%npar, nbat, MPI_INTEGER, npar_array, recv, disp, &
                     MPI_INTEGER, mpicomm, ierr )

!Gather the history information from each processor
disp(1:nproc) = 0
stmker = 0
tempn  = 0
do ip = 1, nproc
   disp(ip) = stmker
   tempn    = tempn + nbat_array(ip) 
   stmker   = tempn*igat
   recv(ip) = nbat_array(ip)*igat
end do
!
! hist_array  - store the processor indices each batch has been to before.
!               e.g., hist_array(1:igat,1) stores the processor indices for the 
!               first igat attempts that batch 1 has been to
!
call MPI_Allgatherv( pref%hist, nbat*igat, MPI_INTEGER, hist_array, recv, disp, &
                     MPI_INTEGER, mpicomm, ierr )
                                      
! Average number of particles each processor expects
aveN = sum(npar_array(1:nbat_tot))/nproc + 1
! Mb0  - initially the number of particles in each batch needs to be assigned
Mb0(1:nbat_tot) = npar_array(1:nbat_tot)

a = -1
left = sum(npar_array)
do while ( left>0 )
   a = a + 1
   Mb = Mb0
   ! Vp  - the amount of vacancy available on each processor; if the particles can
   !       not be all distributed, increase the amount of vacancy 
   Vp = aveN*(1+para)**a
   n_incoming(1:nproc) = 0
   n_outgoing(1:nproc) = 0
   call pref_p( nproc,nv,nbat,nbat_tot,igat,myrank,Mb,Vp,hist_array, &
                 bat_index,bat_start,p,n_outgoing,n_incoming,left)
end do
   
do ip = 1, nv
   if ( p(ip)<=0.or.p(ip)>nproc ) stop 'p(ip) error'
end do

deallocate( npar_array )
deallocate(     Mb     )
deallocate(    Mb0     )
deallocate( hist_array )    
deallocate( bat_index  )
deallocate( bat_start  )

return

end subroutine x2f_pref_p

subroutine pref_p( nproc,nv,nbat,nbat_tot,igat,myrank,Mb,Vp,hist_array, &
                   bat_index,bat_start0,p,n_outgoing,n_incoming,sumMb )
!
! This subroutine implements the algorithm of selecting a pair of batch and 
! processor, so that the particles in the batch can be assign to some processor. 
!
integer, intent(in)    :: nproc, nv, nbat, nbat_tot, igat, myrank, bat_start0(nbat), &
                          hist_array(igat,nbat_tot), bat_index(nbat_tot,2) 
integer, intent(inout) :: Mb(nbat_tot), Vp(nproc), p(nv), n_outgoing(nproc), &
                          n_incoming(nproc), sumMb
integer  :: bat_start(nbat)
integer  :: p_sl, b_sl, ip, ib, temp, tempn, iproc, k_me, i, ipstart, p_sl0, &
            Fp(nproc), np
real     :: Fp_min, Ab, Ab_min, Fp0(nproc), Fp1(nproc)
logical  :: cont, beento

bat_start = bat_start0
k_me  = myrank + 1
sumMb = sum(Mb)
cont  = .true.
p     = 0

Fp = 0
do ip = 1, nproc
   do ib = 1, nbat_tot
      temp  = 1
      do i = 1, igat
         if ( hist_array(i,ib)==ip ) then
            temp  = 0
            Exit
         end if
      end do
      Fp(ip) = Fp(ip) + temp*Mb(ib)
   end do
end do

do while ( cont )   
   Fp_min = huge(1.0)
   p_sl   = 0
   do ip = 1, nproc
      if ( Fp(ip)>0.and.Vp(ip)/=0 ) then
         ! Fp0(ip)  - indicates the number of particles per vacancy that could be
         !            assigned to processor ip 
         Fp0(ip) = real(Fp(ip))/Vp(ip)
         if ( Fp0(ip)<Fp_min ) then
            Fp_min = Fp0(ip)
            p_sl   = ip
         end if
      else
         Fp(ip)  = 0
         Fp0(ip) = 0.0
      end if
   end do
      
   if ( p_sl/=0 ) then
      Ab_min = huge(1.0)
      do ib = 1, nbat_tot
         if ( Mb(ib)/=0 ) then
            temp = 1
            do i = 1, igat
               if ( hist_array(i,ib)==p_sl ) then
                  temp = 0
                  Exit
               end if
            end do
            if ( temp>0 ) then
               ! Ab  - indicate the number of available vacancies per particle in
               !       batch ib
               Ab = real(temp)*Vp(p_sl)/Mb(ib)
               if ( Ab<Ab_min ) then
                  b_sl   = ib
                  Ab_min = Ab
               end if
            end if
         end if
      end do
   
      if ( p_sl<=0.or.p_sl>nproc ) stop 'p_sl error'

      ! Make assignment for particles in selected batch
      ! Update Mb(b_sl)(number of particles in b_sl batch still needs to be assign)
      ! Update Vp(p_sl)(number of vacancies remains on processor p_sl)
      ! Update n_outgoing and n_incoming (needed for MPI calls)
      if ( Mb(b_sl)<=Vp(p_sl) ) then
         np = Mb(b_sl)
         iproc = bat_index(b_sl,1)
         if ( k_me==iproc ) then
            ib = bat_index(b_sl,2) 
            ipstart = bat_start(ib)         
            p(ipstart:ipstart+Mb(b_sl)-1) = p_sl
            bat_start(ib) = bat_start(ib)+Mb(b_sl)
            n_outgoing(p_sl) = n_outgoing(p_sl)+Mb(b_sl)
         elseif ( k_me==p_sl ) then
            n_incoming(iproc) = n_incoming(iproc)+Mb(b_sl)
         end if
         Vp(p_sl) = Vp(p_sl)-Mb(b_sl)
         Mb(b_sl) = 0
      else
         np = Vp(p_sl)
         iproc = bat_index(b_sl,1) 
         if ( k_me==iproc ) then
            ib = bat_index(b_sl,2) 
            ipstart = bat_start(ib)         
            p(ipstart:ipstart+Vp(p_sl)-1) = p_sl 
            bat_start(ib) = bat_start(ib)+Vp(p_sl)
            n_outgoing(p_sl) = n_outgoing(p_sl)+Vp(p_sl)
         elseif ( k_me==p_sl ) then
            n_incoming(iproc) = n_incoming(iproc)+Vp(p_sl)
         end if
         Mb(b_sl) = Mb(b_sl)-Vp(p_sl)
         Vp(p_sl) = 0
      end if
   end if
   
   ! Update Fp(1:nproc)
   do ip = 1, nproc
      if ( Fp(ip)/=0 ) then
         beento = .false.
         do i = 1, igat
            if ( hist_array(i,b_sl)==ip ) then
               beento = .true.
               Exit
            end if
         end do
         if ( .not.beento ) Fp(ip) = Fp(ip) - np
      end if
   end do
   
   if ( sum(Mb)<sumMb ) then
      cont = .true.
   else
      cont = .false.
   end if
   sumMb = sum(Mb)
end do

end subroutine pref_p

subroutine x2f_pref_update( pref, igat, nproc, n_outgoing, myrank )
!
! Update "pref" (pref%nbat, pref%npar, pref%hist)
!
type (pref_type),intent(inout) :: pref
integer, intent(in)   :: igat, nproc, n_outgoing(:), myrank
integer, allocatable  :: nparnew(:), histnew(:,:), btemp2(:), tempb(:), &
                         nparnew1(:), tempsrc(:,:)
integer   :: nbat, np, ibat, nbatnew

np = sum(pref%npar)
allocate( btemp2(np) )
btemp2 = 0

! Count the number of batches we have now
call cntbat(np,nproc,n_outgoing,pref,btemp2,nbat,myrank)

allocate(    nparnew(nbat)     )
allocate( histnew(igat+1,nbat) )

! Update history information
call update(nproc,np,igat,nbat,n_outgoing,btemp2,pref,nparnew,histnew,myrank)
nbatnew = nbat
if (nbat==0) then
   pref%nbat = 0
   npar0 = 0
   return
end if

deallocate( npar0  )
deallocate( hist0  )
deallocate(   b0   )
allocate(   npar0(nbatnew)      )
allocate( hist0(igat+1,nbatnew) )
allocate(      b0(np)        )
pref%npar => npar0
pref%hist => hist0
pref%b    => b0

do ibat = 1, nbatnew
   npar0(ibat) = nparnew(ibat)
   hist0(:,ibat) = histnew(:,ibat)
end do

b0    = btemp2
pref%nbat = nbatnew

deallocate( nparnew )
deallocate( histnew )
deallocate(  btemp2 )

end subroutine x2f_pref_update

subroutine cntbat(np,nproc,n_outgoing,pref,btemp2,nbat,myrank)

implicit none
type (pref_type), intent(in) :: pref
integer, intent(in)          :: np, nproc, n_outgoing(nproc), myrank
integer, intent(inout)       :: nbat, btemp2(np)

integer  :: ibstart, ipstart, ipend, btemp1(np)
integer  :: ip, jp, ib, ib_pre
logical  :: found

nbat = 0;    ibstart = 0
ipstart = 1; ipend   = 0
do ip = 1, nproc+1

   ibstart = 1 + nbat
   if ( ip<nproc+1 ) then
      ipend = ipend + n_outgoing(ip)
   else
      ipend = np
   end if

   if ( ipend>=ipstart ) then
      do jp = ipstart, ipend
         if ( pref%b(jp)>0 ) then
            if ( nbat>=ibstart ) then
               found = .false.
               do ib = ibstart, nbat
                  if ( pref%b(jp)==btemp1(ib) ) then
                     found = .true.
                     ib_pre = ib
                     Exit
                  end if
               end do
               
               if ( .not.found ) then
                  nbat = nbat + 1
                  btemp1(nbat) = pref%b(jp)
                  btemp2(jp) = nbat
               else
                  btemp2(jp) = ib_pre       
               end if               
            else
               nbat = nbat + 1
               btemp1(nbat) = pref%b(jp)
               btemp2(jp) = nbat
            end if
         else
            btemp2(jp) = 0    
         end if
      end do
   end if

   ipstart = ipend + 1
end do

end subroutine cntbat

subroutine update(nproc,np,igat,nbat,n_outgoing,btemp2,pref,nparnew,histnew,myrank)

type (pref_type), intent(in) :: pref
integer, intent(in)    :: nproc, np, igat, nbat, myrank
integer, intent(in)    :: n_outgoing(nproc), btemp2(np)
integer, intent(inout) :: nparnew(nbat), histnew(igat+1,nbat)

integer  :: ipstart, ipend, ip, jp, ibat, ib_pre, ib_temp, itest
logical  :: found

! Update history information
nparnew = 0; 
ipstart = 1; ipend   = 0
do ip = 1, nproc + 1
   if ( ip<nproc+1 ) then
      ipend = ipend + n_outgoing(ip)
   else
      ipend = np
   end if
   ib_pre = 0
   if ( ipend>=ipstart ) then
      do jp = ipstart, ipend
         ibat = btemp2(jp)
         if ( ibat/=0 ) then
            nparnew(ibat) = nparnew(ibat)+1
            if ( ibat>ib_pre ) then
               ib_temp = pref%b(jp)
               itest = 1
               found = .false.
               do while (.not.found.and.itest<=igat)
                  if (ip<pref%hist(itest,ib_temp)) then
                     found = .true.
                     itest = itest - 1
                  else 
                     if ( itest==igat ) then
                        found = .true.
                     else
                        if ( ip<pref%hist(itest+1,ib_temp) ) then
                           found = .true.
                        else
                           itest = itest + 1
                        end if
                     end if
                  end if                  
               end do
               if ( itest==0 ) then
                  histnew(1,ibat)  = ip
                  histnew(2:igat+1,ibat)= pref%hist(1:igat,ib_temp) 
               else
                  histnew(1:itest,ibat) = pref%hist(1:itest,ib_temp)
                  histnew(itest+1,ibat) = ip
                  if (igat>1) then
                     histnew(itest+2:igat+1,ibat) = pref%hist(itest+1:igat,ib_temp)
                  end if
               end if
               ib_pre = ibat
            end if
         end if
      end do
   end if
   ipstart = 1+ipend
end do

end subroutine update

end module x2f_pref