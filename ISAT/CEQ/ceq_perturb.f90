!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ceq_perturb( ns, nsd, nsu, ne, ned, neu, nrc, &
                        zd, cr, BR, E, ifop, lud, &
						eps_el, eps_sp, pert_tol, pert_skip, &
                        zdp, zup, zupper, crp, npert, max_pert, iret)
!------------------
implicit none
integer, intent(in) :: ns, nsd, nsu, ne, ned, neu, nrc, ifop, lud, pert_skip
real(kind(1.d0)), intent(in) :: zd(nsd), cr(nrc), BR(nsu,nrc), E(ns,ne), eps_el, eps_sp, pert_tol
integer, intent(out) :: npert, iret
real(kind(1.d0)), intent(out) :: zdp(nsd), zup(nsu), zupper(nsu), crp(nrc), max_pert
!------------------

!  Generate (possibly) perturbed CE problem
!  Input:
!	ns		- number of species
!	nsd		- number of determined species
!	nsu		- number of undetermined species
!	ne		- number of elements
!	ned		- number of determined elements
!	neu		- number of undetermined elements
!	nrc		- number of reduced constraints
!   zd      - moles of determined species
!   cr      - reduced constraint vector
!   BR      - reduced constraint matrix
!   E       - element matrix
!   ifop    >0 for output
!   lud     - logical unit for diagnostic output
!   eps_el  - relative lower bound on element moles
!   eps_sp  - relative lower bound on species moles
!   pert_tol- largest allowed perturbation (moles/moles of atoms)
!  pert_skip>0 to skip perturbing undetermined species

!  Output:
!   zdp     - perturbed moles of determined species
!   zup     - min-max solution for undetermined species
!   zupper  - upper bound on undetermined species moles
!   crp     - perturbed reduced constraint vector
!   npert   - number of perturbations made
!   max_pert- largest normalized perturbation made
!   iret    =  0  successful operation
!           = -1, non-realizable large perturbation made
!           = -2, failed to determine max-min composition
!           = -3, zero atoms

integer :: j, k
real(kind(1.d0)) :: zdlim, zatoms, &
  cre(neu), sumcre, cref, zed(ned), zeu(neu), ze(ne), zeu_in(neu), zau, zelow, &
  zemax, zumm(nsu), zumin, zulow, zeumax

iret=0         ! anticipate success
npert=0        ! number of perturbations
max_pert=0.d0  ! largest normalized perturbation

if( ifop>=5 ) then
    write(lud,*)'ceq_perturb'
    write(lud,'(a,9i4)')'ne ned neu ns nsd nsu nrc= ', ne, ned, neu, ns, nsd, nsu, nrc
endif

zdp=zd       ! check that determined species are non-negative
zatoms=0.d0  ! estimate of moles of atoms
if( nsd > 0 ) then
  do j = 1, ne
    zatoms=zatoms + abs(dot_product( zd, E(1:nsd,j) ) )
  end do
endif

do j = 1, neu
  zatoms=zatoms + abs( cr(j) )
end do

if( zatoms <= 0.d0 ) then
   if( ifop >=1 ) write(lud,*) 'ceq_perturb: no atoms'
   iret = -3
   return
endif

zdlim=zatoms*pert_tol
do k=1,nsd
    if( zdp(k)<0.d0 ) then
        if( abs(zdp(k))>zdlim   ) then ! significantly negative
            max_pert=max(max_pert, abs(zdp(k))/zatoms)
            npert=npert+1
            if( ifop>=1 ) write(lud,*)'ceq_perturb: negative determined species '
        endif
        zdp(k)=0.d0
    endif
end do

crp=cr
cre=cr(1:neu) ! check that undetermined elements are positive
sumcre=0.d0
do j=1, neu
  sumcre=sumcre+abs(cre(j))
end do
cref=sumcre*pert_tol

do k=1,neu
    if( cre(k)<0.d0 ) then
        if( abs(cre(k))>cref  ) then ! significantly negative
            max_pert=max(max_pert, abs(cre(k))/zatoms )
            npert=npert+1
            if( ifop>=1 ) write(lud,*)'ceq_perturb: negative undetermined element '
        endif
        cre(k)=0.d0
    endif
end do
crp(1:neu)=cre
            
zed=matmul(zdp,E(1:nsd,1:ned) )                !  moles of determined elements
zeu=matmul(zdp,E(1:nsd,ned+1:ne) ) + crp(1:neu)  ! moles of undetermined elements
ze(1:ned)=zed
ze(ned+1:ne)=zeu    ! moles of elements
zeu_in=zeu
zemax=maxval(ze)
zeumax=maxval(zeu)
zau=sum(zeu)       ! moles of atoms in undetermined species

!  impose lower bound on moles of undetermined elements
zelow=max(eps_el*zeumax, eps_el**2*zemax)
do j=1,neu      
   if( zeu(j)<zelow ) then
      npert=npert+1
      max_pert=max(max_pert, (zelow-zeu(j))/zatoms)
      zeu(j)=zelow
   endif
end do

if( ifop>=5 .and. npert>0 ) then
   write(lud,*) 'zeu_in'
   write(lud,'(1p,5e13.4)') zeu_in
   write(lud,*) 'zeu'
   write(lud,'(1p,5e13.4)') zeu
   write(lud,*) 'zeu - zeu_in'
   write(lud,'(1p,5e13.4)') zeu-zeu_in
endif

zupper=0.d0    ! determine upper bound on undetermined species
do j=1,nsu
    zupper(j)=1/maxval(E(nsd+j,ned+1:ne)/zeu)
end do

if( ifop>=5 ) then
   write(lud,*) 'log10(zupper)= '
   write(lud,'(1p,5e13.4)') log10(zupper)
endif

!  determine max-min moles of undetermined species

call ceq_maxmin(nsu,nrc,BR,crp,zumm,zumin,iret)

if( iret<0 ) then
    if( ifop>=1 ) write(lud,*)'ceq_perturb: ceq_maxmin failed, iret=', iret
    iret=-2 
	return 
endif

zup=zumm
do j=1,nsu     ! impose lower limit on undetermined species
   zulow=eps_sp*zupper(j)
   if( zumm(j)<zulow ) then
      zup(j)=zulow
      npert=npert+1
      max_pert=max(max_pert, (zup(j)-zumm(j))/zatoms)
   endif
end do

if( pert_skip > 0 ) return  ! do not modify constraints

!  modify constraints according to the perturbation in undetermined species
crp=crp+matmul(zup-zumm, BR)

if( max_pert > pert_tol ) then
    if( ifop>=1 ) then
        write(lud,*)'ceq_perturb: large perturbation made'
        write(lud,'(a,i5,1p,e13.4)') 'npert, maxpert = ', npert, max_pert
    endif
    iret=-1 
endif

return
end subroutine ceq_perturb