!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module ci_cem_recon

!  Contains subroutine ci_ceq, which determines the 
!  constrained-equilibrium compositions on boundary facets of 
!  the reduced realizable region (and in the interior) by
!  invoking CEQ.

use ci_dat8
use isat_abort_m
use ceq_state_m, only: sys_type, ceq_state
use ci_stats
use ci_utils

implicit none

integer, save :: k_param = 1 ! set =1 to pass non-default parameters to CEQ
                             ! see call to ci_ceq

type :: facet_pointer
   type (facet_type),  pointer :: f
end type facet_pointer

type :: facet_type
   type (sys_type), pointer   :: sys
end type facet_type

logical,                      save :: initialized = .false.
type(facet_pointer), pointer, save :: facet(:)
integer, allocatable,         save :: nsp(:), nrp(:), zp2z(:,:), rp2r(:,:)
integer,                      save :: lu_diag = -2

contains !--------------------------------------------------

subroutine ceq_nrs( nrsp, CSp, rin, hin, press, z_CE, T_CE, stats, info )

  ! Routine to determine constrained-equilibrium composition. This
  ! routine resets CEQ at every call, so use it only if the no. of
  ! constrained species or (nrsp) or the index of constrained species
  ! changes (CSp) changes.

  ! Input:
  !  nrsp  : no. of constrained/represented species
  !  CSp   : index of constrained species
  !  rin   : reduced composition (size nrsp + ne)
  !  hin   : enthalpy
  !  press : pressure

  ! Output:
  !  z_CE  : constrained-equilibrium composition
  !  T_CE  : temperature at z_CE
  !  stats : statistics returned from CEQ
  !  info  : >=0 for success (as returned from CEQ)
  !          < 0  for failure (as returned from CEQ)

    use ceq_system
    use ci_dat8

    implicit none

    integer, intent(in) :: nrsp, CSp(nrsp)
    real(kind(1.d0)), intent(in) :: rin(nrsp+ne), hin, press
    real(kind(1.d0)), intent(out) :: z_CE(ns), T_CE
    integer, intent(out) :: info
    real(kind(1.d0)), intent(out) :: stats(20)
    
    ! local
    integer, parameter :: nbg = 0, diag = -1
    real(kind(1.d0)), parameter:: eps_sp=1.d-9, srat_lim=1.d-9, &
         dec_min=0.5d0, pert_tol=1.d-4, res_tol=1.d-9
  
    integer, save :: init = 0
    integer:: i, luout = 0
    real(kind(1.d0)) :: Bg(1,1), HoR, CEup(ns,ne)
    type (sys_type), pointer :: sys
    
    if ( init /= 0 ) then
       call ceq_sys_rm2(nrsp, 0, sys)
    endif
    
    ! unrepresented element matrix
    CEup = CE
    do i = 1, nrsp
       CEup(CSp(i),1:ne) = 0
    enddo
    
    call ceq_sys_init(ns, ne, nrsp, 0, CEup, CSp, Bg, thermo_ns, luout, &
         diag, sys, info)
    
    call ceq_param_set( sys, diag=diag, eps_sp=eps_sp, srat_lim=srat_lim, & 
         T_low=tbadlo, T_high=tbadhi, dec_min=dec_min, pert_tol=pert_tol, & 
         res_tol=res_tol )
    
    HoR = hin/gascon
    call ceq_state(sys, p_cgs=press, c=rin, HoR=HoR, T_eq=T_CE, N_eq=z_CE, & 
         stats=stats, info=info)
    
    init = 1
    
end subroutine ceq_nrs

subroutine sp_ordering( ncs, nc, CS, nsd, nsu, sp_order )
  
  ! Routine to identify determined and undetermined species using CEQ

  ! Input:
  !  ncs : no. of constraned species = nrs
  !  nc  : no. of constraints = nrs + ne
  !  CS  : index of constrained species

  ! Output:
  !  nsd : no. of determined species
  !  nsu : no. of undetermined species
  !  sp_order : species order, 
  !        first nds determined, followed by nus undetermined
  

  integer, intent(in)  :: ncs, nc, CS(ns)
  integer, intent(out) :: nsd, nsu, sp_order(ns)

  ! local
  integer :: ned, neu, el_order(ne), nrc, nb, iret
  real(kind(1.d0)) :: Bg(1,1), Eout(ns,ne), BR(ns, nc), A(nc, nc), &
       Bout(ns,nc)

  call ceq_red_con( ns, ne, ncs, 0, ne + ncs, CE, CS, Bg, -1, 0, &
       Bout, ned, neu, el_order, nsd, nsu, sp_order, Eout, nrc, nb, BR, &
       A, iret )
       
  
end subroutine sp_ordering
       
subroutine ci_ceq( kfa, r, h, p, guess, z_g, T_g, z, T, & 
     stats, info, ludiag )

!  Determine constrained-equilibrium composition on a facet.

! Input:
!   kfa    index of facet (see note 1 below)
!   r      reduced composition {zr, zue}
!   h      enthalpy (cgs)
!   p      pressure (cgs)
!   guess  if .true., use initial guess (z_g, T_g) (see note 2 below)
!   z_g    initial guess of species specific moles 
!   T_g    initial guess of temperature 
!   ludiag (optional) treatment of diagnostic output (see note 3 below)

! Output:
!   z      specific moles of CE composition
!   T      temperature of CE composition
!   stats  statistics returned from CEQ
!   info   >=0 for success ( = number of Newton iterations)
!          <0 for failure (as returned from CEQ)
!          (see note 4)

! Notes:
! 1/ kfa=0 indicates the interior of the reduced realizable region;
!    for 1 <= kfa <= nrc, the facet is defined by r(kfa)=0.

! 2/ if guess==.false., then z_g and T_g are not referenced.

! 3/ by default, diagnostic output is generated and written to the file
!    ci_ceq.op.  Optionally, on the first call to facet_eq,
!    specify ludiag=-2 to suppress output, or specify ludiag >= 0 to 
!    write diagnostic output to unit ludiag.

! 4/ Treatment of errors: for info in the range -10:-1, the input is invalid
!    and ci_ceq aborts.  For info in the range -17:-11 CEQ failed (possibly
!    with valid input), a diagnostic is generated (unless ludiag=-2) and 
!    ci_ceq returns.  For info = -11 the constraints are unrealizable, and
!    ci_ceq returns.  For info = -18 or -19, T_eq is out of the specified range,
!    and ci_ceq returns.  When facet_ceq returns with info < 0, z and T are
!    set to zero.  

integer,    intent(in)  :: kfa
logical,    intent(in)  :: guess
real(k_dp), intent(in)  :: r(nrc), h, p, z_g(ns), T_g
integer,    intent(out) :: info
real(k_dp), intent(out) :: z(ns), T, stats(20)
integer, optional, intent(in) :: ludiag

integer    :: i
real(k_dp) :: HoR, zp(ns), rp(nrc), z_gp(ns)

! stats
call routine_start(i_ci_ceq)

! on first call, initialize data structures and output file
! for diagnostics

if( .not.initialized ) then
   if( present(ludiag) ) lu_diag = ludiag
   call ci_ceq_init
endif

!  check index of facet

if( kfa < 0  .or.  kfa > nrc  ) call isat_abort( 'ci_ceq', 1, &
                                mess='invalid kfa = ', isv = kfa )

!  initialize facet, if not already initialaized
if( nsp(kfa) == 0 ) call ci_ceq_init_facet(kfa)

! form rp = {zr, zue} for positive species and unrepresented elements
do i = 1, nrp(kfa)
   rp(i) = r( rp2r(i,kfa) )
end do

HoR = h/gascon
z   = 0.d0  !  null values returned in case of failure
T   = 0.d0

!  perform CEQ calculation
if( .not.guess ) then
!  call ceq_state without initial guess
   call ceq_state( facet(kfa)%f%sys, c=rp(1:nrp(kfa)), p_cgs=p, HoR=HoR,   &
                   N_eq=zp, T_eq=T, stats=stats, info=info )
                   
else    
!  call ceq_state with initial guess
   do i = 1, nsp(kfa)
      z_gp(i) = z_g(zp2z(i,kfa))  ! initial guess of positive species
   end do
             
   call ceq_state( facet(kfa)%f%sys, c=rp(1:nrp(kfa)), p_cgs=p, HoR=HoR, &
                   N_g=z_gp, T_g=T_g,   &
                   N_eq=zp, T_eq=T, stats=stats, info=info )
endif
   
!  treat failures
if( info < 0 ) then
   if( info >= -10 ) then
      call isat_abort( 'ci_ceq', 2, mess='CEQ failure, info = ', isv = info ) 
      
   elseif( info <= -17 ) then
      if( lu_diag >= 0 ) write(lu_diag,'(a,2i4)') 'CEQ failure, kfa, info = ', &
                                                   kfa, info
   endif
   !return XXX varun 
endif
    
! success: map positive species to all species
do i = 1, nsp(kfa)
   z( zp2z(i,kfa) ) = zp(i)
end do

call routine_stop(i_ci_ceq)
                        
return
end subroutine ci_ceq

subroutine ci_ceq_init  !-----------------------------------------------

!  Allocate facet arrays and pointers, and set lu_diag

integer :: i

allocate( facet(0:nrc) )

do i = 0, nrc
   allocate( facet(i)%f )
end do

allocate( nsp(0:nrc) )
allocate( nrp(0:nrc) )
allocate( zp2z(ns,0:nrc) )
allocate( rp2r(nrc,0:nrc) )

nsp = 0  ! indicate facet not initialized

if( lu_diag == -1 ) then
   call isat_lu( lu_diag )
   open( lu_diag, file = 'ci_ceq.op' )
endif

initialized = .true.

return

end subroutine ci_ceq_init

subroutine ci_ceq_init_facet( kfa )  !-----------------------------------

!  Perform initialization of the CEQ SYS for facet kfa

!  Conceptually, the input is: 
!      kfa, ns, ne, nrs, CS, CE, thermo_ns

!  On facet kfa, the species and unrepresented elements which are
!  positive are determined and, based on these, the following 
!  corresponding quantities are defined:
!      nsp, nep, nrsp, CSp, CEp, thermo_p
!  Note that CEp is for unrepresented species only.

use ci_dat, only: luout
use ceq_system
integer, intent(in) :: kfa

integer    :: i, j, s_pos, s_rep, c_pos, isz, nep, nrsp, ep2e(ne), &
              CSp(ns), diag = 0, iret, z2zp(ns)
real(k_dp) :: BA(ns+1,nrc+2), CEp(ns,ne), thermo_p(ns,15), Bg(1,1), &
              eps_sp, srat_lim, dec_min, pert_tol, res_tol, &
              T_low, T_high, T_tol

!  The matrix BA is formed as follows
!  Initially BA(1:ns,1:nrc) = BBF. 
!  Rows BA(i,1:nrc) are set to zero for all species i which are zero on the facet.

!  BA(i,s_pos) is an indicator for species i being positive.    (s_pos = nrc+1) 
!  BA(i,s_rep) is an indicator for species i being represented. (s_rep = nrc+2) 

!  BA(c_pos,j) is an indicator for column j being positive      (c_pos = ns+1)

!  Initialize BA

s_pos = nrc+1
s_rep = nrc+2
c_pos = ns+1

BA              = 0.d0
BA(1:ns,s_pos)  = 1.d0    ! all species initially positive
BA(c_pos,1:nrc) = 1.d0    ! all reduced compositions initially positive

BA(1:ns,nrs+1:nrc) = CE      ! elements

do i = 1, nrs  !  represented species
   BA( CS(i), i )        = 1.d0
   BA( CS(i), s_rep )    = 1.d0
   BA( CS(i), nrs+1:nrc) = 0.d0 ! set represented elements to zero 
end do

if( lu_diag >= 0 ) then
   write(lu_diag,'(a,i4)')'Initialization for facet kfa = ', kfa
   write(lu_diag,*)'BA prior to modification'
   write(lu_diag,*)' '
   do i = 1, ns+1
      write(lu_diag,'(40i3)') i, (nint(BA(i,j)),j=1,nrc+2) 
   end do
endif

!  Modify BA based on zero represented composition on facet kfa
!  (no modification for interior, kfa=0) 


if( kfa >=1  .and.  kfa <= nrs ) then 
!  treat case of a represented species being zero

   isz             = CS(kfa)  !  index of zero species
   BA(1:ns+1,kfa)  = 0.d0     ! set column to zero
   BA(isz,1:nrc+1) = 0.d0     ! set row to zero
   
elseif( kfa > nrs  .and.  kfa <= nrc ) then
!  treat case of unrepresented element being zero

   do i = 1, ns  !  loop over species
      if( BA(i,kfa) > 0.d0  .and.  BA(i,s_rep) == 0.d0 ) then
         BA(i,1:nrc+1) = 0.d0  ! unrepresented species is zero
      endif
   end do

elseif( kfa < 0  .or.  kfa > nrc ) then 
   call isat_abort('ci_ceq_init_facet', 1,  &
                    mess='bad value of kfa = ', isv = kfa )  
endif

!  check for zero unrepresented elements
do i = nrs+1, nrc
   if( maxval( BA(1:ns,i) ) == 0.d0 ) BA(c_pos,i) = 0.d0
end do

if( lu_diag >= 0 ) then
   write(lu_diag,*)' '
   write(lu_diag,*)'BA after to modification'
   do i = 1, ns+1
      write(lu_diag,'(40i3)') i, (nint(BA(i,j)),j=1,nrc+2) 
   end do
endif

!  Modification to BA is complete: extract needed information

!  Construct mappings from positive to all

j = 0  ! zp2z - positive species to all species
do i = 1, ns
   if( BA(i,s_pos) > 0.d0 ) then
      j           = j+1
      zp2z(j,kfa) = i
      z2zp(i)     = j
   endif
end do
nsp(kfa) = j

nep = 0  ! ep2e - positive unrepresented elements to unrepresented elements
do i = 1, ne
   if( BA(c_pos,nrs+i) > 0.d0 ) then
      nep       = nep+1
      ep2e(nep) = i
   endif
end do

j = 0  ! rp2r - positive reduced composition to reduced composition
do i = 1, nrc
   if( BA(c_pos,i) > 0.d0 ) then
      j       = j+1
      rp2r(j,kfa) = i
   endif
end do
nrp(kfa) = j

nrsp = 0  !  positive represented species
do i = 1, nrs
   if( BA(c_pos,i) > 0.d0 ) then
      nrsp = nrsp+1
      CSp(nrsp) = CS( z2zp(i) )
   endif
end do

do i = 1, nsp(kfa) 
   do j = 1, nep   !   unrepresented element matrix
      CEp(i,j) = BA( zp2z(i,kfa), nrs+ep2e(j) )
   end do
                   !    thermo
   thermo_p(i,:) = thermo_ns( zp2z(i,kfa), : )
end do

if( lu_diag >= 0 ) call ci_ceq_op

!  initialize CEQ system for facet

call ceq_sys_init( nsp(kfa), nep, nrsp, 0, CEp(1:nsp(kfa),1:nep), CSp, Bg, &
                  thermo_p(1:nsp(kfa),1:15), luout, diag, &
                  facet(kfa)%f%sys, iret )
                   
if( iret < 0 ) call isat_abort('ci_ceq_init_facet', 2,  &
                    mess='ceq_sys_init failed, iret = ', isv = iret )
  
!  CEQ defaults                       
diag     = 0
eps_sp   = 1.d-9
srat_lim = 1.d-9
dec_min  = 0.5d0
pert_tol = 1.d-4
res_tol  = 1.d-9 
T_tol    = 1.d-6

! changes to defaults
!XX diag     = 5 !XXX
res_tol  = 1.d-10 !XXX 
T_tol    = 1.d-7  !XXX

T_low    = tbadlo  !  make allowed range consistent
T_high   = tbadhi

!  as necessary, change parameters here

call ceq_param_set( facet(kfa)%f%sys, diag=diag, eps_sp=eps_sp, &
                    srat_lim=srat_lim, T_low=T_low,T_high=T_high, &
                    dec_min=dec_min, pert_tol=pert_tol, res_tol=res_tol )

return

contains !-------------------------------------------------

subroutine ci_ceq_op

write(lu_diag,*)' '
write(lu_diag,*)'kfa, nsp, nep, nrp, nrsp = '
write(lu_diag,'(10i5)')kfa, nsp(kfa), nep, nrp(kfa), nrsp 

write(lu_diag,*)' '
write(lu_diag,*)'zp2z = '
write(lu_diag,'(20i5)') (zp2z(i,kfa),i=1,nsp(kfa))

write(lu_diag,*)' '
write(lu_diag,*)'ep2e = '
write(lu_diag,'(20i5)') (ep2e(i),i=1,nep)

write(lu_diag,*)' '
write(lu_diag,*)'rp2r = '
write(lu_diag,'(20i5)') (rp2r(i,kfa),i=1,nrp(kfa))

write(lu_diag,*)' '
write(lu_diag,*)'CSp = '
write(lu_diag,'(20i5)') (CSp(i),i=1,nrsp)

write(lu_diag,*)' '
write(lu_diag,*)'CEp = '
do i = 1, nsp(kfa)
   write(lu_diag,'(20i5)') (nint(CEp(i,j)),j=1,nep)
end do
write(lu_diag,*)' '

end subroutine ci_ceq_op

end subroutine ci_ceq_init_facet

end module ci_cem_recon
