!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module ceq_system

   use ceq_state_m
   implicit none

contains

   subroutine ceq_sys_init( ns, ne, ncs, ng, &
        Ein, CS, Bg, thermo_in, lu_op, diag, sys, iret )

!  Specify a constrained equilibrium system.
!  This routine must be called (to generate the data structure  sys)
!  prior to calling ceq_state.
	
!   S.B. Pope 6/30/03

! Input:
!   ns    - number of species
!   ne    - number of elements
!   ncs   - number of constrained species
!   ng    - number of general linear constraints
!   Ein   - element matrix (ns x ne).  A molecule of species k contains
!           Ein(k,j) atoms of element j.
!   CS    - array containing indexes of constrained species (n_cs x 1)
!   Bg    - general linear constraint matrix (ns x ng)
!   thermo_in  - thermodynamic data for species (ns x 15)
!              - see below for details
!   lu_op - logical unit for output 
!   diag  - level of diagnostic output (0=none, 1=severe errors, ...5=full)
!   (The values of lu_op and diag are stored in sys, and used in calls to ceq_state.
!    These values can be changed by, for example, call ceq_param_set( lu_op=3, diag=4 ).)

! Output:
!   sys  - data structure (of type sys_type)
!   iret < 0 for failure

! Unconstrained equilibrium:
!   For unconstrained equilibrium calculations (i.e., in which the only 
!   constraint is on the elements), set ncs=0, ng=0, CS=[] and Bg=[].

! Details of thermo_in:
!   As in Chemkin, the thermodynamic properties of each species are given
!   in terms of 7 non-dimensional coefficients, a(1:7).  Different values 
!   of a(1:7) are given for two temperature ranges.  T* [K] denotes the
!   upper limit of the lower temperature range, and the lower limit of the 
!   upper temperature range.  The contents of the array thermo_in are: 
!   thermo_in(k,1)    = T* for species k
!   thermo_in(k,2:8)  = a(1:7) for species k in the lower temperature range
!   thermo_in(k,9:15) = a(1:7) for species k in the upper temperature range
!   The relevant thermodynamic variables are:
!   T   - temperature [K]
!   H   - molar specific enthalpy [J/kmol]
!   S   - molar specific entropy [J/(kmol K)]
!   R   - universal gas constant [J/(kmol K)]
!   Then for each species, H and S are given by:
!   H/(RT) = a(6)/T + sum_{n=1}^{n=5} a(n) T^{n-1}/n
!   S/R = a(1)log(T) + a(7) + sum_{n=2}^{n=5} a(n) T^{n-1}/(n-1)

!  Notes:
!    The values of lu_op and diag are stored in sys, and used in subsequent calls to ceq_state.
!    These values can be changed by, for example, call ceq_param_set( lu_op=3, diag=4 ).


!------------------------
  implicit none
  integer, intent(in)      :: ns, ne, ncs, ng, CS(ncs), lu_op, diag
  real(k_dp), intent(in)   :: Ein(ns,ne), Bg(ns,ng), thermo_in(ns,15)
  integer, intent(out)     :: iret
  real(k_dp) :: BR(ns,ne+ncs+ng), A(ne+ncs+ng,ne+ncs+ng)
  type (sys_type), pointer :: sys
!------------------------

  integer    :: nc, j, k
  real(k_wp) :: cpu_redcon0, cpu_redcon
    

   allocate(sys)
   call ceq_param_def( sys )  

!  check input --------------------------------------------
   if( ns >=1 ) then
      sys%ns = ns
   else
	  if( diag >=1 ) write(lu_op,*)'ceq_sys_init: bad input, ns= ', ns
	  iret = -1
	  return
   endif

   if( ne >=1 ) then
      sys%ne = ne
   else
	  if( diag >=1 ) write(lu_op,*)'ceq_sys_init: bad input, ne= ', ne
	  iret = -2
	  return
   endif

   if( ncs >=0 ) then
      sys%ncs = ncs
	  if( ncs >= 1 ) then
	     allocate( sys%CS(ncs) )
         sys%CS = CS(1:ncs)
	  endif
   else
	  if( diag >=1 ) write(lu_op,*)'ceq_sys_init: bad input, ncs= ', ncs
	  iret = -3
	  return
   endif

   if( ng >=0 ) then
      sys%ng = ng
	  if( ng >= 1 ) then
	     allocate( sys%Bg(ns,ng) )
         sys%Bg = Bg(1:ns,1:ng)
	  endif
   else
	  if( diag >=1 ) write(lu_op,*)'ceq_sys_init: bad input, ng= ', ng
	  iret = -4 
	  return
   endif

!  allocate arrays
   nc = ne + ncs + ng
   sys%nc = nc
   allocate( sys%el_order(ne) )
   allocate( sys%sp_order(ns) )
   allocate( sys%Ein(ns,ne) )
   allocate( sys%E(ns,ne) )
   allocate( sys%B(ns,nc) )
   allocate( sys%thermo(ns,15) )

   sys%Ein = Ein(1:ns,1:ne)

! form basic and reduced constraint equations

   call ceq_red_con( ns, ne, ncs, ng, nc, Ein, CS, Bg, diag, lu_op, &
	  sys%B, sys%ned, sys%neu, sys%el_order, sys%nsd, sys%nsu, sys%sp_order, &
	  sys%E, sys%nrc, sys%nb, BR, A, iret)

   if( iret < 0 ) return  !  failure

   allocate( sys%BR(sys%nsu,sys%nrc) )
   allocate( sys%A(sys%nb,nc) )

   sys%BR = BR(1:sys%nsu,1:sys%nrc)
   sys%A  = A(1:sys%nb,1:nc)

!re-order thermo
   call ceq_reorder( ns, 15, thermo_in, sys%sp_order, sys%thermo )

   sys%lu   = lu_op
   sys%diag = diag
   sys%initialized = .true.  ! indicate successful initialization

   return
   end subroutine ceq_sys_init  !-------------------------------------------------------------

   subroutine ceq_sys_rm( sys )

!  Remove data structure sys

   type (sys_type), pointer :: sys

   nullify( sys%CS, sys%sp_order, sys%el_order, sys%Ein, sys%E, sys%Bg, sys%B, &
            sys%BR, sys%A, sys%thermo ) 

   nullify(sys)

   return
   end subroutine ceq_sys_rm  !--------------------------------------------------------------

!-------------------------------------------
    subroutine ceq_sys_rm2( ncs, ng, sys )

!  Laniu free the memory
     integer :: ncs, ng
    type (sys_type), pointer :: sys
    
!    if(allocated(sys%CS)) deallocate(sys%CS) 
!    if(allocated(sys%sp_order)) deallocate(sys%sp_order)
!    if(allocated(sys%el_order)) deallocate(sys%el_order) 
!    if(allocated(sys%Ein)) deallocate(sys%Ein)
!    if(allocated(sys%E)) deallocate(sys%E)
!    if(allocated(sys%Bg)) deallocate(sys%Bg)
!    if(allocated(sys%B)) deallocate(sys%B)
!    if(allocated(sys%BR)) deallocate(sys%BR)
!    if(allocated(sys%A)) deallocate(sys%A)
!    if(allocated(sys%thermo)) deallocate(sys%thermo)

       if (ncs>0) deallocate(sys%CS) 
       deallocate(sys%sp_order)
       deallocate(sys%el_order) 
       deallocate(sys%Ein)
       deallocate(sys%E)
       if (ng>0)  deallocate(sys%Bg)
       deallocate(sys%B)
       deallocate(sys%BR)
      deallocate(sys%A)
      deallocate(sys%thermo)
       deallocate(sys)

   return
   end subroutine ceq_sys_rm2  !--------------------------------------------------------------



   subroutine ceq_param_set( sys, diag, lu, T_low, T_high, frac_zm, T_tol, ds_inc, &
        ds_dec, res_tol, ires_fac, logy_lim, srat_lim, err_huge, dec_min, &
		eps_el, eps_sp, pert_tol, pert_skip )

! Reset values of parameters in sys: apart from sys, all arguments are optional.
! If sys%diag=5, the values of the parameters are output.

   type (sys_type) :: sys
   integer,    intent(in), optional :: diag, lu, pert_skip
   real(k_wp), intent(in), optional :: T_low, T_high, frac_zm, T_tol, ds_inc, ds_dec, &
               res_tol, ires_fac, logy_lim, srat_lim, err_huge, dec_min, &
			   eps_el, eps_sp, pert_tol
   integer :: luu

   if( present(diag) )     sys%diag     = diag
   if( present(lu) )       sys%lu       = lu
   if( present(T_low) )    sys%T_low    = T_low
   if( present(T_high) )   sys%T_high   = T_high
   if( present(frac_zm) )  sys%frac_zm  = frac_zm
   if( present(T_tol) )    sys%T_tol    = T_tol
   if( present(ds_inc) )   sys%ds_inc   = ds_inc
   if( present(ds_dec) )   sys%ds_dec   = ds_dec
   if( present(res_tol) )  sys%res_tol  = res_tol
   if( present(ires_fac) ) sys%ires_fac = ires_fac
   if( present(logy_lim) ) sys%logy_lim = logy_lim
   if( present(srat_lim) ) sys%srat_lim = srat_lim
   if( present(err_huge) ) sys%err_huge = err_huge
   if( present(dec_min) )  sys%dec_min  = dec_min
   if( present(eps_el) )   sys%eps_el   = eps_el
   if( present(eps_sp) )   sys%eps_sp   = eps_sp
   if( present(pert_tol) ) sys%pert_tol = pert_tol
   if( present(pert_skip)) sys%pert_skip=pert_skip

   if( sys%diag < 5 ) return

   luu = sys%lu
   write(luu,*)' '
   write(luu,*)' Parameters output from ceq_param_set '
   write(luu,*)' '

   write(luu, '(a,i4)' ) 'diag      = ', sys%diag
   write(luu, '(a,i4)' ) 'lu        = ', sys%lu  
   write(luu, '(a,i4)' ) 'pert_skip = ', sys%pert_skip 

   write(luu, '(a,1p,e13.4)' ) 'T_low    = ', sys%T_low
   write(luu, '(a,1p,e13.4)' ) 'T_high   = ', sys%T_high
   write(luu, '(a,1p,e13.4)' ) 'frac_zm  = ', sys%frac_zm
   write(luu, '(a,1p,e13.4)' ) 'T_tol    = ', sys%T_tol
   write(luu, '(a,1p,e13.4)' ) 'ds_inc   = ', sys%ds_inc 
   write(luu, '(a,1p,e13.4)' ) 'ds_dec   = ', sys%ds_dec 
   write(luu, '(a,1p,e13.4)' ) 'res_tol  = ', sys%res_tol
   write(luu, '(a,1p,e13.4)' ) 'ires_fac = ', sys%ires_fac
   write(luu, '(a,1p,e13.4)' ) 'logy_lim = ', sys%logy_lim
   write(luu, '(a,1p,e13.4)' ) 'srat_lim = ', sys%srat_lim
   write(luu, '(a,1p,e13.4)' ) 'err_huge = ', sys%err_huge
   write(luu, '(a,1p,e13.4)' ) 'dec_min  = ', sys%dec_min
   write(luu, '(a,1p,e13.4)' ) 'eps_el   = ', sys%eps_el
   write(luu, '(a,1p,e13.4)' ) 'eps_sp   = ', sys%eps_sp
   write(luu, '(a,1p,e13.4)' ) 'pert_tol = ', sys%pert_tol

   return
   end subroutine ceq_param_set  !---------------------------------------------------------------

   subroutine ceq_param_def( sys )

!  Reset parameters in sys to their default settings.

   type (sys_type), pointer :: sys

   	sys%initialized = .false. ! = .true. when sys has been initialized

    sys%diag     = 1       ! >0 for diagnostics
    sys%lu       = 0       ! logical unit number for diagnostics
	sys%T_low    = 250.d0  ! lowest allowed temperature
	sys%T_high   = 5000.d0 ! highest allowed temperature
	sys%frac_zm  = 1.d-1   ! fraction of zm used in initial guess
	sys%T_tol    = 1.d-6   ! convergence tolerance on log(T)
    sys%ds_inc   = 1.4d0   ! factor by which ds is increased after success (ceq_fixed_T)
    sys%ds_dec   = 0.25d0  ! factor by which ds is decreased after failure (ceq_fixed_T)
	sys%ds_min   = 1.d-15  ! smallest allowed time step
    sys%res_tol  = 1.d-9   ! convergence tolerance for residual (ceq_fixed_T)
    sys%ires_fac = 1.d2    ! factor by which the irreducible residual can exceed res_tol
	sys%logy_lim = 120.d0  ! upper limit on log(y)  (to prevent overflow in ceq_y)
	sys%srat_lim = 1.d-9   ! lower limit on singular-value ratio (ceq_Sinv)
	sys%err_huge = 1.d6    ! upper limit on error (ceq_newt)
	sys%dec_min  = 0.5d0   ! minimum acceptable decrease in residual (ceq_newt)
	sys%eps_el   = 1.d-9   ! relative lower bound on element moles (ceq_perturb)
	sys%eps_sp   = 1.d-9   ! relative lower bound on species moles (ceq_perturb)
	sys%pert_tol = 1.d-4   ! largest allowed normalized perturbation
	sys%pert_skip= 0       ! set =1 to skip perturbing undetermined species


   return
   end subroutine ceq_param_def

end module ceq_system
