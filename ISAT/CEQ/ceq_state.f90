!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module ceq_state_m

!  CEQ Release 1.0
!  Copyright, 2003, Stephen B. Pope;  All rights reserved.

   use ceq_solve
   implicit none

contains

subroutine ceq_state( sys, c, N, p_atm, p_Pa, p_cgs, T, HoR, N_h, T_h, &
                      N_g, T_g, N_eq, T_eq, HoR_eq, stats, info )

!  Subroutine ceq_state determines the constrained equilibrium state of 
!  an ideal gas mixture consisting of ns species, either at fixed pressure and
!  temperature (p,T), or at fixed pressure and enthalpy (p,H).

!  The nc equality constraints are written:  B' * N = c, where B is the
!  ns x nc basic constraint matrix, N is the ns-vector of species moles,
!  and c is the nc-vector of constraint values.

!  Usage:  
!     The subroutine ceq_sys_init must be called prior to this routine to
!     initialize the data structure  sys. 

!     The calling routine must contain the statements:
!        use ceq_state
!        type (sys_type), pointer :: SYS  
!     where SYS is the actual name corresponding to the dummy argument sys.

!     All arguments, except for the first (sys), are "optional".  Hence they
!     should be referred to by name, e.g.,
!        call ceq_state( SYS, N=N, p_atm=1.d0, T=1500.d0, N_eq=NCEQ, info=INFO )
!     Nevertheless, some "optional" arguments must be provided, including info.

!     All real variables are real(kind(1.d0)) (i.e., double precision).

!     On return from ceq_state, info should be checked. (info >= 0 for success)

!  Input:
!    sys  - data structure of type sys_type created by subroutine ceq_sys_init.
!           Required first argument

!    (Either c or N must be specified.)
!    c    - values of the nc basic constraints (real(nc))
!    N    - moles of species used to calculate c as : c = B' * N (real(ns))

!    (Either p_atm or p_Pa or p_cgs must be specified.)
!    p_atm - pressure in standard atmospheres (real)
!    p_Pa  - pressure in Pascals (SI units) (real)
!    p_cgs - pressure in cgs units (dynes/cm^2) (real)

!    (For a fixed (p,T) equilibrium calculation, specify T: do not specify HoR, N_h or T_h.)
!    T - temperature (K) for fixed-temperature problem

!    (For a fixed (p,H) equilibrium calculation, specify either HoR or N_h and T_h: 
!     do not specify T.)
!    HoR the fixed value of H/R [moles K], where H=enthalpy, R=universal gas constant.
!    N_h species moles used to calculate H  (real(ns))
!    T_h temperature used to calculate H as:  H/R = sum( N_h h(T_h)/R ), where
!        h(T_h)/R [which has dimensions K] is the molar specific species enthalpy.

!    (Initial guesses are not needed, and should not be specified unless they
!     are good guesses.)
!    N_g   - initial guess for species moles
!    T_g   - initial guess for temperature (for fixed (p,H) only)

!  Output:
!    N_eq     - moles of species in equilibrium mixture (real(ns))
!    T_eq     - temperature of the equilibrium mixture (real)
!    HoR_eq   - value of H/R [moles K] for the equilibrium mixture (real)
!    stats    - numerical information about the solution - details below (real(20))
!    info     - information about the solution:
!             > 0, info = the number of Newton iterations performed
!             < 0, an error occurred - see details below

!  Error conditions:  
!    info = -1  sys is not associated (ceq_sys_init must be called first)
!    info = -2  sys has not been initialized (ceq_sys_init must be called first)
!    info = -3  the pressure has been specified more than once
!    info = -4  the pressure has not been specified
!    info = -5  the pressure is non-positive
!    info = -6  the specified temperature is outside the range (T_low, T_high)
!    info = -7  neither T nor enthalpy specified
!    info = -8  the guess for temperature T_g is outside the range (T_low, T_high)
!    info = -9  neither c nor N has not been specified
!    info = -10 the constraints (and determined species) are zero
!    info = -11 the constraints are not realizable
!    info = -12 (ceq_maxmin failed)
!    info = -13 (ceq_ming failed)
!    info = -14 (ceq_h2T failed)
!    info = -15 (ceq_lamg_init failed)
!    info = -16 (ceq_fixed_T failed)
!    info = -17 (ceq_fixed_h failed)
!    info = -18 T_eq is greater than T_high
!    info = -19 T_eq is less than T_low

!  Values of info between -12 and -17 are internal failures (which could be caused by erroneous
!  input).  Other values are caused by erroneous input.  A temperature outside the range
!  (T_low, T_high) may be caused by erroneous input (e.g., of T or HoR).  If the correct
!  temperature is outside this range, the values of T_low and T_high can be changed prior to
!  calling ceq_state by calling ceq_param_set, e.g., call ceq_param_set( sys, T_low = 200.d0 ).
!  To diagnose an error condition, the case can be repeated with diagnostics turned on,
!  by: call ceq_param_set( sys, diag=5 ).

! stats:
!    stats(1) = temp_its   - number of temperature iterations
!    stats(2) = time_steps - number of psuedo-time steps
!    stats(3) = newt_calls - number of Newton calls (iteration not performed if initial residual small)
!    stats(4) = newt_its   - number of Newton iterations
!    stats(5) = nyeval     - number of temperature iterations
!    stats(6) = err        - last residual error
!    stats(7) = npert      - number of quantities perturbed
!    stats(8) = max_pert   - maximum value of perturbation
!    stats(9) = zmm        - minimum maxmin composition

   type (sys_type),   pointer :: sys
   real(k_wp), intent(in),  optional :: c(sys%nc), N(sys%ns), p_atm, p_Pa, p_cgs, T, &
                                        HoR, N_h(sys%ns), T_h, N_g(sys%ns), T_g
   real(k_wp), intent(out), optional :: N_eq(sys%ns), T_eq, HoR_eq, stats(20)
   integer,    intent(out), optional :: info

   type (state_type), pointer :: state
   integer    :: nb, nc, ns, nsd, nsu, nrc, npert, iret, i, lu
   real(k_wp) :: p, T0, errc, max_pert, Nd, Nu_atoms, zumin, zumax, cb(sys%nc), cmod(sys%nb), &
                 zd(sys%nsd), err, cr_norm, cb_norm, res, z_low, ceq_norm, res_tol, err2
   real(k_wp), dimension(sys%ns)  :: z0, z1, h, zeq
   real(k_wp), dimension(sys%nsu) :: zu, zu0, zm, zupper, atoms, zg, gu, xu, rhs, gu0, zu00
   real(k_wp), dimension(sys%nrc) :: cr, lam0
   real(k_wp), parameter :: p0_Pa=1.01325d5 ! standard atmosphere [Pa]
   logical :: fixed_T, fail, diag

   if( .not.present(info) ) then
      write(0,*)'ceq_state: argument INFO must be present; stopping'
	  stop
   endif

   info  = 0  ! anticipate success
   if( present(stats) ) stats = 0.d0
   err = 0.d0

!  check that sys has been initialized  =====================================
   if( .not.associated(sys) ) then
	  info = -1
	  return
   endif

   if( .not.sys%initialized ) then
	  info = -2
	  return
   endif

!  obtain indexes and allocate state

   nb  = sys%nb
   nc  = sys%nc
   nrc = sys%nrc
   ns  = sys%ns
   nsd = sys%nsd
   nsu = sys%nsu

   call ceq_state_init( sys, state )

!  determine if errors are to be reported
   if( sys%diag >= 1 ) then
      diag = .true.
   else
      diag = .false.
   endif
   lu = sys%lu

!  determine pressure  p  in standard atmospheres  !=================
   if( present(p_atm) ) then
      p = p_atm
	  if( present(p_Pa) .or. present(p_cgs) ) then
	     if( diag ) write(lu,*) 'ceq_state: multiple input pressures'
		 info = -3
		 return
	  endif

   elseif( present(p_Pa) ) then
      p = p_Pa / p0_Pa 
	  if( present(p_atm) .or. present(p_cgs) ) then
	     if( diag ) write(lu,*) 'ceq_state: multiple input pressures'
		 info = -3
		 return
	  endif

   elseif( present(p_cgs) ) then
      p = p_cgs * 0.1d0 / p0_Pa  
	  if( present(p_atm) .or. present(p_Pa) ) then
	     if( diag ) write(lu,*) 'ceq_state: multiple input pressures'
		 info = -3
		 return
	  endif
   else
   	  if( diag ) write(lu,*) 'ceq_state: no pressure specified'
	  info = -4
	  return
   endif

   if( p <= 0.d0 ) then
   	  if( diag ) write(lu,'(a,1p,9e13.4)') 'ceq_state: pressure must be strictly positive: p = ', p
	  info = -5
	  return
   endif

   state%p = p

!  determine T or HoR  =============================================

   fixed_T = .false.
   if( present(T) ) then
      state%T = T
	  fixed_T = .true.

	  if( T < sys%T_low  .or.  T > sys%T_high ) then
	     if( diag ) write(lu,'(a,1p,9e13.4)') &
		      'ceq_state: temperature out of range: T, T_low, T_high = ', T, sys%T_low, sys%T_high
	     info = -6
		 return
      endif

   elseif( present(HoR) ) then
	  state%h = HoR

   elseif( present(N_h)  .and.  present(T_h) ) then
      call ceq_reorder( sys%ns, 1, N_h, sys%sp_order, z0 )
	  call ceq_h( sys%ns, T_h, sys%thermo, h )
	  state%h = sum( z0 * h ) * T_h  ! HoR

   else
	  if( diag ) write(lu,*) 'ceq_state: neither T nor HoR nor (N_h and T_h) specified'
	  info = -7
	  return
   endif

!  determine T0 for initial guess, and evaluate gu0 = gu(T0) ===============

   if( fixed_T ) then
      T0 = state%T  ! specified fixed T

   elseif( present(T_g) ) then
      T0 = T_g  ! guess provided

	  if( T_g < sys%T_low  .or.  T_g > sys%T_high ) then
	     if( diag ) write(lu,'(a,1p,9e13.4)') 'ceq_state: T_g out of range: T_g, T_low, T_high = ', &
                                 T_g, sys%T_low, sys%T_high
		 info = -8
	     return
      endif

   else
      T0 = sqrt( sys%T_low * sys%T_high )
	  T0 = max( T0, 0.1d0 * sys%T_high )
   endif

   call ceq_g( nsu, T0, p, sys%thermo(nsd+1:ns,:), gu )

!  form the basic constraint vector  =========================

   if( present(c) ) then
      cb(1:nc) = c

   elseif( present(N) ) then
      cb(1:nc) = matmul( N, sys%B )

   else
   	  if( diag ) write(lu,*) 'ceq_state: neither c nor N specified'
	  info = -9
	  return
   endif

!  form modified and reduced constraints  !===================

   cmod(1:nb) = matmul( sys%A(1:nb,1:nc) ,cb )   ! form modified constraint vector
   zd(1:nsd)  = cmod(1:nsd)                ! determined species

!  treat the special case of no undetermined species

   if( nsu == 0 ) then
      zeq = zd
	  if( .not.fixed_T ) then
	     call ceq_h2T(ns,zeq,state%h,sys%T_low,sys%T_high,sys%thermo,state%T,iret)
		 if( iret < 0 ) then
		    info = -17+iret
			return
		 endif
	  endif
	  
	  go to 500
   endif
    
   cr(1:nrc)  = cmod(nsd+1:nb)             ! reduced constraint vector
   cr_norm    = ceq_norm( nrc, cr )        ! |cr|

   if( cr_norm <= 0.d0 ) then
   ! SBP added 4/9/2009
      if( cr_norm == 0.d0  .and.  nsd > 0  .and.  sum( zd(1:nsd) ) > 0.d0 ) then
      !  only determined species
      zeq = 0.d0
      zeq(1:nsd) = zd(1:nsd)
      if( .not.fixed_T ) then
	     call ceq_h2T(ns,zeq,state%h,sys%T_low,sys%T_high,sys%thermo,state%T,iret)
		 if( iret < 0 ) then
		    info = -17+iret
			return
		 endif
	   endif
	  
	   go to 500
   endif
   ! SBP end of added
   
   	  if( diag ) write(lu,*) 'ceq_state: constraints are zero'
	  info = -10  ! SBP bug fix 4/8/2009
	  return
   endif

!  use initial guess N_g if provided  !====================

   if( present(N_g) ) then
      call ceq_reorder( sys%ns, 1, N_g, sys%sp_order, z0 )
	  zu0(1:nsu) = z0(nsd+1:ns)  ! guessed undetermined species

	  cb(1:nrc) = matmul( zu0(1:nsu), sys%BR ) ! reduced c.v. based on N_g
	  cb_norm   = ceq_norm( nrc, cb )

	  if( cb_norm == 0.d0 ) then
	     if( sys%diag >=2 ) write(lu,*)'ceq_state: N_g yields zero constraint'
	     go to 50
	  endif  

	  res       = ceq_norm( nrc, (cb(1:nrc)/cb_norm - cr/cr_norm) )  !  residual 

	  if( sys%diag >= 2 ) write(lu,'(a,1p,2e13.4)')'ceq_state: initial guess residual = ', res

	  if( res > sys%res_tol ) then  !  adjust initial guess zu0, store in zm

	     z_low = sum( cr(1:sys%neu) ) * 1.d-15

	     call ceq_min_pert( nsu, nrc, sys%BR, cr, zu0, z_low, zm, info )
		 if( info /= 0 ) then  !  ceq_min_pert failed
            if( sys%diag >= 2 ) write(lu,*)'ceq_state: ceq_min_pert failed, info = ', info
			go to 50
         endif

		 zu0       = zm

		 if( sys%diag >= 2 ) then
	        cb(1:nrc) = matmul( zu0(1:nsu), sys%BR ) ! reduced c.v. based on N_g
	        cb_norm   = ceq_norm( nrc, cb )
	        res       = ceq_norm( nrc, (cb(1:nrc)/cb_norm - cr/cr_norm) )
		    write(0,'(a,1p,9e13.4)')  &
		        'ceq_state: after ceq_min_pert: res, min(zu0) = ', res, minval(zu0)
		 endif
      endif

!     accept initial guess (zu0); skip max-min and min_g

	  state%zd = zd
      state%cr = cr

	  go to 100
   endif

50    continue  !  proceed without using initial guess N_g

!  perturb if necessary  =====================================
   
   call ceq_perturb( ns, nsd, nsu, sys%ne, sys%ned, sys%neu, nrc, &
                     zd, cr, sys%BR, sys%E, sys%diag, sys%lu,  &
					 sys%eps_el, sys%eps_sp, sys%pert_tol, sys%pert_skip, &
                     state%zd, zm, zupper, state%cr, npert, max_pert, iret )

   if( iret == -1 ) then 
   	  if( diag ) write(lu,*) 'ceq_state: non-realizable constraint'
	  info = -11
	  return

   elseif( iret == -2 ) then 
   	  if( diag ) write(lu,*) 'ceq_state: ceq_max-min failed'
	  info = -12
	  return
   endif

   if( present(stats) ) then
      stats(7) = npert
	  stats(8) = max_pert
	  stats(9) = minval(zm)
   endif


!  determine min_g composition based on T0  !===========

   call ceq_ming( nsu, nrc, sys%BR, state%cr, gu, zg, iret )
    
   if( iret < 0 ) then
   	  if( diag ) write(lu,*) 'ceq_state: ceq_ming failed'
	  info = -13
	  return
   endif

!  form initial guess zu0   !============================

   zu0 = zg + sys%frac_zm * (zm-zg)

100 continue  !  skip to here if zu0 based on N_g accepted

!  determine normalization-condition vector Q !=========
!     q = Q' * Xu -1 = 0
   Nd       = sum( state%zd )               ! total moles of determined species
   atoms    = sum( sys%E(nsd+1:ns,:), 2 )   ! atoms in undetermined species
   Nu_atoms = sum( state%cr(1:sys%neu) )    ! total moles of atoms in undetermined species
   state%Q  = 1.d0 + atoms * Nd / Nu_atoms  ! consistency condition vector

!  re-estimate T0 and re-evaluate gu if required !======

   if( .not.fixed_T  .and.  .not.present(T_g) ) then
	   z1(1:nsd)    = state%zd
	   z1(nsd+1:ns) = zu0
	   call ceq_h2T( ns, z1, state%h, sys%T_low, sys%T_high, sys%thermo, T0, iret )

	   if( iret == -3 ) then
   	       if( diag ) write(lu,*) 'ceq_state: ceq_h2T failed'
	       info = -14
		   return
	   endif

	   call ceq_g( nsu, T0, p, sys%thermo(nsd+1:ns,:), gu )  ! set gu based on T0
   endif

   state%gu = gu

!  set lam0 and gu0 to satisfy constraints  ==============

   call ceq_lamg_init( sys, state, nsu, nrc, zu0, gu, lam0, gu0, iret )

   if( iret < 0 ) then
      if( diag ) write(lu,*) 'ceq_state: ceq_lamg_init failed, info = ', info
	  info = -15
	  return
   endif

!===============  perform equilibrium calculation  ==================
   
   if( fixed_T ) then         ! fixed (p,T)

      call ceq_fixed_T( sys, state, nsu, nrc, gu, gu0, lam0, 0.d0, sys%res_tol, &
	                    state%lam, state%zu, err, fail )

	  if( fail ) then
	     info = -16
		 if( diag ) write(lu,*) 'ceq_state: ceq_fixed_T failed'
		 return
	  endif

   else                       ! fixed (p,h)

      call ceq_fixed_h( sys, state, nsu, nrc, T0, state%h, gu0, lam0, sys%res_tol, &
                           state%lam, state%zu, state%T, err, iret )

	  if( iret == -1  .or.  iret == -2 ) then
	     info = -17
		 if( diag ) write(lu,*) 'ceq_state: ceq_fixed_h failed'
		 return
      elseif( iret == -3 ) then
	     info = -18
         if( diag ) write(lu,*) 'ceq_state: T_eq is greater than T_high'
		 return
      elseif( iret == -4 ) then
	     info = -19
         if( diag ) write(lu,*) 'ceq_state: T_eq is less than T_low'
		 return
	  endif

   endif

!  determine required output  =======================================

   zeq = (/ state%zd , state%zu /)  !  equilibrium moles (re-ordered)

500	   continue  !  jump to here if there are no undetermined species

   if( present(T_eq) ) T_eq = state%T

   if( present(HoR_eq) ) then
      if( fixed_T ) then
         call ceq_h( ns, state%T, sys%thermo, h )
         state%h = sum( zeq * h ) * state%T
	  endif

	  HoR_eq = state%h
   endif

   if( present(N_eq) ) then
      do i = 1, ns          ! re-order species
         N_eq( sys%sp_order(i) ) = zeq(i)
      end do
   endif

   if( nsu > 0 ) then
      if( present(stats) ) then
         stats(1) = state%temp_its
         stats(2) = state%time_steps
         stats(3) = state%newt_calls
         stats(4) = state%newt_its
         stats(5) = state%nyeval
	     stats(6) = err
      endif

      info = state%nyeval
   endif
   
   call ceq_state_rm( state )

   return


end subroutine ceq_state  !-----------------------------------------------------


subroutine ceq_state_init( sys, state )
!  allocate state
   
   type (sys_type),   pointer :: sys
   type (state_type), pointer :: state

   allocate( state )
   allocate( state%zd(sys%nsd) )
   allocate( state%zu(sys%nsu) )
   allocate( state%lam(sys%nrc) )
   allocate( state%cr(sys%nrc) )
   allocate( state%Q(sys%nsu) )
   allocate( state%gu(sys%nsu) )

   
   state%temp_its   = 0 ! number of temperature iterations performed
   state%time_steps = 0 ! number of pseudo time steps
   state%newt_calls = 0 ! number of Newton solves attempted
   state%newt_its   = 0 ! number of Newton iterations performed
   state%nyeval     = 0 ! number of y evaluations

   return

end subroutine ceq_state_init  !------------------------------------------------

subroutine ceq_state_rm( state )
!  deallocate state
   
   type (state_type), pointer :: state

   deallocate( state%zd, state%zu, state%lam, state%cr, state%Q, state%gu )
   deallocate( state )

   return
end subroutine ceq_state_rm  !------------------------------------------------

subroutine ceq_lamg_init( sys, state, nsu, nrc, zu0, gu, lam0, gu0, info )

!  given gu and initial guess zu0, determine lam0 and gu0 to satisfy constraints.

   type (sys_type),   pointer :: sys
   type (state_type), pointer :: state

   integer, intent(in)     :: nsu, nrc
   real(k_wp), intent(in)  :: zu0(nsu), gu(nsu)
   real(k_wp), intent(out) :: lam0(nrc), gu0(nsu)
   integer,    intent(out) :: info

   integer    :: i
   real(k_wp) :: xu(nsu), rhs(nsu), res, ceq_norm, BRw(nsu,nrc), wt(nsu)
   real(k_wp), parameter :: pow = 0.0d0  !  weight proportinal to X(i)**pow

!  set lam0 to the least=squares solution to:  log(xu) = -gu + BR * lam  !====

   xu = zu0 / ( sum(state%zd) + sum(zu0) )  ! mole fractions

   do i = 1, nsu
      rhs(i) = gu(i) + log( xu(i) )
	  wt(i)  = xu(i)**pow
   end do

   do i = 1, nrc
      BRw(:,i) = sys%BR(1:nsu,i) * wt(1:nsu)
   end do

   call ceq_lss( nsu, nrc, BRw, rhs, lam0, info )

!  set gu0 such that  xu = exp( -gu0 + BR * lam0 ) 

   gu0 = matmul( sys%BR, lam0 )
   do i = 1, nsu
      gu0(i) = gu0(i) - log( xu(i) )
   end do

   return

end subroutine ceq_lamg_init  !-------------------------------------------------

end module ceq_state_m
