!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module ceq_solve

use ceq_types
implicit none

contains

subroutine ceq_fixed_T( sys, state, nsu, nrc, gu, gu0, lam0, res0, res_tol, lam, zu, err, fail )

! Use the converging Gibbs function method to determine the Lagrange
! multipliers lam corresponding to constrained equilibrium at fixed p and T.

! Input:
!   gu     - Gibbs function (at s=1)
!   gu0    - initial (s=0) value of Gibbs function
!   lam0   - initial (s=0) value of Lagrange multipliers
!   res0   - initial (s=0) residual, based on gu0 and lam0
!   res_tol- tolerance for residual

! Output:
!   lam    - final (s=1) value of Lagrange multipliers
!   zu     - moles of undetermined species
!   err    - the residual error (may be larger than res_tol)
!   fail   - true if failed to obtain solution

   type (sys_type),   pointer :: sys
   type (state_type), pointer :: state
   integer, intent(in)        :: nsu, nrc
   real(k_wp), intent(in)     :: gu(nsu), gu0(nsu), lam0(nrc), res0, res_tol
   real(k_wp), intent(out)    :: lam(nrc), zu(nsu), err
   logical, intent(out)       :: fail

   integer    :: ns, nsd
   real(k_wp) :: rtol, dsmin, s, ds, err_i, send, gus(nsu), dlamds(nrc), &
                 lam1(nrc), lam2(nrc), err_tol, v(nrc), y(nsu), r, q, dguds(nsu), &
				 newt_tol, ceq_norm

   ns    = sys%ns
   nsd   = sys%nsd
   dsmin = sys%ds_min   ! smallest allowed time step
   fail  = .false.      ! anticipate success

!  if required, check initial residuals
   if( sys%diag >=3 ) then
	 call ceq_residual( sys, state, nsu, nrc, gu0, lam0, q, r, fail )

	 if( fail ) then
	    write(sys%lu,*)'ceq_fixed_T: ceq_residual failed'
		return
     endif
	    
	 write(sys%lu,'(a,1pe13.4)') 'residual q = ', q
	 write(sys%lu,'(a,1pe13.4)') 'residual r = ', r
   endif

!  initialize prior to time steps
   s     = 0.d0     ! time at start of step
   send  = 1.d0     ! time at end of step
   lam   = lam0     ! value of lam(s)
   err_i = 0.d0     ! irreducible error
   dguds = gu - gu0 ! d(gu)/ds
   err_tol = max( res0 + res_tol , 1.05 * res0 ) 
   err   = 0.d0     ! error on last accepted step

   do while( s < 1.d0 )   !  integrate from s=0 to s=1
      ds = send - s
	  state%time_steps = state%time_steps + 1
      if( sys%diag >=3 )  &
         write(sys%lu,'(4x,a,1p,2e10.2,2i5,1p,e10.2)') &
		      'time: s, ds, steps, newts, err = ', s, ds, state%time_steps, state%newt_its, err

      if( ds < dsmin ) then
	     fail = .true.
         if( sys%diag >= 1 ) write(sys%lu,*) 'ceq_fixed_T: failure, time step collapsed'
		 return
	  endif
    
!  take explicit Euler step (from s to send)
      gus = gu - (1-s) * dguds
      call ceq_rate( sys, state, nsu, nrc, lam, gus, dguds, state%Q, dlamds, fail )

      if( fail ) then  ! overflow in evaluating y
	     send = s + ds * sys%ds_dec
		 cycle         ! repeat step with smaller ds 
	  endif
    
      lam1 = lam + ds * dlamds  ! provisionally accept time step
     
!  perform Newton iteration to convergence
	  gus     = gu - (1-send) * dguds  ! gu(send)

      newt_tol = res_tol
	  call ceq_newt( sys, state, nsu, nrc, gus, lam1, newt_tol, lam2, y, err, err_i, fail )

      if( fail  .or.  err > err_tol ) then  !  Newton did not converge to acceptable tolerance
         send = s + ds * sys%ds_dec 
         cycle         !  decrease ds and repeat time step
      endif
	
!  acceptable convergence: set acceptable tolerance for next step	   
      err_tol = max( err + res_tol , 1.05 * err ) 

      lam  = lam2  ! Newton success,  accept step
      s    = send  ! increment time
      send = min( send + ds * sys%ds_inc,  1.d0 )  ! increase time step

   end do

!  converged solution obtained at s=1

   if( sys%diag >=3 ) write(sys%lu,'(4x,a,1p,2e10.2,2i5,1p,e10.2)') &
		      'time: s, ds, steps, newts, err = ', s, ds, state%time_steps, state%newt_its, err

   zu  = y * y
   v   = matmul( zu, sys%BR )
   zu  = zu * ceq_norm( nrc, state%cr ) / ceq_norm( nrc, v )

   return

end subroutine ceq_fixed_T  !----------------------------------------------------

subroutine ceq_fixed_h( sys, state, nsu, nrc, T0, h0, gu0, lam0, res_tol, &
                           lam, zu, T, err, iret )

!  Determine the CE composition at fixed h and p 
!   S.B. Pope 7/5/03

! Input
!   T0      - initial guess for T [K]. 
!   h0      - the value of H/R [K]
!   gu0     - initial Gibbs function
!   lam0    - initial value of lam 
!   res_tol - tolerance for residual

! Output
!   lam     - Lagrange vector
!   zu      - moles of undetermined species
!   T       - equilibrium temperature
!   err     - residual  (may be larger than res_tol)
!   iret    = 0 for success
!           = -1  ceq_fixed_T failed
!           = -2  ceq_cp_eff failed
!           = -3  T > T_high
!           = -4  T < T_low

   type (sys_type),   pointer :: sys
   type (state_type), pointer :: state
   integer, intent(in)     :: nsu, nrc
   real(k_wp), intent(in)  :: T0, h0, gu0(nsu), lam0(nrc), res_tol
   real(k_wp), intent(out) :: lam(nrc), zu(nsu), T, err
   integer,    intent(out) :: iret

   integer    :: ns, nsd, ifop
   real(k_wp) :: T_low, T_high, dlogT_tol, Tlo, Thi, dT, gu(nsu), &
                 lam1(nrc), cp_eff, h(sys%ns), Tn, hlo, hhi, h1, guold(nsu), &
				 z(sys%ns), res0
   logical    :: fail

   ifop      = sys%diag
   T_low     = sys%T_low
   T_high    = sys%T_high
   dlogT_tol = sys%T_tol
   iret      = 0

   ns  = sys%ns
   nsd = sys%nsd

   Tlo = -1.d30     ! lowest temperature at which h has been evaluated
   Thi =  1.d30     ! highest temperature at which h has been evaluated

   T    = T0         ! set starting values
   dT   = T0
   gu   = gu0
   lam  = lam0
   res0 = 0.d0

   do while( abs(dT) >= T * dlogT_tol ) !  iterate over temperature 
      state%temp_its = state%temp_its + 1

!  determine equilibrium composition at current T
      guold      = gu
	  call ceq_g( nsu, T, state%p, sys%thermo(nsd+1:ns,:), gu)

      call ceq_fixed_T( sys, state, nsu, nrc, gu, guold, lam, res0, res_tol, lam1, zu, err, fail )

	  if( fail ) then
	     if( sys%diag >=1 ) write(sys%lu,*) 'ceq_fixed_h: ceq_fixed_T failed'
		 iret = -1
		 return
	  endif

	  lam  = lam1
	  res0 = err

      call ceq_cp_eff( sys, state, ns, nsd, nsu, nrc, zu, lam, T, cp_eff, fail )

	  if( fail ) then
	     if( sys%diag >=1 ) write(sys%lu,*) 'ceq_fixed_h: ceq_cp_eff failed'
		 iret = -2
		 return
	  endif	
    
! predict dT based on Newton's method 
      z = (/ state%zd , zu /)
	  call ceq_h( ns, T, sys%thermo, h )   ! obtain species h/(RT)
	  h1 = T * sum( z * h ) ! mixture h/R 
      dT = ( h0 - h1 ) / ( cp_eff * sum(z) ) ! predicted dT based on constant cp_eff
    
! check that T is within limits
      if( T == T_high .and. dT > 0 ) then
         if( sys%diag >= 1 ) write( sys%lu, '(a,1p,2e13.4)' ) &
		     'ceq_fixed_h: T > T_high ', T, T_high
         iret=-3    ! fatal error: T > T_high
         return 
      endif
    
      if( T == T_low .and. dT < 0 ) then 
         if( sys%diag >= 1 ) write( sys%lu, '(a,1p,2e13.4)' ) &
		     'ceq_fixed_h: T < T_low ', T, T_low
         iret=-4    ! fatal error: T < T_low
         return 
      endif
    
!    ensure that Tn is within limits  
      Tn = T+dT  
      Tn = min( Tn, T_high )
      Tn = max( Tn, T_low )
    
      if( dT > 0 ) then! use linear interpolation instead if Tn is closer to known bound
         Tlo   = T
         hlo   = h1
         if( Tn >.5 * (Tlo+Thi) ) then
            Tn = Tlo + (Thi-Tlo) * (h0-hlo) / (hhi-hlo)
         endif
      else
         Thi   = T
         hhi   = h1
        if( Tn < .5 * (Tlo+Thi) ) then
           Tn = Tlo + (Thi-Tlo) * (h0-hlo) / (hhi-hlo)
        endif
      endif

	  dT = Tn - T   ! temperature increase for next iteration

 	  if( sys%diag >= 2 ) write(sys%lu,'(a,1p,3e11.4,3i5)') &
	     'Temp: T, Tn, dT, temps, steps, newts = ', T, Tn, dT, state%temp_its, state%time_steps, state%newt_its

      T  = Tn       ! temperature for next iteration  
   end do

   return

end subroutine ceq_fixed_h  !----------------------------------------------------


subroutine ceq_newt( sys, state, nsu, nrc, gu, lam0, err_tol, lam, y, err, err_i, fail )

! Solve for lam by Newton iteration at time s.

! Input:
!   gu      - Gibbs function
!   lam0    - initial guess for lam
!   err_tol - error tolerance

! Output:
!   lam     - solution
!   y
!   err     - error in solution (may be greater than err_tol)
!   err_i   - the irreducible error
!   fail    = .true.  if failed to obtain any solution

   type (sys_type),   pointer :: sys
   type (state_type), pointer :: state

   integer, intent(in)        :: nsu, nrc
   real(k_wp), intent(in)     :: gu(nsu), lam0(nrc), err_tol

   logical, intent(out)    :: fail
   real(k_wp), intent(out) :: lam(nrc), y(nsu), err, err_i

   integer    :: lwork, info, n_small, i
   real(k_wp) :: err_last, HH(nsu,nrc), v(nrc), w(nrc), lam_v(nrc), lam_w(nrc),  &
                 P(nrc), VTw(nrc), q, r, work(20*(nsu+nrc)), Sig(nrc), U(nsu,nrc), &
				 VT(nrc,nrc), Sinv(nrc), RR(nsu), dlam(nrc), vnorm, cnorm, v_star, w_star, &
				 d_al, d_alp, c(nrc), err_best, y_best(nsu), lam_best(nrc), ceq_norm, denom, &
				 HHcopy(nsu,nrc)
				 

   lam  = lam0
   fail = .false.  !  anticipate success
   state%newt_calls = state%newt_calls + 1

   err      = 2.d0 * sys%err_huge  !  current error
   err_last = err                  !  error on last iteration
   err_best = err                  !  smallest error so far
   err_i    = 0.d0
   c        = state%cr
   cnorm    = ceq_norm( nrc, c )
   Sig      = 1.d0

   do !  iterate until exit due to (a) ceq_y failure (b) err < err_tol (c) no longer converging

      call ceq_y( sys, state, nsu, nrc, gu, lam, y, fail )
	  if( fail )  then  ! (a)  y too large - quit
	     err = 2.d0 * sys%err_huge 
		 if( sys%diag >=5 ) write(sys%lu,*) 'ceq_newt: ceq_y failed'
	     exit  !  previous solution may exist
	  endif

      RR   = y * state%Q
      q    = 1.d0 - sum( RR * y )  !  error in mole fraction normalization
	  call ceq_HHform( nsu, nrc, y, sys%BR, HH )
    
      v     = matmul( y, HH )  !  H'*y
	  vnorm = ceq_norm( nrc, v )
      r     = ceq_norm( nrc,  (v/vnorm - c/cnorm) )  !  residual

	  if( sys%diag >=4 ) write(sys%lu,'(8x,a,1p,3e11.4,i5)') &
	            'Newt: q, r, Srat, newts = ', q, r, Sig(nrc)/Sig(1), state%newt_its
    
	  err = r 	  
	  if( err < err_best ) then  !  record best solution
	     err_best = err
		 lam_best = lam
		 y_best   = y
	  endif
	   
      if( err + abs(q) < err_tol ) return  !  (b) convergence achieved, return with latest solution

	  if( err > sys%dec_min * err_last ) exit  ! (c)  no longer converging
      err_last = err
    
!  proceed with Newton iteration
      state%newt_its = state%newt_its + 1

!  form SVD:  HH = U * diag(Sig) * VT
      lwork = size( work )
      HHcopy = HH
      call dgesvd( 'S', 'A', nsu, nrc, HHcopy(1:nsu,1:nrc), nsu, Sig(1:nrc), &
	     U(1:nsu,1:nrc), nsu, VT(1:nrc,1:nrc), nrc, work(1:lwork), lwork, info )

      if( info /= 0 ) then
	     fail = .true. 
	     if( sys%diag >= 1 ) write(sys%lu,*)'ceq_newt: dgesvd failed, info= ', info
         return
      endif

      call ceq_Sinv( nrc, Sig, Sinv, sys%srat_lim, n_small )  ! pseudo-inverse

!  solve for lam_v:   H * lam_v = U * Sig * VT = Y;  lam_v = V * Sig^{-1} * U' * Y

      lam_v = matmul( y, U )       ! = U' * Y
	  lam_v = Sinv * lam_v         ! = Sig^{-1} * U' * Y
	  lam_v = matmul( lam_v, VT )  ! = V * Sig^{-1} * U' * Y

!  form w and solve for lam_w:  H' * H lam_w = V * Sig^2 * V' lam_w = w 
!  lam_w =  V * Sig^{-2} * V' * w

      w = c * (vnorm/cnorm) - v

      VTw   = matmul( VT, w )      ! = V' * w
	  lam_w = Sinv**2 * VTw        ! = Sig^{-2} * V' * w
	  lam_w = matmul( lam_w, VT )  ! = V * Sig^{-2} * V' * w

!  form P, v_star and w_star

      P = matmul( RR, HH )
	  v_star = sum( P * lam_v )
	  w_star = sum( P * lam_w )

	  denom  = max( v_star+w_star, 0.5d0 )

	  d_al   = (q-w_star)/denom
	  d_alp  = (q+v_star)/denom

	  dlam = d_al * lam_v + d_alp * lam_w
      lam  = lam + dlam
    
!  calculate the irreducible error
      if( n_small == 0 ) then
          err_i = 0.d0
      else
          err_i = ceq_norm( n_small, VTw(nrc-n_small+1:nrc) ) / vnorm
      endif
   end do  !  end of Newton iterations ------

!  complete convergence not achieved

   if( err_best > sys%err_huge ) then  !  no solution obtained: return with fail = .true.
      fail = .true.
	  if( sys%diag >= 5 ) write(sys%lu,*) 'ceq_newt: no solution obtained'
	  return  
   endif

!  return with best solution

   fail = .false.
   if( err_best < err ) then 
      err = err_best
      lam = lam_best
      y   = y_best
   endif

   return

end subroutine ceq_newt  !----------------------------------------------------



subroutine ceq_rate( sys, state, nsu, nrc, lam, gu, dgudt, Q, dlamdt, fail )

!  return d(lam)/dt corresponding to given d(gu)/dt  (for t=s or t=T)

   type (sys_type),  pointer :: sys
   type (state_type),pointer :: state
   integer, intent(in)       :: nsu, nrc
   real(k_wp), intent(in)    :: lam(nrc), gu(nsu), dgudt(nsu), Q(nsu)
   logical, intent(out)      :: fail
   real(k_wp), intent(out)   :: dlamdt(nrc)
 
   integer    :: n_small, lwork, info, i
   real(k_wp) :: y(nsu), gdot(nsu), RR(nsu), HH(nsu,nrc), U(nsu,nrc), VT(nrc,nrc), &
                 Sig(nrc), Sinv(nrc), work(20*(nsu+nrc)), HHinv(nrc,nsu), &
				 lamdot_g(nrc), lamdot_y(nrc), alpha, HHcopy(nsu,nrc)

   call ceq_y( sys, state, nsu, nrc, gu, lam, y, fail )

   if( fail ) then  !  failure
      dlamdt = 0.d0
	  if( sys%diag >= 5 ) write(sys%lu,*)'ceq_rate: ceq_y failed' 
      return
   endif
   fail = .false.   !  success

   gdot = y * dgudt
   RR   = y * Q

   call ceq_HHform( nsu, nrc, y, sys%BR, HH )

! HH = U * diag(Sig) * VT  is the economy-size SVD of HH

   lwork = size( work )

   HHcopy = HH
   call dgesvd( 'S', 'A', nsu, nrc, HHcopy(1:nsu,1:nrc), nsu, Sig(1:nrc), &
	    U(1:nsu,1:nrc), nsu, VT(1:nrc,1:nrc), nrc, work(1:lwork), lwork, info )

   if( info /= 0 ) then 
      fail = .true.
	  if( sys%diag >= 1 ) write(sys%lu,*)'ceq_rate: dgesvd failed, info= ', info
      return
   endif

! HHinv = V * Sinv * U' = ( U * Sinv * VT )'

   call ceq_Sinv( nrc, Sig, Sinv, sys%srat_lim, n_small )

   do i = 1, nrc  !  store Sinv * VT in VT
     VT(i,:) = Sinv(i) * VT(i,:)
   end do

   HHinv = transpose( matmul( U, VT ) )  ! = V * Sinv * U'

   lamdot_g = matmul(HHinv, gdot)
   lamdot_y = matmul(HHinv, y)

   alpha    = sum(RR * (gdot-matmul(HH,lamdot_g)) ) / &
              sum(RR * matmul(HH,lamdot_y) )
   dlamdt   = lamdot_g + alpha * lamdot_y

   return
end subroutine ceq_rate  !----------------------------------------------------

subroutine ceq_HHform( nsu, nrc, y, BR, HH )
!  form the matrix HH = diag(y) * BR

   integer, intent(in)     :: nsu, nrc
   real(k_wp), intent(in)  :: y(nsu), BR(nsu,nrc)
   real(k_wp), intent(out) :: HH(nsu,nrc)

   integer                 :: j

   do j = 1, nrc
      HH(:,j) = BR(:,j) * y
   end do

   return

end subroutine ceq_HHform  !----------------------------------------------------

subroutine ceq_Sinv( m, S, Sinv, srat_lim, n_small )

! given the m-vector of singular values, S, which are in decreasing order,
! return the m-vector of pseudo-inverses, Sinv, and the number n_small of
! small singular values.  The j-th singular value is deemed to be small
! if S(j)/S(1) < srat_lim.

   integer, intent(in)     :: m
   real(k_wp), intent(in)  :: S(m), srat_lim
   integer, intent(out)    :: n_small
   real(k_wp), intent(out) :: Sinv(m)

   integer    :: j
   real(k_wp) :: slim

   slim = srat_lim*S(1)

   do j = 1, m
      if( S(j) > slim ) then
        Sinv(j) = 1.d0/S(j)
        n_small   = m-j
      else
	    Sinv(j) = 0.d0
      endif
   end do

   return
end subroutine ceq_Sinv  !----------------------------------------------------

subroutine ceq_y( sys, state, nsu, nrc, gu, lam, y, fail )

! evaluate y = exp( [-gu + BR*lam]/2 )

   type (sys_type),  pointer :: sys
   type (state_type),pointer :: state
   integer, intent(in)       :: nsu, nrc
   real(k_wp), intent(in)    :: gu(nsu), lam(nrc)
   real(k_wp), intent(out)   :: y(nsu)
   logical, intent(out)      :: fail

   integer                   :: i

   state%nyeval = state%nyeval + 1

   y = (-gu + matmul(sys%BR,lam) )*.5  !  form log(y)

   if( maxval(y) > sys%logy_lim ) then  !  failure due to incipient overflow
      fail = .true.
	  if( sys%diag >= 5 ) write(sys%lu,*) 'ceq_y failed: max(y) = ', maxval(y)
	  return
   endif
   fail = .false.

   do i = 1, nsu
      y(i) = exp( y(i) )
   end do

   return

end subroutine ceq_y  !----------------------------------------------------

subroutine ceq_cp_eff( sys, state, ns, nsd, nsu, nrc, zu, lam, T, cp_eff, fail )

!  Evaluate the effective specific heat.
!     cp_eff = (dh/dT)/R at constant constraints and p.

   type (sys_type),   pointer :: sys
   type (state_type), pointer :: state
   integer,    intent(in)     :: ns, nsd, nsu, nrc
   real(k_wp), intent(in)     :: zu(nsu), lam(nrc), T 
   real(k_wp), intent(out)    :: cp_eff
   logical,    intent(out)    :: fail

   real(k_wp) :: gu(nsu), dgudT(nsu), dlamdT(nrc), cpor(ns), hort(ns), &
                 z(ns), Nbar, Xd(nsd), Xu(nsu), X(ns), dXudT(nsu), &
				 dXddT(nsd), dXdT(ns), w(nrc), dNdT

!  determine d(lam)/dT from dg/dT
   call ceq_g( nsu, T, state%p, sys%thermo(nsd+1:ns,:) ,gu )
   call ceq_dgdT( nsu, T,       sys%thermo(nsd+1:ns,:), dgudT )

   call ceq_rate( sys, state, nsu, nrc, lam, gu, dgudT, state%Q, dlamdT, fail )

   if( fail ) then
      if( sys%diag >= 1 ) write(sys%lu,*)'ceq_cp_eff: ceq_rate failed'
	  return
   endif

   call ceq_cp( ns, T, sys%thermo ,cpor )  !  cp and h based on all species
   call ceq_h(  ns, T, sys%thermo ,hort )
 
   z(1:ns)    = (/ state%zd(1:nsd) , zu(1:nsu) /)  !  species moles

   Nbar = sum(z)          ! total moles
   Xd   = state%zd / Nbar ! mole fractions
   Xu   =       zu / Nbar
   X    = (/ Xd , Xu /)

!  rates of change of moles and mole fractions
   dXudT  = (-dgudT + matmul(sys%BR, dlamdT) ) * Xu
   w      = matmul(Xu, sys%BR)
   dNdT   = -Nbar * sum(w * matmul(dXudT, sys%BR) ) / sum(w*w)
   dXddT  = -Xd * dNdT / Nbar
   dXdT   = (/ dXddT , dXudT /)
   cp_eff = sum(cpor*X) + T * sum(hort*dXdT)

   return
end subroutine ceq_cp_eff  !----------------------------------------------------

subroutine ceq_residual( sys, state, nsu, nrc, gu, lam, q, r, fail )  

!  compute the residuals: for diagnostics only

   implicit none
   type (sys_type), pointer   :: sys
   type (state_type), pointer :: state
   integer, intent(in)        :: nsu, nrc
   real(k_wp), intent(in)     :: gu(nsu), lam(nrc)
   real(k_wp), intent(out)    :: q, r
   logical,    intent(out)    :: fail

   real(k_wp) :: y(nsu), RR(nsu), HH(nsu,nrc), v(nrc), cu(nrc), ceq_norm

   call ceq_y( sys, state, nsu, nrc, gu, lam, y, fail )
   if( fail ) return

   RR   = y * state%Q
   q    = 1.d0 - sum( RR * y )  !  error in mole fraction normalization
   call ceq_HHform( nsu, nrc, y, sys%BR, HH )
    
   v  = matmul(y,HH)  !  H'*y
   cu = state%cr / ceq_norm( nrc, state%cr )
   r  = ceq_norm( nrc,  v / ceq_norm(nrc,v) - cu )

   return
end subroutine ceq_residual

end module ceq_solve
