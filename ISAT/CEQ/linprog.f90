!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module linprog

!  The subroutine  lp  solves linear programming problems using the simplex method.
!  See  lp  (below) for detailed descriptions.

!  S.B. Pope 10/26/2006
!  Additions made 3/25/2009: added diagnostic output and optional arguments:
!    indic_pos, lhvars, t_z_stop, t_z_accept, t_pivot, t_neg, t_aof.
!  If these arguments are not used, then the result is indentical to the original version. 

implicit none
integer, parameter, private :: k_prec   = kind(1.d0)  ! kind of reals
integer, parameter, private :: iter_def = 10          ! iter_max=iter_def*(nx+nc)
real(k_prec), private       :: tol_def  = 1.d-9       ! default tolerance
real(k_prec), private       :: tol                    ! tolerance
real(k_prec), private       :: tol_z_stop, tol_z_accept, tol_pivot, tol_neg, &
                               tol_aof      

! set lp_diag=0 for no diagnostic output; set lp_diag = 1, 2, 3 for increasing o/p
integer, private            :: lp_diag   = 0 

! if call_diag > 0, then o/p (with lp_diag=3) on call==call_op only
integer, private            :: call_diag = 0 
integer, private            :: lu = 91
logical, private            :: lu_set=.false.

contains

subroutine lp( nx, nle, nge, neq, ale, age, aeq, ble, bge, beq, f, x, info,  & 
               iter_lim, toler, check, feasible_only, stats, indic_pos, lhvars, &
               t_z_stop, t_z_accept, t_pivot, t_neg, t_aof, lp_d, lu_d  )

!  Solve linear program.

!  Determine the nx-vector x which minimize g = f^T * x, subject to
!    the nle inequality constraints   Ale * x <= ble, 
!    the nge inequality constraints   Age * x >= bge, and
!    the neq equality constraints     Aeq * x  = beq,
!  where the inequalities apply to each component. 
!  The right-hand sides (ble, bge and beq) need not be positive.

!  On return:
!	info =  0 if solution is obtained 
!	info = -1 if g is unbounded
!	info = -2 if there is no feasible solution
!	info = -3 iteration limit exceeded in attempting to obtain feasible solution
!	info = -4 iteration limit exceeded in attempting to optimize solution
!   info = -5 constraint residuals exceed  check (requires that check is present)
!   info >  0 invalid input

!  Optional arguments

!   iter_lim       - upper limit on the number of iterations allowed
!   toler          - small tolerance
!   check          - check that constraint rediduals are less than  check 
!   feasible_only  = .true. if only feasible solution required (i.e., no minimization)
!   stats          - statistics returned of optimization performed
!                  = (/ nx, nle, nge, neq, nc, n_le, n_ge,  &
!                       info, iter_max, iter_feas, iter_opt, tol, err(1:5) /)
!                    err = (/ minval(x), abs residual in LE, GE, EQ, maxval(err(1:4)) /)
!   indic_pos      - indic_pos(i)=1 if x(i) > 0 for all optimal solutions; 
!                  - indic_pos(i)=0 otherwise, i.e., if x(i)=0 for all optimal solutions.
!   lhvars         - indices of left-hand variables
!   t_z_stop, t_z_accept, t_pivot, t_neg, t_aof - specific tolerances
!   lp_d           - set lp_diag=lp_d to control diagnostic output
!   lu_d           - logical unit to be used for diagnostic output (else unit 91 is used)


   integer, intent(in)       :: nx   ! number of variables 
   integer, intent(in)       :: nle  ! number of (LE) inequality constraints
   integer, intent(in)       :: nge  ! number of (GE) inequality constraints    
   integer, intent(in)       :: neq  ! number of (EQ) equality constraints

   real(k_prec), intent(in)  :: ale(nle,nx), age(nge,nx), aeq(neq,nx)  ! constraint matrices
   real(k_prec), intent(in)  :: ble(nle), bge(nge), beq(neq)           ! constraint vectors
   real(k_prec), intent(in)  :: f(nx)  ! objective function vector

   real(k_prec), intent(out) :: x(nx)  ! solution (or 0 if no solution)
   integer,      intent(out) :: info   ! status of solution

   integer,      intent(in),  optional :: iter_lim, lp_d, lu_d
   real(k_prec), intent(in),  optional :: toler
   real(k_prec), intent(in),  optional :: t_z_stop, t_z_accept, t_pivot, t_neg, t_aof
   real(k_prec), intent(in),  optional :: check
   logical,      intent(in),  optional :: feasible_only
   real(k_prec), intent(out), optional :: stats(20)
   integer,      intent(out), optional :: indic_pos(nx)
   integer,      intent(out), optional :: lhvars(nle+nge+neq)

!-----------------------

   real(k_prec), allocatable :: tab(:,:)     ! tableau
   real(k_prec) :: err(5)                    ! residual errors
   integer      :: rhv(nx)                   ! indexes of right-hand (RH) variables
   integer      :: lhv(nle+nge+neq)          ! indexes of left-hand (LH) variables

   real(k_prec), save :: calls = 0.d0 
   integer      :: nc, i, j, ii, n_le, n_ge, iter_max, iter_feas, iter_opt
   logical      :: feas_only

!------- specify diagnostic output  

   calls = calls + 1.d0
   
   if( call_diag > 0 ) then
      lp_diag = 0
      if( nint(calls) == call_diag ) lp_diag = 3
   endif
   
   if( present(lp_d) ) lp_diag = lp_d
   
   if( present( lu_d ) .and. .not.lu_set ) then
      lu     = lu_d
      lu_set = .true.
   endif
 
!-------- set output variables to preliminary settings
  
   info  = 0
   x     = 0.d0

   if( present( stats ) )   stats     = 0.d0
   if( present(indic_pos) ) indic_pos = 0

   nc   = nle + nge + neq  !  total number of constraints

!   check the input
   if( nx <= 0  .or.  nle < 0  .or.  nge < 0  .or.  neq < 0 ) then
      info = 1
      return
   endif

!--------  set tol, iter_max and feas_only
   if( present( toler ) ) then  !  error tolerance

      if( toler <= 0.d0 ) then
         info = 2
         return
      endif
      tol = toler

   else
      tol = tol_def
   endif
   
!  set specific tolerances based on tol
   tol_z_stop   = tol
   tol_z_accept = tol 
   tol_pivot    = tol  
   tol_neg      = huge(0.d0)
   tol_aof      = tol
   
!  over-ride if user-specified
  
   if( present( t_z_stop ) )   tol_z_stop   = t_z_stop
   if( present( t_z_accept ) ) tol_z_accept = t_z_accept
   if( present( t_pivot ) )    tol_pivot    = t_pivot
   if( present( t_neg ) )      tol_neg      = t_neg
   if( present( t_aof ) )      tol_aof      = t_aof

   if( present( iter_lim ) ) then  !  limit on iterations

      if( iter_lim < 0 ) then
         info = 3
         return
      endif

      iter_max = iter_lim
   else
      iter_max = iter_def * ( nx + nc )
   endif

   if( present( feasible_only ) ) then  !  only feasible solution required?
      feas_only = feasible_only
   else
      feas_only = .false.  !  by default, perform optimization
   endif

!------------  set tableau

   allocate( tab(0:nc+1,0:nx) )
   tab = 0.d0
   tab(0,1:nx) = -f  !  objective function ( -ve for minimization) row=0 

!  Change sign of equations if necessary in order to make RHS non-negative.
!  If sign is changed, LE becomes GE and vice versa.
!  Order equations in tableau: LE, GE, EQ.

   n_le = 0  !  count new n_le
   n_ge = 0  !  count new n_ge

   do ii = 1, nle  !  original LE constraints

     if( ble(ii) >= 0.d0 ) then  ! LE with no sign change
        n_le        =  n_le + 1
        i           =  n_le
        tab(i,0)    =  ble(ii)
        tab(i,1:nx) = -ale(ii,1:nx)
     else                       ! GE with sign change
        n_ge        =  n_ge + 1
        i           =  nle + nge + 1 - n_ge
        tab(i,0)    = -ble(ii)
        tab(i,1:nx) =  ale(ii,1:nx)
     endif

   end do

   do ii = 1, nge  !  original GE constraints

     if( bge(ii) <= 0.d0 ) then  ! LE with sign change
        n_le        =  n_le + 1
        i           =  n_le
        tab(i,0)    = -bge(ii)
        tab(i,1:nx) =  age(ii,1:nx)
     else                       ! GE with no sign change
        n_ge        =  n_ge + 1
        i           =  nle + nge + 1 - n_ge
        tab(i,0)    =  bge(ii)
        tab(i,1:nx) = -age(ii,1:nx)
     endif

   end do 
   
   do ii = 1, neq  !  equality EQ constraints
      i  = nle + nge + ii

      if( beq(ii) >= 0.d0 ) then  !  EQ with no sign change
         tab(i,0)    =  beq(ii)
         tab(i,1:nx) = -aeq(ii,1:nx)   
      else                        !  EQ with sign change
         tab(i,0)    = -beq(ii)
         tab(i,1:nx) =  aeq(ii,1:nx)   
      endif
   end do
   
!---------  attempt to solve LP   

   
   if( lp_diag > 0 ) then
      rhv = (/ (i,i=1,nx) /)
      lhv = (/ (i,i=1,nc) /)
      call tab_op( calls, -3, 0, tab, nc, nx, rhv, lhv )
   endif

   call lp_simplex( tab, nc, nx, n_le, n_ge, neq, info, rhv, lhv,  &
                    feas_only, iter_max, iter_feas, iter_opt )
                    
   if( info /=0 .and. lp_diag > 0 ) then
      write(0,*)'linprog failed: info, calls = ', info, nint(calls)
      call tab_op( calls, -7, info, tab, nc, nx, rhv, lhv )
   endif
                    
   if( lp_diag > 0 ) call tab_op( calls, -1, info, tab, nc, nx, rhv, lhv )

!---------  obtain solution from LH variables     

   do j = 1, nc  
      i = lhv(j)
      if( i <= nx ) x(i) = tab(j,0)
   end do
   
   if( present( lhvars ) ) lhvars = lhv
    
!  determine components of x which are essentially positive  
   if( info == 0  .and.  present(indic_pos) ) &
       call lp_pos( tab, nc, nx, rhv, lhv, indic_pos )
      
   if( present( check ) )  then  !  check constraint equations
      call lp_check( nx, nle, nge, neq, ale, age, aeq, ble, bge, beq, x, err )
      if( info == 0  .and.  err(5) > check ) info = -5
   else
      err = 0.d0
   endif
                               
   if( present( stats ) ) then  ! set stats
      stats(1:7)   = (/ nx, nle, nge, neq, nc, n_le, n_ge /)*1.d0
      stats(8:11)  = (/ info, iter_max, iter_feas, iter_opt /)*1.d0
      stats(12:17) = (/ tol, err /) 
   endif
   
   deallocate( tab )
   return
end subroutine lp

subroutine lp_pos( tab, nc, nx, rhv, lhv, indic_pos ) !---------------------------

!  set indic_pos(i) = 1 for all components of x which are not essentially zero,
!           i.e., those that can be positive in the optimal feasible solution.
!  Note: no tolerances are used.
!  In case of failure(i.e., the positive variables cannot be determined unambiguously)
!     this is flagged by indic_pos(1) < 0.

  real(k_prec), intent(in)  :: tab(0:nc+1,0:nx)
  integer,      intent(in)  :: nc, nx, rhv(nx), lhv(nc)
  integer,      intent(out) :: indic_pos(nx)
  
  integer :: i, j, nz_lhv, rhs_pos(nx), lhs_pos(nc), rhs_eqn(nx), &
             tries, changes, n_pos
  
  indic_pos = 0
  
  ! provisionally, set rhs_pos(j)=1 only if the j-th RHS var
  ! does not affect objective function
  
  rhs_pos = 0
  do j = 1, nx
     if( tab(0,j) == 0.d0 ) rhs_pos(j) = 1
  end do
  
  ! provisionally, set lhs_pos(i)=1,  This is correct if the i-th LHS variable
  ! in the optimal feasible solution obtained is positive.
  lhs_pos = 1
  
  ! Treatment of LHS variables which are zero in the optimal feasible solution.
  ! These may be non-zero if their equation contains a positive coefficient
  ! corresponding to a RHS var which does not affect the objective function.
  
  ! Note: this treatment is not fool-proof.  If more than one LHS var
  !       is essentially zero, then there may be incompatible requirement
  !       on the RHS vars deemed positive.  Thus some essentially zero 
  !       RHS vars may not be identified.
 
  ! loop over modifications (of _pos from 1 to 0) until there is no further change.
  do tries = 1, nx+nc  
     changes = 0
     
     do i = 1, nc  !  loop over LHS variable
        if( tab(i,0) /= 0.d0 ) cycle  ! positive LHS var
        if( lhs_pos(i) == 0  ) cycle  ! already determined to be essentially zero
        !  zero LHS set potentially positive
        !  is there a positive coeff of a pos RHS var?
        n_pos = 0
        do j = 1, nx       ! loop over pos RHS
           if( rhs_pos(j) == 0 ) cycle
           if( tab(i,j) > 0.d0 ) then
              n_pos = 1
              exit
           endif
        end do
           
        if( n_pos == 0 ) then ! no positive coeffs
           lhs_pos(i) = 0 ! LHS var is essentially zero
           changes    = 1
              
           do j = 1, nx  ! rhs_pos(j)=1 only if tab(i,j)=0.d0
              if( rhs_pos(j) == 1  .and.  tab(i,j) < 0.d0 ) rhs_pos(j) = 0
           end do
        end if
        
     end do  !  end of loop over LHS vars
     if( changes == 0 ) exit ! all done
  end do
  
  if( changes > 0 ) then  !  failed
     if(  lp_diag > 0 ) then 
        write(0,*)'lp_pos failed:, changes = ', changes
        write(0,*)lhs_pos
        write(0,*)rhs_pos
     endif
     
     indic_pos(1) = -1  ! indicate failure
     return
  endif
  
  !  check for unambiguous determination
  rhs_eqn = 0 ! number of times an RHS var enters an lhs_pos=0 
              ! equation with a non-zero coefficient 
  do i = 1, nc
     if( lhs_pos(i) == 0 ) then
        do j = 1, nx
           if( rhs_pos(j)==1  .and.  tab(i,j) /=0.d0 )  &
               rhs_eqn(j) = rhs_eqn(j) + 1
        end do
     endif
  end do
  
  if( maxval(rhs_eqn) > 1 ) then
     if(  lp_diag > 0 ) write(0,*)'*** lp_pos ***, ambiguous: ', maxval(rhs_eqn)
     indic_pos(1) = -1
     return
  endif
           
  !  success: form indic_pos
  
  do j = 1, nc  ! left-hand variables
     i = lhv(j)
     if( i <= nx  .and.  lhs_pos(j) == 1 ) indic_pos(i) = 1
  end do
  
  do j = 1, nx  ! right-hand variables
     i = rhv(j)
     if( i <= nx  .and.  rhs_pos(j) == 1 ) indic_pos(i) = 1
  end do
  
  return
end subroutine lp_pos

subroutine tab_op( calls, k, info, tab, nc, nx, rhv, lhv ) !---------------------------

!  output tableau

  real(k_prec), intent(in)    :: calls, tab(0:nc+1,0:nx)
  integer, intent(in)         ::k, info, nc, nx, rhv(nx), lhv(nc)
  
  integer :: i, j

  if( .not.lu_set ) then
     open( lu, file = 'lp_diag.op' )
     lu_set = .true.
  endif
  
  write(lu,*)' '
  write(lu,'(a,3i6)')'tab_op for calls, k, info = ', nint(calls), k, info
  write(lu,*)' '

  write(lu,'(13x,100i8)')      (j,j=1,nx)
  write(lu,'(13x,100i8)')      (rhv(j),j=1,nx)
  write(lu,'(2i4,1p,100e8.0)') 0, 0, (tab(0,j),j=0,nx)
  do i = 1, nc
     write(lu,'(2i4,1p,100e8.0)') i, lhv(i), (tab(i,j),j=0,nx)
  end do
  write(lu,'(2i4,1p,100e8.0)') nc+1, 0, (tab(nc+1,j),j=0,nx)
    
  return
end subroutine tab_op

subroutine lp_simplex( tab, nc, nx, nle, nge, neq, info, rhv, lhv,  &
                       feas_only, iter_max, iter_feas, iter_opt )

  real(k_prec), intent(inout) :: tab(0:nc+1,0:nx)
  integer, intent(in)         :: nc, nx, nle, nge, neq, iter_max
  integer, intent(inout)      :: info, rhv(nx), lhv(nle+nge+neq)
  logical, intent(in)         :: feas_only 
  integer, intent(out)        :: iter_feas, iter_opt
      
  integer      :: i, j, ip, jp, ncol, col(nx)
  real(k_prec) :: tmax

  info      = 0  ! anticipate success
  iter_feas = 0  ! iterations to obtain feasible solution
  iter_opt  = 0  ! iterations to optimize solution
   
!  col  - list of columns (active RH variables) to be considered for exchange   [col]
!  ncol - number of columns (active RH variables) to be considered for exchange [ncol]
!  Initially the n original variables (x) are, in order, the RH variables.
!  All n columns are active.
  
  ncol = nx        !  number of columns (active RH variables)
  do j = 1, nx
    col(j) = j     ! index of active columns
    rhv(j) = j     ! index for RH variables
  end do

!  Initially the nc LH variables are, in order:
!     nle slack variables
!     nge slack variables (with associated auxilliary variables)
!     neq auxilliary variables

  do i = 1, nc   
    lhv(i) = nx+i  ! set indexes of LH variables
  end do
  
  if( nge+neq /= 0 ) &
    call lp_feasible( tab, nc, nx, nle, nge, neq, rhv, lhv, col, ncol, info,  &
                      iter_max, iter_feas  )
                      
  if( lp_diag > 0 ) call tab_op( 0.d0, -2, info, tab, nc, nx, rhv, lhv )
            
  if( feas_only  .or.  info /= 0 ) return
     
  do !----  loop over exchanges to optimize
      
!  find the largest element in the z-row (i=1)
    call lp_max( tab, nc, nx, 0, col, ncol, jp, tmax )
      
    if( tmax <= tol_z_stop ) then  !  largest element is essentially zero
      info=0                !  objective function cannot be increased
      return                !  success
    endif

!  find pivot:  given RH variable jp, pivot is LH variable ip  
    call lp_pivot( tab, nc, nx, ip, jp)
      
    if( ip == 0 ) then
      info = -1  ! unbounded.  No limit on increase in RH variable
      if( tmax < tol_z_accept ) info = 0 ! accept solution
      if( lp_diag >= 2  .and.  info /=0 ) call tab_op( -8.d0, iter_opt, info, tab, nc, nx, rhv, lhv )
      return
    endif

! exchange RH (jp) and LH (ip) variables
    call lp_exchange( tab, nc, nx, 0, ip, jp)
    j       = rhv(jp)
    rhv(jp) = lhv(ip)
    lhv(ip) = j
    
    iter_opt = iter_opt + 1
    if( lp_diag >= 2 ) call tab_op( -4.d0, iter_opt, info, tab, nc, nx, rhv, lhv )

    if( iter_opt > iter_max ) then
       info = -4
       return
    endif
  end do 
end subroutine lp_simplex
         
subroutine lp_feasible( tab, nc, nx, nle, nge, neq, rhv, lhv, col, ncol, info,  &
                        iter_max, iter_feas  )
 
!  Attempt to obtain a basic feasible solution by maximizing the
!  auxilliary objective function (AOF).
!  The tableau is modifiec accordingly, and the status of the solution is returned
!  in  info.
     
  real(k_prec), intent(inout) :: tab(0:nc+1,0:nx)
  integer, intent(in)         :: nc, nx, nle, nge, neq, iter_max
  integer, intent(inout)      :: rhv(nx), lhv(nc), col(ncol), ncol, info
  integer, intent(out)        :: iter_feas

  integer      :: i, j, ip, jp, ii, jj
  real(k_prec) :: aof, tmax
  logical      :: exist_aux, ge_swapped(nc)

!  indicate that GE slack variables have not been exchanged
  ge_swapped(1:nge) = .false.

  do  j = 0, nx  !  complete the bottom row of the tableau
                 !  = -sum( rows of nge and neq variables )
    tab(nc+1,j) = - sum( tab(nle+1:nc,j) )
  end do

  iter_feas = 0  !  number of iterations

  do  !--------- loop over exchanges
    if( iter_feas > iter_max ) then
       info = -3
       return
    else
       iter_feas = iter_feas + 1
    endif
    
    if( lp_diag >= 3 ) call tab_op( -5.d0, iter_feas, 0, tab, nc, nx, rhv, lhv )

    aof = tab(nc+1,0)  ! current value of AOF
!  find largest element in the AOF (bottom) row
!  Note: when feasible solution is obtained, bottom row consists of 0 and negative integers.
    call lp_max( tab, nc, nx, nc+1, col, ncol, jp, tmax )  

    if( tmax <= tol_aof .and. aof <= tol_aof ) then 
 !  The auxilliary objective function (AOF) is around zero or less
 !  and cannot be increased. 
  
      if( aof < -tol_aof )then
 ! AOF is strictly negative -- no feasible solution
        info = -2  
        return
      endif
 
 ! AOF is essentially zero ( |AOF| < tol )
 ! Hence, feasible solution found
 ! Clean out auxilliary variables
              
      exist_aux = .false.  !  assume no auxilliary variables
        
      do ip = nle + nge + 1, nc      ! loop over EQ auxilliary variables
        if( lhv(ip) /= ip+nx ) cycle ! aux var has already been exchanged
          call lp_max_abs( tab, nc, nx, ip, col, ncol, jp, tmax )
          if(tmax > tol_aof) then
            exist_aux = .true.       ! exchange required
            exit
          endif                         
      end do

      if( .not.exist_aux ) then
        do  i = nle+1, nle+nge              !  loop over GE aux vars 
          if( .not.ge_swapped(i-nle) ) then !  if not exchanged...
            tab(i,0:nx) = -tab(i,0:nx)      !  ...change sign of row
          endif
        end do
        exit  !  done: no more auxilliary variables
      endif
               
    else      
! tmax is greater than tol  (assuming AOF < tol)

      if( aof > tol_aof  .and.  lp_diag > 0 ) then
         write(0,*)'linprog: positive aof = ', aof
      endif
      
      call lp_pivot( tab, nc, nx, ip, jp )  !  find pivot (jp)
      if( ip == 0 ) then
        info = -1  ! no pivot - AOF is unbounded
        return     ! no feasible solution
      endif
    
    endif
 
!  exchange LH variable ip with RH variable jp
    call lp_exchange( tab, nc, nx, 1, ip, jp )  

    if( lhv(ip) >= nx+nle+nge+1) then
! ip corresponds to an EQ auxilliary variable. 
! It has become a RH variable and does not need to be considered further.
! It has been exchanged to jp, so remove the jp column from consideration. 
! Find j such that col(j)=jp   
      do  j = 1, ncol
        if( col(j) == jp ) exit  ! j found
      end do
      ncol = ncol-1  !  reduce length of list
!  remove col(j) by moving rest of list up
      do jj = j, ncol
        col(jj) = col(jj+1)
      end do

    else
      
      ii = lhv(ip)-nle-nx
      if( ii >= 1 ) then 
!  ip corresponds to a GE slack variable
        if( .not.ge_swapped(ii) )then
          ge_swapped(ii) = .true.  !  indicate that it has been exchanged
!  account for implicit artificial variable
          tab(nc+1,jp)   =   tab(nc+1,jp)+1.d0 
          tab(0:nc+1,jp) = - tab(0:nc+1,jp)
        endif
      endif
        
    endif
      
    jj      = rhv(jp)
    rhv(jp) = lhv(ip)  !  make exchange in indices
    lhv(ip) = jj
      
  end do  

  return
end subroutine lp_feasible

subroutine lp_max( tab, nc, nx, i, col, ncol, jp, tmax )

! col(j), j=1:ncol are the indexes of the ncol active right-hand
! variables.
! i is the index of a left-hand variable (or objective function)
! return jp = col(j) (1 <= j <= ncol) and tmax such that
! tmax = tab(i,jp)  is the largest element of  tab
! corresponding to i and col(:)

  real(k_prec), intent(in)  :: tab(0:nc+1,0:nx)
  integer, intent(in)       :: nc, nx, i, col(ncol), ncol
  integer, intent(out)      :: jp
  real(k_prec), intent(out) :: tmax

  integer :: j
 
  if( ncol <= 0 ) then
    tmax = 0.d0  ! return if no right-hand variables
    jp   = 0
    return
  endif
      
  tmax = -huge(1.d0)
  do j = 1, ncol
    if( tab(i,col(j)) <= tmax ) cycle
    tmax = tab(i,col(j))
    jp   = col(j)       
  end do

  return
end subroutine lp_max
      
      
subroutine lp_max_abs( tab, nc, nx, i, col, ncol, jp, tmax )

! col(j), j=1:ncol are the indexes of the ncol active right-hand
! variables.
! i is the index of a left-hand variable (or objective function)
! return jp = col(j) (1 <= j <= ncol) and tmax such that
! tmax = tab(i,jp)  is the largest element of  tab (in absolute magnitude)
! corresponding to i and col(:)

  real(k_prec), intent(in)  :: tab(0:nc+1,0:nx)
  integer, intent(in)       :: nc, nx, i, col(ncol), ncol
  integer, intent(out)      :: jp
  real(k_prec), intent(out) :: tmax

  integer :: j
 
  if( ncol <= 0 ) then
    tmax = 0.d0  ! return if no right-hand variables
    jp   = 0
    return
  endif
      
  tmax = -huge(1.d0)
  do j = 1, ncol
    if( abs( tab(i,col(j)) ) <= tmax ) cycle
    tmax = abs( tab(i,col(j)) )
    jp   = col(j)       
  end do

  tmax   = tab(i,jp)
  return
end subroutine lp_max_abs

subroutine lp_pivot( tab, nc, nx, ip, jp )
      
!  Given the right-hand variable  jp, identify the pivot ip.
!  If there is no valid pivot, return ip = 0.

  real(k_prec), intent(in) :: tab(0:nc+1,0:nx)
  integer, intent(in)      :: nc, nx, jp
  integer, intent(out)     :: ip

  integer :: i,j
  real(k_prec) :: dxmin, dxi, disc
      
  ip    = 0          
  dxmin = huge(1.d0)
      
  do i = 1, nc 
    if( tab(i,jp)  >= -tol_pivot ) cycle !  skip non-negative components
        
    dxi = -tab(i,0) / tab(i,jp) ! limiting value of RH variable jp, 
                                !  imposed by LH variable i 
          
    if( dxi < dxmin ) then
            ip    = i  ! i is limiting and therefore pivot
            dxmin = dxi
            
    else if ( dxi == dxmin ) then  !  treat degeneracy: Bland's rule
      do j = 1, nx  
        disc = tab(ip,j)/tab(ip,jp) - tab(i,j)/tab(i,jp)
        
        if( disc > 0.d0 ) then
           exit
        elseif( disc < 0.d0 ) then
           ip = i
           exit
        endif
      end do

    endif
                
  end do
  
  return
end subroutine lp_pivot
      
subroutine lp_exchange( tab, nc, nx, if_aux, ip, jp)
    
!  Update the tableau based on exchanging the LH variable ip
!  with the RH variable jp.
!  If the auxilliary problem is being solved, set if_aux=1
!  so that the final row is included.  Otherwise, set if_aux=0.

  real(k_prec), intent(inout) :: tab(0:nc+1,0:nx)
  integer, intent(in)         :: nc, nx, if_aux, ip, jp
  
  integer      :: i
  real(k_prec) :: piv_inv, piv_col(0:nc+1), tab0
      
  piv_inv         = 1.d0 / tab(ip,jp)         ! inverse of pivot
  piv_col(0:nc+1) = tab(0:nc+1,jp) * piv_inv  !  new pivot column
  
!  Adjust piv_col if necessary to prevent negative left-hand variables.
!  This may be necessary if the true pivot is not selected because it is smaller
!  than the tolerance. 
!  To guarantee non-negative left-hand variables, set tol_neg = 0.d0  
!  To suppress, set tol_neg = huge(0.d0).  SBP 3/5/09

  do i = 1, nc
     if( i==ip ) cycle
     tab0 = tab(i,0) - tab(ip,0) * piv_col(i) 
     if( tab0 < -tol_neg ) then
        tab(i,jp)  = 0.d0  !  deemed to be zero within roud-off
        piv_col(i) = 0.d0
        tab(i,0)   = 0.d0  !  set potentially negative LHS variable to zero
     endif
  end do
      
  do i = 0, nc + if_aux   !  loop over all rows (except pivot row)  
    if( i /= ip ) tab(i,0:nx) = tab(i,0:nx) - tab(ip,0:nx) * piv_col(i)             
  end do

  tab(0:nc+1,jp) =  piv_col(0:nc+1)       ! set pivot column
  tab(ip,0:nx)   = -tab(ip,0:nx)*piv_inv  ! set pivot row
  tab(ip,jp)     =  piv_inv               ! set pivot element
      
  return
end subroutine lp_exchange   

subroutine lp_check( nx, nle, nge, neq, ale, age, aeq, ble, bge, beq, x, err )

!  measure residuals in constraints

!  return in err(1:5)
!     err(1) = error in x(i) >=0
!     err(2) = error in LE constraints
!     err(3) = error in GE constraints
!     err(4) = error in EQ constraints
!     err(5) = max( err(1:4) )

   integer, intent(in)       :: nx, nle, nge, neq

   real(k_prec), intent(in)  :: ale(nle,nx), age(nge,nx), aeq(neq,nx)  
   real(k_prec), intent(in)  :: ble(nle), bge(nge), beq(neq), x(nx)    
   real(k_prec), intent(out) :: err(5) 
   
   real(k_prec) :: zle(nle), zge(nge), zeq(neq)

   err    = 0.d0
   err(1) = max( -minval(x), 0.d0 )  !  max error in x(i) >= 0
   
   if( nle > 0 ) then  !  error in LE constraint
      zle(1:nle) = ble(1:nle) - matmul( ale(1:nle,1:nx) , x(1:nx) )
      err(2) = max( 0.d0,  -minval( zle(1:nle) ) )
   endif
   
   if( nge > 0 ) then  !  error in GE constraint
      zge(1:nge) = bge(1:nge) - matmul( age(1:nge,1:nx) , x(1:nx) )
      err(3) = max( 0.d0,  maxval( zge(1:nge) ) )
   endif
   
   if( neq > 0 ) then  !  error in EQ constraint
      zeq(1:neq) = beq(1:neq) - matmul( aeq(1:neq,1:nx) , x(1:nx) )
      err(4) = max( maxval(zeq), -minval( zeq(1:neq) ) )
   endif
   
   err(5) = maxval( err(1:4) )
   
end subroutine lp_check

end module linprog
