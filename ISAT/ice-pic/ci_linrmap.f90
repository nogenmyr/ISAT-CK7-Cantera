!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ci_linrmap( z, t, dz, dt, Tz, h, p, zmap, TR, A_ODE, info )

  ! Determine the linearised reaction mapping R(z+dz,t+dt), given z, t, dz and dt.
  ! Use Prof. Van Loan's solution using block matrices for both +ve and -ve dt values.
  
  ! input:
  !   z : initial point
  !   t : initial time
  !  dz : finite increment in z
  !  dt : finite increment in t
  !  Tz : temperature at z
  !   h : enthalpy
  !   p : pressure

  ! output:
  !  zmap : linearised reaction mapping
  !    TR : temperature at zmap
  ! A_ODE : gradient matrix
  !  info >  0, success
  !       =  0, clean successful retrieve
  !       =  1, within S_tol, so using saved quantities, values reset
  !       =  2, called ci_ice_chem_map, initialized, successful
  !       =  3, initial linear mapping failed, ci_linrmap recalled, values reset
  !       =  4, some species were negative |z| < |z_ntol|, reset to 0.d0
  !       <  0, failure
  !       = -1, Mapping R(z+dz,t+dt) failed, due to ci_ice_chem_map failure
  !       = -2, T/h out of range in temphy

  use ci_utils
  use ci_dat8
  use ci_sens
  use ci_cksubs
  use ci_ice_cksubs
  use ci_stats

  implicit none

  real(k_dp), intent(inout) :: z(ns), t, dz(ns), dt, Tz ! reset for info = 1, 3
  real(k_dp), intent(in)    :: h, p
  real(k_dp), intent(out)   :: zmap(ns), TR, A_ODE(ns+1,ns+1)
  integer, intent(out)      :: info

  ! local variables
  integer    :: i, flag, check = 1, t_flag
  real(k_dp) :: z_ntol = -1.d-30 ! Negative species concentration tolerance.
  real(k_dp) :: S_tol  = 0.d0   ! Acceptable error tolerance in S, | S - (So + J*dz) | 
  real(k_dp) :: S_err, zo(ns), To, A(ns,ns), S(ns)
  real(k_dp) :: J(ns,ns), C(ns+1, ns+1), eCdt(ns+1, ns+1)
  real(k_dp) :: dzs(ns), zr(nsp1,nsp4), maxel

  ! saved quantities
  real(k_dp), save  ::  t_pre, initialised = 0
  real(k_dp), save, allocatable :: zo_pre(:), z_pre(:), A_ODE_pre(:, :), C_pre(:, :)

  ! stats
  call routine_start(i_ci_linrmap)
  
  info = 0 ! anticipate success

  if(initialised == 1) then
     if(maxval( abs( z(1:ns) - z_pre(1:ns) ) ) == 0.d0 .and. (t == t_pre)) then
        ! copy saved values
        zo(1:ns) = zo_pre(1:ns)
        A(1:ns,1:ns) = A_ODE_pre(1:ns, 1:ns)
        A_ODE(1:ns+1, 1:ns+1) = A_ODE_pre(1:ns+1, 1:ns+1)
        C(1:ns+1, 1:ns+1) = C_pre(1:ns+1, 1:ns+1)
     else
        ! compute S
        call ciS( z, Tz, p, S )

        ! compute error using the saved Jacobian, J
        ! S_err = norm(S - (S_pre + J_pre*(z - z_pre)))

        S_err = norm(S - (C_pre(1:ns, ns + 1) + matmul(C_pre(1:ns, 1:ns), (z - z_pre))))
        
        if(S_err <= S_tol) then
           !ci_linrmap: using old saved values
           info = 1

           ! copy values
           zo(1:ns) = zo_pre(1:ns)
           A(1:ns,1:ns) = A_ODE_pre(1:ns, 1:ns)
           A_ODE(1:ns+1, 1:ns+1) = A_ODE_pre(1:ns+1, 1:ns+1)
           C(1:ns+1, 1:ns+1) = C_pre(1:ns+1, 1:ns+1)

           ! reset z, dz, t and dt
           dz = z - z_pre + dz
           dt = t - t_pre + dt
           z = z_pre
           t = t_pre
        else
           call initialize
           if(info == -1) return
        endif
     endif
  else
     call initialize
     if(info == -1) return
  endif
  
  call get_zmap

  if(maxval(zmap) > 1.d0 .or. minval(zmap) < z_ntol) then
     ! failure, reset values and recalculate

     z = z + dz
     t = t + dt
     if(t < 0.d0) then
        info = -1
        return
     endif
     dz = 0.d0
     dt = 0.d0
     call temphz( h, z, Tz, check, t_flag )
     if( t_flag /= 0 ) then
        flag = -2
        return
     endif
     call initialize
     if(info == -1) return
     call get_zmap
     info = 3
  endif

  ! make sure all the components of zmap are positive
  do i=1,ns
     if(zmap(i) < 0.d0) then 
        zmap(i) = 0.d0
        info = 4
     endif
  end do


  ! Enforce normalization
  zmap = zmap / dot_product( zmap, amolwt )

  ! Compute TR at zmap
  call temphz( h, zmap, TR, check, t_flag )
  if( t_flag /= 0 ) then
     flag = -2
     return
  endif

  call routine_stop(i_ci_linrmap)

  return

contains
  
  subroutine initialize

    ! check if memory has been allocated
    if(initialised == 0) then
       initialised = 1
       allocate(z_pre(ns))
       allocate(zo_pre(ns))
       allocate(A_ODE_pre(ns+1,ns+1))
       allocate(C_pre(ns+1,ns+1))
       z_pre = 0
    end if

    ! calculate new values and save
    info = 2

    ! reaction rate: S
    call ciS( z, Tz, p, S )
    
    ! zo = R(z,t) and sensitivity matrix: A
    call ci_ice_chem_map( 1, t, p, z, h, To, zo, A_ODE, flag)
    if( flag > 0 ) then
       info = -1
    endif
    A(1:ns, 1:ns) = A_ODE(1:ns, 1:ns)
    
    ! Jacobian, J = dS/dz
    call jacobian( t, ns, z, p, Tz, J )
    
    ! Form the special matrix C = [J S; 0 0] (ns+1) x (ns+1)
    C = 0.d0
    C(1:ns, 1:ns) = J
    C(1:ns, ns+1) = S

    ! save the quantities
    t_pre = t
    z_pre(1:ns) = z(1:ns)
    zo_pre(1:ns) = zo(1:ns)
    A_ODE_pre(1:ns+1,1:ns+1) = A_ODE(1:ns+1, 1:ns+1)
    C_pre(1:ns+1, 1:ns+1) = C(1:ns+1, 1:ns+1)

  end subroutine initialize
  
  subroutine get_zmap

    ! Prof. Van Loan's solution using block matrix exponential.
    ! C = [J, S,; 0, 0] expm(C*dt) = eCdt
    if(dt == 0.d0) then
       zmap = zo + matmul(A, dz)
    else
       call expm( ns+1, C*dt, eCdt, maxel )
       dzs = matmul(eCdt(1:ns, 1:ns), dz) + eCdt(1:ns, ns+1)
       zmap = zo + matmul( A, dzs )
    endif
  end subroutine get_zmap
  
end subroutine ci_linrmap
