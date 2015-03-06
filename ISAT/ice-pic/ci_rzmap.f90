!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ci_rzmap( r, dt, h, p, rmap, rzmap )

  ! Determine the reduced mapping r(dt) by solving only the reduced
  ! equations, and compare it to the mapping obtained by solving the
  ! full set of ODEs (both using linearized equations)
  
  ! input:
  !   r : initial reduced representation
  !  dt : time step
  !   h : enthalpy
  !   p : pressure

  ! output:
  !  rmap : reaction mapping in reduced space
  ! rzmap : reaction mapping reduced from full mapping

  use ci_utils
  use ci_dat8
  use ci_sens
  use ci_cksubs
  use ci_stats
  use ci_cem_recon

  implicit none

  real(k_dp), intent(in)  :: r(nrc), dt, h, p
  real(k_dp), intent(out) :: rmap(nrc), rzmap(nrc)

  ! local variables
  integer :: info, luout, ipiv(nrc)
  real(k_dp) :: work(nrc), stats(20), maxel
  real(k_dp) :: z_CE(ns), T_CE, S(ns), J(ns,ns), T_CEM(ns+1,nrc+2), BTTI(nrc,nrc)
  real(k_dp) :: Cz(ns+1, ns+1), Cr(nrc+1, nrc+1), eCzdt(ns+1, ns+1), eCrdt(nrc+1, nrc+1)

  ! saved quantities
  real(k_dp), save  ::  initialised = 0
  real(k_dp), save, allocatable :: r_pre(:), Cr_pre(:, :), Cz_pre(:,:)

  if(initialised == 1) then
     if( maxval(abs(r_pre - r)) == 0.d0 ) then
        Cz(ns+1, ns+1) = Cz_pre(ns+1, ns+1)
        Cr(nrc+1, nrc+1) = Cr_pre(nrc+1, nrc+1)
     else
        call initialize
     endif
  else
     call initialize
  endif

 ! exponentials
  call expm( nrc+1, Cr*dt, eCrdt, maxel )
  call expm( ns+1, Cz*dt, eCzdt, maxel )
  
  ! compute rmap and zmap
  rmap = r + eCrdt(1:nrc, nrc+1)
  rzmap = matmul( BBT, z_CE + eCzdt(1:ns, ns+1) )

contains

  subroutine initialize

    ! check if memory has been allocated
    if(initialised == 0) then
       initialised = 1
       allocate(r_pre(nrc))
       allocate(Cz_pre(ns+1,ns+1))
       allocate(Cr_pre(nrc+1,nrc+1))
    end if

    ! constrained-equilibrium composition, z_CE
    call ci_ceq( 0, r, h, p, .false., z_CE, T_CE, z_CE, T_CE, stats, info, luout )

    ! reaction rate, S
    call ciS( z_CE, T_CE, p, S )
    
    ! Jacobian, J = dS/dz
    call jacobian( 0.d0, ns, z_CE, p, T_CE, J )

    ! CEM Tangent Vectors, T_CEM
    call ci_cem_tan(z_CE, h, T_CE, p, thermo_ns, ns, nrc, BB, T_CEM, info) 

    ! Cr = [B'*J*A_CEM B'*S; 0 0] (nrc+1) x (nrc+1)
    Cr = 0.d0
    ! BTTI = inv(B'*T)
    BTTI = matmul( BBT, T_CEM(1:ns,1:nrc) )
    call dgetrf( nrc, nrc, BTTI, nrc, ipiv, info )
    call dgetri( nrc, BTTI, nrc, ipiv, work, nrc, info )
    call print_var(BTTI, "BTTI", 23)
    Cr(1:nrc, 1:nrc) = matmul(matmul(matmul(BBT,J),T_CEM(1:ns,1:nrc)),BTTI)
    Cr(1:nrc, nrc+1) = matmul(BBT, S)
    
    ! Cz = [J S; 0 0] (ns+1) x (ns+1)
    Cz = 0.d0
    Cz(1:ns, 1:ns) = J
    Cz(1:ns, ns+1) = S

    ! save the quantities
    r_pre(1:nrc) = r(1:nrc)
    Cz_pre(ns+1, ns+1) = Cz(ns+1, ns+1)
    Cr_pre(nrc+1, nrc+1) = Cr(nrc+1, nrc+1)

  end subroutine initialize
 
end subroutine ci_rzmap
