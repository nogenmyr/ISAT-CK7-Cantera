!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module ceq_types   !  define data types
  implicit none

  integer, parameter :: k_sp = kind(1.e0), k_dp = kind(1.d0), k_wp = kind(1.d0)

!---------------------- sys_type ------------------------------------------
  type :: sys_type

	integer :: ne   ! number of elements
	integer :: ned  ! number of determined elements
	integer :: neu  ! number of undetermined elements
    integer :: ns   ! number of species
	integer :: nsd  ! number of determined species    
	integer :: nsu  ! number of undetermined species
	integer :: ncs  ! number of constrained species
	integer :: ng   ! number of general linear constraints
	integer :: nc   ! number of basic constraints =ne+ncs+ng
	integer :: nrc  ! number of reduced constraints
	integer :: nb   ! number of independent constraints =nsd+nrc

    integer, pointer :: CS(:)       ! indexes of constrained species (n_cs x 1)
	integer, pointer :: sp_order(:) ! the k-th ordered species is sp_order(k) (ns x 1)
	integer, pointer :: el_order(:) ! the j-th ordered element is el_order(j) (ne x 1)

	real(k_wp), pointer :: Ein(:,:) ! element matrix (ns x ne).  A molecule of species
!                                     k contains Ein(k,j) atoms of element j.
    real(k_wp), pointer :: E(:,:)   ! the re-ordered element matrix (ns x ne)
	real(k_wp), pointer :: Bg(:,:)  ! general linear constraint matrix (ns x ng)
	real(k_wp), pointer :: B(:,:)   ! basic constraint matrix (ns x nc)
	real(k_wp), pointer :: BR(:,:)  ! the reduced constraint matrix (nsu x nrc)
	real(k_wp), pointer :: A(:,:)   ! modified constraint transformation matrix (nb x nc)
	real(k_wp), pointer :: thermo(:,:)  ! thermodynamic coefficients (re-ordered) (ns x 15)

	logical :: initialized != .false. ! = .true. when sys has been initialized

! numerical parameters
! (when parameters are altered, ceq_param_set and ceq_param_def must be updated.)
    integer    :: diag     != 1       ! >0 for diagnostics
    integer    :: lu       != 0       ! logical unit number for diagnostics
	real(k_wp) :: T_low    != 250.d0  ! lowest allowed temperature
	real(k_wp) :: T_high   != 5000.d0 ! highest allowed temperature
	real(k_wp) :: frac_zm  != 1.d-1   ! fraction of zm used in initial guess
	real(k_wp) :: T_tol    != 1.d-6   ! convergence tolerance on log(T)
    real(k_wp) :: ds_inc   != 1.4d0   ! factor by which ds is increased after success (ceq_fixed_T)
    real(k_wp) :: ds_dec   != 0.25d0  ! factor by which ds is decreased after failure (ceq_fixed_T)
	real(k_wp) :: ds_min   != 1.d-15  ! smallest allowed time step
    real(k_wp) :: res_tol  != 1.d-9   ! convergence tolerance for residual (ceq_fixed_T)
    real(k_wp) :: ires_fac != 1.d2    ! factor by which the irreducible residual can exceed res_tol
	real(k_wp) :: logy_lim != 120.d0  ! upper limit on log(y)  (to prevent overflow in ceq_y)
	real(k_wp) :: srat_lim != 1.d-9   ! lower limit on singular-value ratio (ceq_Sinv)
	real(k_wp) :: err_huge != 1.d6    ! upper limit on error (ceq_newt)
	real(k_wp) :: dec_min  != 0.5d0   ! minimum acceptable decrease in residual (ceq_newt)
	real(k_wp) :: eps_el   != 1.d-9   ! relative lower bound on element moles (ceq_perturb)
	real(k_wp) :: eps_sp   != 1.d-9   ! relative lower bound on species moles (ceq_perturb)
	real(k_wp) :: pert_tol != 1.d-4   ! largest allowed normalized perturbation
	integer    :: pert_skip != 0       ! set =1 to skip perturbing undetermined species

  end type sys_type

!---------------------- state_type -----------------------------------------
  type :: state_type  ! properties in the constrained equilibrium state
    real(k_wp)          :: p       ! pressure
    real(k_wp)          :: T       ! temperature
    real(k_wp)          :: h       ! enthalpy
    real(k_wp), pointer :: zd(:)   ! moles of determined species (nsd)
    real(k_wp), pointer :: zu(:)   ! moles of undetermined species (nsu)
    real(k_wp), pointer :: lam(:)  ! Lagrange multipliers (nrc)
    real(k_wp), pointer :: cr(:)   ! reduced constraint vector
    real(k_wp), pointer :: Q(:)    ! consistency vector Q' * Xu = 1 (nsu)
    real(k_wp), pointer :: gu(:)   ! Gibbs functions of undetermined species (nsu)

    integer :: temp_its   != 0 ! number of temperature iterations performed
	integer :: time_steps != 0 ! number of pseudo time steps
    integer :: newt_calls != 0 ! number of Newton solves attempted
    integer :: newt_its   != 0 ! number of Newton iterations performed
    integer :: nyeval     != 0 ! number of y evaluations

  end type state_type

end module ceq_types
