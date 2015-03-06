!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ci_recon( rin, press, hin, z_DR, T_DR )
  ! routine to perform species reconstruction for modes 8 and 9
  
  use ci_8
  use ci_dat8
  use ci_cem_recon
  implicit none

  real(k_dp), intent(in) :: rin(nrc), press, hin
  real(k_dp), intent(out):: z_DR(ns), T_DR

  ! local variables
  integer    :: info, lu, k_facet, n_iters, n_integs, method
  real(k_dp) :: r_e(nrc), r_g(nrc), z_e_CE(ns), z_g(ns), tau_ICE, stats(20)
  real(k_dp) :: A_ICE_ODE(ns+1, ns+1), r_normal(nrc), T_e_ICE, T_g

  if( modeci == 8 ) then
     call ci_ceq( 0, rin, hin, press, .false., z_DR, T_DR, &
          z_DR, T_DR,  stats, info, lu )
  elseif( modeci == 9 ) then
     call ci_ice_recon( rin, hin, press, z_e_CE, r_g, z_g, T_g, tau_ICE,  &
          r_e, z_DR, T_DR, A_ICE_ODE, r_normal, k_facet,  &
          n_iters, n_integs, method, info )
  else
     call isat_abort( 'ci_recon', 1, mess = 'bad modeci = ', isv=modeci )
  endif
  return

end subroutine ci_recon
