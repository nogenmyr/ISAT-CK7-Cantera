!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ci_test_dpt

  use ci_dat, only: modeci
  use ci_dat8
  use ci_cem_recon
  use ci_cksubs
  use ci_utils
  use ci_ice_cksubs
  use ci_dpt_dr

  implicit none

  integer, parameter :: nr_max = 100000
  integer, parameter :: nr_sel = 200
  logical, parameter :: estimate_A = .true.
  integer :: i, lu_in, lu_out
  real(kind(1.d0)) :: x(nsp1), xdt(nsp1), y(ns), tmp
  real(kind(1.d0)) :: hs, h, dpt(3), dpt_fc(3), dpt_dr(3)
  real(kind(1.d0)) :: ZUN( nus, nr_sel ), YT( ne, ne ), XT(ne, nus)

  ! new arrays using ZZ^T
  integer ::  j, k
  real(kind(1.d0)) :: zui(nus), EDD(ne,nus), RHS(ne), HH(ne,ne)

  ! arrays for the LS soln for we via dgelss
  integer :: rank, info_ls, lwork, nm
  real(kind(1.d0)) :: ZUE(nr_sel, ne), SZU(nr_sel,1), SE(ne)
  real(kind(1.d0)), allocatable :: work_ls(:)
  
  ! matrices for SVD of ZUN
  real(kind(1.d0)) :: UU(nus,nus), SS(nus), VVT(1,1), ZS(ne,ne)
  real(kind(1.d0)) :: work(nus*nus+20*nus)

  integer :: info, info_svd, ipiv(ne+1)
  real(kind(1.d0)) :: cp, cpr, hr, h_fc, h_dr
  real(kind(1.d0)) :: rin(nrc), dpt_ce(3), xr(nrc+1), hs_n
  real(kind(1.d0)) :: z_CE(ns), x_CE(nsp1), T_CE, stats(20), hs_ce
  real(kind(1.d0)) :: err_rho_ce(nr_max), err_temp_ce(nr_max), err_rho_dr(nr_max), err_temp_dr(nr_max)
  real(kind(1.d0)) :: err_rho_a(nr_max), err_temp_a(nr_max)
  real(kind(1.d0)) :: err_rho_ce_stats(3), err_temp_ce_stats(3)
  real(kind(1.d0)) :: err_rho_dr_stats(3), err_temp_dr_stats(3)
  real(kind(1.d0)) :: err_rho_a_stats(3), err_temp_a_stats(3)

  call isat_lu( lu_in )
  call isat_lu( lu_out )
  open( lu_in,  file='rxn_map.op' )
  open( lu_out, file='dpta_errors.op', position='append')

  DD = 0.d0
  if(estimate_A) then
     ! compute the special matrix A
     do i = 1, nr_sel
        read(lu_in,*, end=100, err=100) k, x, xdt, h, tmp, dpt

        rin(1:nrc) = matmul(BBT, x(1:ns))
        call ci_ceq( 0, rin, h, prc, .false., z_CE, T_CE, &
             z_CE, T_CE, stats, info)

        ZUN(:,i) = z_CE(US)
        zui(1:nus) = z_CE(US)
        do j = 1, nus
           DD(:,j) = DD(:,j) + zui(j)*zui
        enddo
        ZUE(i,:) = rin(nrs+1:nrc)
        SZU(i,1) = sum(z_CE(US))
     enddo

     ! SVD of DD = UU* SS* VVT
     call dgesvd( 'S', 'N', nus, nus, DD, nus, SS, UU, nus, VVT, 1, &
          work, nus*nus+20*nus, info_svd)
     
     ! extract first ne columns of U
     XT = transpose(UU(:,1:ne))
     
     ! Form YT = XT*CEu
     YT = matmul(XT,CEu)
     
     ! Solve YT*PT = XT
     call dgesv( ne, nus, YT, ne, ipiv, XT, ne, info )
     
     ! update PT
     PT = XT
     
     ! set we
     call set_we
  endif

  err_rho_ce = 0.d0
  err_temp_ce = 0.d0
  err_rho_dr = 0.d0
  err_temp_dr = 0.d0

  rewind(lu_in)
  do i = 1, nr_max
     read(lu_in,*, end=100, err=100) k, x, xdt, h, tmp, dpt

     dpt_fc = 0.d0
     dpt_fc(2) = prc

     modeci = 7
     call ci_dens_temp( nsp1, x, dpt_fc )

     ! do reduction 
     rin(1:nrc) = matmul(BBT, x(1:ns))

     ! approximate dens, temp
     dpt_dr = 0.d0
     dpt_dr(2) = prc

     call h2hs_n( h, href_n, rin(1:nrc), nrc, hs_n )
     xr(1:nrc) = rin(1:nrc)
     xr(nrc+1) = hs_n

     modeci = 8
     call ci_dens_temp( nrc+1, xr, dpt_dr )

     ! compute dens temp using ceq reconstrcution
     call ci_ceq( 0, rin, h, prc, .false., z_CE, T_CE, &
          z_CE, T_CE, stats, info)

     call h2hs( h, href, z_CE, ns, hs_ce )
     x_CE(1:ns) = z_CE(1:ns)
     x_CE(nsp1) = hs_ce

     dpt_ce = 0.d0
     dpt_ce(2) = prc
     
     modeci = 7
     call ci_dens_temp( nsp1, x_CE, dpt_ce )

     ! compute errors
     err_rho_ce(i)  = abs((dpt_ce(1) - dpt_fc(1))/dpt_fc(1))*100
     err_temp_ce(i) = abs((dpt_ce(3) - dpt_fc(3))/dpt_fc(3))*100
     err_rho_dr(i)  = abs((dpt_dr(1) - dpt_fc(1))/dpt_fc(1))*100
     err_temp_dr(i) = abs((dpt_dr(3) - dpt_fc(3))/dpt_fc(3))*100
     err_rho_a(i)  = abs((dpt_dr(1) - dpt_ce(1))/dpt_ce(1))*100
     err_temp_a(i) = abs((dpt_dr(3) - dpt_ce(3))/dpt_ce(3))*100

     !write(lu_out,'(i8,1p,100e16.2)') i, dpt_fc(1), dpt_ce(1), dpt_dr(1), dpt_fc(3), dpt_ce(3), dpt_dr(3)
  enddo
  
100 continue

  ! compute error statistics
  err_rho_ce_stats(1)  = sum(err_rho_ce)/(i-1)   ! mean
  err_rho_ce_stats(2)  = sqrt(sum((err_rho_ce - err_rho_ce_stats(1))**2)/(i-1)) ! std
  err_rho_ce_stats(3)  = sqrt(dot_product(err_rho_ce,err_rho_ce)/nr_max)  ! rms error

  err_temp_ce_stats(1) = sum(err_temp_ce)/(i-1)  ! mean
  err_temp_ce_stats(2) = sqrt(sum((err_temp_ce - err_temp_ce_stats(1))**2)/(i-1)) ! std
  err_temp_ce_stats(3) = sqrt(dot_product(err_temp_ce,err_temp_ce)/nr_max) ! rms error

  err_rho_dr_stats(1)  = sum(err_rho_dr)/(i-1)   ! mean
  err_rho_dr_stats(2)  = sqrt(sum((err_rho_dr - err_rho_dr_stats(1))**2)/(i-1)) ! std
  err_rho_dr_stats(3)  = sqrt(dot_product(err_rho_dr,err_rho_dr)/nr_max) ! rms error

  err_temp_dr_stats(1) = sum(err_temp_dr)/(i-1)  ! mean
  err_temp_dr_stats(2) = sqrt(sum((err_temp_dr - err_temp_dr_stats(1))**2)/(i-1)) ! std
  err_temp_dr_stats(3) = sqrt(dot_product(err_temp_dr,err_temp_dr)/nr_max) ! rms error

  err_rho_a_stats(1)  = sum(err_rho_a)/(i-1)   ! mean
  err_rho_a_stats(2)  = sqrt(sum((err_rho_a - err_rho_a_stats(1))**2)/(i-1)) ! std
  err_rho_a_stats(3)  = sqrt(dot_product(err_rho_a,err_rho_a)/nr_max) ! rms error

  err_temp_a_stats(1) = sum(err_temp_a)/(i-1)  ! mean
  err_temp_a_stats(2) = sqrt(sum((err_temp_a - err_temp_a_stats(1))**2)/(i-1)) ! std
  err_temp_a_stats(3) = sqrt(dot_product(err_temp_a,err_temp_a)/nr_max) ! rms error

  write(lu_out,'(i3,i8,1p,100e16.8)') nrs, i-1, err_rho_ce_stats, &
       err_rho_dr_stats, err_rho_a_stats, err_temp_ce_stats, &
       err_temp_dr_stats, err_temp_a_stats

end subroutine ci_test_dpt
