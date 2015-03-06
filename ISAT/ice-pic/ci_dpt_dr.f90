!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module ci_dpt_dr

! TEMPERATURE AND DENSITY COMPUTATION

! Module for maintaining data related to approximation of temperature
! and density in the DR modes.

  use ci_prec
  use ci_dat
  use ci_dat6
  use ci_dat8
  use ci_utils
  implicit none

  integer, save    :: n_recon    ! no. of reconstruction points stored in DD
  logical, save    :: updated    ! true if PT successfully updated
  integer, save    :: ludpt      ! logical unit for reading PT matrix
  character(30)    :: fname      ! file to store the PT matrix
  character(30)    :: blank, head, tail

  ! NOTE: The file "dta.in" is used to store the special matrix PT. On
  ! rerun if the file "dta.in" exists, then the PT matrix will be read
  ! from this file for initialization.

  real(k_dp), allocatable :: DD(:,:), PT(:,:), cPT(:,:), nPT(:,:), we(:)
  ! DD(nus, nus)   special matrix to store the zu values
  ! PT(1:ne,1:nus) special matrix to compute enthalpy for the reduced modes
  !                cPT -> current PT
  !                nPT -> new PT
  ! we(1:ne)       special vector to approximate sum(zu) to compute density

  real(k_dp), save :: cur_error, new_error

contains

!=========================================================================

subroutine ci_dpt_alloc

  use ci_dat8

  allocate( we(ne) )
  allocate( PT(ne, nus) )
  allocate( cPT(ne, nus) )
  allocate( nPT(ne, nus) )
  allocate( DD(nus, nus) )

  n_recon = 0
  DD = 0.d0
  updated = .false.

  ! dta_#.in input file
  blank = repeat(' ',30)
  head  = blank
  head  = 'dta'
  tail  = blank   
  tail  = 'in'
  call isat_file_name( head, -1, myrank, tail, fname )
  call isat_lu(ludpt)

end subroutine ci_dpt_alloc

!=========================================================================

subroutine set_we

  implicit none
  integer :: i

  ! reset the special vector we
  do i=1,ne
     we(i) = sum(PT(i,:))
  enddo

end subroutine set_we

!=========================================================================
subroutine ci_dpt_init

  ! Initialize PT and we arrays
  ! PT = pseduo-inverse(CEu^T)
  ! we = sum(PT)

  use ci_dat8
  implicit none

  ! SVD computation array
  logical :: exists = .false.
  integer :: i, info_svd
  integer :: nspin, nein, nrsin, rCS(nrs)
  real(k_dp), allocatable :: CEuT(:,:), XU(:,:), XS(:), XVT(:,:)
  real(k_dp), allocatable :: XSinvT(:,:), EUT(:,:), work(:)

  call ci_dpt_alloc
  
  if(dtlog==1) write(dta_lu,*) 'ci_dpt_dr: initializing dta matrix'

  ! check if the "dta_#.in" file exists
  inquire( file = fname, exist = exists )
  
  if( exists ) then
     ! read PT matrix for the file "dta.in"
     open( ludpt, file = fname )
     read( ludpt, * ) nspin, nein, nrsin
     if ( nspin == ns .and. nein == ne .and. nrsin == nrs ) then
        read( ludpt, * ) rCS
        if( maxval(abs(rCS - CS)) == 0 ) then
           read( ludpt, * ) PT
           if(dtlog==1) write(dta_lu,*) 'ci_dpt_dr: info read from', fname, 'successfully.'
        else
           if(dtlog==1) then
              write(dta_lu,*) 'ci_dpt_dr: reading ', fname, ' failed.'
              write(dta_lu,'(a,3i3)') ' The represented species do not match with the saved species.'
           endif
        endif
     else
        if( dtlog==1 ) then
           write(dta_lu,*) 'ci_dpt_dr: reading ', fname, ' failed.'
           write(dta_lu,'(a,3i3)') ' File corrupt and/or incompatible: ns, ne, nrs /= ', nspin, nein, nrsin
        endif
        exists = .false.
     endif
     close( ludpt )
  endif
  if( .not.exists ) then
     ! compute pseudo-inverse transpose of CEuT
     allocate(CEuT(ne,nus))
     allocate(XU(ne,ne))
     allocate(XS(ne))
     allocate(XVT(nus,nus))
     allocate(XSinvT(ne,nus))
     allocate(work(nus*ne+10*(nus+ne)))
     
     CEuT = transpose(CEu)
     
     ! compute SVD of CEuT = XU * XS * XVT
     call dgesvd('A','A', ne, nus, CEuT, ne, XS, XU, ne, &
          XVT, nus, work, nus*ne+10*(nus+ne), info_svd )
     
     ! compute the pseudo-inverse transpose of XS
     XSinvT = 0.d0
     do i=1,ne
        XSinvT(i,i) = 1.d0/max(XS(i),XS(1)*1.d-10)
     enddo
     
     ! PT = pseudo-inverse transpose of CEuT
     PT = matmul(matmul(XU,XSinvT),XVT)
  endif

  cPT = PT
  nPT = PT

  ! set the special vector we
  call set_we

  if(allocated(CEuT))       deallocate(CEuT)
  if(allocated(XU))         deallocate(XU)
  if(allocated(XS))         deallocate(XS)
  if(allocated(XVT))        deallocate(XVT)
  if(allocated(XSinvT))     deallocate(XSinvT)
  if(allocated(work))       deallocate(work)

  return
end subroutine ci_dpt_init

!=========================================================================

subroutine update_nPT( info )

  ! Special routine to compute the new PT matrix using DD

  ! Output: info
  ! =  0: successful
  ! = -1: dgesvd failed
  ! = -2: dgesv failed

  use ci_dat8
  implicit none

  integer, intent(out) :: info

  ! Local variables
  integer :: info_svd, info_sv, ipiv(ne+1)
  real(kind(1.d0)) :: UU(nus,nus), SS(nus), VVT(1,1)
  real(kind(1.d0)) :: work(nus*nus+20*nus), DD_copy(nus, nus)
  real(kind(1.d0)) :: YT( ne, ne ), XT(ne, nus)

  DD_copy = DD
  ! SVD of DD = UU* SS* VVT
  call dgesvd( 'S', 'N', nus, nus, DD_copy, nus, SS, UU, nus, VVT, 1, &
       work, nus*nus+20*nus, info_svd)
  
  if( info_svd /= 0) then
     info = -1
     return
  endif

  ! Extract first ne columns of U
  XT = transpose(UU(:,1:ne))
  
  ! Form YT = XT*CEu
  YT = matmul(XT,CEu)
  
  ! Solve YT*PT = XT
  call dgesv( ne, nus, YT, ne, ipiv, XT, ne, info_sv )
  
  if (info_sv /= 0) then
     info = -2
     return
  endif

  ! set nPT (new PT)
  nPT = XT
  info = 0

  return

end subroutine update_nPT

!=========================================================================

subroutine reset_dpt_dr

  ! subroutine to reset PT, cPT to nPT
  PT  = nPT
  cPT = nPT
  call set_we

  ! dump the matrix in ludpt
  open( ludpt, file=fname )
  rewind( ludpt )
  write(ludpt, *) ns, ne, nrs
  write(ludpt, *) CS(1:nrs)
  write(ludpt, *) PT

  ! print information on screen
  write(0,*) 'ci_dpt_dr: dta.in file has been created on rank', myrank

end subroutine reset_dpt_dr

!=========================================================================

subroutine use_nPT

  ! use nPT for approximation
  PT  = nPT
  call set_we

end subroutine use_nPT

!=========================================================================

subroutine use_cPT

  ! use cPT for approximation
  PT  = cPT
  call set_we

end subroutine use_cPT

!=========================================================================

subroutine store_zu(zu)

  ! store zu values in the special matrix DD

  implicit none
  real(k_dp), intent(in) :: zu(nus)
  real(k_dp) :: stats(10)
  integer :: i

  do i = 1, nus
     DD(:,i) = DD(:,i) + zu(i)*zu
  enddo
  n_recon = n_recon + 1

  ! check if PT needs to be updated
  if( mod(n_recon, ufreq) == 0 .and. autou == 1 .and. .not.updated ) then
     call cidpt_dr_update( tatol, updated, stats )
  endif
end subroutine store_zu

!=========================================================================

subroutine ha_dr( temp, r, ha )

  ! routine to return approximated enthalpy
  implicit none

  real(k_dp), intent(in)  :: temp, r(nrc)
  real(k_dp), intent(out) :: ha

  real(k_dp) :: hi(ns)

  ! enthalpies per unit mole
  call cthml( temp,   hi, gas )

  ha = dot_product( r(1:nrs), hi(CS) ) +  &
       dot_product( matmul(PT, hi(US)), r(nrs+1:nrc) )

  return
end subroutine ha_dr

!=========================================================================

subroutine cpa_dr( temp, r, cpa )

  ! routine to return approximated specific heat
  implicit none

  real(k_dp), intent(in)  :: temp, r(nrc)
  real(k_dp), intent(out) :: cpa

  real(k_dp) :: cpi(ns)

  ! specific heats per unit mole
  call ctcpml( temp,   cpi, gas )

  cpa = dot_product( r(1:nrs), cpi(CS) ) +  &
       dot_product( matmul(PT, cpi(US)), r(nrs+1:nrc) )

  return
end subroutine cpa_dr

!=========================================================================

subroutine reset_dpt_errors

  cur_error = 0.d0
  new_error = 0.d0

end subroutine reset_dpt_errors

!=========================================================================
subroutine temp_error_leaves( replace, nx, x, nf, nh, i_wrk, r_wrk, fa, ga, ha )

  use ci_prec
  use ci_dat8
  use ci_utils

  integer,    intent(in)  :: nx, nf, nh, i_wrk(3)
  integer,    intent(out) :: replace(3)
  
  real(k_dp) :: x(nx), r_wrk(nx+nf+11), fa(nf), ga(nf,nx), ha(nh)
  real(k_dp) :: xscale(nx), fscale(nf), etolsq, prop(10), xx(nx)
  real(k_dp) :: ff(nf), ct(nrc+1), dpt(3), T_DR, Ta_c, Ta_n, press

  integer    :: nclipt, nclipt_log

  replace = 0
  xscale  = r_wrk(1:nx)
  fscale  = r_wrk(nx+1:nx+nf)
  etolsq  = r_wrk(nx+nf+1)
  prop    = r_wrk(nx+nf+2:nx+nf+12)

  xx = x*xscale    ! unscale x
  ff = fa*fscale   ! unscale f
  ct = ff(1:nrc+1) ! {r(dt), hs_n(dt)}
  Ta_c = ff(nf)
  T_DR = ha(1)

  ! set pressure
  if( mode_pdt == 1 .or. mode_pdt == 2 ) then
     press = xx(nrc+2)
  else
     press = prc
  endif

  ! update cur_error
  cur_error = cur_error + abs(Ta_c - T_DR)/T_DR

  ! use nPT to recompute Ta
  call use_nPT
  ! since new PT matrix may be ill-conditioned, clip temp.
  nclipt = clipt
  nclipt_log = clipt_log
  clipt = 1       ! clip the temperature
  clipt_log = 0   ! disable temperature clipping log

  dpt(2) = press
  call cidpt_dr( ct, dpt )
  Ta_n = dpt(3)

  ! reset clipt values
  clipt = nclipt           
  clipt_log = nclipt_log   

  new_error = new_error + abs(Ta_n - T_DR)/T_DR

  ! reset back PT to cPT
  call use_cPT

  return
end subroutine temp_error_leaves
!=========================================================================
end module ci_dpt_dr
