!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

!  This file contains all the Chemistry Interface routines for the Dimension Reduction (DR) modes, namely:

!  cinit_dr : ciinit routine for the DR modes
 
!  cirxn_dr( t, c0, ct, dpt ) : cirxn routine for the DR modes
   
!  cicmp_dr( c, krep, comp )  : cicomp routine for the DR modes

!  cicomp_sr( ncv, c, krep, nfullv, comp, cname ) : species reconstruction

!  cidpt_dr( cc, dpt ) : routine to return approximated T and rho

!  cidpt_dr_update(etol, success, stats_dpt) : routine to reset and
!                 update the temp. and density values

!=====================================================================

subroutine cinit_dr
  ! Initialization routine for modes 8 and 9
  use ci_8
  call cinit8

end subroutine cinit_dr

!=====================================================================

subroutine cirxn_dr( t, c0, ct, dpt )

  ! Interface routine to cirxn8
  use ci_dat
  use ci_8
  implicit none
  
  real(k_dp), intent(in)    :: t, c0(nc)
  real(k_dp), intent(inout) :: dpt(3)
  real(k_dp), intent(out)   :: ct(nc)
  
  call cirxn8( t, c0, ct, dpt )
end subroutine cirxn_dr

!=====================================================================
subroutine cicmp_dr( c, krep, comp )

  ! Interface routine to cicmp8
  use ci_dat
  use ci_8
  implicit none
  
  real(k_dp), intent(in)    :: c(nc)
  real(k_dp), intent(inout) :: comp(nfull)
  integer, intent(in)       :: krep
 
  call cicmp8( c, krep, comp )
end subroutine cicmp_dr

!=====================================================================

subroutine ciz2ze_dr( ncv, z, nev, ze )
  ! rotine to determine specific moles of the elements from specific moles of species

! input:
!  ncv    - number of reduced composition variables (integer)
!  z      - composition vector in specific moles (length ncv, double)
!  nev    - number of elements
      
! output:
!  ze     - specific moles of the elements (ename)
  
  use ci_dat
  use ci_dat8
  implicit none

  integer, intent(in)        :: ncv, nev
  real(k_dp), intent(in)     :: z(ncv)
  real(k_dp), intent(out)    :: ze(nev)

  real(k_dp) :: zer(nev)

  zer = matmul(transpose(CEr), z(1:nrs))
  ze = zer + z(nrs+1:nrs+nev)

  return
end subroutine ciz2ze_dr

!=====================================================================

subroutine cicomp_sr( ncv, c, krep, nfullv, comp, cname )

  !   Chemistry interface routine to return full composition, comp,
  !   corresponding to reduced composition, c, for modes 8 and 9 by
  !   performing a species reconstruction (sr) c = {r, hs_n}; comp =
  !   {x/y/z, dens, temp, press, h}

  ! input:
  !      ncv     - number of reduced composition variables (integer)
  !              - If ncv=0, then only cname is returned: other input ignored.
  !      c       - composition vector (length ncv, double)
  !      krep    - type of representation required
  !          = 1 - express species as mole fractions
  !          = 2 - express species as mass fractions
  !          = 3 - express species as specific mole numbers
  !      nfullv  = number of full composition variables ( = ns + 4 )
  !      comp(1) - pressure in Chemkin units	

  ! output:
  !      comp    - full composition vector (length ns+4, double)
  !	  comp(i)    = species(i), i=1, ns
  !         comp(ns+1) = density (Chemkin units)
  !         comp(ns+1) = temperature (K)
  !         comp(ns+2) = pressure (Chemkin units)
  !         comp(ns+3) = enthalpy (Chemkin units)
  !      cname   - names of composition variables (nfullv, character*(*))

  use ci_2
  use ci_6
  use ci_dat8
  use ci_cksubs

  implicit none
  external ci_recon

  integer, intent(in)        :: ncv, krep, nfullv
  real(k_dp), intent(in)     :: c(ncv)
  real(k_dp), intent(inout)  :: comp(nfullv)
  character*(*), intent(out) :: cname(nfullv)

  ! local variables
  real(k_dp) :: press, hs_n, h, z_DR(ns), T_DR
  real(k_dp) :: y(ns), dpt(3)

  cname = cmpsym_f

  if( ncv == 0 ) return

  ! check if nfullv == ns + 4
  if( nfullv /= ns + 4 ) then
     call isat_abort( 'cicomp_sr', 1, mess = 'nfullv /= ns + 4 = ', isv=ns+4 )
  endif

  press      = comp(1)
  dpt(2)     = press
  hs_n       = c(ncv)

  if( modeci == 8  .or.  modeci == 9 ) then 
     call hs2h_n( hs_n, href_n, c(1:nrc), nrc, h ) ! h = enthalpy
     call ci_recon( c(1:nrc), press, h, z_DR, T_DR )
  else
     call isat_abort( 'cicomp_sr', 1, mess = 'bad modeci = ', isv=modeci )
  endif

  ! convert to right units

  if( krep == 1 ) then
     comp(1:ns) = z_DR(1:ns) / sum( z_DR(1:ns) )
  elseif( krep == 2 ) then
     comp(1:ns) = y(1:ns)
  elseif( krep == 3 ) then
     comp(1:ns) = z_DR(1:ns)
  endif

  ! set density and temperature returned from cidpt_dr
  call cidpt_dr( c, dpt )

  comp(ns+1) = dpt(1)
  comp(ns+2) = dpt(3)
  comp(ns+3) = dpt(2)
  comp(ns+4) = h

  return
end subroutine cicomp_sr

!=====================================================================
subroutine cidpt_dr( cc, dpt )
  ! routine to return approximated density and temperature

  use ci_dat8
  use ci_ice_cksubs
  use ci_stats
  use ci_dpt_dr
  implicit none
  
  real(k_dp), intent(in)    :: cc(nc)
  real(k_dp), intent(inout) :: dpt(3)

  !local 
  real(k_dp) :: temp, hs_n, h, sz

  call routine_start(i_cidpt_dr)

  ! compute approximated temperature
  hs_n = cc(nc)
  call hs2h_n( hs_n, href_n, cc(1:nrc), nrc, h )
  call temphr( h, cc(1:nrc), temp )
  dpt(3) = temp

  ! compute approximated density
  sz = sum(cc(1:nrs)) + dot_product(we, cc(nrs+1:nrc))
  dpt(1) = dpt(2)/(sz*gascon*temp)

  call routine_stop(i_cidpt_dr)

end subroutine cidpt_dr

!=========================================================================
subroutine cidpt_dr_update(etol, success, stats_dpt)
  ! routine to reset and update dens. and temp. values

  ! input 
  !  etol : error tolerance
  
  ! output:
  ! success: logical
  !   true:  the dens. temp. values will be updated
  !   false: the dens. temp. values will not be updated
  ! stats: statistics
  !   stats(1) = current mean relative error in approximate temp.
  !   stats(2) = mean relative error after temp. values are updated
  !   stats(3) = no. of leaves in the ISAT table
  !   stats(4) = no. of ISAT special calls made for computing the error

  ! Note: In the current implementation:
  ! 1/ If current error < etol, nothing is done. success = .true.
  ! 2/ If current error > etol and new error < etol, then all the
  ! exisiting leaves in the ISAT table are deleted, and dens. and
  ! temp. values are updated in the subsequent calls to ci_dens_temp.

  use ci_dat 
  use ci_dat6
  use ci_dpt_dr
  implicit none
  
  real(k_dp), intent(in)  :: etol
  logical, intent(out)    :: success
  real(k_dp), intent(out) :: stats_dpt(10)

  ! local vairables
  integer    :: mode, info(100), iusr(1), dpt_info
  real(k_dp) :: x(1), rusr(1), f(1), dfdx(1,1), hvar(1), rinfo(50), stats(100)
  integer, save :: deferred = 0

  if(modeci /= 8 .and. modeci /= 9) return

  ! log message
  if(dtlog==1) then
     if( n_recon < 50 ) then
        write(dta_lu,*) 'cidpt_dr_update: stored zu values (less than 50) = ', n_recon
        deferred = deferred + 1
        if( deferred > 20 ) then
           write(dta_lu,*) 'cidpt_dr_update: deferred 20 times, returning with success = .true.'
           success = .true.
        else
           write(dta_lu,*) 'cidpt_dr_update: deferring the update, count = ', deferred
        endif
        return
     endif
     write(dta_lu,*) 'cidpt_dr_update: attempting to update using ', n_recon, ' stored zu values'
  endif

  ! reset dpt errors to zero
  call reset_dpt_errors
  
  ! initialize values for isat special call
  mode  = 24
  x     = 0.d0
  iusr  = 0
  rusr  = 0.d0
  info  = 0
  rinfo = 0
  stats = 0.d0
  
  info(80) = -12345 ! call_flag
  info(81) = 100    ! max_call
  info(82) = 2      ! leaf_stride

  ! compute nPT
  call update_nPT(dpt_info)
  if (dpt_info /= 0) then
     write(dta_lu,*)'ci_dpt_dr: updating PT failed, rank, dpt_info = ', myrank, dpt_info
     return
  endif

  call isatab( 1, mode, 1, x, 1, 0, 1, temp_error_leaves, &
       iusr, rusr, info, rinfo, f, dfdx, hvar, stats )

  cur_error = cur_error/stats(2) ! mean error
  new_error = new_error/stats(2) ! mean new error

  stats_dpt = 0.d0
  stats_dpt(1) = cur_error
  stats_dpt(2) = new_error
  stats_dpt(3) = stats(1)
  stats_dpt(4) = stats(2)

  if(dtlog==1) then
     write(dta_lu,*) 'cidpt_dr_update: current error =', cur_error
     write(dta_lu,*) 'cidpt_dr_update: error achieved (after resetting) =', new_error
  endif

  if( cur_error <= etol ) then
     if(dtlog==1) write(dta_lu,*) 'cidpt_dr_update: doing nothing, current error <= etol'
     if(dtlog==1) write(dta_lu,*) 'cidpt_dr_update: done!'
     success = .true.

  elseif( new_error <= etol ) then
     if(dtlog==1) write(dta_lu,*) 'cidpt_dr_update: resetting dta matrix, new error <= etol'
     call reset_dpt_dr

     ! delete all the leaves
     mode  = 21
     x     = 0.d0
     iusr  = 0
     rusr  = 0.d0
     info  = 0
     rinfo = 0
     stats = 0.d0
     call isatab( 1, mode, 1, x, 1, 0, 1, temp_error_leaves, &
       iusr, rusr, info, rinfo, f, dfdx, hvar, stats )

     success = .true.
  else
     success = .false.
     if(dtlog==1) write(dta_lu,*) 'cidpt_dr_update: updating matrix failed, error > etol'
  endif

  return
end subroutine cidpt_dr_update
!=========================================================================
