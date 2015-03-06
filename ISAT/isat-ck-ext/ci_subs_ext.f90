!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

!  This file contains the extended Chemistry Interface routines, namely:

!   cistrm_name( cstrm, ncv, c, dpt )                            scistrm_name

!   cisize_full_ext( nfull_ext )                                 cisize_full_ext

!   cicomp_ext( ncv, c, krep, nfull_ext, comp_ext, cname_ext )   scicomp_ext

!   ci_mix_frac( mode, es_names, nc, z, workspace, mf )          sci_mix_frac

!   ci_dens_temp( ncv, cc, dpt )                                 sci_dens_temp

!   ci_dens_temp_update( tatol, success, stats )                 sci_dens_temp_update

!==================================================================================

subroutine cistrm_name(cstrm, ncv, c, dpt)
!
!  chemistry interface routine to return composition of specified stream.
!
!  input:
!      cstrm  - name of the stream (character)
!      ncv    - number of composition variables (integer)
!  output:
!      c      - composition vector for stream istrm (length ncv, double)
!      dpt    - density, pressure and temperature of stream istrm (length 3, double)

  use ci_2
  use ci_prec
  use ci_dat
  use streams_mod
  implicit none
  
  character*(*),intent(in):: cstrm
  integer, intent(in)     :: ncv
  real(k_dp), intent(out) :: c(ncv), dpt(3)
  integer                 :: istrm
  integer                 :: cistrmIdx
  
  istrm = cistrmIdx(cstrm)
  
  if(modeci == 2 .or. modeci == 3) then
     c(1) = StrmDat%Cvalue(istrm,1)
     call ci_dpt2( c(1), dpt )
  else
     call cistrm( istrm, ncv, c, dpt )
  end if
  
end subroutine cistrm_name

!=====================================================================

subroutine cisize_full_ext( nfull_ext )

  !     routine to retrieve size of the full extended composition

  !     input:      - none
  !     output:
  !      nfull_ext - size of the full extended composition

  use ci_dat
  use ci_dat6
  implicit none

  integer, intent(out) :: nfull_ext

  if ( modeci == 8 .or. modeci == 9 ) then
     nfull_ext = ns + 4
  else
     nfull_ext = nfull
  endif

  return
end subroutine cisize_full_ext

!=====================================================================

subroutine cicomp_ext( ncv, c, krep, nfull_ext, comp_ext, cname_ext )

  ! Chemistry interface routine to return the full extended
  ! composition, comp_ext, given the reduced composition, c.

  ! input:
  !      ncv         - number of reduced composition variables (integer)
  !                  - If ncv=0, then only cname is returned: other input ignored.
  !      c           - composition vector (length ncv, double)
  !      krep        - type of representation required
  !          = 1     - express species as mole fractions
  !          = 2     - express species as mass fractions
  !          = 3     - express species as specific mole numbers
  !      nfull_ext   - number of full extended composition variables (value can be
  !                    obtained from the cisize_full_ext routine)
  !      comp_ext(1) - pressure in Chemkin units

  ! output:
  !      comp_ext    - full extended composition vector (length nfull_ext, double)
  !      cname_ext   - names of composition variables (nfull_ext, character*(*))

  ! Notes:
  !  1/ For modes 1, 2, 3, 6 and 7 comp_ext = comp
  !  2/ For modes 8 and 9, comp_ext is computed by performing species
  ! reconstruction. Given c = {r, hs_n}; comp_ext = {z, dens, T, p, h}

  use ci_prec
  use ci_dat
  implicit none
  external cicomp_sr

  integer, intent(in)        :: ncv, krep, nfull_ext
  real(k_dp), intent(in)     :: c(ncv)
  real(k_dp), intent(inout)  :: comp_ext(nfull_ext)
  character*(*), intent(out) :: cname_ext(nfull_ext)

  if( modeci == 8 .or. modeci == 9 ) then
     call cicomp_sr( ncv, c, krep, nfull_ext, comp_ext, cname_ext )
  else
     call cicomp( ncv, c, krep, nfull_ext, comp_ext, cname_ext )
  end if

end subroutine cicomp_ext

!=====================================================================

subroutine ci_mix_frac( mode, es_names, ncv, z, workspace, mf )
  ! Routine to compute mixture fraction (mf)

  ! Input:
  !  mode = -1 - initialization
  !       =  1 - computation
  !  es_names = names of the elements and streams considered for
  !             computing the mixture fraction, character array (size = 3)
  !             1/ First character is names of the elements
  !             2/ Second is the name of the base stream (str1)
  !             3/ Third is the name of the reference stream (str2)
  !             [ mf  = (w^T*(ze - ze1))/(w^T*(ze2 - ze1)) ]
  !  ncv - number of components
  !    z - compact composition (in species specific moles)

  ! InOut:
  !  workspace : real, size = 2*ncv
  !   Input (mode = -1) :
  !       nec   - number of elements considered
  !       wc    - weights for the mixture fraction, size(nec)
  !   Output :
  !       wfull - size(ne)
  !       ze1   - ze for the base stream - size(ne)
  !       beta  - 1/wfull^T*(ze2 - ze1) (constant factor)

  ! Output:
  !   mf : mixture fraction

  use ci_dat
  use ci_dat6
  use streams_mod
  implicit none

  integer, intent(in) :: mode, ncv
  real(k_dp), intent(in)  :: z(ncv)
  real(k_dp), intent(inout) :: workspace(2*ncv)
  real(k_dp), intent(out) :: mf
  character*(*), intent(in) :: es_names(3)

  ! local variables
  integer :: i, index, nec, str1, str2
  logical :: error
  real(k_dp) :: beta, wfull(ne), ze(ne), ze1(ne), ze2(ne)
  real(k_dp) :: wc(int(workspace(1)))
  character(16) :: elems

  if(mode == -1) then
     ! initialization
     nec = workspace(1)
     wc = workspace(2:1+nec)
     elems = es_names(1)
     wfull = 0.d0

     do i = 1, nec
        call ctcomp(elems(i:i), index, gas)
        wfull(index) = wc(i)
     enddo

     ! get streams info
     str1 = streams_index(es_names(2), error)
     str2 = streams_index(es_names(3), error)

     call ciz2ze( nc, ccstrm(:,str1), ne, ze1 )
     call ciz2ze( nc, ccstrm(:,str2), ne, ze2 )

     beta = 1.d0/dot_product(wfull, ze2 - ze1)

     ! save the values in workspace
     workspace = 0.d0
     workspace(1:ne) = wfull(1:ne)
     workspace(ne+1:2*ne) = ze1(1:ne)
     workspace(2*ne+1) = beta

     call ciz2ze( nc, z, ne, ze )
     mf = dot_product(wfull, ze - ze1)*beta
     return
  else if(mode == 1) then
     ! copy values from workspace

     wfull(1:ne) = workspace(1:ne)
     ze1(1:ne) = workspace(ne+1:2*ne)
     beta = workspace(2*ne+1)

     call ciz2ze( nc, z, ne, ze )
     mf = dot_product(wfull, ze - ze1)*beta
     return

  else
     write(0,*) ' ci_mix_frac: unrecognized mode = ', mode
  endif

end subroutine ci_mix_frac

!=====================================================================
subroutine ci_dens_temp( ncv, cc, dpt )
  ! Routine to compute density and temperature for a given compact
  ! composition. In the Dimension Reduction modes 8 and 9, this
  ! routine returns the approximated density and temperature.

  ! input:
  !  ncv - number of components
  !    c - compact composition
  !        for modes 6 to 9, c = z (specific moles)
  ! dpt(2) = pressure
  
  ! output:
  ! dpt(1) = density
  ! dpt(3) = temperature

  use ci_dat
  use ci_prec
  use ci_stats

  implicit none
  external cidpt_dr

  integer, intent(in)       :: ncv
  real(k_dp), intent(in)    :: cc(ncv)
  real(k_dp), intent(inout) :: dpt(3)

  call routine_start(i_ci_dens_temp)

  if ( modeci == 1 ) then
     call ci_dpt1( ncv, cc, dpt )
  elseif ( modeci == 2 .or. modeci == 3 ) then
     call ci_dpt23( ncv, cc, dpt )
  elseif ( modeci == 6 .or. modeci == 7 ) then
     call ci_dpt67( ncv, cc, dpt )
  elseif ( modeci == 8 .or. modeci == 9 ) then
     call cidpt_dr( cc, dpt )
  endif
  
  call routine_stop(i_ci_dens_temp)

end subroutine ci_dens_temp

!=====================================================================
subroutine ci_dens_temp_update( tatol, success, stats )

  ! This routine can be used to improve the temperature and density
  ! approximations returned by the ci_dens_temp routine. For a
  ! specified value of error tolerance, tatol, if the mean relative
  ! error in temperature can be made less than the specified error
  ! tolerance value, then the dens. & temp. values are reset and
  ! updated in the subsequent calls of ci_dens_temp routine.

  ! Instead of calling this routine, one can also set the option to
  ! auto-update the density and temperature values by setting the
  ! following parameters in the ci_ext.nml namelist file:

  ! autou      binary variable for auto-updating dens. and temp values.
  !            To enable auto-update, set autou = 1
  !             [default value, autou = 0 (do not auto-update)]
  ! ufreq      frequency with which to attempt updating the temp. values
  !             [default value, ufreq =  50 ]
  ! tatol      mean relative error tolerance in the approx. temp. Ta
  !             [default value, tatol = 1.d-2 (1% error)]

  ! If this routine is invoked, the specified input value of tatol
  ! overrides the value specified in the ci_ext.nml namelist file.

  ! Note: In the current implementation:
  ! 1/ If current error < etol, nothing is done. success = .true.
  ! 2/ If current error > etol and new error < etol, then all the
  ! existing ISAT leaves are deleted, and the dens. and temp. values
  ! are updated in the subsequent calls of ci_dens_temp. So this
  ! routine should be used prudently, as it deletes the previous
  ! additions from the ISAT table.

  ! input: 
  !  tatol = error tolerance ( suggested value, tatol = 1.d-2 )
  !          if the mean relative error in Ta < tatol, update dens. &
  !          temp. values
  
  ! output:
  !  success: logical
  !   true  = the dens. & temp. values will be updated in the
  !           subsequent calls to ci_dens_temp
  !   false = the dens. & temp. values will not be updated, 
  !           (current error > tatol)

  !  stats: statistics
  !   stats(1) = current mean relative error in approximate temp.
  !   stats(2) = mean relative error after temp. values are updated
  !   stats(3) = no. of leaves in the ISAT table
  !   stats(4) = no. of ISAT special calls made for computing the error

  use ci_dat
  use ci_prec
  use ci_stats
  implicit none
  external cidpt_dr_update

  real(k_dp), intent(in)  :: tatol
  logical, intent(out)    :: success
  real(k_dp), intent(out) :: stats(10)

  if(modeci /= 8 .and. modeci /= 9) then
     success = .true.
     return
  else
     call cidpt_dr_update( tatol, success, stats )
  endif
  
end subroutine ci_dens_temp_update
!=====================================================================
subroutine scistrm_name(cstrm, ncv, sc, sdpt)
    use ci_prec
    implicit none

    character*(*),intent(in):: cstrm
    integer, intent(in)     :: ncv
    real(k_sp), intent(out) :: sc(ncv), sdpt(3)

    ! local variables
    real(k_dp)              :: c(ncv), dpt(3)
    
    ! call cistrm_name
    call cistrm_name(cstrm, ncv, c, dpt)
    
    ! convert to single precision
    sc   = c
    sdpt = dpt
    
end subroutine scistrm_name

!=====================================================================

subroutine scicomp_ext( ncv, sc, krep, nfull_ext, scomp_ext, cname_ext )

  use ci_prec
  implicit none

  integer, intent(in)        :: ncv, krep, nfull_ext
  real(k_sp), intent(in)     :: sc(ncv)
  real(k_sp), intent(inout)  :: scomp_ext(nfull_ext)
  character*(*), intent(out) :: cname_ext(nfull_ext)

  ! local variables
  real(k_dp) :: c(ncv), comp_ext(nfull_ext)

  ! convert to double precision
  c = sc
  comp_ext(1) = scomp_ext(1)

  ! call cicomp_ext
  call cicomp_ext( ncv, c, krep, nfull_ext, comp_ext, cname_ext )

  ! convert back
  scomp_ext = comp_ext

  return
end subroutine scicomp_ext

!=====================================================================

subroutine sci_mix_frac( mode, es_names, ncv, sz, sworkspace, smf )

  use ci_prec
  implicit none

  integer, intent(in) :: mode, ncv
  real(k_sp), intent(in)  :: sz(ncv)
  real(k_sp), intent(inout) :: sworkspace(2*ncv)
  real(k_sp), intent(out) :: smf
  character*(*), intent(in) :: es_names(3)

  ! local variables
  real(k_dp) :: z(ncv), workspace(2*ncv), mf

  ! convert to double precision
  z = sz
  workspace = sworkspace

  ! call ci_mix_frac
  call ci_mix_frac( mode, es_names, ncv, z, workspace, mf )

  ! convert back
  sworkspace = workspace
  smf = mf

  return
end subroutine sci_mix_frac
  
!=====================================================================
subroutine sci_dens_temp( ncv, sc, sdpt )

  use ci_dat
  use ci_prec
  implicit none

  integer, intent(in)       :: ncv
  real(k_sp), intent(in)    :: sc(ncv)
  real(k_sp), intent(inout) :: sdpt(3)

  ! local variables
  real(k_dp) :: cc(ncv), dpt(3)

  ! convert to double precision
  cc = sc
  dpt(2) = sdpt(2)

  ! call ci_dens_temp
  call ci_dens_temp( ncv, cc, dpt )

  ! convert back
  sdpt = dpt
  
  return
end subroutine sci_dens_temp

!=====================================================================
subroutine sci_dens_temp_update( statol, success, sstats )

  use ci_dat
  use ci_prec
  use ci_stats
  implicit none

  real(k_sp), intent(in)  :: statol
  logical, intent(out)    :: success
  real(k_sp), intent(out) :: sstats(10)

  ! local variables
  real(k_dp) :: tatol, stats(10)

  tatol = statol
  call ci_dens_temp_update( tatol, success, stats)
  sstats = stats

  return
end subroutine sci_dens_temp_update
!=====================================================================
