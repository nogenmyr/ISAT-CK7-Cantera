!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module isat_change_m

   use isat_defaults
   use isat_types
   use isat_subs
   implicit none

   integer :: i_signal = -12345
   real    :: r_signal = -12345.

contains

subroutine isat_change( table, info, rinfo, nchanges )

!  Change the values of some parameters set through info and rinfo.
!  Set info(80)=1 to suppress output.

   implicit none

   type (table_type), pointer :: table
   integer, intent(in)        :: info(l_info)
   real(k_xf), intent(in)     :: rinfo(l_rinfo)
   integer, intent(out)       :: nchanges  ! number of changes

   integer    :: lu_log, ichange(l_info), ibefore(l_info), rchange(l_rinfo), &
                 my_ichout, my_isatop, id

   real(k_xf)    :: rbefore(l_rinfo), etol_rat, etol_ratsq, &
                    len_rat

   lu_log = -1
   if( table%write_log  .and.  info(80)==0 ) lu_log = table%lu_log

!---------  make changes based on info  ---------------------------------

   ichange = 0
   ibefore = 0

   call isat_ichange( table%if_g,          2, info, lu_log, ichange, ibefore )
   call isat_ichange( table%if_h,          3, info, lu_log, ichange, ibefore )
   call isat_ichange( table%pr_mru,        4, info, lu_log, ichange, ibefore )
   call isat_ichange( table%pr_mfu,        5, info, lu_log, ichange, ibefore )
   call isat_ichange( table%ret_min,       6, info, lu_log, ichange, ibefore )
   call isat_ichange( table%ret_max,       7, info, lu_log, ichange, ibefore )
   call isat_ichange( table%grow_min,      8, info, lu_log, ichange, ibefore )
   call isat_ichange( table%grow_max,      9, info, lu_log, ichange, ibefore )

   call isat_ichange( table%ichout,       11, info, lu_log, ichange, ibefore )
   call isat_ichange( table%isatop,       12, info, lu_log, ichange, ibefore )

   call isat_ichange( table%ifull,        21, info, lu_log, ichange, ibefore )
   call isat_ichange( table%istats,       22, info, lu_log, ichange, ibefore )
   call isat_ichange( table%no_pr,        23, info, lu_log, ichange, ibefore )
   call isat_ichange( table%no_sr,        24, info, lu_log, ichange, ibefore )
   call isat_ichange( table%no_grow,      25, info, lu_log, ichange, ibefore )
   call isat_ichange( table%no_add,       26, info, lu_log, ichange, ibefore )
   call isat_ichange( table%no_dev,       27, info, lu_log, ichange, ibefore )
   call isat_ichange( table%idites,       28, info, lu_log, ichange, ibefore )
   call isat_ichange( table%defer,        29, info, lu_log, ichange, ibefore )

   call isat_ichange( table%grow_mode,    35, info, lu_log, ichange, ibefore )
   call isat_ichange( table%mode_con,     36, info, lu_log, ichange, ibefore )
   call isat_ichange( table%mode_eoi,     37, info, lu_log, ichange, ibefore )
   call isat_ichange( table%mode_shrink,  38, info, lu_log, ichange, ibefore )

   call isat_ichange( table%degen_r,      39, info, lu_log, ichange, ibefore )
   call isat_ichange( table%degen_g,      40, info, lu_log, ichange, ibefore )
   call isat_ichange( table%degen_s,      41, info, lu_log, ichange, ibefore )
   call isat_ichange( table%degen_c,      42, info, lu_log, ichange, ibefore )
   call isat_ichange( table%ecc_act,      43, info, lu_log, ichange, ibefore )

   call isat_ichange( table%ad_look,      44, info, lu_log, ichange, ibefore )
   call isat_ichange( table%ad_use,       45, info, lu_log, ichange, ibefore )
   call isat_ichange( table%ad_shrink,    46, info, lu_log, ichange, ibefore )

   call isat_ichange( table%aff_na,       48, info, lu_log, ichange, ibefore )
   call isat_ichange( table%aff_lim,      49, info, lu_log, ichange, ibefore )

   call isat_ichange( table%sr_mode,      50, info, lu_log, ichange, ibefore )
   call isat_ichange( table%gr_mode,      51, info, lu_log, ichange, ibefore )
   call isat_ichange( table%pair_cover,   52, info, lu_log, ichange, ibefore )
   call isat_ichange( table%cut_add,      53, info, lu_log, ichange, ibefore )
   call isat_ichange( table%cut_rem,      54, info, lu_log, ichange, ibefore )
   call isat_ichange( table%cut_upd,      55, info, lu_log, ichange, ibefore )

   call isat_ichange( table%icheck,       60, info, lu_log, ichange, ibefore )
   call isat_ichange( table%n_spool,      62, info, lu_log, ichange, ibefore )
   call isat_ichange( table%n_add_op,     63, info, lu_log, ichange, ibefore )

!---------  make changes based on rinfo  ---------------------------------

   rchange = 0
   rbefore = 0

   call isat_rchange( table%etola,         1, rinfo, lu_log, rchange, rbefore )
   call isat_rchange( table%etolc,         3, rinfo, lu_log, rchange, rbefore )
   call isat_rchange( table%eoa_lim,       4, rinfo, lu_log, rchange, rbefore )
   call isat_rchange( table%eoi_lim,       5, rinfo, lu_log, rchange, rbefore )
   call isat_rchange( table%ret_frac,      6, rinfo, lu_log, rchange, rbefore )
   call isat_rchange( table%grow_frac,     7, rinfo, lu_log, rchange, rbefore )
   call isat_rchange( table%stomby,        8, rinfo, lu_log, rchange, rbefore )
   call isat_rchange( table%outinc,        9, rinfo, lu_log, rchange, rbefore )
   call isat_rchange( table%chk_inc,      10, rinfo, lu_log, rchange, rbefore )
   call isat_rchange( table%chk_full,     11, rinfo, lu_log, rchange, rbefore )
   call isat_rchange( table%aff_adds,     12, rinfo, lu_log, rchange, rbefore )
   call isat_rchange( table%aff_inc,      13, rinfo, lu_log, rchange, rbefore )
   call isat_rchange( table%aff_thresh,   14, rinfo, lu_log, rchange, rbefore )
   call isat_rchange( table%aff_ratio,    15, rinfo, lu_log, rchange, rbefore )
   call isat_rchange( table%ad_theta,     16, rinfo, lu_log, rchange, rbefore )
   call isat_rchange( table%pair_sep,     17, rinfo, lu_log, rchange, rbefore )
   call isat_rchange( table%slow_prog_thresh, 18, rinfo, lu_log, rchange, rbefore )
   call isat_rchange( table%ret_frac_slow,    19, rinfo, lu_log, rchange, rbefore )
   call isat_rchange( table%grow_frac_slow,   20, rinfo, lu_log, rchange, rbefore )
   call isat_rchange( table%ecc_freq,         21, rinfo, lu_log, rchange, rbefore )

!---------  calculate nchanges,  return if zero  --------------------------

   nchanges = sum(ichange) + sum(rchange)
   if( nchanges == 0 ) return

!----------  address consequences of changes to info  ---------------------

   if( ichange(11) == 1 ) then   ! ichout has changed ---------------------

      my_ichout = 0       ! new value of my_ichout: =1 for checkpoint file
      if( table%myrank == 0 ) then
         if( table%ichout > 0  ) my_ichout = 1
      else
         if( table%ichout == 2 ) my_ichout = 1
      endif

   if( my_ichout == 0 ) then  !  set table%nextch
	     table%nextch = huge( table%nextch )
	  elseif( table%my_ichout == 0 ) then
         table%nextch = 0  !  force checkpointing on next call
      endif

	  table%my_ichout = my_ichout
   endif

   if( ichange(12) == 1 ) then   !  isatop has changed

      my_isatop = 0       ! new value of my_isatop: =1 for op file for this process
      if( table%myrank == 0 ) then
         if( table%isatop > 0  ) my_isatop = 1
      else
         if( table%isatop == 2 ) my_isatop = 1
      endif

	  if( my_isatop == 0  .and.  table%my_isatop == 1 ) then
	     close( table%lu_op )
		 table%lu_op   = 0
		 table%isatop_everon = .true.
	  elseif( my_isatop == 1  .and.  table%my_isatop == 0 ) then
         call isat_lu( table%lu_op )
		 if( table%isatop_everon ) then  !  append to existing file
            open( table%lu_op , file = table%isat_op  , status = 'old', &
               action = 'write', position = 'append' )
         else                      !  create new file
            open( table%lu_op , file = table%isat_op  , status = 'replace', &
               action = 'write' )
         endif
		 table%nextop = 0  !  force output on next call
	  endif

	  table%my_isatop = my_isatop
   endif
   
   if( ichange(49) == 1 ) call isat_na_max( table%aff_lim, table%nx, table%na_max )

   if( ichange(60) == 1 ) then  !  icheck has changed
      call id_set_check( table%idlist,     table%icheck )
	  call sell_set_check( table%seoa,     table%icheck )
	  call sell_set_check( table%speoa,    table%icheck )
	  call sell_set_check( table%seoi,     table%icheck )
	  call sell_set_check( table%speoi,    table%icheck )
	  call bt_param_set(   table%eoa_bt,   table%icheck ) 
	  call ll_param_set(   table%eoa_mru,  table%icheck ) 
	  call ll_param_set(   table%eoa_mru,  table%icheck ) 
	  call ebt_param_set(  table%peoa_ebt, check=table%icheck ) 
	  call ebt_param_set(  table%peoi_ebt, check=table%icheck ) 
   endif

!----------  address consequences of changes to rinfo  ---------------------

   if( rchange(8) == 1 ) then  !  stomby  has changed --------------

      if( table%stomby > 0.d0 ) then  !  based on specified storage
	     call isat_leaves( table%stomby, table%nx, table%nf, table%nh, table%na,  &
                  table%leaves, nint(table%stats(14)), table%maxleaves )
	  else                            !  specified number of leaves
	     table%maxleaves = -table%stomby
	  endif

	  if( table%maxleaves > table%leaves ) then
	     table%full = .false.
      else
	     table%maxleaves = table%leaves
	     table%full      = .true.
      endif
	  table%stats(13)    = table%maxleaves 
   endif

   if( rchange(3) == 1 ) then  !  etolc  has changed --------------
      table%etol_consq = table%etolc**2
!        Note: if  etolc  is decreased, no changes are made to the EOA.
!        To re-impose  the etolc  criterion  requires an SVD for each leaf.
!        Consider introducing this as an option.
   endif

   etol_rat   = 1.d0
   len_rat    = 1.d0

   if( rchange(1) == 1 ) then  !  etola has changed  ---------------
      etol_rat    = table%etola / rbefore(1)  !  etola_new / etola_old
	  table%etolr = table%etolr * etol_rat    !  change etolr in proportion
	  len_rat     = sqrt( etol_rat )
	  etol_ratsq  = etol_rat**2

	  do id = 1, table%n_pt  !  reset etolsq for all leaves
         if( associated( table%leaf_pt(id)%leaf ) )  &
		    table%leaf_pt(id)%leaf%etolsq = table%leaf_pt(id)%leaf%etolsq * etol_ratsq
	  end do
   endif

   if( len_rat < 1.d0 ) then  !  EOA/EOI are to be shrunk
      table%stats(15:19) = table%stats(15:19) * len_rat
	  table%stats(72)    = 0.d0  !  reset errmax
! shrink all leaf ELLs by factor len_rat
	  call sell_stretch( table%seoa,  len_rat ) 
	  call sell_stretch( table%seoi,  len_rat ) 
	  call sell_stretch( table%speoi, len_rat ) 
	  call sell_stretch( table%speoa, len_rat ) 

	  if( info(81) == 1  .or.  info(81) == 2 ) then  !  rebuilt EBTs
         call ebt_rebuild( table%peoa_ebt, info(81), table%pair_cover )
         call ebt_rebuild( table%peoi_ebt, info(81), table%pair_cover )
	  endif
   endif
   
   if( rchange(21) == 1 ) then     ! ecc_freq  has changed
      if( table%ecc_freq <= -1.d0 ) then
         table%ecc_q = huge(1.d0)  !  turn off ECC
      else
         table%ecc_q = 0.d0        !  force ECC on next query
      endif
   endif

   return

end subroutine isat_change

!===========================================================================

subroutine isat_ichange( iset, i, info, lu_log, ichange, ibefore )

! Attempt to change the value of  iset  to  info(i).
! The new value is inew=info(i), or if info(i)=0 then inew=info_def(i).
! Make no change if  inew=iset.
! Make no change if  info(i)  violates specified bounds.
! If imax>0, info(i) must satisfy  imin <= info(i) <= imax.
! If imax=0, info(i) must be non-negative.
! (No check performed for imax < 0.)
! If lu_log >=0, record changes and errors on lu_log.
! ichange(i) is set to 0 or 1 if there is no change or change, respectively, to iset.
! ibefore(i) is set to the initial value of iset.

! Note: if info(i)=0, then iset is set to the default.
! Thus, variables with non-zero defaults cannot be reset to zero (in one call).
 
   implicit none
   integer, intent(in)    :: i, info(i), lu_log 
   integer, intent(inout) :: iset, ichange(i), ibefore(i)

   integer :: imin, imax, inew
   logical :: valid

   ichange(i) = 0
   ibefore(i) = iset

   if( info(i) == i_signal ) return    ! no change for signalled value
   inew = info(i)
   if( inew == 0 ) inew = info_def(i)  !  use default value
   if( inew == iset ) return           !  no change

   valid = .true.  !  check validity of new value
   imin  = info_min(i)
   imax  = info_max(i)

   if( imax > 0 ) then
      if( inew < imin  .or.  inew > imax ) valid = .false.
   elseif( imax == 0 ) then
      if( inew <0 ) valid = .false.
   endif

   if( valid ) then
      if( lu_log >=0 ) write(lu_log,'(a,i3,a,i9,a,i9)') 'isat_ichange: info(i) for i=', &
	                                         i, ' changed from ', iset, ' to ', inew
      iset       = inew
	  ichange(i) = 1
   elseif( lu_log >=0 ) then
      write(lu_log,'(a,i3,a,i9)') 'isat_ichange: ERROR: attempt to change info(i) for i=',&
		                           i, ' to the invalid value ', inew
   endif

   return

end subroutine isat_ichange

!===========================================================================

subroutine isat_rchange( rset, i, rinfo, lu_log, rchange, rbefore )

! Attempt to change the value of  rset  to  rinfo(i).
! The new value is rnew=rinfo(i), or if rinfo(i)=0. then rnew=rinfo_def(i).
! Make no change if  rnew=rset.
! Make no change if  rinfo(i)  is less than  rinfo_min(i).
! If lu_log >=0, record changes and errors on lu_log.
! rchange(i) is set to 0 or 1 if there is no change or change, respectively, to rset.
! rbefore(i) is set to the initial value of rset.

! Note: if rinfo(i)=0., then rset is set to the default.
! Thus, variables with non-zero defaults cannot be reset to zero (in one call).
 
   implicit none
   integer, intent(in)       :: i, lu_log 
   real(k_xf), intent(in)    :: rinfo(i)

   integer, intent(inout)    :: rchange(i)
   real(k_xf), intent(inout) :: rset, rbefore(i)

   real(k_xf) :: rnew

   rchange(i) = 0
   rbefore(i) = rset

   if( rinfo(i) == r_signal ) return     ! no change for signalled value
   rnew = rinfo(i)
   if( rnew == 0. ) rnew = rinfo_def(i)  !  use default value
   if( rnew == rset ) return             !  no change

   if( rnew < rinfo_min(i) ) then
      if( lu_log >=0 ) write(lu_log,'(a,i3,a,1p,e13.4,a,1p,e13.4)')  &
	            'isat_rchange: ERROR: attempt to change rinfo(i) for i=',&
		         i, ' to rnew = ', rnew, ' which is less than rinfo_min(i) = ', rinfo_min(i)
	  return
   endif
      
   if( lu_log >=0 ) write(lu_log,'(a,i3,a,1p,e13.4,a,1p,e13.4)') &
          'isat_rchange: rinfo(i) for i=', i, ' changed from ', rset, ' to ', rnew
   rset       = rnew
   rchange(i) = 1
   return

end subroutine isat_rchange

end module isat_change_m
