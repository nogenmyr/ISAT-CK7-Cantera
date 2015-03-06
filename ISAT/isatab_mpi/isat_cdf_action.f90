!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module isat_cdf_action

contains

subroutine isat_cdf_act( table, mode, info, rinfo )

!  Perform required actions on CDF's

!  mode = 25 - act on cdf_err  (the CDF of the interpolation error)

!  info(80) = 1 - initialize or re-initialize CDF's using default parameters
!           = 2 - initialize or re-initialize CDF's using parameters from info and rinfo
!           = 3 - force output of CDF's
!           = 4 - turn CDF formation OFF
!           = 5 - turn CDF formation back ON

   use isat_types
   use isat_cdf

   implicit none

   type (table_type), pointer :: table
   integer, intent(in)        :: mode, info(l_info)
   real(k_xf), intent(in)     :: rinfo(l_rinfo)

   integer :: action
   integer, parameter :: nbin = 10000, op_inc=1000
   real(k_xf)    :: x_lower, x_upper
   character(30) :: blank, head, tail, name1

if( mode < 25  .or.  mode > 25 ) then  !  check mode
   if( table%write_log ) write(table%lu_log,'(a,i8)') &
      'isat_cdf_act: Warning, called for invalid mode = ', mode
   return
endif

action = info(80)                      ! check action = info(80)
if( action < 1  .or.  action > 5 ) then
   if( table%write_log ) write(table%lu_log,'(a,i8)') &
       'isat_cdf_act: Warning, called for invalid action: info(80) = ', action
   return
endif
      
blank  = repeat(' ',30)

if( mode == 25 ) then  !  cdf_err  ----------------------------------------------------

   if( action == 1  .or.  action ==2 ) then  !  initialize (or re-initialize)

      if( table%cdf_err%initialized ) call isat_cdf_kill( table%cdf_err )
      if( associated(table%cdf_err) ) deallocate(table%cdf_err)

	  head  = blank
	  head  = 'isat_err'
	  tail  = blank
	  tail  = 'cdf'
	  call isat_file_name( head, table%idtab, table%idproc, tail, name1 )
   endif

   if( action == 1 ) then   !  initialize with default parameters

      x_lower = table%etola * 1.e-3
      x_upper = table%etola * 1.e3
	  call isat_cdf_init( nbin, x_lower, x_upper, .true., name1, op_inc, table%cdf_err )

	  if( table%write_log ) write(table%lu_log,'(a)') &
	    'cdf_err (re-)initialized with default parameters' 

   elseif( action == 2 ) then  !  initialize with parameters from info/rinfo

   	  if( table%write_log ) write(table%lu_log,'(a,2i8,1p,2e13.4)') &
	    'cdf_err about to be (re-)initialized with nbin, op_inc, x_lower, x_upper = ', &
		info(81), info(84), rinfo(41), rinfo(44)  

      call isat_cdf_init( info(81), rinfo(41), rinfo(44), .true., name1, info(84), table%cdf_err )
   
   elseif( action == 3 ) then  !  force output

      if( table%cdf_err%initialized ) then
         if( table%cdf_err%samples > 1.d0 ) then 
	        call isat_cdf_op( table%cdf_err )
		    if( table%write_log ) write(table%lu_log,'(a)') &
			     'isat_cdf_act: cdf_err has been output'
	     else
	     	 if( table%write_log ) write(table%lu_log,'(a)') &
		          'isat_cdf_act: Warning - cdf_err is empty, cannot output'
	     endif
      else
	     if( table%write_log ) write(table%lu_log,'(a)') &
		  'isat_cdf_act: Warning - cdf_err has not been initialized, cannot output'
	  endif

   elseif( action == 4 ) then !  turn OFF

      if( table%cdf_err%initialized ) then
	     table%cdf_err%on = .false.
		 if( table%write_log ) write(table%lu_log,'(a)') &
			   'isat_cdf_act: cdf_err has been turned off'
      else
	     if( table%write_log ) write(table%lu_log,'(a)') &
		  'isat_cdf_act: Warning - cdf_err has not been initialized, cannot turn off'
	  endif
      
   elseif( action == 5 ) then !  turn ON

      if( table%cdf_err%initialized ) then
	     table%cdf_err%on = .true.
		 if( table%write_log ) write(table%lu_log,'(a)') &
			   'isat_cdf_act: cdf_err has been turned on'
      else
	     if( table%write_log ) write(table%lu_log,'(a)') &
		  'isat_cdf_act: Warning - cdf_err has not been initialized, cannot turn on'
	  endif

   endif

endif

return
end subroutine isat_cdf_act

end module isat_cdf_action