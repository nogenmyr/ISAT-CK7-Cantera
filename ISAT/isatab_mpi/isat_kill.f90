!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module isat_kill

use isat_types
use isat_cdf
 
contains

!===========================================================================

subroutine isat_table_kill( table, first_table )

!  kill table

   implicit none

   type (table_type), pointer :: table, first_table

   type (table_type), pointer :: this_table

   integer :: lu_log, idtab, id
   logical :: write_log

!--- check that table exists
   if( .not. associated(table)  ) then
      write(lu_err,*)'isat_table_kill: warning; tried to kill non-existent table'
      return
   endif

!--- report start of killing
   idtab     = table%idtab
   write_log = table%write_log
   if( write_log ) then
      lu_log = table%lu_log
      write(lu_log,"(a)")' '
      write(lu_log,"(a,3i8)") &
            'killing table: idtab, leaves=', &
             idtab, table%leaves

      write(lu_log,"(a)")' '
   else
      lu_log = 0
   endif

!---  remove table from list of tables
   if( associated(table,first_table) ) then
      if( associated(table%next_table) ) then
         first_table => table%next_table
      else
         nullify(first_table)
      endif
   else
      this_table => first_table
      do 
         if( .not.associated(this_table%next_table) ) then
             call isat_abort('isat_table_kill',1, &
			     mess='table not in list',isv=table%idtab )
         elseif( .not.associated(this_table%next_table,table) ) then
	       this_table => this_table%next_table
	       cycle
!  this_table is the predecessor of table
	     elseif( associated(table%next_table) ) then
	        this_table%next_table => table%next_table
	        exit
	     else
	        nullify( this_table%next_table )
	       exit
	     endif
      end do
   endif

!--- kill data structures
   call ebt_purge(   table%peoa_ebt )
   call ebt_purge(   table%peoi_ebt )
   call ebt_destroy( table%peoa_ebt )
   call ebt_destroy( table%peoi_ebt )
   call ll_destroy(  table%eoa_mru )
   call ll_destroy(  table%eoa_mfu )
   call bt_purge(    table%eoa_bt )
   call bt_destroy(  table%eoa_bt )

   call sell_destroy( table%seoa )
   call sell_destroy( table%speoa )
   call sell_destroy( table%seoi )
   call sell_destroy( table%speoi )
   call id_destroy(   table%idlist )

!---  kill leaves
   do id = 1, table%n_pt
      if( associated( table%leaf_pt(id)%leaf ) ) &
	     call isat_leaf_kill( table%leaf_pt(id)%leaf )  
   end do

!--- kill query
   call isat_query_kill( table%query )

!--- kill CDF's
   if( table%cdf_err%initialized  ) call isat_cdf_kill( table%cdf_err  )
    
!----  close all files (except lu_log which is closed below)
   if( table%lu_op > 0 ) close( table%lu_op )

!--- deallocate table 
   deallocate( table%leaf_pt, table%ua, table%xscale, table%fscale )
   deallocate( table%xsci, table%fsci, table%g2gs, table%gs2g, table%p )
   deallocate( table%dsq, table%ids )
   if( associated(table%spool) ) deallocate( table%spool )
			   
   deallocate( table )

!--- report success
   if( write_log) then
      write(lu_log,"(a,3i10)") &
            'completed killing table: idtab', idtab
      close( lu_log )
   endif

   return

end subroutine isat_table_kill  !------------------------------------------

subroutine isat_leaf_kill( leaf )

! kill leaf and deallocate all associated storage

   implicit none

   type (leaf_type), pointer :: leaf

   nullify( leaf%table )
   deallocate( leaf%xfh, leaf%g )
   deallocate( leaf)

   return
end subroutine isat_leaf_kill  !-------------------------------------------------------

subroutine isat_query_kill( query )

   implicit none
   type(query_type), pointer :: query

   nullify( query%table )
   deallocate( query%xs, query%xa, query%fs )
   deallocate( query )

   return
end subroutine isat_query_kill

!===========================================================================
 
end module isat_kill
