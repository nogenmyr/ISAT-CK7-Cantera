!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module isat_io

use isat_types
use isat_cdf
use isat_subs
 
contains

!===========================================================================

subroutine isat_table_write( table )

!  write table

   implicit none

   type (table_type), pointer :: table
   integer :: lu, n_leaves, id

!--- check that table and leaves exist
   if( .not. associated(table)  ) &
      call isat_abort('isat_table_write',1,mess='table not associated')

   if( table%leaves <= 0 ) &
      call isat_abort('isat_table_write',2,mess='empty, no leaves=' )

!--- check consistency of datastructures
	if( table%icheck >= 2 ) then
	   call isat_integrity( table, table%lu_log, id ) 
	   if( id /= 0 ) call isat_abort('isat_table_write', 3, isv = id, &
            mess = 'inconsistent data structures, info = ' )
    endif

!--- open file
   call isat_lu( lu )
   table%lu_dat = lu
   open( lu, file = table%isat_tab, &
         form = 'unformatted', err = 100 )

!--- write idtab and leaves
   write(lu) table%idtab, table%leaves

!--- write parameters that must match when file is read
   write(lu) table%nx, table%nf, table%nh, &
             table%iscale, table%isat_vers, table%n_spool

!--- write error tolerance
   write(lu) table%etola

!--- write affine space
   write(lu) table%na
   write(lu) table%ua(:,:)

   if( table%iscale /= 0 ) write(lu) table%xscale(:), table%fscale(:)

!---  write stats
   write(lu) table%stats
   write(lu) table%i_hist
   if( table%i_hist > 0 ) write(lu) table%stats_hist(:,1:table%i_hist)

!--- spooled output
   write(lu) table%i_spool
   if( table%i_spool > 0 ) write(lu) table%spool(1:100,1:table%i_spool)

!--- other quantities needed for perfect restart
   write(lu) table%deferred

!--- write n_pt and leaves
   write(lu) table%n_pt
   do id = 1, table%n_pt
      if( associated( table%leaf_pt(id)%leaf ) )  &
         call isat_leaf_write( table%leaf_pt(id)%leaf )
   end do

!--- write all data structures
   call id_write( table%idlist, lu )

   call sell_write( table%seoa, lu )
   call sell_write( table%speoa, lu )
   call sell_write( table%seoi, lu )
   call sell_write( table%speoi, lu )

   call bt_write( table%eoa_bt, lu, n_leaves )
   call ll_write( table%eoa_mru, lu )
   call ll_write( table%eoa_mfu, lu )
   call ebt_write( table%peoa_ebt, lu, n_leaves ) 
   call ebt_write( table%peoi_ebt, lu, n_leaves ) 

   if( table%cdf_error > 0 ) then
      call isat_cdf_write( lu, table%cdf_err )
   else
      write(lu) .false.
   endif

!--- report success
   if( table%write_log) write(table%lu_log,"(a,4i8)") &
         'table written: idtab, leaves, grows=', &
	 table%idtab, table%leaves, nint( table%stats(44) )
   close( lu )

   return

100  call isat_abort('isat_table_write',4,mess='error opening file', chv= table%isat_tab )

end subroutine isat_table_write

!============================================================================

subroutine isat_leaf_write( leaf )

   implicit none

   type (leaf_type), pointer :: leaf
   integer :: lu

   lu = leaf%table%lu_dat

   write(lu) leaf%id
   write(lu) leaf%etolsq, leaf%xfh
   write(lu) leaf%g, leaf%props

   return
end subroutine isat_leaf_write

!===========================================================================

subroutine isat_table_read( table )

!  read table
   use isat_subs
   implicit none

   type (table_type), pointer :: table
   type (leaf_type),  pointer :: leaf

   integer    :: lu, lu_log, i, idtab, leaves, nx, nf, nh, na, &
                 iscale, isat_version, n_spool, n_pt, id, n_leaves
   real(k_xf) :: xscale(table%nx), fscale(table%nf)

!--- check that table exists
   if( .not. associated(table)  )  &
      call isat_abort('isat_table_read',1,mess='table not associated' )

!--- report start of reading
   if( table%write_log) then
      lu_log = table%lu_log
      write(lu_log,"(a)")' '
      write(lu_log,"(a,i10)")'starting to read table, idtab=', &
                                table%idtab
      write(lu_log,"(a)")' '
   else
      lu_log = 0
   endif

!--- open file
   call isat_lu( lu )
   table%lu_dat = lu
   open( lu, file = table%isat_tab, status='old', action='read', &
         form = 'unformatted', err = 100 )

!--- read idtab and leaves
   read(lu,err=110,end=110) idtab, leaves
   if( idtab /= table%idtab ) then
      write(lu_log,"(a,2i6)")'isat_table_read: WARNING idtab mismatch', &
                              idtab, table%idtab
   endif
   if( leaves > table%maxleaves ) then
   
      if( .false. ) then  !SBP 9/5/2010 - changed from fatal error to warning
         call isat_abort('isat_table_read',2,mess='leaves > maxleaves', &
			                ivar = (/ leaves, table%maxleaves /) )
      else
         write(lu_log,"(a,2i6)")'isat_table_read: WARNING: leaves > maxleaves', &
                                 leaves, table%maxleaves
         table%maxleaves = leaves
         table%full = .true.
      endif
      
   elseif( leaves == table%maxleaves ) then
      table%full = .true.
   endif
   table%leaves = leaves

!--- read parameters that must match with existing data
   read(lu,err=120,end=120) nx, nf, nh, iscale, isat_version, n_spool 
   call check_i( 'nx    ',    table%nx    ,      nx    )
   call check_i( 'nf    ',    table%nf    ,      nf    )
   call check_i( 'nh    ',    table%nh    ,      nh    )
   call check_i( 'iscale',    table%iscale,      iscale)
   call check_i( 'isat_vers', table%isat_vers,   isat_version  )
   call check_i( 'n_spool',   table%n_spool  ,   n_spool  )

!--- read error tolerance
   read(lu,err=125,end=125) table%etola

   read(lu,err=130,end=130) na
   read(lu,err=140,end=140) table%ua(:,:)

   if( na /= table%na ) then
      table%nga = (na*(na+1))/2
	  table%na  = na
   endif

   if( table%iscale /= 0 ) then
      read(lu,err=150,end=150) xscale(:), fscale(:)
      do i = 1, nx
         call check_d( 'xscale', table%xscale(i), xscale(i) )
      end do
      do i = 1, nf
         call check_d( 'fscale', table%fscale(i), fscale(i) )
      end do
   endif

!----  read stats
   read(lu,err=160,end=160) table%stats
   read(lu,err=170,end=170) table%i_hist
   if( table%i_hist > 0 ) read(lu,err=180,end=180) table%stats_hist(:,1:table%i_hist)
   table%i_hist = -table%i_hist  ! set negative to indicate that stats_hist to be o/p

!--- spooled output
   read(lu,err=190,end=190) table%i_spool
   if( table%i_spool > 0 ) read(lu,err=200,end=200 ) table%spool(1:100,1:table%i_spool)

!-- other quantities needed for perfect restart
   read(lu,err=210,end=210 ) table%deferred

!--- read n_pt 
   read(lu) n_pt
   if( n_pt > table%n_pt ) then

!--- if necessary, re-initialize  table%leaf_pt
      do id = 1, table%n_pt
	     nullify( table%leaf_pt(id)%leaf )
	  end do
	  deallocate( table%leaf_pt )

	  allocate( table%leaf_pt( n_pt ) )
	  table%n_pt = n_pt
	  do id = 1, n_pt
         nullify( table%leaf_pt(id)%leaf )
	  end do

   endif

!--- read leaves
   do i = 1, leaves
      nullify( leaf )
      call isat_leaf_init( table, leaf )
      call isat_leaf_read( leaf, lu )
	  table%leaf_pt(leaf%id)%leaf => leaf
   end do

!--- read all data structures
   call id_read( table%idlist, lu )

   call sell_read( table%seoa, lu, table%idlist )
   call sell_read( table%speoa, lu, table%idlist  )
   call sell_read( table%seoi, lu, table%idlist  )
   call sell_read( table%speoi, lu, table%idlist  )

   call bt_read( table%seoa, table%eoa_bt, table%idlist, lu, n_leaves )
   call ll_read( table%eoa_mru, lu )
   call ll_read( table%eoa_mfu, lu )

   call ebt_read( table%speoa, table%peoa_ebt, table%idlist, lu, n_leaves ) 
   call ebt_read( table%speoi, table%peoi_ebt, table%idlist, lu, n_leaves ) 

   if( table%cdf_error == 1 ) call isat_cdf_read( lu, table%cdf_err )

!--- report success
   if( table%write_log) then
      write(lu_log,"(a,3i8)") 'table read:    idtab, trees, leaves=', &
	     table%idtab, table%leaves
      write(lu_log,"(a)")' '
      call isat_flush(lu_log)
      close( lu )
   endif

   return
100  call isat_abort('isat_table_read',3, mess='error opening file', chv= table%isat_tab )
110  call isat_abort('isat_table_read',4, mess='error reading idtab, leaves' )
120  call isat_abort('isat_table_read',5, mess='error reading nx' )
125  call isat_abort('isat_table_read',6, mess='error reading etola' )
130  call isat_abort('isat_table_read',7, mess='error reading na' )
140  call isat_abort('isat_table_read',8, mess='error reading ua' )
150  call isat_abort('isat_table_read',9, mess='error reading xscale(:), fscale(:)' )
160  call isat_abort('isat_table_read',10, mess='error reading stats' )
170  call isat_abort('isat_table_read',11, mess='error reading i_hist' )
180  call isat_abort('isat_table_read',12,mess='error reading stats_hist' )
190  call isat_abort('isat_table_read',13,mess='error reading i_spool' )
200  call isat_abort('isat_table_read',14,mess='error reading spool' )
210  call isat_abort('isat_table_read',15,mess='error reading deferred' )

contains

   subroutine check_i( name, i1, i2 )
   character(6) name
   integer i1, i2
   if( i1 == i2 ) return
   write(lu_err,*)'isat_table_read: mismatch for:  ', name, i1, i2
   call isat_abort('isat_table_read',9 )

   end subroutine check_i

   subroutine check_d( name, d1, d2 )
   character(6) name
   real(k_xf) d1, d2
   if( d1 == d2 ) return
   write(lu_err,*)'isat_table_read: mismatch for:  ', name, d1, d2
   call isat_abort('isat_table_read',10 )
   end subroutine check_d

end subroutine isat_table_read

!============================================================================

subroutine isat_leaf_read( leaf, lu )

   implicit none

   type (leaf_type), pointer :: leaf
   integer, intent(in)       :: lu

   read(lu,err=201,end=201) leaf%id
   read(lu,err=202,end=202) leaf%etolsq, leaf%xfh
   read(lu,err=203,end=203) leaf%g, leaf%props

   return

201   call isat_abort( 'isat_leaf_read', 2, mess='error reading  id' )
202   call isat_abort( 'isat_leaf_read', 3, mess='error reading  etolsq, x, f, h' )
203   call isat_abort( 'isat_leaf_read', 4, mess='error reading  g, props' )

end subroutine isat_leaf_read

!============================================================================

end module isat_io
