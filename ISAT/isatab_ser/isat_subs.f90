!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module isat_subs

!  This module contains the core routines used for retrieving, growing and 
!  adding in isatab.  The routines are:

!   RETRIEVING
! isat_primary_retrieve( query )
! isat_secondary_retrieve( query )
! isat_get_fgh( query, need_fs, need, fa, ga, ha )

!   ERROR CHECKING AND CORRECTION
! isat_ecc( query, accurate, errsq )

!   GROWING
! isat_grow( query ) 
! isat_leaf_degenerate( table, id )
! isat_leaf_grow( query, id, grown, inclusive, degenerate )
! isat_eoa_err( leaf, fs1, fs2, accurate, errsq )

!   ADDING
! isat_add( table, xs, fs, ga, ha, modret )
! isat_leaf_delete( id, table )
! isat_leaf_set( table, leaf, xs, fs, ga, ha )
! isat_leaf_init( table, leaf )
! isat_leaf_add( table, leaf, xs )  
! isat_leaf_geom_init( leaf, nx, nf, ng, xs, g_eoa, g_eoi ) 
! select_pts( nout, nin, x, x_lim, jpass )
! leaf_pt_double( table ) 

!   AFFINE SPACE
! isat_affine( table, sell, changed ) 
! isat_aff_rebuild( table, sell, spell, ebt )

!   TEST ROUTINES 
! isat_integrity( table, lu, info )
! isat_table_plot( table )   
! isat_leaf_plot( nx, c, gg, k )

!----------------------------------------------------------------------------------

use isat_kill
use sell_ll
use sell_ebt
use sell_bt

implicit none

private isat_leaf_set, isat_leaf_add, isat_leaf_geom_init, select_pts, &
        leaf_pt_double

contains  !------------------------------------------------------------------------

!=========  RETRIEVING=============================================================

subroutine isat_primary_retrieve( query )

!  Attempt to find leaf whose EOA covers  xs = query%xs.
!  If successful, id = query%id > 0 is set to the ID of the leaf.
!  If unsuccessful, query%id is negative, and query%xa is evaluated.
!  If query%id <= -2, then id_pl=-query%id is the ID of primary leaf.

! 1/ BT is traversed to identify primary leaf (unless no_bt > 0).
! 2/ If  xs  is in EOA of primary leaf: done.
! 3/ Form  xa , projection of xs onto affine space.
! 4/ Test first  pr_mru  (untested) entries in MRU cache.
! 5/ Test first  pr_mfu  (untested) entries in MFU cache.

! Note: pr_mru=-1 for exhaustive search; pr_mru<=-2 to suppress search;
!       similarly for pr_mfu

   type (query_type), pointer :: query  ! on entry query%xs is defined; on exit
!                                query%id is set; and, if id is negative, 
!                                query%xa is set

   type (table_type), pointer :: table  ! pointer to table

   integer :: id, nx, na, depth, mru_tests, mfu_tests, &
              eoa_tests, peoa_tests, id_pl
   logical :: in

   table    => query%table
   query%id = -1
   if( table%leaves < 1 ) return

   nx  = table%nx
   na  = table%na

search_methods: do  !  only loop once through BT, MRU, MFU

!--------  BT traverse
   if( table%no_bt == 0 ) then
      call bt_traverse( table%eoa_bt, query%xs, id, depth )
      id_pl = id

      table%stats(25) = table%stats(25) + 1.d0
      table%stats(27) = table%stats(27) + depth

      call sell_ell_pt_in( nx, table%seoa%ngeom, table%seoa%ell_pt(id)%ell%geom, &
                           query%xs, in )

      if( in ) then
         table%stats(26) = table%stats(26) + 1.d0
         exit search_methods
      endif
   endif

!--------- preliminaries for MRU and MFU caches

!  map query point onto affine space
   query%xa(1:na) = matmul( query%xs-table%ua(:,nx+1) , table%ua(:,1:na) )  

   if( table%pr_mru >=0 ) then  !  determine maximum number of MRU tests
      mru_tests = table%pr_mru
   elseif( table%pr_mru == -1 ) then
      mru_tests = table%leaves
   else
      mru_tests = 0
      if( table%pr_mfu <= -2 ) return  ! no MRU or MFU testing 
   endif

   if( table%pr_mfu >=0 ) then  !  determine maximum number of MFU tests
      mfu_tests = table%pr_mfu
   elseif( table%pr_mfu == -1 ) then
      mfu_tests = table%leaves
   else
      mfu_tests = 0
   endif

!  indicate that all EOAs which might be tested are currently untested...
   call ll_set_tested( table%eoa_mru, table%seoa, mru_tests+1, .false. )
   call ll_set_tested( table%eoa_mfu, table%seoa, mru_tests+mfu_tests+1, .false. )

!  ...except that the primary leaf may have been tested
   if( table%no_bt == 0 ) table%seoa%ell_pt(id)%ell%tested = .true.
   
   id = -1   !  indicate BT unsuccessful 

!---------  MRU cache

   if( mru_tests > 0 ) then

      if( na == nx ) then  ! use ll_query
         call ll_query( table%eoa_mru, table%seoa, query%xs, mru_tests, &
		                -1, id, eoa_tests )
         peoa_tests = 0

      else                 ! use ll_query_pair		 		 
		 call ll_query_pair( table%eoa_mru, table%speoa, table%seoa,  &
		                     query%xa(1:na), query%xs,  mru_tests, mru_tests, &
							 -1, id, peoa_tests, eoa_tests )
      endif

      table%stats(68) = table%stats(68) + eoa_tests
      table%stats(70) = table%stats(70) + peoa_tests
      table%stats(28) = table%stats(28) + 1.d0

      if( id > 0 ) then
         table%stats(29) = table%stats(29) + 1.d0
         table%stats(30) = table%stats(30) + peoa_tests
	     table%stats(69) = table%stats(69) + 1.d0
         exit search_methods
      endif

   endif

!---------  MFU cache

   if( mfu_tests > 0 ) then

      if( na == nx ) then  ! use ll_query  
         call ll_query( table%eoa_mfu, table%seoa, query%xs, mfu_tests,  &
		                -1, id, eoa_tests )
         peoa_tests = 0
      else  ! use ll_query_pair
	 		 
		 call ll_query_pair( table%eoa_mfu, table%speoa, table%seoa,  &
		                     query%xa(1:na), query%xs, mfu_tests, mfu_tests, &
							 -1, id, peoa_tests, eoa_tests )
      endif


      table%stats(68) = table%stats(68) + eoa_tests
      table%stats(70) = table%stats(70) + peoa_tests
      table%stats(31) = table%stats(31) + 1.d0

      if( id > 0 ) then
         table%stats(32) = table%stats(32) + 1.d0
         table%stats(33) = table%stats(33) + peoa_tests
	  	 table%stats(69) = table%stats(69) + 1.d0
         exit search_methods
      endif

   endif

   query%id = -id_pl
   return !  failed to find EOA covering xs

end do search_methods

! success
   query%id = id  !  leaf whose EOA covers query point found 
   table%leaf_pt(id)%leaf%props(2) = table%stats(1)  !  query of last retrieve

   return
end subroutine isat_primary_retrieve  !--------------------------------------------

subroutine isat_secondary_retrieve( query )

!  Attempt to find a leaf whose EOA covers  xs = query%xs.
!  If successful, query%id is set to the (positive) ID of the leaf.
!  If unsuccessful, query%id is set to  -1.
!  On entry, query%xs and query%xa are defined.

   type (query_type), pointer :: query ! query%xs and query%xa are defined on entry.

   type (table_type), pointer :: table ! pointer to table

   type(be_type), pointer :: be
   integer    :: id, max_f_tests, max_be_tests, f_tests, a_tests, be_tests
   real(k_xf) :: cpu_f, cpu_sr, n_f, n_eoa, sr_success, sr_total, ret_frac

   table => query%table
   nullify( be )

   max_be_tests = 3*table%leaves  ! unlimited BE testing

! determine maximum allowed number of EOA tests (max_f_tests)

   cpu_f  = table%stats(91)   !  CPU for f evals
   cpu_sr = table%stats(85)   !  CPU for SR (only accounts for successful SR)
   
   if( query%slow_prog ) then
      ret_frac = table%ret_frac_slow   
   else
      ret_frac = table%ret_frac
   endif

   if( ret_frac > 0.d0  .and.  cpu_sr > 0.1d0 .and. cpu_f > 0.1d0 ) then
!     limit so that CPU time is no more than ret_fac * cpu_f
      sr_success = table%stats(3)                   ! number of successful SRs
	  sr_total   = table%stats(1) - table%stats(2)  !  total number of SR

      n_eoa  = table%stats(36)   !  EOA tests in SR
	  n_eoa  = n_eoa * sr_success / sr_total  !  estimate of EOA tests in successful SR
	  cpu_sr = cpu_sr / n_eoa    !  CPU per EOA test

      n_f    = table%stats(77)   !  f evals
      cpu_f  = cpu_f / n_f       !  CPU per f eval

      max_f_tests = ret_frac * cpu_f / cpu_sr
	  
      if( table%ret_max > 0 ) max_f_tests = min( max_f_tests, table%ret_max )
      max_f_tests = max( max_f_tests, table%ret_min )

   elseif( table%ret_max < 0 ) then
      max_f_tests = table%leaves  !  effectively unlimited

   else
      if( table%ret_max > 0 ) then
         max_f_tests = table%ret_max  !  set to ret_max
	  else
         max_f_tests = table%leaves   !  set unlimited
	  endif
   endif

   max_f_tests = min( max_f_tests, table%leaves )  ! upper bound on possible tests

   table%stats(38) = max_f_tests

   call ebt_query_pair( table%seoa, query%xs, table%peoa_ebt, query%xa, table%sr_mode, &
                        max_f_tests, max_be_tests, be, f_tests, a_tests, be_tests ) 

   table%stats(34) = table%stats(34) + 1.d0
   table%stats(36) = table%stats(36) + f_tests
   table%stats(37) = table%stats(37) + be_tests
   table%stats(68) = table%stats(68) + f_tests
   table%stats(70) = table%stats(70) + a_tests 

   if( .not.associated(be) ) then  !  unsuccessful
      query%id = -1
	  if( f_tests >= max_f_tests )  then   ! unsuccessful incomplete search
	     table%stats(39) = table%stats(39) + 1.d0 
      endif
      return
   endif

   table%stats(35) = table%stats(35) + 1.d0
   table%stats(69) = table%stats(69) + 1.d0

   id       = be%id
   query%id = id

   table%leaf_pt(id)%leaf%props(2) = table%stats(1)  !  query of last retrieve

   return
end subroutine isat_secondary_retrieve  !------------------------------------------

subroutine isat_get_fgh( query, need_fs, need, fa, ga, ha )

!  Selectively, obtain from leaf values of fs, fa, ga and ha.

   type (query_type), pointer :: query    !from which to retrieve
   integer, intent(in)        :: need_fs  ! = 1 if fs is to be set
   integer, intent(in)        :: need(3)  ! need = (1,1,1) if (fa,ga,ha),
!                                           respectively, are to be set
   real(k_xf)                 :: fa(query%table%nf)
   real(k_xf)                 :: ga(query%table%nf, query%table%nx)
   real(k_xf)                 :: ha(query%table%nh)
   real(k_xf)                 :: dxs(query%table%nx)

!  If  fs  is to be determined, query%xs must be set on entry.
!  fa  is determined from  fs.  Hence, on entry, either query%fs 
!  must be set, or need_fs=1 is required.

   type (table_type), pointer :: table
   type (leaf_type),  pointer :: leaf

   table => query%table
   leaf  => table%leaf_pt( query%id )%leaf

!  set  fs  if required
   if( need_fs == 1 ) then
      dxs = query%xs - leaf%xfh(1:table%nx)
      ! check if it is a direct hit i.e., query = leaf
      if( maxval(abs(dxs)) == 0.d0 ) then
         query%fs =  leaf%xfh(table%nxp1:table%jfe)
      else
         query%fs =  leaf%xfh(table%nxp1:table%jfe) + matmul(leaf%g, dxs)
      endif
   endif

!  set  fa  if required
   if( need(1) == 1 ) then
      if( table%iscale == 0 ) then
	     fa = query%fs
	  else
	     fa = query%fs * table%fscale
	  endif
   endif

!  set  ga  if required
   if( need(2) == 1 ) then  
      if( table%iscale == 0 ) then
         ga = leaf%g		
      else
         ga = leaf%g * leaf%table%gs2g 	
      endif
   endif

!  set  ha  if required
   if( need(3) == 1  .and.  table%nh > 0 ) ha = leaf%xfh(table%jhs:table%jhe)

   return
end subroutine isat_get_fgh  !-----------------------------------------------------

!=================  ERROR CHECKING AND CORRECTION (ECC) ===============================

subroutine isat_ecc( query, accurate, errsq )

!  Error checking and correction.
!  This routine is invoked occassionally on successful retrieves.
!  The error is measured, and if it exceeds the tolerance, then the EOA is shrunk.

   use isat_rnu
   type (query_type), pointer :: query    ! query%xs and query%xa are defined on entry.
   logical, intent(in)        :: accurate ! true if error less than tolerance
   real(k_xf), intent(in)     :: errsq    ! square of error
   
   type (table_type), pointer :: table ! pointer to table
   type (ell_type),   pointer :: eoa   ! EOA
   type (ell_type),   pointer :: pell  ! PEOA

   integer    :: id, nx, ng, na, nga, ngeom, info, j
   real(k_xf) :: x_shrink(query%table%nx), prob, cpu_q, cpu_f, ranu, r_eoa(2)
   real(k_xf), parameter :: shrink = 0.99d0  !  attenuation of shrink point

   table => query%table

! specify prob of ECC testing, then next query on which testing is to be performed

   if( table%ecc_freq < 0.d0 ) then
      prob = -table%ecc_freq  ! specified probability of testing
   
   elseif( table%stats(4) > 0.d0 ) then
!  set ECC test probability so that ECC testing frequency is the 
!     fraction ecc_freq of the growing frequency    
      prob  = table%ecc_freq * table%stats(4) / table%stats(1)
!  attenuate probability if already over target
      if( table%stats(4) > 100.d0  .and.  table%stats(79) > table%ecc_freq * table%stats(4) )  &
         prob = prob * (table%ecc_freq * table%stats(4) / table%stats(79) )**2
   else
      prob  = 0.1d0 * table%ecc_freq  ! insufficient grows for estimate
   endif

! reset  table%ecc_q   
   call rnu( ranu )  
   ranu = -log(1.d-5+ranu)/prob
   table%ecc_q = table%stats(1) + ranu
   
   table%stats(53) = prob
   table%stats(79) = table%stats(79) + 1.d0  ! increment counter of acc. tests
   
   if( accurate ) return  !  all done if accurate
   
!  shrink EOA based on inaccurate query point

   id    =  query%id
   eoa   => table%seoa%ell_pt(id)%ell 

   nx    = table%nx
   ng    = table%ng
   na    = table%na
   nga   = table%nga 
   ngeom = table%seoa%ngeom 
   r_eoa = eoa%geom(ngeom-1:ngeom) ! store radii 
   
!  attenuate shrink point so that query point is just outside shrunk EOA
   x_shrink(1:nx) = shrink*query%xs(1:nx) + (1.d0-shrink)*eoa%geom(1:nx)

   call ell_pt_shrink( nx, eoa%geom(1:nx), eoa%geom(nx+1:nx+ng), x_shrink(1:nx), &
	                      1, info )
	                       
   call sell_ell_rad_set( eoa )           !  update ELL radii 
   table%stats(16:17) = table%stats(16:17) + (eoa%geom(ngeom-1:ngeom) - r_eoa)
   if( eoa%geom(ngeom) > table%stats(15) ) table%stats(15) = eoa%geom(ngeom)

   pell => table%speoa%ell_pt(id)%ell  ! reform PEOA
   call ell_aff_pr( nx, eoa%geom(1:nx), eoa%geom(nx+1:nx+ng), na, &
	                 table%ua(:,nx+1), table%ua(:,1:na),  &
	                 pell%geom(1:na), pell%geom(na+1:na+nga) )

   call sell_ell_rad_set( pell ) 

   call ebt_update( table%peoa_ebt, table%peoa_ebt%be_pt(id)%be, &
                    table%pair_cover, table%cut_upd )
      
   table%stats(54) = table%stats(54) + 1.d0  !  increment shrink counter
              
return 
end subroutine isat_ecc

!======================= GROWING =================================+================

subroutine isat_grow( query ) 

!  Attempt to grow EOAs and shrink EOIs based on grow-point xs.
!  On entry, query%[xs,xa,fs] are set.
!  On return, query%id is the ID of a leaf which has been inclusively grown,
!  or query%id <= 0 if there is none.

!  Note: if a leaf is identified for degeneration, the degeneration is delayed 
!        until the search no longer references the EOI.
   
   type (query_type), pointer :: query   ! query

   type (table_type), pointer :: table   !  table 

   type (be_type), pointer :: be
   integer    :: f_tests, a_tests, be_tests, f_total, a_total, be_total, &
                 n_update, max_be_tests, max_f_tests, n_eoi, id, id_degen, &
                 gr_max, n_grows, mode_eoi
			     
   logical    :: grown, inclusive, degenerate
   real(k_xf) :: cpu_q, cpu_f, cpu_gr, cpu_be, n_q, n_f, n_gr, n_be, err_rel, grow_frac
   real       :: cpu_1, cpu_2
   
   query%id = 0
   table    => query%table
   if( table%leaves < 1 ) return ! quick return if no leaves

   n_eoi    =  table%seoi%n_ell
   if( n_eoi < 1 ) return  !  quick return if no EOIs

   nullify(be)
   max_be_tests = 3*n_eoi  ! unlimited BE testing
   max_f_tests  = 3*n_eoi  ! unlimited EOI testing
   
   call cpu_time( cpu_1 )
   table%stats(97) = table%stats(97) + (cpu_1 - query%cpu_start)  !  query time to here 

! determine maximum allowed number of EOI grows (gr_max)

   cpu_q  = table%stats(82)             !  CPU for queries
   cpu_f  = table%stats(91)             !  CPU for f evals
   cpu_gr = sum( table%stats(98:100) )  !  CPU for growing 
   cpu_be = table%stats(98)             !  CPU for searching during growing
   
   n_q    = table%stats(1)   !  number of queries
   n_f    = table%stats(77)  !  number of f evals
   n_gr   = table%stats(44)  !  number of EOA/EOI grows
   n_be   = table%stats(52)  !  number of BE tests during growing
   
   grow_frac    = table%grow_frac
   mode_eoi     = table%mode_eoi

   if( query%slow_prog ) then  ! special treatment for slow progress
      grow_frac      = table%grow_frac_slow
      table%mode_eoi = 3       !  supress EOI shrinking (re-set below)
   endif

   if( grow_frac > 0.d0  .and.  n_gr > 1.d0  .and. n_be > 1.d0 .and. &
       cpu_gr > 0.1d0  .and.  cpu_f > 0.1d0  .and.  cpu_be > 0.1d0 ) then
!     limit so that CPU time is no more than grow_frac * cpu_f

      cpu_f  = cpu_f  / n_f     !  CPU per f eval
	  cpu_gr = cpu_gr / n_gr    !  CPU per EOA/EOI grow
	  cpu_be = cpu_be / n_be    !  CPU per BE test

      gr_max       = nint( grow_frac * cpu_f / cpu_gr )
      max_be_tests = nint( grow_frac * cpu_f / cpu_be )
	  
      if( table%grow_max > 0.d0 ) then
         gr_max       = min( gr_max,       table%grow_max )
         max_be_tests = min( max_be_tests, table%grow_max )
      endif
      
      gr_max = max( gr_max, table%grow_min )
   else
      if( table%grow_max > 0 ) then
         gr_max = table%grow_max  ! set to grow_max
	  else
	     gr_max = n_eoi+1    ! set unlimited
	  endif
   endif

   table%stats(50) = gr_max

   n_grows  = 0
   f_total  = 0
   a_total  = 0
   be_total = 0
   id_degen = -1

   do  !  loop over EOIs containing the grow point

      call ebt_query_pair( table%seoi, query%xs, table%peoi_ebt, query%xa, &
	                       table%gr_mode, max_f_tests, max_be_tests, be, &
						   f_tests, a_tests, be_tests )

      f_total  = f_total  + f_tests
      a_total  = a_total  + a_tests
      be_total = be_total + be_tests
      max_be_tests = max_be_tests - be_tests
      
      call cpu_time( cpu_2 )
      table%stats(98) = table%stats(98) + (cpu_2 - cpu_1) ! searching time
      cpu_1 = cpu_2
      
	  if( .not.associated(be) ) exit  !  all EOIs have been grown 

      id = be%id
	  call isat_leaf_grow( query, id, grown, inclusive, degenerate, err_rel )

      call cpu_time( cpu_1 )
      table%stats(99) = table%stats(99) + (cpu_1 - cpu_2) ! growing time
      
      if( inclusive  .and.  query%id <= 0 ) then
	     query%id = be%id
      endif

	  if( degenerate ) then
!  degenerate previously-identified leaf (id_degen) (if any):
!  then mark this leaf (id) for future degeneration
         if( id_degen > 0 ) call isat_leaf_degenerate( table, id_degen )
		 id_degen = id
	  endif

      table%stats(41) = table%stats(41) + 1.d0
      
	  if( grown ) then
	     n_grows = n_grows + 1
	     table%stats(44) = table%stats(44) + 1.d0
	  endif

!  test for max grows and max BE tests
      if( n_grows >= gr_max  .or.  max_be_tests <= 0 ) then  
	     table%stats(51) = table%stats(51) + 1.d0  !  incomplete grow event
	     exit
	  endif

   end do

   if( id_degen > 0 ) call isat_leaf_degenerate( table, id_degen )

   call ebt_mark_update( table%peoa_ebt, table%pair_cover, table%cut_upd, n_update )

   call cpu_time( cpu_2 )
   table%stats(100) = table%stats(100) + (cpu_2 - cpu_1) ! updating time
   
   table%stats(49) = table%stats(49) + n_update

   table%stats(40) = table%stats(40) + 1.d0
   table%stats(42) = table%stats(42) + f_total
   table%stats(43) = table%stats(43) + a_total
   table%stats(47) = be_total / (2.d0 * n_eoi - 1.d0) ! fractioon of BEs tested
   table%stats(52) = table%stats(52) + be_total

   if( query%id > 0 ) then
      table%stats(45) = table%stats(45) + 1.d0  !  inclusive
   elseif( n_grows > 0 ) then 
      table%stats(46) = table%stats(46) + 1.d0  !  exclusive
!  else - no EOA/EOI modifications: don't count as either
   endif

   table%mode_eoi = mode_eoi  !  restore original value
   return
end subroutine isat_grow  !--------------------------------------------------------

subroutine isat_leaf_degenerate( table, id )

!  Degenerate leaf with ID=id in table.
!  PEOI is removed from EBT on PEOIs.
!  PEOI and EOI are destroyed.

   type (table_type), pointer :: table
   integer, intent(in)        :: id

   integer :: ngeom

   if( id < 1  .or.  .not.associated( table%seoi%ell_pt(id)%ell ) ) &
      call isat_abort('isat_leaf_degenerate',1, mess='bad leaf ', isv=id )

   ngeom = table%seoi%ngeom 

   table%leaf_pt(id)%leaf%props(6) = table%stats(1) ! query on which leaf degenerated
   table%stats(14)    = table%stats(14) + 1.d0      ! number of degenerated leaves
   table%stats(18:19) = table%stats(18:19) - &
                        table%seoi%ell_pt(id)%ell%geom(ngeom-1:ngeom) ! sum of radii

   call ebt_remove( table%peoi_ebt%be_pt(id)%be, table%pair_cover, 0 )
   call sell_ell_destroy( table%speoi, id )
   call sell_ell_destroy( table%seoi,  id )

   return
end subroutine isat_leaf_degenerate  !---------------------------------------------

subroutine isat_leaf_grow( query, id, grown, inclusive, degenerate, err_rel )

!  Given the grow point xs (which is in the EOI of the leaf with ID=id)
!  at which f(xs)=fs, "grow" the leaf.  
!  That is, the EOA is possibly grown, and the EOI is possibly shrunk.
!  The grow is inclusive if xs is covered by the grown EOA.
!  The leaf may be marked for degenerated (i.e., the EOI is to be destroyed).
!  Note that, if retrieving is incomplete, xs may be in the EOA.

   type (query_type), pointer :: query      !  query%[xs,fs] are set on entry
   integer, intent(in)        :: id         !  ID of leaf
   logical, intent(out)       :: grown      ! .true. if leaf is modified
   logical, intent(out)       :: inclusive  ! .true. if final EOA covers xs
   logical, intent(out)       :: degenerate ! .true. if leaf is to be degenerated
   real(k_xf), intent(out)    :: err_rel    !  error, relative to tolerance

!  Parameters controlling operations

!  mode_con - action to be taken in the event of a conflict (i.e., if the
!             shrunk EOI does not cover the grown EOA
!           = 1 - restore original EOA; degenerate leaf
!           = 2 - shrink EOA to be covered by EOI; degenerate leaf
!           = 3 - shrink EOA to be covered by EOI; do not degenerate leaf
!           = 4 - do nothing (conflict not resolved)

!  mode_eoi - method used to define shrink point, in shrinking the EOI
!           = 0 - conservative shrinking
!           = 1 - shrink based on xs if inaccurate
!           = 2 - shrink based on err ~ r**2
!           = 3 - no shrinking

!  mode_shrink - type of ellipsoid shrinking to be used
!           = 1 - maximun-volume algorithm       (E_v)
!           = 2 - maximum near content algorithm (E_n)
!           = 3 - conservative algorithm         (E_c)

   type (table_type), pointer :: table      !  table to which leaf belongs
   type(leaf_type), pointer   :: leaf
   type(ell_type),  pointer   :: eoa, eoi, pell

! The grow point  xs  is extended by the factor  extend  so that a repeat query
! with the same value of  xs  results in a retrieve.
   real(k_xf), parameter      :: extend = 1.d0 + 1.d-5

   integer    :: mode_con, mode_eoi, mode_shrink, nx, na, ng, nga, ngeom,  &
                 n_shrink, info, n_mark
   real(k_xf) :: eoa_geom(query%table%seoa%ngeom), dxs(query%table%nx), &
                 x0(query%table%nx), xse(query%table%nx), xs_fac, &
				 fl(query%table%nf), errsq, r_eoa(2), r_eoi(2)
   logical    :: shrink_eoi, grow_eoa, in, accurate, conflict

   grown      = .false.
   inclusive  = .false.
   degenerate = .false.
   conflict   = .false.
   err_rel    = huge(0.d0)

   table => query%table

   if( table%icheck > 0 ) then  !  check input
      if( .not.associated(table) ) then
         call isat_abort('isat_leaf_grow', 1, mess= 'table not associated' )

      elseif( .not.associated(table%leaf_pt(id)%leaf) ) then
         call isat_abort('isat_leaf_grow', 2, mess= 'LEAF not associated, ID = ', &
                         isv=id )

      elseif( .not.associated(table%seoa%ell_pt(id)%ell) ) then
         call isat_abort('isat_leaf_grow', 3, mess='EOA not associated, ID = ', &
                         isv=id )
	  endif  
   endif

! if ELL is degenerate (i.e., without an EOI), no growth is possible
   if( .not.associated(table%seoi%ell_pt(id)%ell) ) return  

   mode_con    = table%mode_con   !  set local variables
   mode_eoi    = table%mode_eoi
   mode_shrink = table%mode_shrink

   nx    = table%nx 
   na    = table%na 
   ng    = table%ng 
   nga   = table%nga 
   ngeom = table%seoa%ngeom 

   leaf  => table%leaf_pt(id)%leaf  !  set local pointers
   eoa   => table%seoa%ell_pt(id)%ell 
   eoi   => table%seoi%ell_pt(id)%ell

   if( table%icheck > 0 ) then  !  check that  xs  is in EOI
      call sell_ell_pt_in( nx, ngeom, eoi%geom(1:ngeom), query%xs, in ) 
      if( .not.in ) &
         call isat_abort('isat_leaf_grow', 4, mess = 'xs not covered by EOI' )
   endif

!  form extended grow point

   x0   = eoa%geom(1:nx)    !  leaf center
   dxs  = query%xs - x0     !  displacement of xs from leaf center
   xse  = x0 + dxs * extend !  extended grow point

!  return if xse is in EOA: this can occurs if retrieving is incomplete
   call sell_ell_pt_in( nx, ngeom, eoa%geom, xse(1:nx), in )  
   if( in ) then
      inclusive =.true.
	  return  
   endif

!  make copy of initial EOA in case re-setting is required
   if( mode_con == 1 ) eoa_geom(1:ngeom) = eoa%geom(1:ngeom)  
   r_eoa = eoa%geom(ngeom-1:ngeom)  
   r_eoi = eoi%geom(ngeom-1:ngeom)  

!----- evaluate error, and grow EOA if accurate

!  form linear approximation to fs
   fl  = leaf%xfh(table%nxp1:table%jfe) + matmul( leaf%g, dxs(1:nx) ) 

   call isat_eoa_err( leaf, fl, query%fs, accurate, errsq ) ! evaluate error

   if( accurate ) then  ! accurate case

! check that error in piece-wise constant approximation is not too large
      fl = leaf%xfh(table%nxp1:table%jfe)  - query%fs
	  if( sum( fl*fl ) > table%etol_consq ) return

!  grow EOA
      call ell_pt_modify( nx, x0, eoa%geom(nx+1:nx+ng), xse(1:nx) )  !  grow EOA
      grow_eoa = .true.                    ! indicate EOA grown
	  leaf%props(3) = table%stats(1)       ! query of last EOA grow
	  leaf%props(8) = leaf%props(8) + 1.d0 ! number of EOA grows

   else  !  inaccurate
      grow_eoa = .false.   ! indicate EOA not grown
   endif

!----- start of shrinking EOI

   err_rel    = sqrt( errsq / leaf%etolsq )  !  error relative to tolerance
   xs_fac     = 0.d0     ! shrink point is xs_fac * xs  (xs_fac=0 for mode_eoi=3)

   if( err_rel > 1.d0 ) then  ! inaccurate case  

      if(     mode_eoi == 0 ) then
	     xs_fac = max( 0.25d0 , 1.d0 / err_rel**0.25d0 )
      elseif( mode_eoi == 1 ) then
	     xs_fac = 1.d0
      elseif( mode_eoi == 2 ) then
	     xs_fac = 1.d0 / sqrt( err_rel )
      endif

   else  !  accurate case

      if(     mode_eoi == 0 ) then
	     if( err_rel > 0.25d0 ) xs_fac = 1.d0 / err_rel
      elseif( mode_eoi == 1 ) then
	     xs_fac = 0.d0
      elseif( mode_eoi == 2 ) then
	     xs_fac = 1.d0 / max( sqrt( err_rel ), 1.d-30 )
      endif

   endif

   shrink_eoi = .false.  !  indicate not yet shrunk

   if( xs_fac > 0.d0 ) then   !  shrink EOI
      xse = x0 + xs_fac * dxs !  shrink point
      call ell_pt_shrink( nx, x0(1:nx), eoi%geom(nx+1:nx+ng), xse(1:nx), &
	                      mode_shrink, info )

!  Note that info >0 indicates that the shrink point is not covered by EOI
	  if( info == 0 ) then  
	     shrink_eoi    = .true.               ! indicate EOI shrunk
	     leaf%props(4) = table%stats(1)       ! query of last EOI shrink
	     leaf%props(9) = leaf%props(9) + 1.d0 ! number of EOI shrinks
	  endif
   endif

!----------  identify and treat conflict (EOA not covered by EOI)
   n_shrink = 0  
   if( grow_eoa  .or. shrink_eoi ) then  ! possible conflict
      grown = .true.

      if( mode_con /=4 ) then  ! detect conflict, and shrink EOA (if conflict)
         call ell_pair_shrink( nx, eoi%geom(nx+1:nx+ng), eoa%geom(nx+1:nx+ng), &
	                           n_shrink )

      else  !  mode_con==4       detect conflict if degen_c > 0; do not shrink EOA
	     eoa_geom(1:ngeom) = eoa%geom(1:ngeom)
	     if( table%degen_c > 0 ) then  
            call ell_pair_shrink( nx, eoi%geom(nx+1:nx+ng), eoa_geom(nx+1:nx+ng), &
	                              n_shrink )
	     else
	        n_shrink = 0
	     endif
      endif
   endif

   if( n_shrink > 0 ) then  !  there is a conflict
      conflict = .true.
	  leaf%props(5)  = table%stats(1)        ! query of last EOA/EOI conflict
	  leaf%props(10) = leaf%props(10) + 1.d0 ! number of EOA/EOI conflicts
   endif

!  determine whether leaf is to be degenerated
   if( conflict  .and.  mode_con < 3 ) degenerate = .true.

   if( (table%degen_r > 0  .and.  leaf%props(7)  >= table%degen_r)  .or.  &
       (table%degen_g > 0  .and.  leaf%props(8)  >= table%degen_g)  .or.  &
       (table%degen_s > 0  .and.  leaf%props(9)  >= table%degen_s)  .or.  &
       (table%degen_c > 0  .and.  leaf%props(10) >= table%degen_c)   )    &
	      degenerate = .true.

!  reset PEOI and EOI radii unless EOI has not changed, or will be degenerated
   if( shrink_eoi  .and.  .not.degenerate ) then
       call sell_ell_rad_set( eoi )
	   table%stats(18:19) = table%stats(18:19) + (eoi%geom(ngeom-1:ngeom) - r_eoi)

	   pell => table%speoi%ell_pt(id)%ell  ! reform PEOI
	   call ell_aff_pr( nx, eoi%geom(1:nx), eoi%geom(nx+1:nx+ng), na, &
	                    table%ua(:,nx+1), table%ua(:,1:na),  &
	                    pell%geom(1:na), pell%geom(na+1:na+nga) )

       call sell_ell_rad_set( pell )

	   call ebt_mark_set( table%peoi_ebt, table%peoi_ebt%be_pt(id)%be, n_mark )
       table%stats(48) = table%stats(48) + n_mark
   endif

!  reset EOA geometry
   if( n_shrink > 0  .and.  mode_con == 1 ) then
      eoa%geom(1:ngeom) = eoa_geom(1:ngeom)  !  restore original EOA

   elseif( grow_eoa  .or.  (n_shrink >0  .and.  mode_con /= 4) ) then
!  reset PEOA and EOA radii unless EOA has not changed

	  call sell_ell_rad_set( eoa )           !  update ELL radii 
	  table%stats(16:17) = table%stats(16:17) + (eoa%geom(ngeom-1:ngeom) - r_eoa)
	  if( eoa%geom(ngeom) > table%stats(15) ) table%stats(15) = eoa%geom(ngeom)

	  
	  pell => table%speoa%ell_pt(id)%ell  ! reform PEOA
	  call ell_aff_pr( nx, eoa%geom(1:nx), eoa%geom(nx+1:nx+ng), na, &
	                    table%ua(:,nx+1), table%ua(:,1:na),  &
	                    pell%geom(1:na), pell%geom(na+1:na+nga) )

      call sell_ell_rad_set( pell ) 

	  call ebt_mark_set( table%peoa_ebt, table%peoa_ebt%be_pt(id)%be, n_mark )
      table%stats(48) = table%stats(48) + n_mark
   endif

!-----------  Determine type of grow:  inclusive = .false. set above

   if( grow_eoa ) then
      if( n_shrink == 0  .or.  mode_con == 4 ) then
         inclusive = .true.  !  EOA has not been shrunk to resolve conflict

      elseif( mode_con > 1 ) then
         call sell_ell_pt_in( nx, ngeom, eoa%geom, query%xs(1:nx), in )
		 inclusive = in
      endif
   endif
   
   return
end subroutine isat_leaf_grow   !--------------------------------------------------

subroutine isat_eoa_err( leaf, fs1, fs2, accurate, errsq )

!  Measure error between two scaled values of f (fs1 and fs2)
!  relative to error tolerance of leaf.

   type (leaf_type), pointer  :: leaf	    	! leaf from which etol obtained
   real(k_xf), intent(in)     :: fs1(:), fs2(:)	! values of f to be tested
   logical   , intent(out)    :: accurate   	! true if err < etol
   real(k_xf), intent(out)    :: errsq         	! square of error

   real(k_xf)                 :: df(size(fs1))

   df    = fs1 - fs2
   errsq = dot_product( df, df )

   if(  errsq < leaf%etolsq ) then
      accurate = .true.
   else
      accurate = .false.
   endif

   leaf%table%stats(73) = leaf%etolsq

   return
end subroutine isat_eoa_err

!=================== ADDING =======================================================

subroutine isat_add( table, xs, fs, ga, ha, modret )

!  add leaf to table

   type (table_type), pointer :: table    ! pointer to table  
   real(k_xf), intent(in)     :: xs(:)    ! position of query point (scaled)
   real(k_xf), intent(in)     :: fs(:)    ! f(x) (scaled)
   real(k_xf), intent(in)     :: ga(:,:)  ! g(x)
   real(k_xf), intent(in)     :: ha(:)    ! h(x)
   integer,    intent(out)    :: modret   ! =1 for add, =2 for delete/add

   type (leaf_type), pointer  :: newleaf
   integer  :: id

!  if table is full, delete last leaf in MFU list before adding new leaf

   if( table%full ) then
      modret = 2
	  id   = table%eoa_mfu%last  !  least frequently used leaf
	  call isat_leaf_delete( id, table )

   else  !  add to not-full table
      modret = 1
   endif

!  add leaf

   call isat_leaf_set( table, newleaf, xs, fs, ga, ha )
   call isat_leaf_add( table, newleaf, xs )  

   table%leaves = table%leaves + 1	! one more leaf in table

!  update maxleaves
   if( table%stomby > 0.d0 ) then
      call isat_leaves( table%stomby, table%nx, table%nf, table%nh, table%na,  &
                        table%leaves, nint(table%stats(14)), table%maxleaves )
      table%stats(13) = table%maxleaves
   endif
                     
   if( table%leaves >= table%maxleaves ) table%full = .true.

   return
end subroutine isat_add  !---------------------------------------------------------

subroutine isat_leaves( stomby, nx, nf, nh, na, leaves, degen, max_leaves )

   real(k_xf), intent(in)  :: stomby ! storage in megabytes
   integer,    intent(in)  :: nx     ! dimension of full space 
   integer,    intent(in)  :: nf     ! dimension of f 
   integer,    intent(in)  :: nh     ! dimension of h
   integer,    intent(in)  :: na     ! dimension of affine space 
   integer,    intent(in)  :: leaves ! current number of leaves (degenerate and non-degenerate)
   integer,    intent(in)  :: degen  ! current number of degenerate leaves
   integer,    intent(out) :: max_leaves  !  maximum number of leaves
   
   integer :: weoa, wpeoa, weoi, wpeoi, words, max_words
     
! 8-byte words per item  
   weoa  = (nx*(nx+1))/2 + (nx*nf*k_g)/k_dp + 3*nx + nf + nh + 30 ! leaf and ELL
   wpeoa = na*(na+1) + 3*na + 16  ! PEOA and BEs
   weoi  = (nx*(nx+1))/2 + nx +2  ! ELL
   wpeoi = wpeoa
   
   words      = (weoa+wpeoa)*(leaves-degen) + (weoi+wpeoi)*degen !  words current used
   max_words  = stomby*2**(20-4)     ! max number of 16-bit words
   max_leaves = leaves + (max_words - words) / (weoa+wpeoa+weoi+wpeoi) 
   
   return
end subroutine isat_leaves   !-----------------------------------------------------
  
subroutine isat_leaf_delete( id, table )

! delete leaf of ID=id from table

   integer, intent(in)         :: id       ! ID
   type (table_type), pointer  :: table	   ! table

   type (leaf_type), pointer   :: leaf     ! leaf
   integer    :: nri, nro
   real(k_xf) :: rad(2)
   logical    :: degenerate

!  determine whether leaf is degenerate
   if( associated(table%peoi_ebt%be_pt(id)%be) ) then
      degenerate = .false.
   else
      degenerate = .true.
   endif

!  subtract radii from sum
   nri = table%nx+table%ng+1  ! index of r_in
   nro = nri+1                ! index of r_out

   rad = table%seoa%ell_pt(id)%ell%geom(nri:nro)  ! EOA radii
   table%stats(16:17) = table%stats(16:17) - rad

   if( .not.degenerate ) then
      rad = table%seoi%ell_pt(id)%ell%geom(nri:nro)  ! EOI radii
      table%stats(18:19) = table%stats(18:19) - rad
   endif
  
   leaf => table%leaf_pt(id)%leaf

!  remove EOA from all data structures

   call bt_remove( table%eoa_bt, id )
   call ll_remove( table%eoa_mru, id )
   call ll_remove( table%eoa_mfu, id )
   call ebt_remove( table%peoa_ebt%be_pt(id)%be, table%pair_cover, table%cut_rem )

!  destroy EOA ELLs

   call sell_ell_destroy( table%seoa,  id )
   call sell_ell_destroy( table%speoa, id )

!  treat EOI similarly if it exists
   if( .not.degenerate ) then
      call ebt_remove( table%peoi_ebt%be_pt(id)%be, table%pair_cover, table%cut_rem )
	  call sell_ell_destroy( table%seoi,  id )
      call sell_ell_destroy( table%speoi, id )
   else
      table%stats(14) = table%stats(14) - 1.d0  !  one less degenerate leaf
   endif

   call id_return( table%idlist, id )  !  return ID

!  remove leaf

   call isat_leaf_kill( leaf )
   table%leaves    = table%leaves - 1
   table%stats(12) = table%leaves

   return
end subroutine isat_leaf_delete   !------------------------------------------------

subroutine isat_leaf_set( table, leaf, xs, fs, ga, ha )

!  allocate leaf (by call to isat_leaf_init) in table and set values
!  of xs, fs, ga, ha

   type (table_type), pointer :: table  ! table to which leaf belongs
   type (leaf_type),  pointer :: leaf   ! leaf to be set
   real(k_xf), intent(in)     :: xs(:)  ! (scaled)
   real(k_xf), intent(in)     :: fs(:)  ! (scaled)
   real(k_xf), intent(in)     :: ga(:,:), ha(:)
   integer :: nx, nf, nh

   call isat_leaf_init( table, leaf )  !  allocate leaf

!  set values

   leaf%etolsq = ( table%etola + table%etolr * sqrt( dot_product(fs,fs) ) )**2

   nx = table%nx
   nf = table%nf
   nh = table%nh

   leaf%xfh(1:nx)                             = xs
   leaf%xfh(table%nxp1:table%jfe)             = fs
   leaf%g                                     = ga * table%g2gs
   if( nh > 0 ) leaf%xfh(table%jhs:table%jhe) = ha

   leaf%props    = 0.d0
   leaf%props(1) = table%stats(1)  ! query on which leaf added

   return
end subroutine isat_leaf_set !--------------------------------------------------------
 
subroutine isat_leaf_init( table, leaf )

!  allocate leaf in table

   type (table_type), pointer :: table
   type (leaf_type),  pointer :: leaf 

   integer    :: nx, nf

   nx = table%nx
   nf = table%nf

!  set or nullify pointers

   allocate( leaf )
   leaf%table => table

!  allocate storage

   allocate( leaf%xfh(table%jhe) )
   allocate( leaf%g(nf,nx) )

   return
end subroutine isat_leaf_init  !----------------------------------------------------
 
subroutine isat_leaf_add( table, leaf, xs )  

! Add new leaf and associated ELLs to table and all datastructures.
! The geometry of the ELLs is set via isat_leaf_geom_init

   type (table_type), pointer  :: table	   ! table
   type (leaf_type), pointer   :: leaf     ! leaf
   real(k_xf), intent(in)      :: xs(:)    ! position of query point (scaled)

   real(k_xf) :: g_eoa(table%ng),   g_eoi(table%ng), rad(2)
   real(k_xf) :: g_peoa(table%nga), g_peoi(table%nga), pxs(table%na)

   integer :: id, nx, nf, ng, na, nri, nro

   nx = table%nx
   nf = table%nf
   ng = table%ng  ! Cholesky dimension
   na = table%na  ! dimension of affine space
   nri = nx+ng+1  ! index of r_in
   nro = nx+ng+2  ! index of r_out

   call id_get( table%idlist, id )  !  get ID

   if( id > table%n_pt ) call leaf_pt_double( table ) 

   leaf%id = id
   table%leaf_pt(id)%leaf => leaf

! set geometry of EOA and EOI
   call isat_leaf_geom_init( leaf, nx, nf, ng, xs, g_eoa, g_eoi ) 

! project EOA and EOI onto affine space (centers pxs are the same)
   call ell_aff_pr( nx, xs, g_eoa, na, table%ua(:,nx+1), table%ua(:,1:na), &
                    pxs, g_peoa )
   call ell_aff_pr( nx, xs, g_eoi, na, table%ua(:,nx+1), table%ua(:,1:na), &
                    pxs, g_peoi )

! create all ELLs
   call sell_ell_create( table%seoa,  id )
   call sell_ell_create( table%speoa, id )
   call sell_ell_create( table%seoi,  id )
   call sell_ell_create( table%speoi, id )

! set all Ells' geometry
   call sell_ell_geom_set( table%seoa%ell_pt(id)%ell,  xs,  g_eoa )
   call sell_ell_geom_set( table%seoi%ell_pt(id)%ell,  xs,  g_eoi )
   call sell_ell_geom_set( table%speoa%ell_pt(id)%ell, pxs, g_peoa )
   call sell_ell_geom_set( table%speoi%ell_pt(id)%ell, pxs, g_peoi )

   rad = table%seoa%ell_pt(id)%ell%geom(nri:nro)  ! EOA radii
   table%stats(55:56) = rad
   table%stats(16:17) = table%stats(16:17) + rad
   if( rad(2) > table%stats(15) ) table%stats(15) = rad(2)

   rad = table%seoi%ell_pt(id)%ell%geom(nri:nro)  ! EOI radii
   table%stats(57:58) = rad
   table%stats(18:19) = table%stats(18:19) + rad

! add ELLs to datastructures
   call bt_add( table%eoa_bt, id )

   call ll_add( table%eoa_mru, id, 1 )
   call ll_add( table%eoa_mfu, id, 2 )

   call ebt_add( table%peoa_ebt, id, .true., table%pair_cover, table%cut_add )
   call ebt_add( table%peoi_ebt, id, .true., table%pair_cover, table%cut_add )

   return
end subroutine isat_leaf_add   !------------------------------------------------------

subroutine isat_leaf_geom_init( leaf, nx, nf, ng, xs, g_eoa, g_eoi )    

!  Initialize the EOA and EOI of leaf.

! Specification of EOA

!  The (scaled) gradient matrix has the SVD:   g = df/dx = u sv vt.
!  The radius of the sphere of accuracy in f-space is  etol.
!  The inverse of g maps this sphere in f-space into an ellipsoid in x-space.
!  (For nf < nx, this ellipsoid is infinite in nx-nf directions.)
!  This corresponds to the region of accuracy for the constant approximation.
!  In the i-th principal direction in x-space (for i <= nf), the half-length
!  of the principal axis is:  r(i) = etol / sv(i)   where sv(i) is the i-th
!  singular value of g.  However, if sv(i) is very small, r(i) is very large,
!  and the quadratic term in the expansion for error may be significant.
!  Hence, the specified half-length is limited to be no greater than r0_eoa,
!  defines as:  r0_eoa = eoa_lim * sqrt( etol ).

! Specification of EOI

!  Provisionally, the EOI is set to the ball B of radius 
!  r0_eoi = eoi_lim * sqrt( etol ).  This is taken to be the EOI if
!  table%ad_look is less than one, or table%ad_theta is non-positive.
!  Otherwise, up to table%ad_look existing tabulation points are found
!  within B; and up to table%ad_use of these are used to construct the
!  EOI so that it does not cover these points.  This is done by
!  ell_pts_uncover which uses the parameter table%ad_theta, and returns
!  the parameter phi, which is the factor by which the intermediate ellipsoid 
!  E_n is shrunk to form the EOI.  (If table%ad_shrink=1, this shrinking is 
!  reversed.)  If necessary, the EOA is then shrunk so that it is covered by 
!  the EOI.

   type (leaf_type), pointer  :: leaf
   integer, intent(in)        :: nx, nf, ng
   real(k_xf), intent(in)     :: xs(:)    ! position of query point (scaled)
   real(k_xf), intent(out)    :: g_eoa(ng), g_eoi(ng)


   type (table_type), pointer :: table
   real(k_xf)   :: gc(nf,nx), sv(min(nx,nf)), uu(1,1), vT(nx,nx), &
                   etol, radius, lam(nx), v(nx,nx), radinv, r0_eoa, svmin
   real(k_xf)   :: work(10*(nx+nf))
   real(k_xf)   :: r0_eoi, phi, dsq_lim, r_eoa_max, r_eoi_min
   integer      :: INFO, i, j, npts, npts_use, cp_tests, in_tests, &
                   jpass, n_shrink, ad_look

   table => leaf%table

!--------- set EOA ----------------------------------------------------------

   gc     = leaf%g  !  form SVD of g = uu * sv * vT:  only sv and vT are used 

   call dgesvd( 'N', 'A', nf, nx, gc, nf, sv, uu, 1, vT, nx, &
	                 work, 10*(nx+nf), INFO )
  
   if( INFO /= 0 ) call isat_abort( 'isat_leaf_init', 1, &
                                      mess = 'svd failed, INFO=', isv = INFO )
                                      
   etol   = sqrt( leaf%etolsq )           ! error tolerance for leaf
   r0_eoa = table%eoa_lim * sqrt( etol )  ! upper limit on EOA semi-axes
   svmin  = etol / r0_eoa                 ! corresponding value of sv

! set principal semi-axes of EOA
   do i = 1, nx
      if( i <= nf ) then
		 radius = etol / max( sv(i), svmin ) 
      else
         radius = r0_eoa
      endif
      lam(i)  = 1. / radius
   end do

   r_eoa_max = radius  ! largest principal semi-axis of EOA

!  EOA matrix is A = L * L^T = V * lam^2 * V^T
!  Set  v = (vT)^T

   v = transpose( vT )
   call ell_eig2chol( nx, v, lam, g_eoa )

!  The EOA may be modified below if it not covered by the EOI.

!--------- set EOI ----------------------------------------------------------

   table%stats(59) = 0.d0
   table%stats(63) = 0.d0

   r0_eoi = table%eoi_lim * sqrt(etol)  ! upper limit on EOI semi-axes

!  Provisionally, set EOI to ball of radius r0_eoi

!  A = I / r0_eoi^2;   L = I / r0_eoi

   v      = 0.d0  !  form L (unpacked)
   radinv = 1.d0 / r0_eoi
   do i = 1, nx
      v(i,i) = radinv
   end do

! pack L into EOI chol
   call ell_full2low( nx, v, g_eoi )

!  Consider shrinking EOI based on existing tabulated points

   ad_look = min( table%ad_look, table%ad_use )  ! number of points to look for

   if( ad_look < 1  .or.  table%ad_theta <= 0.d0 ) return  !  do not shrink EOI

! Find up to ad_look points covered by initial EOI (ball of radius r0_eoi).
! npts points are found with IDs table%ids and square distances table%dsq from
!    the center of the EOI 

   call bt_ball_pts( table%eoa_bt, xs, r0_eoi, table%ad_look, npts, &
                     table%ids, table%dsq, cp_tests, in_tests ) 

   if( npts == 0 ) return !  no points to use

   if( npts <= table%ad_use ) then  !  use all points
      npts_use = npts
      dsq_lim  = maxval( table%dsq(1:npts) )  ! square distance to furthest point

   else  !  select the table%ad_use closest points 
      npts_use = table%ad_use
	  call select_pts( table%ad_use, npts, table%dsq, dsq_lim, jpass )
   endif

!  Copy the npts_use selected points into p.  Points not selected
!  are indicated by dsq(i) < 0.

   j = 0
   do i = 1, npts  
	  if( table%dsq(i) < 0.d0 ) cycle  ! indicates point not selected
	  j = j + 1
      table%p(:,j) = table%seoa%ell_pt( table%ids(i) )%ell%geom(1:nx)
	  if( j == npts_use ) exit
   end do

   call ell_pts_uncover( nx, xs, r0_eoi, table%ad_theta, npts_use, table%p,  &
                         g_eoi, phi, r_eoi_min )

   if( table%ad_shrink == 1 ) then  ! reverse shrinking of EOI
      g_eoi     = g_eoi * phi
	  r_eoi_min = r_eoi_min / phi 
   endif 

!---------  shrink EOA if necessary (not necesary if EOI has not been shrunk)

   if( r_eoa_max <= r_eoi_min ) then
      n_shrink = 0
   else
      call ell_pair_shrink( nx, g_eoi, g_eoa, n_shrink )
   endif

   table%stats(59) = sqrt( dsq_lim)
   table%stats(63) = phi
   table%stats(60) = table%stats(60) + npts_use
   table%stats(61) = table%stats(61) + cp_tests
   table%stats(62) = table%stats(62) + in_tests
   table%stats(64) = table%stats(64) + n_shrink

   return
end subroutine isat_leaf_geom_init  !------------------------------------------------------------

subroutine select_pts( nout, nin, x, x_lim, jpass )

!  Identify the nout (or more) smallest components of x(1:nin), x(i) >= 0.
!  Components not identified are set negative.
!  In non-pathalogical cases, and for sufficiently large npass (set below),
!  exactly (nin-nout) components are set negative.  In general, at most
!  (nin-nout) components are set negative

   integer,    intent(in)    :: nout   !  required number of components
   integer,    intent(in)    :: nin    !  number of components of x
   real(k_xf), intent(inout) :: x(nin) !  components, x(i) >= 0
   real(k_xf), intent(out)   :: x_lim  !  selected components 0 <= x(i) <= x_lim
   integer,    intent(out)   :: jpass  !  number of passes performed

   integer, parameter    :: npass = 3     !  maximum number of passes
   integer, parameter    :: nbins = 100   !  number of bins
   real(k_xf), parameter :: small = 1.d-6 !  small parameter to avoid rounding 
                                          !  problems

   integer    :: ipass, bin(nbins), i, j, jless
   real(k_xf) :: xlow, xhigh, dx, dxi

   jpass = 0
   xlow  = minval( x )
   xhigh = maxval( x )
   x_lim = xhigh

   if( nin <= nout ) return  !  all components identified

   passes: do ipass = 1, npass  !  loop over passes

      dx    = (xhigh-xlow)/(nbins-2)
	  if( dx <= 0.d0 ) return  !  all components coincident
	  jpass = ipass

      xlow  = xlow  - small*dx
      xhigh = xhigh + small*dx
      dx    = (xhigh-xlow)/(nbins-2)
      dxi   = 1.d0 / dx
      bin   = 0  

! x(i) is in bin j if  xlow+(j-2)*dx <= x(i) < xlow+(j-1)*dx
! if x(i) < xlow, then j=1;  if x(i) > xhigh, then j=nbins

	  do i = 1, nin  !  count entries in bins
	     j = 2 + dxi * ( x(i) - xlow )

		 if( j < 2 ) then
		    bin(1) = bin(1) + 1
		 elseif( j >= nbins ) then
		    bin(nbins) = bin(nbins) + 1
		 else
		    bin(j) = bin(j) + 1
		 endif
      end do

	  jless = 1        !  x_lim is in bin jless+1
	  do j = 2, nbins  !  cumulative sum of bin counts
         bin(j) = bin(j-1) + bin(j)
		 if( bin(j) < nout ) jless = j
	  end do

      xhigh = xlow + dx * jless   !  bounds on x_lim
	  xlow  = xhigh - dx

! x_lim found if bin(jless+1) has exactly one entry
	  if( bin(jless+1)-bin(jless) == 1  .or.  ipass == npass ) exit passes
      
   end do  passes

   x_lim = xhigh

   j = 0
   do i = 1, nin
      if( x(i) < x_lim ) then
		 j = j + 1
	  else
		 x(i) = -x(i)
	  endif
   end do

   if( j < nout ) call isat_abort('select_pts', 1, mess='select_pts, error: j, nout = ', &
                                  ivar=(/j,nout/) )
   return

end subroutine select_pts  !-------------------------------------------------------

subroutine leaf_pt_double( table ) 

!  double size of array table%leaf_pt

      type (table_type),   pointer :: table

      type (leaf_pointer), pointer :: pt_temp(:)
	  integer                      :: id

!  allocate temporary array pt_temp of twice size
   	  allocate( pt_temp(2*table%n_pt) )  

	  do id = 1, 2*table%n_pt
	     nullify( pt_temp(id)%leaf )
      end do

!  copy leaf_pt%pt to pt_temp
	  do id = 1, table%n_pt
	     if( associated( table%leaf_pt(id)%leaf ) ) pt_temp(id)%leaf => table%leaf_pt(id)%leaf
		 nullify( table%leaf_pt(id)%leaf )
	  end do

!  double size of table%ell_pt
	  deallocate( table%leaf_pt )

      table%n_pt = 2 * table%n_pt
	  allocate( table%leaf_pt(table%n_pt) )  

	  do id = 1, table%n_pt
	     
	     if( associated( pt_temp(id)%leaf ) ) then
		    table%leaf_pt(id)%leaf => pt_temp(id)%leaf
	     else
		    nullify( table%leaf_pt(id)%leaf )
		 endif

		 nullify(pt_temp(id)%leaf) 
	  end do

	  deallocate( pt_temp )

	  return

end subroutine leaf_pt_double   !-----------------------------------------------

!============  AFFINE SPACE ==========================================================

subroutine isat_affine( table, sell, nx, changed )  

!  Re-define affine space (if necessary).

!  The current affine space A has dimension na.
!  The new affine space (if redefined) has dimension na_new.
!  For table%aff_na < 0,  na_new = nx
!  For table%aff_na = 0,  na_new  is determined as described below
!  For table%aff_na > 0,  na_new = table%aff_na <= nx.

!  The covariance W of the existing ELL centers is formed.
!  The square roots of the normalized singular values sig(i) of W
!  are formed: sig(i) represent the normalized rms deviation in the i-th
!  principal direction.  The normalization is such that sig(1)=1.
!  The normalized rms deviation from the existing na-dimensional affine space,
!  sig_star, is evaluated.  This satisfies  sig_star >= sig(na+1).

!  For table%aff_na=0, the value of na_new is set such that sig(na_new+1)
!  is less than table%aff_thresh.  This value of na_new is then modified
!  if need be in order to satisfy  na_new >= na, and na_new <= table%na_max.

!  If na_new /= na, then the affine space is re-defined, and this is
!  indicated by  changed=/true.

!  If na_new == na, then the affine space is re-defined only if doing so
!  would decrease sig_star by at least a factor of table%aff_ratio, i.e.,
!  if sig_star > sig(na+1) * table%aff_ratio.

   type (table_type), pointer :: table    ! table
   type (sell_type),  pointer :: sell     ! SELL  - set of ellipsoids   
   integer, intent(in)        :: nx       ! dimension of x-space 
                                          ! (for dimension statements) 
   logical, intent(out)       :: changed  ! true if affine space redefined

   integer   :: i, lwork, info, na, na_new, nan
                                
   real(k_d) :: sig_star, r_in, &
                cm(nx), sv(nx), sig(nx), lam(nx), vt(1), &
                w(nx,nx), u(nx,nx), ww(nx,nx), &
				g( (nx*(nx+1))/2 ), gp( (nx*(nx+1))/2 )
			    

   real(k_d) :: work(10*nx*(nx+10))
				
   lwork = 10*nx*(nx+10)   !  ...for SVD

   changed = .false.   ! ...until found otherwise

!  return if no table or sell
   if( .not.associated(table)  .or.  .not.associated(sell) ) call isat_abort(  &
         'isat_affine',1, mess='table, sell or spell not associated' )
                                  
!  check that nx is consistent
   if( nx /= table%nx ) call isat_abort( 'isat_affine', 2,  &
      mess= 'nx, table%nx inconsistent ', ivar=(/nx, table%nx/) )

!  return if no leaves		
   if( sell%n_ell < 1 )  return 

!  return if  na=nx  is required and already set
   na = table%na
   if( table%aff_na < 0  .and.  na == nx ) return

!-----  form sig(i), the normalized rms deviations in the principal directions

!  form mean and covariance W of ELL centers
   call sell_cov( sell, cm, w )

!  form SVD of covariance W = u * sv * vt
   ww = w
   call dgesvd( 'A', 'N', nx, nx, ww, nx, sv, u, nx, vt, 1, work, lwork, info )

   if( info /= 0 ) call isat_abort( 'isat_affine', 3,  &
                                     mess= 'SVD failed, info = ', isv=info)

!  form normalized rms sig(i);  
!  and detrmine na_new: the value of na such that sig(na+1) < table%thresh
   na_new = nx

   do i = 1, nx
      sig(i) = sqrt( sv(i) / sv(1) )
	  lam(i) = 1.0d0 / max( sig(i), table%aff_thresh * 1.d-6 )
	  if( na_new == nx  .and.  sig(i) < table%aff_thresh ) na_new = i-1
   end do

!----- form sig_star, the normalized deviation from the existing affine space

   sig_star = 0.d0
   if( na < nx ) then

!     form ellipsoid E = { x | x^T sv(1)*W^{-1} x <= 1 }  
!     corresponding to W (normalized)

      call ell_eig2chol( nx, u, lam, g )

!  project E onto the orthogonal complement of the affine space to form Ep

      lam = 0.d0
      nan = nx - na
      call ell_aff_pr( nx, lam, g, nan, lam, table%ua(:,na+1:nx), sv, gp )

!  sig_star is the outer radius of Ep

      call ell_radii( nan, gp, r_in, sig_star )
   endif
   
!-----------  specify na_new

   if( table%aff_na < 0 ) then  !  fix na_new = nx
	  na_new = nx

   elseif( table%aff_na > 0 ) then  !  fix na_new = aff_na
      na_new = min( table%aff_na, nx-1 )

   else  !  impose limits on na_new
      na_new = max( na, min( na_new, table%na_max ) )
   endif
   
!----------  determine whether affine space is to be re-defined

   if( na_new /= na ) then  !  yes, if na changed
      changed = .true.
      
   elseif( na_new < nx ) then  !  yes, if sig_star can be reduced by a factor of at least aff_ratio
      if( sig_star >= sig(na+1) * table%aff_ratio ) changed = .true.
   endif

!---------  re-define affine space
   if( changed ) then  
      na               = na_new
      table%na         = na
	  table%nga        = (na*(na+1))/2
	  table%ua(:,1:nx) = u
	  table%ua(:,nx+1) = cm

	  table%stats(68:70) = 0.d0 

      if( table%write_log ) write(table%lu_log,'(a,2i9,1p,3e13.4)') &
	      'New affine space: n_ell, na, sig_star, sig(na:na+1) = ',  &
		  sell%n_ell, na, sig_star, sig(na), sig(min(na+1,nx))
   endif

   table%stats(65) = na * 1.d0
   table%stats(66) = sig_star
   table%stats(67) = sig(min(na+1,table%nx)) 
   
   return
end subroutine isat_affine  !---------------------------------------------------------------

subroutine isat_na_max( aff_lim, nx, na_max )

!  Set the upper limit  na_max  on the allowed dimension of the affine space.

  integer, intent(in)  :: aff_lim  ! parameter controlling na_max -- see if-statement below
  integer, intent(in)  :: nx       ! dimension of full space
  integer, intent(out) :: na_max   ! upper limit on dimension of the affine space
  
  if( aff_lim > 0 ) then
     na_max = min( aff_lim, nx )
     
  elseif( aff_lim == -1 ) then
     na_max = nx
     
  elseif( aff_lim == -2 ) then
     na_max = max( nint(sqrt(1.d0*nx)), 1 )
     
  elseif( aff_lim == -3 ) then
     na_max = max( nint(log(1.d0*nx)/log(2.d0)), 1 )
     
  else
     call isat_abort( 'isat_na_max', 1, mess= 'invalid value of aff_lim ', &
                      ivar=(/aff_lim, nx, na_max/)  )
  endif

end subroutine isat_na_max  !---------------------------------------------------------------

subroutine isat_aff_rebuild( table, sell, spell, ebt )  !  rebuild SPELL

!  Given a set of ellipsoids (SELL) and an affine space, generate the set
!  of ellipsoids (SPELL) obtained by projecting SELL onto the affine space.
!  Then build an EBT for SPELL.  
!  Existing data structures (EBT, SPELL) are first destroyed.

   type (table_type), pointer :: table
   type (sell_type),  pointer :: sell   !  SELL (unchanged)
   type (sell_type),  pointer :: spell  !  SPELL (corresponding to SELL)
   type (ebt_type),   pointer :: ebt    !  EBT on SPELL

   type (ell_type),  pointer  :: ell
   integer   :: nx, ng, na, n_pell_0, n_pell, id, i, ids_assigned, id_max
   integer   :: id_array(table%n_pt)
   real(k_d) :: c(table%na), g(table%nga)
   
   if( .not.associated(table)  .or.  .not.associated(sell)  .or.  &
       .not.associated(spell)  .or.  .not.associated(ebt) )  &
          call isat_abort( 'isat_aff_rebuild', 1, mess=  &
                           'table, sell, spell or ebt not associated' )

   n_pell_0 = spell%n_ell      !  initial number of PELLs
   if( n_pell_0 == 0 ) return  !  nothing to do

   nx = table%nx
   ng = table%ng
   na = table%na

!  Destroy EBT and SPELL 
   call ebt_purge( ebt )
   call ebt_destroy( ebt )
   call sell_destroy( spell )

!  Re-initialize  EBT and SPELL
   call sell_create( spell, na, check=table%icheck, idlist=table%idlist )
   call ebt_initialize( spell, ebt, check=table%icheck, idlist=table%idlist,  &
                        pair_sep_qual=table%pair_sep )

!  loop over extant ELLs in SELL to form and add PELLs to SPELL
   n_pell = 0
   do id  = 1, sell%n_pt
      if( .not.associated(sell%ell_pt(id)%ell) ) cycle

      n_pell = n_pell + 1
      ell    =>  sell%ell_pt(id)%ell

!  form geometry (c,g) of PELL by projecting ELL
      call ell_aff_pr( nx, ell%geom(1:nx), ell%geom(nx+1:nx+ng), &
	                   na, table%ua(:,nx+1), table%ua(:,1:na), c, g )

	  call sell_ell_create( spell, id )
	  ell => spell%ell_pt(id)%ell
	  call sell_ell_geom_set( ell, c, g )

   end do

   if( n_pell /= n_pell_0 ) call isat_abort( 'isat_aff_rebuild', 1, mess=  &
                    'bad count of PELLs ', ivar=(/n_pell_0, n_pell/) )

!  Build EBT on SPELL by adding PELLs at random

   call id_rand_array( table%idlist, id_array )
   call id_stats( table%idlist, ids_assigned, id_max )

   do i = 1, ids_assigned  !  add PELLs; do not update BEs
      id = id_array(i)
      if( .not.associated(sell%ell_pt(id)%ell) ) cycle
      call ebt_add( ebt, id, .false., table%pair_cover, -1 )
   end do

   do id  = 1, sell%n_pt  !  make all nodes for updating
      if( associated( ebt%be_pt(id)%be ) ) &
	     call ebt_mark_set( ebt, ebt%be_pt(id)%be, i )
   end do

   call ebt_mark_update( ebt, table%pair_cover, 0, i )  ! update

   return
end subroutine isat_aff_rebuild   !--------------------------------------------

!============== TEST ROUTINES ===============================================

subroutine isat_integrity( table, lu, info )

!  Check the integrity of the data structures in table.

type (table_type), pointer :: table
integer, intent(in)        :: lu   ! logical unit for error messages
!                                    < 0 to suppress error messages
integer, intent(out)       :: info ! =0 for consistent datastructures 
!                                    /0 indicates inconsistency

   integer :: id_status, ids_assigned, id_max, id, ierr, n_eoa, n_eoi
   logical :: op

   if( lu >= 0 ) then  !  output required?
      op = .true.
   else
      op = .false.
   endif

    if( op ) write(lu,*)'isat_integrity: starting...'
!  check association of all data structures

    info = 0
	if( .not.associated( table%idlist ) )   info = 1
    if( .not.associated( table%seoa ) )     info = 2
    if( .not.associated( table%speoa ) )    info = 3
    if( .not.associated( table%seoi ) )     info = 4
    if( .not.associated( table%speoi ) )    info = 5
    if( .not.associated( table%eoa_bt ) )   info = 6
    if( .not.associated( table%eoa_mru ) )  info = 7
    if( .not.associated( table%eoa_mfu ) )  info = 8
    if( .not.associated( table%peoa_ebt ) ) info = 9
    if( .not.associated( table%peoi_ebt ) ) info = 10

	if( info /= 0 ) then
	   if( op ) write(lu,*)'structure not associated ', info
	   return
    endif

!  check id list
	call id_query( table%idlist, -1, id_status, ids_assigned, id_max )

	if( id_status /= 0 ) then
	   info = id_status
	   if( op ) write(lu,*)'idlist inconsistent ', info
	   return
    endif

!  check all leaves and pointers
    n_eoa = 0
	n_eoi = 0

	do id = 1, table%n_pt
       call id_query( table%idlist, id, id_status, ids_assigned, id_max )

	   if( id_status < 1  .or.  id_status > 4 ) then
	      info = 11
	      if( op ) write(lu,*)'id, (invalid) id_status = ', id, id_status

	   elseif( id_status == 2  .or.  id_status == 3 ) then ! not currently assigned
	      ierr = 0
	      if( associated( table%leaf_pt(id)%leaf ) )      ierr = 12
	      if( associated( table%seoa%ell_pt(id)%ell ) )   ierr = 13
	      if( associated( table%seoi%ell_pt(id)%ell ) )   ierr = 14
	      if( associated( table%speoa%ell_pt(id)%ell ) )  ierr = 15
	      if( associated( table%speoi%ell_pt(id)%ell ) )  ierr = 16

	      if( associated( table%eoa_bt%nd_pt(id)%nd ) )   ierr = 17
	      if( associated( table%peoa_ebt%be_pt(id)%be ) ) ierr = 18
	      if( associated( table%peoi_ebt%be_pt(id)%be ) ) ierr = 19

          info = max( info, ierr )
		  if( ierr /= 0  .and. op )  &
             write(lu,*)'association problem A, id, info = ', id, info

		else  !  currently assigned
          n_eoa = n_eoa + 1
		  ierr = 0
		  if( table%leaf_pt(id)%leaf%id /= id )            ierr = 21
		  if( table%seoa%ell_pt(id)%ell%ell_id /= id )     ierr = 22
		  if( table%speoa%ell_pt(id)%ell%ell_id /= id )    ierr = 23

		  if( associated( table%seoi%ell_pt(id)%ell ) ) then  !  EOI exists
		     n_eoi = n_eoi + 1
		     if( table%seoi%ell_pt(id)%ell%ell_id /= id )  ierr = 25
		     if( table%speoi%ell_pt(id)%ell%ell_id /= id ) ierr = 26
		     if( table%leaf_pt(id)%leaf%props(6) /= 0 )    ierr = 27
		  else  !  degenerate
		     if( associated( table%seoi%ell_pt(id)%ell ) ) ierr = 28
		     if( associated( table%speoi%ell_pt(id)%ell) ) ierr = 29
		     if( table%leaf_pt(id)%leaf%props(6) <= 0 )    ierr = 30
		  endif

          info = max( info, ierr )
          if( ierr /= 0 .and. op )  &
             write(lu,*)'association problem B, id, info = ', id, info
      endif
   end do

!  check counts

   ierr = 0
   if( table%leaves /= n_eoa )      ierr = 31
   if( table%seoa%n_ell /= n_eoa )  ierr = 32
   if( table%speoa%n_ell /= n_eoa ) ierr = 33
   if( table%seoi%n_ell /= n_eoi )  ierr = 34
   if( table%speoi%n_ell /= n_eoi ) ierr = 35
   if( ids_assigned /= n_eoa )      ierr = 36

   info = max( info, ierr )
   if( ierr /= 0 .and. op )  &
      write(lu,*)'count problem, ierr = ', ierr

!  check EBTs: PEOA_EBT

   do id = 1, table%n_pt
      call id_query( table%idlist, id, id_status, ids_assigned, id_max )
      
	  ierr = 0
	  if( id_status == 1 ) then !  EOA exists
         if( .not.associated(table%peoa_ebt%be_pt(id)%be) ) then
		    ierr = 41
		 else  !  EOA exists and associated
		    if( table%peoa_ebt%be_pt(id)%be%id /= id ) ierr = 42
		    if( .not.associated( table%peoa_ebt%be_pt(id)%be%geom, &
			                 table%speoa%ell_pt(id)%ell%geom ) )ierr = 43
		 endif		  
	  else  !  no EOA
         if( associated(table%peoa_ebt%be_pt(id)%be) ) ierr = 44
	  endif

   info = max( info, ierr )
   if( ierr /= 0 .and. op )  &
        write(lu,*)'PEOA_EBT problem, id, info = ', id, info

   end do

!  check EBTs: PEOI_EBT

   do id = 1, table%n_pt
      call id_query( table%idlist, id, id_status, ids_assigned, id_max )
      
	  ierr = 0
	  if( id_status == 1  .and.  table%leaf_pt(id)%leaf%props(6) == 0 ) then !  EOI exists
         if( .not.associated(table%peoi_ebt%be_pt(id)%be) ) then
		    ierr = 51
		 else  !  EOI exists and associated
		    if( table%peoi_ebt%be_pt(id)%be%id /= id ) ierr = 52
		    if( .not.associated( table%peoi_ebt%be_pt(id)%be%geom, &
			                 table%speoi%ell_pt(id)%ell%geom ) )ierr = 53
		 endif		  
	   else  !  no EOI
         if( associated(table%peoi_ebt%be_pt(id)%be) ) ierr = 54
	   endif

       info = max( info, ierr )
       if( ierr /= 0 .and. op )  &
          write(lu,*)'PEOI_EBT problem, id, info = ', id, info

   end do

   if( op ) write(lu,*)'isat_integrity: completed, info = ', info
   
return
end subroutine isat_integrity  !---------------------------------------------------

subroutine isat_xf_range( table, nx, nf, lu )

!  Determine the range of leaf values of x and f in table.

type (table_type), pointer :: table
integer, intent(in)        :: nx   ! dimension of x
integer, intent(in)        :: nf   ! dimension of f
integer, intent(in)        :: lu   ! logical unit for error messages
!                                    < 0 to suppress error messages

type(leaf_type), pointer :: leaf
integer    :: ids_assigned, id_max, id, id_count, i, j, nsamp, id_array(table%n_pt),  &
              lwork, info, rank, iwrk1(1)
real(k_dp) :: xf_min(nx+nf), xf_max(nx+nf), xf_mean(nx+nf), g(nf,nx), f_min(nf), &
              f_max(nf), fp(nf)
real(k_dp), allocatable :: a(:,:), b(:,:), sv(:), work(:)
integer,    allocatable :: iwork(:)

   logical :: op

   if( lu >= 0 ) then  !  output required?
      op = .true.
   else
      op = .false.
   endif
   
   if( table%leaves < 1 ) return  !  return if no leaves
   
   if( nx /= table%nx  .or.  nf /= table%nf ) &  !  check consistency
      call isat_abort( 'isat_xf_range:', 1, mess= 'inconsistent nx, nf ', &
                        ivar=(/nx, table%nx, nf, table%nf/)  )
   
   call id_stats( table%idlist, ids_assigned, id_max )
   
   xf_min  =  huge(0.d0)  !  find max and min of x and f
   xf_max  = -huge(0.d0)
   xf_mean = 0.d0
   
   id_count = 0
   do id = 1, id_max
      if( .not.associated(table%leaf_pt(id)%leaf) ) cycle
      id_count = id_count + 1
      leaf => table%leaf_pt(id)%leaf
      
      do j = 1, nx+nf
         xf_min(j)  = min( xf_min(j), leaf%xfh(j) )
         xf_max(j)  = max( xf_max(j), leaf%xfh(j) )
         xf_mean(j) = xf_mean(j) + leaf%xfh(j)
      end do
   end do
   
   xf_mean = xf_mean / id_count
   
!  form LMSE estimate fh(x) of f(x)-mean(f)
!  least squares solution to: A Y = B, where:
!    rows of A are samples of  x - x_mean
!    rows of B are samples of  f - f_mean
!    Y(j,k) = df_k / dx_j

if( table%leaves < nx ) then  !  too few samples
   f_min = 0.d0
   f_max = 0.d0
else

   nsamp = min( 1000, table%leaves )
   
   call id_rand_array( table%idlist, id_array )
   allocate( a(nsamp,nx) )
   allocate( b(nsamp,nf) ) 
   allocate( sv(nsamp+nx+nf) )

   do i = 1, nsamp
      id   = id_array(i)
      leaf => table%leaf_pt(id)%leaf
      a(i,:) = leaf%xfh(1:nx) - xf_mean(1:nx)  !  rows of A are x vectors
      b(i,:) = leaf%xfh(nx+1:nx+nf) - xf_mean(nx+1:nx+nf)  !  rows of B are f vectors
   end do
   
   lwork = -1  !  pre-call to determone lwork
   !  XXX SBP minor bug fix 1/8/09
   call dgelsd( nsamp, nx, nf, a, nsamp, b, nsamp, sv, 1.d-9, rank, fp, lwork, iwrk1, info )
   
   lwork = 2*fp(1)
   allocate( work(lwork) )
   allocate( iwork(lwork) )
   
   call dgelsd( nsamp, nx, nf, a, nsamp, b, nsamp, sv, 1.d-9, rank, work, lwork, iwork, info )
  
   if( info /= 0 ) call isat_abort( 'isat_xf_range:', 3, &
                                    mess= 'dgelsd info = ', isv = info )

!  Find range of deviation of f from  LMSE
   g = transpose( b(1:nx,1:nf) )
   f_min  =  huge(0.d0)  
   f_max  = -huge(0.d0)

   do i = 1, nsamp
      id   = id_array(i)
      leaf => table%leaf_pt(id)%leaf
      fp   = leaf%xfh(nx+1:nx+nf)-xf_mean(nx+1:nx+nf)  &
             - matmul( g, leaf%xfh(1:nx)-xf_mean(1:nx) )
      do j = 1, nf
         f_min(j) = min( f_min(j), fp(j) )
         f_max(j) = max( f_max(j), fp(j) )
      end do
   end do
endif
      
   if( op ) then
      write(lu,*)' '
      write(lu,'(a)')'     j    x_min(j)    x_max(j)    x_diff(j)   x_mean(j)'
      do j = 1, nx
         write(lu,'(i6,1p,9e12.2)') j, xf_min(j), xf_max(j), xf_max(j)-xf_min(j), xf_mean(j)
      end do 
      
      write(lu,*)'  '
      write(lu,'(a)')'     j    f_min(j)    f_max(j)    f_diff(j)   fp_diff(j)    f_mean(j)'
      do j = 1, nf
         write(lu,'(i6,1p,9e12.2)') j, xf_min(nx+j), xf_max(nx+j),  &
                xf_max(nx+j)-xf_min(nx+j), f_max(j)-f_min(j), xf_mean(nx+j)
      end do     
   endif
   
!laniu 
    if(allocated(a))   deallocate(a)
    if(allocated(b))   deallocate(b)
    if(allocated(sv))  deallocate(sv)
    if(allocated(work)) deallocate(work)
    if(allocated(iwork)) deallocate(iwork)
       

   return

end subroutine isat_xf_range  !------------------------------------------------------

subroutine isat_table_plot( table )

!  write out all EOAs and EOIs for plotting

   implicit none
   type (table_type), pointer :: table  !  table to which EOA belongs

   type (ell_type),  pointer  :: eoa, eoi

   integer :: id, id_status, ids_assigned, id_max, nx, ng

   call id_query( table%idlist, 1, id_status, ids_assigned, id_max )

   if( id_status < 0 ) call isat_abort( 'isat_table_plot', 1, mess = &
                                    'table%idlist not associated' )

   nx = table%nx
   ng = table%ng

   do id = 1, id_max
	  eoa => table%seoa%ell_pt(id)%ell
	  if( .not.associated(table%seoi%ell_pt(id)%ell) ) then
!  plot degenerate EOA blue
         call isat_leaf_plot( nx, eoa%geom(1:nx), eoa%geom(nx+1:nx+ng), 3 )
      else
!  plot EOA green and EOI red
	     eoi => table%seoi%ell_pt(id)%ell
         call isat_leaf_plot( nx, eoa%geom(1:nx), eoa%geom(nx+1:nx+ng), 2 )
         call isat_leaf_plot( nx, eoi%geom(1:nx), eoi%geom(nx+1:nx+ng), 1 )
      endif

   end do

   return
end subroutine isat_table_plot  !-----------------------------------------------------------

subroutine isat_leaf_plot( nx, c, gg, k )

!  write out ELL for plotting

   integer, intent(in)   :: nx
   real(k_d), intent(in) :: c(nx)                !  center
   real(k_d), intent(in) :: gg( (nx*(nx+1))/2 )  !  packed Cholesky
   integer, intent(in)   :: k  ! color

   real(k_d) :: u(nx,nx), lam(nx)

   if( nx /= 2 ) call isat_abort( 'isat_leaf_plot', 1, mess = 'called with nx = ', isv=nx )

   call ell_chol2eig( nx, gg, u, lam )

   write(41,'(1p,2e20.10)') 4.d0, 1.d0*k, c, lam, u

   return
end subroutine isat_leaf_plot  !-----------------------------------------------------------

end module isat_subs
