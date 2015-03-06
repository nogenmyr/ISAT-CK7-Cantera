!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module sell_ll

!  Module for searching sets of ellipsoids using linked lists (LLs).

!  S.B. Pope  1/15/2006

!  Data structures
!     ll_type - for linked list LL

!  Subroutines
!   ll_initialize       initialize empty LL, including allocating ll
!   ll_destroy          destroy LL, including deallocating ll
!   ll_add              add ELL with ID=id at start or end of list
!   ll_remove           remove ELL with ID=id from LL
!   ll_mru_update       update ELL in most-recently-used LL
!   ll_mfu_update       update ELL in most-frequently-used LL
!   ll_mfu_cum          determine number of ELLs accounting for fraction of uses
!   ll_mfu_check        check consistency of MFU list
!   ll_set_tested       set entries as tested or untested
!   ll_query            attempt to find ELL which covers point x
!   ll_query_list       attempt to find list of ELLs which covers point x
!   ll_query_pair       attempt to find ELL which covers point x based on projected ELLs
!   ll_query_pair_list  attempt to find list of ELLs which covers point x based on projected ELLs
!   ll_write            write LL for checkpointing
!   ll_read             read LL
!   ll_param_set        set parameters
!   ll_param_get        get parameters

! ll_initialize( ll, check, n_pt, idlist )
! ll_destroy( ll )
! ll_add( ll, id, place )
! ll_remove( ll, id )
! ll_mru_update( ll, id )
! ll_mfu_update( ll, sell, id )
! ll_mfu_cum( ll, sell, used, k )
! ll_mfu_check( ll, sell, status )
! ll_set_tested( ll, sell, nset, truth )
! ll_query( ll, sell, x, max_tests, id0, id, tests )
! ll_query_list( ll, sell, x, max_tests, id0, n_id_max, n_id, id, tests )
! ll_query_pair( ll, sell_a, sell_b, xa, xb, max_a_tests, max_b_tests, id0, &
!                id, a_tests, b_tests )
! ll_query_pair_list( ll_a, sell_a, sell_b, xa, xb, max_a_tests, max_b_tests, &
!                     id0, n_id_max, n_id, id, a_tests, b_tests  )
! ll_write( ll, lu )
! ll_read( ll, lu, idlist )
! ll_param_set( ll, check )
! ll_param_get( ll, check )


!**********  consistency tests should be removed once testing complete  *****************

use sell_m
implicit none

type :: ll_type
   integer          :: first, last          !  id of first and last ell in linked list
   integer, pointer :: before(:), after(:)  !  id of ell before and after this ell
   integer          :: n_pt                 !  dimension of pointers
   type (id_list_type), pointer :: idlist   !  data structure used in id_list to generate IDs

! Optional:  set through ll_initialize or ll_param_set
   integer :: check                !  level of checking to be performed

end type ll_type

contains  !===================================================================

subroutine ll_initialize( ll, check, n_pt, idlist )

!  Allocate and initialize (empty) linked list ll.  Optionally set check.

   type (ll_type), pointer                :: ll      ! (out) data structure for linked list.

!  Optional
   type (id_list_type), pointer, optional :: idlist  ! (in) data structure used in id_list to generate IDs
   integer, intent(in), optional          :: check   ! level of checking  =0 minimal; =1 some; =2 maximal
   integer, intent(in), optional          :: n_pt    ! initial size of pointer array

   allocate( ll )
   ll%first = 0 ! null
   ll%last  = 0 ! null

   if( present(idlist) ) then
      ll%idlist => idlist
   else
      nullify( ll%idlist)
   endif

   if( present(check) ) then
      call ll_param_set( ll, check )
   else
      ll%check = 1
   endif

   if( present(n_pt) ) then
      if( n_pt < 1 ) then
	     write(0,*)'ll_initialize: invalid n_pt = ', n_pt
		 stop
	  endif
      ll%n_pt  = n_pt
   else
      ll%n_pt = n_pt_0
   endif

   allocate( ll%before(ll%n_pt) )
   allocate( ll%after(ll%n_pt) )

   ll%before(1:ll%n_pt) = 0  ! null
   ll%after( 1:ll%n_pt) = 0  ! null

   return
end subroutine ll_initialize  !------------------------------------------------------------

subroutine ll_destroy( ll )

!  Destroy and deallocate linked list ll.

   type (ll_type), pointer :: ll    ! (in) data structure for linked list.

   if( .not.associated(ll) ) then
      write(0,*) 'll_destroy: ll not associated; stopping'
	  stop
   endif

   deallocate( ll%before )
   deallocate( ll%after )
   deallocate( ll )

   return
end subroutine ll_destroy  !------------------------------------------------------------

subroutine ll_add( ll, id, place )

!  Add ELL with ID=id to linked list ll at start (for place = 1) or at end (for place = 2)

   type (ll_type), pointer :: ll     !  data structure for linked list. 
   integer, intent(in)     :: id     !  ID of ell
   integer, intent(in)     :: place  !  =1 to add at start,  =2 to add at end

   integer                 :: n_pt, id_status, ids_assigned, id_max 
   integer, allocatable    :: temp(:)

   if( .not.associated(ll) ) then
      write(0,*) 'll_add: ll not associated; stopping'
	  stop
   endif

   if( ll%check > 0 ) then
      if( id < 1 ) then  !  check input
         write(0,*) 'll_add: invalid id = ', id
	     stop
      endif

      if( ll%check > 1  .and.  associated( ll%idlist ) ) then
	     call id_query( ll%idlist, id, id_status, ids_assigned, id_max )
		 if( id_status /= 1 ) then
		    write(0,*) 'll_add: id has not been assigned, id = ', id
	        stop
	     endif
	  endif
   endif

   if( id > ll%n_pt ) then  !  need to increase size of arrays
      n_pt = max( id, 2*ll%n_pt )
	  allocate( temp(n_pt) )

	  temp(1:ll%n_pt)      = ll%before
	  temp(ll%n_pt+1:n_pt) = 0
	  deallocate( ll%before )
	  allocate( ll%before(n_pt) )
	  ll%before(1:n_pt) = temp(1:n_pt)

	  temp(1:ll%n_pt)      = ll%after
	  deallocate( ll%after )
	  allocate( ll%after(n_pt) )
	  ll%after(1:n_pt) = temp(1:n_pt)

      deallocate( temp )
	  ll%n_pt = n_pt
   endif

!  add ELL

   if( ll%first < 1 ) then  !  first ell to be added to empty list
      ll%first = id
      ll%last  = id
	  ll%before(id) = 0 ! null
	  ll%after( id) = 0 ! null
	  return
   endif

   if( place == 1 ) then  !  add to start
      ll%before(id) = 0  !  null
	  ll%after( id) = ll%first
	  ll%before(ll%first) = id
	  ll%first = id

   elseif( place == 2 ) then  !  add to end
      ll%after( id) = 0  !  null
	  ll%before(id) = ll%last
	  ll%after(ll%last) = id
	  ll%last = id

   else
      write(0,*) 'll_add: invalid place = ', place
	  stop
   endif

   return
end subroutine ll_add  !------------------------------------------------------------------

subroutine ll_remove( ll, id )

!  Remove ELL with ID=id from linked list ll

   type (ll_type), pointer :: ll  !  data structure for linked list. 
   integer, intent(in)     :: id  !  ID of ell to be removed

   integer                 :: id_status, ids_assigned, id_max 

   if( .not.associated(ll) ) then
      write(0,*) 'll_remove: ll not associated; stopping'
	  stop
   endif

   if( ll%check > 0 ) then
      if( id < 1 ) then  !  check input
         write(0,*) 'll_remove: invalid id = ', id
	     stop
      endif

      if( id > ll%n_pt ) then  
         write(0,*) 'll_remove: invalid id, n_pt = ', id, ll%n_pt
	     stop
      endif

      if( ll%check > 1  .and.  associated( ll%idlist ) ) then
	     call id_query( ll%idlist, id, id_status, ids_assigned, id_max )
		 if( id_status /= 1 ) then
		    write(0,*) 'll_remove: id has not been assigned, id = ', id
	        stop
	     endif
	  endif
   endif

   if( ll%before(id) > 0 ) then    !  ell is not first
      if( ll%after(id) > 0 ) then  !  ell is not last (nor first)
	     ll%before( ll%after(id) ) = ll%before(id)
		 ll%after(  ll%before(id)) = ll%after(id)
	  else                         !  ell is last (but not first)
	     ll%after(  ll%before(id)) = 0
		 ll%last                   = ll%before(id)
	  endif
   else                            !  ell is first
      if( ll%after(id) > 0 ) then  !  ell is not last (but is first)
	     ll%before( ll%after(id) ) = 0
		 ll%first                  = ll%after(id)
	  else                         !  ell is last (and first)
	     ll%first = 0
		 ll%last  = 0
	  endif
   endif

   ll%before(id) = 0
   ll%after(id)  = 0

   return
end subroutine ll_remove  !------------------------------------------------------------------

subroutine ll_mru_update( ll, id )

!  Update most-recently-used linked list by moving ELL to the front

   type (ll_type), pointer :: ll  !  data structure for linked list. 
   integer, intent(in)     :: id  !  ID of ELL

   integer                 :: id_status, ids_assigned, id_max 

   if( .not.associated(ll) ) then
      write(0,*) 'll_mru_update: ll not associated; stopping'
	  stop
   endif

   if( ll%check > 0 ) then
      if( id < 1 ) then  !  check input
         write(0,*) 'll_mru_update: invalid id = ', id
	     stop
      endif

      if( ll%check > 1  .and.  associated( ll%idlist ) ) then
	     call id_query( ll%idlist, id, id_status, ids_assigned, id_max )
		 if( id_status /= 1 ) then
		    write(0,*) 'll_mru_update: id has not been assigned, id = ', id
	        stop
	     endif
	  endif
   endif

   if( ll%before(id) < 1 ) return    !  already at front
   call ll_remove( ll, id )          !  remove from current position
   call ll_add( ll, id, 1 )          !  add to beginning

   return

end subroutine ll_mru_update  !------------------------------------------------------------------

subroutine ll_mfu_update( ll, sell, id )

!  Update most-frequently-used linked list by promoting ELL as necessary

   type (ll_type), pointer   :: ll   !  data structure for linked list. 
   type (sell_type), pointer :: sell !  SELL to which ELLs belong
   integer, intent(in)       :: id   !  ID of ELL to be promoted

   integer                   :: id_mf, id_lf, id_status, ids_assigned, id_max 
   real(k_d)                 :: used

   if( .not.associated(ll) ) then
      write(0,*) 'll_mfu_update: ll not associated; stopping'
	  stop
   endif

   if( ll%check > 0 ) then
      if( id < 1 ) then  !  check input
         write(0,*) 'll_mfu_update: invalid id = ', id
	     stop
      endif

      if( ll%check > 1  .and.  associated( ll%idlist ) ) then
	     call id_query( ll%idlist, id, id_status, ids_assigned, id_max )
		 if( id_status /= 1 ) then
		    write(0,*) 'll_mfu_update: id has not been assigned, id = ', id
	        stop
	     endif
	  endif
   endif

!  traverse list from id towards front until more-frequently-used entry (id_mf) is encountered

   used  = sell%ell_pt(id)%ell%used
   id_mf = id

   do
      if( ll%before(id_mf) < 1 ) then  !  ELL to be moved to front of list
	     call ll_mru_update( ll, id )
		 return
      endif

	  id_mf = ll%before(id_mf)
	  if( used <= sell%ell_pt( id_mf )%ell%used ) exit  ! id_mf identified
   end do

   if( id_mf == ll%before(id) ) return    !  no re-ordering needed

   call ll_remove( ll, id )          !  remove from current position

! add id after id_mf

   id_lf            = ll%after(id_mf)
   ll%after(id_mf)  = id
   ll%before(id)    = id_mf
   ll%after(id)     = id_lf
   ll%before(id_lf) = id

   return

end subroutine ll_mfu_update  !------------------------------------------------------------------

subroutine ll_mfu_cum( ll, sell, used, k )

!  Return the smallest index k such that the cumulative number of uses of leaves
!  in the MFU linked list up to and including k is at least equal to  used.

   type (ll_type), pointer   :: ll     !  data structure for linked list. 
   type (sell_type), pointer :: sell   !  SELL to which ELLs belong
   real(k_d), intent(in)     :: used   !  specified
   integer, intent(out)      :: k      !  index in linked list

   integer                   :: id
   real(k_d)                 :: cum

   if( .not.associated(ll) ) then
      write(0,*) 'll_mfu_cum: ll not associated; stopping'
	  stop
   endif

   id  = ll%first
   k   = 0
   cum = 0.d0

   do
      k   = k + 1
	  cum = cum + sell%ell_pt(id)%ell%used
	  if( cum >= used  .or.   ll%after(id) < 1 ) then
		 return
	  endif
	  id = ll%after(id)
   end do

end subroutine ll_mfu_cum  !------------------------------------------------------------------

subroutine ll_mfu_check( ll, sell, status )

!  Check that MFU list is consistent

   type (ll_type), pointer   :: ll     !  data structure for linked list. 
   type (sell_type), pointer :: sell   !  SELL to which ELLs belong
   integer, intent(out)      :: status ! = 0 if consistent

   integer                   :: id

   if( .not.associated(ll) ) then
      write(0,*) 'll_mfu_check: ll not associated; stopping'
	  stop
   endif

   status = -1
   id = ll%first

   do
      if( ll%after(id) < 1 ) exit  !  traverse completed
	  if( sell%ell_pt(id)%ell%used < sell%ell_pt( ll%after(id) )%ell%used ) return  !  inconsistent
	  id = ll%after(id)
   end do

   status = 0

   return

end subroutine ll_mfu_check  !------------------------------------------------------------------

subroutine ll_set_tested( ll, sell, nset, truth )

!  Traverse nset entries in the linked list LL, setting  ell%tested = truth.

   type (ll_type), pointer   :: ll     ! LL
   type (sell_type), pointer :: sell   ! SELL
   integer,   intent(in)     :: nset   ! maximum number of ELL to be set
   logical,   intent(in)     :: truth  !  value of  tested  to be set

   type (ell_type), pointer  :: ell
   integer                   :: sets, id
 
   if( nset <= 0 ) return  !  quick return if no sets allowed

   if( .not.associated(ll) ) then
      write(0,*) 'll_set_tested: ll not associated; stopping'
	  stop
   endif

   if( ll%check > 0 ) then
      if( .not.associated(sell) ) then
         write(0,*) 'll_set_tested: sell not associated; stopping'
	     stop
      endif
   endif

   if( ll%first  < 1 ) return  !  empty linked list

   id   = ll%first
   sets = 0

!  loop until covering ellipsoid found or until end of list
   do
      ell => sell%ell_pt(id)%ell

	  ell%tested = truth
	  sets       = sets + 1

	  if( sets == nset ) exit   !   nset exceeded

	  if( ll%after(id) < 1 ) exit  !  at end of list
	  
	  id = ll%after(id)
   end do

   return
end subroutine ll_set_tested  !--------------------------------------------------------------

subroutine ll_query( ll, sell, x, max_tests, id0, id, tests )

!  Given a query point x, search the SELL, using the linked list LL, to attempt to find an
!  ELL (id) which covers x. Perform no more than max_tests tests in this attempt.
!  There may be 0, 1 or more ELLs satisfying the query.  
!  On entry, if id0 is less than 1, a fresh search is started, and the
!  ID of the first ELL (if any) found satisfying the query is returned in id.  If no such
!  ELL is found, then on return id is negative.  After a successful call (id>0), a
!  subsequent call with id0 set to the (positive) value of id from the previous call 
!  resumes the search.  Only those ELLs with ell%tested==.false. are tested, and then
!  ell%tested is set to .true..  The values of ell%tested should be set by an previous 
!  call to ll_set_tested.

   type (ll_type), pointer   :: ll         ! LL
   type (sell_type), pointer :: sell       ! SELL
   real(k_d), intent(in)     :: x(sell%nx) ! location of query point
   integer,   intent(in)     :: max_tests  ! maximum number of ELL tests
   integer,   intent(in)     :: id0        ! ID of starting ELL (or <1 for fresh search)
   integer,   intent(out)    :: id         ! ID of ELL satisfying query (or <0 for none)
   integer,   intent(out)    :: tests      ! number of ELL tests performed

   type (ell_type), pointer  :: ell
   integer                   :: id_status, ids_assigned, id_max
   logical                   :: in
 
   id    = -1 
   tests = 0
   if( max_tests <= 0 ) return  !  quick return if no tests allowed

   if( .not.associated(ll) ) then
      write(0,*) 'll_query: ll not associated; stopping'
	  stop
   endif

   if( ll%check > 0 ) then
      if( .not.associated(sell) ) then
         write(0,*) 'll_query: sell not associated; stopping'
	     stop
      endif

	  if( ll%check > 1  .and.  associated( ll%idlist )  .and.  id0 > 0 ) then
	     call id_query( ll%idlist, id0, id_status, ids_assigned, id_max )
		 if( id_status /= 1 ) then
		    write(0,*) 'll_query: id0 is not an assigned ID, id0 = ', id0
	        stop
	     endif
	  endif 
   endif

   if( ll%first  < 1 ) return  !  empty linked list

!  determine id at which to start
   if( id0 < 1 ) then  !  fresh start 
      id = ll%first

   elseif( id0 <= ll%n_pt ) then     ! potentially valid id0

      if( ll%after(id0) > 0 ) then
         id = ll%after(id0)          !  start from ll%after(id0)

	  elseif( id0 == ll%last ) then  ! already at end of list
	     id = -1
		 return

      else  !  inconsistency in linked list
	    write(0,*) 'll_query: inconsistency, id0, ll%last = ', id0, ll%last
		stop
      endif

   else  !  invalid id0
      write(0,*) 'll_query: invalid id0 > ll%n_pt ', id0, ll%n_pt
	  stop
   endif

!  loop until covering ellipsoid found or until end of list
   do
      ell => sell%ell_pt(id)%ell

	  if( ell%ell_id /= id ) then  !  check consistency   ***
	     write(0,*) 'll_query: inconsistency, ell%id /= id ', ell%ell_id, id
		 stop
      endif

      if( .not.ell%tested ) then
	     call sell_ell_pt_in( sell%nx, sell%ngeom, ell%geom(1:sell%ngeom), x, in )
	     tests      = tests + 1
		 ell%tested = .true.
	  else
	     in = .false.
      endif

	  if( in ) exit  !  covering ellipsoid found

	  if( tests == max_tests ) then
	     id = -1     !   max_tests exceeded
		 exit
	  endif

	  if( ll%after(id) < 1 ) then  !  at end of list
	  
	     if( id == ll%last ) then  ! check for consistency  ***
	        id = -1
		    exit
         else  !  inconsistency in linked list
	        write(0,*) 'll_query: inconsistency, id, ll%last = ', id, ll%last
		    stop
         endif
      endif

	  id = ll%after(id)
   end do

   nullify(ell)

   return
end subroutine ll_query  !--------------------------------------------------------------

subroutine ll_query_list( ll, sell, x, max_tests, id0, n_id_max, n_id, id, tests )

!  Given a query point x, search the SELL, using the linked list LL, to attempt to find 
!  up to n_id_max ELLs which covers x. Perform no more than max_tests tests in this attempt.
!  There may be 0, 1 or more ELLs found to satisfy the query.  
!  The number of covering ELLs found is returned in n_id and their IDs in id(:).
!  On entry, if id0 is less than 1, a fresh search is started.  After a successful call (n_id>0), a
!  subsequent call with id0 set to id(n_id) from the previous call resumes the search.

   type (ll_type), pointer   :: ll           ! LL
   type (sell_type), pointer :: sell         ! SELL
   real(k_d), intent(in)     :: x(sell%nx)   ! location of query point
   integer,   intent(in)     :: max_tests    ! maximum number of ELL tests
   integer,   intent(in)     :: id0          ! ID of starting ELL (or <1 for fresh search)
   integer,   intent(in)     :: n_id_max     ! maximum nimber of IDs to be returned
   integer,   intent(out)    :: n_id         ! number of covering ELLs found
   integer,   intent(out)    :: id(n_id_max) ! array of IDs of covering Ells
   integer,   intent(out)    :: tests        ! number of ELL tests performed

   type (ell_type), pointer  :: ell
   integer                   :: idc, id_status, ids_assigned, id_max
   logical                   :: in
 
   n_id  = 0
   id    = 0
   tests = 0

   if( max_tests <= 0 ) return  !  quick return if no tests allowed

   if( .not.associated(ll) ) then
      write(0,*) 'll_query_list: ll not associated; stopping'
	  stop
   endif

   if( ll%check > 0 ) then
      if( .not.associated(sell) ) then
         write(0,*) 'll_query_list: sell not associated; stopping'
	     stop
      endif

	  if( ll%check > 1  .and.  associated( ll%idlist )  .and.  id0 > 0 ) then
	     call id_query( ll%idlist, id0, id_status, ids_assigned, id_max )
		 if( id_status /= 1 ) then
		    write(0,*) 'll_query_list: id0 is not an assigned ID, id0 = ', id0
	        stop
	     endif
	  endif

	  if( n_id_max < 1 ) then
	     write(0,*) 'll_query_list: invalid n_id_max = ', n_id_max
	     stop
	  endif
   endif


   if( ll%first  < 1 ) return  !  empty linked list

!  determine id (idc) at which to start
   if( id0 < 1 ) then  !  fresh start 
      idc = ll%first

   elseif( id0 <= ll%n_pt ) then     ! potentially valid id0

      if( ll%after(id0) > 0 ) then
         idc = ll%after(id0)          !  start from ll%after(id0)

	  elseif( id0 == ll%last ) then  ! already at end of list
		 return

      else  !  inconsistency in linked list
	    write(0,*) 'll_query_list: inconsistency, id0, ll%last = ', id0, ll%last
		stop
      endif

   else  !  invalid id0
      write(0,*) 'll_query_list: invalid id0 > ll%n_pt ', id0, ll%n_pt
	  stop
   endif

!  loop until covering ellipsoid found or until end of list
   do
      ell => sell%ell_pt(idc)%ell

	  if( ell%ell_id /= idc ) then  !  check consistency  ***
	     write(0,*) 'll_query_list: inconsistency, ell%id /= idc ', ell%ell_id, idc
		 stop
      endif

	  call sell_ell_pt_in( sell%nx, sell%ngeom, ell%geom(1:sell%ngeom), x, in )
	  tests = tests + 1

	  if( in ) then  !  add ID to list
	     n_id     = n_id + 1
		 id(n_id) = idc
		 if( n_id == n_id_max ) exit !  n_id_max IDs found
	  endif

	  if( tests >= max_tests ) exit  !  max. tests reached
	  if( ll%after(idc) < 1 )  exit  !  at end of list
	  idc = ll%after(idc)
   end do

   nullify(ell)

   return
end subroutine ll_query_list  !--------------------------------------------------------------

subroutine ll_query_pair( ll, sell_a, sell_b, xa, xb, max_a_tests, max_b_tests, id0, &
                          id, a_tests, b_tests )

!  xb is a query point in the same space as the set of ellipsoids SELL_B.
!  xa is the projection of xb onto the (usually lower-dimensional) space containing the 
!  set of ellipsoids SELL_A, which are the projections of those in SELL_B.  
!  Thus, if xb is covered by an ellipsoid in SELL_B, then xa is covered by the corresponding
!  ellipsoid in SELL_A.

!  Given a query point xb, use the linked list LL, to attempt to find an
!  ELL in SELL_B which covers xb.  Tests are performed first to see if xa is covered by the
!  ellipsoid in SELL_A, and only then is the corresponding ellipsoid in SELL_B tested.
!  Perform no more than max_a_tests tests in SELL_A, and no more than max_b_tests in SELL_B.
!  There may be 0, 1 or more ELLs satisfying the query.  
!  On entry, if id0 is less than 1, a fresh search is started, and the
!  ID of the first ELL (if any) found satisfying the query is returned in id.  If no such
!  ELL is found, then on return id is negative.  After a successful call (id>0), a
!  subsequent call with id0 set to the (positive) value of id from the previous call 
!  resumes the search.  Only those ELLs with ell%tested==.false. are tested, and then
!  ell%tested is set to .true..  The values of ell%tested should be set by an previous 
!  call to ll_set_tested.

   type (ll_type), pointer   :: ll              ! LL
   type (sell_type), pointer :: sell_a, sell_b  ! SELL_A, SELL_B
   real(k_d), intent(in)     :: xa(sell_a%nx)   ! location of query point in SELL_A-space
   real(k_d), intent(in)     :: xb(sell_b%nx)   ! location of query point in SELL_B-space
   integer,   intent(in)     :: max_a_tests     ! maximum number of ELL tests in SELL_A
   integer,   intent(in)     :: max_b_tests     ! maximum number of ELL tests in SELL_B
   integer,   intent(in)     :: id0             ! ID of starting ELL (or <1 for fresh search)
   integer,   intent(out)    :: id              ! ID of ELL satisfying query (or <0 for none)
   integer,   intent(out)    :: a_tests         ! number of ELL tests performed in SELL_A
   integer,   intent(out)    :: b_tests         ! number of ELL tests performed in SELL_B

   type (ell_type), pointer  :: ell_a, ell_b
   integer                   :: id_status, ids_assigned, id_max
   logical                   :: in
 
   id      = -1  !  set null values
   a_tests = 0
   b_tests = 0
   if( max_a_tests <= 0  .or.  max_b_tests <= 0 ) return  !  quick return if no tests allowed

   if( .not.associated(ll) ) then
      write(0,*) 'll_query_pair: ll not associated; stopping'
	  stop
   endif

   if( ll%check > 0 ) then
      if( .not.associated(sell_a) ) then
         write(0,*) 'll_query_pair: sell_a not associated; stopping'
	     stop
      endif

      if( .not.associated(sell_b) ) then
         write(0,*) 'll_query_pair: sell_b not associated; stopping'
	     stop
      endif

	  if( ll%check > 1  .and.  associated( ll%idlist )  .and.  id0 > 0 ) then
	     call id_query( ll%idlist, id0, id_status, ids_assigned, id_max )
		 if( id_status /= 1 ) then
		    write(0,*) 'll_query_pair: id0 is not an assigned ID, id0 = ', id0
	        stop
	     endif
	  endif 
   endif

   if( ll%first  < 1 ) return  !  empty linked list

!  determine id at which to start
   if( id0 < 1 ) then  !  fresh start 
      id = ll%first

   elseif( id0 <= ll%n_pt ) then     ! potentially valid id0

      if( ll%after(id0) > 0 ) then
         id = ll%after(id0)          !  start from ll%after(id0)

	  elseif( id0 == ll%last ) then  ! already at end of list
		 return

      else  !  inconsistency in linked list
	    write(0,*) 'll_query_pair: inconsistency, id0, ll%last = ', id0, ll%last
		stop
      endif

   else  !  invalid id0
      write(0,*) 'll_query_pair: invalid id0 > ll%n_pt ', id0, ll%n_pt
	  stop
   endif

!  loop until covering ellipsoid found or until end of list
   do
      ell_a => sell_a%ell_pt(id)%ell
      ell_b => sell_b%ell_pt(id)%ell

	  if( ell_a%ell_id /= id ) then  !  check consistency  ***
	     write(0,*) 'll_query_pair: inconsistency, ell_a%id /= id ', ell_a%ell_id, id
		 stop
      endif

! test for xa covered by ELL_A

      if( .not.ell_b%tested ) then

	     call sell_ell_pt_in( sell_a%nx, sell_a%ngeom, ell_a%geom(1:sell_a%ngeom), xa, in )
	     a_tests      = a_tests + 1
		 ell_b%tested = .true.

	     if( in ) then  ! in ELL_A: test for xb covered by ELL_B
		    call sell_ell_pt_in( sell_b%nx, sell_b%ngeom, ell_b%geom(1:sell_b%ngeom), xb, in )
	        b_tests = b_tests + 1
	     endif
	  else
         in = .false.
      endif

	  if( in ) exit  !  covering ellipsoid found

	  if( a_tests >= max_a_tests  .or.  b_tests >= max_b_tests ) then
	     id = -1     !   max_tests exceeded
		 exit
	  endif

	  if( ll%after(id) < 1 ) then  !  at end of list
	  
	     if( id == ll%last ) then  ! check for consistency  ***
	        id = -1
		    exit
         else  !  inconsistency in linked list
	        write(0,*) 'll_query_pair: inconsistency, id, ll%last = ', id, ll%last
		    stop
         endif
      endif

	  id = ll%after(id)
   end do

   return
end subroutine ll_query_pair  !--------------------------------------------------------------

subroutine ll_query_pair_list( ll_a, sell_a, sell_b, xa, xb, max_a_tests, max_b_tests, &
                               id0, n_id_max, n_id, id, a_tests, b_tests  )

!  xb is a query point in the same space as the set of ellipsoids SELL_B.
!  xa is the projection of xb onto the (usually lower-dimensional) space containing the 
!  set of ellipsoids SELL_A, which are the projections of those in SELL_B.  
!  Thus, if xb is covered by an ellipsoid in SELL_B, then xa is covered by the corresponding
!  ellipsoid in SELL_A.

!  Given the query point xb, use the linked list LL, to attempt to find up to n_id_max
!  ELLs in SELL_B which cover xb.  Tests are performed first to see if xa is covered by the
!  ellipsoid in SELL_A, and only then is the corresponding ellipsoid in SELL_B tested.
!  Perform no more than max_a_tests tests in SELL_A, and no more than max_b_tests in SELL_B.
!  There may be 0, 1 or more ELLs found to satisfy the query.  
!  The number of covering ELLs found is returned in n_id and their IDs in id(:).
!  On entry, if id0 is less than 1, a fresh search is started.  
!  After a successful call (n_id>0), a subsequent call with id0 set to id(n_id) from the 
!  previous call resumes the search.  Only those ELLs with ell%tested==.false. are tested, and then
!  ell%tested is set to .true..  The values of ell%tested should be set by an previous 
!  call to ll_set_tested.

   type (ll_type), pointer   :: ll_a            ! LL for SELL_A
   type (sell_type), pointer :: sell_a, sell_b  ! SELL_A, SELL_B
   real(k_d), intent(in)     :: xa(sell_a%nx)   ! location of query point in SELL_A-space
   real(k_d), intent(in)     :: xb(sell_b%nx)   ! location of query point in SELL_B-space
   integer,   intent(in)     :: max_a_tests     ! maximum number of ELL tests in SELL_A
   integer,   intent(in)     :: max_b_tests     ! maximum number of ELL tests in SELL_B
   integer,   intent(in)     :: id0             ! ID of starting ELL (or <1 for fresh search)
   integer,   intent(in)     :: n_id_max        ! maximum number of IDs to be returned
   integer,   intent(out)    :: n_id            ! number of covering ELLs found
   integer,   intent(out)    :: id(n_id_max)    ! array of IDs of covering Ells
   integer,   intent(out)    :: a_tests         ! number of ELL tests performed in SELL_A
   integer,   intent(out)    :: b_tests         ! number of ELL tests performed in SELL_B


   type (ell_type), pointer  :: ell_a, ell_b
   integer                   :: idc, id_status, ids_assigned, id_max
   logical                   :: in
 
   n_id    = 0  !  set null values
   a_tests = 0
   b_tests = 0
   id      = 0

   if( max_a_tests <= 0  .or.  max_b_tests <= 0 ) return  !  quick return if no tests allowed

   if( .not.associated(ll_a) ) then
      write(0,*) 'll_query_pair_list: ll_a not associated; stopping'
	  stop
   endif

   if( ll_a%check > 0 ) then
      if( .not.associated(sell_a) ) then
         write(0,*) 'll_query_pair_list: sell_a not associated; stopping'
	     stop
      endif

      if( .not.associated(sell_b) ) then
         write(0,*) 'll_query_pair_list: sell_b not associated; stopping'
	     stop
      endif

	  if( ll_a%check > 1  .and.  associated( ll_a%idlist )  .and.  id0 > 0 ) then
	     call id_query( ll_a%idlist, id0, id_status, ids_assigned, id_max )
		 if( id_status /= 1 ) then
		    write(0,*) 'll_query_pair_list: id0 is not an assigned ID, id0 = ', id0
	        stop
	     endif
	  endif

	  if( n_id_max < 1 ) then
	     write(0,*) 'll_query_pair_list: invalid n_id_max = ', n_id_max
	     stop
	  endif
   endif

   if( ll_a%first  < 1 ) return  !  empty linked list

!  determine ID (idc) at which to start
   if( id0 < 1 ) then  !  fresh start 
      idc = ll_a%first

   elseif( id0 <= ll_a%n_pt ) then     ! potentially valid id0

      if( ll_a%after(id0) > 0 ) then
         idc = ll_a%after(id0)          !  start from ll_a%after(id0)

	  elseif( id0 == ll_a%last ) then  ! already at end of list
		 return

      else  !  inconsistency in linked list
	    write(0,*) 'll_query_pair_list: inconsistency, id0, ll_a%last = ', id0, ll_a%last
		stop
      endif

   else  !  invalid id0
      write(0,*) 'll_query_pair_list: invalid id0 > ll_a%n_pt ', id0, ll_a%n_pt
	  stop
   endif

!  loop until covering ellipsoid found or until end of list
   do
      ell_a => sell_a%ell_pt(idc)%ell
      ell_b => sell_b%ell_pt(idc)%ell

	  if( ell_a%ell_id /= idc ) then  !  check consistency ***
	     write(0,*) 'll_query: inconsistency, ell_a%id /= idc ', ell_a%ell_id, idc
		 stop
      endif

      if( .not.ell_b%tested ) then
        ! test ELL_A
	     call sell_ell_pt_in( sell_a%nx, sell_a%ngeom, ell_a%geom(1:sell_a%ngeom), xa, in )
	     a_tests = a_tests + 1
		 ell_b%tested = .true.

	     if( in ) then  !  in ELL_A, test ELL_B
	        call sell_ell_pt_in( sell_b%nx, sell_b%ngeom, ell_b%geom(1:sell_b%ngeom), xb, in )
	        b_tests = b_tests + 1
	     endif
	  else
         in = .false.
      endif

	  if( in ) then  !  add ID to list
	     n_id     = n_id + 1
		 id(n_id) = idc
		 if( n_id == n_id_max ) exit !  n_id_max IDs found
	  endif

	  if( a_tests >= max_a_tests  .or.  b_tests >= max_b_tests ) exit  !  max. tests reached
	  if( ll_a%after(idc) < 1 ) exit    !  at end of list
	  idc = ll_a%after(idc)
   end do

   return
end subroutine ll_query_pair_list  !--------------------------------------------------------------

subroutine ll_write( ll, lu )

!  Write LL to logical unit LU for checkpointing

   type (ll_type), pointer  :: ll      ! (out) data structure for linked list.
   integer, intent(in)      :: lu

   logical :: opened

   if( .not.associated(ll) ) then
      write(0,*)'ll_write: ll not associated'
	  stop
   endif

   if( lu < 1 ) then
      write(0,*)'ll_write: invalid lu = ', lu
	  stop
   endif

   inquire( lu, opened=opened )
   if( .not.opened ) then
      write(0,*)'ll_write: LU not opened '
	  stop
   endif

   write(lu, err=100) ll%n_pt, ll%check
   write(lu, err=105) ll%first, ll%last
   write(lu, err=110) ll%before(1:ll%n_pt), ll%after(1:ll%n_pt)

   return

100   write(0,*)'ll_write: error writing n_pt'
      stop

105   write(0,*)'ll_write: error writing first'
      stop

110   write(0,*)'ll_write: error writing before'
      stop

end subroutine ll_write  !------------------------------------------------------------

subroutine ll_read( ll, lu, idlist )

!  Initialize and read LL from logical unit LU.

   type (ll_type), pointer                :: ll  
   integer, intent(in)                    :: lu
   type (id_list_type), pointer, optional :: idlist

   integer :: n_pt, check
   logical :: opened

   if( lu < 1 ) then
      write(0,*)'ll_read: invalid lu = ', lu
	  stop
   endif

   inquire( lu, opened=opened )
   if( .not.opened ) then
      write(0,*)'ll_read: LU not opened '
	  stop
   endif

   read(lu, err=100) n_pt, check

   if( present(idlist) ) then
      call ll_initialize( ll, check=check, n_pt=n_pt, idlist=idlist )
   else
      call ll_initialize( ll, check=check, n_pt=n_pt )
   endif

   read(lu, err=105) ll%first, ll%last
   read(lu, err=110) ll%before(1:ll%n_pt), ll%after(1:ll%n_pt)

   return

100   write(0,*)'ll_read: error reading n_pt'
      stop

105   write(0,*)'ll_read: error reading first'
      stop

110   write(0,*)'ll_read: error reading before'
      stop

end subroutine ll_read  !------------------------------------------------------------

subroutine ll_param_set( ll, check ) 

! set parameters

   type (ll_type), pointer       :: ll     ! data structure for linked list. 
   integer, intent(in), optional :: check  ! level of checking [0,1,2]

   if( .not.associated(ll) ) then
      write(0,*) 'll_param_set: ll not associated; stopping'
	  stop
   endif

   if( present(check) ) then
      if( check < 0  .or.  check >2 ) then
	     write(0,*) 'll_param_set: invalid  check = ', check
	     stop
	  endif

      ll%check = check
   endif

   return
end subroutine ll_param_set  !--------------------------------------------------------------

subroutine ll_param_get( ll, check ) 

! get parameters

   type (ll_type), pointer        :: ll     ! data structure for linked list. 
   integer, intent(out), optional :: check  ! level of checking 

   if( .not.associated(ll) ) then
      write(0,*) 'll_param_get: ll not associated; stopping'
	  stop
   endif

   if( present(check) ) check = ll%check

   return
end subroutine ll_param_get  !--------------------------------------------------------------


end module sell_ll
 
