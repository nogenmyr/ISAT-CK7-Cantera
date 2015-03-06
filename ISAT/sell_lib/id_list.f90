!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module id_list

!  Data structure and procedures for generating unique ID numbers.

!  S.B. Pope 1/12/06, 5/14/06

!  ID numbers are the set of positive intergers.  
!  Each ID number is either "assigned" or "available."

!  Usage:
!    In the calling routine, declare  idlist  of type (id_list_type), pointer.

! call id_init(idlist,check,kseed)  to initialize -- initially all strictly positive
!                                integers are available as ID numbers.
! call id_get(idlist,id)      to assign a new ID number, id.
! call id_return(idlist,id)   to return an ID number, no longer assigned.
! call id_destroy( idlist )   to destroy the data structure
! call id_stats( idlist, ids_assigned, id_max )  to return ids_assigned and id_max
! call id_query( idlist, id, id_status, ids_assigned, id_max ) to query id and
!                             check consistency of idlist
! call id_rand( idlist, id, ids_assigned ) to return an assigned ID number 
!                             selected at random
! call id_rand_array( idlist, id_array ) to return an array of assigned
!                             ID numders selected at random
! call id_set_check( idlist, check ) to change checking level
! call id_print( idlist, lu ) to write the data structure (for diagnostic purposes)
! call id_write( idlist, lu ) to write idlist for checkpointing
! call id_read( idlist, lu )  to read idlist from checkpoint file

!  id_max is the largest ID ever assigned.
!  The value of id returned by id_get is between 1 and id_max+1 (inclusive).

! Implementation:
!  idlist%first points to the first item in a linked list of type id_item_type.
!  The last item in the linked list has  item%id = id_max+1  and  item%next nullified.
!  For all other items in the linked list (if any),  item%id  is an ID once assigned,
!  but now available (through a call to id_return).  
!  Thus every value of item%id is an available ID.
!  Note that the last ID in the list is the largest, but the others are not ordered.

! ISAT random number generator:
!  In order to assure repeatability and machine independence when used by ISATAB,
!  the random numbers used are from the ISAT routine rnu contained in the module isat_rnu.
!  The random number sequence used by rnu by setting kseed>0 in the call to id_init;
!  otherwise, the current seed values at the time of calling id_init are used.
!  For use without ISAT, set isat_random=.false. immediately below, and comment out
!  isat-related lines, which are indicated by ! ISAT.

   implicit none
   logical, parameter, private :: isat_random = .true.  

   type :: id_list_type
      type(id_item_type), pointer  :: first  !  first item in linked list
      integer :: check        !  level of checking: =0, minimal; =1, maximal
	  integer :: id_max       !  largest ID number ever assigned
      integer :: ids_assigned !  number of ID numbers currently assigned
      integer :: is1, is2     !  random number seeds
   end type id_list_type

   type :: id_item_type
      integer :: id                        ! ID number which could be assigned
	  type(id_item_type), pointer :: next  ! next item in linked list
   end type id_item_type

contains  !====================================================================

subroutine id_init( idlist, check, kseed )

!  allocate and initialize idlist

   use isat_rnu  ! ISAT

   type (id_list_type), pointer  :: idlist  !  sidlist to be initialized
   integer, intent(in), optional :: check   !  level of checking [0,1]
   integer, intent(in), optional :: kseed   !  for kseed > 0, set rnu seeds to 
                                            !  the kseed-th seed

   type (id_item_type),pointer   :: first   !  first entry in linked list
   integer :: ist1, ist2
   
   allocate(idlist)
   idlist%id_max       = 0
   idlist%ids_assigned = 0

   if( present(check) ) then
      call id_set_check( idlist, check )
   else
	  idlist%check = 1  ! by default, checking of input is performed
   endif
   
   if( isat_random ) then  !  set random number seeds
      call rnuget( idlist%is1, idlist%is2 )         !  ISAT
      if( present(kseed) ) then
         if( kseed > 0 ) then
            call rnuget( ist1, ist2 )  ! store current seeds               ! ISAT
            call rnused( kseed )       ! re-set seed                       ! ISAT
            call rnuget( idlist%is1, idlist%is2 )  ! get seeds for ID    ! ISAT
            call rnuput( ist1, ist2 )  ! restore original seeds            ! ISAT
         endif
      endif   
   endif

   allocate(first)
   first%id = 1
   nullify(first%next)

   idlist%first => first

   return

end subroutine id_init  !------------------------------------------------------

subroutine id_get( idlist, id )

!  return first available id from idlist

   type (id_list_type), pointer  :: idlist 
   integer, intent(out)          :: id

   type (id_item_type), pointer  :: id_temp

   if( .not.associated(idlist) ) then
      write(0,*)'id_get: error, idlist not associated'
	  stop
   endif

   id = idlist%first%id  !  return first ID in list

   if( .not.associated(idlist%first%next) ) then
!    if  idlist%first  is sole entry in list, increment id
      idlist%id_max   = idlist%first%id
      idlist%first%id = idlist%first%id + 1  
   else
!    otherwise make idlist%first point to next entry, then deallocate original idlist%first

	    id_temp       => idlist%first
		idlist%first  => idlist%first%next
		deallocate(id_temp)
   endif 
   
   idlist%ids_assigned = idlist%ids_assigned + 1
   return  

end subroutine id_get  !------------------------------------------------------

subroutine id_return( idlist, id )

!  return id to idlist (i.e., make id available)

   type (id_list_type), pointer :: idlist   
   integer, intent(in)          :: id      
   
   type (id_item_type), pointer, save :: id_new 

   if( .not.associated(idlist) ) then
      write(0,*)'id_return: error, idlist not associated'
	  stop
   endif

   if( idlist%check > 0 ) then  !  check value of id

      id_new => idlist%first  !  check that id is not already in list
      do
         if( id_new%id == id ) then
            write(0,*)'id_return: error ID already in list, stopping; id = ', id
		    stop
         endif
	     if( .not.associated( id_new%next ) ) exit
	     id_new => id_new%next
      end do

	  if( id < 1  .or.  id >= id_new%id ) then
            write(0,*)'id_return: error ID out of range, stopping; id, id_max = ', id, id_new%id-1
		    stop
	  endif
      nullify( id_new )
   endif

! add new item to the beginning of the list
   allocate( id_new )
   id_new%id    = id
   id_new%next  => idlist%first
   idlist%first => id_new

   idlist%ids_assigned = idlist%ids_assigned - 1

   return

end subroutine id_return  !------------------------------------------------------

subroutine id_destroy( idlist )

!  destroy and deallocate data associated with idlist

   type (id_list_type), pointer  :: idlist 

   type (id_item_type), pointer  :: id_this, id_next

   if( .not.associated(idlist) ) then
      write(0,*)'id_destroy: error, idlist not associated'
	  stop
   endif

   id_this => idlist%first
   do
      if( .not.associated(id_this%next) ) then
	     deallocate(id_this)
		 exit
	  endif

	  id_next => id_this%next
	  deallocate( id_this )
	  id_this => id_next
   end do

   deallocate( idlist )

   return  

end subroutine id_destroy  !------------------------------------------------------

subroutine id_stats( idlist, ids_assigned, id_max )

! Return:
!     ids_assigned - number if IDs currently assigned (or -1 if idlist not associated)
!     id_max       - maximum id ever assigned

   type (id_list_type), pointer  :: idlist 
   integer, intent(out)          :: ids_assigned, id_max

   if( associated(idlist) ) then  
      ids_assigned = idlist%ids_assigned
      id_max       = idlist%id_max
   else
	  ids_assigned = -1
	  id_max       = -1
   endif

   return  

end subroutine id_stats  !------------------------------------------------------

subroutine id_query( idlist, id, id_status, ids_assigned, id_max )

!  Query and check ID list.
!  There are two modes of operation:

!  Mode 1, for id >= 0
!    id = value of ID to be queried (input)

!    id_status - status of id (output)
!              =  0 - id=0 (i.e., invalid)
!              =  1 - currently assigned
!              =  2 - once assigned, but returned
!              =  3 - never assigned
!              = -2 - idlist not associated

!  Mode 2, for id < 0
!    id = -1 - check consistency of list
!    id = -2 - check consistency of list and write message if inconsistent
!    id = -3 - check consistency of list and write message

!    id_status - status of list (output)
!              =  0 - consistent
!              = -1 - inconsistent
!              = -2 - idlist not associated

!  Both modes
!     ids_assigned - number if IDs currently assigned
!     id_max       - maximum id ever assigned

   type (id_list_type), pointer  :: idlist  !  first entry in linked list
   integer, intent(in)           :: id
   integer, intent(out)          :: id_status, ids_assigned, id_max

   type (id_item_type), pointer  :: id_temp
   integer                       :: n_items, j, k
   integer, allocatable          :: i(:)

   if( .not.associated(idlist) ) then  !  check for idlist not being associated
      id_status    = -2
      ids_assigned = 0
      id_max       = 0
	  if( id < -1 ) write(0,*)'id_query: idlist not associated'
	  return
   endif

   ids_assigned = idlist%ids_assigned
   id_max       = idlist%id_max

   if( id == 0 ) then
      id_status = 0
	  return
   elseif( id > id_max ) then
      id_status = 3
	  return
   endif

!  traverse list, looking for id
   id_status = 1   !  assume assigned
   id_temp   => idlist%first
   n_items   = 0

   do
      n_items = n_items + 1
	  if( id_temp%id == id ) then
	     id_status = 2
		 exit
	  endif

	  if( .not.associated( id_temp%next ) ) exit
	  id_temp => id_temp%next
   end do

   if( id >= 0 ) return !  all done for mode 1

! mode 2: check consistency: 
   if( n_items /= 1+id_max-ids_assigned ) then
      write(0,*)'id_query: inconsistent n_items, id_max, ids_assigned = ', &
                                        n_items, id_max, ids_assigned
      stop
   endif

!put list into array i
   id_status = 0
   allocate( i(n_items) )

   j = 0
   id_temp=> idlist%first
   do 
      j = j + 1
      i(j) = id_temp%id
	  if( .not.associated( id_temp%next ) ) exit
	  id_temp => id_temp%next
   end do

   do j = 1, n_items  !  check that IDs are positive
      if( i(j) < 1 ) then
	     id_status = -1
		 if( id < -1 ) write(0,*)'id_query: inconsistent, i(j) < 1; j, i(j) = ', j, i(j)
	  endif

	  do k = 1, j-1 !  check that IDs are distinct
	     if( i(j) == i(k) ) then
		    id_status = -1
		    if( id < -1 ) write(0,*)'id_query: inconsistent, i(j) = i(k); ', &
			                        'j, k, i(j), i(k) = ', j, k, i(j), i(k) 
		 endif
	  end do
   end do


   do j = 1, n_items - 1  !  check that IDs are less than i(n_items)
      if( i(j) >= i(n_items) ) then
	     id_status = -1
		 if( id < -1 ) write(0,*)'id_query: inconsistent, i(j) >= i(n_items); ', &
		                         'j, i(j), i(n_items) = ', j, i(j), i(n_items)
	  endif
   end do

   if( id == -3  .and.  id_status == 0 ) write(0,*) 'id_query: consistent', &
              ' ids_assigned, id_max, n_items = ', ids_assigned, id_max, n_items

   deallocate(i)

   return  

end subroutine id_query  !------------------------------------------------------

subroutine id_rand( idlist, id, ids_assigned )

!  Return an assigned ID selected at random (uniformly) from the list idlist.

   use isat_rnu  ! ISAT
   type (id_list_type), pointer  :: idlist  !  first entry in linked list
   integer, intent(out)          :: id      !  assigned ID (selected at random)
   integer, intent(out)          :: ids_assigned ! total number of IDs assigned

   type (id_item_type), pointer  :: id_temp
   integer                       :: id_max, n_less, n_less_last, ist1, ist2
   real                          :: rand

   if( .not.associated(idlist) ) then  !  empty list
	  id           = 0
      ids_assigned = 0
	  return
   endif

   id_max       = idlist%id_max
   ids_assigned = idlist%ids_assigned

!  determine index of ID selected at random
   if( isat_random ) then
      call rnuget( ist1, ist2 )              ! ISAT
      call rnuput( idlist%is1, idlist%is2 )  ! ISAT
      
      call rnu( rand )                       ! ISAT
      
      call rnuget( idlist%is1, idlist%is2 )  ! ISAT
      call rnuput( ist1, ist2 )              ! ISAT
   else
      call random_number( rand )
   endif

   id = int( rand * ids_assigned ) + 1
   id = min( id, ids_assigned )

!  traverse list to find id-th assigned ID
   n_less_last = 0
   do  ! loop until n_less does not change
      
      id_temp => idlist%first
	  n_less      = 0
      do  !  cont n_less
         if( id_temp%id <= id ) then
		    n_less = n_less + 1
		 endif

	     if( .not.associated( id_temp%next ) ) exit
	     id_temp => id_temp%next
      end do

	  if( n_less == n_less_last ) exit

	  id           = id + n_less - n_less_last
      n_less_last = n_less
   end do

   return  

end subroutine id_rand  !------------------------------------------------------

subroutine id_rand_array( idlist, id_array )

!  Return in id_array, all assigned IDs from the list idlist in a random order.

   use isat_rnu  ! ISAT
   type (id_list_type), pointer  :: idlist  
   integer, intent(out)          :: id_array(idlist%ids_assigned) 

   type (id_item_type), pointer  :: id_temp
   integer                       :: idt(idlist%id_max), i, j, id_count, idj, ist1, ist2
   real                          :: rand(idlist%ids_assigned)
 
   if( .not.associated(idlist) ) then  !  empty list
      write(0,*)'id_rand_array: idlist not associated'
	  stop
   endif

   if( idlist%ids_assigned < 1 ) then  !  empty list
      write(0,*)'id_rand_array: ids_assigned < 1 ', idlist%ids_assigned
	  stop
   endif

!  For assigned IDs,   set idt(id) = id
!  For unassigned IDs, set idt(id) = -1

   do i = 1, idlist%id_max
      idt(i) = i  !  assume ID assigned
   end do

   id_count =  idlist%id_max
   id_temp  => idlist%first  
   do
      if( .not.associated(id_temp%next)  ) exit  !  last ID never assigned
	  idt( id_temp%id ) = -1  ! for unassigned IDs, set id_array to -1
	  id_count = id_count - 1
      id_temp  => id_temp%next
   end do

   if( id_count /= idlist%ids_assigned ) then
      write(0,*)'id_rand_array: count error ', id_count, idlist%ids_assigned
	  stop
   endif

!  copy assigned IDs into id_array

   j = 0
   do i = 1, idlist%id_max
      if( idt(i) > 0 ) then
	     j = j + 1
		 id_array(j) = idt(i)
       endif
   end do

!  randomly re-order

   if( isat_random ) then
      call rnuget( ist1, ist2 )              ! ISAT
      call rnuput( idlist%is1, idlist%is2 )  ! ISAT
      
      call rnu( rand )                       ! ISAT
      
      call rnuget( idlist%is1, idlist%is2 )  ! ISAT
      call rnuput( ist1, ist2 )              ! ISAT   
   else
      call random_number( rand )
   endif

   do i = 1, id_count - 1
      j = i+ (id_count+1-i) * rand(i)
	  j = min( j, id_count )  !  select j between i and id_count
	  idj = id_array(j)       !  commute id_array(i) and id_array(j)
	  id_array(j) = id_array(i)
	  id_array(i) = idj
   end do

   return  

end subroutine id_rand_array  !------------------------------------------------------

subroutine id_set_check( idlist, check )

!  set checking check: check=0 - minimal checking, check=1 - maximal checking
!  (check=2 is changed to idlist%check=1)

   type (id_list_type), pointer :: idlist  !  first entry in linked list
   integer, intent(in)          :: check

   if( .not.associated(idlist) ) then
      write(0,*)'id_set_check: error, idlist not associated'
	  stop
   endif

   if( check < 0  .or.  check  > 2 ) then
      write(0,*) 'id_set_check: check outside valid range [0,1], check= ', check
	  stop
   endif

   idlist%check = min( check, 1 ) 

   return  

end subroutine id_set_check  !------------------------------------------------------

subroutine id_print( idlist, lu )

!  write list of IDs on unit lu

   type (id_list_type), pointer  :: idlist  !  pointer to start of list
   integer, intent(in)           :: lu

   type (id_item_type), pointer  :: id
   integer :: n
   
   if( .not.associated(idlist) ) then
      write(lu,*)'id_print:  idlist not associated'
	  return
   endif

   id => idlist%first
   n   = 0

   do
      n = n + 1
      write(lu,'(2i10)') n, id%id

     if( .not.associated(id%next) ) exit

     id => id%next
   end do

   return

end subroutine id_print  !------------------------------------------------

subroutine id_write( idlist, lu )  !  write idlist to file for checkpointing

   type (id_list_type), pointer  :: idlist  !  pointer to start of list
   integer, intent(in)           :: lu      !  logical unit for output (must be open)

   type (id_item_type), pointer  :: id_temp
   integer                       :: n_items
   logical                       :: opened
   
   if( .not.associated(idlist) ) then  !  check input
      write(0,*)'id_write:  idlist not associated'
	  stop
   endif

   if( lu < 1 ) then
      write(0,*)'id_write:  lu < 1, ', lu
	  stop
   endif

!  check file status
   inquire( lu, err=100, opened=opened )

   if( .not.opened ) then
      write(0,*)'id_write:  file not opened'
	  stop
   endif

   write( lu, err=105 ) idlist%check

!  traverse list to count entries
   n_items = 0
   id_temp => idlist%first

   do  
      n_items = n_items + 1
	  if( .not.associated( id_temp%next ) ) exit
	  id_temp => id_temp%next
   end do

   write( lu, err=110 ) n_items, idlist%id_max, idlist%ids_assigned, idlist%is1, idlist%is2

   id_temp => idlist%first  !  traverse list and write id
   do  
      write( lu, err=120 ) id_temp%id
	  if( .not.associated( id_temp%next ) ) exit
	  id_temp => id_temp%next
   end do

   return

100   continue
      write(0,*)'id_write:  error making file inquiry'
	  stop

105   continue
      write(lu,*)'id_write:  error writing check'
	  stop

110   continue
      write(lu,*)'id_write:  error writing n_items'
	  stop

120   continue
      write(0,*)'id_write:  error writing id, check'
	  stop

end subroutine id_write  !--------------------------------------------------------------

subroutine id_read( idlist, lu )  !  read idlist from checkpoint file

!  idlist should be allocated prior to call

   type (id_list_type), pointer  :: idlist  !  pointer to start of list
   integer, intent(in)           :: lu      !  logical unit for reading (must be open)

   type (id_item_type), pointer  :: id_temp, id_last
   integer :: n_items, k, check
   logical :: opened
   
   if( lu < 1 ) then
      write(0,*)'id_read:  lu < 1, ', lu
	  stop
   endif

!  check file status
   inquire( lu, err=100, opened=opened )

   if( .not.opened ) then
      write(0,*)'id_read:  file not opened'
	  stop
   endif

   read( lu, err=105 ) check
   call id_init( idlist, check )

!  traverse list to count entries
   
   read( lu, err=110 ) n_items, idlist%id_max, idlist%ids_assigned, idlist%is1, idlist%is2
   if( n_items < 1 ) then
      write(0,*)'id_read:  n_items < 1 ', n_items
	  stop
   endif

   do k = 1, n_items  !  loop over list
      allocate( id_temp )
      read( lu, err=120 ) id_temp%id
	  nullify( id_temp%next )

	  if( k == 1 ) then
	     idlist%first => id_temp
	  else
	     id_last%next => id_temp
	  endif

	  id_last => id_temp
	  nullify( id_temp )
   end do

   return

100   continue
      write(0,*)'id_read:  error making file inquiry'
	  stop

105   continue
      write(0,*)'id_read:  error reading check'
	  stop

110   continue
      write(0,*)'id_read:  error reading n_items'
	  stop

120   continue
      write(0,*)'id_read:  error writing id, check'
	  stop

end subroutine id_read  !---------------------------------------------------------------

end module id_list 
