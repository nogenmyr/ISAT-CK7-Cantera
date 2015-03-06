!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the x2f_mpi      !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module sort_tools

contains ! ----------------------------------------------------------------

subroutine bucket_sort( nlist, ld, nelem, list, bucket_assigns, &
      nbucket, bucket_totals, list_origplaces, orig_slots )

! Given a list of vectors in the form of a rank-2 array (1:nelem,1:nlist),
! together with an array of bucket assignments (1:nlist) associating each
! vector with a "bucket number" (integer value between 1 and nbucket):
! sort the list (but not the assignments) by bucket number, then return the
! sorted list along with an array orig_slots giving the original locations.
! Only the first nelem (<= ld) components of each list item are relocated.
! The input list may be long, so bucket_totals is an additional input array
! telling bucket_sort how many vectors to expect in each bucket.

! Note #1: the sort is done in place, in one pass.
! Note #2: the contents of each bucket are not further sorted.

! When a bucket number outside of the range 1:nbucket is encountered, the
! associated vector is placed in a "bad bucket" at the end of the list.
! Therefore, sum( bucket_totals ) may be less than nlist.

implicit none

integer, intent(in)             :: nlist, ld, nelem
real(kind(1.d0)), intent(inout) :: list(ld,nlist)
integer, intent(in)             :: bucket_assigns(nlist)
integer, intent(in)             :: nbucket, bucket_totals(nbucket)
integer, intent(inout)          :: orig_slots(nlist)
integer, intent(inout)          :: list_origplaces(nlist)

! Local variables:
! bucket_end_pts(j) is the ending point of bucket j in final, sorted list
! bucket_cur_pts(j) starts at bucket_end_pts(j-1)+1, moves up during sort
! dimension is nbucket+1 due to extra, "bad" bucket at end of list

integer :: bucket_cur_pts(nbucket+1), bucket_end_pts(nbucket+1)
integer :: total_all_good, marker, j, current_total, i, jdest, jup
integer :: chainstart, listplace, nextmove, backchain
logical :: done_bucket(nbucket+1), found_start, found_end
real(kind(1.d0)) :: temp_vec(nelem)
integer :: temp_place

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The idea of the algorithm is to scan the list for the first item A that's
! in the wrong bucket for its position.  Go to the correct bucket for A and
! find a second item B that's also in the wrong place.  Proceed to plan out
! an entire chain of displacements A -> B -> C -> ... ending at some item Z
! (say) which is the last item needed to complete bucket j.  Thus, if there
! are three items--including A--that are out of place in bucket j, the
! chain must begin at A and return to bucket j three times, ending at A.
! It is natural and efficient to store this plan in the orig_slots array by
! setting orig_slots(B) = A, etc.
!
! To sort the items, first copy Z into temp_vec.  Then execute the chain
! backwards from Z: Y -> Z, X -> Y, ... A -> B.  Finally, put the contents
! of temp_vec into slot A.  Increment j by 1.  If the new bucket j is not
! yet complete, start a new chain, and so on, until all buckets are done.
!
! Clearly this is an order-N algorithm, because (nlist + 2*nbucket) is the
! upper limit of moves.  Generally this upper limit will never be reached,
! as items that are already in the right location don't move.  Moreover,
! all the 2*nbucket copies to/from temp_vec are seldom required, as some
! buckets k > j may be finished off in the process of completing bucket j.
! Obviously the algorithm works well when nbucket << nlist, but it also
! works well for large values of nbucket, since the probability of a given
! chain returning to bucket j declines as nbucket grows.  In fact, it can
! be shown that the length of chains always tends to be of order nlist,
! independent of nbucket, so only a few big chains ever need to be built.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

total_all_good = sum( bucket_totals(1:nbucket) )

if ( total_all_good > nlist ) then
   write(0,*) 'bucket_sort: bad bucket_totals array, sum exceeds nlist'
   return
end if

! set some boundaries for the ultimate, sorted list
! empty buckets will have bucket_cur_pts(j) > bucket_end_pts(j)
! later this will also become true for completed buckets
! for empty and (eventually) completed buckets, done_bucket(j) = .true.

done_bucket = .false.
marker = 1
do j = 1, nbucket + 1
   if ( j < nbucket + 1 ) then
      current_total = bucket_totals(j)
   else
      current_total = nlist - total_all_good
   end if
   if ( current_total == 0 ) done_bucket(j) = .true.
   bucket_cur_pts(j) = marker
   marker = marker + current_total
   bucket_end_pts(j) = marker - 1
end do

orig_slots = (/ ( i, i = 1, nlist ) /)

do j = 1, nbucket + 1

   ! begin to plan the unique swap chain that starts in bucket j:
   ! consider only places in bucket j that aren't yet done;
   ! whole buckets can and will be skipped as the algorithm proceeds

   ! scan for the next item out of place, if any, and start chain there
   found_start = .false.
   do while ( .not. ( found_start .or. done_bucket(j) ) )

      chainstart = bucket_cur_pts(j)
      jdest = bucket_assigns(chainstart)
      if ( ( jdest < 1 ) .or. ( jdest > nbucket ) ) jdest = nbucket + 1

      if ( jdest /= j ) found_start = .true.

      bucket_cur_pts(j) = bucket_cur_pts(j) + 1
      if ( bucket_cur_pts(j) > bucket_end_pts(j) ) done_bucket(j) = .true.

   end do

   if ( found_start ) then

      ! continue to map out the swap chain until bucket j is exhausted
      listplace = chainstart
      found_end = .false.
      do while ( .not. found_end )

         nextmove = bucket_cur_pts(jdest)

         if ( nextmove > bucket_end_pts(jdest) ) then
            if ( jdest == j ) then
               found_end = .true.
            else
               write(0,*) 'bucket_sort: can''t move into bucket', jdest
               return
            end if

         else
            jup = bucket_assigns(nextmove)
            if ( ( jup < 1 ) .or. ( jup > nbucket ) ) jup = nbucket + 1

            if ( jup == jdest ) then  ! pass it by; it's in the right spot
               bucket_cur_pts(jdest) = bucket_cur_pts(jdest) + 1

            else  ! add a link to the chain
               orig_slots(nextmove) = listplace
               listplace = nextmove
               bucket_cur_pts(jdest) = bucket_cur_pts(jdest) + 1
               if ( bucket_cur_pts(jdest) > bucket_end_pts(jdest) ) then
                  done_bucket(jdest) = .true.
               end if
               jdest = jup

            end if

         end if

      end do

      ! execute the swap chain
      temp_vec(1:nelem) = list(1:nelem,chainstart)
      temp_place        = list_origplaces(chainstart)
      list_origplaces(chainstart) = list_origplaces(listplace)
      list(1:nelem,chainstart) = list(1:nelem,listplace)
      orig_slots(chainstart) = listplace
      backchain = orig_slots(listplace)
      do while ( backchain /= chainstart )
         list(1:nelem,listplace) = list(1:nelem,backchain)
         list_origplaces(listplace) = list_origplaces(backchain)
         listplace = backchain
         backchain = orig_slots(listplace)
      end do
      list(1:nelem,listplace) = temp_vec(1:nelem)
      list_origplaces(listplace) = temp_place

   end if

end do

end subroutine bucket_sort


! --------------------------------------------------------------------

subroutine pigeonhole_sort( nlist, ld, nelem, list, assigns, orig_slots )

! Given a list of vectors in the form of a rank-2 array (1:nelem,1:nlist),
! together with an array of pigeonhole assignments (1:nlist) associating
! each vector with a UNIQUE pigeonhole number (integer value between 1 and
! nlist): sort the list (but not the assignments) by bucket number and
! return the sorted list.  Only the first nelem (<= ld) components of each
! list item are relocated.  The original locations, if desired, can be
! deduced from the assigns array.  Each pigeonhole should get one and only
! one item in the list assigned to it.  This sort is useful for restoring
! the original order of a previously sorted list.

! Note: the sort is done in place, in one pass.

implicit none

integer, intent(in)             :: nlist, ld, nelem
real(kind(1.d0)), intent(inout) :: list(ld,nlist)
integer, intent(in)             :: assigns(nlist)
integer, intent(out)            :: orig_slots(nlist)

integer          :: i, j, thisone, moveto, backchain
real(kind(1.d0)) :: temp_vec(nelem)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The pigeonhole sort is a special case of the bucket sort; therefore, this
! routine uses the same algorithm as bucket_sort above.  The differences
! are that the code is simpler and the argument list is shorter, due to the
! one-to-one correspondence between list items and buckets.  Technically,
! a pigeonhole sort may have more buckets than list items; however, that
! case is not dealt with here.  It could be treated through a separate
! routine mapping the unique pigeonhole assignments onto the range 1:nlist.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

orig_slots = (/ ( i, i = 1, nlist ) /)

do j = 1, nlist

   if ( assigns(orig_slots(j)) /= j ) then

      thisone = j
      moveto = assigns(thisone)
      do while ( moveto /= j )
         orig_slots(moveto) = thisone
         thisone = moveto
         moveto = assigns(thisone)
      end do

      ! execute the swap chain
      temp_vec(1:nelem) = list(1:nelem,j)
      list(1:nelem,j) = list(1:nelem,thisone)
      orig_slots(j) = thisone
      backchain = orig_slots(thisone)
      do while ( backchain /= j )
         list(1:nelem,thisone) = list(1:nelem,backchain)
         thisone = backchain
         backchain = orig_slots(thisone)
      end do
      list(1:nelem,thisone) = temp_vec(1:nelem)

   end if

end do

end subroutine pigeonhole_sort


! --------------------------------------------------------------------

! The following subroutine for performing a bucket sort works correctly,
! but the algorithm in bucket_sort above is more efficient and is thus
! preferred.  The subroutine is mainly included for completeness, in case 
! some future situation suggests it as a superior approach.

! --------------------------------------------------------------------

subroutine marble_sort( nlist, ld, nelem, list, bucket_assigns, &
      nbucket, bucket_totals, orig_slots)

! Given a list of vectors in the form of a rank-2 array (1:nelem,1:nlist),
! together with an array of bucket assignments (1:nlist) associating each
! vector with a unique bucket number (integer value between 1 and nbucket):
! sort the list (but not the assignments) by bucket number, then return the
! sorted list along with an array orig_slots giving the original locations.
! Only the first nelem (<= ld) components of each list item are relocated.
! The input list may be long, so bucket_totals is an additional input array
! telling bucket_sort how many vectors to expect in each bucket.

! Note #1: the sort is done in place, in one pass.
! Note #2: the contents of each bucket are not further sorted.

! When a bucket number outside of the range 1:nbucket is encountered, the
! associated vector is placed in a "bad bucket" at the end of the list.
! Therefore, sum( bucket_totals ) may be less than nlist.

implicit none

integer, intent(in)             :: nlist, ld, nelem
real(kind(1.d0)), intent(inout) :: list(ld,nlist)
integer, intent(in)             :: bucket_assigns(nlist)
integer, intent(in)             :: nbucket, bucket_totals(nbucket)
integer, intent(out)            :: orig_slots(nlist)

! Local variables:
! bucket_end_pts(j) is the ending point of bucket j in final, sorted list
! bucket_cur_pts(j+1) starts at bucket_end_pts(j)+1, moves up during sort
! dimension is nbucket+1 due to extra, "bad" bucket at end of list

integer :: bucket_cur_pts(nbucket+1), bucket_end_pts(nbucket+1)
integer :: total_all_good, marker, j, current_total, i, jdest, jup
integer :: chainstart, listplace, nextmove, backchain
logical :: done_bucket(nbucket+1), found_start, found_end
real(kind(1.d0)) :: temp_vec(nelem)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Picture of the algorithm: imagine a row of marbles that you wish to sort
! by color.  There are M distinct colors in this row, but you only care
! about sorting some number nbucket (= M, maybe) of these colors.  In the
! end, you want your particular colors to appear in a certain order along
! the row.  You've counted how many marbles there are in each color, so
! you know in advance which color is correct for each position in the row.
! The colors you don't care about will be assigned to the end of the row.
!
! What do you do?
!
! Start at one end and pick up the first marble that's the wrong color for
! its position in the row.  Determine the range of positions in which it
! belongs.  With your other hand, pick up a wrong-color marble from that
! range and replace it with the first marble.  Then look at this second
! marble, determine its correct range, etc., alternating hands as you go.
! Every so often the destination turns out to be the empty place where
! you picked up the first marble.  When that occurs, fill it with the one
! in your hand and restart the procedure, scanning down the row to find
! the next marble that's out of place.  Continue until all are sorted.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

total_all_good = sum( bucket_totals(1:nbucket) )

if ( total_all_good > nlist ) then
   write(0,*) 'marble_sort: bad bucket_totals array, sum exceeds nlist'
   return
end if

! set some boundaries for the ultimate, sorted list
! empty buckets will have bucket_cur_pts(j) > bucket_end_pts(j)
! later this will also become true for completed buckets
! for empty and (eventually) completed buckets, done_bucket(j) = .true.

done_bucket = .false.
marker = 1
do j = 1, nbucket + 1
   if ( j < nbucket + 1 ) then
      current_total = bucket_totals(j)
   else
      current_total = nlist - total_all_good
   end if
   if ( current_total == 0 ) done_bucket(j) = .true.
   bucket_cur_pts(j) = marker
   marker = marker + current_total
   bucket_end_pts(j) = marker - 1
end do

orig_slots = (/ ( i, i = 1, nlist ) /)

! Start sort at listplace = 1, in the first nonempty (post-sort) bucket...
! sort procedure is as follows: [A] check bucket assignment of list item;
! [B] if item is already in the right bucket, increment listplace; else,
! (1) "pick up" item from listplace by copying it into the temporary array
! (2) find next out-of-place item in the destination bucket, pick it up too
! (3) "put down" first copied item in correct bucket, get next destination
! (4) repeat (2) & (3) until destination is once again listplace (bucket j)
! (5) put down the last copied item into listplace, completing a swap cycle
! [C] increment listplace and restart procedure at [A] until list is done
! ...progress is guaranteed because on each trip through a "do while" loop,
! one item is either put in, or confirmed to be in, the right location, AND
! a bucket position counter is incremented--proving the sort is order-N

do j = 1, nbucket + 1

   ! begin to plan the unique swap chain that starts in bucket j:
   ! consider only places in bucket j that aren't yet done;
   ! whole buckets can and will be skipped as the algorithm proceeds

   ! scan for the next item out of place, if any, and start chain there
   found_start = .false.
   do while ( .not. ( found_start .or. done_bucket(j) ) )

      chainstart = bucket_cur_pts(j)
      jdest = bucket_assigns(chainstart)
      bucket_cur_pts(j) = bucket_cur_pts(j) + 1
      if ( ( jdest < 1 ) .or. ( jdest > nbucket ) ) jdest = nbucket + 1

      if ( jdest /= j ) found_start = .true.
      if ( bucket_cur_pts(j) > bucket_end_pts(j) ) done_bucket(j) = .true.

   end do

   if ( found_start ) then

      ! continue to map out the swap chain until bucket j is exhausted
      listplace = chainstart
      found_end = .false.
      do while ( .not. found_end )

         nextmove = bucket_cur_pts(jdest)

         if ( nextmove > bucket_end_pts(jdest) ) then
            if ( jdest == j ) then
               found_end = .true.
            else
               write(0,*) 'marble_sort: can''t move into bucket', jdest
               return
            end if

         else
            jup = bucket_assigns(nextmove)
            if ( ( jup < 1 ) .or. ( jup > nbucket ) ) jup = nbucket + 1

            if ( jup == jdest ) then  ! pass it by; it's in the right spot
               bucket_cur_pts(jdest) = bucket_cur_pts(jdest) + 1

            else  ! add a link to the chain
               orig_slots(nextmove) = listplace
               listplace = nextmove
               bucket_cur_pts(jdest) = bucket_cur_pts(jdest) + 1
               if ( bucket_cur_pts(jdest) > bucket_end_pts(jdest) ) then
                  done_bucket(jdest) = .true.
               end if
               jdest = jup

            end if

         end if

      end do

      ! execute the swap chain
      temp_vec(1:nelem) = list(1:nelem,chainstart)
      list(1:nelem,chainstart) = list(1:nelem,listplace)
      orig_slots(chainstart) = listplace
      backchain = orig_slots(listplace)
      do while (backchain /= chainstart)
         list(1:nelem,listplace) = list(1:nelem,backchain)
         listplace = backchain
         backchain = orig_slots(listplace)
      end do
      list(1:nelem,listplace) = temp_vec(1:nelem)

   end if

end do

end subroutine marble_sort

end module sort_tools