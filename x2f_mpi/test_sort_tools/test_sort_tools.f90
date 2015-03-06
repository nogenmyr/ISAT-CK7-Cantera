!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the x2f_mpi      !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

program sort_test

use sort_tools

implicit none

integer, parameter   :: nlist = 100, ld = 1, nelem = 1, nbucket = 25
real(kind(1.d0))     :: list(ld,nlist), listsave(ld,nlist)
integer              :: bucket_assigns(nlist)
integer              :: bucket_totals(nbucket+1)
integer              :: list_origplaces(nlist)
integer              :: orig_slots(nlist), seedsize, i, buck, buckend
integer, allocatable :: seednum(:)
real                 :: x1
character(9)         :: header
real(kind(1.d0))     :: tempitem(nelem)
integer              :: thisone, nextone
integer              :: sorted_slots(nlist)

bucket_totals = 0
call random_seed( size = seedsize )
allocate( seednum(seedsize) )
seednum = (/ ( 4550 * i, i = 1, seedsize ) /)
call random_seed( put = seednum(1:seedsize) )
do i = 1, nlist
   call random_number( x1 )
   list(1,i) = x1
   buck = int( 10.0 * x1 )
   buck = ( buck + 1 ) * ( buck + 3 ) / 4
   ! want buckets 0 and 30 to be out of range 1:nbucket == 1:25
   ! also want lots of buckets with zero totals for testing purposes
   bucket_assigns(i) = buck
   if ( ( buck > 0 ) .and. ( buck <= nbucket ) ) then
      bucket_totals(buck) = bucket_totals(buck) + 1
   else
      bucket_totals(nbucket+1) = bucket_totals(nbucket+1) + 1
   end if
end do

listsave = list
call bucket_sort( nlist, ld, nelem, list, bucket_assigns, &
      nbucket, bucket_totals, list_origplaces, orig_slots)
!call marble_sort( nlist, nelem, nelem, list, bucket_assigns, &
!      nbucket, bucket_totals, orig_slots)

buckend = 0
do i = 1, nbucket + 1
   buck = buckend + 1
   buckend = buckend + bucket_totals(i)
   if ( i <= nbucket ) then
      write (header,"('bucket',i3)") i
   else
      write (header,"('leftovers')")
   end if
   if ( bucket_totals(i) == 0 ) then
      print "('bucket',i3,':  total =  0')", i
   else
      print "(a9,':  total =',i3,',  min =',f12.9,          &
              ',  max =',f12.9)", header, bucket_totals(i), &
              minval( list(1:nelem,buck:buckend), dim=2 ),  &
              maxval( list(1:nelem,buck:buckend), dim=2 )
   end if
end do

call pigeonhole_sort( nlist, ld, nelem, list, orig_slots, sorted_slots)

print "('error sum from restoring list order was',e23.15)", sum( abs ( &
        list(1:nelem,1:nlist) - listsave(1:nelem,1:nlist) ) )

end program sort_test

