!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine isat_lu( lun )

!  return number (between luf and lul) of an avialable logical unit
!  call isat_lu_set( luf, lul ) to change the values of luf and lul.

   use isat_lu_m
   use isat_abort_m

   implicit none

   integer, intent(out) :: lun
   logical              :: usedio
   integer              :: lus

   if( lu < luf  .or.  lu > lul ) lu = luf

   lun = lu
   lus = lu

   do
      inquire( lun, opened = usedio )
      if( .not.  usedio ) exit
      lun = lun + 1

      if( lun == lus ) call isat_abort('isat_lu',1, &
	     mess='no unopened logical units between luf and lul', &
		 ivar=(/ luf, lul /) )

      if( lun > lul ) lun = luf
   end do

   lu  = lun + 1

   return
end subroutine isat_lu  !------------------------------------------------

subroutine isat_flush(lu)

!  routine to simulate standard unix flush.
!  for lu > 0, and if if_flush is .true., isat_flush flushes the buffer
!      of logical unit lu.
!  by default, if_flush is .false.
!  a call with lu = -1 sets if_flush = .true.
!  a call with lu = -2 sets if_flush = .false.

!  flush operation:
!  1/	determine file name
!  2/	close file
!  3/	open file
!
    use isat_abort_m

	implicit none

	integer, intent(in) :: lu
	integer             :: recl
	logical, save       :: if_flush = .false.
	logical             :: exist, opened
	character(100)      :: flname

	if( lu == -1 ) then
	   if_flush = .true.
	   return
	elseif( .not.if_flush ) then
	   return
	elseif( lu == -2 ) then
	   if_flush = .false.
	elseif( lu < 1 ) then
	   call isat_abort('isat_flush',1,mess='bad value, lu = ', isv=lu)
	endif

!  flush buffer of logical unit lu

	inquire( lu, name = flname, exist = exist, opened = opened, recl = recl, err=100 )
	if( .not.exist  .or.  .not.opened ) return
	close( lu )
	open ( lu, file = flname, position = 'append', action = 'write', &
	        recl = recl, err = 200 )
	return

100   call isat_abort('isat_flush',2,mess='inquire error for lu = ', isv=lu)
200   call isat_abort('isat_flush',3,mess='open error for lu = ', isv=lu)

end subroutine isat_flush  !--------------------------------------------

subroutine isat_file_name( head, n, p, tail, name )

!  construct the character variable name:
!
!     name = head_n_p.tail    (for n>=0 and p>=0)
!
!  for n<0, '_n' is omitted
!  for p<0, '_p' is omitted.
 

   implicit none

   character(30), intent(in) :: head, tail
   integer, intent(in)       :: n, p
   character(30), intent(out):: name

   character(30):: num_n, num_p, dot, under, blank, name1, name2

   blank = repeat(' ',30)
   dot   = blank
   dot   = '.'
   under = blank
   under = '_'

   write(num_n,"(i30)") n
   write(num_p,"(i30)") p
   
   name2 = head

   if( n >= 0 ) then
      call concat( name2, under, name1 )  ! name1 = head_
	  call concat( name1, num_n, name2 )  ! name2 = head_n
   endif

   if( p >= 0 ) then
      call concat( name2, under, name1 )  ! name1 = head_n_
	  call concat( name1, num_p, name2 )  ! name2 = head_n_p
   endif

   call concat( name2, dot,  name1)       ! name1 = head_n_p.
   call concat( name1, tail, name )       ! name  = head_n_p.tail

   return

contains

   subroutine concat( head, tail, joined )

!  concatinate head and tail to form joined

   character(30), intent(in)   :: head, tail
   character(30), intent(out)  :: joined

   character(30) :: he, ta
   integer :: lh, lt

   joined = repeat(' ',30)

   he = adjustl(head)
   ta = adjustl(tail)

   lh = len_trim(he)
   lt = len_trim(ta)

   joined(1:lh)         = he
   joined(lh+1:lh+lt)   = ta

   return
   end subroutine concat

end subroutine isat_file_name