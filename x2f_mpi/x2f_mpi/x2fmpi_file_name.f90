!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the x2f_mpi      !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine x2fmpi_file_name( head, n, p, tail, name )

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

end subroutine x2fmpi_file_name