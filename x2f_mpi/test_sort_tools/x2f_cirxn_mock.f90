!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the x2f_mpi      !
!  source directory.)                                                 !
!---------------------------------------------------------------------!


subroutine x2f_cirxn_mock( ldxf, nx, x, f, k_pos, info_xf, rinfo_xf )

implicit none
integer, intent(in)             :: nx, ldxf, k_pos, info_xf(*)
real(kind(1.d0)), intent(in)    :: x(nx), rinfo_xf(*)
real(kind(1.d0)), intent(inout) :: f(ldxf)

real :: temp

f = 0.0
if ( info_xf(1)==1 ) then
   return
elseif ( info_xf(1)==2 ) then
   f(2) = x(1)
   f(1) = x(2)  
   if ( ldxf>=3 ) f(3:ldxf) = x(1) + 1.0 
elseif ( info_xf(1)==0 ) then
   call random_number( temp )
   if ( temp<rinfo_xf(1) ) then
      f(2) = x(1)
      f(1) = x(2)
      if ( ldxf>=3 ) f(3:ldxf) = x(1) + 1.0 
   else
      f(1:nx) = x(1:nx)
      f(k_pos) = -x(k_pos)
   end if
end if

end subroutine 

