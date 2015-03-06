!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the x2f_mpi      !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine x2f_setmethod_mock( method, info_xf, rinfo_xf )
implicit none
integer, intent(in)    :: method
integer, intent(inout) :: info_xf(*)
real(kind(1.d0)), intent(inout)  :: rinfo_xf(*)

if ( method==1 ) then
   info_xf(1)  = 0
   rinfo_xf(1) = 0.3
elseif ( method==2 ) then
   info_xf(1)  = 0
   rinfo_xf(1) = 0.6
elseif ( method==3 ) then
   info_xf(1)  = 2
end if

end subroutine
