!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the x2f_mpi      !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine x2f_setmethod( method, info_xf, rinfo_xf )

implicit none
integer, intent(in)     :: method
integer, intent(inout)  :: info_xf(*)
real(kind(1.d0)),intent(inout) :: rinfo_xf(*)

integer :: info(100)
real(kind(1.d0)) :: rinfo(50), stats(100)

info     = -12345
rinfo    = -12345.
info(80) = 1

if ( method==1 ) then
   !Invoke Primary Retrieve
   info(23)    = 0
   info(24:27) = 1
elseif ( method==2 ) then
   !Invoke Secondary Retrieve
   info(23:24) = 0
   info(25:27) = 1
else
   !Invoke Direct Evaluation
   info(23:27) = 0
end if
call cisat( 1, info, rinfo, stats )

end subroutine x2f_setmethod
