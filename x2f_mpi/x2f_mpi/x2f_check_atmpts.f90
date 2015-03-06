!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the x2f_mpi      !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine x2f_check_atmpts(gat,n_spec1,info,igat,nproc)
!
! Check to see if the number of attempts is valid or not!
!
integer, intent(in)  :: gat, n_spec1, info(gat,n_spec1), igat, nproc

integer  :: nolocal, npref, i

npref = 0
do i = igat, gat
   if ( info(i,2)==4 ) then
      nolocal = info(igat,3)
      if ( npref==0 ) then 
         npref = npref + 1
      elseif ( npref>0.and.info(i,6)==0.and.info(i-1,2)==4 ) then
         npref = npref + 1
      else
         Exit
      end if
   end if
end do
if ( (nolocal==0.and.npref>nproc).or.(nolocal==1.and.npref>nproc-1) ) then
   print *, 'npref & nproc:', npref, nproc
   stop 'Err: The number of PREF attempts is larger than processor involved!'
end if

end subroutine x2f_check_atmpts
   

