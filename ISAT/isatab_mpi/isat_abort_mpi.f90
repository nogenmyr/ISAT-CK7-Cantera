!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module isat_abort_m

integer, save :: lu_err = 0

contains

subroutine isat_abort( sub, loc, mess, chv, isv, ivar, rsv, rvar )

!  Aborts ISATAB because of error: prints diagnostic and stops

   use mpi
   use isat_prec
   implicit none

   character(*), intent(in)           :: sub     ! name of calling routine
   integer, intent(in)                :: loc     ! location number
   character(*), intent(in), optional :: mess    ! message
   character(*), intent(in), optional :: chv     ! character variable
   integer, intent(in), optional      :: isv, ivar(:) ! integer variables
   real(k_xf), intent(in), optional   :: rsv, rvar(:) ! real variables

   integer :: val, ierr
   logical :: flag

   write(lu_err,*)' '
   write(lu_err,*)'********** ISAT_ABORT  **************'
   write(lu_err,*)' '
   write(lu_err,*)'routine  = ', sub
   write(lu_err,*)'location = ', loc
   if( present(mess) ) write(lu_err,*)'message  = ', mess
   if( present(chv ) ) write(lu_err,1) chv 
   if( present(isv ) ) write(lu_err,2) isv 
   if( present(ivar) ) write(lu_err,2) ivar
   if( present(rsv ) ) write(lu_err,3) rsv  
   if( present(rvar) ) write(lu_err,3) rvar 
   write(lu_err,*)' '
   write(lu_err,*)'****** END ISAT_ABORT  **************'
   write(lu_err,*)' '

   call MPI_INITIALIZED( flag, ierr )
   if( flag ) then
      call MPI_ABORT( MPI_COMM_WORLD, val, ierr )
   else
      stop
   endif

1     format((a))
2     format((8i10))
3     format((1p,5e13.4))
end subroutine isat_abort

end module isat_abort_m
