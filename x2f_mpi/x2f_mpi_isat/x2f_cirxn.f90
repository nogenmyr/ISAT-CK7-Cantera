!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the x2f_mpi      !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine x2f_cirxn( ldxf, nx, x, f, k_pos, info_xf, rinfo_xf )

!  x2f routine for Chemistry Interface (CI) 
!  x = { c(1:ncv), dt, press  }
!  f = { c(1:ncv), dens, press, temp}

implicit none
integer, parameter      :: k_dp = kind(1.d0)

integer,    intent(in)    :: ldxf, nx, k_pos, info_xf(*)
real(k_dp), intent(in)    :: x(nx), rinfo_xf(*)
real(k_dp), intent(out)   :: f(ldxf)

integer       :: ncv
real(k_dp)    :: dpt(3), dt

if ( info_xf(1)==-1 ) then
   ! measure message passing time, do not evaluate function
   f(1:ldxf) = x(1)
else 
   !---------  (try to) evaluate f(x) -------
   dpt    = 0.d0
   ncv    = nx - 2
   dt     = x(ncv+1)
   dpt(2) = x(ncv+2)

   call cirxn( dt, ncv, x(1:ncv), f(1:ncv), dpt )
   
   if ( dpt(1) > 0.d0 ) then
      f(ncv+1) = dpt(1)
      f(ncv+2) = dpt(2)
      f(ncv+3) = dpt(3)  
   else
      f(1:nx) = x(1:nx)
      f(k_pos) = -x(k_pos)
      if ( ldxf>nx ) then
         f(nx+1:ldxf) = 0.0
      end if
   end if
end if

end subroutine x2f_cirxn
