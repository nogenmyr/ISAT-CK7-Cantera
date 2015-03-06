!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module isat_val

   use isat_rnu
   use isat_abort_m

   implicit none

   integer, parameter :: nin = 5, nk = 2*nin, coderet = 9632468, &
                         keyin5 = 255255

   contains

   subroutine isat_lic( luout, retcode )

   implicit none

   integer, intent(in)  :: luout
   integer, intent(out) :: retcode
   integer :: v(8), lu, i, keyin(nk), keyout(nk), iy, im, id, &
              idate, days_left, expire, exit_status, ip(4)
   logical :: exist
   logical, save :: validated = .false.
   real(kind(1.d0)) :: ipstart, ipend, iphost

   character(len=128) :: &
   comm0 = 'cp isat.key isatkey.loc1', &
   comm0w= 'COPY /Y isat.key isatkey.loc1 1>NUL:', &
   comm1 = 'cp ${ISATKEY} isatkey.loc1', &
   comm1w= 'COPY %ISATKEY% isatkey.loc1 1>NUL:', &
   comm2 = 'grep `hostname` /etc/hosts | sed ''s/\./ /g'' > isatkey.loc2', &
   comm2a= 'grep `hostname` /etc/hosts | sed ''s/\\./ /g'' > isatkey.loc2', &
   comm2w= 'ipconfig | grep ''IP A'' | grep -v ''n IP'' | sed ''s/.*:/ /g'' | sed ''s/\./ /g'' > isatkey.loc2', &
   comm3 = 'rm -f isatkey.loc*', &
   comm3w= 'DEL /Q isatkey.loc*'
   character(len=6) :: systype

   interface
      function isat_sys( command )
         implicit none
         character(128), intent(in) :: command
         integer isat_sys
      end function isat_sys
   end interface

   retcode = coderet
   validated = .true.  !laniu 08-23-07 check ISAT keyXXXX
   if( validated ) return

!  set system-dependent commands
   call isat_systype( systype )

   if( systype == 'WIN_NT' ) then
      comm0 = comm0w
      comm1 = comm1w
      comm2 = comm2w
      comm3 = comm3w
   elseif( systype == 'AIX___' )  then
      comm2 = comm2a
   endif

! use local file isat.key if it exists,
! otherwise use environment variable ISATKEY

   inquire( file = 'isat.key', exist = exist ) 
   if( exist ) then
      exit_status = isat_sys(comm0)
      if( exit_status /=0 ) call isat_abort( 'isat_license', 0, &
         mess = 'Failed to make copy of isat.key' )
   else
      exit_status = isat_sys(comm1)
      if( exit_status /=0 ) call isat_abort( 'isat_license', 1, &
         mess = 'License key not found: check setting of ISATKEY' )
   endif

   call isat_lu( lu )
   open( lu, file='isatkey.loc1', err=100 )

   do i = 1, nk
      read(lu,*,err=120,end=140) keyin(i)
   end do
   close(lu)

   call lic_key( keyin, keyout )

   do i = 1, nk
      if( keyin(i) /= keyout(i) ) call isat_abort( 'isat_license', &
          2, mess='License key is not valid.')
   end do

   call date_and_time( values=v )
   iy = v(1)
   im = v(2)
   id = v(3)
   idate = milldays( iy, im, id )

   im = keyin(1)/10000
   id = (keyin(1)-10000*im)/100
   iy = (keyin(1)-10000*im-100*id) + 2000
   expire = milldays( iy, im, id )

   days_left = expire - idate

   if( days_left < 0 ) then
      call isat_abort( 'isat_license', 4, &
         mess='License has expired.' )
   elseif( days_left < 60 ) then
      if( luout >= 0 ) write(luout,200) days_left
200   format('***WARNING*** License expires in approximately ',i2,' days.' )
   endif

   if( keyin(5) /= keyin5 ) then
      exit_status = isat_sys(comm2)
      if( exit_status /=0 ) call isat_abort( 'isat_license', 3, &
         mess = 'Cannot determine IP address' )

      open( lu, file='isatkey.loc2', err=160 )
90    continue
      read(lu,*,end=180,err=90) ip
      close(lu)

      ipstart = keyin(2) * 1000000.d0 + keyin(3)
      ipend   = keyin(4) * 1000000.d0 + keyin(5)
      iphost  = 1.d9 * ip(1) + 1.d6 * ip(2) + 1.d3 * ip(3) + ip(4)

      if( iphost+.5 < ipstart  .or.  iphost-.5 > ipend ) then
         if( luout >= 0 ) then
            write(lu_err,'(e25.13)') ipstart
            write(lu_err,'(e25.13)') iphost 
            write(lu_err,'(e25.13)') ipend   
         endif
         call isat_abort( 'isat_license', 5, mess='Invalid IP address' )
      endif
   endif

   exit_status = isat_sys(comm3)

   validated = .true.
   if( luout >= 0 ) write(luout,*)'License validated'

   return

100   call isat_abort( 'isat_license', 6, &
               mess='Error opening ISATKEY.' )
120   call isat_abort( 'isat_license', 7, &
               mess='Error reading ISATKEY.' )
140   call isat_abort( 'isat_license', 8, &
               mess='Hit end of file  ISATKEY.' )
160   call isat_abort( 'isat_license', 9, &
               mess='Error opening IP file.' )
180   call isat_abort( 'isat_license', 10, &
               mess='Error reading IP file.' )

   end subroutine isat_lic

   integer function milldays( iy, im, id )

   implicit none
   integer :: iy, im, id

   milldays = id + 31* ( (im-1) + 12*(iy-2000) )

   return
   end function milldays

   subroutine lic_key( keyin, keyout ) 

   implicit none

   integer, intent(in)  :: keyin(:)
   integer, intent(out) :: keyout(:)

   integer :: is1, is2, ig1, ig2, i, ii
   integer :: jx(nk)

   if( size(keyin) /= nk  .or.  size(keyout) /= nk ) &
       call isat_abort( 'isat_license', 11, mess='Invalid call' )

   call rnuget( ig1, ig2 )

   do i = 1, nin

      is1 = keyin(i)
      is2 = 2147483560 - is1

      call rnuput( is1, is2 )
      do ii = 1, i
         call rnu( jx(i:nk) )
      end do
      keyout(i) = keyin(i)
   end do

   call rnuput( ig1, ig2 )

   keyout(nin+1:nk) = jx(1:nk-nin)

   return
   end subroutine lic_key

end module isat_val
