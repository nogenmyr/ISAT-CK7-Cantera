!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the x2f_mpi      !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module pasr_multi_subs

use MPI
use isat_rnu

implicit none

integer, parameter :: k_sp = kind(1.e0)

!------------------------------------------------------------------------------
type :: pasr_type


   real(k_sp), pointer :: f(:,:), fin(:,:), flstrm(:), xstrm(:), temp(:)
   real(k_sp)          :: t, dt, prout, prpar
   character(30)       :: op, op1
   integer             :: gpasr, cond, ip1, ip2, nstr

end type pasr_type

contains  !===========================================================

subroutine pasr_init( n, ncomp, nstr, gpasr, pasr )

!  initialize data structure pasr

integer, intent(in) :: n, ncomp, nstr, gpasr
type (pasr_type)    :: pasr

integer :: lus

character(30) :: blank, head, tail

  allocate( pasr%f(n,ncomp) )
  allocate( pasr%fin(ncomp,nstr) )
  allocate( pasr%flstrm(nstr) )
  allocate( pasr%xstrm(nstr) )
  allocate( pasr%temp(n) )

  pasr%xstrm = 0.
  pasr%t     = 0.
  pasr%dt    = 0.
  pasr%prout = 0.
  pasr%prpar = 0.
  pasr%gpasr = gpasr
  pasr%nstr  = nstr

  call rnused( gpasr )
  call rnuget( pasr%ip1, pasr%ip2 )  
  
  blank = repeat(' ',30)
  head  = blank
  head  = 'pasr'
  tail  = blank
  tail  = 'out'

  call isat_lu( lus )
  call isat_file_name( head, gpasr, -1, tail, pasr%op )
  open( lus, file = pasr%op, status = "replace" )
  close( lus )
  
  head = 'pdata'
  call isat_lu( lus )
  call isat_file_name( head, 0, -1, tail, pasr%op1 )

  return
end subroutine pasr_init

!===========================================================

    subroutine mix( dt, tmix, n, ncomp, ldf, f )

!  advance particle properties by mixing for time dt.
!  For odd i, particle i relaxes to particle i+1 on the time scale tmix.

!  input:
!       dt      - time interval
!       tmix	- reaction time scale
!       n	- number of particles
!       ncomp	- number of composition variables
!       ldf	- leading dimension of f
!	f	- composition

!  output:
!	f	- new composition

	implicit none
	integer, intent(in)             :: n, ncomp, ldf
	real(kind(1.e0)), intent(in)    :: dt, tmix
	real(kind(1.e0)), intent(inout) :: f(ldf,ncomp)

	integer          :: k, i, j
        real(kind(1.e0)) :: decay, dfj
        real(kind(1.d0)) :: exparg

!  check that n is even

        if( mod(n,2) /= 0 ) then
           write(0,*)'mix: called with odd n = ', n
           stop
        endif

        exparg   = 2. * dt / tmix
        decay    = 0.5 * ( 1. - exp( -exparg ) )

        do k = 1, ncomp
        do i = 1, n-1, 2
           j        = i + 1
           dfj      = decay * ( f(i,k) - f(j,k) )
           f(i,k)   = f(i,k) - dfj
	   f(j,k)   = f(j,k) + dfj
	end do
	end do

    return
    end subroutine mix

!============================================================================

    subroutine pickpr( npick, n, ncomp, ldf, f )

!  randomly select pairs and place them at the end of the array.
!  with probability 1/2, particles in a pair are commuted.

!  input:
!	npick	- number of pairs to be picked
!	n	- number of particles
!	ncomp	- number of compositions
!	ldf	- leading dimension of f
!	f	- array of particle compositions

!  output:
!	f	- reordered array

!  work
!	work(n)

	implicit none

	integer, intent(in) ::  npick, n, ncomp, ldf
    real(kind(1.e0)), intent(inout) :: f(ldf,ncomp)

	integer :: nh, ipick, iplast, ip, i, j, k, jl, il, jt
    real(kind(1.e0)) :: work(n), temp

!  check that n is even

        if( mod(n,2) /= 0 ) then
           write(0,*)'pickpr: called with odd n=', n
           stop
        endif

        nh = n / 2

!  check npick

        if( npick == 0  .or.  npick == nh ) then
           return
        elseif( npick < 0  .or.  npick  >  nh ) then
           write(0,*)'pickpr: called with bad npick=', npick
           stop
        endif

        call rnu( work(1:n) )

!  loop over pairs

        do ipick = 1, npick

!  select at random pair (i,j=i+1)

           iplast = nh + 1 - ipick
           ip     = 1 + work(ipick) * iplast
           j      = 2 * ip
           i      = j - 1

           jl     = 2 * iplast
           il     = jl - 1

!  commute pair or not at random

           if( work(ipick+npick) > 0.5 ) then
              jt = jl
              jl = il
              il = jt
           endif

!  commute (i,j) and (il,jl)

           do k = 1, ncomp
              temp     = f(i ,k)
              f(i,k)   = f(il,k)
              f(il,k)  = temp

              temp     = f(j ,k)
              f(j,k)   = f(jl,k)
	      f(jl,k)  = temp
	   end do
	end do

    return
    end subroutine pickpr

!============================================================================

     integer function inflow( pasr )

!  return index of stream of next inflowing particle

!  input:
!	pasr - data structure for reactor

	implicit none
    type (pasr_type) :: pasr

	integer :: i, infl
	real(kind(1.e0)) :: flsum, xmax

!  determine inflow

    flsum = 0.
    xmax  = 0.

    do i = 1, pasr%nstr
       pasr%xstrm(i) = pasr%xstrm(i) + pasr%flstrm(i)
       flsum    = flsum + pasr%flstrm(i)

       if( pasr%xstrm(i) > xmax ) then
          infl  = i
          xmax  = pasr%xstrm(i)  !  current max.
       endif
	end do

!  check flsum

    if( flsum <= 0. ) then
       write(0,*)'inflow: flsum=', flsum
       stop
    endif

    pasr%xstrm(infl) = pasr%xstrm(infl) - flsum
    inflow = infl

    return
    end function inflow
    
 !============================================================================
    
    subroutine set_nufm_fac( iproc, nufm_fac )
    !
    ! Set nonuniform factor, nufm_fac
    !
    implicit none
    integer, intent(in)    :: iproc
    integer, intent(inout) :: nufm_fac
    
    nufm_fac = iproc+1
    
    end subroutine set_nufm_fac
 
 !============================================================================
    
    subroutine rxn( mpicomm, press, ncomp, n, nl_pasr, pasr, stomby )
	
	! Advance particle properties by reaction for time dt.
	!
	! The procedure here is 
	!     [QT] -> pack xf array -> call ISAT_MP -> unpack xf array
	!
	! To use x2f_mpi effectively, the user needs to set qt(T.or.F), mode and
	! info_imp array. For detailed information about mode and info_imp, please
	! refer to subroutine ISAT_MP. 
	!
	! Input:
	! mpicomm - MPI communicator
	!  press  - pressure in Chemkin units
	!  ncomp  - number of composition variables
	!    n    - number of particles in each reactor
	! nl_pasr - number of reactors on each processor
	!  pasr   - data structure for reactor
    !
	! Parameters set by user:
    !   qt    - use qt or not
    !  mode   - integer determining the action to be taken(see ISAT_MP)
    !info_imp - integer array dimensioned 10 controlling x2f_mpi performance
    !
	! Output:
	!  pasr   - updated data after reaction step
	! stats_rxn_#.out - timing statistics 
	!
	
	implicit none
		
	real, parameter                 :: op_inc = 1.02
	real(kind(1.e0)), intent(in)    :: press, stomby
	integer, intent(in)             :: mpicomm, ncomp, n, nl_pasr
	type(pasr_type), intent(inout)  :: pasr(nl_pasr)
	
	integer                         :: n_tot, nv_un, n_qt, nx, nf, ldxf, lpasr, i
	integer                         :: p(n*nl_pasr), info_imp(10)
	real(kind(1.e0))                :: cc(ncomp), ct(ncomp), dpt(3), dt
	real(kind(1.d0)), allocatable   :: xf(:,:)

    character(30)                   :: head, tail
    character(30), save             :: filename
    
    logical, save                   :: qt 
    integer, save                   :: istep=1, next_op=1, mode, adapt, nSR2, lu_save
    integer                         :: lu, iproc, ierr, wall_rate, wall_max, &
                                       wall_b, wall_e
    real                            :: cpu_b, cpu_e, cpu_qt, cpu_isatmp, &
                                       wall_qt, wall_isatmp
    real, save                      :: wall_cum=0.0, cpu_cum=0.0
    logical                         :: exist
    
    data qt    / .false. /
    data mode  /   1     /
    data nSR2  /   1     /
    data adapt /   0     /
    
    namelist / rxnnml / qt, mode, adapt, nSR2

    if ( istep==1 ) then
       !--- read namelist (if the file pasr.nml exists)

       inquire( file = 'rxn.nml', exist = exist )
       if ( exist ) then
          call isat_lu( lu )
          open( lu, file = 'rxn.nml', position ='rewind',action ='read')
          read( lu, nml = rxnnml )
          close( lu )
       endif
    
       ! Generate the different names for output files on different processors
	   call mpi_comm_rank( mpicomm, iproc, ierr )
	   head = 'stats_rxn'
	   tail = 'out'
	   call isat_file_name(head,iproc,-1,tail,filename)
	   call isat_lu( lu_save )
       open( lu_save, file=filename, status='replace', action = 'write' )
	   
	   ! Initialize ISAT
	   cc(1:ncomp) = pasr(1)%f(1,1:ncomp)
	   dt          = pasr(1)%dt
	   dpt(2)      = press
	   call scirxn( dt, ncomp, cc, ct, dpt )
	   
	   ! Initialize ISAT_MP
	   info_imp = 0
	   if ( qt ) info_imp(1) = 1
	   info_imp(3) = adapt
	   info_imp(4) = 0
	   
	end if
	
	call MPI_Barrier( MPI_COMM_WORLD, ierr )
    n_tot = n*nl_pasr
    ! Measure wall clock time spent in QT mode
    call system_clock( count=wall_b )  
    ! Measure CPU time spent in QT mode
    call cpu_time( cpu_b )
	p(1:n_tot) = 0
	n_qt = 0
	
	if ( qt ) then    ! invoke QT mode if qt=.true.
	   nv_un = 0      ! counter used to count the number of unevaluated particles
	   call x2f_setmethod(2)   ! set method for QT mode: 
	                           ! 1 - PR(primary retrieve)
	                           ! 2 - PR+SR(primary and secondary retrieve)
	   do lpasr = 1, nl_pasr
	      do i = 1, n
	         cc(1:ncomp) = pasr(lpasr)%f(i,1:ncomp)
	         dt          = pasr(lpasr)%dt
	         dpt(2)      = press
	         call scirxn( dt, ncomp, cc, ct, dpt )
	         if ( dpt(1)>0.0 ) then   ! based on the value of returned density,
	                                  ! determine if the particle has been evaluted
	                                  ! or not
	            pasr(lpasr)%f(i,1:ncomp) = ct(1:ncomp)
	            pasr(lpasr)%temp(i)      = dpt(3)
	            p(i+n*(lpasr-1))         = 1
	            n_qt = n_qt + 1
	         else
	            nv_un = nv_un + 1
	         end if
	      end do
	   end do
	else
	   nv_un = n_tot
	end if
	
	call MPI_Barrier( MPI_COMM_WORLD, ierr )
	call cpu_time( cpu_e )
    call system_clock( count=wall_e, count_rate=wall_rate, count_max=wall_max )
    if ( wall_e<wall_b ) wall_e=wall_e+wall_max
    wall_qt = real(wall_e-wall_b)/real(wall_rate)
    cpu_qt  = cpu_e-cpu_b
	   
	nx   = ncomp + 2
	nf   = ncomp + 3
	ldxf = max(nx,nf)

	allocate( xf(ldxf,nv_un) )	
    ! Pack xf array
	nv_un = 0
	do lpasr = 1, nl_pasr
	   do i = 1, n
	      if ( p(i+n*(lpasr-1))==0 ) then
	         nv_un = nv_un + 1
	         xf(1:ncomp,nv_un) = pasr(lpasr)%f(i,1:ncomp)
	         xf(ncomp+1,nv_un) = pasr(lpasr)%dt
	         xf(ncomp+2,nv_un) = press
	      end if
	   end do
	end do
		
	call ISAT_MP( nx,nf,ldxf,nv_un,xf,mode,nSR2,info_imp,istep )
	            
	! Unpack xf array
	nv_un = 0
	do lpasr = 1, nl_pasr
	   do i = 1, n
	      if ( p(i+n*(lpasr-1))==0 ) then
	         nv_un = nv_un + 1
	         pasr(lpasr)%f(i,1:ncomp) = xf(1:ncomp,nv_un)
	         pasr(lpasr)%temp(i)      = xf(ncomp+3,nv_un)
	      end if
	   end do
	end do
	
	deallocate( xf )
	call MPI_Barrier( MPI_COMM_WORLD, ierr )
	call cpu_time( cpu_b )
	call system_clock( count=wall_b, count_rate=wall_rate, count_max=wall_max )
	if ( wall_b<wall_e ) wall_b = wall_b+wall_max
	wall_isatmp = real(wall_b-wall_e)/real(wall_rate)
    cpu_isatmp  = cpu_b-cpu_e
    wall_cum = wall_cum + wall_isatmp + wall_qt
    cpu_cum  = cpu_cum + cpu_isatmp + cpu_qt
    
    if ( istep==next_op ) then
       next_op = max( istep+1.0, istep*op_inc )
       next_op = istep+1.0
       !call isat_lu( lu )
       !open(lu, file=filename, position='append')
       write(lu_save,'(3i8, 100e16.6)') &
             istep, n_tot, n_qt, cpu_qt, wall_qt, cpu_isatmp, wall_isatmp, cpu_cum, wall_cum
    end if
    
    istep = istep + 1
	
	end subroutine 
	
!============================================================================
	
	subroutine rxn_bk( mpicomm, press, ncomp, n, nl_pasr, pasr )
	
	! Advance particle properties by reaction for time dt.
    !
    ! This reaction subroutine considers the case which has large number of 
    ! particles during each call. After qt step, the particle array is first 
    ! divided into several blocks; then by calling ISAT_MP for each block, all of
    ! the particles get evaluated. 
	! The procedure here is  
	!     [QT] -> divide particle array into blocks -> loop over blocks
	! For each block
	!     pack xf array -> call ISAT_MP -> unpack xf array
	!
    ! To use x2f_mpi effectively, the user needs to set qt(T.or.F), bksize(which 
    ! control the size of each block), mode and info_imp array. For detailed 
    ! information about mode and info_imp, please refer to subroutine ISAT_MP.
    !
	! Input:
    ! mpicomm - MPI communicator
	!  press  - pressure in Chemkin units
	!  ncomp  - number of composition variables
	!    n    - number of particles in each reactor
	! nl_pasr - number of reactors on each processor
	!  pasr   - data structure for reactor
	!
    ! Parameters set by user:
    !   qt    - use qt or not
    !  bksize - block size in megabytes allowed for using ISAT_MP 
    !  mode   - integer determining the action to be taken(see ISAT_MP)
    !info_imp - integer array dimensioned 10 controlling x2f_mpi performance
    !
	! Output:
	!  pasr   - updated data after reaction step
	!stats_rxn_#.out - timing statistics
	!
	
	implicit none
		
	real, parameter                 :: op_inc = 1.0
	real(kind(1.e0)), intent(in)    :: press
	integer, intent(in)             :: mpicomm, ncomp, n, nl_pasr
	type(pasr_type), intent(inout)  :: pasr(nl_pasr)
	
	integer                         :: n_tot, nv_un, n_qt, nx, nf, ldxf, lpasr, i 
	integer                         :: p(n*nl_pasr), info_imp(10)
	real(kind(1.e0))                :: cc(ncomp), ct(ncomp), dpt(3), dt
	real(kind(1.d0)), allocatable   :: xf(:,:)
	
	integer                         :: nb, nbmax, np, np_lastB, nv_tmp, ib
	integer                         :: np_tmp, i0, lpasr0
	
	character(30)                   :: head, tail
    character(30), save             :: filename
    
    logical, save                   :: qt = .true.
    real, save                      :: bksize = 1.0
    integer, save                   :: istep=1, next_op=1, mode=3
    integer                         :: lu, ierr, iproc, wall_rate, wall_max, &
                                       wall_b, wall_e
    real                            :: cpu_b, cpu_e, cpu_qt, cpu_isatmp, &
                                       wall_qt, wall_isatmp
                                       
    if ( istep==1 ) then
       ! Generate the different names for output file on different processors
	   call mpi_comm_rank( mpicomm, iproc, ierr )
	   head = 'stats_rxn'
	   tail = 'out'
	   call isat_file_name(head,iproc,-1,tail,filename)
	   
	   ! Initialize ISAT
	   cc(1:ncomp) = pasr(1)%f(1,1:ncomp)
	   dt          = pasr(1)%dt
	   dpt(2)      = press
	   call scirxn( dt, ncomp, cc, ct, dpt )
	end if

	n_tot = n*nl_pasr
    call system_clock( count=wall_b )
    call cpu_time( cpu_b )
	p(1:n_tot) = 0
	n_qt = 0
	
	if ( qt ) then    ! invoke QT mode if qt=.true.
	   nv_un = 0      ! counter used to count the number of unevaluated particles
	   call x2f_setmethod(2)   ! set method for QT mode: 
	                           ! 1 - PR(primary retrieve)
	                           ! 2 - PR+SR(primary and secondary retrieve)  
	   do lpasr = 1, nl_pasr
	      do i = 1, n
	         cc(1:ncomp) = pasr(lpasr)%f(i,1:ncomp)
	         dt          = pasr(lpasr)%dt
	         dpt(2)      = press
	         call scirxn( dt, ncomp, cc, ct, dpt )
	         if ( dpt(1)>0.0 ) then   ! based on the value of returned density,
	                                  ! determine if the particle has been evaluted
	                                  ! or not  
	            pasr(lpasr)%f(i,1:ncomp) = ct(1:ncomp)
	            pasr(lpasr)%temp(i)      = dpt(3)
	            p(i+n*(lpasr-1))         = 1
	            n_qt = n_qt + 1
	         else
	            nv_un = nv_un + 1
	         end if
	      end do
	   end do
	else
	   nv_un = n_tot
	end if
	
	call cpu_time( cpu_e )
    call system_clock( count=wall_e, count_rate=wall_rate, count_max=wall_max )
    if ( wall_e<wall_b ) wall_e=wall_e+wall_max
    wall_qt = real(wall_e-wall_b)/real(wall_rate)
    cpu_qt  = cpu_e-cpu_b

	nx   = ncomp + 2
	nf   = ncomp + 3
	ldxf = max(nx,nf)
	
	call ISAT_MP_bksize( mpicomm, nv_un, ldxf, bksize, nb, nbmax, np )
	print *, 'np, nb, nbmax, np_lastb:', np, nb, nbmax, nv_un-(nb-1)*np

	allocate( xf(ldxf,np) )	
	i      = 1
	i0     = 1
	lpasr  = 1
	lpasr0 = 1
	do ib = 1, nbmax
	   if ( ib<nb ) then
	      nv_tmp = np
	   elseif ( ib==nb ) then
	      nv_tmp = nv_un - (nb-1)*np 
	   else
	      nv_tmp = 0
	   end if 
	   
	   np_tmp = 0
	   do while( np_tmp<nv_tmp.and.lpasr<=nl_pasr )
	      if ( i>n ) then
	         i     = 1
	         lpasr = lpasr + 1
	      else
	         if ( p(i+n*(lpasr-1))==0 ) then
	            np_tmp = np_tmp + 1
	            xf(1:ncomp,np_tmp) = pasr(lpasr)%f(i,1:ncomp)
	            xf(ncomp+1,np_tmp) = pasr(lpasr)%dt
	            xf(ncomp+2,np_tmp) = press
	         end if
	         i = i + 1
	      end if
	   end do
	   
       info_imp = 0
       info_imp(1) = 1
	   call ISAT_MP( nx,nf,ldxf,nv_tmp,xf,mode,info_imp )
	   
	   np_tmp = 0
	   do while( np_tmp<nv_tmp.and.lpasr0<=nl_pasr )
	      if ( i0>n ) then
	         i0     = 1
	         lpasr0 = lpasr0 + 1
	      else
	         if ( p(i0+n*(lpasr0-1))==0 ) then
	            np_tmp = np_tmp + 1
	            pasr(lpasr0)%f(i0,1:ncomp) = xf(1:ncomp,np_tmp) 
	            pasr(lpasr0)%temp(i0)      = xf(ncomp+3,np_tmp)
	         end if
	         i0 = i0 + 1
	      end if
	   end do
	end do
	
	deallocate( xf )
	
	call cpu_time( cpu_b )
	call system_clock( count=wall_b, count_rate=wall_rate, count_max=wall_max )
	if ( wall_b<wall_e ) wall_b = wall_b+wall_max
	wall_isatmp = real(wall_b-wall_e)/real(wall_rate)
    cpu_isatmp  = cpu_b-cpu_e
    
    if ( istep==next_op ) then
       next_op = max( istep+1.0, istep*op_inc )
       call isat_lu( lu )
       open(lu, file=filename, position='append')
       write(lu,'(3i8, 4e16.6)') &
             istep, n_tot, n_qt, cpu_qt, wall_qt, cpu_isatmp, wall_isatmp
       close(lu)
    end if
    
    istep = istep + 1
	
	end subroutine 

end module pasr_multi_subs
