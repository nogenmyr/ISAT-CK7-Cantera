!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the x2f_mpi      !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

program pasr_multi
!
! PaSR - partially-stirred reactor
!
! Perform multiple PaSR calculations to test x2f_mpi and ISAT-CK. This version 
! set up to mimic premixed CH4 flame calculation.
!
! There are a total of ng_pasr reactors, which have n_cond distinctly different
! operating conditions.(ng_pasr is required to be an integer multiple of n_cond.) 
! Thus there are ng_pasr/n_cond statistically identical reactors with each of the
! operating conditions.
!
! The operating conditions are defined by the compositions and flow rates of the
! nstr inflowing streams.
!
! There are n_proc processes, where it is required that ng_pasr be an integer 
! multiple of n_proc. In each process there are nl_pasr = ng_pasr/nproc reactors.

! Each reactor consists of n particles. Particles flow in and out at a rate 
! determined by the residence time tres. There are nstr inflowing streams, with 
! relative mass flow rates flstrm. The initial condition for all particles is set
! to the composition of the first stream. The particles react; and they mix on a 
! timescale tmix with partners which are randomly re-assigned on a timescale tpair.
!
! Physical parameters:
!	tres	- residence time (s)
!	tmix	- mixing time scale (s)
!	tpair	- pairing time scale (s)
!	anres	- length of run (number of residence times)
!	flstrm  - relative flow rates of streams
!	pr_atm	- pressure (atm.)
!
! The largest time step dtmax is determined by the specified constant cdtrp. 
! The smallest time step dtmin is smaller than dtmax by the factor dtfac(i.e., 
! dtmin = dtmax * dtfac). The time step used is random,uniformly distributed 
! between dtmin and dtmax.
!
! Numerical parameters:
!	n	    - number of particles (must be even) on each reactor
!             can be uniform/non-uniform
! nufm_fac  - nonuniform factor
!           = 1, uniform number of particles on each reactor, n=n0
!           > 1, reactors on some processor have more particles, n=n0*nufm_fac
!	cdtrp	- time-step control (must be less than 0.25)
!	dtfac	- min-to-max substep ratio
!	cdtmix	- sub-step control, based on tmix (less than 0.5)
!
! ISAT parameters:
!   etola   - ISAT error tolerance 
!   stomby  - ISAT storage(Mbyte)
!  const_dt - constant or varying time step
!
! Other parameters:
!   rstr    - start from scratch or restart from existing data
!   op_inc  - output increment
!  savefrqc - frequency(e.g, every 1000 steps) to save particle information
!
! Pre-processing:
!   ckinterp.exe must be run with chem.inp and therm.dat as input files, to 
!   produce chem.bin and chem.out
!
! Input files:
!   chem.bin   - required by ISAT-CK
!	streams.in - required by ISAT-CK (provides value of nstr)
!   pasr.nml, ci.nml and isat.nml are used if they exist.
!
! Output files
!	pasr_#.out - t, T(1), T_mean  from pasr with gpasr #
!  pdata_#.out - particle information used as restart file
!stats_rxt_#.out - timing statistics output from reaction subroutine
!  isat_mp.op  - the strategy used in x2f_mpi output from subroutine ISAT_MP 

use MPI
use pasr_multi_subs
implicit none

type (pasr_type), allocatable :: pasr(:)
integer, allocatable :: g_cond(:)
real, allocatable    :: ftemp(:),fstr(:,:),dfin(:),fin(:,:,:),flstrm(:,:),tempf(:)	                    
integer :: nproc, iproc, ng_pasr, n_cond, ci_info(20), info(100), gpasr, ncomp, &
           nfull, nstr, ns, nl_pasr, kf_gpasr, kl_gpasr, mode_cond, ist, lpasr, &
           cond, i, istep, nprout, npick, istrm, nprpar, ifst, k, ii, luerr, &
           luout, lu_op, lu, ierror, next_op, init_stream, const_dt, nufm_fac, &
		   n0, n, savefrqc, rcd, mpicomm, kase, j
real(kind(1.e0)) :: ci_rinfo(20), rinfo(50), dpt(3), patm, pr_atm, press, tres, &
                    tmix, tpair, cdtrp, dtsub, dtmax, dtmin, dtfac, dtave, &
                    cdtmix, tend, anstep, dt, t, tmean, tmax, etola, op_inc, &
                    anres, stomby, o2_mean, h2o_mean, co_mean
logical          :: exist, rstr            
integer :: dmi_out,dmi_in,status1       	

!default parameters for premixed CH4 combustion
data tres, tmix, tpair    / 1.e-2, 1.0e-3, 1.0e-3 / 
data cdtrp, cdtmix, dtfac / 0.1,   0.040,  1.0    /
data anres    /  1.0   /
data etola    / 1.e-5  / 
data stomby   / 1.e0   /
data const_dt /   1    /
data n0       /  5000  /
data nufm_fac /   1    /
data rstr     / .true. /
data op_inc   /  1.02  /
data savefrqc /  250  /
data pr_atm, patm / 1.0, 1.01325e6 /
data luerr, luout / 0, 6 /
data ng_pasr / 8 /

namelist / pasrnml / tres, tmix, tpair, cdtrp, cdtmix, dtfac, anres, etola, &
                     stomby, const_dt, n0, rstr, ng_pasr, nufm_fac
                     
!----------------  initialization  ---------------------------------------

!--- read namelist (if the file pasr.nml exists)

inquire( file = 'pasr.nml', exist = exist )
if ( exist ) then
   call isat_lu( lu )
   open( lu, file = 'pasr.nml', position ='rewind',action ='read')
   read( lu, nml = pasrnml )
   close( lu )
endif
	
!======  test parameters,Some parameters are read from input file  =======

!ng_pasr   = 8    ! total number of reactors (over all processors) ! Now set in namelist
n_cond    = 1    ! number of distinct conditions
mode_cond = 1    ! =1 assign conditions to reactors in blocks; =2 wrap

!--- initialize MPI

call mpi_init( ierror )
if ( ierror /= 0 ) then
   write(luerr,*) 'pasr_multi: error initializing MPI', ierror
   stop
endif
mpicomm = MPI_COMM_WORLD
call mpi_comm_size( mpicomm, nproc, ierror )
call mpi_comm_rank( mpicomm, iproc, ierror )

!INITIALIZE ALL CONNECTIONS
do j = 0,nproc-1
  if (iproc == j) then
    do i = j+1, nproc-1
      call mpi_send(dmi_out,1,MPI_INTEGER,i,123+i,mpicomm,ierror)
      call mpi_recv(dmi_in ,1,MPI_INTEGER,i,456+i,mpicomm,status1,ierror)
    enddo
    print *, iproc, 'mpi initialized...'
  elseif (iproc > j) then
    call mpi_recv(dmi_in ,1,MPI_INTEGER,j,123+iproc,mpicomm,status1,ierror)
    call mpi_send(dmi_out,1,MPI_INTEGER,j,456+iproc,mpicomm,ierror)
  endif
  call mpi_barrier(mpicomm,ierror)
enddo

write(0,*) 'Starting process: iproc, nproc = ', iproc, nproc
  
if ( mod(ng_pasr,n_cond) /= 0 ) then
   write(luerr,*)'pasr_multi: bad ng_pasr,n_cond = ', ng_pasr,n_cond
   stop
endif

if ( mod(ng_pasr,nproc) /= 0 ) then
   write(luerr,*)'pasr_multi: bad ng_pasr,nproc = ', ng_pasr,nproc
   stop
endif
print *, 'pasr', iproc, 'started...'

!--- initialize chemistry interface

	ci_info  = 0
	ci_rinfo = 0.d0
	info     = 0
	rinfo    = 0.d0

	ci_info(1)  = 1             ! constant pressure
	ci_info(2)  = const_dt      ! constant time step
	ci_info(14) = 1             ! output values of all parameters
	ci_info(16) = 0             ! logical unit number for error output

	ci_rinfo(10) = 1.d-4        ! reference time interval
	
	ci_info(8)  = 1

        ! VH, 11/28/2012: mpi_nml = 1, use single isat nml file for all the cores
        info(67) = 1                

    kase = 5
    if ( kase==4 ) then   
       rinfo(1)    = etola         ! specified absolute error tolerance
       rinfo(9)    = stomby        ! storage in megabytes allowed for ISAT table
    elseif( kase==5 ) then
       info(12)    = 2
       rinfo(1)    = etola
       rinfo(8)    = stomby
    end if
	call sciparam( ci_info, ci_rinfo, info, rinfo )
    call ciinit( ncomp, nfull, nstr )
    ns  = nfull - 4

print *, 'chemistry in', iproc, 'initialized...'

!--- misc. initializations

  press   = pr_atm * patm
  rcd     = 1
  next_op = 1

!--- allocate arrays 

  nl_pasr    = ng_pasr / nproc      ! number of reactors per processor 
  kf_gpasr = iproc * nl_pasr + 1    ! global index of first reactor on the process
  kl_gpasr = kf_gpasr + nl_pasr - 1 ! global index of last reactor on the process

  allocate( pasr(nl_pasr) )
  allocate( ftemp(ncomp) )
  allocate( fstr(ncomp,nstr) )
  allocate( g_cond(ng_pasr) )
  allocate( dfin(ncomp) )
  allocate( fin(ncomp,nstr,n_cond) )
  allocate( flstrm(nstr,n_cond) )
  
!-----  allocate an operating condition to each (global) pasr...

  if( mode_cond == 1 ) then
!  ...in blocks...
     do gpasr = 1, ng_pasr
        g_cond(gpasr) = 1 + (gpasr-1) /  (ng_pasr / n_cond)
     end do
  else
!  ...or wrapped
     do gpasr = 1, ng_pasr
        g_cond(gpasr) = 1 + mod( gpasr-1, n_cond )
     end do
  endif

!-----  get stream compositions

  do ist = 1, nstr
     call scistrm( ist, ncomp, dfin, dpt )
	 fstr(:,ist) = dfin(:)
  end do

!-----  specify operating conditions -- case dependent

  do cond = 1, n_cond
     fin(1:ncomp,1:nstr,cond)  = fstr(1:ncomp,1:nstr)
  end do
  
  if ( n_cond>1 ) then
     write(0,*) 'Coded for only 1 condition: n_cond=',n_cond
     stop
  end if
!---- initialize reactors

print *, 'streams in', iproc, 'initialized...'

  gpasr = kf_gpasr - 1
  
  allocate( tempf(ncomp) )
  if ( iproc==0 ) then
     n = nufm_fac * n0
  else
     n = n0
  end if
  !call set_nufm_fac( iproc, nufm_fac )

  if ( rstr==.false. ) then
     do lpasr = 1, nl_pasr
        gpasr = gpasr + 1
        call pasr_init( n, ncomp, nstr, gpasr, pasr(lpasr) )

        cond = g_cond(gpasr)
        pasr(lpasr)%fin(1:ncomp,1:nstr) = fin(1:ncomp,1:nstr,cond)
	    !pasr(lpasr)%flstrm(1:nstr) = flstrm(1:nstr,cond)
	    pasr(lpasr)%flstrm(1:nstr) = (/0.05, 0.85, 0.1/)  
               
        init_stream = 1
	    do i = 1, n  !  set initial condition to stream 2
	       pasr(lpasr)%f(i,1:ncomp)   = pasr(lpasr)%fin(1:ncomp,init_stream)
	    end do
     end do
  elseif ( rstr==.true. ) then
     do j = 0, nproc-1
      if (j==iproc) then
       do lpasr = 1, nl_pasr
         gpasr = gpasr + 1
         call pasr_init( n, ncomp, nstr, gpasr, pasr(lpasr) )
        
         cond = g_cond(gpasr)
         cond = 1
         pasr(lpasr)%fin(1:ncomp,1:nstr) = fin(1:ncomp,1:nstr,cond)
	    !pasr(lpasr)%flstrm(1:nstr) = flstrm(1:nstr,cond)
	    pasr(lpasr)%flstrm(1:nstr) = (/0.05, 0.85, 0.1/)  
        
         call isat_lu( lu_op )
         open( lu_op, file = pasr(lpasr)%op1, form='unformatted', status='old' )
         !do i = 1, n
         !   read(lu_op) pasr(lpasr)%f(i,1:ncomp)
         !end do
         do i = 1, n0
           read(lu_op) tempf
           do ii = 1, nufm_fac
             pasr(lpasr)%f((i-1)*nufm_fac+ii,1:ncomp) = tempf(:)
           end do
        end do        
        close( lu_op )
       end do
       print *, j, 'initialized...'
      end if
      call MPI_Barrier(mpicomm, ierror)
     end do
  end if       
        
!--- set time step parameters

    dtmax  = min( cdtrp, 0.25) * min( tres, tpair )
    dtsub  = min( cdtmix, 0.5) * tmix
	dtmax  = min( dtmax, dtsub )

	dtmin  = min( dtfac, 1.0 ) * dtmax
	dtave  = 0.5 * ( dtmax + dtmin )

    tend   = anres * tres
    anstep = tend / dtave

!--- write values of parameters to standard output (processor 0 only)
  
    if( iproc == 0  .and.  luout >= 0 ) then
	    write(luout,1 ) nproc
	    write(luout,2 ) iproc
	    write(luout,3 ) ng_pasr
	    write(luout,4 ) nl_pasr
	    write(luout,5 ) n_cond
	    write(luout,6 ) g_cond
        write(luout,7 ) ns
        write(luout,8 ) nstr
        write(luout,10) n
        write(luout,12) tres
        write(luout,14) tmix
        write(luout,16) tpair
        write(luout,17) dtmin, dtmax
        write(luout,18) anstep
        write(luout,19) const_dt
        write(luout,20) etola
        write(luout,21) stomby

1       format(' Number of processes               =', i8 )
2       format(' Rank of this process              =', i8 )
3       format(' Total number of reactors          =', i8 )
4       format(' Number of reactors this process   =', i8 )
5       format(' Number of distinct conditions     =', i8 )
6       format(' Index of conditions of all react. =', (20i4) )
7       format(' Number of species                 =', i8 )
8       format(' Number of streams                 =', i8 )
10      format(' Number of particles               =', i8 )
12      format(' Residence time, tres (seconds)    =', 1p,e13.4)
14      format(' Mixing time, tmix (seconds)       =', 1p,e13.4)
16      format(' Pairing time, tpair (seconds)     =', 1p,e13.4)
17      format(' Min. and max. time step (seconds) =', 1p,2e13.4)
18      format(' Number of steps                   =', 1p,e13.4)
19      format(' 1:Constant time step;0:varying step=',i8)
20      format(' error tolerance set for ISAT      =', 1p,e13.4)
21      format(' Storage set for ISAT (Mbyte)      =', 1p,e13.4)
    endif

!===============  start of time steps  =================================

    istep  = 0
	dt     = dtave
    t      = -dt

write(0,*) pasr(1)%t, tend

time_steps: do while (pasr(1)%t < tend .and. istep<floor(anstep))  
  
  istep  = mod( istep + 1, 1000000 )
  
  if (iproc.eq.0) write(luout,*) 'istep', istep

  reactor_loop: do lpasr = 1, nl_pasr  !  loop over all reactors 

  call rnuput( pasr(lpasr)%ip1, pasr(lpasr)%ip2 ) !re-set random number sequence
 
!--- set the time step

	if( mod(istep,2) == 0 ) then
	   dt = 2.*dtave - pasr(lpasr)%dt
	else
	   if( dtmin == dtmax ) then
	      dt = dtmin
	   else
	      call rnu( dt )
	      dt     = dtmin + dt * ( dtmax - dtmin )
	   endif
	endif
	
	pasr(lpasr)%dt = dt
    pasr(lpasr)%t  = pasr(lpasr)%t + dt

!-------  outflow/ inflow -------------------------------------------

    pasr(lpasr)%prout  = pasr(lpasr)%prout + 0.5 * n * dt / tres
    nprout             = pasr(lpasr)%prout
    pasr(lpasr)%prout  = pasr(lpasr)%prout - nprout

    if( nprout > 0 ) then

!  select 2 * nprout pairs of particles at random; put at end of f array

       npick  = 2 * nprout
       call pickpr( npick, n, ncomp, n, pasr(lpasr)%f )

!  set alternate particles to inflowing properties

       i = n + 2
       do ii = 1, nprout
          i     = i - 2
          istrm = inflow( pasr(lpasr) )
	      pasr(lpasr)%f(i,:) = pasr(lpasr)%fin(:,istrm)
	   end do

	endif

!-------  pairing  ------------------------------------------------

    pasr(lpasr)%prpar  = pasr(lpasr)%prpar + 0.5 * n * dt / tpair
    nprpar = pasr(lpasr)%prpar

    if( nprpar >=2 ) then

       pasr(lpasr)%prpar  = pasr(lpasr)%prpar - nprpar

!  select nprpar pairs of particles at random; put at end of f array

       call pickpr( nprpar, n, ncomp, n, pasr(lpasr)%f )

!  rotate particle properties

       ifst = n - 2 * nprpar + 2
       ftemp(1:ncomp)    = pasr(lpasr)%f(n,1:ncomp)

       do i = n, ifst+2, -2
	      pasr(lpasr)%f(i,1:ncomp)   = pasr(lpasr)%f(i-2,1:ncomp)
	   end do

	   pasr(lpasr)%f(ifst,1:ncomp)= ftemp(1:ncomp)

    endif

!---- mixing step  --------------------------------------------------

    call mix( dt, tmix, n, ncomp, n, pasr(lpasr)%f )

	call rnuget( pasr(lpasr)%ip1, pasr(lpasr)%ip2 )  ! get random number sequence

  end do reactor_loop  

!================= reaction  ================================================= 
   
  call rxn( mpicomm, press, ncomp, n, nl_pasr, pasr, stomby )
  !call rxn_bk( mpicomm, press, ncomp, n, nl_pasr, pasr )
  
!-----  output
  if ( istep == next_op ) then
     next_op = max( istep+1.d0, istep * op_inc )
     do lpasr = 1, nl_pasr
        tmean = sum( pasr(lpasr)%temp(1:n) ) / float(n)
        tmax  = maxval( pasr(lpasr)%temp(1:n) )
        o2_mean  = sum( pasr(lpasr)%f(1:n,1) ) / float(n)
        h2o_mean = sum( pasr(lpasr)%f(1:n,3) ) / float(n)
        co_mean  = sum( pasr(lpasr)%f(1:n,6) ) / float(n)
	    cond  = g_cond(pasr(lpasr)%gpasr)
	    call isat_lu( lu_op )
        open( lu_op, file = pasr(lpasr)%op, position = "append" )
	    write(lu_op,'(1p,1i8,10e14.6,6i4)') &
		    istep, pasr(lpasr)%t, tmean, pasr(lpasr)%temp(1), tmax, o2_mean, h2o_mean, co_mean, &
            pasr(lpasr)%f(1,1), pasr(lpasr)%f(1,3), pasr(lpasr)%f(1,6), &
            nproc, iproc, ng_pasr, nl_pasr, lpasr, cond
	    close( lu_op )
     end do
  endif
  
  if ( istep==rcd*savefrqc ) then
     rcd = rcd + 1
     do lpasr = 1, nl_pasr
        gpasr = gpasr + 1        
        call isat_lu( lu_op )
        open( lu_op, file = pasr(lpasr)%op1, form='unformatted', position='rewind' )
        do i = 1, n
           write(lu_op) pasr(lpasr)%f(i,1:ncomp)
        end do        
        close( lu_op )
     end do
  end if
   
  !write(luout,*) 'end of istep', istep
  !if (istep==50) stop
end do time_steps  

write(0,*)'pasr_multi: stopped, iproc, t, istep = ', iproc, pasr(1)%t, istep
call mpi_finalize( ierror )
stop

end program pasr_multi
