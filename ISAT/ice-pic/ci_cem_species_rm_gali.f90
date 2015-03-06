!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ci_cem_species_rm_gali
  ! Greedy With Local Imporvements (GALI). Run the Greedy Algorithm
  ! 'ci_cem_species_rm' routine before running this

  use ci_dat8
  use ci_stats
  use ci_utils
  use ci_cksubs
  use ci_ice_cksubs
  use ci_rmap
  use ci_cem_recon
  implicit none

  ! local variables
  integer, parameter :: nr_max = 5000
  integer, parameter :: nr_in  = 500
  integer, parameter :: f_err  = 500
  integer, parameter :: nrs_max = 25
  integer            :: lu_in, lu_ce, lu_ge, lu_sp, lu_er, lu_it, lu_ti
  integer            :: nx_FC, nf_FC, info, spidx, nswaps, niter
  real(kind(1.d0))   :: dt, press
  real(kind(1.d0))   :: cpu_start, cpu_stop, cpu_tic, cpu_tac

  ! full composition (FC)
  logical :: kerr, changed
  character(16) :: symb(ns), repsp
  real(kind(1.d0)) :: FC(nr_max, 2*nsp1), zr_ref
  real(kind(1.d0)) :: hin, hsin, z_FC(ns), zr_FC(ns), x_CE(nsp1)
  real(kind(1.d0)), allocatable :: x_FC(:), f_FC(:)
  real(kind(1.d0)) :: zr(nsp1,nsp4), zr_CE(ns)

  ! RCCE variables
  real(kind(1.d0)) :: z_CE(ns), T_CE, stats(20)

  ! reduced composition (RC)
  integer :: i, j, k, s, m, num_rs, minl(1), rep(ns), unrep(ns), det(ns)
  integer :: nsd, nsu, sp_order(ns), srep(ns), crep(ns), trep(ns), tmp
  real(kind(1.d0)) :: gerr(ns), err_z(ns), min_err, new_min_err, error, ngerr

  write(0,*) 'Species selection with GALI'

  ! Get species names
  call ctsyms(  6, symb, kerr, gas )

  ! Read x = z(0), hs(0), f = z(dt), hs(dt)
  nx_FC = nsp1
  nf_FC = nsp1

  ! XXX Assuming constant pressure and dt mode XXX
  dt = dtc
  press = prc

  ! allocate arrays for FC
  allocate( x_FC(nx_FC) )
  allocate( f_FC(nf_FC) )

  call isat_lu( lu_in )
  call isat_lu( lu_ce )
  call isat_lu( lu_sp )
  call isat_lu( lu_er )
  call isat_lu( lu_ge )
  call isat_lu( lu_it )
  call isat_lu( lu_ti )

  open( lu_in, file='fc_xe_fe.op' )
  open( lu_ce, file='cem_sp_order_rm_5000.op' )
  open( lu_ge, file='cem_sp_error_rm_5000.op' )
  open( lu_sp, file='cem_sp_order_rm_gali_500.op' )
  open( lu_er, file='cem_sp_error_rm_gali_500.op' )
  open( lu_ti, file='cem_sp_times_rm_gali_500.op' )
  open( lu_it, file='cem_sp_error_rm_gali_iters_500.op' )

  ! Read FC data in FC array
  do nr = 1, nr_max
     read(lu_in,*, end=100, err=100) FC(nr,:)
  enddo

  ! compute zr_ref
  zr_ref = 0.d0
  do m = 1, nr_max
     zr_ref = zr_ref + dot_product(FC(m,nsp1+1:nsp1+ns),FC(m,nsp1+1:nsp1+ns))
  enddo
  zr_ref = sqrt(zr_ref/nr_max)

  ! cpu time
  call cpu_time( cpu_start )

  num_rs = 0
  srep = 0

  ! Represented species
  do i = 1, min(ns, nrs_max)
     read(lu_ce,'(A)', end=100, err=100) repsp
     call ctcomp(repsp, spidx, gas)
     srep(i) = spidx
  enddo

  do i = 1, min(ns - ne - 1, nrs_max)
     ! read greedy spo errors
     read(lu_ge,*) tmp, gerr(i)
     write(0, '(i3,1p,10e26.17)') i, gerr(i)

     ! compute the error
     ! call compute_error( nr_max, i, srep(1:i), error )
     ! gerr(i) = sqrt(error/nr_max)/zr_ref
     ! write(lu_ge, '(i3,1p,10e26.17)') i, gerr(i)
  enddo

  do num_rs = 1, min(ns - ne - 1, nrs_max)
     
     if(num_rs < 3) then
        write(lu_sp,'(100a16)') symb(srep(1:num_rs))
        write(lu_er,'(i3,1p,10e26.17)') num_rs, gerr(num_rs)
        cycle
     endif

     rep = srep
     crep = srep
     trep = srep

     ! print inital set of species
     write(0,'(i3,a2,100a16)') num_rs, '', symb(rep(1:num_rs))

     ! initial error
     call compute_error( nr_in, num_rs, rep(1:num_rs), min_err )

     write(f_err,'(1p,e16.8,a2,100a16)') min_err,'',symb(rep(1:num_rs)),'**'

     ! Find the undetermined unrepresented species
     call sp_ordering( num_rs, num_rs + ne, rep(1:num_rs), nsd, nsu, sp_order )
     unrep(1:nsu) = sp_order(nsd+1:ns)

     k = 1 
     changed = .false.
     niter = 0

     write(lu_it,'(i3)') num_rs
     write(lu_it,'(2i3,1p,e26.17)') 0, 0, gerr(num_rs)

     species: do
        nswaps = 0
        niter = niter + 1
        do i = 2, num_rs
           ! leave the last species in the first run if unchanged, already checked
           if( k == 1 .and. .not.changed .and. i == num_rs ) cycle

           err_z = 0.d0
           do s = 1, nsu
              trep(i) = unrep(s)
              call compute_error( nr_in, num_rs, trep(1:num_rs), err_z(s) )
              write(f_err,'(1p,e16.8,a2,100a16)') err_z(s),'',symb(trep(1:num_rs))
           enddo

           ! check if the error improved
           new_min_err = minval(err_z(1:nsu))
           if(new_min_err < min_err) then
              nswaps = nswaps + 1
              changed = .true.
              minl = minloc(err_z(1:nsu))
              rep(i) = unrep(minl(1))
              min_err = new_min_err

              write(f_err,'(1p,e16.8,a2,100a16)') min_err,'',symb(rep(1:num_rs)),'<--<<'

              ! find the undetermined unrepresented species
              call sp_ordering( num_rs, num_rs + ne, rep(1:num_rs), nsd, nsu, sp_order )
              unrep(1:nsu) = sp_order(nsd+1:ns)
           endif

           trep = rep
        enddo

        call compute_error( nr_max, num_rs, rep(1:num_rs), error )
        ngerr = sqrt(error/nr_max)/zr_ref

        ! Exit if the error increased
        if(ngerr < gerr(num_rs)) then
           write(lu_it,'(2i3,1p,e26.17)') niter, nswaps, ngerr
        else
           write(lu_it,'(2i3,1p,e26.17)') niter, 0, gerr(num_rs)
           exit species
        endif

        ! Exit if the species ordering didn't change
        if( maxval(abs(rep(1:num_rs) - crep(1:num_rs))) == 0 ) then
           exit species
        endif

        ! reset crep to current rep
        crep = rep
        k = 0
     enddo species

     ! accept new species ordering if the error improved
     if(ngerr < gerr(num_rs)) then
        write(0,'(i3,a2,100a16)') num_rs, '', symb(rep(1:num_rs))
        write(lu_sp,'(100a16)') symb(rep(1:num_rs))
        write(lu_er,'(i3,1p,10e26.17)') num_rs, sqrt(error/nr_max)/zr_ref
     else
        write(0,'(i3,a16,100a16)') num_rs, 'no improvement'
        write(lu_sp,'(100a16)') symb(srep(1:num_rs))
        write(lu_er,'(i3,1p,10e26.17)') num_rs, gerr(num_rs)
     endif

     ! compute cputime
     call cpu_time( cpu_stop )
     write( lu_ti,'(i8,1p,3e13.4)' ) num_rs, (cpu_stop - cpu_start)

  enddo

100 continue
  return

contains

  subroutine compute_error(nr_pts, nrsin, repin, n_err)

    integer, intent(in) :: nr_pts, nrsin, repin(1:nrsin)
    real(kind(1.d0)), intent(out) :: n_err
    real(kind(1.d0)), allocatable :: rin(:), BT(:,:)

    allocate( rin(nrsin+ne) )
    allocate( BT(nrsin+ne,ns) )
    call BT_red( nrsin, repin(1:nrsin), BT )

    n_err = 0.d0
    do nr = 1, nr_pts
       ! Initialize z_FC = x_FC(1:ns)
       z_FC = FC(nr, 1:ns)
       zr_FC = FC(nr, nsp1+1:nsp1+ns)
       hsin = FC(nr, nsp1)
       call hs2h( hsin, href, z_FC, ns, hin )

       rin = matmul(BT, z_FC)
       call routine_start(i_ceq_nrs)
       call ceq_nrs( nrsin, repin(1:nrsin), rin, hin, press, z_CE, T_CE, stats, info )
       call routine_stop(i_ceq_nrs)

       x_CE(1:ns) = z_CE(1:ns)
       x_CE(nsp1) = T_CE

       ! Compute Reaction Mapping
       call routine_start(i0_rmap2)
       call rmap2( 0, nsp1, nsp3, x_CE, press, dt, zr )
       call routine_stop(i0_rmap2)

       zr_CE(1:ns) = zr(1:ns,1)

       ! Compute Errors
       n_err = n_err + dot_product(zr_CE - zr_FC, zr_CE - zr_FC)
    enddo
    if(allocated(BT))   deallocate(BT)
    if(allocated(rin))  deallocate(rin)

  end subroutine compute_error

end subroutine ci_cem_species_rm_gali
