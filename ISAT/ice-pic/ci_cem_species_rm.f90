!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ci_cem_species_rm
  ! Species selection based on reaction mapping error

  use ci_dat8
  use ci_stats
  use ci_utils
  use ci_cksubs
  use ci_ice_cksubs
  use ci_rmap
  use ci_cem_recon
  implicit none

  ! local variables
  integer, parameter :: nr_max = 2500
  integer            :: lu_in, lu_ce, lu_pl, lu_ck, lu_ct, lu_er, lu_ti
  integer            :: nx_FC, nf_FC, info
  real(kind(1.d0))   :: dt, press
  real(kind(1.d0))   :: cpu_start, cpu_stop, cpu_tic, cpu_tac

  ! full composition (FC)
  logical :: kerr
  character(16) :: symb(ns)
  real(kind(1.d0)) :: FC(nr_max, 2*nsp1), zr_ref
  real(kind(1.d0)) :: hin, hsin, z_FC(ns), zr_FC(ns), x_CE(nsp1)
  real(kind(1.d0)), allocatable :: x_FC(:), f_FC(:)
  real(kind(1.d0)) :: zr(nsp1,nsp4), zr_CE(ns)

  ! RCCE variables
  real(kind(1.d0)) :: z_CE(ns), T_CE, stats(20)

  ! reduced composition (RC)
  integer :: i, j, k, s, m, num_rs, minl(1), rep(ns), unrep(ns), det(ns)
  integer :: nsd, nsu, sp_order(ns), selected(ns)
  real(kind(1.d0)) :: err_z(ns), err2_z(ns)
  real(kind(1.d0)), allocatable :: rin(:), BT(:,:)
  
  write(0,*) 'Computing species order based on reaction mapping error.'

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
  call isat_lu( lu_pl )
  call isat_lu( lu_ck )
  call isat_lu( lu_ct )
  call isat_lu( lu_er )
  call isat_lu( lu_ti )

  open( lu_in, file='fc_xe_fe.op' )
  open( lu_ck, file='chemkin_sp_order.op' )
  open( lu_ce, file='cem_sp_order_rm_2500.op' )
  open( lu_er, file='cem_sp_error_rm_2500.op' )
  open( lu_ti, file='cem_sp_times_rm_2500.op' )
  open( lu_pl, file='cem_uusp_errors_rm_2500.op' )
  open( lu_ct, file='cem_sp_cputime_rm_2500.op')

  ! Read FC data in FC array
  do m = 1, nr_max
     read(lu_in,*, end=100, err=100) FC(m,:)
  enddo

  ! compute zr_ref
  zr_ref = 0.d0
  do m = 1, nr_max
     zr_ref = zr_ref + dot_product(FC(m,nsp1+1:nsp1+ns),FC(m,nsp1+1:nsp1+ns))
  enddo
  zr_ref = sqrt(zr_ref/nr_max)

  ! Write chemkin species order to file
  do i = 1,ns
     write(lu_ck,*) symb(i)
  enddo

  ! cpu time
  call cpu_time( cpu_start )

  ! selected: species index 1 if selected, 0 if not selected
  selected = 0

  ! Represented species
  ! You may fix some of the major species here if needed like
  ! num_rs = 4
  ! rep(1:num_rs) = (/7,3,9,5/) !CH4, O2, CO2 and H20 from chem.out
  num_rs = 0
  rep = 0

  ! Mark selected species if any
  do i = 1, num_rs
     selected(rep(i)) = 1
  enddo

  ! Find the undetermined unrepresented species
  call sp_ordering( num_rs, num_rs + ne, rep(1:num_rs), nsd, nsu, sp_order )
  unrep(1:nsu) = sp_order(nsd+1:ns)

  ! Print fixed represented species if any
  if ( num_rs > 0 ) then
     do i = 1, num_rs
        write(lu_ce,*) symb(rep(i))
        write(0,*) symb(rep(i))
     enddo
  endif

  species: do
     call cpu_time( cpu_tic )

     num_rs = num_rs + 1
     err_z = 0.d0

     allocate( rin(num_rs+ne) )
     allocate( BT(num_rs+ne,ns) )

     do s = 1, nsu
        rep(num_rs) = unrep(s)
        call BT_red( num_rs, rep(1:num_rs), BT )

        do m = 1, nr_max
           ! Initialize z_FC and zr_FC
           z_FC = FC(m, 1:ns)
           zr_FC = FC(m, nsp1+1:nsp1+ns)
           hsin = FC(m, nsp1)
           call hs2h( hsin, href, z_FC, ns, hin )

           rin = matmul(BT, z_FC)
           call routine_start(i_ceq_nrs)
           call ceq_nrs( num_rs, rep(1:num_rs), rin, hin, press, z_CE, T_CE, stats, info )
           call routine_stop(i_ceq_nrs)
           x_CE(1:ns) = z_CE(1:ns)
           x_CE(nsp1) = T_CE

           ! Compute Reaction Mapping
           call routine_start(i0_rmap2)
           call rmap2( 0, nsp1, nsp3, x_CE, press, dt, zr )
           call routine_stop(i0_rmap2)

           zr_CE(1:ns) = zr(1:ns,1)

           ! Compute Errors
           err_z(s) = err_z(s) + dot_product(zr_CE - zr_FC, zr_CE - zr_FC)
        enddo
     enddo

     call cpu_time( cpu_tac )
     cpu_tac = cpu_tac - cpu_tic
     write( lu_ct,'(2i8,1p,3e13.4)' ) num_rs, nsu, 1.e6*cpu_tac, 1.e6*cpu_tac/(nr_max*nsu)
     
     write(lu_pl,'(i3)') nsu
     do i = 1, nsu
        ! Compute 2-norm errors
        err2_z(i)  = sqrt(err_z(i)/nr_max)/zr_ref
        write(lu_pl,'(i3,1p,10e26.17)') unrep(i), err2_z(i)
     enddo

     ! minimum error location (based on err_z)
     minl = minloc(err2_z(1:nsu))
     rep(num_rs) = unrep(minl(1))
     selected(rep(num_rs)) = 1
     
     write(lu_ce,*) symb(rep(num_rs))
     write(0,*) symb(rep(num_rs))
     write(lu_er,'(i3,1p,10e26.17)') num_rs, err2_z(minl(1))

     ! Set the selected species index to zero
     unrep(minl(1)) = 0
     
     ! find the undetermined unrepresented species
     call sp_ordering( num_rs, num_rs + ne, rep(1:num_rs), nsd, nsu, sp_order )

     ! Exit if no undetermined unrepresented species left
     if( nsu == 0 ) then
        ! Write the remaining undetermined species
        do i = 1, ns - num_rs
           if(unrep(i) /= 0) then
              write(lu_ce,*) symb(unrep(i))
              write(0,*) symb(unrep(i))
              num_rs = num_rs + 1
           endif
        enddo
        ! Write the determined species
        do i = 1, ns - num_rs
           if(det(i) /= 0) then
              write(lu_ce,*) symb(det(i))
              write(0,*) symb(det(i))
              num_rs = num_rs + 1
           endif
        enddo

        exit species
     endif

     ! Note the remianing undetermined and determined species index
     unrep = 0
     unrep(1:nsu) = sp_order(nsd+1:ns)

     ! Identify species determined using atom conservation
     det = 0
     k = 1
     do i = 1,nsd
        if(selected(sp_order(i)) == 0) then
           det(k) = sp_order(i)
           k = k + 1
        endif
     end do

     if(allocated(rin))  deallocate(rin)
     if(allocated(BT))   deallocate(BT)

     call cpu_time( cpu_stop )
     write( lu_ti,'(i8,1p,3e13.4)' ) num_rs, (cpu_stop - cpu_start)

  enddo species

  ! compute cputime
  call cpu_time( cpu_stop )
  cpu_stop = cpu_stop - cpu_start
  write( lu_ct,'(2i8,1p,3e13.4)' ) ns, 0, 1.e6*cpu_stop, 1.e6*cpu_stop/(nr_max*ns)
 
100 continue
  return

end subroutine ci_cem_species_rm
