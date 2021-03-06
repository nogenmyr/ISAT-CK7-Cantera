!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ci_ice_species_rm
  ! Species selection based on reaction mapping error

  use ci_8
  use ci_dat8
  use ci_utils
  use ci_cksubs
  use ci_ice_cksubs
  use ci_rmap
  use ci_cem_recon
  implicit none

  ! local variables
  integer, parameter :: nr_max = 50
  integer            :: lu_in, lu_ce, lu_pl, lu_ck
  integer            :: nx_FC, nf_FC
  real(kind(1.d0))   :: dt, press

  ! full composition (FC)
  logical :: kerr
  character(16) :: symb(ns)
  real(kind(1.d0)) :: FC(nr_max, 2*nsp1), zr_ref
  real(kind(1.d0)) :: hin, hsin, z_FC(ns), zr_FC(ns), x_ICE(nsp1)
  real(kind(1.d0)), allocatable :: x_FC(:), f_FC(:)
  real(kind(1.d0)) :: zr(nsp1,nsp4), zr_ICE(ns)

  ! ICE-PIC variables
  real(kind(1.d0)) :: z_CE(ns), z_ICE(ns), T_ICE, tau_ICE, z_g(ns), &
       T_g, A_ICE_ODE(nsp1,nsp1)
  integer :: k_facet, n_iters, n_integs, method, info_ice

  ! reduced composition (RC)
  integer :: i, k, s, m, num_rs, minl(1), rep(ns), unrep(ns), det(ns)
  integer :: nsd, nsu, sp_order(ns), selected(ns)
  real(kind(1.d0)) :: err_z(ns), err2_z(ns)
  real(kind(1.d0)), allocatable :: rin(:), BT(:,:)
  real(kind(1.d0)), allocatable :: r_g(:), r_ICE(:), r_normal(:)
  
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

  open( lu_in, file='fc_xe_fe.op' )
  open( lu_ck, file='chemkin_sp_order.op' )
  open( lu_ce, file='ice_sp_order_rm_50.op' )
  open( lu_pl, file='ice_uusp_errors_rm_50.op' )

  ! Read FC data in FC array
  do m = 1, nr_max
     read(lu_in,*, end=100, err=100) FC(m,:)
  enddo

  ! compute z_ref
  zr_ref = 0.d0
  do m = 1, nr_max
     zr_ref = zr_ref + dot_product(FC(m,nsp1+1:nsp1+ns), FC(m,nsp1+1:nsp1+ns))
  enddo
  zr_ref = sqrt(zr_ref/nr_max)

  ! Write chemkin species order to file
  do i = 1,ns
     write(lu_ck,*) symb(i)
  enddo

  ! selected: species index 1 if selected, 0 if not selected
  selected = 0
  rep = 0

  ! Represented species
  ! You may fix some of the major species here if needed like
  num_rs = 1
  rep(1:num_rs) = (/7/) !check chem.out
  !num_rs = 0
  !rep = 0

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

     num_rs = num_rs + 1
     err_z = 0.d0

     allocate( rin(num_rs+ne) )
     allocate( r_g(num_rs+ne) )
     allocate( r_ICE(num_rs+ne) )
     allocate( r_normal(num_rs+ne) )
     allocate( BT(num_rs+ne,ns) )

     do s = 1, nsu
        rep(num_rs) = unrep(s)
        call BT_red( num_rs, rep(1:num_rs), BT )
        call ci_reset9( num_rs, rep(1:num_rs), BT )

        do m = 1, nr_max
           ! Initialize z_FC = x_FC(1:ns)
           z_FC = FC(m, 1:ns)
           zr_FC = FC(m, nsp1+1:nsp1+ns)
           hsin = FC(m, nsp1)
           call hs2h( hsin, href, z_FC, ns, hin )

           rin = matmul(BT, z_FC)
           call ci_ice_recon( rin, hin, press, z_CE, r_g, z_g, T_g,  &
                tau_ICE, r_ICE, z_ICE, T_ICE, A_ICE_ODE, r_normal, &
                k_facet, n_iters, n_integs, method, info_ice )

           write(0,'(a,5i4)') 'nrs, s, nr, mth, info =', num_rs, s, m, &
                method, info_ice
           
           x_ICE(1:ns) = z_ICE(1:ns)
           x_ICE(nsp1) = T_ICE
           
           ! Compute Reaction Mapping
           call rmap2( 0, nsp1, nsp3, x_ICE, press, dt, zr )
           zr_ICE(1:ns) = zr(1:ns,1)

           ! Compute Errors
           err_z(s) = err_z(s) + dot_product(zr_ICE - zr_FC, zr_ICE - zr_FC)
        enddo
     enddo
     
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
     
     write(lu_ce,*) symb(unrep(minl(1)))
     write(0,*) symb(unrep(minl(1)))
     
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
     if(allocated(r_g))  deallocate(r_g)
     if(allocated(r_ICE))  deallocate(r_ICE)
     if(allocated(r_normal))  deallocate(r_normal)
     if(allocated(BT))   deallocate(BT)

  enddo species
 
100 continue
  return

end subroutine ci_ice_species_rm
