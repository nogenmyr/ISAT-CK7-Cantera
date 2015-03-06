!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

      subroutine cirmap_dr_ice_pic( need, nx, x, nf, nh, iusr, rusr,
     1     f, dfdx, hvar )

!     Determine the reaction mapping (f), its gradient (dfdx), and additional 
!     variables (hvar) as a function of the ISAT tabulation parameters (x). 
!     This version performs dimension reduction using ICE-PIC.

!     For ICE-PIC (info(47)=de_nearby=1), f, g and h are evaluated not exactly at x, 
!     but at a nearby location, xe.  On entry, rusr(2) contains the tolerance  tol_re 
!     which is used to enforce  |x(1:nrc)-xe(1:nrc)| < tol_re.

!  On return from cirmap_dr_ice_pic:
!     x          contains xe 
!     rusr(1:nf) contains f(xe)
!     f          contains an approximation to f(x). 
!     
!  definitions:
!     x = r(0), h_sn(0), p, dt  for mode_pdt = 1, general case,   nx=nc+2
!     x = r(0), h_sn(0), p      for mode_pdt = 2, const dt,       nx=nc+1
!     x = r(0), h_sn(0), dt     for mode_pdt = 3, const p,        nx=nc+1
!     x = r(0), h_sn(0)         for mode_pdt = 4, const p and dt, nx=nc
!     r(0): represented species + elements in unrepresented species
!     nc= nrs+ne+1
!     f  = r(dt), h_sn(dt), T(dt);   nf = nc+1

!  input:
!     need(1)	= 1 - f is needed (otherwise need(1) = 0 )
!     need(2)	= 1 - df/dx is needed
!     need(3)	= 1 - hvar is needed
!     nx	- number of parameters
!     x	        - parameters 
!     nf	- number of dependent variables (nf=nsp2)
!     nh	- number of additional variables
!     rusr(1)   - length of rusr, must be at least 2
!     rusr(2)   = tol_re

!  output:
!     x     -  xe: the point where reconstruction is performed
!     f	    - variables; an approximation to f(x)
!     dfdx  - df/dx
!     hvar  - additional variable (hvar(1) = T_DR(dt))
!             for diagnostic purposes, if nh >= 3, 
!                  hvar(2) = info_ice
!                  hvar(3) = k_facet
!     rusr  - rusr(1:nf) contains f(xe)

      use ci_dat
      use ci_dat8
      use ci_rmap
      use ci_stats
      use ci_utils
      use ci_ice_cksubs
      use ci_dpt_dr

      implicit none
      
      integer, intent(in)     :: need(3), nx, nf, nh, iusr(1)
      real(k_dp), intent(inout)  :: x(nx), rusr(*)
      real(k_dp), intent(out) :: f(nf), dfdx(nf,nx), hvar(max(nh,1))
      real(k_dp)              :: press, dt
      real(k_dp), save        :: calls=0.d0
      
!     local arrays

      real(k_dp), allocatable :: ROE(:,:), SR(:,:), RM(:,:), RED(:,:), 
     1     APP(:,:)
      real(k_dp), allocatable :: A(:,:), A_ICE(:,:)
      
      real(k_dp) :: hs_e_n, rin(nrs+ne), hin, hsin_n, hs_n, y(ns),
     1              cp1, enth1(ns), tol_re, tmp(ns), fnew(nf), zu(nus)

      real(k_dp) :: hi(ns), hri(nrs), cpi(ns), hua(ne), cpa, dpt(3)
           
!     saved quantities from  ci_ice_recon  
      real(k_dp),  save :: h_e, T_CE, T_DR, T_a, T_e_ICE, tau_ICE, T_g, 
     1      hsin_nsave, press_save, dt_save
      
      real(k_dp),  save, allocatable ::  xe(:), r_e(:), r_g(:), z_CE(:),
     1     z_e_ICE(:), z_g(:),  B_red_full(:, :), A_ICE_ODE(:,:),
     2     zin(:), zr(:,:), rin_save(:), A_CEM(:,:) , r_normal(:)
      
      integer, save :: k_facet, ns_pos, nx_save, lu_dice, 
     1     info_ice, method
      integer, save, allocatable :: index_pos(:)
      
      integer :: ifst=0, calls_diag = -1, n_iters, n_integs
      integer :: check = 1, info, i
      logical :: match, call_ICE_ODE
      integer, parameter :: ice_diag = -1 ! >0 for ICE diagnostics      

!----------------------start of execution  -------------------------------
      if( ifst == 0 ) then
         ifst = 1
         allocate( xe(nx) ) 
         allocate( r_e(nrs+ne) )
         allocate( r_g(nrs+ne) )
         allocate( r_normal(nrs+ne) )
         allocate( z_CE(ns)  )
         allocate( z_e_ICE(ns) )
         allocate( z_g(ns) )
         allocate( B_red_full(ns, nrc-1) )
         allocate( A_ICE_ODE(ns+1, ns+1) )
         allocate(index_pos(ns)) 
         allocate(zin(nsp1))
         allocate(zr(nsp1,nsp4)) 
         allocate(rin_save(nrs+ne))
         allocate(A_CEM(ns+3,nc+2)) 
         
         if( ice_diag > 0 ) then
            call isat_lu( lu_dice )
            open( lu_dice, file='dice.op')
         endif
         
      endif   
      
      calls = calls + 1.d0

      if( need(1) == 0  .and.  need(2) == 0 ) return
      
!------check input  and  determine  press  and  dt-----------------------
      
      if( nf /= nc+1 ) call isat_abort('cirmap_dr_ice_pic', 1,
     1     mess = 'bad value of nf', ivar = (/ nf, nc+1 /) )
      if( nh < 0 ) call isat_abort('cirmap_dr_ice_pic', 2,
     1     mess = 'bad value of nh', isv = nh )
      if( nint( rusr(1) ) < 2 ) call isat_abort('cirmap_dr_ice_pic', 3,
     1     mess = 'bad value of rusr(1)', rsv = rusr(1) )
           
      if( mode_pdt == 1 ) then  !  variable  p  and  dt
         if( nx /= nc+2 ) call bad_nxDR_ICE( nx, nc+2 )
         press = x(nc+1)
         dt    = x(nc+2)
      elseif( mode_pdt == 2 ) then !  variable  p, const  dt
         if( nx /= nc+1 ) call bad_nxDR_ICE( nx, nc+1 )
         press = x(nc+1)
         dt    = dtc
      elseif( mode_pdt == 3 ) then !  const  p,  variable  dt
         if( nx /= nc+1 ) call bad_nxDR_ICE( nx, nc+1 )
         press = prc
         dt    = x(nc+1)
      elseif( mode_pdt == 4 ) then !  const  p,  and  dt
         if( nx /= nc ) call bad_nxDR_ICE( nx, nc )
         press = prc
         dt    = dtc
      else
         call isat_abort('cirmap_dr_ice_pic', 4, 
     1        mess = 'bad value of mode_pdt', isv = mode_pdt )
      endif

!--------  enforce realizability  ----------------------------------
      do i = 1, nrs+ne
         rin(i)  = max( x(i), 0.d0 )
      end do

! enforce normalization
      rin = rin/dot_product(rin, amolwt_n)
      
      hsin_n = x(nrs+ne+1)
      call hs2h_n( hsin_n, href_n, rin, nrc, hin)
          
!----------------determine whether or not to call the reconstructed and DDASAC subroutine---------------------------
      if( need(2) == 0 ) then   ! gradient not needed
         call_ICE_ODE = .true.
!     update the saved quantities 
         rin_save   = rin
         hsin_nsave = hsin_n
         nx_save    = nx
         dt_save    = dt
         press_save = press   
         
      else                      !  gradient needed 
         match = .true.
         if( nx /= nx_save  .or.  press /= press_save  
     1        .or.  dt /= dt_save ) match= .false.
         do i = 1, nrs+ne
            if( abs(rin(i)- rin_save(i)) >1.d-12) match =.false. 
         end do
         if( abs(hsin_n- hsin_nsave) >1.d-5)  match =.false. 
         
         
         if( match ) then    
            call_ICE_ODE = .false.
         else
            write(lu_err,*)'cirmap_dr_ice_pic: match failed'
            call_ICE_ODE = .true.
         endif

      endif
!----------------------------------------------------------------------------------
      
      
!--------determine constrained species mass fraction and temperature-------------
      if (call_ICE_ODE) then 
         tol_re = rusr(2)
        
         call routine_start(i_ci_ice_recon)

         call ci_ice_recon( rin, hin, press, z_CE, r_g, z_g, T_g, 
     1        tau_ICE, r_e, z_e_ICE, T_e_ICE, A_ICE_ODE, r_normal, 
     2        k_facet, n_iters, n_integs, method, info_ice )

         call routine_stop(i_ci_ice_recon)

         call temphz( hin, z_CE, T_CE, check, info )
         h_e = hin

         ! if ci_ice_recon failed z_CE is the best approximation
         if( info_ice < 0 ) then
            z_e_ICE = z_CE
            r_e = matmul( BBT, z_CE )
            T_e_ICE = T_CE
            if( ice_diag > 0 ) then
               write(lu_dice,*) nint(calls), 
     1              'ci_ice_recon: all methods failed'
               call print_var(rin, "rin", lu_dice)
               call print_var(z_CE, "z_CE", lu_dice)
               call print_var(z_g, "zg", lu_dice)
               call print_var(z_e_ICE, "z_ICE", lu_dice)
            endif
            write(0,*) nint(calls), 'ci_ice_recon: all methods failed'
         endif
        
!     The real reduced composition reconstruction build on       
         xe(1:nrc)= r_e
         call h2hs_n(h_e, href_n, r_e(1:nrs+ne), nrs+ne, hs_e_n )
         xe(nrc+1)=hs_e_n
!     ---------------------------------------------       
         
         zin(1:ns) = z_e_ICE(1:ns)
         zin(nsp1) = T_e_ICE

!---- store the zu values for temperature approximation
         zu(1:nus) = zin(US)
         call store_zu(zu)

!---- determine the reaction mapping; call rmap2 -------------
         
         call rmap2( 0, nsp1, nsp3, zin, press, dt, zr )

      endif
      
      if (.not.call_ICE_ODE .and. need(2) /=0)  
     1     call rmap2( 1, nsp1, nsp3, zin, press, dt, zr )
      
!---- successful integration:  transform results  ---------------
      do i = 1, ns              !  enforce realizability 
         zr(i,1)  = max( zr(i,1), 0.d0 )
      end do

      rusr(1:nrs+ne) = matmul(BBT(1:nrs+ne, 1:ns), zr(1:ns,1))
      
      call h2hs_n(hin,href_n,rusr(1:nrs+ne), nrs+ne, hs_n )
      
      T_DR = zr(nsp1,1) ! temperature
      rusr(nrs+ne+1) = hs_n
      rusr(nrs+ne+2) = T_DR

!---- set hvar(1) to T_DR(dt)
      hvar(1) = T_DR

!---- approximate temperate if needed
         dpt(2) = press
         call cidpt_dr( rusr(1:nrc+1), dpt )
         T_a = dpt(3)
         rusr(nrs+ne+2) = T_a ! approximated temperature
      
!---- using the constant approximation for {r, h} ---------------
      f(1:nrs+ne+2)= rusr(1:nrs+ne+2) 	

!---- diagnostic output
      if( nh >= 3 )  then
         hvar(2) = info_ice
         hvar(3) = k_facet
      endif

      if( need(2) == 0 ) then
         x = xe
         return  
      endif            
      
!  nc = (nrs+ne) + 1 = nrc + 1
! ------- compute the sensitivity matrix dfdx  --------------
!  The sensitivity matrix involves 4 parts:
!  1: Recovery of enthalpy: ROE(nc+2, nc+2)
!      INPUT : x  =  {r, zu^e, hs^n, p, dt}
!      OUTPUT: y2 =  {r, zu^e, h, p, dt}
!  2: Species reconstruction: SR(ns+3, nc+2)
!      INPUT : y2 =  {r, zu^e, h, p, dt}
!      OUTPUT: y3 =  {z^ICE, T^ICE, p, dt}
!  3. Reaction mapping: RM(ns+1,ns+3) (A_ODE)
!      INPUT : y3 =  {z^ICE, T^ICE, p, dt}
!      OUTPUT: y4 =  {z^R, T^R}
!  4. Reduction : RED(nc+1,ns+1)
!      INPUT : y4 =  {z^R, T^R}
!      OUTPUT: f  =  {r^R, zu^eR, hs^nR, T^R}
!  5. Approximation: APP(nc+1,nc+1)
!      INPUT:  f  =  {r^R, zu^eR, hs^nR, T^R}
!      OUTPUT: fa =  {r^R, zu^eR, hs^nR, T^a}

!  The final sensitivity matrix  A = RED * RM * SR * ROE
!------------------------------------------------------------

!Part 1: Recovery of enthalpy
      allocate(ROE(nc+2,nc+2))
      ROE=0.d0
      
      do i=1, nc+2
         ROE(i,i)=1.d0
      enddo
      ROE(nc,1:nc-1)=href_n(1:nc-1)

!Part 2: Species reconstruction
      allocate(A_ICE(ns+1,nrc+2))
      A_ICE=0.d0

      allocate(SR(ns+3,nc+2))
      SR=0.d0
           
      ! d{z^ICE, T^ICE}/d{r, zu^e, h, p}

      if( info_ice < 0 ) then
       ! for unsuccessful cases use CEM_Tan
       call ci_cem_tan(z_CE(1:ns), hin, T_CE, press, thermo_ns, ns,
     1        nrc, BB, A_ICE, info)
      else
         if(k_facet <= nrc) then
            r_normal = 0.d0
            r_normal(k_facet) = -1.d0
         endif
         if( info_ice > 1 ) k_facet = 0 !attractive

       call ci_ice_tan(z_g(1:ns), T_g, z_e_ICE, T_e_ICE, hin, press,
     1        A_ICE_ODE, BB, k_facet, r_normal, A_ICE, info)
      endif

      if (info < 0) then  
         call isat_abort( 'cirmap_dr_ice_pic', 1,
     1        mess=' failure in calculating tangent space ' )
      endif

      SR(1:ns+1,1:nc+1) = A_ICE

      ! dp/dp
      SR(ns+2,nc+1) = 1.d0
         
      ! dt = dt
      SR(ns+3,nc+2) = 1.d0

!Part 3: Reaction mapping: A_ODE
      allocate(RM(ns+1,ns+3))
      RM = 0.d0

      RM = zr(1:nsp1,2:nsp4)

!Part 4: Reduction
      allocate(RED(nc+1,ns+1))
      RED = 0.d0
      
      ! {rR, zu^eR} = B^T * {z^R}
      RED(1:nrc,1:ns) = transpose(BB)

      ! hs^nR
      call phi2y( zr(1:ns,1), amolwt, ns, y )
      call ctcpbs( zr(nsp1,1), y,   cp1, gas )
      call cthml( zr(nsp1,1),   enth1, gas )
      
      tmp = matmul(href_n(1:nrc), RED(1:nrc,1:ns))
      do i = 1,ns
         RED(nc,i) = enth1(i) - tmp(i)
      end do

      RED(nc,ns+1) = cp1

      ! T^R = T^R
      RED(nc+1,ns+1) = 1.d0

! Part 5: Approximation
      allocate(APP(nc+1,nc+1))
      APP = 0.d0

      ! fa(1:nc) = f(1:nc)
      do i=1, nc
         APP(i,i)=1.d0
      enddo

      ! d{T^a}/d{r^R} = 0.d0 
      ! => APP(nc+1,1:nrs) = 0.d0
     
      ! d{T^a}/d{zu^eR} = 1/cpa*(href_n(nrs+1:nrc) - AT*h^u)^T
      call cpa_dr( T_a, f(1:nrc), cpa )

      call cthml( T_a,   hi, gas )
      hua = matmul(PT, hi(US))
      hri = hi(CS)

      do i=1,nrs
         APP(nc+1,i) = 1.d0/cpa*(href_n(i) - hri(i))
      enddo

      do i=1,ne
         APP(nc+1,nrs+i) = 1.d0/cpa*(href_n(nrs+i) - hua(i))
      enddo

      ! d{T^a}/d{hs^nR} = 1/cpa
      APP(nc+1,nc) = 1.d0/cpa

      ! d{fa}/d{T^R}  = 0.d0 => APP(:,nc+1) = 0.d0

! DONE: Construct full matrix A
      allocate(A(nc+1,nc+2))
      A= matmul(APP, matmul(RED, matmul(RM, matmul(SR, ROE))))

      if( mode_pdt == 1 ) then
         dfdx(1:nc+1,1:nc+2) = A(1:nc+1,1:nc+2) ! r, h_n, p, dt
      elseif( mode_pdt == 2 ) then
         dfdx(1:nc+1,1:nc+1) = A(1:nc+1,1:nc+1) ! r, h_n, p
      elseif( mode_pdt == 3 ) then
         dfdx(1:nc+1,1:nc) = A(1:nc+1,1:nc) ! r, h_n,
         dfdx(1:nc+1,nc+1) = A(1:nc+1,nc+2) !          dt
      else
         dfdx(1:nc+1,1:nc) = A(1:nc+1,1:nc) ! phi, h,
      endif
      
      if( allocated(ROE))       deallocate(ROE)
      if( allocated(SR))        deallocate(SR)
      if( allocated(RM))        deallocate(RM)
      if( allocated(A_ICE))     deallocate(A_ICE)
      if( allocated(RED))       deallocate(RED)
      if( allocated(A))         deallocate(A)

!---- set hvar(1) to T_DR(dt)
      hvar(1) = T_DR

!---- diagnostic output
      if( nh >= 3 )  then
         hvar(2) = info_ice
         hvar(3) = k_facet
      endif
       
!---- return linear approximation to f(x)
      fnew = f + matmul( dfdx, x-xe )
      if( minval(fnew) >= 0.d0 ) then
         f = fnew
      endif

!---- enforce normalization
      f(1:nrc) = f(1:nrc)/dot_product(f(1:nrc), amolwt_n)

      x = xe
      return	
       
!===============  internal subroutines   =====================================

      contains
      
      subroutine bad_nxDR_ICE( n1, n2 ) !-----------------------------------
      
      integer :: n1, n2
      
      call isat_abort('cirmap_dr_ice_pic', 3, mess = 'bad value of nx',
     1     ivar = (/ n1, n2 /) )
      
      end subroutine bad_nxDR_ICE
      
      end subroutine cirmap_dr_ice_pic
