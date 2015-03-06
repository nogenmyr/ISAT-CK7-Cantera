!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

      subroutine cirmap_dr_rcce( need, nx, x, nf, nh, iusr, rusr,
     1	                  f, dfdx, hvar )

!  Determine the reaction mapping (f), its gradient (dfdx), and additional 
!  variables (hvar) as a function of the ISAT tabulation parameters (x). 
!  This version is for the reduced representation provided by RCCE combined with ISAT.

!  definitions:
!   x = r(0), h_sn(0), p, dt  for mode_pdt = 1, general case,   nx=nc+2
!   x = r(0), h_sn(0), p      for mode_pdt = 2, const dt,       nx=nc+1
!   x = r(0), h_sn(0), dt     for mode_pdt = 3, const p,        nx=nc+1
!   x = r(0), h_sn(0)         for mode_pdt = 4, const p and dt, nx=nc
!   r(0): represented species + elements in unrepresented species
!   nc= nrs+ne+1
!   f  = r(dt), h_sn(dt), T(dt);   nf = nc+1

!   hvar - none

!  input:
!       need(1)	= 1 - f is needed (otherwise need(1) = 0 )
!       need(2)	= 1 - df/dx is needed
!       need(3)	= 1 - hvar is needed
!	nx	- number of parameters
!	x	- parameters 
!	nf	- number of dependent variables  (nf=nsp2)
!	nh	- number of additional variables (nh=1)
!  output:
!	f	- variables
!	dfdx	- df/dx
!       hvar    - additional variable (hvar(1) = T_DR(dt))
!
!  Note: It is assumed that the chemistry (and hence nx and nf) are the same on each call.
!        It is assumed that a call with need=(/1,0,0/) is made prior to a call with need(2)>0.

      use ci_dat
      use ci_dat8
      use ci_rmap
      use ci_cem_recon
      use ci_utils
      use ci_stats
      use ci_dpt_dr

      implicit none

      integer, intent(in)     :: need(3), nx, nf, nh, iusr(1)
      real(k_dp), intent(in)  :: x(nx), rusr(1)
      real(k_dp), intent(out) :: f(nf), dfdx(nf,nx), hvar(nh)
      real(k_dp)              :: press, dt
      
! local arrays
      logical :: repeat_call
      real(k_dp), allocatable :: ROE(:,:), SR(:,:), RM(:,:), RED(:,:), 
     1     APP(:,:)
      real(k_dp), allocatable :: A(:,:), A_CEM(:,:)
      
      integer    :: index_pos(ns), info, i
      real(k_dp) :: rin(nrs+ne), hin, hsin_n, hs_n
      real(k_dp) :: y(ns), cp1, enth1(ns), tmp(ns), CBg(1,1), 
     1     zu(nus), stats(20)
      
      real(k_dp) :: hi(ns), hri(nrs), cpi(ns), hua(ne), cpa, dpt(3)

! saved quantities
      integer, save    :: nx_last = -2
      real(k_dp), save :: T_CE, T_a, T_DR, hsin_n_last, calls=0.d0  
      real(k_dp), save, allocatable ::  x_last(:), f_last(:), 
     1            rin_last(:), zin(:),zr(:,:),  phi_CE(:)     
      
!----------------------  start of execution  -------------------------------
      if( nx_last == -2 ) then
         nx_last = -1
         allocate( x_last(nx) )
         allocate( f_last(nf) )
         allocate( rin_last(nrs+ne) )
         allocate( phi_CE(ns) )
         allocate( zin(nsp1) )
         allocate( zr(nsp1,nsp4) ) 
      endif   

      calls = calls + 1.d0
      if( need(1) == 0  .and.  need(2) == 0 ) return

!------ check input  and  determine  press  and  dt-----------------------

	if( nf /= nc+1 ) call isat_abort('cirmap_dr_rcce', 1,
     1	      mess = 'bad value of nf', ivar = (/ nf, nc+1 /) )
	if( nh < 0 ) call isat_abort('cirmap_dr_rcce', 2,
     1	      mess = 'bad value of nh', isv = nh )
     
	if( mode_pdt == 1 ) then   !  variable  p  and  dt
	   if( nx /= nc+2 ) call bad_nxDR( nx, nc+2 )
     	   press = x(nc+1)
	   dt    = x(nc+2)
	elseif( mode_pdt == 2 ) then   !  variable  p, const  dt
	   if( nx /= nc+1 ) call bad_nxDR( nx, nc+1 )
     	   press = x(nc+1)
	   dt    = dtc
	elseif( mode_pdt == 3 ) then   !  const  p,  variable  dt
	   if( nx /= nc+1 ) call bad_nxDR( nx, nc+1 )
	   press = prc
     	   dt    = x(nc+1)
	elseif( mode_pdt == 4 ) then   !  const  p,  and  dt
	   if( nx /= nc ) call bad_nxDR( nx, nc )
	   press = prc
     	   dt    = dtc
      else
	   call isat_abort('cirmap_dr_rcce', 4, mess 
     1 	         = 'bad value of mode_pdt', isv = mode_pdt )
      endif
     	
!---- determine if x is the same as on the last call
      
      repeat_call = .false.
      if( nx == nx_last ) then
         if( maxval( abs( x(1:nx) - x_last(1:nx) ) ) == 0.d0 ) then        
            repeat_call = .true.
            f(1:nf)     = f_last(1:nf)
            if( need(2) == 0 ) return !  only f required
              
!---- retrieve information from last call: note that zin and zr are saved           
            rin    = rin_last
            hsin_n = hsin_n_last          
         endif
      endif

      if( .not.repeat_call ) then  
!---- perform species reconstruction (if not saved)
         index_pos = (/(i, i=1,ns)/)
	   do i = 1, nrs+ne          !  enforce realizability 
	      rin(i)  = max( x(i), 0.d0 )
	   end do
           
         hsin_n = x(nrs+ne+1)  ! nominal sensible enthalpy         
         call hs2h_n( hsin_n, href_n, rin, nrc, hin)  ! hin=enthalpy
         
!---- species reconstruction by RCCE
         call ci_ceq( 0, rin, hin, press, .false., phi_CE, T_CE, 
     1                      phi_CE, T_CE,  
     2                      stats, info )
       	
         zin(1:ns) = phi_CE(1:ns)
         zin(nsp1) = T_CE
     
         if( info < 0 ) call isat_abort( 'cirmap_dr_rcce', 1, mess=
     1	  'Constrained equilibrium calculation failed, info = ',
     2      isv = info )

!---- store the zu values for temperature approximation
         zu(1:nus) = zin(US)
         call store_zu(zu)

!---- perform the reaction mapping
         call rmap2( 0, nsp1, nsp3, zin, press, dt, zr ) 

!---- extract f from solution zr  ---------------
         do i = 1, ns       !  enforce realizability 
            zr(i,1)  = max( zr(i,1), 0.d0 )
         end do

         f(1:nrs+ne) = matmul( BBT(1:nrs+ne, 1:ns), zr(1:ns,1) )
            
         call h2hs_n( hin, href_n, f(1:nrs+ne), nrs+ne, hs_n )  

         T_DR = zr(nsp1,1)  ! temperature
         f(nrs+ne+1) = hs_n ! nominal sensible enthalpy
         f(nrs+ne+2) = T_DR

!---- set hvar(1) to T_DR(dt)
         hvar(1) = T_DR

!---- approximate temperate if needed
         dpt(2) = press
         call cidpt_dr( f(1:nrc+1), dpt )
         T_a = dpt(3)
         f(nrs+ne+2) = T_a ! approximated temperature

!---- save information for subsequent call 
         nx_last      = nx
         rin_last     = rin
         hsin_n_last  = hsin_n
         x_last(1:nx) = x(1:nx)
         f_last(1:nf) = f(1:nf)  
      endif 

! Gradient matrix          		
      if( need(2) == 0 ) return
      
! Get mapping gradient
      call rmap2( 1, nsp1, nsp3, zin, press, dt, zr ) 
      do i = 1, ns              !  enforce realizability 
         zr(i,1)  = max( zr(i,1), 0.d0 )
      end do


!  nc = (nrs+ne) + 1 = nrc + 1
! ------- compute the sensitivity matrix dfdx  --------------
!  The sensitivity matrix involves 4 parts:
!  1: Recovery of enthalpy: ROE(nc+2, nc+2)
!      INPUT : x  =  {r, zu^e, hs^n, p, dt}
!      OUTPUT: y2 =  {r, zu^e, h, p, dt}
!  2: Species reconstruction: SR(ns+3, nc+2)
!      INPUT : y2 =  {r, zu^e, h, p, dt}
!      OUTPUT: y3 =  {z^CE, T^CE, p, dt}
!  3. Reaction mapping: RM(ns+1,ns+3) (A_ODE)
!      INPUT : y3 =  {z^CE, T^CE, p, dt}
!      OUTPUT: y4 =  {z^R, T^R}
!  4. Reduction : RED(nc+1,ns+1)
!      INPUT : y4 =  {z^R, T^R}
!      OUTPUT: f  =  {r^R, zu^eR, hs^nR, T^R}
!  5. Approximation: APP(nc+1,nc+1)
!      INPUT:  f  =  {r^R, zu^eR, hs^nR, T^R}
!      OUTPUT: fa =  {r^R, zu^eR, hs^nR, T^a}

!  The final sensitivity matrix  A = APP * RED * RM * SR * ROE
!------------------------------------------------------------

!Part 1: Recovery of enthalpy
      allocate(ROE(nc+2,nc+2))
      ROE=0.d0

      do i=1, nc+2
         ROE(i,i)=1.d0
      enddo
      ROE(nc,1:nc-1)=href_n(1:nc-1)

!Part 2: Species reconstruction
      allocate(A_CEM(ns+1,nrc+2))
      A_CEM=0.d0

      allocate(SR(ns+3,nc+2))
      SR=0.d0
           
      ! d{z^CE, T^CE}/d{r, zu^e, h, p}
      call ci_cem_tan(phi_CE(1:ns), hin, T_CE, press, thermo_ns, 
     1     ns, nrc, BB, A_CEM, info) 

      if (info < 0) then  
         call isat_abort( 'cirmap_dr_rcce', 1,
     1        mess=' failure in calculating tangent space ' )
      endif

      SR(1:ns+1,1:nc+1) = A_CEM

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
      hua = matmul( PT, hi(US) )
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
      A = matmul(APP, matmul(RED, matmul(RM, matmul(SR, ROE))))

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
      if( allocated(A_CEM))     deallocate(A_CEM)
      if( allocated(RED))       deallocate(RED)
      if( allocated(A))         deallocate(A)

!---- set hvar(1) to T_DR(dt)
      hvar(1) = T_DR
   
      return

!===============  internal subroutines   =====================================

      contains

      subroutine bad_nxDR( n1, n2 )

      integer :: n1, n2

      call isat_abort('cirmap_dr_rcce', 3, mess = 'bad value of nx',
     1         ivar = (/ n1, n2 /) )
      
      end subroutine bad_nxDR

      end subroutine cirmap_dr_rcce
