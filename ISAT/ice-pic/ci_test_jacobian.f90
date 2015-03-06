!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ci_test_jacobian( nx, nf, kr_out )
  ! Test the jacobian J = dS/dz evaluated at constant enthalpy

  use ci_dat8
  use ci_cksubs
  use ci_utils
  use ci_cem_recon

  implicit none
  
  integer, intent(in) :: nx, nf, kr_out

  ! local variables
  integer, parameter  :: np=60      ! number of values of s (ns >=1)
  real(kind(1.d0))    :: tau=1.d-7  ! some value of tau
  real(kind(1.d0))    :: s0=1.d-30   ! initial value of s (s0<=s<=smax)
  real(kind(1.d0))    :: smax=1.d0  ! initial value of s (s0<=s<=smax)
  integer             :: nr_min=1, nr_max = 50 ! first and last record to be treated

  integer    :: i, k, info, lui, luz, luJ, luT, luS, logz, luout = 0

  real(k_dp) :: incz, ampz, s, sref, xa(nx), xb(nx), zc(ns), y(ns), yc(ns), ss(np)
  real(k_dp) :: p, hsin_n, ha, hb, za_CE(ns), zb_CE(ns), Ta_CE, Tb_CE, Tc, CBg(1,1), index_pos(ns), ceq_stats(20)
  real(k_dp) :: dz(ns), J(ns, ns), Sr(ns), Srn(ns), errJ(np), errT(np), errS(np)

  call isat_lu( lui )
  call isat_lu( luz )
  call isat_lu( luJ )
  call isat_lu( luT )
  call isat_lu( luS )

  open( lui, file='adds.op' )
  open( luz, file='dzpts.op' )
  open( luJ, file='err_J.op' )
  open( luS, file='err_S.op' )
  open( luT, file='err_T.op' )

  ! set p to prc
  p = prc

  ! log scale or linear scale
  logz = 1

  if(logz==1) then
     ampz  = (smax/s0)**(1.d0/(np-1))
     do i = 1, np
        ss(i)  = s0 * ampz**(i-1)
     enddo
  else
     incz = (smax - s0)/(np - 1)
     do i = 1, np
        ss(i)  = s0 + incz*(i-1)
     enddo
  endif

  ! store s and dt values in lmap_st.op
  write( luz,'(1p,200e25.16)' ) ss

  do nr = 1, nr_max
     read(lui,*,end=100,err=100) xa, xb
     if( nr < nr_min ) cycle

     index_pos = (/(i, i=1,ns)/)

     hsin_n = xa(nrs+ne+1)
     call hs2h_n( hsin_n, href_n, xa(1:nrc), nrc, ha)

     ! compute za_CE = CE(xa)
     call ci_ceq( 0, xa(1:nrc), ha, p, .false., za_CE, Ta_CE, & 
          za_CE, Ta_CE, ceq_stats, info, luout)

     hsin_n = xb(nrs+ne+1)
     call hs2h_n( hsin_n, href_n, xb(1:nrc), nrc, hb)

     ! compute zb_CE = CE(xb)
     call ci_ceq( 0, xb(1:nrc), hb, p, .false., zb_CE, Tb_CE, & 
          zb_CE, Tb_CE, ceq_stats, info, luout)
     
     call phi2y( za_CE, amolwt, ns, yc )
     call temphy( ha, yc, Ta_CE )

     ! S at za_CE
     call ciS( za_CE, Ta_CE, p, Sr )

     ! Jacobian, J = dS/dz
     call jacobian( tau, ns, za_CE, p, Ta_CE, J )

     do i = 1, np
        s = ss(i)

        zc = za_CE + s*(zb_CE - za_CE)
        call phi2y( zc, amolwt, ns, yc )
        call temphy( ha, yc, Tc )
        dz = s*(zb_CE - za_CE)

        ! S at zc
        call ciS( zc, Tc, p, Srn )
        
        ! compute error in J, T and S
        errJ(i) = norm(Srn - (Sr + matmul(J, dz))) ! absolute error
        errT(i) = abs(Tc - Ta_CE)/Ta_CE ! relative error
        errS(i) = norm(Srn - Sr) ! absolute error

     enddo

     write( luJ,'(1p,200e25.16)' ) errJ
     write( luT,'(1p,200e25.16)' ) errT
     write( luS,'(1p,200e25.16)' ) errS

  enddo

100 continue
  return

end subroutine ci_test_jacobian
