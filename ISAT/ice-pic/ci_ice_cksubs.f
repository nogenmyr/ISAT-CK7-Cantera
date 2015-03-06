!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

	module ci_ice_cksubs

!  This module conatins various subroutines that evaluate thermochemical
!  properties.  Some call Chemkin routines.
!  Subroutines:

!  ci_chem_ddassac(map,z0,press,dt,Tdt,zt,dzdz0,iflag_ODE)
!  ci_chem_ddassac5(nsp1, tstart, z, zp, tout, ...

	use ci_dat
	use ci_dat8
	use isat_abort_m
	use ci_cksubs
	use ci_utils
	use ci_prec

	implicit none

	contains

!==========================================================================

      subroutine ci_chem_ddassac(map,z0,press,dt,Tdt,zt,dzdz0,iflag_ODE)
!--------------------------------------------------------------------------
!  Determine the reaction mapping zt, its gradient dzdz0, 
!     given z0,dt, pressure .
!  definitions:
!     const pressure and dt
!     z0, zt = {phi, h} : ns specific mole numbers and enthalpy (ns+1)      
!     dzdz0 : sensitivity matrix  (fixed enthalpy)  (ns+1)*(ns+1)
!---------------------------------------------------------------------------
!  Input:
!       map	= 0  - Tdt, zt  are needed 
!          	= 1  - Tdt, zt, dzdz0 are needed
!       z0	     - initial composition  
!       press	 - pressure (cm-g-s), constant
!       dt       - time step (sec)
!
!  Output:
!	  Tdt	      - final temperature
!	  z	      -  reaction mapping composition
!       dzdz0     - sensitivity matrix
!       iflag_ODE - integration succesuful or not
!---------------------------------------------------------------------------
   
	implicit real(kind(1.d0)) (a-h, o-z), integer (i-n)
     
      integer, intent(in)          :: map
	real(kind(1.d0)), intent(in) :: z0(ns+1), press, dt
	real(kind(1.d0)), intent(out):: Tdt, zt(ns+1), dzdz0(ns+1,ns+1)


!local array
  
      integer          :: ii, jj, info(18), ipar(ns+2)
      real(kind(1.d0)) :: rpar(ns+2)
	real(kind(1.d0)) :: zstart(ns+1)
      real(kind(1.d0)) :: z(ns+1,ns+2),  zp(ns+1,ns+2),
     1                    zst(ns+1,ns+2),zpst(ns+1,ns+2)  
	real(kind(1.d0)) :: rtol(ns+1, ns+2),  atol(ns+1,ns+2),
     1                    rtolt(ns+1, ns+2), atolt(ns+1,ns+2)   

      real(kind(1.d0)) :: y(ns), enth0(ns),enth1(ns) 

	integer, save                 :: liwdas, lrwdas
     	integer, allocatable, save    :: iwork(:)
     	real(kind(1.d0)), allocatable, save :: rwork(:)

      data ifst/0/
      data ifjac, nonneg, knewt/ 1, 0, 0 /      
	data lu / -1 /

      iflag_ODE=0
!---------  on first call allocate work arrays for ddasac  --------------------
     
	if( ifst == 0 ) then
	   ifst = 1
	   liwdas = 20 + (ns+1) * (ns+2)
	   lrwdas = 40+10*(ns+1)*(ns+2)+(ns+1)**2
	   allocate( iwork(liwdas) )
	   allocate( rwork(lrwdas) )
	endif
!--------- open file for ddasac diagnostic output  -----------------

        if( lu .eq. -1 ) then
           call isat_lu( lu )
           open( lu , file = 'ddasac.op' )
       endif

!--------  set sensitivities to the identity  ---------------------
c
        if( map .eq. 1 ) then
         do  ii   = 1, ns+1
           do  jj   = 2, ns+2
           z(ii,jj)   = 0.d0
	     enddo
           z(ii,ii+1) = 1.d0
	   enddo
        endif	

	   nover = 0
	   ifxh0 = 0
	   tleft = dt

!--------copy z0 to z(:,1) and zstart
        call phi2y(z0,amolwt,ns,y)
         h=z0(ns+1)
	  call temphy(h,y,temp)
	 do  ii   = 1, ns
	  z(ii,1)    = z0(ii)
	 enddo
        z(ns+1,1)=temp
	
10	 continue
      do  ii   = 1, ns+1
	z(ii,1)    = max( z(ii,1), 0.d0 )
	zstart(ii) = z(ii,1)
      enddo
      
      zp = 0.d0 ! SBP XXX
      
!    so the last component in z and zstart is temperature instead of enthalpy.

!---------------  initial Cp and enthalpies  ----------------------

	temp    = z(ns+1,1)
	
	call phi2y( z, amolwt, ns, y )
	call ctcpbs( temp, y,   cp0, gas )
	call cthml( temp,   enth0, gas )

!------------  prepare for call to ddasac  ------------------------------

	do  ii = 1, 18
	iwork(ii) = 0
	info(ii) = 0
      enddo
c
	do  ii = 1, 44
	rwork(ii) = 0.d0
      enddo
c
        info(2)  = 1
        info(3)  = 0
        info(5)  = ifjac
        info(10) = nonneg
        info(11) = 1
        info(12) = map * (ns + 1)
        info(14) = 1
        info(18) = knewt
c
        rwork(3) = 1.d0
c
        tstart  = 0.d0
        tout    = tleft
c

!  fix initial time step size h0?
!
	if( ifxh0 .eq. 1 ) then
	   info(8)  = 1
	   rwork(3) = 1.1 * tout
	   ifxh0    = 0
	endif
!
!  set tolerances
c
	do  ii   = 1, ns+1
	atol(ii,1) = atolc * phiref
	rtol(ii,1) = rtolc
      enddo
!
        atol(ns+1,1) = atolc * tempref
c
c  additional info for sensitivities ( map = 1 )
c
        if( map .eq. 1 ) then
          do ii = 1, ns+1
           ipar(ii) = ii
           rpar(ii) = z(ii,1)
           do jj = 2, ns+2
           atol(ii,jj) = atol(ii,1)
           rtol(ii,jj) = rtol(ii,1)
           enddo
	    enddo
        endif  
!---------  set pressure (set to rpar(nspar) in cires5 etc.)
         ipar(ns+2) = ns+2
	   rpar(ns+2) = press

          idid=0
c------------  loop over tries at integration  --------------------
c

110     call ci_chem_ddassac5(ns+1, tstart, z, zp, tout,
     1               info, rtol, atol, idid,
     1               rwork, lrwdas, iwork, liwdas, rpar,
     1               ipar, ns, zst, zpst, rtolt, atolt)
c
c---------------  diagnostic information  (ddasac.op)  ------------
c

              !   write(*,*) tstart, idid, zst(ns+1,1)
!---if all going well but not finished, keep going
           if(idid.eq.-1) then
	        info(1)=1
		 goto 110
		 
	   elseif (idid .eq.-6) then
	        info(11)=1
		info(1)=0
		goto 110	 
           elseif (idid .lt.0) then
	      write(*,*)'ddasac failure', idid  
              iflag_ODE=1
              return
	   endif   
1000          format(1p,100e16.8)
1002          format(1p, 100e16.8)
1003          format(1i8,100e16.8)
1004          format(400e16.8)

c
c----------------------  successful integration  ---------------
c copy z to zt	
           temp=z(ns+1,1)
	 do   ii =1,ns 
	   zt(ii)=z(ii,1)
	 enddo
!ren           zt(ns+1)=h
         Tdt=temp	
	call phi2y( zt, amolwt, ns, y )
	   
	call cthbms( Tdt, y,   h, gas )
	zt(ns+1)=h 

	if( map .eq. 0 ) return

c  set density	
	call phi2y( z, amolwt, ns, y )
	call ctrhoy( press, temp, y,   dens, gas )	
c
c------  convert from sensitivities w.r.t. temperature to w.r.t. enthalpy ----
c
	do 120 i     = 1, ns+1
	dzdh         = z(i,1+ns+1) / cp0
	z(i,1+ns+1)  = dzdh
c
	do 121 j     = 1, ns
	z(i,1+j)     = z(i,1+j) - enth0(j) * dzdh
121     enddo
120	enddo
c
c------  convert sensitivity of temperature to enthalpy --------------------
c
	call ctcpbs( temp, y,   cp1, gas )
	call cthml( temp,   enth1, gas )
c
	do 140 j    = 1, ns+1
	hdfdz       = 0.
	do 150 i    = 1, ns
	hdfdz       = hdfdz + enth1(i) * z(i,1+j)
150     enddo
	z(ns+1,1+j) = cp1 * z(ns+1,1+j) + hdfdz
140      enddo
c copy the array of z(1:ns+1,2:ns+2) to the sensitivity matrix dzdz0
       do 123 i=1,ns+1
	 do 124 j=1,ns+1
         dzdz0(i,j)=z(i,j+1)
124      enddo
123      enddo
      return
      
      end subroutine ci_chem_ddassac

!==========================================================================

      	subroutine ci_chem_ddassac5(nsp1, tstart, z, zp, tout, 
     1                     info, rtol, atol, idid, 
     1                     rwork, lrwdas, iwork, liwdas, rpar, 
     1                     ipar, ns, ztst, zptst, rtolt, atolt)
c
c  routine to call ddasac with arrays of leading dimension nsp1,
c    then assign the values to arrays with leading dimension ns + 1.
c
      
	implicit real(kind(1.d0)) (a-h, o-z), integer (i-n)
c
	integer info, iwork, idid
	dimension z(nsp1,nsp1+1), zp(nsp1,nsp1+1), iwork(liwdas),
     1	  rwork(lrwdas), info(18), 
     1	  rtol(nsp1,nsp1+1), atol(nsp1,nsp1+1), ipar(nsp1+1),
     1	  rpar(nsp1+1)
c
	dimension ztst(ns+1,ns+2), zptst(ns+1,ns+2), 
     1            rtolt(ns+1,ns+2), atolt(ns+1,ns+2)
c
      external cires5, ciesub5, cijac5,cibsub5 !, cijac5adf
	data lud/ 0 /
c
	do i = 1, ns+1
	do j = 1, ns+2
	ztst(i,j) = z(i,j)
	zptst(i,j) = zp(i,j)
	atolt(i,j) = atol(i,j)
	rtolt(i,j) = rtol(i,j)
	enddo
	enddo

!------------------------------------------------------------------------------       
       if(info(5).eq.1) then
	call ddasac5( tstart, tout, ns+1, ztst, zptst, rtolt, atolt, 
     1                info, rwork, lrwdas, iwork, liwdas, rpar, ipar, 
     2                idid, lud, 0, cires5, ciesub5,cibsub5)!,cijac5adf
       else
        call ddasac5( tstart, tout, ns+1, ztst, zptst, rtolt, atolt, 
     1                info, rwork, lrwdas, iwork, liwdas, rpar, ipar, 
     2                idid, lud, 0, cires5, ciesub5, cijac5, cibsub5 )
       endif
!------------------------------------------------------------------------------

      if (idid==-33) then
         call isat_abort( 'ci_chem_ddassac5', 1,
     1           mess='DDASAC5 failed, idid = -33: t0, t1, t1-t0= ',
     2           rvar = (/tstart, tout, tout-tstart/) )

      endif
       
	do i = 1, ns+1
	do j = 1, ns+2
	z(i,j) = ztst(i,j)
	zp(i,j) = zptst(i,j)
	enddo
	enddo
c
	return
	end subroutine ci_chem_ddassac5

!=====================================================================	
	subroutine temphr( h, r, T )

!  Routine to determine temperature given enthalpy and reduced representation
!
!  input:
!	h  - specific enthalpy of gas mixture
!	r  - the reduced representation
!  output:
!	t  - temperature (K)

	use ci_dat8
	use ci_dpt_dr
	implicit none

	real(k_dp), intent(in)  :: h, r(ns)
	real(k_dp), intent(out) :: T

	! Local variables
	real(k_dp) :: za(ns), zua(nus)
	real(k_dp) :: hlo, hhi, ya(ns)
	integer    :: i, check, info

	! compute zua
	zua = matmul(transpose(PT), r(nrs+1:nrc))

	! set za
	za = 0.d0
	za(CS) = r(1:nrs)
	za(US) = zua(1:nus)

	! call temphz with za
	check = 1
	info = 0
	call temphz( h, za, T, check, info )
	
	if(info /= 0) then
	 if( clipt == 1 ) then
	  if( clipt_log == 1 ) then
	   if( info == 1 ) then
	 write(lu_err,*) myrank, 'temphr: temperature clipped, T > tbadhi'
	   else
	 write(lu_err,*) myrank, 'temphr: temperature clipped, T < tbadlo'
	    endif
	  endif
	  return
	  else
	      call phi2y( za, amolwt, ns, ya )
	      call cthbms( tbadlo, ya,   hlo, gas )
	      call cthbms( tbadhi, ya,   hhi, gas )

	      write(lu_err,*)'------------------------------------'
	      write(lu_err,*)'temphr: temperature out of range'
	      write(lu_err,*)'temphr: tbadlo, tbadhi', tbadlo, tbadhi
	      write(lu_err,*)'temphr: h, hlo, hhi', h, hlo, hhi
	      write(lu_err,*)'species mass fractions'
	      do i = 1, nrc
		 write(lu_err,600) i, r(i)
 600		 format(1x,i5,1p,e13.4)
	      end do
	      
	      write(lu_err,*)'------------------------------------'
	      write(lu_err,*)'Please try one of the following options to
     1 fix this issue:'
	  write(lu_err,*)'(1) Change the represented species'
	  write(lu_err,*)'(2) Increase the number of represented species'
	  write(lu_err,*)'(3) Clip the temperature values between
     1 [tbadlo, tbadhi]'
	      write(lu_err,*)'    by setting clipt = 1 in ci_ext.nml'
	      
	      call isat_abort( 'temphr:', 1, mess='h outside range')
	   endif
	endif
	return

	end subroutine temphr

!=====================================================================

	function exp_safe( x, n )
	
! For each component of x, return min( exp(x), huge(x) )
	
	real(k_dp), dimension(n) :: exp_safe
	integer   , intent(in)   :: n
	real(k_dp), intent(in)   :: x(n)
	real(k_dp), parameter    :: log_huge = log( huge(1.d0))
	integer                  :: i
	
	if( maxval( abs(x) ) < log_huge ) then
	   exp_safe = exp(x)
	   return
	endif
	
	do i = 1, n
	
	  if( x(i) >= log_huge ) then
	     exp_safe(i) = huge(1.d0)
	     	  
	  else
	     exp_safe(i) = exp( x(i) )
	  endif
	  
	end do
	
	return
	end function exp_safe

!=====================================================================

	subroutine set_href_n( k_mode )
	
!  Set the nominal enthalpy of formation for the reduced composition.
!  For represented species:  href_n(i) = href(CS(i)) for i = 1, nrs
!  For elements of unrepresented species ( href_n(nrs+1:nrc) ), treatment depends on k_mode:
!  k_mode  = -2 - set href_n(nrs+1:nrc) = 0
!  k_mode >= -1 - href_n(nrs+1:nrc) is taken as the minimum-norm, least-squares solution to
!     w.*CEu * href_n = w.*href_u, where w is a weighting vector:
!  k_mode = -1     - w = 1
!  k_mode = 0      - w = sum of streams' unrepresented species specific moles
!  k_mode = 1-nstr - w = unrepresented species specific moles of stream k_mode
!

	use ci_dat
	use ci_dat6
	use ci_dat8
	implicit none

	integer, intent(in) :: k_mode
	integer    :: mode, info, i, ist
	real(k_dp) :: w(nus), A(nus,ne), b(nus), S(ne), U(nus,ne),
     1              VT(ne,ne), work(10+10*nus*ne), thresh, Sinv(ne,ne)

	if( k_mode == -2 ) then
	   href_n(nrs+1:nrc) = 0.d0   
	   return
	   
	elseif( k_mode == -1 ) then
	   w = 1.d0		! equal weighting
	   
	elseif( k_mode == 0 ) then !sum of stream compositions
	   do i = 1, nus
	      w(i) = sum( ccstrm_f( US(i),1:nstr) )
	   end do
	   
	elseif( k_mode > 0  .and.  k_mode <= nstr ) then ! composition of stream kwt
	   w = ccstrm_f(US,k_mode)
         
	else
	   call isat_abort('ci_8:set_href_n', 1, 
     1                 mess='bad value of k_mode = ', isv=k_mode )
	endif
      
	do i = 1, nus		!  set up equation as A*x=b
	   A(i,:) = w(i)*CEu(i,:)
	   b(i)   = w(i)*href(nrs+i)
	end do
	
	call dgesvd( 'S', 'S', nus, ne, A, nus, S, U, nus, VT, ne, 
     1              work, 10+10*nus*ne, info )
      
	if( info /= 0 ) call isat_abort('ci_8:set_href_n', 2, 
     1                 mess='dgesvd failed, info = ', isv=info )
     
	Sinv = 0.d0		!  form psuedo-inverse
	thresh = max( 1.d-9, 1.d-6*S(1) )  
	do i = 1, ne
	   if( S(i) > thresh ) Sinv(i,i) = 1.d0 /S(i)
	end do
	
	href_n(nrs+1:nrc) = matmul( b, matmul( U, matmul( Sinv, VT ) ) )
	
	return
	end subroutine set_href_n


	end module ci_ice_cksubs
!=====================================================================
