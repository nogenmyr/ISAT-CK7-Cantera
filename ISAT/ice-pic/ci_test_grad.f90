!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ci_test_grad( nx, nf, kr_out )

! Test the accuracy of the mapping gradient.
! Read (adds.op) values of x, {xa, xb}.
! Determine: fa and ga.
! For xc(s) = xa + s*(xb-xa), (for 0<s<1) evaluate fc(s), and the linear approximation
!  fcl(s) = fa + ga*(xc-xa).  
! Evaluate the ISAT error errf = |(fcl-fc)./fref| and output on grad_e.op.
! For the selected record, kr_out, evaluate the gradients (fc-fa)/s and (fcl-fa)/s.
! Output to grad.op the max abs scaled errors.

! Input:
!   nx - number of components of x
!   nf - number of components of f
!   kr_out - record number for which output to be generated

use ci_8, only: amolwt_n, nrc
use ci_dat,  only: modeci
use ci_dat6, only: phiref, hsfref, tempref
use ci_utils

implicit none
integer, intent(in) :: nx, nf, kr_out

integer, parameter  :: ns=60      ! number of values of s (ns >=1)
real(kind(1.d0))    :: s0=1.d-6   ! initial value of s (s0<=s<=smax)
real(kind(1.d0))    :: smax=1.d0   ! initial value of s (s0<=s<=smax)
integer             :: nr_min=1, nr_max = 50 ! first and last record to be treated

integer :: i, iref, nr, need(3), nh, iusr(1), ice_stats, j, jmax, nr_errmax, &
                    lu1, lu2, lu3, lu4, k_facet_a, k_facet_c, log, usedefault
real(kind(1.d0)) :: inc, amp, xa(nx), xb(nx), xe(nx), xc(nx), xd(nx), hvar(3), dfdx(nf,nx), s, &
                    fa(nf), fc(nf), fe(nf), fcl(nf), dfds(nf), dfdsl(nf), err(nf), &
                    fref(nf), errf(ns), errmax, errmax_max, sref, ss(ns)

iref = (ns+1)/2
iref = 2           ! index of s at which error is evaluated 
iref = ns !XXX

! log scale or linear scale
log = 1

if(log==1) then
   amp  = (smax/s0)**(1.d0/(ns-1))
   sref = s0*amp**(iref-1)  !  value of s at which error is evaluated
else
   inc = (smax - s0)/(ns - 1)
   sref = s0 + inc*(iref-1)
endif

if(modeci==9) then 
   nh = 3
else
   nh = 0
endif

errmax_max = -1.d0

fref       = phiref
fref(nf-1) = hsfref
fref(nf)   = tempref

write(0,*) 'nr, ice, k_a, k_c, jmax, Tc, errf, errmax'
call isat_lu( lu1 )
call isat_lu( lu2 )
call isat_lu( lu3 )
call isat_lu( lu4 )
open( lu1, file='adds.op' )
open( lu2, file='grad.op' )
open( lu3, file='grad_e.op' )
open( lu4, file='grad_s.op' )
if(log==1) then
   write(lu2,'(1p,200e25.16)') fref, s0, amp, sref, fref(4:nf), fref, fref
else
   write(lu2,'(1p,200e25.16)') fref, s0, inc, sref, fref(4:nf), fref, fref
endif

! set default x
xd(1:nrc) = 1.0/sum(amolwt_n)
xd(nrc+1) = hsfref
usedefault = 0
! write(0,*) xd


do nr = 1, nr_max
   read(lu1,*,end=100,err=100) xa, xb
   if( nr < nr_min ) cycle
   
   if(usedefault == 1) xa = xd

! evaluate fa and ga 
   need    = 0
   need(1) = 1
   xe      = xa
   hvar(1) = -1
   
   if( modeci == 7 ) then
      call cirmap1( need, nx, xe, nf, nh, iusr, fe, fa, dfdx, hvar )
      need = 1
      call cirmap1( need, nx, xe, nf, nh, iusr, fe, fa, dfdx, hvar )
   
   elseif( modeci == 8 ) then
      call cirmap_dr_rcce( need, nx, xe, nf, nh, iusr, fe, fa, dfdx, hvar )
      need = 1
      call cirmap_dr_rcce( need, nx, xe, nf, nh, iusr, fe, fa, dfdx, hvar )
      
   elseif( modeci == 9 ) then
      fe(1) = nf
      fe(2) = -huge(1.d0)  !  tol_re
      call cirmap_dr_ice_pic( need, nx, xe, nf, nh, iusr, fe, fa, dfdx, hvar )
      need = 1
      xe = xa
      fe(1) = nf
      fe(2) = -huge(1.d0)  !  tol_re
      call cirmap_dr_ice_pic( need, nx, xe, nf, nh, iusr, fe, fa, dfdx, hvar )
      ! Varun: we need x, f(x) and df/dx at same x, so use x = xe
      xa = xe
      fa = fe
   endif
   
   ice_stats = nint(hvar(2))
   k_facet_a = nint(hvar(3))
   need(2)   = 0

   do i = 1, ns
      if(log==1) then
         s = s0 * amp**(i-1)
      else
         s = s0 + inc*(i-1)
      endif

      ss(i) = s
      xc = xa + s*(xb-xa)
      xe = xc

      if( modeci == 7 ) then
         call cirmap1( need, nx, xe, nf, nh, iusr, fe, fc, dfdx, hvar )
         
      elseif( modeci == 8 ) then
         call cirmap_dr_rcce( need, nx, xe, nf, nh, iusr, fe, fc, dfdx, hvar )
      
      elseif( modeci == 9 ) then
         fe(1) = nf
         fe(2) = -huge(1.d0)  !  tol_re
         call cirmap_dr_ice_pic( need, nx, xe, nf, nh, iusr, fe, fc, dfdx, hvar )
         ! Varun: calculate fcl at xe and define new ss
         xc = xe
         ss(i) = norm(xc - xa)/norm(xb - xa)
      endif
      
      ice_stats = nint(hvar(2))
      k_facet_c = nint(hvar(3))
   
      fcl     = fa + matmul( dfdx, xc-xa )
      errf(i) = norm( (fcl-fc)/fref )
      
      dfds  = (fc-fa) /s
      dfdsl = (fcl-fa)/s
      
      if( nr == kr_out ) write(lu2,'(1p,200e25.16)') fa, fc, dfds, dfdsl

      write(0,'(5i4,1p,20e13.4)') nr-1, i, ice_stats, k_facet_a, k_facet_c, jmax, &
           fc(nf), errf(i), errmax
     
      if( i==iref ) then
         errmax = -1.d0
         do j = 1, nf
            err(j) = abs( dfds(j) - dfdsl(j) ) / fref(j) ! max abs fractional error
            if( err(j) > errmax ) then
               errmax = err(j)
               jmax = j
            endif
         end do
              
         if( errf(i) > errmax_max ) then
            errmax_max = errf(i)
            nr_errmax     = nr
         endif
      endif
            
   end do
   write(lu4,'(1p,200e25.16)') ss
   write(lu3,'(1p,200e25.16)') errf

end do

100  continue
   
write(0,'(a,i6,1p,e13.4)')'maximum error: record, error = ', nr_errmax, errmax_max
   
return
end subroutine ci_test_grad
