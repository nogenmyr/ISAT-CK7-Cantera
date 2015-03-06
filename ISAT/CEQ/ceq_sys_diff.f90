!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module ceq_diff

  use ceq_types
  implicit none

contains

  subroutine ceq_sys_diff( sys1, sys2, luout, info, diag, thresh )

    ! Compare two CEQ systems and report any differences.

    ! Input:
    !   sys1 - first system to be compare
    !   sys2 - second system to be compare
    !   luout - logical unit number on which to report any differences
    
    ! Optional:
    !   diag = 0 (default) - minimum output (first difference encountered)
    !   diag = 1 - full output (all differences)
    !   thresh - threshold for differences in reals (thresh >=0.)

    ! Output:
    !   info  = 0 - sys1 and sys2 are identical
    !   info /= 0 - sys1 and sys2 differ

    ! Variables not compared: diag, lu, thermo 

    implicit none

    type (sys_type), pointer      :: sys1, sys2
    integer, intent(in)           :: luout
    integer, intent(in), optional :: diag
    real(kind(1.d0)), intent(in), optional :: thresh
    integer, intent(out)          :: info
    logical                       :: min_diag = .true.
    real(kind(1.d0))              :: eps = 0.d0

    info  = 0
    if( present(diag) )   then
       if( diag /=0 ) min_diag = .false.
    endif
    
    if( present(thresh) ) eps   = thresh

    if(.not.associated(sys1)) then
       write (luout,*)'ceq_sys_diff: Error :: sys1 not associated!'
       info = -2
       return 
    endif

    if(.not.associated(sys2)) then
       write (luout,*)'ceq_sys_diff: Error :: sys2 not associated!'
       info = -2
       return 
    endif

    if(.not.sys1%initialized) then
       write (luout,*)'ceq_sys_diff: Warning :: sys1 not initialized yet!'
       info = -2
    endif

    if(.not.sys2%initialized) then
       write (luout,*)'ceq_sys_diff: Warning :: sys2 not initialized yet!'
       info = -2
    endif

    if(sys1%ne /= sys2%ne) then
       write (luout,100)'ceq_sys_diff:','different ne => ','sys1 ne =', sys1%ne,'sys2 ne =', sys2%ne
       info = -1
       if(min_diag) return 
    endif

    if(sys1%ned /= sys2%ned) then
       write (luout,100)'ceq_sys_diff:','different ned => ','sys1.ned =', sys1%ned,'sys2.ned =', sys2%ned
       info = -1
       if(min_diag) return
    endif

    if(sys1%neu /= sys2%neu) then
       write (luout,100)'ceq_sys_diff:','different neu => ','sys1.neu =', sys1%neu,'sys2.neu =', sys2%neu
       info = -1
       if(min_diag) return
    endif

    if(sys1%ns /= sys2%ns) then
       write (luout,100)'ceq_sys_diff:','different ns => ','sys1.ns =', sys1%ns,'sys2.ns =', sys2%ns
       info = -1
       if(min_diag) return
    endif

    if(sys1%nsd /= sys2%nsd) then
       write (luout,100)'ceq_sys_diff:','different nsd => ','sys1.nsd =', sys1%nsd,'sys2.nsd =', sys2%nsd
       info = -1
       if(min_diag) return
    endif

    if(sys1%nsu /= sys2%nsu) then
       write (luout,100)'ceq_sys_diff:','different nsu => ','sys1.nsu =', sys1%nsu,'sys2.nsu =', sys2%nsu
       info = -1
       if(min_diag) return 
    endif

    if(sys1%ncs /= sys2%ncs) then
       write (luout,100)'ceq_sys_diff:','different ncs => ','sys1.ncs =', sys1%ncs,'sys2.ncs =', sys2%ncs
       info = -1
       if(min_diag) return 
    endif

    if(sys1%ng /= sys2%ng) then
       write (luout,100)'ceq_sys_diff:','different ng => ','sys1.ng =', sys1%ng,'sys2.ng =', sys2%ng
       info = -1
       if(min_diag) return 
    endif

    if(sys1%nc /= sys2%nc) then
       write (luout,100)'ceq_sys_diff:','different nc => ','sys1.nc =', sys1%nc,'sys2.nc =', sys2%nc
       info = -1
       if(min_diag) return 
    endif

    if(sys1%nrc /= sys2%nrc) then
       write (luout,100)'ceq_sys_diff:','different nrc => ','sys1.nrc =', sys1%nrc,'sys2.nrc =', sys2%nrc
       info = -1
       if(min_diag) return 
    endif

    if(sys1%nb /= sys2%nb) then
       write (luout,100)'ceq_sys_diff:','different nb => ','sys1.nb =', sys1%nb,'sys2.nb =', sys2%nb
       info = -1
       if(min_diag) return 
    endif

    if(maxval(abs(sys1%CS - sys2%CS)) > eps) then
       write (luout,*)'ceq_sys_diff: different CS'
       write (luout,300)'sys1.CS =', sys1%CS
       write (luout,300)'sys2.CS =', sys2%CS
       info = -1
       if(min_diag) return 
    endif

    if(maxval(abs(sys1%sp_order - sys2%sp_order)) > eps) then
       write (luout,*)'ceq_sys_diff: different sp_order'
       write (luout,300)'sys1.sp_order = ', sys1%sp_order
       write (luout,300)'sys2.sp_order = ', sys2%sp_order
       info = -1
       if(min_diag) return 
    endif

    if(maxval(abs(sys1%el_order - sys2%el_order)) > eps) then
       write (luout,*)'ceq_sys_diff: different el_order'
       write (luout,300)'sys1.el_order =', sys1%el_order
       write (luout,300)'sys2.el_order =', sys2%el_order
       info = -1
       if(min_diag) return 
    endif

    if(maxval(abs(sys1%Ein - sys2%Ein)) > eps) then
       write (luout,*)'ceq_sys_diff: different Ein'
       write (luout,400)'sys1.Ein =', sys1%Ein
       write (luout,400)'sys2.Ein =', sys2%Ein
       info = -1
       if(min_diag) return 
    endif

    if(maxval(abs(sys1%E - sys2%E)) > eps) then
       write (luout,*)'ceq_sys_diff: different E'
       write (luout,400)'sys1.E =', sys1%E
       write (luout,400)'sys2.E =', sys2%E
       info = -1
       if(min_diag) return 
    endif

    if(maxval(abs(sys1%Bg - sys2%Bg)) > eps) then
       write (luout,*)'ceq_sys_diff: different Bg'
       write (luout,400)'sys1.Bg =', sys1%Bg
       write (luout,400)'sys2.Bg =', sys2%Bg
       info = -1
       if(min_diag) return 
    endif

    if(maxval(abs(sys1%B - sys2%B)) > eps) then
       write (luout,*)'ceq_sys_diff: different B'
       write (luout,400)'sys1.B =', sys1%B
       write (luout,400)'sys2.B =', sys2%B
       info = -1
       if(min_diag) return 
    endif

    if(maxval(abs(sys1%BR - sys2%BR)) > eps) then
       write (luout,*)'ceq_sys_diff: different BR'
       write (luout,400)'sys1.BR =', sys1%BR
       write (luout,400)'sys2.BR =', sys2%BR
       info = -1
       if(min_diag) return 
    endif

    if(maxval(abs(sys1%A - sys2%A)) > eps) then
       write (luout,*)'ceq_sys_diff: different A'
       write (luout,400)'sys1.A =', sys1%A
       write (luout,400)'sys2.A =', sys2%A
       info = -1
       if(min_diag) return 
    endif
    
    if(maxval(abs(sys1%thermo - sys2%thermo)) > eps) then
       write (luout,*)'ceq_sys_diff: different thermo'
       write (luout,400)'sys1.thermo =', sys1%thermo
       write (luout,400)'sys2.thermo =', sys2%thermo
       info = -1
       if(min_diag) return 
    endif

    if(sys1%diag /= sys2%diag) then
       write (luout,200)'ceq_sys_diff:','diag => ','sys1.diag =', sys1%diag,'sys2.diag =', sys2%diag
       info = -1
       if(min_diag) return 
    endif
    
    if(sys1%lu /= sys2%lu) then
       write (luout,200)'ceq_sys_diff:','lu => ','sys1.lu =', sys1%lu,'sys2.lu =', sys2%lu
       info = -1
       if(min_diag) return 
    endif

    if(sys1%T_low /= sys2%T_low) then
       write (luout,200)'ceq_sys_diff:','T_low => ','sys1.T_low =', sys1%T_low,'sys2.T_low =', sys2%T_low
       info = -1
       if(min_diag) return 
    endif

    if(sys1%T_high /= sys2%T_high) then
       write (luout,200)'ceq_sys_diff:','different T_high => ','sys1.T_high =', sys1%T_high,'sys2.T_high =', sys2%T_high
       info = -1
       if(min_diag) return 
    endif

    if(sys1%frac_zm /= sys2%frac_zm) then
       write (luout,200)'ceq_sys_diff:','different frac_zm => ','sys1.frac_zm =', sys1%frac_zm,'sys2.frac_zm =', sys2%frac_zm
       info = -1
       if(min_diag) return 
    endif

    if(sys1%T_tol /= sys2%T_tol) then
       write (luout,200)'ceq_sys_diff:','different T_tol => ','sys1.T_tol =', sys1%T_tol,'sys2.T_tol =', sys2%T_tol
       info = -1
       if(min_diag) return 
    endif

    if(sys1%ds_inc /= sys2%ds_inc) then
       write (luout,200)'ceq_sys_diff:','different ds_inc => ','sys1.ds_inc =', sys1%ds_inc,'sys2.ds_inc =', sys2%ds_inc
       info = -1
       if(min_diag) return 
    endif

    if(sys1%ds_dec /= sys2%ds_dec) then
       write (luout,200)'ceq_sys_diff:','different ds_dec => ','sys1.ds_dec =', sys1%ds_dec,'sys2.ds_dec =', sys2%ds_dec
       info = -1
       if(min_diag) return 
    endif

    if(sys1%ds_min /= sys2%ds_min) then
       write (luout,200)'ceq_sys_diff:','different ds_min => ','sys1.ds_min =', sys1%ds_min,'sys2.ds_min =', sys2%ds_min
       info = -1
       if(min_diag) return 
    endif

    if(sys1%res_tol /= sys2%res_tol) then
       write (luout,200)'ceq_sys_diff:','different res_tol => ','sys1.res_tol =', sys1%res_tol,'sys2.res_tol =', sys2%res_tol
       info = -1
       if(min_diag) return 
    endif

    if(sys1%ires_fac /= sys2%ires_fac) then
       write (luout,200)'ceq_sys_diff:','different ires_fac => ','sys1.ires_fac =', sys1%ires_fac,'sys2.ires_fac =', sys2%ires_fac
       info = -1
       if(min_diag) return 
    endif

    if(sys1%logy_lim /= sys2%logy_lim) then
       write (luout,200)'ceq_sys_diff:','different logy_lim => ','sys1.logy_lim =', sys1%logy_lim,'sys2.logy_lim =', sys2%logy_lim
       info = -1
       if(min_diag) return 
    endif

    if(sys1%srat_lim /= sys2%srat_lim) then
       write (luout,200)'ceq_sys_diff:','different srat_lim => ','sys1.srat_lim =', sys1%srat_lim,'sys2.srat_lim =', sys2%srat_lim
       info = -1
       if(min_diag) return 
    endif

    if(sys1%err_huge /= sys2%err_huge) then
       write (luout,200)'ceq_sys_diff:','different err_huge => ','sys1.err_huge =', sys1%err_huge,'sys2.err_huge =', sys2%err_huge
       info = -1
       if(min_diag) return 
    endif

    if(sys1%dec_min /= sys2%dec_min) then
       write (luout,200)'ceq_sys_diff:','different dec_min => ','sys1.dec_min =', sys1%dec_min,'sys2.dec_min =', sys2%dec_min
       info = -1
       if(min_diag) return 
    endif

    if(sys1%eps_el /= sys2%eps_el) then
       write (luout,200)'ceq_sys_diff:','different eps_el => ','sys1.eps_el =', sys1%eps_el,'sys2.eps_el =', sys2%eps_el
       info = -1
       if(min_diag) return 
    endif

    if(sys1%eps_sp /= sys2%eps_sp) then
       write (luout,200)'ceq_sys_diff:','different eps_sp => ','sys1.eps_sp =', sys1%eps_sp,'sys2.eps_sp =', sys2%eps_sp
       info = -1
       if(min_diag) return 
    endif

    if(sys1%pert_tol /= sys2%pert_tol) then
       write (luout,200)'ceq_sys_diff:','different pert_tol => ','sys1.pert_tol =', sys1%pert_tol,'sys2.pert_tol =', sys2%pert_tol
       info = -1
       if(min_diag) return 
    endif

    if(sys1%pert_skip /= sys2%pert_skip) then
       write (luout,200)'ceq_sys_diff:','different pert_skip => ','sys1.pert_skip =', sys1%pert_skip,'sys2.pert_skip =', &
            sys2%pert_skip
       info = -1
       if(min_diag) return 
    endif

    if(info /= 0) write (luout,*)'ceq_sys_diff: sys1 and sys2 differ.'

100 format(1x,a,a20,1p,a10,i8,' ; ',a10,i8)
200 format(1x,a,a20,1p,a10,1p,e13.4,' ; ',a10,1p,e13.4)
300 format(1x,a,(10i8))
400 format(1x,a,1p,(6e13.4))

  end subroutine ceq_sys_diff

end module ceq_diff
