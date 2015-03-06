!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

	module ci_dat6

!  data required for modeci = 6 and 7

	use ci_prec
	use ceq_types, only: sys_type
	use cantera
	
	implicit none
	save

!--- thermochemical variables

	integer :: ns, nsin, ne, nr, nxx, mode_pdt,
     2     nsp1, nsp2, nsp3, nsp4, iorv

	real(k_dp) :: phiref, tempref, pressref, hsref, dtref,
     1	              hsfref, gascon, prc, dtc, patm, cpmref,
     2	              temp_std=298.15d0

      type(phase_t) ::  gas


	real(k_dp), allocatable   :: amolwt(:), dev(:,:), devt(:,:), 
     1	              href(:), rxnip(:), thermo_ns(:,:)

	logical :: const_pr, const_dt, user_rate,
     1		   radiation
     
      type(sys_type), pointer :: sys_uc  ! for unconstrained equilibrium

!--- numerical parameters

	integer    :: kreal, ichdas, ickcorr, m_sens, njacs, q_pade, 
     1              exp_m, lu_sens, myrank, nprocs, idproc,
     2              kreal_h, kreal_t, m_jac, m_ckwyp

	integer    :: ci_info_n(20), info_n(100)

	real(k_dp) :: tlow, thigh, tbadlo, tbadhi, temtol,
     1	              atolc, rtolc, atols, rtols, tolneg,
     2	              treal, errfac, pewarn, tewarn,
     3	              rmap1t, rmt1_op, rmap2t, rmt2_op,
	4                  sens_lim

	real(k_dp) :: ci_rinfo_n(20), rinfo_n(50)

	logical    :: ciinit_called = .false. , ciparam_called=.false., 
     1              isatab_called = .false.

!--- variables passed through to ISATAB

	integer :: ichin, ichout, ntree, isatop, noisat, idites,
     1	           ifull, kcpv, mess_pass, mpi_log

	real(k_dp) :: errtol, elpmax, elp0max, stomby, outinc

!---  other

	real(k_dp), allocatable :: xref(:)

	contains

!-----------------------------------------------------------------

	subroutine ci6_alloc

!  allocate arrays:  requires  ns, ne
!  set variables from Cantera

      implicit none

	nsp1 = ns+1
	nsp2 = ns+2
	nsp3 = ns+3
	nsp4 = ns+4

	allocate( amolwt(ns) )
	allocate( dev(ns,ne) )
	allocate( devt(ne,ns) )
	allocate( href(ns) )
	allocate( rxnip(nsp3) )
	allocate( thermo_ns(ns,15) )

	end subroutine ci6_alloc
	
!==========================================================================
      subroutine ct_thermo( ns, ctthermo, gas)

!  Return in ctthermo the thermo data for each species.
!  The data are:  switch temperature; 7 coeffs in lower range; 7 coeffs in upper range

c*****precision > double
      use cantera
      use canteraAccess
      use isat_abort_m
      implicit double precision (a-h, o-z), integer (i-n)
      type(phase_t) :: gas
      dimension ctthermo(ns,15)

      ierr = nasadata(gas%thermo_id, ctthermo)

      if (ierr .ne. 0) then
          call isat_abort('ct_thermo', 1,
     1	          mess='Failed to retrieve nasa polynoms' )
      endif

        return
	end subroutine ct_thermo

!=====================================================================


	end module ci_dat6
