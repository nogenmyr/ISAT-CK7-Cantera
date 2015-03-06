!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

	module ci_dat8

! This module contains data required for modeci = 8 and modeci = 9

!===================================================================================
!  The following is a summary of the notation used in the code for modeci=8, 9.

! FULL REPRESENTATION

! ns                    number of species (ci_dat6)
! ne                    number of elements (ci_dat6)
! z(1:ns)               specific moles of species, in order given in the Chemkin input file
! amolwt(1:ns)          molecular weights (Chemkin order, ci_dat6)
! thermo_ns(1:ns,1:15)  thhermo data (Chemkin order, ci_dat6)
! CE(1:ns,1:ne)         element matrix
! ze(1:ne)              specific moles of elements in all species (ze = CE^T*z)
!                       

! REPRESENTED AND UNREPRESENTED SPECIES

! indic_rs(1:ns)  indicator of represented species 
!                 (1->represented, 0->unrepresented)
! nrs             number of represented species:   nrs = sum( indic_rs )
! nus             number of unrepresented species: nus = ns - nrs

! zr(1:nrs)       specific moles of represented species (in same order as in z)
! zu(1:nus)       specific moles of unrepresented species (in same order as in z)

! CS(1:nrs)       species index of the represented species:   zr(irs) = z( CS(irs) )
! US(1:nus)       species index of the unrepresented species: zu(ius) = z( US(ius) )
! S2RS(1:ns)      reduced index of represented species: z(is) = zr( S2RS(is) ) for
!                 represented species; S2RS(is) = 0 for unrepresented species

! CEr(1:nrs,1:ne) element matrix for represented species
!                 =diag(indic_rs) * CE with zero rows omitted
! CEu(1:nus,1:ne)  element matrix for unrepresented species
!                 =diag(1-indic_rs) * CE with zero rows omitted

! zer(1:ne)       specific moles of elements in represented species
!                 zer = CEr^T * zr
! zeu(1:ne)       specific moles of elements in unrepresented species
!                 zeu = CEu^T * zu

! REDUCED COMPOSITIONS

! nrc             number of reduced compositions:  nrc = nrs + ne
! r(1:nrc)        reduced composition (in terms of unrepresented elements)
!                 r = {zr,zeu}
! ra(1:nrc)       alternative reduced composition (in terms of all elements)
!                 ra = {zr,ze} = r + {0, zer}

! amolwt_n(1:nrc) molecular weights corresponding to reduced species and elements
!                 r is realizable if (a) r(i) >=0, (b) dot_product( r, amolwt_n ) = 1

! BB(1:ns,1:nrc)  reduction matrix, defined such that:  r = BB^T * z
! BBT(1:nrc,1:ns) = BB^T

! BBF(1:ns,1:nrc) alternative reduction matrix, defined such that:  ra = BBF^T * z
! BBFT(1:nrc,1:ns) = BBF^T

! FACETS

! kfa             index of facet on which r(kfa)=0, r(i)>0, for i/=kfa
!                 kfa=0 denotes the interior, r(i)>0, all i

! DENSITY AND TEMPERATURE APPROXIMATION

! dta_lu     log file for density temperature approximations
! clipt_log  enable/disable temperature clipping log messages [default 1]

! Parameters used for density and temperature approximation in ci_dpt_dr module.
! The following parameters can be set in the ci_ext.nml namelist file.

! autou      binary variable for auto-updating dens. and temp values.
!            To enable auto-update, set autou = 1
!             [default value, autou = 0 (do not auto-update)]
! ufreq      frequency with which to attempt updating the temp. values
!             [default value, ufreq =  50 ]
! tatol      mean relative error tolerance in the approx. temp. Ta
!             [default value, tatol = 1.d-2 (1% error)]
! clipt      binary, clip temperature values between [tbadlo, tbadhi]
!            To enable set clipt = 1 [default value, clipt = 0]
! dtlog      write dens. temp. approximations attempts to dta.log file
!            To enable dtlog = 1 [default values. dtlog = 1]


!====================================================================================
!         nrs:   number of represented species
!         rname: name of the represented species
!         nus:   number of unrepresented species
!         nrc:   nrs+ne, number of represented compositions
!         CE :   element matrix of all the species
!         CEr:   element matrix of the represented species,  dimension   [nrs, ne]   
!         CEu:   element matrix of the unrepresented species, dimension [nus, ne] 
!         BB, dimension  ns*nrc:  =[CS CE_tmp], where in CE_tmp the corresponding 
!                                   rows for represented species are zero 
!         BBT:   transpose of BB
!         BBF:   =[CS  CE],  dimension  ns*nrc:
!         BBFT:  transpose of BBF 
!         CS:    species index of the represented species, dimension nrs: r(irs) = z( CS(irs) )
!         S2RS:  reduced index of represented species, dimension ns: z(is) = r( S2RS(is) ) for
!                represented species; S2RS(is) = 0 for unrepresented species
!         indic_rs:  indicator for represented species, 
!                    =1 for represented species, =0 for unrepresented s[ecies
!         i_first_us: index of first unrepresented species
!
!         [Rspace, Rtemp]= qr(BBF): Rspace and Rtemp are matrices obtained from QR 
!                                   factorization of matrix BBF 
!
!         href_n:  dimension nrc: nominal enthalpy of formation  
!                  for represented species,  enthalpy of formation is taken from thermal data

!         amolwt_n: dimension nrc:  nominal molecular weight
!                   for represented species,  taken  to be its molecular weight 
!                   for elements in unrepresented species, taken to be the atomic weight 
!                   amolwt_n is needed when transforming the species specific moles to mass fractions, 
!                   see subroutine cicmp8 

!         mol_n, dimension nrc:   nominal mole parameter
!                   for represented species equal to 1 
!                   for elements in unrepresented species, a nominal value is given, see 
!                   subroutine ci_init8 for more details.
!                   mol_n is needed when transorming the species specific moles to mole fraction
!                   (see subroutine cicmp8) as well as to evaluate the density after ISAT retrieve, 
!                   see subroutine  cirxn8

	use ci_prec
	use ci_dat6
	implicit none
	save

!--- thermochemical variables

	integer :: nrs, nus, nrc      ! number of represented species
	integer :: i_first_us         ! index of first unrepresented species

	integer, parameter    :: k_href_n = 0 ! nominal enthalpy computation mode (set_href_n)

	real(k_dp), parameter ::  epsilon_rank=1.d-9, epsilon=1.d-8,
     1                          epsilon_b=1.d-8, dec_fac =0.95
     
	real(k_dp), parameter :: X_tiny = 1.d-200 ! smallest mole fraction
	real(k_dp), parameter :: z_tiny = 1.d-30 ! smallest positive z
	real(k_dp), parameter :: tol_lp = 1.d-9	! tolerance used in LP
      
 ! CE point is taken as ICE point if any reduced specific mole is less than tol_rceb
	real(k_dp), parameter :: tol_rceb = 1.d-10
      
 ! At boundary, species taken to be zero if specific mole is less than tol_zlinb1
	real(k_dp), parameter :: tol_zlinb1 = 1.d-18
      
 ! At boundary, if other tests fail to identify zero species, take smallest species
 !   to be zero if it is less than tol_zlinb2
	real(k_dp), parameter :: tol_zlinb2 = 1.d-15
      
 ! CE point is taken as ICE point if |z_CE - z_EQ| < tol_eq
	real(k_dp), parameter :: tol_eq = 1.d-10
      
 ! Component of vector kT is deemed to be zero if less than tol_kT
	real(k_dp), parameter :: tol_kT = 1.d-8
      
 !  Criteria for extrapolating to the boundary
	real(k_dp), parameter :: tol_extra_1 = 1.d-4
	real(k_dp), parameter :: tol_extra_r = 5.d-2
	real(k_dp), parameter :: tol_extra_t = 5d0
	real(k_dp), parameter :: tol_extra_a = 1d-30
       
 ! tolerance on CEQ mole fractions
!XXX    real(k_dp), parameter :: tol_dX = huge(1.d0) 
	real(k_dp), parameter :: tol_dX = 1.e-8 ! XXXXXXXXXXXX
     
!XXX    integer,parameter :: npoint_max = 10,  iter_max=8
	integer,parameter :: npoint_max = 10,  iter_max=20
	
! dens. and temp. approximations
	integer    :: dta_lu          ! lu for dta log file
	integer    :: autou  = 0      ! by default do not auto-update
	integer    :: ufreq  = 50     ! check for update every 50th species recon.
	real(k_dp) :: tatol  = 1.d-2  ! error tolerance for Ta
	integer    :: clipt = 0       ! by default do not clip temperature values
	integer    :: dtlog = 1	      ! by default enable desn. temp. approx. log
	integer    :: clipt_log = 1   ! enable/disable temp clipping log messages

	real(k_dp), allocatable :: CE(:,:), CEr(:,:), CEu(:,:),
     1	        BB(:,:), BBT(:,:), mol_n(:), amolwt_n(:), href_n(:),
     2          BBF(:,:), BBFT(:,:), Rspace(:,:), Rtemp(:,:),
     3          r_unit(:)
     
	character(16), allocatable :: rname(:)
     	integer, allocatable       :: CS(:), US(:), S2RS(:), indic_rs(:)

! stream info
     	real(k_dp), allocatable :: ccstrm_f(:,:) ! stream info
	character(16), allocatable :: cmpsym_f(:) ! species names
     	
!  Quantities set in subroutine determined
	integer                 :: nzud, nrud
	integer, allocatable    :: i_zud(:), i_rud(:), k_zdet(:), 
     1                           urud_to_ur(:)
	real(k_dp), allocatable :: r2zdet(:,:)
     	
     	real(k_dp) :: dice(100) ! ICE diagnostic information

! dice( 1) - call number
! dice( 2) - outcome, ICE_stats(1)
! dice( 3) - T_CE, constrained equil temp
! dice( 4) - r_CE_min
! dice( 5) - dz_eq, z-distance from CE to full equilibrium
! dice( 6) - dzdt_CE, |dzdt| at CE
! dice( 7) - npoint, number of pre-image points
! dice( 8) - SVA, attenuation factor in first unrepresented composition
! dice( 9) - |dr|, r-residual
! dice(10) - |dg_pre|, g-residual
! dice(11) - cos_last, cos angle between last 2 PIC segments
! dice(12) - r_pre_min, min(r) at pre-image point
! dice( 1) - call number
! dice( 1) - call number
! dice( 1) - call number
! dice( 1) - call number

	contains

!-----------------------------------------------------------------

	subroutine ci8_alloc

!  allocate arrays:  requires  ns, ne
	nrc = nrs+ne
	nus = ns-nrs
	
	allocate( rname(nrs) ) 
	allocate( CS(nrs) )
	allocate( US(nus) )
	allocate( S2RS(ns) )
	allocate( indic_rs(ns) )
	allocate( mol_n(nrc) )
	allocate( amolwt_n(nrc) )
	allocate( href_n(nrc) )
	allocate( CE(ns, ne) )
	allocate( CEr(nrs, ne))
	allocate( CEu(nus, ne))
	allocate( BB(ns, nrc) )
	allocate( BBT(nrc, ns) )
	allocate( BBF(ns, nrc) )
	allocate( BBFT(nrc, ns) )
	allocate( Rspace(ns, nrc) )
	allocate( Rtemp(nrc,nrc)) ![Rspace, Rtemp]= qr(BBF)

	end subroutine ci8_alloc

!-----------------------------------------------------------------

	subroutine ci8_param_init

	use ci_dat
	use ci_dat6
	implicit none
	
	logical :: exists
	integer :: lu
        character(30) :: blank, head, tail, name
	namelist / ci_extnml / autou, tatol, ufreq, clipt, dtlog
	
!-----------  read changes from namelist  -----------------------------

	blank = repeat(' ',30)!  generate file name:  ci_ext_P.nml
	head  = blank
	head  = 'ci_ext'
	tail  = blank   
	tail  = 'nml'
	call isat_file_name( head, -1, idproc, tail, name )

	inquire( file = name, exist = exists )

	if(op_rank) write(luout,*) ' '
	if( exists ) then
	   call isat_lu( lu )
	   open( lu, file = name, position = 'rewind',
     1	         action = 'read' )
	   read( lu, nml = ci_extnml )
	   if(op_rank) write(luout,*) 'Variables have been read from ', 
     1                                 name
	else
	   if(op_rank) write(luout,*) 'No variables have been read from ', 
     1                                name
	endif

	if(op_rank) write( luout,* ) ' '
	if(op_rank) write( luout, nml = ci_extnml )

! open dta log file
	blank = repeat(' ',30)!  generate file name:  ci_ext_P.nml
	head  = blank
	head  = 'dta'
	tail  = blank   
	tail  = 'log'
	call isat_file_name( head, -1, myrank, tail, name )

	call isat_lu(dta_lu)
	open( dta_lu, file = name, status='replace', action='write' )
	
	end subroutine ci8_param_init
!-----------------------------------------------------------------
	end module ci_dat8
