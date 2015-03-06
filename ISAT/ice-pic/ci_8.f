!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

	module ci_8

!  This module, contains all the ci interfaces for ISAT-DR:  
!    modeci = 8 (ISAT/RCCE) and 
!    modeci = 9 (ISAT/ICE-PIC).
!
!  It consists of the following subroutines:
!  cinit8: CI initialization routine for modeci=8 and modeci=9
!
!  cirxn8: CI routine to return composition c(t) resulting from
!       reaction for a time t from the initial composition c(0) for modeci=8 and modeci=9.

!  cicmp8: CI routine to  return composition  given compact representation for modeci=8 and modeci=9.
!----------------------------------------------------------------------------------------------------
!  DEVELOPMENT NOTE: replace stops by calls to isat_abort
!----------------------------------------------------------------------------------------------------

	use ci_dat
	use ci_dat8
	use ci_ck
	use isat_val, only: isat_lic, coderet
	use ci_cksubs
	use ci_utils
	use streams_mod
	use ci_dpt_dr
	contains

!=========================================================================

        subroutine cinit8

!  chemistry interface initialization routine for  modeci = 8 and modeci = 9
!   which use Chemkin.

!  files:
!       streams.in - input file of stream information.
!BEGEXTRACT
!
!      Specification of the file  streams.in  for modeci = 8 and 9
!
!         1st record - modeci (8 or 9)
!         2nd record - nstr, nsin  
!                      nstr - number of streams (integer)
!                      nsin - number of non-zero species (integer)
!         subsequent - k, p, T, c(1), c(2),...,c(nsin)
!                      k = 1 if stream is as stated (integer)
!                      k = 2 if stream is an equilibrium mixture with the
!                            same elemental composition, pressure and
!                            enthalpy as that stated.
!                      p    - pressure in atm  (real)
!                      T    - temperature in K (real)
!                      c(i) - composition in relative volume/mole units (real)
!         subsequent  - nrs,  number of represented species
!         subsequent  - name of the represented species, one line for each species.

!ENDEXTRACT

!  comment: 
!        c(i) does not need to be normalized
!       
	implicit none

	integer :: retcode

	if(op_rank) write(luout,*)' '

	if( myrank == 0 ) then
	   call isat_lic( luout, retcode )
	   if( retcode /= coderet ) call isat_abort('cinit8',0, 
     1	       mess='Invalid license code' )
	endif

!----------------  set defaults, and read changes from ci.nml

!  set all parameters to zero if they have not been set in ciparam
	if( .not.ciparam_called ) then
	   ci_info_n  = 0
	   ci_rinfo_n = 0.d0
	   info_n     = 0
	   rinfo_n    = 0.d0
	endif

! SBP XXX why are these turned off?  
! The following is particular to DR, no accuracy check 
      if ( modeci == 8 .or. modeci == 9 ) then  !  SBP XXX would be be here otherwise?
         info_n(28) = 0       ! idites
         rinfo_n(21) = -2.d0  ! ecc_r=freq
      endif

!This is particular to ICE-PIC, the reconstructed point is not exactly the given point x  
   
      if ( modeci == 9 ) info_n(47) = 1  !  de_nearby
               
	call ci_info_set
	call ci8_param_init

!--------------  read number of streams and non-zero species  ----------------

	if( luin > 0 ) then
	   if( strm_format == 1 ) then
	      read(luin,*,err=140,end=142) nstr, nsin
	      if( nstr < 0 ) call isat_abort( 'cinit8',3,
     1	                  mess='bad value of  nstr', isv = nstr )
	      if(op_rank) write(luout,*)'number of streams, nstr = ', nstr
	   else
	      nstr = StrmDat%nStream
	   endif
     	else
	   nstr = 0
	   nsin = 0
	endif
     	   
!---------------  perform all other initialization

	call ci_init8

	return

!-------------  error conditions

140	call isat_abort('cinit8',1, mess='error reading nstr and nsin' )
142	call isat_abort('cinit8',2, mess='hit end reading nstr and nsin' )

	end subroutine cinit8

!=========================================================================

	subroutine ci_init8

!  routine to get basic data from Chemkin and to initialize for modeci = 8 and modeci=9

!  NOTE:  reference values are obtained from: ci.nml; or if not set
!         there, from streams; or if not from there, set here.

	use ceq_state_m
	use ci_ice_cksubs, only: set_href_n

	implicit none

	logical       :: kerr
	integer       :: ev(ns,ne), rv(ns,nr), kstr(nstr),
     1	                kkindex, kpress, ist, is, iii,jjj,kkk
	real(k_dp)    :: phiavg(ns), cpr, pstr(nstr),
     1	               tstr(nstr), ystr(ns,nstr), hstr(nstr), 
     2	               rhostr(nstr),  tempcp, ewt(ne), hs_n,
     3                 phi_eq(ns), evtemp(ns, ne)
   
	integer    ::  i, j, qr_info, info
	real(k_dp) :: qr_work(ns+10)
	real(k_dp), allocatable :: qr_tau(:), BBFtemp(:,:)

! Initialized to full representation, then copy to reduced representation
	integer       :: nc_f, nfull_f

	real(k_dp), allocatable :: dptstr_f(:,:)

	nc_f    = nsp1
	nfull_f = nsp4
  
	allocate( ccstrm_f(nc_f, nstr) )
	allocate( dptstr_f(3, nstr)    )
	allocate( cmpsym_f(nfull_f)  )

	call ctnu( ns,   rv, gas )  
	ev = dev
	evtemp = dev

!-----  define names of full composition variables

	cmpsym_f(1:ns) = snames(1:ns)
	cmpsym_f(ns+1) = symb_dens
	cmpsym_f(ns+2) = symb_temp
	cmpsym_f(ns+3) = symb_press
	cmpsym_f(ns+4) = symb_enth

!  print out species info

	if(op_rank) call op_species8

	if( nstr == 0 ) go to 700

!----- treatment with streams -------------
	if ( strm_format == 1 ) then
	   if( nsin < 1  .or.  nsin > ns ) call isat_abort( 'ci_init8', 1,
     1	     mess = 'bad value of nsin = ', isv = nsin )

	   ccstrm_f = 0.d0

	   do ist = 1, nstr
	      read(luin,*,err=144,end=146)kstr(ist),pstr(ist),tstr(ist),
     1	       (ccstrm_f(is,ist),is=1,nsin)
     	   end do
	else
	   kstr = StrmDat%kstrm
	   pstr = StrmDat%Cvalue(ns+1,:)
	   tstr = StrmDat%Cvalue(ns+2,:)
	   do ist = 1,nstr
	      if(StrmDat%kmole(ist) == 0) then
		 call y2phi( StrmDat%Cvalue(1:ns,ist), amolwt, ns, 
     1                  ccstrm_f(1:ns,ist))
	      else
		 ccstrm_f(1:ns,ist) = StrmDat%Cvalue(1:ns,ist)
	      endif
	   enddo
	endif


! read in the represented species
	if ( strm_format == 1 ) then
	   read(luin,*,err=148,end=150) nrs
       
	   if ( nrs<=0 .or. nrs>=ns-ne )  call isat_abort( 'ci_init8', 1,
     1	     mess = 'bad value of nrs = ', isv = nrs )
	   
	else
	   nrs = StrmDat%nSp
	endif

	nc   = nrs+ne+1
	nfull= nc+3

! allocate the necessary arrays needed in modeci=8 and modeci=9       
       	call ci_alloc   ! allocate ccstrm, dptstr, cflstr and cmpsym
       	call ci8_alloc  ! allocate rname, CS, mol_n, amolwt_n, CE, CEr, CEu, BB 
       	CS=0
       	US=0
       	S2RS=0
       	indic_rs = 0
       	mol_n=0
       	amolwt_n=0
       	href_n=0     
       	CE=0
       	CEr=0
       	CEu=0
       	BB=0
       	BBT=0
       	BBF=0
       	BBFT=0
       	Rspace=0
       	Rtemp=0
      
       	
! read in the represented species name and determine its corresponding index
	if ( strm_format == 1 ) then
	   do  iii=1, nrs
	      read(luin,*,err=152,end=154) rname(iii)
	   enddo
	else
	   do  iii=1, nrs
	      rname(iii) = StrmDat%Spnames(iii)
	   enddo
	endif
         
        do  iii=1, nrs
            call ctcomp(rname(iii), kkindex, gas)
            CS(iii)=kkindex
            indic_rs(kkindex) = 1 ! indicator of represented species
        enddo

	! sort the represented species
	call sort(CS)
	do iii=1,nrs
	   S2RS(CS(iii)) = iii
	enddo
	
	j = 0
	do i = 1, ns
           if( indic_rs(i) == 0 ) then
              j = j + 1
              US(j) = i  !  j-th unrepresented species is species US(j)
           endif	
	end do

        do i = 1, ns
           if( indic_rs(i) == 0 ) then
              i_first_us = i
              exit
           endif
        end do
        
        
! form the CE, CEr, CEu, and BB
	CE(1:ns, 1:ne)=dev(1:ns,1:ne)
         
         do iii=1,nrs
           CEr(iii, 1:ne) =dev(CS(iii), 1:ne)
           evtemp(CS(iii), 1:ne)=0 ! obtain the element vector for the unrepresented speceis
           BB(CS(iii), iii)=1
           BBF(CS(iii), iii)=1
           cmpsym(iii)=cmpsym_f(CS(iii))
           amolwt_n(iii)=amolwt(CS(iii))
           mol_n(iii)=1
           href_n(iii)=href(CS(iii))       ! The elmenents (norminal species) enthalpy of formation is set to be zero
         enddo
! VH: prepend element names with symb_elem='E', eg. EO to distinguish from species O
	 do iii=nrs+1,nrs+ne
	    cmpsym(iii) = symb_elem
	    cmpsym(iii)(2:16) = enames(iii-nrs)(1:15)
	 enddo

         cmpsym(nrs+ne+1) = symb_dens
	 cmpsym(nrs+ne+2) = symb_temp
	 cmpsym(nrs+ne+3) = symb_press
	 cmpsym(nrs+ne+4) = symb_enth

!  assign mole_n, amolwt_n  
         call ctawt(  ewt, gas)
         amolwt_n(nrs+1:nrs+ne)=ewt(1:ne)
           
          BB(1:ns, nrs+1: nrs+ne)= evtemp(1:ns, 1:ne)
          BBF(1:ns, nrs+1: nrs+ne)= CE(1:ns, 1:ne)
       	  BBT=transpose(BB) 
       	  BBFT=transpose(BBF) 
       	  jjj=0;
       	  do iii=1,ns
       	    if (maxval(evtemp(iii,1:ne))>0) then
       	      jjj=jjj+1 
       	      CEu(jjj,1:ne)=evtemp(iii, 1:ne)
            endif
       	  enddo
       	  
 !   form the represented space and unrepresented space
        allocate ( qr_tau(nrc) )
       	allocate ( BBFtemp(ns, nrc) )
       	    	    
        BBFtemp(1:ns,1:nrc)=BBF(1:ns,1:nrc)
        qr_tau=0
        qr_work=0
        qr_info=0
        
        call dgeqrf(ns,nrc,BBFtemp,ns,qr_tau,qr_work,nrc, qr_info )  ! [Q,R] = qr([wt pt])
       
        if(qr_info /= 0 ) then
         write(0,*) 'ci_8: dgeqrf failed, qr_info = ', qr_info
         stop
        endif
        
       Rtemp=0
       do i=1, nrc
        Rtemp(i, i:nrc)=BBFtemp(i,i:nrc)
       enddo
        
      call dorgqr(ns,nrc,nrc,BBFtemp,ns,qr_tau,qr_work,nrc,qr_info)
      if(qr_info /= 0 ) then
       write(0,*) 'ci_8: dorgqr failed, qr_info = ', qr_info
       stop
      endif
        
      Rspace=BBFtemp(1:ns,1: nrc)       
! check whether the representation matrix is singular or not.
       if(abs(Rtemp(nrc,nrc)) <1.d-6) then
        write(0,*) 'ci_8, the representation matrix in reduced 
     1              representation is singular'
        write(0,*) (Rtemp(i,i), i=1, nrc)
        write(0,*) 'Bad choice of represented species'
!DR-08 This part need to be further improved in the future. 
!It should provide more informative diagonostic information.         
        stop  ! SBP XXX should call isat_abort
       endif 

! Initialize the ci_dpt_dr module 
! CEu must be set before this call
!-----------------------------
	call ci_dpt_init
!-----------------------------

!  the nominal mole for elements is tricky, for inert species, equal to its corresponding element vector component       	  
       	 
       	 do iii=1, ne
       	    jjj=0
       	    kkk=0
       	   do is=1, nus 
	       if(CEu(is, iii)>0) then
	        kkk=kkk+CEu(is, iii)
	        jjj=jjj+1
	       endif
	     enddo
	     if(jjj==1) then
	      mol_n(iii+nrs)=kkk*1.
	     else
	      mol_n(iii+nrs)= kkk*1./(jjj*ne*0.5)
	     endif 
	   enddo   
     	  
!-----  derived stream properties  -----------------------------------------

!  normalize to form specific mole numbers, enthalpy and density
! (note that hstr is enthalpy, not sensible enthalpy)

	kpress = 0
	do ist    = 1, nstr

	   pstr(ist)     = pstr(ist) * patm
	   if( pstr(ist) .ne. pstr(1) ) kpress = max( 1, kpress )

	   call phi2y( ccstrm_f(1:ns,ist), amolwt, ns, ystr(:,ist) )
	   call y2phi( ystr(:,ist), amolwt, ns, ccstrm_f(1:ns,ist) )
	   call cthbms( tstr(ist), ystr(:,ist),   hstr(ist), gas )

!  determine equilibrium if required: kstr=2 for fixed (p,h); kstr>=3 for fixed (p,T)

	   if( kstr(ist) >= 2 ) then
	      if(op_rank) then
	      write(luout,*)' '
	      write(luout,'(a,a,i3,a,i2)')
     1                      'Starting equilibrium calculation ',
     1	                    'for stream',ist, ', kstr = ', kstr(ist)
	      endif

	      if( kstr(ist) == 2 ) then
		 call ceq_state( sys_uc, N=ccstrm_f(1:ns,ist),
     1                      p_cgs=pstr(ist), N_h=ccstrm_f(1:ns,ist),
     2                      T_h=tstr(ist), N_eq=phi_eq, T_eq=tstr(ist),
     3                      info=info )
	      else
		 call ceq_state( sys_uc, N=ccstrm_f(1:ns,ist),
     1                      p_cgs=pstr(ist), T=tstr(ist),
     2                      N_eq=phi_eq, info=info )
	      endif

	      if(op_rank) then
	      write(luout,*)' '
	      write(luout,'(a,f8.2,a)')'Equilibrium temperature ='
     1                               ,tstr(ist), ' (K)'
	      endif
	      call phi2y( phi_eq, amolwt, ns, ystr(:,ist) )
	      call y2phi( ystr(:,ist), amolwt, ns, ccstrm_f(1:ns,ist) )
	      call cthbms( tstr(ist), ystr(:,ist),hstr(ist),gas)
	   endif

	   call ctrhoy( pstr(ist), tstr(ist), ystr(:,ist),  
     1		rhostr(ist), gas )

!  set the sensible enthalpy contained in the represented species hs_n= h- href_n * c
	   ccstrm(1:nrs+ne,ist)=matmul(BBT, ccstrm_f(1:ns, ist)) 
	  
      end do
      
      call set_href_n( k_href_n ) ! set the nominal enthalpies of formation
   
      do ist    = 1, nstr

	   call h2hs_n(hstr(ist),href_n,ccstrm(1:nrs+ne,ist), 
     1               nrs+ne, hs_n )
         
         ccstrm_f(nc_f,ist)= hs_n
         ccstrm(nc,ist)= ccstrm_f(nc_f,ist)
         
         dptstr_f(1,ist) = rhostr(ist)
	   dptstr_f(2,ist) = pstr(ist)
	   dptstr_f(3,ist) = tstr(ist)
	   dptstr(1:3,ist)=dptstr_f(1:3,ist)
	   
	end do
	
!  obtain  tempref and pressref from streams

	if( tempref  == 0.d0 ) tempref  = maxval( tstr )
	if( pressref == 0.d0 ) pressref = maxval( pstr )

!  check assumption of constant pressure

	if( const_pr  .and.  kpress /= 0 ) then
	   if(op_rank) then
	      write(luout,*)' '
	      write(luout,*)'===WARNING===: const_pr = .true., but',
     1	                 ' streams have different pressures.'
	      write(luout,*)'Using variable pressure'
	   endif
	   const_pr = .false.
	endif
      
       
!  print stream info
	if(op_rank) call op_streams8

!---------  end of streams considerations

700	continue

!----------- define reference values
!  For isatab:  x = phi(0), hs(0), press,    dt,   and these are scaled by:
!                   phiref, hsref, pressref, dtref
!               f = phi(dt), hs(dt), temp(dt),   which are scaled by
!                   phiref,  hsfref, tempref
!  The reference values can be set in ci.nml

	if( phiref   == 0.d0 ) phiref   = 0.1d0
	if( tempref  == 0.d0 ) tempref  = 1000.d0
	if( pressref == 0.d0 ) pressref = patm

!  define reference specific heat at  tempref   with
!      mole fraction = 1/ns for all species -------

	phiavg = 1.d0/float(ns)

	tempcp = min( tempref, tbadhi )
	call ctcpbs( tempcp, phiavg,   cpr, gas)
	if( hsref == 0.d0 ) hsref = cpr * tempref

	if(op_rank) then
	   write(luout,*)' '
	   write(luout,*)'Reference values: '
	   write(luout,600) phiref
	   write(luout,601) tempref
	   write(luout,602) pressref
	   write(luout,603) hsref
	endif

	if( hsfref == 0.d0 ) then
	   hsfref = hsref
	else
	   if(op_rank) write(luout,604) hsfref
	endif

	if( dtref /= 0.d0 ) then
	   if(op_rank) write(luout,605) dtref
	endif
	   
600	format('Species:     phiref   = ', 1p,e12.4)
601	format('Temperature: tempref  = ', 1p,e12.4)
602	format('Pressure:    pressref = ', 1p,e12.4)
603	format('enth.: hsref    = ', 1p,e12.4)
604	format('enth.: hsfref   = ', 1p,e12.4)
605	format('Time interv: dtref    = ', 1p,e12.4)

!---------------------  set mode_pdt and nxx

	if(op_rank) write(luout,*)' '

	if( .not.const_pr ) then
	   if(op_rank) write(luout,*)'Pressure is taken to be variable.'
	   prc = -huge(1.d0)
	   if( .not.const_dt ) then
	      mode_pdt = 1
	      nxx      = nc +2 
	   else
	      mode_pdt = 2
	      nxx      = nc+1
	   endif
	else
	   if( nstr > 0 ) then
	      prc = pstr(1)
	   else
	      prc = pressref
	   endif 

	   if( .not.const_dt ) then
	      mode_pdt = 3
	      nxx      = nc+1
	   else
	      mode_pdt = 4
	      nxx      = nc
	   endif
	   if(op_rank) write(luout,710) prc
710	   format('Pressure is taken to be constant, = ', 1p,e13.4)
	endif

	if(op_rank) write(luout,*)' '

	if( const_dt ) then
	   if(op_rank) write(luout,*)'Time step is taken to be constant.'
	else
	   if(op_rank) write(luout,*)'Time step is taken to be variable.'
	endif

	if(op_rank) write(luout,*)' '

	if( radiation ) then
	   if(op_rank) write(luout,*)'Enthalpy equation includes ',
     1	                 'radiative heat loss.'
	else
	   if(op_rank) write(luout,*)'Enthalpy equation does not ',
     1	                 'include radiative heat loss.'
	endif

	if(op_rank) write(luout,*)' '
	if(op_rank) write(luout,*)
     1	  'Number of independent ISAT variables, nx = ', nxx


	if(op_rank) call isat_flush(luout)
	
	if (allocated(qr_tau))    deallocate(qr_tau)
        if (allocated(BBFtemp))   deallocate(BBFtemp)
        if (allocated(dptstr_f )) deallocate(dptstr_f)

	return

!---------  error conditions

140	call isat_abort(  'ci_init8', 2, mess = 'error reading nstr' )
142	call isat_abort(  'ci_init8', 3, mess = 'hit end reading nstr' )
144	call isat_abort(  'ci_init8', 4, mess =
     1	   ' error reading stream data. ist=', isv = ist )
146	call isat_abort(  'ci_init8', 5, mess =
     1	   ' hit end reading stream data. ist=', isv = ist )
148	call isat_abort(  'ci_init8', 6, mess =' error reading nrs' )
150	call isat_abort(  'ci_init8', 7, mess =' hit end reading nrs' )
152	call isat_abort(  'ci_init8', 8, mess =
     1     ' error reading represented species' )
154	call isat_abort(  'ci_init8', 9, mess =' hit end reading rname' )
	return

!------------  internal subroutines  ------------------------------------

	contains
	
	subroutine op_species8 !---  print information about species

	integer :: is, ie, ir

	write(luout,*)' '
	write(luout,*)'Species, mol. wt., and their elemental composition'
	write(luout,*)' '
	do is = 1, ns
	   write(luout,600)is,snames(is),amolwt(is),(ev(is,ie),ie=1,ne)
	end do
	call isat_flush(luout)

	if( iorv == 1 ) then
	   write(luout,*)' '
	   write(luout,*)'Reaction vectors (transposed) '
	   write(luout,*)'(set  iorv=0  in  ci.nml  to suppress ',
     1	                 'this output.)'
	   write(luout,*)' '
	   write(luout,610)(is, is=0, ns)
	   do ir = 1, nr
	      write(luout,610)ir,(rv(is,ir),is=1,ns)
	   end do

	   write(luout,*)' '
	   call isat_flush(luout)
600	   format(1x,i4,2x,a16,f8.4,20i3)
610	   format(1x,25i3,/,(4x,24i3))
	endif

	end subroutine op_species8 

	subroutine op_streams8 !---  print information about streams

	integer :: is

	write(luout,*)' '
	write(luout,*)'Output from ci_init8: stream information ',
     1	              '(Chemkin units)'
	write(luout,*)' '
	write(luout,600)(pstr(ist),ist=1,nstr)
	write(luout,610)(tstr(ist),ist=1,nstr)
	write(luout,620)(rhostr(ist),ist=1,nstr)
	write(luout,630)(hstr(ist),ist=1,nstr)
	write(luout,635)(ccstrm_f(nc_f,ist),ist=1,nstr)
	write(luout,*)' '
	write(luout,*)'Specific mole numbers'
	write(luout,*)' '
	do is = 1, ns
	   write(luout,640)is,snames(is),(ccstrm_f(is,ist),ist=1,nstr)
	end do
	write(luout,*)' '
	write(luout,*)'Mass fractions '
	write(luout,*)' '
	do is = 1, ns
	   write(luout,640)is,snames(is),(ystr(is,ist),ist=1,nstr)
	end do
	write(luout,*)' '
	write(luout,*)'Reduced representation by DR '
	write(luout,*)' '
	write(luout,*)'Specific mole numbers and enthalpy'
	write(luout,*)' '
	do is = 1, nrs+ne
	   write(luout,650)is,cmpsym(is),(ccstrm(is,ist),ist=1,nstr)
	end do	
         write(luout,650)is,cmpsym(nc+3),(ccstrm(nc,ist),ist=1,nstr)
600	format('Press.  ',1p,20e12.4)
610	format('Temp    ',1p,20e12.4)
620	format('Dens.   ',1p,20e12.4)
630	format('Enth    ',1p,20e12.4)
635	format('Enth',1p,20e12.4)
640	format(i4,2x,a16,2x,     1p,20e12.4)
650     format(i4,2x,a16,2x,     1p,20e12.4)
	end subroutine op_streams8

	end subroutine ci_init8
!=========================================================================

	subroutine cicmp8( cc, mode, comp )

!  given compact representation, return composition:

!  input:
!	cc    - compact representation, size (nc = nrs+ne+1)
!	mode  = 1 - express species as mole fractions
!	      = 2 - express species as mass fractions
!	      = 3 - express species as specific mole numbers
!	comp(1) - pressure in Chemkin units
!  output:
!	comp - composition consisting of:
!	       ( size = nrs+ne+4; nrc = nrs+ne ) 
!	comp(1:nrc) = species concentrations
!	comp(nrc+1) = density (g/cc)
!	comp(nrc+2) = temperature (K)
!	comp(nrc+3) = pressure (atm)
!	comp(nrc+4) = enthalpy (ergs/mol)

	implicit none

	real(k_dp), intent(in)    :: cc(nc)
	real(k_dp), intent(inout) :: comp(nfull)
	integer, intent(in)       :: mode
	integer :: ii
	real(k_dp) :: hs_n, h, y(nrc), x(nrc), dpt(3)

! set species and make conversions
	if( mode == 1 ) then
	   do ii=1,nrc
	      x(ii)=cc(ii)/mol_n(ii)
	   enddo
	   comp(1:nrc) = x(1:nrc) / sum( x(1:nrc) )
	   
	elseif( mode == 2 ) then
	   do ii=1,nrc
	      y(ii)=cc(ii)*amolwt_n(ii)
	   enddo
	   comp(1:nrc) = y(1:nrc) / sum( y(1:nrc) )
	   
	elseif( mode == 3 ) then
	   comp(1:nrc) = cc(1:nrc)
	endif

! set dpt returned from cidpt_dr
	dpt(2) = comp(1)
	call cidpt_dr( cc, dpt )

	comp(nrc+1) = dpt(1)
	comp(nrc+2) = dpt(3)
	comp(nrc+3) = dpt(2)
	
! set enthalpy
	hs_n       = cc(nc) 
	call hs2h_n( hs_n, href_n, cc(1:nrc), nrc, h )  ! h = enthalpy
	comp(nrc+4) = h

	return
	end subroutine cicmp8


!=========================================================================

	subroutine cirxn8( t, c0, ct, dpt )

!  chemistry interface routine to return composition c(t) resulting from
!       reaction for a time t from the initial composition c(0).  Also
!       returned are density, pressure and temperature.
!       This version for DR/ISAT

!  input:
!       t     - time, duration of reaction (real)
!       c0    - initial composition vector (r, h) (real)
!               r: the specific mole number of  nrs represented species 
!                  + specific mole number of ne element in the unrepresented species
!	dpt(2) - pressure
!  output:
!       ct    - final composition vector (real)
!       dpt   - density, pressure and temperature (real)

	implicit none

	real(k_dp), intent(in)    :: t, c0(nc)
	real(k_dp), intent(inout) :: dpt(3)
	real(k_dp), intent(out)   :: ct(nc)

	real(k_dp), allocatable, save :: rinfo(:), rusr(:)
	integer, save    :: ifst=0, idtab, info(100), nx, nf, nh, mode

	real(k_dp) :: x(nxx), f(nc+1), press, dfdx(nc+1,nxx), hvar(1), 
     1              stats(100), temp_a, ftemp(nrs+ne), sz

	integer    :: iusr(1) ,iii
	
	external cirmap_dr_rcce, cirmap_dr_ice_pic

!--------------  start of execution  -------------------------------

	press = dpt(2)

	if( ifst /= 0 ) go to 60
	ifst = 1

!===============  initialization on first call =======================

!----------  check press and dt,  set dtc  -----------------------

	if( press <= 0.d0 ) call isat_abort('cirxn8', 1,
     1	   mess='bad pressure', rsv = press )

	if( t <= 0.d0 ) call isat_abort('cirxn8', 2,
     1	   mess='bad dt', rsv = press )

     	if( const_dt ) then
	   dtc = t
	else
	   dtc = -1.d0
	endif

	if( dtref <= 0.d0 ) dtref = t   !  SBP  8/20/04

!----------  check nxx set correctly

	if( mode_pdt == 1 ) then
	   if( nxx /= nc+2 ) call isat_abort('cirxn8', 3,
     1	       mess='bad nxx: nxx, nc+2, mode_pdt=',
     2	       ivar = (/ nxx, nc+2, mode_pdt/) )
	elseif( mode_pdt == 2  .or.   mode_pdt == 3 ) then
	   if( nxx /= nc+1 ) call isat_abort('cirxn8', 4,
     1	       mess='bad nxx: nxx, nc+1, mode_pdt=',
     2	       ivar = (/ nxx, nc+1, mode_pdt/) )
	else
	   if( nxx /= nc ) call isat_abort('cirxn8', 5,
     1	       mess='bad nxx: nxx, nc, mode_pdt=',
     2	       ivar = (/ nxx, nc, mode_pdt/) )
     	endif

!-----------  set isat controls  ----------------------------------

	idtab  = 1
	mode   = 0
	nx     = nxx
	nf     = nc+1
	nh     = 1

	info     = info_n
	info( 1) = 1            !  iscale
	!info(12) set using ciparam
	!info(12) = 2            !  isatop
	info(13) = lu_err       !  error output
	info(66) = ci_info_n(15)          !  mpi_uniq

!  allocate space for f to be returned for modeci=9
	allocate( rusr(nrc+2) ) 
	rusr = nrc+2

	allocate( rinfo(50+nxx+nc+1) )
	rinfo = 0.d0
	rinfo(1:50) = rinfo_n

!  scaling for x

	rinfo(51:50+nrs+ne) = phiref
	rinfo(50+nc)  = hsref

	if( mode_pdt == 1 ) then
	   rinfo(50+nc+1)  = pressref 
	   !SBP 4/26/02:  
	   rinfo(50+nc+2)  = dtref
	elseif( mode_pdt == 2 ) then
	   rinfo(50+nc+1)  = pressref 
	elseif( mode_pdt == 3 ) then
		!SBP 4/26/02:
	   rinfo(50+nc+1)  = dtref
	endif

	allocate( xref(nx) )
	xref = rinfo(51:50+nx)
!  scaling for f

	rinfo(50+nx+1:50+nx+nrs+ne) = phiref 

	rinfo(50+nx+nc) = hsfref
	rinfo(50+nx+nc+1) = tempref

	isatab_called = .true.

!==========  end of initialization  ============================

60	continue

!-----------  check constant press and dt  ----------------------

	if( const_pr ) then
	   if( press /= prc ) call isat_abort('cirxn8', 6,
     1	      mess='pressure not constant', rvar=(/ press, prc/) )
     	endif

	if( const_dt ) then
	   if( t /= dtc ) call isat_abort('cirxn8', 7,
     1	      mess='dt not constant', rvar=(/ t, dtc/) )
     	endif

!-----------  set parameters for isat ---------------------------

!  x = { r(1:nrs+ne), h, p, t }  (at time 0),  mode_pdt = 1
!  x = { r(1:nrs+ne), h, p }     (at time 0),  mode_pdt = 2
!  x = { r(1:nrs+ne), h, t }     (at time 0),  mode_pdt = 3
!  x = { r(1:nrs+ne), h }        (at time 0),  mode_pdt = 4

!  f = { r(1:nrs+ne), h, T }     (at time t)

	x(1:nc) = c0

	if( mode_pdt == 1 ) then
	   x(nc+1) = press
	   x(nc+2) = t
	elseif( mode_pdt == 2 ) then
	   x(nc+1) = press
	elseif( mode_pdt == 3 ) then
	   x(nc+1) = t
	endif

!-----------  enforce realizability

! DR-08 enforce realizability may be needed  
! similar to modeci=6 call cireal( 2, treal, x(1:ns), ismall, vsmall )
! DR-08 a simple realizability check, if negative, force to be zero 	
	x(1:nc-1) = max(x(1:nc-1), 0.d0)

!--- call isatab
! if modeci ==8 call RCCE  and if modeci ==9 call ICE-PIC
	
	if ( modeci == 8 ) then
	   call isatab( idtab, mode, nx, x, nf, nh, 1, cirmap_dr_rcce,
     1	   iusr, rusr, info, rinfo, f, dfdx, hvar, stats )
	elseif( modeci == 9 ) then 
	   rusr(1) = nrc+2
	   rusr(2) = -huge(1.d0) !  set tol_re  XXXXXXXXXXXX
        call isatab( idtab, mode, nx, x, nf, nh, 1, cirmap_dr_ice_pic,
     1	   iusr, rusr, info, rinfo, f, dfdx, hvar, stats )
      endif
      

!--- test for non-retrieve
	if( nint( stats(9) ) == 8 ) then
	   dpt(1) = -1.d0  !  signal that f has not been returned
	   return
	endif


!  enforce realizability may be needed similar to modeci=6
!  call cireal( 2, treal, f(1:ns), ismall, vsmall )
      
!  DR-08 a simple realizability check, if negative, force to be zero 	

	ct(1:nc) = f(1:nc)
	ct(1:nc-1) = max(ct(1:nc-1), 0.d0)
	temp_a = f(nc+1)

	! set approximated temp and dens values
	dpt(3) = temp_a
	sz = sum(ct(1:nrs)) + dot_product(we, ct(nrs+1:nrc))
	dpt(1) = dpt(2)/(sz*gascon*temp_a)
	
	return
	end subroutine cirxn8

!=========================================================================
	subroutine ci_reset9( num_rs, rep, BT ) 

	use ci_dat8
	use ci_ice_cksubs
	implicit none
	
	integer, intent(in) :: num_rs, rep(num_rs)
	real(k_dp), intent(in) :: BT(num_rs+ne, ns)
	integer    :: i, j, CEtemp(ns, ne)
	real(k_dp) :: ewt(ne)
	
	nrs = num_rs
	nrc = nrs + ne
	nus = ns - nrs
	
	deallocate( amolwt_n )
	deallocate( href_n )
	deallocate( CS )
	deallocate( US )
	deallocate( CEu )
	deallocate( CEr )
	deallocate( BB )
	deallocate( BBT )
	deallocate( BBF )
	deallocate( BBFT )
	
	allocate( amolwt_n(nrc) )
	allocate( href_n(nrc) )
	allocate( CS(nrs) )
	allocate( US(nus) )
	allocate( CEr(nrs, ne))
	allocate( CEu(nus, ne))
	allocate( BB(ns, nrc) )
	allocate( BBT(nrc, ns) )
	allocate( BBF(ns, nrc) )
	allocate( BBFT(nrc, ns) )
	
	CS = rep(1:num_rs)
	indic_rs = 0
	S2RS = 0
	do i = 1,nrs
	   S2RS(CS(i)) = i
	   indic_rs(CS(i)) = 1
	   amolwt_n(i) = amolwt(CS(i))
	enddo
	call ctawt(  ewt)
	amolwt_n(nrs+1:nrs+ne)=ewt(1:ne)
	
	j = 0
	do i = 1, ns
	   if( indic_rs(i) == 0 ) then
	      j = j + 1
	      US(j) = i	!  j-th unrepresented species is species US(j)
	   endif
	end do


	CEtemp = CE
	do i = 1,nrs
	   CEr(i,1:ne) = CE(CS(i),1:ne)
	   CEtemp(CS(i),1:ne) = 0
	end do

	j = 0
	do i=1,ns
	   if (maxval(CEtemp(i,1:ne))>0) then
	      j=j+1 
	      CEu(j,1:ne)=CEtemp(i, 1:ne)
	   endif
	enddo

	call set_href_n( k_href_n )	! set the nominal enthalpies of formation
	BBT = BT
	BB = transpose(BT)
	BBF(1:ns, 1:nrs) = BB(1:ns, 1:nrs)
	BBF(1:ns, nrs+1:nrs+ne) = CE(1:ns, 1:ne)
	BBFT = transpose(BBF)
	
	end subroutine ci_reset9

	end module ci_8
