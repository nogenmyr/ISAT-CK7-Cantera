!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module isat_init

   use isat_cdf_action
   use isat_lu_m
   use isat_defaults
   use isat_change_m
   use isat_subs
   use isat_io

   private :: info_n, rinfo_n

   integer,    save, dimension(l_info)  :: info_n
   real(k_xf), save, dimension(l_rinfo) :: rinfo_n
           
contains

!====================================================================== 

   subroutine isat_table_init( idtab, nx, nf, nh, info, rinfo, table )

!  allocate table and initialize all sub-objects

   implicit none

   integer, intent(in)        :: idtab, nx, nf, nh, info(l_info)
   real(k_xf), intent(in)     :: rinfo(l_rinfo+nx+nf)
   type (table_type), pointer :: table

   integer       :: i, j, wleaf, lu_log, na, icheck, &
		            idproc, nchanges, n_pt_0=100, def_set
   real(k_xf)    :: stomby
   real          :: cputime
   character(30) :: blank, head, tail

!----------  check version of isatab  -------------------------------------

	  if( info(61) > 0 ) then  
	     if( info(61) /= isat_vers ) then
		    call isat_abort('isatab', -1, mess='Wrong version', ivar=(/info(61), isat_vers/) )
		 endif
	  endif

!----------  allocate table and set pointers  ----------------------------

   allocate(table)

   nullify(table%next_table)
   table%leaves   = 0 

!------  check scalar input arguments  ------------------------------------

   if( idtab < 0 ) call isat_abort( 'isat_table_init', 1, &
      mess = 'idtab must be >= 0: idtab = ', isv = idtab )

   if( nx <= 0 ) call isat_abort( 'isat_table_init', 2, &
      mess = 'nx must be > 0: nx = ', isv = nx )

   if( nf <= 0 ) call isat_abort( 'isat_table_init', 3, &
      mess = 'nf must be > 0: nf = ', isv = nf )

   if( nh < 0 ) call isat_abort( 'isat_table_init', 4, &
      mess = 'nh must be >= 0: nh = ', isv = nh )

!---------  set permanent quantities from arguments  ------------------------

   table%idtab  = idtab
   table%nx     = nx
   table%nf     = nf
   table%nh     = nh

   table%nxp1   = nx + 1
   table%nx2    = 2 * nx
   table%jfe    = nx + nf
   table%jhs    = nx + nf + 1
   table%jhe    = nx + nf + nh

!----------  set defaults for info -----------------------------

   def_set = info(15) !  default settings
   call isat_info_defaults( def_set )

!------  check values from info and set corresponding table variables -------
!        valueof info used is stored in  info_n

   info_n = info

   call isat_iset( table%iscale,        1 )
   call isat_iset( table%if_g,          2 )
   call isat_iset( table%if_h,          3 )
   call isat_iset( table%pr_mru,        4 )
   call isat_iset( table%pr_mfu,        5 )
   call isat_iset( table%ret_min,       6 )
   call isat_iset( table%ret_max,       7 )
   call isat_iset( table%grow_min,      8 )
   call isat_iset( table%grow_max,      9 )
   call isat_iset( table%ichin,        10 )
   call isat_iset( table%ichout,       11 )
   call isat_iset( table%isatop,       12 )

   call isat_iset( table%ich0,         19 )
   call isat_iset( table%stats0,       20 )
   call isat_iset( table%ifull,        21 )
   call isat_iset( table%istats,       22 )
   call isat_iset( table%no_pr,        23 )
   call isat_iset( table%no_sr,        24 )
   call isat_iset( table%no_grow,      25 )
   call isat_iset( table%no_add,       26 )
   call isat_iset( table%no_dev,       27 )
   call isat_iset( table%idites,       28 )
   call isat_iset( table%defer,        29 )
   call isat_iset( table%no_bt,        30 )

   call isat_iset( table%grow_mode,    35 )
   call isat_iset( table%mode_con,     36 )
   call isat_iset( table%mode_eoi,     37 )
   call isat_iset( table%mode_shrink,  38 )

   call isat_iset( table%degen_r,      39 )
   call isat_iset( table%degen_g,      40 )
   call isat_iset( table%degen_s,      41 )
   call isat_iset( table%degen_c,      42 )
   call isat_iset( table%ecc_act,      43 )

   call isat_iset( table%ad_look,      44 )
   call isat_iset( table%ad_use,       45 )
   call isat_iset( table%ad_shrink,    46 )
   
   call isat_iset( table%de_nearby,    47 )

   call isat_iset( table%aff_na,       48 )
   call isat_iset( table%aff_lim,      49 )

   call isat_iset( table%sr_mode,      50 )
   call isat_iset( table%gr_mode,      51 )
   call isat_iset( table%pair_cover,   52 )
   call isat_iset( table%cut_add,      53 )
   call isat_iset( table%cut_rem,      54 )
   call isat_iset( table%cut_upd,      55 )

   call isat_iset( table%icheck,       60 )
   call isat_iset( table%isat_vers,    61 )
   call isat_iset( table%n_spool,      62 )
   call isat_iset( table%n_add_op,     63 )

   call isat_iset( table%mpi_log,      65 )
   call isat_iset( table%mpi_uniq,     66 )
   call isat_iset( table%mpi_nml,      67 )

   call isat_iset( table%cdf_error,    70 )


!---------  check and store quantities from  rinfo  --------------------------

   rinfo_n = rinfo
   call isat_rset( table%etola,         1 )
   call isat_rset( table%etolr,         2 )
   call isat_rset( table%etolc,         3 )
   call isat_rset( table%eoa_lim,       4 )
   call isat_rset( table%eoi_lim,       5 )
   call isat_rset( table%ret_frac,      6 )
   call isat_rset( table%grow_frac,     7 )
   call isat_rset( table%stomby,        8 )
   call isat_rset( table%outinc,        9 )
   call isat_rset( table%chk_inc,      10 )
   call isat_rset( table%chk_full,     11 )
   call isat_rset( table%aff_adds,     12 )
   call isat_rset( table%aff_inc,      13 )
   call isat_rset( table%aff_thresh,   14 )
   call isat_rset( table%aff_ratio,    15 )
   call isat_rset( table%ad_theta,     16 )
   call isat_rset( table%pair_sep,     17 )
   call isat_rset( table%slow_prog_thresh, 18 )
   call isat_rset( table%ret_frac_slow,    19 )
   call isat_rset( table%grow_frac_slow,   20 )
   call isat_rset( table%ecc_freq,         21 )

!---------  initialize MPI  --------------------------------------------

   call isat_mpi_rank( table%myrank, table%nprocs )

!---------  determine files needed  ------------------------------------

   table%write_log = .true.  ! .true.  for log file for this process
   table%my_isatop = 0       ! =1 for op file for this process
   table%my_ichout = 0       ! =1 for checkpoint file

   if( table%myrank == 0 ) then
      if( table%isatop > 0 ) table%my_isatop = 1
      if( table%ichout > 0 ) table%my_ichout = 1
   else
      if( table%isatop == 2 ) table%my_isatop = 1
	  if( table%ichout == 2 ) table%my_ichout = 1
	  if( table%mpi_log == 1 ) table%write_log = .false.
   endif

!----------  create file names ---------------------------------------------

   idproc = table%myrank  ! = rank for unique names; = -1 otherwise
   if( table%mpi_uniq == 1 ) idproc = -1  
   if( table%nprocs == 1 ) then
      idproc = -1
      if( table%mpi_uniq == 2 ) idproc = 0
   endif
   table%idproc = idproc

   blank = repeat(' ',30)
   head  = blank
   head  = 'isat'

   tail = blank   !  generate file name:  isat_#_P.log
   tail = 'log'
   call isat_file_name( head, idtab, idproc, tail, table%isat_log )

   tail = blank   !  generate file name:  isat_#_P.op
   tail = 'op'
   call isat_file_name( head, idtab, idproc, tail, table%isat_op  )

   tail = blank   !  generate file name:  isat_#_P.op
   tail = 'leaf'
   call isat_file_name( head, idtab, idproc, tail, table%isat_leaf  )

   tail  = blank      !  isat_tab = 'isat_n_p.tab'
   tail  = 'tab'
! VH - 05/11/2012 - add option to save table with specified id
   if( info(72) == 0 ) then
      call isat_file_name( head, idtab, idproc, tail, table%isat_tab )
   else
      call isat_file_name( head, idtab, info(72), tail, table%isat_tab )
      write(0,*) 'isat_tab: rank, partition id =', idproc, info(72)
   endif

   tail  = blank       !  isat_nml = 'isat_n_p.nml'
   tail  = 'nml'
! VH - 11/28/2012 - use only one input file if mpi_nml /= 0 
   if( table%mpi_nml == 0 ) then
      call isat_file_name( head, idtab, idproc, tail, table%isat_nml )
   else
      call isat_file_name( head, idtab, -1, tail, table%isat_nml )
   endif

   tail  = blank       !  isat_err = 'isat_n_p.nml'
   tail  = 'err'
   call isat_file_name( head, idtab, idproc, tail, table%isat_err )

!----------  open files as needed ---------------------------------

   if( info(13) < 0 ) then  !  error file to be created
      call isat_lu( lu_err )
      table%lu_err = lu_err    !  open error file on logical unit lu_err
      open( lu_err, file = table%isat_err , status = 'replace', &
            action = 'write' )
   else
      table%lu_err = info(13)
   endif

   if( table%write_log ) then  !  log file to be written
      call isat_lu( lu_log )
      table%lu_log = lu_log    !  open log file on logical unit lu_log
      open( lu_log, file = table%isat_log , status = 'replace', &
            action = 'write' )
   else
      lu_log         = 0
      table%lu_log   = 0
   endif

   if( table%my_isatop == 1 ) then !  isat_op = 'isat_n_p.op'
      call isat_lu( table%lu_op )
      open( table%lu_op , file = table%isat_op  , status = 'replace', &
               action = 'write', recl = 4000 )
      table%isatop_everon = .true.

   else
      table%lu_op   = 0
      table%lu_leaf = 0
      table%isatop_everon = .false.
   endif

   call isat_lu( table%lu_leaf ) 
   table%leaf_everon = .false.

!---------  start output to log file  ---------------------------------

   if( table%write_log ) then  
      write(lu_log,*)'ISATAB Version ', isatab_version
      write(lu_log,"(/a,a,i6)")'isat_table_init: initializing table ', &
                            'for idtab=', idtab
   endif

!---  reset info_n and rinfo_n and table% constants from isat_#_P.nml  ----

   call isat_const(table)

!---------  check and store quantities from  rinfo  --------------------------

   table%etol_consq = table%etolc**2

!------------  output quantities from info and rinfo  ------------------------

   if( table%write_log ) then
      write(lu_log,4) table%idtab, table%nprocs, table%myrank
      write(lu_log,5) nx, nf, nh
	  write(lu_log,'(/,a)') 'Values of info(1:70) to be used:'
      write(lu_log,7) info_n(1:70)
	  write(lu_log,'(/,a)') 'Values of rinfo(1:25) to be used:'
      write(lu_log,9) rinfo_n(1:25)
   endif

4  format(/'idtab, nprocs, myrank=', 3i6 )                                  
5  format(/'nx, nf, nh=', 3i6 )                                  
7  format(/(5i8)) 
9  format(/(1p,5e12.3))

!---------------  store scale factors --------------------------------------

   allocate( table%xscale(nx) )
   allocate( table%fscale(nf) )
   allocate( table%xsci(nx) )
   allocate( table%fsci(nf) )
   allocate( table%g2gs(nf,nx) )
   allocate( table%gs2g(nf,nx) )
   allocate( table%p(nx,table%ad_use) )
   allocate( table%dsq(table%ad_look) )
   allocate( table%ids(table%ad_look) )

   if( table%iscale == 0 ) then
      table%xscale = 1.
      table%fscale = 1.
      table%xsci   = 1.
      table%fsci   = 1.
      table%g2gs   = 1.
      table%gs2g   = 1.
   else
      table%xscale = rinfo(l_rinfo+1:l_rinfo+nx)
      table%fscale = rinfo(l_rinfo+1+nx:l_rinfo+nx+nf)

      if( table%write_log ) then
         write(lu_log,10) rinfo(l_rinfo+1:l_rinfo+nx)
         write(lu_log,11) rinfo(l_rinfo+1+nx:l_rinfo+nx+nf)
      endif

10    format(/'xscale=',1p,5e12.3,/,('       ',1p,5e12.3))
11    format(/'fscale=',1p,5e12.3,/,('       ',1p,5e12.3))

      if( minval(table%xscale) <= 0. ) call isat_abort( &
	     'isat_table_init', 1, rvar = table%xscale, &
         mess = 'xscale must be positive; xscale = ' )

      if( minval(table%fscale) <= 0. ) call isat_abort( &
	     'isat_table_init', 1, rvar = table%fscale, &
         mess = 'fscale must be positive; fscale = ' )

      table%xsci = 1. / table%xscale
      table%fsci = 1. / table%fscale

      do i   = 1, nf
      do j   = 1, nx
         table%g2gs(i,j)  = table%xscale(j) / table%fscale(i)
         table%gs2g(i,j)  = table%fscale(i) / table%xscale(j)
      end do
      end do
   endif

!------------ set affine space -------------------------------------------
   
   call isat_na_max( table%aff_lim, nx, table%na_max )  !  set na_max, upper limit on na

   if( table%aff_na == 0 ) then
      na = 1
   elseif( table%aff_na < 0 ) then
      na = nx
   elseif( table%aff_na >= nx ) then
      table%aff_na = -1
	  na = nx
   else
	  na = table%aff_na
   endif

   table%na = na
   allocate( table%ua(nx,nx+1) )
   table%ua = 0.d0

   do i = 1, nx
      table%ua(i,i) = 1.d0
   end do

!------------ initialize ISAT5 data structures  --------------------------

   table%ng  = (nx*(nx+1))/2
   table%nga = (na*(na+1))/2

   icheck = table%icheck
   call id_init( table%idlist, check=icheck )

   if( table%ichin == 0 ) then  !  do not create here if they will be read
      call sell_create( table%seoa,  nx, check=icheck, idlist=table%idlist )
      call sell_create( table%speoa, na, check=icheck, idlist=table%idlist )
      call sell_create( table%seoi,  nx, check=icheck, idlist=table%idlist )
      call sell_create( table%speoi, na, check=icheck, idlist=table%idlist )

      call bt_initialize( table%seoa, table%eoa_bt, check=icheck, idlist=table%idlist )

      call ebt_initialize( table%speoa, table%peoa_ebt, check=icheck, idlist=table%idlist, pair_sep_qual=table%pair_sep )
      call ebt_initialize( table%speoi, table%peoi_ebt, check=icheck, idlist=table%idlist, pair_sep_qual=table%pair_sep )

      call ll_initialize( table%eoa_mru, check=icheck, idlist=table%idlist )
      call ll_initialize( table%eoa_mfu, check=icheck, idlist=table%idlist )
   endif

   allocate( table%leaf_pt(n_pt_0) )
   table%n_pt = n_pt_0

   do i = 1, n_pt_0
      nullify( table%leaf_pt(i)%leaf )
   end do

!-----------  allocate query ---------------------------------------------

   allocate( table%query )
   allocate( table%query%xs(nx) )
   allocate( table%query%xa(nx) )
   allocate( table%query%fs(nf) )
   table%query%table => table

!----------  set maxleaves for given storage  ----------------------------

!  calculate (approximate) size of one leaf plus associated storage

   stomby = table%stomby
      
   if( stomby > 0.d0 ) then
      call isat_leaves( stomby, nx, nf, nh, table%na_max, 0, 0, wleaf )
   else
      wleaf = nint( -stomby )
   endif
   table%maxleaves = max( wleaf, 1 )

   if( table%write_log ) write(lu_log,"(/a,i8)")'maximum number of table entries = ',&
                             table%maxleaves
  
   if( table%maxleaves < 1 ) call isat_abort( 'isat_table_init', 1, &
      mess = 'Too few leaves.  Increase rinfo(9)=stomby. maxleaves = ', &
      isv = table%maxleaves )

!----------------  allocate CDF and initialize as needed

   allocate( table%cdf_err ) 
   table%cdf_err%initialized  = .false. 

   if( table%cdf_error > 0 ) then
      info_n(80) = 1
      call isat_cdf_act( table, 25, info_n, rinfo_n )
   else
	  table%cdf_err%on  = .false.
   endif

!----------------  initialize remaining quantities

   table%full   = .false.
   table%nextop = 0.
   if( table%isatop == 0 ) table%nextop = huge(1.d0)

   table%i_spool = 0
   if( table%n_spool > 1 ) then
	  allocate( table%spool(200, table%n_spool) )
   else
      nullify( table%spool )
   endif

   table%nextch = table%ich0
   if( table%ichout == 0 ) table%nextch = huge(1.d0)

   call cpu_time( cputime )
   table%cpu_end = cputime

   call system_clock( count = table%wall_0 )  !  wall clock time
   table%wall_last  = table%wall_0
   table%wall_cycle = 0

   table%deferred = 0.d0
   
   table%ecc_q = huge(1.d0)  !  ECC 
   if( table%ecc_freq > 0.d0 ) then  !  based on growing
      table%ecc_q = 1.d1 / table%ecc_freq
   elseif( table%ecc_freq > -1.d0  .and.  table%ecc_freq < 0.d0 ) then  !  based on probability
      table%ecc_q = -1.d0 / table%ecc_freq 
   endif

!----------  initialize stats ----------------------------------------------

   call isat_stats_init( table, 0 )

!----------  read in checkpointed table?   ---------------------------------

   if( table%ichin > 0 ) then
      call isat_table_read( table )  !  read in old table

	  call isat_change( table, info_n, rinfo_n, nchanges )  !  reset parameters if changed

	  if( table%write_log ) then
	     write(table%lu_log,*)' '
	     if( nchanges == 0 ) then
		    write(table%lu_log,'(a)') 'Table read in has identical parameters to those set.'
		 elseif( nchanges == 1 ) then
		    write(table%lu_log,'(a)') 'Table read in has 1 parameter that differs from those set.'
		 else
		    write(table%lu_log,'(a,i2,a)') 'Table read in has', nchanges, &
			                          ' parameters that differ from those set.'
		 endif
		 write(table%lu_log,*)' '
      endif

	  if( table%stats0 == 1 ) call isat_stats_init( table, 0 )

      table%nextch = max( table%nextch , table%chk_inc * table%leaves )
   endif

   if( table%write_log) call isat_flush( lu_log )

   return
end subroutine isat_table_init

!===========================================================================

subroutine isat_iset( iset, i )

! Set the value of iset to info_n(i), or to idef(i) if info_n(i)=0. 
! Reset info_n(i) = iset.
! If imax>0, check that info_n(i) satisfies  imin <= info_n(i) <= imax.
! If imax=0, check that info_n(i) is non-negative.
! (No check performed for imax < 0.)
 
   implicit none
   integer, intent(in)  :: i
   integer, intent(out) :: iset

   integer :: idef, imin, imax

   iset = info_n(i)
   idef = info_def(i)
   imin = info_min(i)
   imax = info_max(i)

   if( imax > 0 ) then
      if( iset < imin  .or.  iset > imax ) call isat_abort('isat_iset', 1, &
	                      mess = 'invalid value of, info(i): i = ', isv = i )
   elseif( imax == 0 ) then
      if( iset <0 ) call isat_abort('isat_iset', 2, &
	                      mess = 'negative value of, info(i): i = ', isv = i ) 
   endif
   
   if( iset == 0 ) iset = idef
   info_n(i) = iset
   
   return

end subroutine isat_iset

!===========================================================================

subroutine isat_rset( rset, i )

! Set the value of rset to rinfo_n(i), or to rinfo_def if rinfo_n(i)=0. 
! Reset rinfo_n(i) = rset.
! Check that rinfo_n(i) is non-negative.
 
   implicit none
   integer, intent(in)     :: i
   real(k_xf), intent(out) :: rset

   rset = rinfo_n(i)
   if( rset == 0. ) rset = rinfo_def(i)

   if( rset < rinfo_min(i) ) call isat_abort('isat_rset', 1, &
	                      mess = 'invalid value of, rinfo(i): i = ', isv = i )

   rinfo_n(i) = rset
   
   return

end subroutine isat_rset


!===========================================================================

subroutine isat_stats_init( table, mode )

!  Initialize stats.
!  mode = 0 - initialize
!  mode = 1 - re-initialize (except for stats(1))
!  mode = 2 - re-initialize and rewind isat_#.op
!  (stats_hist is reset)

   implicit none
   
   type (table_type), pointer :: table
   integer, intent(in)        :: mode

   if( mode /= 1 ) table%stats(1) = 0.d0
   table%stats(2:100) = 0.d0
   table%stats(13)    = table%maxleaves

   table%stats_hist   = 0.d0
   table%i_hist       = 0

   if( mode /=2 ) return

   if( table%my_isatop == 0 ) return
   close( table%lu_op )
   open( table%lu_op , file = table%isat_op  , status = 'replace', &
               action = 'write', recl = 2000 )

   return
end subroutine isat_stats_init

!===========================================================================

subroutine isat_const( table )

!  This routine has three functions:
! 1/ if the file  table%isat_nml exists,  set info and rinfo (and
!    other) variables specified in that file.
! 2/ set variables in table%mpi
! 3/ set logical unit limits (luf, lul) and set  flush  on if required

implicit none
type (table_type), pointer :: table

!----------- namelist variables

integer    :: iscale, if_g, if_h, pr_mru, pr_mfu, ret_min, &
              ret_max, grow_min, grow_max, ichin, ichout, isatop, &
			  ich0, stats0, ifull, istats, no_pr, no_sr, no_grow, no_add, &
			  no_dev, idites, defer, no_bt, grow_mode,  mode_con, mode_eoi, &
			  mode_shrink, degen_r, degen_g, degen_s, degen_c, ecc_act, ad_look, ad_use, &
			  ad_shrink, de_nearby, aff_na, aff_lim, sr_mode, gr_mode, &
			  pair_cover, cut_add, cut_rem, cut_upd, icheck, isat_vers, &
			  n_spool, n_add_op, mpi_log, mpi_uniq, mpi_nml, cdf_error

integer    :: ngsb, pend_max, nasb, luf, lul, if_flush, nml_op

real(k_xf) :: etola, etolr, etolc, eoa_lim, eoi_lim, ret_frac, grow_frac, &
              stomby, outinc, chk_inc, chk_full, aff_adds, aff_inc, &
			  aff_thresh, aff_ratio, ad_theta, pair_sep, slow_prog_thresh, &
			  ret_frac_slow, grow_frac_slow, ecc_freq

!----------- other variables

integer    :: lu
logical    :: exist
               
namelist / isat_nml /  nml_op, &
              if_g, if_h, pr_mru, pr_mfu, ret_min, &
              ret_max, grow_min, grow_max,  &
			  ich0, stats0, ifull, istats, no_pr, no_sr, no_grow, no_add, &
			  no_dev, idites, defer, no_bt, grow_mode,  mode_con, mode_eoi, &
			  mode_shrink, degen_r, degen_g, degen_s, degen_c, ecc_act, ad_look, ad_use, &
			  ad_shrink, de_nearby, aff_na, aff_lim, sr_mode, gr_mode, &
			  pair_cover, cut_add, cut_rem, cut_upd, icheck, n_spool, n_add_op,  &
			  etola, etolr, etolc, eoa_lim, eoi_lim, ret_frac, grow_frac, &
			  stomby, outinc, chk_inc, chk_full, aff_adds, aff_inc, &
			  aff_thresh, aff_ratio, ad_theta, pair_sep, slow_prog_thresh, &
			  ret_frac_slow, grow_frac_slow, ecc_freq

!---------  set  mpi  and  flush  constants to default values

   ngsb         = 10    ! size of the grow-send buffer
   pend_max     = 1000  ! number of pending messages allowed
   nasb         = 10    ! size of add-send buffer

   luf          = 60    ! first logical unit number in range to be used
   lul          = 99    ! last logical unit number in range to be used
   if_flush     =  0    ! set to 1 to turn on isat_flush
   nml_op       =  0    ! set /=0 to produce output on isat_log.#

!-----  copy from info_n, rinfo_n to local variables  -------------

   call iset1( iscale,        1 )
   call iset1( if_g,          2 )
   call iset1( if_h,          3 )
   call iset1( pr_mru,        4 )
   call iset1( pr_mfu,        5 )
   call iset1( ret_min,       6 )
   call iset1( ret_max,       7 )
   call iset1( grow_min,      8 )
   call iset1( grow_max,      9 )
   call iset1( ichin,        10 )
   call iset1( ichout,       11 )
   call iset1( isatop,       12 )

   call iset1( ich0,         19 )
   call iset1( stats0,       20 )
   call iset1( ifull,        21 )
   call iset1( istats,       22 )
   call iset1( no_pr,        23 )
   call iset1( no_sr,        24 )
   call iset1( no_grow,      25 )
   call iset1( no_add,       26 )
   call iset1( no_dev,       27 )
   call iset1( idites,       28 )
   call iset1( defer,        29 )
   call iset1( no_bt,        30 )

   call iset1( grow_mode,    35 )
   call iset1( mode_con,     36 )
   call iset1( mode_eoi,     37 )
   call iset1( mode_shrink,  38 )

   call iset1( degen_r,      39 )
   call iset1( degen_g,      40 )
   call iset1( degen_s,      41 )
   call iset1( degen_c,      42 )
   call iset1( ecc_act,      43 )

   call iset1( ad_look,      44 )
   call iset1( ad_use,       45 )
   call iset1( ad_shrink,    46 )
   
   call iset1( de_nearby,    47 )

   call iset1( aff_na,       48 )
   call iset1( aff_lim,      49 )

   call iset1( sr_mode,      50 )
   call iset1( gr_mode,      51 )
   call iset1( pair_cover,   52 )
   call iset1( cut_add,      53 )
   call iset1( cut_rem,      54 )
   call iset1( cut_upd,      55 )

   call iset1( icheck,       60 )
   call iset1( isat_vers,    61 )
   call iset1( n_spool,      62 )
   call iset1( n_add_op,     63 )

   call iset1( mpi_log,      65 )
   call iset1( mpi_uniq,     66 )
   call iset1( mpi_nml,      67 )

   call iset1( cdf_error,    70 )

   call rset1( etola,         1 )
   call rset1( etolr,         2 )
   call rset1( etolc,         3 )
   call rset1( eoa_lim,       4 )
   call rset1( eoi_lim,       5 )
   call rset1( ret_frac,      6 )
   call rset1( grow_frac,     7 )
   call rset1( stomby,        8 )
   call rset1( outinc,        9 )
   call rset1( chk_inc,      10 )
   call rset1( chk_full,     11 )
   call rset1( aff_adds,     12 )
   call rset1( aff_inc,      13 )
   call rset1( aff_thresh,   14 )
   call rset1( aff_ratio,    15 )
   call rset1( ad_theta,     16 )
   call rset1( pair_sep,     17 )
   call rset1( slow_prog_thresh, 18 )
   call rset1( ret_frac_slow,    19 )
   call rset1( grow_frac_slow,   20 )
   call rset1( ecc_freq,         21 )

!---------  if file  isat.nml exists, read namelist

   inquire( file = table%isat_nml, exist = exist )
   if( exist ) then
      call isat_lu( lu )
      open( lu, file = table%isat_nml, position = 'rewind', action = 'read', err = 100 )
      read( lu, nml = isat_nml, err=101 )
      close( lu )

!   report changes
      if( table%write_log .and. nml_op /= 0 ) then
		 	write( table%lu_log, * )'Read from namelist file: ', table%isat_nml
            write( table%lu_log, * )' '
            write( table%lu_log, * )'------ Values of isatab parameters:&
                                  & after modification by isat_nml -----'
            write( table%lu_log, isat_nml )
            write( table%lu_log, * )'----------------------------------&
                                  &----------------------------------'
      endif

   elseif( table%write_log ) then
	  write( table%lu_log, '(a)' )'Namelist file not found: ', table%isat_nml
   endif

!--------  treat  flush

   call isat_lu_set( luf, lul )
   if( if_flush == 1 ) call isat_flush( -1 )

!-----  copy from local variables to table% , info_n, rinfo_n  -------------

   call iset2( iscale,       table%iscale,        1 )
   call iset2( if_g,         table%if_g,          2 )
   call iset2( if_h,         table%if_h,          3 )
   call iset2( pr_mru,       table%pr_mru,        4 )
   call iset2( pr_mfu,       table%pr_mfu,        5 )
   call iset2( ret_min,      table%ret_min,       6 )
   call iset2( ret_max,      table%ret_max,       7 )
   call iset2( grow_min,     table%grow_min,      8 )
   call iset2( grow_max,     table%grow_max,      9 )
   call iset2( ichin,        table%ichin,        10 )
   call iset2( ichout,       table%ichout,       11 )
   call iset2( isatop,       table%isatop,       12 )

   call iset2( ich0,         table%ich0,         19 )
   call iset2( stats0,       table%stats0,       20 )
   call iset2( ifull,        table%ifull,        21 )
   call iset2( istats,       table%istats,       22 )
   call iset2( no_pr,        table%no_pr,        23 )
   call iset2( no_sr,        table%no_sr,        24 )
   call iset2( no_grow,      table%no_grow,      25 )
   call iset2( no_add,       table%no_add,       26 )
   call iset2( no_dev,       table%no_dev,       27 )
   call iset2( idites,       table%idites,       28 )
   call iset2( defer,        table%defer,        29 )
   call iset2( no_bt,        table%no_bt,        30 )

   call iset2( grow_mode,    table%grow_mode,    35 )
   call iset2( mode_con,     table%mode_con,     36 )
   call iset2( mode_eoi,     table%mode_eoi,     37 )
   call iset2( mode_shrink,  table%mode_shrink,  38 )

   call iset2( degen_r,      table%degen_r,      39 )
   call iset2( degen_g,      table%degen_g,      40 )
   call iset2( degen_s,      table%degen_s,      41 )
   call iset2( degen_c,      table%degen_c,      42 )
   call iset2( ecc_act,      table%ecc_act,      43 )

   call iset2( ad_look,      table%ad_look,      44 )
   call iset2( ad_use,       table%ad_use,       45 )
   call iset2( ad_shrink,    table%ad_shrink,    46 )
   
   call iset2( de_nearby,    table%de_nearby,    47 )

   call iset2( aff_na,       table%aff_na,       48 )
   call iset2( aff_lim,      table%aff_lim,      49 )

   call iset2( sr_mode,      table%sr_mode,      50 )
   call iset2( gr_mode,      table%gr_mode,      51 )
   call iset2( pair_cover,   table%pair_cover,   52 )
   call iset2( cut_add,      table%cut_add,      53 )
   call iset2( cut_rem,      table%cut_rem,      54 )
   call iset2( cut_upd,      table%cut_upd,      55 )

   call iset2( icheck,       table%icheck,       60 )
   call iset2( isat_vers,    table%isat_vers,    61 )
   call iset2( n_spool,      table%n_spool,      62 )
   call iset2( n_add_op,     table%n_add_op,     63 )

   call iset2( mpi_log,      table%mpi_log,      65 )
   call iset2( mpi_uniq,     table%mpi_uniq,     66 )
   call iset2( mpi_nml,      table%mpi_nml,      67 )

   call iset2( cdf_error,    table%cdf_error,    70 )

   call rset2( etola,        table%etola,         1 )
   call rset2( etolr,        table%etolr,         2 )
   call rset2( etolc,        table%etolc,         3 )
   call rset2( eoa_lim,      table%eoa_lim,       4 )
   call rset2( eoi_lim,      table%eoi_lim,       5 )
   call rset2( ret_frac,     table%ret_frac,      6 )
   call rset2( grow_frac,    table%grow_frac,     7 )
   call rset2( stomby,       table%stomby,        8 )
   call rset2( outinc,       table%outinc,        9 )
   call rset2( chk_inc,      table%chk_inc,      10 )
   call rset2( chk_full,     table%chk_full,     11 )
   call rset2( aff_adds,     table%aff_adds,     12 )
   call rset2( aff_inc,      table%aff_inc,      13 )
   call rset2( aff_thresh,   table%aff_thresh,   14 )
   call rset2( aff_ratio,    table%aff_ratio,    15 )
   call rset2( ad_theta,     table%ad_theta,     16 )
   call rset2( pair_sep,     table%pair_sep,     17 )
   call rset2( slow_prog_thresh, table%slow_prog_thresh, 18 )
   call rset2( ret_frac_slow,    table%ret_frac_slow,    19 )
   call rset2( grow_frac_slow,   table%grow_frac_slow,   20 )
   call rset2( ecc_freq,      table%ecc_freq,    21 )

   return

100   call isat_abort('isat_const', 1, mess='error opening file',chv=table%isat_nml )
101   call isat_abort('isat_const', 2, mess='error reading file',chv=table%isat_nml )

contains  !-------------------------------------

   subroutine iset1( ivar, i )
   integer :: ivar, i

      ivar = info_n(i)

   return
   end subroutine iset1
   
   subroutine iset2( ivar, tvar, i )
   integer :: ivar, tvar, i

      tvar      = ivar
      info_n(i) = ivar

   return
   end subroutine iset2

   subroutine rset1( rvar, i )
   integer    :: i
   real(k_xf) :: rvar

      rvar = rinfo_n(i)

   return
   end subroutine rset1  
   
   subroutine rset2( rvar, tvar, i )
   integer    :: i
   real(k_xf) :: rvar, tvar

      tvar       = rvar
      rinfo_n(i) = rvar

   return
   end subroutine rset2

end subroutine isat_const

!===========================================================================

end module isat_init
