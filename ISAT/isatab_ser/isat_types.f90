!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module isat_types
!  This module contains most of the data structures used in isatab.
!------------------------------------------------------------------------------

use isat_abort_m  ! this is included so that it is available by use association
use isat_cdf
use sell_ll
use sell_ebt
use sell_bt

implicit none

character(3)       :: isatab_version = '5.1'
integer            :: isat_vers      = 51
integer, parameter :: n_hist = 1000  ! number of records of stats history
integer, parameter :: l_info=100, l_rinfo=50  ! dimensions of info/rinfo

!------------------------------------------------------------------------------
type :: table_type

! The tables are stored in a linked list.  
! This table has  table%leaves  leaves, each of which has a unique positive integer ID.
! Collectively, the leaves are SELLs, and they have various data structures defined on them.
! In isat_table_init, all fixed data are stored in this structure.

   type (table_type),pointer  :: next_table ! pointer to next table

   type (cdf_type),  pointer  :: cdf_err    !  CDF or retrieve error

   integer                    :: leaves     ! total number of leaves

!--- dimensions and indexes
   integer :: nx, nf, nh
   integer :: na         ! dimension of affine space
   integer :: ng         ! Cholesky storage = (nx*(nx+1))/2
   integer :: nga        ! Cholesky storage = (na*(na+1))/2
   integer :: n_pt       ! dimension of leaf_pt(:)
   integer :: nxp1, nx2, jfe, jhs, jhe

!--- pointer to leaves
   type (id_list_type), pointer :: idlist     ! ID list for leaves
   type (leaf_pointer), pointer :: leaf_pt(:) ! pointer to leaves

!--- affine space used to define projected ELLs (ellipsoids)
   real(k_xf), pointer          :: ua(:,:)    ! affine space basis vectors

!--- SELLs (sets of ellipsoids) defined on leaves
   type (sell_type),    pointer :: seoa       ! SELL of EOAs
   type (sell_type),    pointer :: speoa      ! SELL of PEOAs
   type (sell_type),    pointer :: seoi       ! SELL of EOIs
   type (sell_type),    pointer :: speoi      ! SELL of PEOIs

!--- data structures defined on SELLs

   type (bt_type),      pointer :: eoa_bt     ! BT for EOAs
   type (ll_type),      pointer :: eoa_mru    ! MRU linked list of leaves
   type (ll_type),      pointer :: eoa_mfu    ! MFU linked list of leaves
   type (ebt_type),     pointer :: peoa_ebt   ! EBT for PEOAs
   type (ebt_type),     pointer :: peoi_ebt   ! EBT for PEOIs

!--- data structure for query information

   type (query_type),   pointer :: query

!--- quantities taken from arguments and info

   integer    :: idtab, iscale, if_g, if_h, pr_mru, pr_mfu, ret_min, &
                 ret_max, grow_min, grow_max, ichin, ichout, isatop, &
                 ich0, stats0, ifull, istats, no_pr, no_sr, no_grow, no_add, &
				 no_dev, idites, defer, no_bt, grow_mode,  mode_con, mode_eoi, &
				 mode_shrink, degen_r,degen_g, degen_s, degen_c, ecc_act, ad_look, ad_use, &
				 ad_shrink, de_nearby, aff_na, aff_lim, sr_mode, gr_mode, &
				 pair_cover, cut_add, cut_rem, cut_upd, icheck, isat_vers, &
				 n_spool, n_add_op, mpi_log, mpi_uniq, cdf_error
				 
!--- quantities taken from rinfo

   real(k_xf) :: etola, etolr, etolc, eoa_lim, eoi_lim, ret_frac, grow_frac, &
                 stomby, outinc, chk_inc, chk_full, aff_adds, aff_inc, &
				 aff_thresh, aff_ratio, ad_theta, pair_sep, slow_prog_thresh, &
				 ret_frac_slow, grow_frac_slow, ecc_freq

   real(k_xf), pointer :: xscale(:), fscale(:)

!--- fixed quantities set in isat_table_init

   real(k_xf), pointer :: xsci(:), fsci(:), g2gs(:,:), gs2g(:,:)
   real(k_xf)          :: etol_consq, r0_eoa, r0_eoi
   character(30)       :: isat_op, isat_log, isat_tab, isat_nml, isat_err, &
                          isat_leaf
   integer             :: lu_op, lu_log, lu_dat, lu_err, lu_leaf, &
                          maxleaves,  nprocs,  myrank, idproc, na_max, &
						  my_isatop, my_ichout, wall_0, wall_cycle, wall_last
   logical             :: write_log, isatop_everon, leaf_everon

!---  quantities that change during execution

   real(k_xf), pointer :: p(:,:), dsq(:), spool(:,:)
   integer, pointer    :: ids(:)

   real(k_xf)  :: stats(100), stats_hist(100,n_hist), deferred(3), nextop, &
                  nextch, cpu_end, ecc_q
   logical     :: full
   integer     :: i_hist, i_spool

end type table_type


!------------------------------------------------------------------------------
type :: leaf_type
   integer                        :: id      ! ID of leaf  ! ISAT5
   type (table_type), pointer     :: table
   real(k_xf)           :: etolsq
   real(k_xf),  pointer :: xfh(:)	 ! x and f are scaled
   real(k_g),   pointer :: g(:,:)   
   real(k_xf)           :: props(10) ! properties -- see isatab_leaf_props

!  Note: to reduce the number of pointers to arrays, x, f and h are stored
!  together as: x = xfh(1:nx), f = xfh(nxp1:jfe), h = xfh(jhs:jse), where
!  nxp1=nx+1, jfe = nx+nf, jhs = jfe+1, jhe = nx + nf + nh.
end type leaf_type
!------------------------------------------------------------------------------
type :: leaf_pointer
   type (leaf_type),  pointer :: leaf
end type leaf_pointer

!------------------------------------------------------------------------------
type :: query_type
   type (table_type), pointer :: table  ! table to which query pertains
   integer                    :: id     ! ID of leaf resolving query (or <0 if none)
   real(k_xf),  pointer       :: xs(:)	! query point (scaled)
   real(k_xf),  pointer       :: xa(:)	! projected query (scaled) point
   real(k_xf),  pointer       :: fs(:)	! f(x) (scaled)
   real                       :: cpu_start  !  CPU time at the start of the query
   logical                    :: slow_prog  !  true for slow progress
end type query_type

!------------------------------------------------------------------------------

end module isat_types