!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

!  Version 5.0

!=================================================================================

!BEGEXTRACT

recursive subroutine isatab( idtab, mode, nx, x, nf, nh, nhd, usrfgh, iusr, rusr, &
                             info, rinfo, fa, ga, ha, stats )

!  ISATAB - In Situ Adaptive TABulation
!  Fortran 90/MPI version.  S.B. Pope,  8/11/2000
!              Version 4.0  S.B. Pope, 10/04/2003
!              Version 5.0  S.B. Pope, 10/07/2006
! 
!  All rights reserved by Ithaca Combustion Enterprise, LLC.
!
!==========================================================================
!
!-------  OVERVIEW:
!
!  There is a user-supplied subroutine, usrfgh, which (on request)
!  returns the values of the user-defined function f(x), the 
!  derivatives g(x) = df(x)/dx, and (optionally) a second user-defined
!  function h(x).
!
!  x, f and h are vectors of length nx, nf and nh respectively.
!  Given x, ISATAB returns a value fa which (with high probability) is
!  within a specified error tolerance of the exact value, f(x). 
!
!  Optionally, piecewise-constant approximations to g(x) and h(x)
!  are returned, without error control.
!
!  Multiple tables:  isatab can be used to construct distinct tables
!     for different functions (i.e., for different, f, g, h and usrfgh),
!     and it can be called recursively.
!
!  All real variables are double precision.
!
!  Method: isatab implements the in situ adaptive tabulation algorithm
!     of Pope (1997), with many augmentations to improve efficiency and
!     error control.
!
!  Relation to previous versions:  version 5.0 is essentially a new code.
!     While the calling sequence to isatab is the same as in version 4.0,
!     many components of mode, info, rinfo and stats are re-defined (see below).
!
!-------  INPUT
!
!	idtab   - unique identifier of the table (idtab >= 0 )
!   mode    - integer determining the action to be taken 
!           = 0 for a "regular query call" (i.e., given x, return fa)
!           > 0 for "special calls" (described below)
!	nx	    - number of components of x  ( nx >= 1 )
!	x	    - components of x
!	nf	    - number of components of f  ( nf >= 1 )
!	nh	    - number of components of h  ( nh >= 0 )
!	nhd	    - dimension of ha ( nhd >= max(1,nh) )
!	usrfgh	- name of the user-supplied subroutine that returns
!		        f(x), g(x) = df(x)/dx, and h(x) (see below)
!	iusr	- user-defined integer array passed to usrfgh
!	rusr	- user-defined real array passed to usrfgh
!	info	- integer array controlling ISATAB operations (see below)
!	rinfo	- real array controlling ISATAB operations (see below)
!
!-------  OUTPUT
!
!	fa	- piecewise-linear approximation to f(x)
!	ga	- piecewise-constant approximation to g(x) (for info(2)=1)
!		  (ga must be dimensioned at least nf*nx)
!	ha	- piecewise-constant approximation to h(x) (for info(3)=1)
!	stats	- statistics of ISATAB performance (see below)
!
!  Note: Depending on  mode  and other settings, some of the input may not be
!  referenced, and some of the output may not be set.
!
!----- SUBROUTINE USRFGH
!
!  The subroutine 'usrfgh' must be provided.  The calling sequence is:
!
!	call usrfgh( need, nx, x, nf, nh, iusr, rusr, fa, ga, ha )
!
!  This routine (on request) returns the values of fa=f(x), ga=g(x) = df/dx
!  and ha=h(x), for the specified value of x.  The array ga contains
!  ga(i,j) = d f(i) / dx(j).  The array rusr must be dimensioned at least one, 
!  and rusr(1) must contain this dimension.  Except for need, the variables 
!  are as defined above. 
!  The integer array  need(3) specifies which functions are required:
!	 need(1) = 1 - fa is needed
!	 need(2) = 1 - ga is needed
!	 need(3) = 1 - ha is needed
!
!  Expert user option: for info(47)=de_nearby=1 (see below), f, g and h are 
!  evaluated not at x, but at a nearby location, xe.  On return from usrfgh:
!    x contains xe; rusr(1:nf) contains f(xe); fa contains an approximation to f(x). 
!
!--------  SCALING
!
!  ISATAB works in terms of the scaled variables  xs(i) = x(i) / xscale(i)
!  and  fs(j) = f(j) / fscale(j).  By default, xscale(i) and fscale(j)
!  are taken to be unity, and in general they must be strictly positive.
!  For best performance, xscale(i) and fscale(j) should be chosen so that 
!  differences in xs(i) and fs(j) are of order unity.  In particular, elements 
!  of the matrix  gs(xs) = d fs(xs) / d xs  should not be large compared to unity.
!
!  ISATAB works best if f is approximately linear in x.  The user should
!  consider whether a transformation of variables (prior to calling
!  isatab) can be used to reduce the level of non-linearity.
!
!-------  LOCAL VARIABLE NAMES
!
!  Most of the parameters input via the arrays  info  and  rinfo  are given 
!  names within the code, i.e., local variable names.  In the description given
!  here, these names (e.g., etola, etolr) are used; and they are defined where
!  the components of info and rinfo are defined.
!
!-------  ERROR TOLERANCE
!
!  The tabulation error is defined as the 2-norm:
! 	err = | fnorm |,   where
!	fnorm(j) =  [ fa(j) - f(j) ] / fscale(j).
!  The tabulation error is deemed acceptable if
!  	err < etol = etola + etolr * | f(j)/fscale(j) |,
!  where  etola = rinfo(1),  and  etolr = rinfo(2)  are the specified error
!  tolerances.  etola must be strictly positive.
!  See note below on the specification of the error tolerance.
!
!------  MPI IMPLEMENTATION
!
!  ISATAB can be used with parallel codes using MPI.  A separate
!  table is used on each processor, so that there is no message passing.
!  The principal issues are with input and output files, and with
!  synchronizing "maintenance" -- see below.
!
!------  CONCEPTS and DATA STRUCTURES
!
!  The basic storage element is a leaf.  Conceptually this consists of:
!     id   - the unique positive ID number of the leaf
!     x    - the value of x stored
!     f=f(x), g=g(x), h=h(x)
!
!     EOA  - ellipsoid of accuracy: an ellipsoid, centered at x, such that, 
!            for all points q in the EOA, the error in the linear approximation
!            fa = f + g * (q-x) is estimated to be less than the error tolerance.
!     EOI  - ellipsoid of inaccuracy: an ellipsoid, centered at x, such that, 
!            for all points q outside the EOI, the error in the linear 
!            approximation  fa = f + g * (q-x) is estimated to be greater than 
!            the error tolerance.
!
!     SELL - data structure for Sets of ELLipsoids, used for the EOAs and EOIs.
!     BT   - binary tree: defined on SELL of EOAs
!     MRU  - linked list of most-recently-used leaves (where "used" means 
!            "retrieved from")
!     MFU  - linked list of most-frequently-used leaves 
!
!     Affine space: an affine space A of dimension na (1 <= na <= nx) is determimed
!            such that the tabulation points x are close to A.  This is used
!            so that searches can be performed more efficiently in na-space
!            rather than in nx-space.
!
!     PEOA, PEOI - the projection of the EOAs and EOIs onto the affine space.
!            These are also stored as SELLs
!     EBT   - ellipsoidal binary trees defined on PEOAs and PEOIs
!
!    A leaf may become "degenerate" in which case it does not have an EOI (nor a PEOI).
!
!------  ALGORITHM
!
!  On a regular call (mode=0), there is a "query" given by the vector x,
!       and the objective is to "resolve" the query, i.e., to return fa,
!       an approximation to f(x).
!
!  The basic algorithm consists of the 5 stages described below.  
!       If  fa  is resolved in any stage, then fa is returned, and
!       the subsequent stages are not performed. 
!
! 1/ Primary Retrieve
! 2/ Secondary retrieve
! 3/ Grow
! 4/ Add or replace 
! 5/ Direct evaluation
!
!  In more detail, the stages are as follows.  Note that many of the parameters
!  set in  info  and  rinfo  control the detailed behaviors of these stages.
!
!  Primary retrieve:
!       a) The BT of EOAs is traversed to identify the primary leaf.
!          If  x  is in this EOA, the query is resolved.
!       b) The MRU is tested for an EOA containing x.
!       c) The MFU is tested for an EOA containing x.
!
!  Secondary retrieve:
!       The EBT of PEOAs is traversed in an attempt to find an EOA containing x.
!
!  Grow:
!       f(x) is evaluated.  The EBT of PEOIs is traversed and if x is in an EOI, 
!       then the EOA/EOI is "modified".  The EOA may be "grown" and/or the
!       EOI may be shrunk.
!
!  If the preceding three stages fail to resolve the query, then, if the table 
!       is not full, an "add" is performed.  Otherwise (depending on the settings)
!       either a "replace" or a "direct evaluation" (DE) is performed.
!
!  Add: g(x) and h(x) are evaluated and a new leaf is added.
!
!  Replace: the least recently used leaf is deleted, and then an add is performed
!
!  Direct evaluation: f(x) (and g and h if required) are evaluated and returned. 
!
!----- CONTROL OF STAGES / UNRESOLVED QUERIES
!  
!  There are "expert" settings (no_pr, no_sr, no_grow, no_add, no_dev, described 
!  below) which can be used to suppress one or more of the 5 stages.  
!  If direct evaluation is suppressed (by setting no_dev=1), it is possible that
!  the query will not be resolved.  Such a query is "unresolved", and, on return,
!  this outcome is indicated by stats(9)==8.  (For a regular call, stats(9) is 
!  always set so that unresolved queries can be identified.) 
!
!----- FILE NAMES
!
!  Several optional output files may be generated, and there is an optional input
!  file.  The files names are of the form  isat_#_P.ext.  Here # denotes the table
!  number (specified in idtab) and P is the MPI rank of the process.  If there is
!  a single process, or if mpi_uniq=1, then the file name is abbreviated to 
!  isat_#.ext (unless mpi_uniq=2).  The extension  ext  distinguished the different
!  files for a given table and process.  
!
!  The files involved are:
!  isat_#_P.nml   - namelist file that can be used to over-ride some 
!                   parameters specified in  info  and  rinfo.
!  isat_#_P.op    - output on the performance of ISATAB, from the array stats.
!  isat_#_P.log   - output file providing a log of the operations
!                   performed.
!  isat_#_P.tab   - data file for writing, reading and checkpointing the table.
!  isat_#_P.err   - file containing error messages (if any).
!  isat_#_P.leaf  - file containing diagnostic information on leaves.  Created by
!                   mode=13
!  isat_err_#_P.cdf - files containing output of CDF of retrieve error (see below) 
!
!-----SLOW PROGRESS
!
!  Some of the specifications of the variables  info  and  rinfo  (described below)
!  are altered if ISAT is deemed to be making "slow progress."  See the explanations
!  of rinfo(18:20).  (The values of ret_frac, grow_frac and mode_eoi are affected.)
!
!-----INFO
!
!  info	- Integer array dimensioned 100.  The contents of info(1:79) specify 
!         control parameters.  For most components of  info  the default value is 
!         0: non-zero default values are shown below in [ ].  If  info(i)  is set
!         to zero, then the default value is used.  On a regular query call 
!         (mode=0), info is referenced only on the first call (for the given table,
!         for which the value of  idtab  is denoted by #).  On special calls with
!         mode=1, the values of some control parameters are reset to the values 
!         specified in info.  Quantities that can be changed in this way are 
!         indicated below either as +info( ) or as -info( ).  (The +info( ) 
!         variables can also be set in the namelist; the -info( ) variables cannot.)
!         The components  info(80:100)  are used on some special calls (see below).
!         As decribed immediately below, the meaning of some  info(i)  are modified
!         with MPI.  The info items providing basic control are listed first, 
!         then those for expert users.
!
!  Scaling

! info( 1) = iscale    = 0, xscale(i) and fscale(j) are taken to be unity
!                      = 1, xscale(i) and fscale(j) are taken from  rinfo

!  Control of actions to be taken

!+info( 2) = if_g      = 1, if the approximation  ga  to df/dx is to be returned
!+info( 3) = if_h      = 1, if the approximation  ha  is to be returned

!  Primary retrieving

!+info( 4) = pr_mru    =  n > 0 test at most n MRU cache entries        [5]
!                      = -1     test all MRU cache entries
!                      = -2     do not use MRU cache
!                      = -3     do not maintain MRU cache
!                      
!+info( 5) = pr_mfu    =  n > 0 test at most n MFU cache entries        [30]
!                      = -1     test all MFU cache entries
!                      = -2     do not use MFU cache
!                      = -3     do not maintain MFU cache

!  Secondary retrieving

!+info( 6) = ret_min   - minimum number of EOA tests in secondary retrieving [10]
!                        (set to -1 to suppress secondary retrieving)
!+info( 7) = ret_max   - maximum number of EOA tests in secondary retrieving [-1]
!                        (set to -1 to impose no limit on secondary retrieving)

!  Growing

!+info( 8) = grow_min  - minimum number of EOA/EOIs to be grown [ 1]
!+info( 9) = grow_max  - maximum number of EOA/EOIs to be grown [-1]
!                        (set to -1 to impose no limit on growing)

!  Checkpointing

! info(10) = ichin     = 0, the ISAT table is created from scratch
!                      = 1, the initial ISAT table is to be read from the file  
!                           isat_#_P.tab
!-info(11) = ichout    = 1, to checkpoint the table occasionally on node 0
!                      = 2, to checkpoint the table occasionally on all nodes

!  Control of I/O

!-info(12) = isatop    = 1, for isat performance output on node 0
!                      = 2, for isat performance output on all nodes 
!                           (on file isat_#_P.op)
! info(13) = lu_err    >=0, logical unit number for error output;
!                      = -1, error output goes to file isat_#_P.err

!---  For expert users ---

!  Checkpointing

!+info(19) = ich0      - checkpoint table first after ich0 adds  [100] 
!                        (or, if full, grows or replaces) 
!                        (requires ichout > 0)
!-info(20) = stats0    = 0, for ichin=1, continue from  stats  read from 
!                           isat_#_P.tab
!                      = 1  for ichin=1, re-initialize  stats

!  Control of actions to be taken

!+info(21) = ifull     = 0, if table becomes full, stop adding leaves
!                      = 1, if table becomes full, replace old leaf with new  
!+info(22) = istats    = 1, return stats as well as fa

!+info(23) = no_pr     = 1  to suppress primary retrieving
!+info(24) = no_sr     = 1  to suppress secondary retrieving
!+info(25) = no_grow   = 1  to suppress growing
!+info(26) = no_add    = 1  to suppress adding and replacing
!+info(27) = no_dev    = 1  to suppress direct evaluation

!+info(28) = idites    > 0, perform accuracy test every  idites  retrieves
!                           (incurs performance penalty, unless idites is large)
!+info(29) = defer     = 1  to defer maintenance activities (affine, stats_op,
!                           checkpoint).  See explanation below.
!+info(30) = no_bt     = 1  to suppress BT in primary retrieve

!  Growing

!+info(35) = grow_mode - assumed region of accuracy used for growing:
!                      = 1, ellipsiodal; =2, Chew modification; =3, conical. [1]

!+info(36) = mode_con  - action to be taken if there is a conflict, i.e.,
!                        if the shrunk EOI does not cover the grown EOA [4]
!                      = 1 - restore original EOA; degenerate leaf
!                      = 2 - shrink EOA to be covered by EOI; degenerate leaf
!                      = 3 - shrink EOA to be covered by EOI; do not degenerate 
!                            leaf
!                      = 4 - do nothing (conflict not resolved)

!+info(37) = mode_eoi  - method used in shrinking EOI
!                      = 0 - conservative shrinking
!                      = 1 - shrink based on x if inaccurate
!                      = 2 - shrink based on err ~ r**2
!                      = 3 - no shrinking
!                        Note: in the case of "slow progess" mode_eoi is (temporarily) set to 3.

!+info(38) = mode_shrink - type of ellipsoid shrinking to be used [3]
!                      = 1 - maximun-volume algorithm       (E_v)
!                      = 2 - maximum near content algorithm (E_n)
!                      = 3 - conservative algorithm         (E_c)

!+info(39) = degen_r   - degenerate leaf after degen_r retrieves          [-1]
!                        (set degen_r = -1 to deactivate)
!+info(40) = degen_g   - degenerate leaf after degen_g grows              [-1]
!                        (set degen_g = -1 to deactivate)
!+info(41) = degen_s   - degenerate leaf after degen_s shrinks            [-1]
!                        (set degen_s = -1 to deactivate)
!+info(42) = degen_c   - degenerate leaf after degen_c conflicts          [-1]
!                        (set degen_c = -1 to deactivate)
!+info(43) = ecc_act   - error checking and control action                [1]
!                      =  1 - return exact f; add error to CDF
!                      =  2 - return approximate f; add error to CDF
!                      = -1 - return exact f; do not add error to CDF
!                      = -2 - return approximate f; do not add error to CDF
!  Adding

!+info(44) = ad_look   - maximum number of points to look for covered by 
!                        initial EOI ball                               [1000]
!+info(45) = ad_use    - maximum number of points to use to set EOI 
!                        (ad_use <= add_look)                           [100]
!+info(46) = ad_shrink = 1 to expand the initial EOI based on phi

!  Direct evaluation

!+info(47) = de_nearby = 0 - subroutine usrfgh evaluates f, g, h at the given point x
!                      = 1 - f, g, h are evaluated at a nearby point, xe

!  Affine space

!+info(48) = aff_na    = 0 - automatically determine dimension na of affine space
!                      < 0 - set na = nx
!                      > 0 - set na = aff_na < nx
!+info(49) = aff_lim   - impose upper limit on na,  na <= na_max, where  [-3]
!                      >  0 - na_max = aff_lim
!                      = -1 - na_max = nx (i.e., no limit) 
!                      = -2 - na_max = sqrt(nx)
!                      = -3 - na_max = log_2(nx) 

!  Control of EBT

!+info(50) = sr_mode   - mode used in EBT for SR   ( 1 <= sr_mode <= 8 ) [2]
!+info(51) = gr_mode   - mode used in EBT for grow ( 1 <= gr_mode <= 8 ) [2]
!+info(52) = pair_cover- mode used for ellipsoid pair covering           [2]
!                      = 1 - spheroid algorithm (no shrinking)
!                      = 2 - covariance algorithm (QL implementation)
!                      = 3 - iterative algorithm
!                      = 4 - spheroid algorithm (with shrinking)
!                      = 5 - covariance algorithm (Cholesky implementation)
!+info(53) = cut_add   - maximum number of iterations in ell_pair_separate when 
!                        forming cutting planes when adding to EBT
!                        (set to -1 to suppress updating of cutting planes of 
!                         antecedent nodes)
!+info(54) = cut_rem   - same as cut_add, but when removing
!+info(55) = cut_upd   - same as cut_add, but when updating

!  Miscellaneous

!+info(60) = icheck    = 0, 1, 2 - level of checking used in sell_m etc.
! info(61) = isat_vers - version of ISATAB (not checked for isat_vers=0)
!+info(62) = n_spool   - for n_spool > 1, spool n_spool records of stats
!+info(63) = n_add_op  - output (on add.op) x for the primary leaf and 
!                        the query on first n_add_op adds
!
! info(65) = mpi_log   = 0, for log file for each processor; =1, for log
!                          file for processor 0 only
! info(66) = mpi_uniq  = 0, use unique file names (isat_n_p.tail)
!                      = 1, use the same file name on each processor (isat_n.tail)
!
! info(70) = cdf_error > 0 to form CDF of query error 
!                        (Requires idites=info(28) > 0  or  ecc_act=info(71) > 0.)
!                        With checkpointing, cdf_error=1 re-starts from the CDF 
!                        which is read in; cdf_error=2, re-starts from scratch.
!
! info(80:100) --- for use with mode > 0
!
!-----RINFO
!
!  rinfo - Real array dimensioned  50+nx+nf.
!          The values of rinfo(1:50) are set on the first call, and can be 
!          changed only through a special call with mode=1.  
!          The components which can be changed with mode=1 are indicated
!          below by +rinfo( ).  Default values, indicated by [ ], are used if 
!          rinfo(i) is set to zero (for i=1:50).
!
! +rinfo( 1) = etola     - specified absolute error tolerance >0       [1.d-4]
!  rinfo( 2) = etolr     - specified relative error tolerance >=0
! +rinfo( 3) = etolc     - EOA is not grown if |f-f0| > etolc          [1.d20]
! +rinfo( 4) = eoa_lim   - the initial outer radius of the EOA, r0_eoa, is limited
!                          to r0_eoa < eoa_lim * sqrt(etol)            [0.1d0]
! +rinfo( 5) = eoi_lim   - the initial outer radius of the EOI, r0_eoi, is limited
!                          to r0_eoi < eoi_lim * sqrt(etol)            [10.d0]
! +rinfo( 6) = ret_frac  - perform secondary retrieving for a CPU time equivalent 
!                          to ret_frac * (CPU time of f evaluation).  This is 
!                          approximate and is limited by ret_min and ret_max.  
!                          See note below.                             [0.5]
!                          Temporarily modified in the case of slow progress.
! +rinfo( 7) = grow_frac - perform growing for a CPU time equivalent to
!                          grow_frac * (CPU time of f evaluation).  This is 
!                          approximate and is limited by grow_min and grow_max.  
!                          See note below.                             [2.0]
!                          Temporarily modified in the case of slow progress.
! +rinfo( 8) = stomby    - storage in megabytes allowed for ISAT table [500.]
!                          (for stomby < 0., max. leaves is set to |stomby|)
! +rinfo( 9) = outinc    - output increment (>=1.0)                    [1.02]
!                        = 1.0 for output on every call
! +rinfo(10) = chk_inc   - checkpoint table when the number of leaves has increased
!                          by the factor chk_inc (>=1.0)               [1.2]
! +rinfo(11) = chk_full  - when table is full, checkpoint when the number of grows
!                          (for ifull=0) or the number of replaces (for ifull=1) 
!                          has increased by the  chk_full * the number of leaves 
!                          (>=1.0)                                     [1.2]
! +rinfo(12) = aff_adds  - first set affine space after aff_adds adds  [10.]
!                          Setting aff_ads negative suppresses re-defining the 
!                          affine space.
! +rinfo(13) = aff_inc   - increment factor for aff_adds ( > 1. )       [1.2]
! +rinfo(14) = aff_thresh -threshold on rms deviation for affine space  [1.d-1]
! +rinfo(15) = aff_ratio - parameter (>1.) for re-defining affine space [2.d0]
!                        - for aff_na=0, re-define affine space if rms deviation
!                          can be reduced by at least a factor of aff_ratio
! +rinfo(16) = ad_theta  - parameter controlling EOI initialization not to cover 
!                          existing tabulation points. 0 < ad_theta <= 1.  
!                          Set negative to suppress                     [0.9]
! +rinfo(17) = pair_sep  - value of pair_sep_qual used in ell_pair_separate (EBT)
!                                                                      [0.9]
! +rinfo(18) = slow_prog_thresh - ISAT is deemed to be making slow progress
!                          if the ratio of (average CPU time per query) /
!                          (direct evaluation time) exceeds slow_prog_thresh
!                                                                      [1.0]
! +rinfo(19) = ret_frac_slow - value of ret_frac temporarily assigned to 
!                          ret_frac  in the case of slow progress      [0.5]
! +rinfo(20) = grow_frac_slow - value of grow_frac temporarily assigned to 
!                          grow_frac  in the case of slow progress     [0.5]
! +rinfo(21) = ecc_freq  - frequency of error checking and correction (ECC) [0.1]
!                           0 < ecc_freq < 1  - frequency of ECC relative to growing
!                          -1 < ecc_freq < 0  - |ecc_freq| = ECC frequency
!                               ecc_freq <= -1  to suppress ECC 
!  rinfo(22:50)          - reserved for future use and with special modes
!
! The following are referenced only on the first call and cannot subsequently 
! be changed.
!
!  rinfo(50+i)          - xscale(i) > 0, i=1,nx (required for iscale=1)
!  rinfo(50+nx+j)       - fscale(j) > 0, j=1,nf (required for iscale=1)
!
!
! Note:  Setting ret_frac >=0  or  grow_frac >=0 leads to non-deterministic 
!        performance because of the use of CPU time in determining the amount of 
!        secondary retrieving and growing.  To obtain deterministic performance, 
!        set ret_frac < 0 and grow_frac < 0, in which case  ret_limit = ret_max  
!        and  grow_limit = grow_max.
!
!-----  NAMELIST
!
!  The initial values of the indicated (+) parameters set by info and rinfo can be
!  changed by setting the value of the local variable in the namelist file  
!  isat_#_P.nml  (e.g., etola = 1.d-5).    The values set in the namelist file 
!  supercede those specified in info and rinfo (except for lu_err which cannot be
!  changed through the namelist file).  Setting  nml_op=1  in the namelist file 
!  causes the values of all quantities in the namelist to be output to the isat log
!  file.  No checking is performed on values set through the namelist, and hence
!  directly specifying  info  and  rinfo  through the arguments is preferable.
!
!-----  DEFERRED MAINTENANCE ( defer=info(29)=1 )
!
!   In normal operation (with defer=0), some maintenance operations occur
!   occasionally on regular calls.  These operations are: (1) possibly re-forming 
!   the affine space; (2) outputting stats; (3) checkpointing the table.  
!   For parallel operation it is best to synchronize these operations.  
!   This is achieved by:
!   1/ setting info(29)=defer=1
!   2/ on return (occasionally) checking whether stats(21:23) are greater than
!      0.d0, indicating that there is some maintenance to perform
!   3/ making special calls (with mode=20) to perform all maintenance, or
!      making special calls (with mode=17, 9, 12) to perform the individual
!      maintenance actions. 
! 
!-----MODE
!
!  For normal operation, set mode=0.  This is called a "regular query call".
!  Calls with mode > 0 are "special calls".  A list of special calls is now
!  given, with details to follow.
!
!  mode = 0  - regular query call (given x return fa)
!
! Change some parameters
!  mode = 1  - change the value of one or more parameters set by info and/or rinfo
!       = 2  - set the table to be "full"
!       = 3  - set the table to be "not full"
!       = 4  - set to "quick mode" 
!
! Obtain and/or write statistics (stats)
!  mode = 6  - return the statistics  stats  and output  stats 
!       = 7  - re-initialize the statistics  stats  (except stats(1)): 
!              continue output to isat_#_P.op
!       = 8  - re-initialize the statistics  stats: rewind isat_#_P.op
!       = 9  - write spooled  stats  output
!       = 10 - return in stats(95) the MFU index corresponding to a cumulative
!              fraction frac = rinfo(41) of uses 
!       = 11 - with F(x) being the CDF of retrieve errors, and for specified 
!              F=rinfo(41), return in stats(96) the value of x such that 
!              F(x)=rinfo(41).
!
! Other output
!  mode = 12 - checkpoint the table
!       = 13 - write leaf properties to isat_#_P.leaf (post-process with leaf.m)
!
! Perform action on data structures
!  mode = 15 - test consistency of isat data structures
!       = 16 - rebuild BT 
!       = 17 - consider forming new affine space
!       = 18 - rebuild the EBT on SPEOA using mode=info(81)
!       = 19 - rebuild the EBT on SPEOI using mode=info(81)
!       = 20 - perform all deferred maintenance
!       = 21 - degenerate all leaves
!       = 22 - delete the ISAT table
!       = 23 - compute ranges of x and f
!       = 24 - call usrfgh with leaf values of (x,f,g,h)
!
! Perform action on CDF of retirve error
!  mode = 25 - act on cdf_err  (the CDF of the interpolation error:
!              see below for details)
!
!  mode = 28 - return in stats the distribution of leaf usage.
!
!----  GENERAL NOTES ON SPECIAL CALLS
!
! On "special action" calls (mode > 0), fa etc. are not returned, and much of the 
! input is not referenced.  
! Initialization is performed by a regular call (mode=0): no special call is needed.
! A special call cannot be made before the table has been initialized.
! Some special calls reference info(80:100) and rinfo(40:50).  
! By defualt, a message is written to the log file when some special actions
! are taken.  Most of these messages can be suppressed by setting info(80)=1.
!
! Further details on Special Calls are given below.
!
!-----  SPECIAL CALL (mode=1)  CHANGE  INFO  AND  RINFO
!
!  If (in a mode=1 call) info(i) is set to the special value  info(i)= -12345  then
!     the existing value of info(i) is not changed.  Similarly if  
!     rinfo(i) = -12345.  Thus, only the components of info and rinfo that are to 
!     be changed need to be set different from -12345.
!
!  The user should be aware of the following when using special calls (mode=1) to 
!  change values of  rinfo: 
!     If etola is changed, then etolr is changed in proportion.
!     If etola is decreased, then all ELLs are shrunk accordingly; and,
!     if info(81) is set to 1 or 2, then all EBTs are updated (using mode=info(81)) .
!     If etolc is changed, existing EOAs are not modified, and hence may not 
!        satisfy the condition |f-f0| < etolc  for all points in the EOA.
!     If parameters related to the affine space are changed, the affine space
!        is not updated.  A separate call with mode=17 performs this update.
!
!  Normally, changes are recorded on the log file. Set info(80) =1 to suppress output.
!
!-----  SPECIAL CALL (mode=4)  QUICK MODE	
!
!   The optimal parameters depend significantly on the nature of the problem, and 
!   they can change as the computation proceeds.  The default parameters are 
!   intended as a reasonable compromise for a variety of cases.  Towards the end 
!   of a computation (and in some other circumstances) the optimal strategy is not 
!   to perform further growing and adding.  The special call with mode=4 sets the 
!   parameters appropriately for this "quick mode".  Specifically, the following 
!   are set: no_grow = 1, no_add = 1; ret_max  = -1, ret_frac = 0.5.
!
!-----  SPECIAL CALL (mode=6)  STATS
!
!   stats is returned.  If info(81)>0, stats is written; if info(81)>1 the spooled
!   output is written.  To complete the output at the end of a run, make this
!   call with info(81)=2.
!
!-----  SPECIAL CALL (mode=12)  CHECKPOINTING
!
! For ichout<=1, when using MPI the table is checkpointed only for process 0.  
! To checkpoint from all processes, set info(81)=1.
!
!-----  SPECIAL CALL (mode=13)  WRITE LEAF PROPERTIES	
!
! Properties are written only from process 0, unless info(81)=1.
! Second and subsequent calls with mode=13 overwrite isat_#_P.leaf unless info(82)=1.
! For k=info(83) > 0, output every k-th leaf.
!
!----  SPECIAL CALL (mode=15)    CHECK CONSISTENCY OF DATA STRUCTURES
!  
!   On return, stats(100) /= 0.d0 indicates inconsistency.
!   By default, output is written to the log file.  
!   Set info(81)=1 for output to the logical unit number given in info(82).
!   Set info(81)=2 to suppress output. 
!
!-----  SPECIAL CALL (mode=20)   PERFORM DEFERRED MAINTENANCE
!
!   The maintenance performed is (1) possibly re-defining the affine space
!   (2) flushing spooled output (3) checkpointing the table.  Item k is performed
!   if  table%deferred(k) > 0.d0, or info(80+k) > 0.  With MPI, the table
!   is checkpointed only on processor 0, unless info(83)=2.
!

!-----  SPECIAL CALL (mode=24)   CALL USRFGH WITH LEAF VALUES OF (x,f,g,h)
!
!  The subroutine usrfgh is called for all (or selected) leaves.  The subroutine 
!  usrfgh can then extract information from the leaf values of (x,f,g,h) and/or
!  it can modify these values.  The subroutine usrfgh used may be different from 
!  that used in a regular call: call_flag=info(80) is passed to usrfhg and can
!  be used to identify the call.  By default all leaves are treated.  For 
!  max_call=info(81)>0, at most max_call calls are made.  For leaf_stride=info(82)
!  >0, the leafs are traversed with stride leaf_stride.  See internal subroutine
!  isat_leaf_usrfgh for further information.
!
!------  SPECIAL CALL (mode=25)    CDF OF RETRIEVE ERROR
!
!  Optionally, the cumulative distribution function (CDF) of the retrive error
!  cdf_err can be formed.  This is controlled through info(70) and mode = 25.  
!  The output file generated isat_err_#_P.cdf can be post-processed using the
!  Matlab script cdf_err.m.  This option requires idites=info(28) > 0.
!
!  Actions can be performed on the CDF through special calls with mode = 25.
!  The action taken is determined by the value of info(80).
!  info(80) = 1 - initialize or re-initialize CDF's using default parameters
!           = 2 - initialize or re-initialize CDF's using parameters from info 
!                 and rinfo
!           = 3 - force output of CDF's
!           = 4 - turn CDF formation OFF
!           = 5 - turn CDF formation back ON
!
!------  SPECIAL CALL (mode=28)    DISTRIBUTION OF LEAF USAGE
!
!  On return, stats(j) equals the number of leaves which have been used no more than
!  fac^j times, where fac=rinfo(40) > 1. is specified.  stats(100) = leaves.
!
!  For mode = 25 and info(80) = 2 (i.e., initializing cdf_err with specified 
!  parameters), the parameters to be used are specified by:
!  info(81)  = nbin    - number of bins used in CDF
!  info(84)  = op_inc  - output increment (output after op_inc additions to the 
!                        CDF)
!  rinfo(41) = x_lower - lower bound on range 
!  rinfo(44) = x_upper - upper bound on range 
!
!  The default values are nbin=10000, op_inc=1000, x_lower=etola/1e3, 
!    x_upper=etola*1e3
!
!  The output files consists of two numbers per line with the following lines:
!      x_min    x_max
!      samples  wt_sum
!      xr(i)    cdf(i)
!  where: 
!       x_min  is the smallest sample value
!       x_max  is the largest sample value
!       sample  is the total number of samples
!       wt_sum  is the sum of the numerical weights of the samples
!       xr(:)  is a vector of values of x (in increasing order)
!       cdf(:) ia a vector containing the corresponding values of F( x(:) ) 
!   (The number of lines may be less than nbin+2.)
!
!-----STATS
!
!  stats(100)	- Real array containing current isat performance statistics.
!                 The Matlab script  isat.m  can be used to postprocess stats 
!                 written to isat+#_P.op

!  Query outcomes
!
!  stats( 1) - queries:    total number of queries
!  stats( 2) - p_ret:	   total number of primary retrieves
!  stats( 3) - s_ret:	   total number of secondary retrieves
!  stats( 4) - gro:        total number of queries resulting in grows
!  stats( 5) - add:        total number of queries resulting in adds
!  stats( 6) - rep:        total number of queries resulting in replaces
!  stats( 7) - dev:        total no. of queries resulting in direct evaluation
!  stats( 8) - unres:      total number of unresolved queries
!  stats( 9) - last:       last action, 2=primary retrieve, 4=grow, etc.

!  Table and EOA statistics

!  stats(12) - leaves:     number of leaves in the table
!  stats(13) - nleaves:    allowed maximum number of leaves in the table
!  stats(14) - degenerate: number of degenerate leaves (i.e., without an EOI).
!  stats(15) - eoa_max:	   largest EOA principal semi-axis
!  stats(16) - sr_eoa_in    sum of EOA inner radii
!  stats(17) - sr_eoa_out   sum of EOA outer radii
!  stats(18) - sr_eoi_in    sum of EOI inner radii
!  stats(19) - sr_eoi_out   sum of EOI outer radii

!  Maintenance deferred

!  stats(21) =1.d0 if there is deferred mainetance on affine space
!  stats(22) =1.d0 if there is deferred mainetance on stats_op
!  stats(23) =1.d0 if there is deferred mainetance on checkpointing

!  Primary retrieving

!  stats(25) = bt_travs     number of BT traverses
!  stats(26) = bt_in        number of BT traverses resulting in retrieves
!  stats(27) = bt_depth     sum of depth of BT traverses

!  stats(28) = mru_at       number of MRU attempts
!  stats(29) = mru_in       number of MRU attempts resulting in retrieves
!  stats(30) = mru_tests_in number of cache (PEOA) tests in MRU when successful

!  stats(31) = mfu_at       number of MFU attempts
!  stats(32) = mfu_in       number of MFU attempts resulting in retrieves
!  stats(33) = mfu_tests_in number of cache (PEOA) tests in MFU when successful

!  Secondary retrieving

!  stats(34) = sr_at        number of SR attempts
!  stats(35) = sr_in        number of SR attempts resulting in retrieves
!  stats(36) = sr_x_tests   number of EOA tests in SR attempts
!  stats(37) = sr_be_tests  number of BE tests in SR attempts
!  stats(38) = sr_x_max     maximum number of allowed EOA tests in SR attempts
!  stats(39) = sr_incomp    number of incomplete SR searches

!  Growing

!  stats(40) = gr_at        number of queries resulting in grow attempts
!  stats(41) = gr_eoi       number of EOIs identified in grow attempts
!  stats(42) = gr_x_tests   number of EOI tests in grow attempts
!  stats(43) = gr_a_tests   number of PEOI tests in grow attempts
!  stats(44) = gr_mod       number of EOA/EOIs modified in grow attempts
!  stats(45) = gr_inc       number of queries resulting in inclusive grows
!  stats(46) = gr_exc       number of queries resulting in exclusive grows
!  stats(47) = gr_pr_be     prob. of BE being tested
!  stats(48) = gr_mark_set  number of nodes set (or re-set) marked for updating
!  stats(49) = gr_marked    number of nodes marked for updating
!  stats(50) = gr_max       maximun number of EOA/EOI grows
!  stats(51) = gr_incomp    number of incomplete EOI searches
!  stats(52) = gr_be        number of BE tests

!  Error checking and correction (ECC) (see also stats(79)=nf_ecc)

!  stats(53) = ecc_prob     probability (per query) of ECC
!  stats(54) = ecc_shrink   number of EOA shrinks in ECC

!  Adding

!  stats(55:64) - adding and EOA/EOI initialization. (Numbers are cumulative)
!  stats(55) = r_eoa_in     inner radius of initial EOA
!  stats(56) = r_eoa_out    outer radius of initial EOA
!  stats(57) = r_eoi_in     inner radius of initial EOI
!  stats(58) = r_eoi_out    outer radius of initial EOI
!  stats(59) = r_eoi_lim    greatest distance of points used to shrink EOI
!  stats(60) = ad_npts      number of points used to shrink EOI
!  stats(61) = ad_cpt       number of cutting-plane tests in search for points
!  stats(62) = ad_int       number of "in" tests in search for points
!  stats(63) = ad_phi       value of phi in ell_pts_uncover
!  stats(64) = ad_shrink    number of EOA principle axes shrunk

!  Affine space

!  stats(65) - na           dimension of affine space
!  stats(66) - sig_star     normalized rms deviation from affine space 
!  stats(67) - sig_nap      principal value: sigma(na+1)/sigma(1) 

!  stats(68:70) pertain to pair searching since the affine space was last set

!  stats(68) - eoa_tests    number of EOA tests 
!  stats(69) - eoa_in       number of positive EOA tests 
!  stats(70) - peoa_tests   number of PEOA tests 

!  Interpolation errors

!  stats(71) - errsq:      error (squared) bewteen retrieve and f (if f evaluated)
!  stats(72) - errmax:     maximum of errsq
!  stats(73) - tolsq:      tolerance (squared) on last EOA test

!  Function evaluation statistics

!  stats(75) - sincef:      number of queries since usrfgh called for f
!  stats(76) - sinceg:      number of queries since usrfgh called for g
!  stats(77) - n_f			number of calls to usrfgh to evaluate f
!  stats(78) - n_g			number of calls to usrfgh to evaluate g
!  stats(79) - nf_ecc       number of f evaluations for ECC
!  CPU times

!  stats(81) - cpu_sec	  cpu seconds for this call (including time in usrfgh)
!  stats(82) - cpu_cum	  cumulative cpu seconds in isatab (for this table, inc.
!                           usrfgh)
!  stats(83) - cpu_out	  cumulative cpu seconds outside isatab
!  stats(84) - cpu_pr	  cumulative cpu seconds on primary retrieves
!  stats(85) - cpu_sr	  cumulative cpu seconds on secondary retrieves
!  stats(86) - cpu_gro    cumulative cpu seconds on grows (excl. time in usrfgh)
!  stats(87) - cpu_add	  cumulative cpu seconds on adds (excl. time in usrfgh)
!  stats(88) - cpu_rep	  cumulative cpu seconds on replaces (excl. time in usrfgh)
!  stats(89) - cpu_dev	  cumulative cpu seconds on direct evaluation (excl. time
!                           in usrfgh)
!  stats(90) - cpu_unres  cumulative cpu seconds on unresolved
!  stats(91) - cpu_f	  cumulative cpu seconds on evaluating f (in usrfgh)
!  stats(92) - cpu_g	  cumulative cpu seconds on evaluating g (in usrfgh)
!  stats(93) - cpu_spec	  cumulative cpu seconds on special actions (e.g., affine)
!  stats(94) - wall_sec   wall-clock time since first call to isatab
!
!  stats(95) - mfu_index  MFU index determined in last call with mode=10
!  stats(96) - x_cdf_err  retrieve error determined in last call with mode=11
!
!  CPU time detail for growing

!  stats(97)  - cpug_prior   cum. cpu secs. prior to call to isat_grow (excluding f-eval)
!  stats(98)  - cpug_search  cum. cpu secs. searching EOIs
!  stats(99)  - cpug_grow    cum. cpu secs. growing EOA/EOI
!  stats(100) - cpug_update  cum. cpu secs. updating EBTs
!
!------------------  NAME CONFLICTS
!
!  All ISATAB names begin with isat_, and hence such names should not be used
!  by the user.
!
!--------  SPECIFICATION OF ERROR TOLERANCE
!
!  The specification of the error tolerance is crucial to the accurate and
!  efficient use of ISATAB.  The following should be considered in each new
!  application of ISATAB.
!  a) Obviously, the smaller the error tolerance, the more accurately
!     fa(x) approximates f(x).
!  b) Specifying too small an error tolerance degrades performance because:
!	i)   a larger table is needed to achieve peak performance
!	ii)  more time is required to build the larger table
!	iii) for given storage (and therefore table size) more queries
!	     require direct evaluation of f(x) rather than retrieves
!  c)  The acceptable error is application dependent.  Tests should be
!      performed to determine the largest acceptable error tolerance.
!
!---------  APPROXIMATION OF DERIVATIVES
!
!  Large errors can result if values of fa(x) from ISATAB are used in
!  divided-difference approximations to derivatives such as, for example,
!	f'(x) = [ f(x+h) - f(x) ] / h
!  This is because ISATAB provides a piecewise-linear approximation to f(x).
!  If  x  and  x+h  lie in different approximation regions, the error
!  in  [ fa(x+h) - fa(x) ] / h  as an approximation to  f'(x)  can be
!  of order  etol/h.  Since h is typically small, this error can be large.
!  Instead, for accuracy and efficiency, the piecewise-constant approximation
!  provided by  ga(x)  should be used.  The error in this approximation
!  can be estimated to be of order  sqrt( etol |f''(x)| ).
!
!--------  NOTES
!
!  1/ If the additional function h(x) is not required, set nh=0, nhd=1.
!
!  2/ For idites=info(28)>=1, accuracy checking is performed.  This can be very
!        inefficient as f(x) is directly evaluated for each test.  However,
!        the performance penalty is small if  idites  is set sufficiently
!        large, specifically  idites >> cpu_f / cpu_pr.
!        
!  3/ The storage (in Mbytes) required for a given table is system-dependent,
!         because the storage required for pointers is system-dependent.
!         Hence, the storage required for a full table may be different from
!         the specified value rinfo(8)=stomby.
!
!  4/ Error output: the surest way to obtain accurate error output is to open
!         a file for writing prior to calling isatab, and to pass the logical
!         unit number in info(13).  Setting info(13)=0 results in the error
!         messages going to the standard output.  If info(13) is set negative,
!         then the error output goes to isat_#_P.err.  However, if an error 
!         occurs before this file is opened, the error message may go to unit 0
!         or to isat_ID_P.err, where ID is the value of idtab of another table. 
! 
!  5/ The routines detsv.f and dgqt.f in ell_lib are from Minpack-2.
!         See CopyrightMINPACK.txt for legal notices.     
!
!ENDEXTRACT
!
!==========================================================================
!
   use isat_init
   use isat_io
   use isat_val
   use isat_subs
   use isat_cdf
   use isat_kill

!---------------------------------------------------------------------------

   implicit none

   integer, intent(in) :: mode, idtab, nx, nf, nh, nhd, info(l_info)
   integer             :: iusr(*)

   real(k_xf), intent(in)  :: x(nx), rinfo(l_rinfo+nx+nf)
   real(k_xf), intent(out) :: fa(nf), ga(nf,nx), ha(nhd), stats(100)
   real(k_xf)              :: rusr(*)

   type (table_type), pointer, save :: first_table
   type (table_type), pointer       :: table, last_table
   type (query_type), pointer       :: query
   type (leaf_type),  pointer       :: leaf 

   real(k_xf), allocatable :: fe(:)

   real(k_xf) :: xs(nx), xe(nx),fs(nf), ft(nf), fts(nf), errsq, n_ret
   real(k_xf) :: used, mga
   real       :: cpu_sec, cpu_end, cpu_out, cpu_start, cpu_before,  &
                 cpu_after, cpu_sub, cpu_spec

   integer    :: need(3), modret, kisat, retcode, i, j, id, id_pl, &
                 nprocs, ichout, grorep, isatop, skip, n_rusr, &
				 nchanges, wall, wall_rate, wall_max, na, need_g, need_h
					
   integer, save :: ifst = 0, myrank, lu_adds = -1
   logical       :: accurate, full0, changed, f_evaluated, checkpoint, op, &
                    acc_test, call_ecc, ret_fex, err_rec
                    
!---------------------------------------------------------------------------

   interface
      subroutine usrfgh( need, nx, x, nf, nh, iusr, rusr, fa, ga, ha )
	     use isat_prec
         integer :: need(3),  nx, nf, nh, iusr(*)
	     real(k_xf) :: x(nx), rusr(*), fa(nf), ga(nf,nx), ha(*)
      end subroutine usrfgh
   end interface

!================  start of execution  ====================================

   stats(9) = 0.d0  !  indicate no action taken (modified below)

!-------  first call?

   if( ifst == 0 ) then
      if( mode > 0 ) return  !  no special actions prior to initialization
      ifst = 1
	  call isat_mpi_rank( myrank, nprocs )  

	  if( info(13) >= 0 ) then   
	     lu_err = info(13)  ! use specified logical unit for error output
	  else
	     lu_err = 0  ! use standard error: reset in isat_init if info(13) < 0
	  endif

      if( myrank == 0 ) then  !  check license on processor 0
         call isat_lic( -1, retcode )  
         if( retcode /= coderet ) call isat_abort('isatab',0, mess='Invalid license code' )
      endif

!  first_table is an empty table, created on the first call, and never used or killed      
      allocate( first_table )
      nullify( first_table%next_table )
      first_table%idtab = -huge(1)
   endif

!-------  search through tables for correct one, or create if non-existent

   table => first_table
   do
      if( table%idtab == idtab ) then
         exit
      elseif( associated(table%next_table) ) then 
         table => table%next_table
         cycle
      endif

!  table not found: check idtab and generate new table

      if( idtab < 0 ) call isat_abort('isatab', 1, isv = idtab, &
           mess = 'idtab must be positive: idtab =' )

      if( mode > 0 ) return !  no special actions prior to initialization of table

!  initialize new table
      last_table=> table  !  currently last table
      call isat_table_init(idtab, nx, nf, nh, info, rinfo, table )
      last_table%next_table => table  !  new table is successor of previously last table
      
	  if( table%icheck >= 2 ) then
	     call isat_integrity( table, table%lu_log, retcode ) 
	     if( retcode /= 0 ) call isat_abort('isatab', 2, isv = retcode, &
              mess = 'inconsistent data structures, info = ' )
      endif


   end do

!-------  table found or created -----------------

   ichout = table%my_ichout  
   lu_err = table%lu_err

!================  start of ISAT operation: exit this loop to return ===========

   isat_operation: do    ! one loop only

   call cpu_time( cpu_start )          ! ISATAB CPU timing starts here
   cpu_out  = cpu_start - table%cpu_end ! CPU out of ISAT since last ISAT call
   cpu_sub  = 0.
   cpu_spec = 0.

!===============  treat special calls ======================================

   if( mode /= 0 ) then
      call isatab_special
	  if( mode == 22 ) return  !  table has been deleted
	  if( mode == 28 ) return  !  do not corrupt stats
	  kisat    = -1  
	  exit isat_operation
   endif

   query    => table%query  ! data structure for query information
   query%id = -1
   query%cpu_start = cpu_start
   query%slow_prog = .false.
   
   need_g   =  table%if_g
   need_h   =  table%if_h   
   f_evaluated = .false.


!--------------  standard ISAT query -----------------------------------------

   table%stats( 1) = table%stats( 1) + 1.d0  ! increment queries
   table%stats(75) = table%stats(75) + 1.d0  ! increment sincef
   table%stats(76) = table%stats(76) + 1.d0  ! increment sinceg

!------------  form xs =  scaled x

   if( table%iscale == 0 ) then
      xs = x
   else
      xs = x * table%xsci
   endif
   query%xs = xs

!------------ mapping of xs onto affine space formed in isat_primary_retrieve if
!             retrieve is unsuccessful
 
   na = table%na

!===========  Stages 1 and 2: primary and secondary retrieving  ==================

   if( table%no_pr == 0 ) then
      call isat_primary_retrieve( query )
      id_pl = abs( query%id ) 
      if( query%id > 0 ) kisat = 2       ! outcome = primary retrieve

   else  ! set  xa  which otherwise is set in isat_primary_retrieve
      query%xa(1:na) = matmul( query%xs-table%ua(:,nx+1) , table%ua(:,1:na) )  
   endif
   
!  determine whether progress is deemed to be slow
   if( table%stats(82) > 10.d0  .and.  table%stats(91) > 1.d0 ) then
      if( table%stats(82)*table%stats(77) >  &
          table%stats(91)*table%stats(1)*table%slow_prog_thresh ) &
             query%slow_prog = .true.
   endif

   if( query%id < 0  .and.  table%no_sr == 0 ) then  ! attempt seconday retrieve
	     call isat_secondary_retrieve( query )
         if( query%id > 0 ) kisat = 3     ! outcome = secondary retrieve
   endif

   if( query%id > 0 ) then   !------------ successful retrieve

   	  need = (/ 1, need_g, need_h /)  !  set  fs and fa, and g and h as needed
	  call isat_get_fgh( query, 1, need, fa, ga, ha )

!----------  perform accuracy testing and/or error checking and correction
!  set: acc_test (perform accuracy testing); call_ecc (perform ECC); 
!       ret_fex (return exact f); err_rec (record error)

      acc_test = .false.
      
      if( table%idites >0 ) then  !  perform accuracy test every  idtest  retrieves
	     n_ret = sum( table%stats(2:3) ) + 1.d0 ! total number of retrieves (PR+SR)
	     if( nint( mod( n_ret, 1.d0*table%idites ) ) == 0 ) then
	        acc_test = .true.
	        err_rec  = .true.
	        call_ecc = .false.
	        ret_fex  = .false.
	     endif
	  endif
	      
      if( table%stats(1) >= table%ecc_q  .and. &
         (table%ecc_freq < 0.d0  .or.  table%no_grow == 0) ) then  !  perform ECC ?
         acc_test = .true.
         call_ecc = .true.
         err_rec  = .false.
	     ret_fex  = .false.
         if( table%ecc_act > 0.d0 )    err_rec = .true.
         if( abs(table%ecc_act) == 1 ) ret_fex = .true.
      endif
         
      if( acc_test ) then      !  perform accuracy testing         
         call isatab_accuracy  !  evaluate:  ft, errsq, accurate
         if( call_ecc ) call isat_ecc( query, accurate, errsq ) 
         if( ret_fex  ) fa(1:nf) = ft(1:nf) 
         
         if( err_rec ) then
            table%stats(71) = errsq   !  square of retrieve error
            table%stats(72) = max( errsq, table%stats(72) )
            if( table%cdf_err%on )  &  !  form CDF of retrieve error
		        call isat_cdf_add1( table%cdf_err, sqrt( errsq ) )
         endif
      endif   
      
!-----------  end of accuracy testing:  complete retrieve and exit

      exit isat_operation
   endif

!===== Determine whether query will be unresolved: if not, evaluate f(x)  =========

   if( table%no_grow ==1  .and.  table%no_dev ==1   .and.  &
       (table%no_add ==1  .or.  (table%full  .and.  table%ifull==0) ) ) then
      kisat = 8   !  outcome = unresolved
      exit isat_operation
   endif

! set fe: on entry to usrfgh, fe contains rusr; 
! for de_nearby=1, on exit fe contains fe /approx f(xe)

   call isatab_set_fe

   need = (/ 1, 0, 0 /)  !  evaluate f
   xe   = x              !  location at which f is evaluated

   call cpu_time( cpu_before )
   call usrfgh( need, nx, xe, nf, nh, iusr, fe, fa, ga, ha ) 
   call cpu_time( cpu_after )
   table%stats(75) = 0.d0                      ! reset sincef to zero
   table%stats(77) = table%stats(77) + 1.d0                    ! increment n_f
   table%stats(91) = table%stats(91) + (cpu_after-cpu_before)  ! increment cpu_f
   cpu_sub         = cpu_after-cpu_before      ! charge time to cpu_f, not to query
   query%cpu_start = query%cpu_start + cpu_sub ! don't count f-eval time to prior

   f_evaluated = .true.
   if( table%de_nearby == 0 ) then
      fs          = fa * table%fsci
      query%fs    = fs
      
   else  ! re-define point for growing and adding as xe (nearby to x)
      xs             = xe * table%xsci
      query%xs       = xs
      query%xa(1:na) = matmul( query%xs-table%ua(:,nx+1) , table%ua(:,1:na) )  
      fs             = fe(1:nf) * table%fsci
      query%fs       = fs
   endif

!=========  Stage 3: growing ======================================================

   if( table%no_grow == 0 ) then

      call isat_grow( query )

	  nullify(leaf)
	  if( query%id > 0 ) then  !  successful inclusive grow
	     leaf => table%leaf_pt(query%id)%leaf
         kisat  = 4   !   outcome = grow

	     need   = (/ 0, need_g, need_h /)   !  get ga and ha if required
	     call isat_get_fgh( query, 0, need, fa, ga, ha )

!  checkpoint if required when table is full (based on grows and replaces)

         if( table%full .and. ichout==1 ) then
	         grorep = table%stats(44) + table%stats(6)
		     if( grorep > table%nextch ) then

			    if( table%defer == 0.d0 ) then
	            call isat_table_write( table )
				else
                   table%deferred(3) = 1.d0
				endif

                table%nextch = grorep + table%chk_full * table%leaves + 1.1
             endif
         endif

         exit isat_operation

      endif
   endif

!===========  Stage 4: adding  ===================================================

   if( table%no_add == 0  .and.  (.not.table%full  .or.  table%ifull==1) ) then
      full0 = table%full  !  indicates value of full before adding
      call isatab_set_fe
      need  = 1
      xe    = x
   
      call cpu_time( cpu_before )
      call usrfgh( need, nx, xe, nf, nh, iusr, fe, fa, ga, ha )
      call cpu_time( cpu_after )
      table%stats(75) = 0.d0                              ! reset sincef to zero
      table%stats(76) = 0.d0                              ! reset sinceg to zero
      table%stats(78) = table%stats(78) + 1.d0                   ! increment n_g
      table%stats(92) = table%stats(92) + (cpu_after-cpu_before) ! increment cpu_fgh
! charge time to cpu_fgh, not to query
      cpu_sub = cpu_sub + (cpu_after-cpu_before) 

      f_evaluated = .true.
      
      if( table%de_nearby == 0 ) then
         fs          = fa * table%fsci
         query%fs    = fs
      
      else  ! re-define point for growing and adding as xe (nearby to x)
         xs             = xe * table%xsci
         query%xs       = xs
         query%xa(1:na) = matmul( query%xs-table%ua(:,nx+1) , table%ua(:,1:na) )  
         fs             = fe(1:nf) * table%fsci
         query%fs       = fs
      endif

!XXX detect and treat failure to obtain sensitivities

     mga = max( maxval(ga), -minval(ga) )

     if( (mga /= 0  .and.  mga == 2*mga) .or. mga > 1.d50  ) then
!detect infinity
        if( need_g == 0  .and.  need_h == 0 ) then
           kisat = 7
        else
           kisat = 8
        endif

        write(0,*)'bad sensitivities detected, kisat = ', kisat
        write(0,*)'mga =', mga
        write(0,*)'need =', need
        write(0,*)'nx, nf, nh, iusr = ', nx, nf, nh, iusr(1)
        write(0,*)'xe =', xe
        write(0,*)'fe =', fe
        write(0,*)'fa =', fa
        write(0,*)'ga =', ga

        exit isat_operation ! skip adding - VH - 05/29/2011
     endif
!XXX

      call isat_add( table, xs, fs, ga, ha, modret )
      kisat = 4 + modret  !  outcome = add (5) or replace (6)

!--------- diagnostic output: x(primary leaf), x(query) --------------------

      if( table%n_add_op > 0 ) then
         if( table%leaves <= table%n_add_op .and. id_pl > 2 ) then
         
            if( lu_adds < 0 ) then
               call isat_lu( lu_adds )
               open( lu_adds, file = 'adds.op', form = 'formatted' )
            endif
                
            write(lu_adds,'(1p,200e26.17)') table%leaf_pt(id_pl)%leaf%xfh(1:nx) * table%xscale, x  
            if( table%leaves == table%n_add_op) close(lu_adds)
         endif
      endif

!----------   possibly re-define affine space

      if( table%aff_adds > 0.d0  .and.   &
	      table%stats(5) + table%stats(6) > table%aff_adds ) then 

	     if( table%defer == 0 ) then
		    call cpu_time( cpu_before ) 
            call isatab_affine
            call cpu_time( cpu_after )
			cpu_sub = cpu_sub +   (cpu_after-cpu_before) ! charge time to cpu_spec
			cpu_spec = cpu_spec + (cpu_after-cpu_before)
		 else
		    table%deferred(1) = 1.d0
		 endif

         table%aff_adds = (table%stats(5) + table%stats(6)) * table%aff_inc
      endif

!---------  checkpointing

      if( ichout == 1 ) then  !  checkpoint table if required
	     checkpoint = .false.

         if( .not.table%full ) then  !  decide based on leaves when not full
            if( table%leaves >= table%nextch ) then
               checkpoint = .true.
               table%nextch = max(1.d0*table%leaves,table%nextch) * table%chk_inc + 1.1
            endif

	     elseif( .not.full0 ) then  !  table has just filled
               checkpoint = .true.
            grorep = table%stats(44) + table%stats(6)
            table%nextch = grorep + table%chk_full * table%leaves + 1.1d0

	     else  !  full,  base on grows plus replaces
            grorep = table%stats(44) + table%stats(6)
		    if( grorep > table%nextch ) then
               checkpoint = .true.
               table%nextch = grorep + table%chk_full * table%leaves + 1.1d0
            endif

         endif

		 if( checkpoint ) then
		    if( table%deferred(3) == 0.d0 ) then
		       call cpu_time( cpu_before ) 
               call isat_table_write( table )
			   call cpu_time( cpu_after )
			   cpu_sub = cpu_sub +   (cpu_after-cpu_before) !charge time to cpu_spec
			   cpu_spec = cpu_spec + (cpu_after-cpu_before)
			else
               table%deferred(3) = 1.d0
			endif
		 endif

      endif

      exit isat_operation
   endif

!===========  Stage 5: Direct evaluation ==========================================
!----------  grow failed:  if table full and ifull=0, return unless dfdx needed

   if( f_evaluated  .and.  need_g ==0  .and.  need_h == 0 ) then
!   f has been evaluated, and g and h are not required:  all done
      kisat = 7   !  outcome = direct evaluation

   elseif( table%no_dev == 0 ) then  !  need to evaluate g and h
      need = (/ 1, need_g, need_h /)  !  obtain g and h 
      if( f_evaluated ) need(1) = 0
      call isatab_set_fe
      xe = x

      call cpu_time( cpu_before )
      call usrfgh( need, nx, xe, nf, nh, iusr, fe, fa, ga, ha )
      call cpu_time( cpu_after )

	  if( need_g == 1 ) then
	     table%stats(76) = 0.d0                      ! reset sinceg to zero
	     table%stats(78) = table%stats(78) + 1.d0    ! increment n_g
	  endif

      table%stats(92) = table%stats(92) + (cpu_after-cpu_before) ! charge to cpu_g
      cpu_sub = cpu_sub + (cpu_after-cpu_before) !charge time to cpu_fgh, not to query
      kisat = 7   !  outcome = direct evaluation

   else
      kisat = 8   !  outcome = unresolved
   endif

   exit isat_operation

   end do isat_operation  !========================================================

   if( kisat > 0  .and.  kisat <=4 ) then !  update MRU/MFU based on leaf used
      table%seoa%ell_pt(query%id)%ell%used  = table%seoa%ell_pt(query%id)%ell%used  + 1.d0
	  table%leaf_pt(query%id)%leaf%props(7) = table%leaf_pt(query%id)%leaf%props(7) + 1.d0 

      if( table%pr_mru /= -3 ) call ll_mru_update( table%eoa_mru, query%id )
      if( table%pr_mfu /= -3 ) call ll_mfu_update( table%eoa_mfu, table%seoa, query%id )

   elseif( kisat == 8 ) then  !  non-retreive
      table%stats( 1) = table%stats( 1) - 1.d0  ! do not count as query
      table%stats(75) = table%stats(75) - 1.d0  ! decrement sincef
      table%stats(76) = table%stats(76) - 1.d0  ! decrement sinceg
	  fa(1) = -huge(1.e0)                       ! flag that fa is not set
   endif

!----------------  deal with stats and return  ------------------------------------

   call cpu_time( cpu_end )
   table%cpu_end   = cpu_end
   cpu_sec         = cpu_end - cpu_start

   table%stats(81) = cpu_sec                    ! CPU in ISAT this call
   table%stats(82) = table%stats(82) + cpu_sec  ! cum. CPU in ISAT
   table%stats(83) = table%stats(83) + cpu_out  ! cum. CPU out of ISAT
   table%stats(93) = table%stats(93) + cpu_spec ! cum. CPU on specials actions

   stats(9)       = kisat  !  outcome this call (always return in stats)
   table%stats(9) = kisat  

   stats(21:23)   = table%deferred
   table%stats(21:23) = stats(21:23)

   if( kisat <= 0 ) then
      table%stats(93) = table%stats(93) + cpu_sec ! charge CPU to specials actions
      return  !  return from special call

   elseif( kisat > 0 ) then
      table%stats(kisat)    = table%stats(kisat) + 1.d0 !total number of this event
!     cum. CPU on this event
      table%stats(82+kisat) = table%stats(82+kisat) + (cpu_sec - cpu_sub) 
   endif

   isatop = 0
   if( table%my_isatop==1 .and.  &
      (table%stats(1) >= table%nextop  .or.  table%i_hist < 0)  ) isatop = 1

   if( table%istats == 1  .or.  isatop == 1 ) then
      table%stats(12)  = table%leaves   ! number of leaves
      if( table%istats == 1 ) stats = table%stats

      if( isatop == 1 ) then
	     call isatab_stats_op  !   output this step
	     table%nextop = table%stats(1) * table%outinc
	  endif

   endif

!Laniu deallocate the fe array
    if(allocated(fe))  deallocate(fe)
   
   return

contains  !========== internal subroutines ========================================

subroutine isatab_special
   implicit none

   op = .false.  !  output to log file?
   if( table%write_log  .and.  info(80) == 0 ) op = .true.

   if( mode == 1 ) then  !  change some parameters
      call isat_change( table, info, rinfo, nchanges )
      if( op ) write(table%lu_log,'(a,i4)') &
		    'isatab: mode=1, number of parameters changed = ', nchanges

   elseif( mode == 2 ) then   !  set table to be full
      if( op ) write(table%lu_log,'(a,1pe20.10)') &
		      'Setting table full (mode=2 call) at query = ', table%stats(1)
	  if( table%full ) return ! already full
	  table%full = .true.
	  if(  ichout==1 ) then
         call isat_table_write( table ) 
!          future checkpointing based on grows and replaces
	     grorep = table%stats(44) + table%stats(6)
         table%nextch = grorep + table%chk_full * table%leaves + 1.1
	  endif

   elseif( mode == 3 ) then    !  set table to be not full
      if( op ) write(table%lu_log,'(a,1pe20.10)') &
		       'Setting table not full (mode=3 call) at query = ', table%stats(1)
	  if( table%leaves < table%maxleaves ) table%full = .false.

   elseif( mode == 4 ) then    !  quick mode
	  table%ret_max  = -1
      table%no_grow  = 1
	  table%no_add   = 1
	  table%ret_frac = 0.5d0
	  if( op ) write(table%lu_log,'(a,1pe20.10)') &
		       'Going to quick mode at query = ', table%stats(1)

   elseif( mode == 6 ) then    !  set stats
      stats     = table%stats
	  stats(9)  = 0 ! (kisat) null call
	  stats(12) = table%leaves
          ! VH: 08/3/2011, output stats only is my_isatop == 1
          if( table%my_isatop==1 ) then
             if( info(81) > 0 ) call isatab_stats_op
             if( info(81) > 1 ) call isatab_spool_op
          endif

   elseif( mode == 7 ) then    !  re-set stats
      if( op ) write(table%lu_log,'(a,1pe20.10)') &
		       'Resetting stats (mode=7 call) at query = ', table%stats(1)
	  call isat_stats_init( table, 1 )

   elseif( mode == 8 ) then    !  re-set stats
      if( op ) write(table%lu_log,'(a,1pe20.10)') &
		       'Resetting stats (mode=8 call) at query = ', table%stats(1)
	  call isat_stats_init( table, 2 )

   elseif( mode == 9 ) then
      call isatab_spool_op   !  flush spooled output

   elseif( mode == 10 ) then  !  return index from MFU
	  used = rinfo(41)*sum( table%stats(2:4) )
	  call ll_mfu_cum( table%eoa_mfu, table%seoa, used, j )
	  stats(95)       = j
	  table%stats(95) = stats(95)

   elseif( mode == 11 ) then  !  return error percentile
	  call isat_cdf_f2x( table%cdf_err, rinfo(41), stats(96) )
	  table%stats(96) = stats(96)

   elseif( mode == 12 ) then   !  checkpoint table
      if( ichout == 1  .or.  myrank == 0  .or.  info(81) == 1 ) &
		    call isat_table_write( table )
	  table%deferred(3) = 0.d0

   elseif( mode == 13 ) then   !  write leaf properties
	  if(  myrank == 0  .or.  info(81) == 1 ) then

         if( info(82) == 1  .and.  table%leaf_everon ) then  ! open file
            open( table%lu_leaf , file = table%isat_leaf  , &
			         status = 'old', action = 'write', position = 'append' )
		 else
            open( table%lu_leaf , file = table%isat_leaf  , &
			      status = 'replace', action = 'write' )
		 endif

         skip = 1
         if( info(83) > 0 ) skip = info(83)
	     call isat_leaf_props
	     close( table%lu_leaf ) 
		 table%leaf_everon = .true.

         if( op ) write(table%lu_log,'(a,1pe20.10,i10)') &
		          'Leaves written: query, leaves = ', table%stats(1), table%leaves
      endif

   elseif( mode == 15 ) then  !  check consistency of data structures
      if( info(81) == 1 ) then  !  set lu
		 lu = info(82)
	   elseif( info(81) == 2 ) then
		 lu = -1
	   else
		  lu = table%lu_log
	   endif

       call isat_integrity( table, lu, retcode )
	   stats(100) = retcode

   elseif( mode == 16 ) then  !  rebuild eoa_bt
      call bt_rebuild( table%eoa_bt, table%idlist, 1 )
      if( op ) write(table%lu_log,'(a,1pe20.10,i10)') &
		       'BT rebuilt: query, leaves = ', table%stats(1), table%leaves

   elseif( mode == 17 ) then  !  possibly re-set affine space
      call isatab_affine
      if( changed  .and.  op ) write(table%lu_log,'(a,1pe20.10,i10)') &
		       'Affine space redefined: query, na = ', table%stats(1), table%na

   elseif( mode == 18 ) then  !  rebuild EBT on PEOA with mode=info(81)
	  if( info(81) < 1  .or.  info(81) > 2 )  &
		 call isat_abort('isatab', 3, isv = info(81), &
            mess = 'mode = 18 invoked with invalid info(81) = ' )

      call ebt_rebuild( table%peoa_ebt, info(81), table%pair_cover )
      if( op ) write(table%lu_log,'(a,1pe20.10,i10)') &
		             'EBT of PEOAs rebuilt: query, PEOAs = ', table%stats(1),  &
		              table%speoa%n_ell

   elseif( mode == 19 ) then  !  rebuild EBT on PEOI with mode=info(81)
	  if( info(81) < 1  .or.  info(81) > 2 )  &
		 call isat_abort('isatab', 4, isv = info(81), &
            mess = 'mode = 19 invoked with invalid info(81) = ' )

      call ebt_rebuild( table%peoi_ebt, info(81), table%pair_cover )
      if( op ) write(table%lu_log,'(a,1pe20.10,i10)') &
		             'EBT of PEOIs rebuilt: query, PEOIs = ', table%stats(1),  &
		              table%speoi%n_ell

   elseif( mode == 20 ) then  ! all deferred maintenance
      if( table%deferred(1) > 0.d0  .or.  info(81) == 1 ) then
         call isatab_affine
         if( changed  .and.  op ) write(table%lu_log,'(a,1pe20.10,i10)') &
		       'Affine space redefined: query, na = ', table%stats(1), table%na
	  endif

	  if( table%deferred(2) > 0.d0  .or.  info(82) == 1 )   &
	     call isatab_spool_op   !  flush spooled output

      if( table%deferred(3) > 0.d0  .or.  info(83) >= 1 ) then
         if( ichout == 1  .or.  myrank == 0  .or.  info(83) == 2 ) &
		    call isat_table_write( table )
	     table%deferred(3) = 0.d0
	  endif

	  if( op ) write(table%lu_log,'(a,1pe20.10)') &
		         'Deferred maintenance performed at query = ', table%stats(1)

   elseif( mode == 21 ) then  !  degenerate all leaves
      i = table%speoi%n_ell
	  do id = 1, table%seoi%n_pt
         if( associated( table%seoi%ell_pt(id)%ell ) )  &
			call isat_leaf_degenerate( table, id )
	   end do
      if( op ) write(table%lu_log,'(a,1pe20.10,i10)') &
		             'Degenerated all leaves: query, PEOIs = ', table%stats(1), i

   elseif( mode == 22  ) then   !  delete table
      call isat_table_kill( table, first_table )
      
   elseif( mode == 23 ) then  !  compute the ranges of x and f
      if( info(81) == 1 ) then  !  set lu
		 lu = info(82)
	   elseif( info(81) == 2 ) then
		 lu = -1
	   else
		  lu = table%lu_log
	   endif

       call isat_xf_range( table, table%nx, table%nf, lu )
       
   elseif( mode == 24  ) then   !  call isat_leaf_usrfgh
      call isat_leaf_usrfgh
      
   elseif( mode == 25 ) then  !  act of CDF
   
	  call isat_cdf_act( table, mode, info, rinfo )

   elseif( mode == 28 ) then
      call isatab_mfu_dist
      
   else
	  if( op ) write(table%lu_log,'(a,i8)') &
		             'isatab: Warning, called with invalid mode = ', mode
   endif

   return
end subroutine isatab_special  !---------------------------------------------------

subroutine isatab_set_fe

! if necessary, allocated the array fe, and then set fe(1:n_rusr) = rusr

   if( .not.allocated(fe) ) then
      n_rusr = nint( rusr(1) )
      if( n_rusr < 1 ) call isat_abort('isatab', 5, isv = n_rusr, &
            mess = 'rusr(1) must be set (>=1) to the length of rusr ' )
      allocate( fe( max(nf,n_rusr) ) )
   endif
   
   fe(1:n_rusr) = rusr(1:n_rusr)
   return
   
end subroutine isatab_set_fe   !---------------------------------------------------

subroutine isatab_accuracy

   need = (/ 1, 0, 0 /)  !  evaluate f
   xe   = x

   call cpu_time( cpu_before )
   call isatab_set_fe
   call usrfgh( need, nx, xe, nf, nh, iusr, fe, ft, ga, ha )
   call cpu_time( cpu_after )
   table%stats(91) = table%stats(91) + (cpu_after-cpu_before)  ! increment cpu_f
   table%stats(77) = table%stats(77) + 1.d0                    ! increment n_f
   cpu_sub = cpu_sub + (cpu_after-cpu_before) ! charge time to cpu_f, not to query

   fts  =  ft * table%fsci
   leaf => table%leaf_pt(query%id)%leaf

   call isat_eoa_err( leaf, query%fs, fts, accurate, errsq )

   return
end subroutine isatab_accuracy  !--------------------------------------------------

subroutine isatab_affine
   implicit none

   call isat_affine( table, table%seoa, table%nx, changed ) 

   if( changed ) then  ! rebuild SPEOA and SPEOI
      call isat_aff_rebuild( table, table%seoa, table%speoa, table%peoa_ebt )
      call isat_aff_rebuild( table, table%seoi, table%speoi, table%peoi_ebt )

	  if( table%icheck >= 2 ) then
	     call isat_integrity( table, table%lu_log, retcode ) 
	     if( retcode /= 0 ) call isat_abort('isatab', 6, isv = retcode, &
              mess = 'inconsistent data structures, info = ' )
      endif

   endif

   table%deferred(1) = 0.d0

   return
end subroutine isatab_affine  !----------------------------------------------------

subroutine isatab_stats_op
   implicit none

!  obtain wall clock time
   call system_clock( count = wall, count_rate = wall_rate, count_max = wall_max )  

   if( wall < table%wall_last ) table%wall_cycle = table%wall_cycle + 1
   table%wall_last = wall

!  wall clock time in seconds
   table%stats(94) = (float( wall - table%wall_0) +  &
                     table%wall_cycle * float(wall_max) ) / float(wall_rate)  


   if( table%i_hist < 0 ) then !  stats_hist to be output
      table%i_hist = -table%i_hist
	  do i = 1, table%i_hist
	     write(table%lu_op,"(1p,100e15.7)") table%stats_hist(1:100,i)
	  end do

   elseif( table%i_hist == n_hist ) then !  compress stats_hist if necessary
       do i = 1, n_hist, 2
	      j = (i+1)/2
		  table%stats_hist(:,j) = table%stats_hist(:,i)
	   end do
       table%i_hist = j
   endif

   table%i_hist = table%i_hist + 1  !  record stats
   table%stats_hist(:,table%i_hist) = table%stats

   if( table%n_spool <= 1 ) then  !  write statistics to lu_op 
      write(table%lu_op,"(1p,100e15.7)") table%stats(1:100)  
      call isat_flush( table%lu_op )

   else                           !  spool output
      if( table%i_spool == table%n_spool ) call isatab_spool_op  ! spool is full
      table%i_spool     = table%i_spool + 1
	  stats(22)         = table%i_spool
	  table%stats(22)   = stats(22)
	  table%deferred(2) = stats(22)
	  table%spool(1:100,table%i_spool) = table%stats(1:100)
   endif

   return
end subroutine isatab_stats_op  !--------------------------------------------------

subroutine isatab_spool_op
   implicit none

   if( table%i_spool < 1 ) return

!  write statistics to lu_op 
   write(table%lu_op,"((1p,100e15.7))") table%spool(1:100,1:table%i_spool)  
   call isat_flush( table%lu_op )
   table%i_spool     = 0
   table%deferred(2) = 0.d0

end subroutine isatab_spool_op !---------------------------------------------------

subroutine isatab_mfu_dist

! return in stats the "use" distribution of leaves
   implicit none
   real(k_xf) :: fac, thresh, used, used_tot
   
   fac = rinfo(40)
   if( fac < 1.d0 ) call isat_abort('isatab_mfu_dist',7,  &
                    mess='rinfo(40) must be greater than 1.', rsv=fac )
   thresh   = 1.5d0
   j        = 1
   id       = table%eoa_mfu%last
   stats    = table%leaves
   stats(1) = 0.d0
   used_tot = 0.d0
   
   leaf_loop: do 
      used = table%seoa%ell_pt(id)%ell%used + 1.d0  ! +1 to include initial add
      used_tot = used_tot + used
    
      do while (used >= thresh )  ! as needed, increase j and until thresh > used
         j = j + 1
         if( j == 100 ) exit leaf_loop
         stats(j) = stats(j-1)
         if( thresh > 1.d8 ) then
            thresh = fac * thresh
         else
            thresh   = max( thresh+1.d0, 0.5d0+nint(fac*thresh) )
         endif
      end do
      
      stats(j) = stats(j) + 1.d0
      if( id == table%eoa_mfu%first ) exit
      id = table%eoa_mfu%before(id)
   
   end do leaf_loop
      
   return

end subroutine isatab_mfu_dist

subroutine isat_leaf_props

!  Write to table%lu_leaf the properties of leaves, 
!  skipping skip leaves

!  Definition of leaf properties, leaf%props(k)

! k=1  query on which leaf added to table
! k=2  query of last retrieve from leaf (or 0 if none)
! k=3  query of last EOA grow (or 0 if none)
! k=4  query of last EOI shrink (or 0 if none)
! k=5  query of last EOA/EOI conflict (or 0 if none)
! k=6  query on which leaf was degenerated (or 0 if not)

! k=7  number of retrieves (and inc. grows) from leaf
! k=8  number of EOA grows
! k=9  number of EOI shrinks
! k=10 number of EOA/EOI conflicts

! k=11 PEOA r_in
! k=12 PEOA r_out
! k=13 PEOI r_in
! k=14 PEOI r_out

   type (ell_type), pointer   :: pell
   integer                    :: id, n_leaves
   real(k_xf)                 :: rad(4)

   if( table%lu_leaf < 1  .or.  table%leaves < 1 ) return

   n_leaves = 0
   do id = 1, table%n_pt
      if( .not.associated( table%leaf_pt(id)%leaf ) ) cycle
	  n_leaves = n_leaves + 1
	  if( mod(n_leaves,skip) /= 0 ) cycle

	  pell => table%speoa%ell_pt(id)%ell  ! PEOA
	  rad(1:2) = pell%geom( pell%sell%ngeom-1:pell%sell%ngeom )

      if( associated( table%speoi%ell_pt(id)%ell ) ) then
	     pell => table%speoi%ell_pt(id)%ell  ! PEOI
	     rad(3:4) = pell%geom( pell%sell%ngeom-1:pell%sell%ngeom )
      else
	     rad(3:4) = 0.d0
	  endif

	  write(table%lu_leaf,"((1p,200e15.7))") table%leaf_pt(id)%leaf%props,  &
	                                         rad, 1.d0*id
   end do

   call isat_flush( table%lu_leaf )

   return
end subroutine isat_leaf_props  !------------------------------------------------

subroutine isat_leaf_usrfgh

!  Loop over leaves and make the following call to usrfgh:
!	  call usrfgh( replace, nx, x, nf, nh, i_wrk, r_wrk, fa, ga, ha )
!
!  call_flag   = info(80) - flag passed to usrfgh in i_wrk(1)
!  max_call    = info(81) - maximum number of calls to be made
!  leaf_stride = info(82) - stride for loop over leaves
!
!  Note that fa and ga are scaled: xscale and fscale are passed to usrfgh in r_wrk.
!  The values of x, fa and ha (stored in the leaf) can be changed in usrfgh.
!  On return, if replace(2)==1, then the leaf value of g is replaced by the
!  returned value of ga.
!
!  Warning: on this special call, input arguments nx, nf, etc. to isatab may not be set correctly:
!  use table%nx, table%nf, etc.
!
!  On return: stats(1) = no. of leaves,  stats(2) = no. of calls to usrfgh

   type (leaf_type), pointer  :: leaf
   integer                    :: call_flag, max_call, leaf_stride, id, n_leaves, &
                                 i_wrk(3), replace(3), nnx, nnf, nxp1, jfe, jhs, &
                                 jhe
   real(k_xf)                 :: r_wrk(table%nx+table%nf+12)
   
   stats(1) = table%leaves
   stats(2) = 0.d0

   if( table%leaves < 1 ) return
   
   nnx  = table%nx
   nnf  = table%nf
   nxp1 = table%nxp1
   jfe  = table%jfe 
   jhs  = table%jhs
   jhe  = table%jhe
   jhs  = min( jhs, jhe )
   
   call_flag   = info(80)
   max_call    = huge(1)
   leaf_stride = 1

   if( info(81) > 0 ) max_call    = info(81)
   if( info(82) > 0 ) leaf_stride = info(82)

   n_leaves = 0
   do id = 1, table%n_pt
      if( .not.associated( table%leaf_pt(id)%leaf ) ) cycle
	  n_leaves = n_leaves + 1
	  if( mod(n_leaves-1,leaf_stride) /= 0 ) cycle
	  if( n_leaves/leaf_stride > max_call-1 ) return

      stats(2) = stats(2) + 1.d0
	  leaf     => table%leaf_pt(id)%leaf

	  i_wrk(1:3)          = (/ call_flag, leaf%id, n_leaves /)
	  r_wrk(1:nnx)         = table%xscale(1:nnx)
	  r_wrk(nxp1:jfe)     = table%fscale(1:nnf)
	  r_wrk(jfe+1)        = leaf%etolsq
	  r_wrk(jfe+2:jfe+12) = leaf%props(1:10)
	  
	  replace = 0
	  ga      = leaf%g(1:nnf,1:nnx)
	  
	  call usrfgh( replace, nnx, leaf%xfh(1:nnx), nnf, table%nh, i_wrk, r_wrk, &
	               leaf%xfh(nxp1:jfe), ga, leaf%xfh(jhs:jhe) )
	               
	  if( replace(2) == 1 ) leaf%g(1:nnf,1:nnx) = ga
	    
   end do

end subroutine isat_leaf_usrfgh  !------------------------------------------------

end subroutine isatab
