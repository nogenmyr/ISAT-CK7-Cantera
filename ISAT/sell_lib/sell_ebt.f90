!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module sell_ebt

!  Module for searching sets of ellipsoids (SELLs) using ellipsoidal
!  binary trees (EBTs).

!  An EBT has nodes and leaves, with a bounding ellipsoid (BE) associated 
!  with each.  The BE associated with a leaf is an ellipsoid (ELL) in the SELL.
!  Each node has two children (1 and 2) and the node's BE
!  covers those of its children's BEs.  The EBT supports querying to determine
!  whether a given point x is covered by one or more of the leaf ELLs.

!  S.B. Pope  2/14/2006, 5/17/2006

!  Data structures
!     ebt_type   - for EBTs
!     be_type    - for bounding ellipsoids (BE)
!     be_pointer - ponter to BE

!================================================

!  Subroutines: public

! be_create       !  allocate and initialize BE to null settings
! be_destroy      !  deallocate and nullify BE

! ebt_initialize  !  initialize EBT
! ebt_destroy     !  destroy EBT
! ebt_add         !  add leaf BE (for ELL) to EBT
! ebt_remove      !  remove leaf BE from EBT
! ebt_purge       !  remove all BEs from EBT
! ebt_mark_set    !  mark nodes of EBT to be updated
! ebt_mark_update !  update marked nodes of EBT
! ebt_update      !  update nodes of EBT based on change to BE
! ebt_rebuild     !  rebuild EBT
! ebt_traverse    !  traverse EBT based on x

! ebt_query       !  find an ELL (if any) containing x
! ebt_query_pair  !  find an ELL (if any) containing x, using affine space
! ebt_query_list  !  return list of ELLs containing x
! ebt_query_pair_list  !  return list of ELLs containing x, using affine space

! ebt_status      !  check EBT for problems
! ebt_param_set   !  set EBT parameters
! ebt_param_get   !  get EBT parameters
! ebt_prop        !  get EBT properties
! ebt_write       !  write EBT to logical unit LU (for checkpointing)
! ebt_read        !  read EBT from logical unit LU

!  Subroutines: private

! be_write          !  write BE to logical unit LU
! be_read           !  create and read BE from logical unit LU
! ebt_query_which   !  decide on which child to visit first in search
! ebt_query_success !  update node props of EBT based on BE
! ebt_be_pt_double  !  double size of pointer array

use sell_m
implicit none

   private ebt_query_which, ebt_query_success, ebt_be_pt_double

   type :: ebt_type  !  Ellipsoidal Binary Tree  (EBT) !--------
      type (sell_type), pointer    :: sell     !  SELL to which EBT belongs
      type(be_type),    pointer    :: root     !  root BE of EBT
	  type(be_pointer), pointer    :: be_pt(:) !  pointer to leaves: be_pt(id) 
	  integer                      :: n_pt     !  dimension of be_pt

      integer                      :: check    ! level of checking [0,1,2]
      type (id_list_type), pointer :: idlist   ! data structure used in id_list to generate IDs

      real(k_d)                    :: prop(10) ! properties of EBT - see ebt_prop

      integer                      :: force_search   ! >0 to force search outside bounding hyperplanes
      real(k_d)                    :: pair_sep_qual  !  parameter (less than 1.0) controlling
!                                        the quality of separating hyperplanes, used as cutting planes.
!                                        Call ebt_param_set to change.
      real(k_d)                    :: det_lim  !  limit on determinant of G allowed in content calculation 
   end type ebt_type

   type :: be_type   !  Bounding Ellipsoid (BE) ----------------
      type (ebt_type),   pointer :: ebt       ! EBT to which BE belongs
      type (be_type),    pointer :: parent    ! be%parent  is the parent of BE
      type (be_pointer), pointer :: child(:)  ! be%child(k)%be  is the k-th child of BE
	  integer                    :: id        ! be%id is the ID of BE's leaf (if BE is a leaf of the EBT)
	                                          ! be%id < 0 indicates that BE is a node
	  integer                    :: k         ! BE is parent's k-th child (k = 1 or 2)

	  real(k_d), pointer         :: geom(:)   ! geometry of BE {c, gg, r_in, r_out} (for leaf, points to ell%geom)
	  real(k_d), pointer         :: p(:)      ! point on cutting hyperplane
	  real(k_d), pointer         :: v(:)      ! cutting plane vector:  s(x) = v^T * (x-p)
	  real(k_d)                  :: smin(2)   ! min(s) for x in children of BE
	  real(k_d)                  :: smax(2)   ! max(s) for x in children of BE
	  real(k_d)                  :: prop(6)   ! properties of BE (see below)
	  integer                    :: k_first   ! first child to be processed
	  logical                    :: done(2)
	  logical                    :: marked    ! indicator for ebt_mark_update
   end type be_type  !------------------------------------------

   type :: be_pointer  
      type (be_type), pointer :: be
   end type be_pointer

!  be%prop:
! 1 - number of leaves in sub-tree defined by BE
! 2 - sum of content (volume) of leaves in sub-tree defined by BE
! 3 - sum of radii of circumscribed balls of leaves in sub-tree defined by BE
! 4 - index (k=1 or 2) of child visited on last successful query
! 5 - total number of visits to child 1 on successful queries
! 6 - total number of visits to child 2 on successful queries

contains  !===================================================================

subroutine be_create( ebt, node_leaf, be )  !  allocate and initialize BE to null settings

   type (ebt_type), pointer :: ebt       !  (in)  EBT to which BE belongs
   integer, intent(in)      :: node_leaf ! =1 for node; =2 for leaf
   type (be_type),  pointer :: be        !  (out) BE

! Note that be%geom is not allocated

   if( .not.associated( ebt ) ) then  ! check that EBT exists
      write(0,*)'be_create; ebt not associated'
	  stop
   endif

   if( node_leaf < 1  .or.  node_leaf > 2 ) then
      write(0,*)'be_create; bad value of node_leaf = ', node_leaf
	  stop
   endif

   allocate( be )
   be%ebt => ebt
   nullify(  be%parent )

   be%id      = -1
   be%k       = -1
   be%k_first = -1
   be%done    = .false.
   be%marked  = .false.
   be%smin    = -1.d0
   be%smax    = -1.d0
   be%prop    =  0.d0

   if( node_leaf == 1 ) then  ! node
      allocate( be%child(2) )
      nullify(  be%child(1)%be )
      nullify(  be%child(2)%be )
      allocate( be%geom( ebt%sell%ngeom ) )
      allocate( be%p(ebt%sell%nx) )
      allocate( be%v(ebt%sell%nx) )
	  be%geom    = -1.d0
	  be%p       = -1.d0
      be%v       = -1.d0
   else                        ! leaf
      nullify( be%child )               
      nullify( be%geom )               
      nullify( be%p )               
      nullify( be%v )               
   endif

   return
end subroutine be_create  !------------------------------------------------

subroutine be_destroy( be )  !   deallocate and nullify BE

   type (be_type), pointer :: be

   if( .not.associated( be ) ) then  ! check that BE exists
      write(0,*)'be_destroy; be not associated'
	  stop
   endif

   if( be%id < 1 ) then  !  node
      deallocate( be%child )
      deallocate( be%geom )
      deallocate( be%p )
      deallocate( be%v )
   else
      nullify( be%geom )
   endif

   deallocate( be )

   return
end subroutine be_destroy  !------------------------------------------------

subroutine be_write( be, lu )  !  write BE to logical unit LU

   type (be_type),  pointer :: be 
   integer, intent(in)      :: lu

   integer :: nx, ngeom

   nx    = be%ebt%sell%nx
   ngeom = be%ebt%sell%ngeom

! Note: be%k_first, be%done(1:2) are not written (a) because they are temporary, and not needed
!       for a perfect restart (b) to ensure repeatability
 
   write(lu, err=100 ) be%id
   write(lu, err=110 ) be%k, be%smin(1:2), be%smax(1:2), be%prop(1:6)

   if( be%id < 0 ) write(lu, err=120) be%geom(1:ngeom), be%p(1:nx), be%v(1:nx)

   return

100  write(0,*)'be_write: error writing id'
     stop

110  write(0,*)'be_write: error writing k'
     stop

120  write(0,*)'be_write: error writing geom'
     stop

end subroutine be_write  !------------------------------------------------

subroutine be_read( ebt, be, lu )  !  create and read BE from logical unit LU

   type (ebt_type), pointer :: ebt       !  (in)  EBT to which BE belongs
   type (be_type),  pointer :: be        !  (out) BE
   integer, intent(in)      :: lu

   integer :: id, node_leaf, nx, ngeom

   nx    = ebt%sell%nx
   ngeom = ebt%sell%ngeom

   read(lu, err=100 ) id

   node_leaf = 1
   if( id > 0 ) node_leaf = 2

   call be_create( ebt, node_leaf, be )

   be%id = id
   read(lu, err = 110 ) be%k, be%smin(1:2), be%smax(1:2), be%prop(1:6)

   if( be%id < 0 ) read(lu, err = 120 ) be%geom(1:ngeom), be%p(1:nx), be%v(1:nx)

   return

100  write(0,*)'be_read: error reading id'
     stop

110  write(0,*)'be_read: error reading k'
     stop

120  write(0,*)'be_read: error reading geom'
     stop

123  write(0,*)'be_read: error reading geom(1)'
     stop

end subroutine be_read  !------------------------------------------------

subroutine ebt_initialize( sell, ebt, n_pt, check, idlist, force_search, pair_sep_qual, det_lim )

!  Allocate and initialize (empty) EBT.  Optionally set check.

   type (sell_type), pointer              :: sell    ! (in) SELL on which EBT is to be formed
   type (ebt_type),  pointer              :: ebt     ! (out) data structure for EBT

!  Optional
   integer, intent(in), optional          :: n_pt    ! initial dimension of pointer array
   type (id_list_type), pointer, optional :: idlist  ! (in) data structure used in id_list to generate IDs

   integer, intent(in), optional          :: check   ! level of checking  =0 minimal; =1 some; =2 maximal
   integer, intent(in), optional          :: force_search   ! >0 to force search outside bounding hyperplanes
   real(k_d), intent(in), optional        :: pair_sep_qual  !  quality parameter (less than 1.0) 
   real(k_d), intent(in), optional        :: det_lim !  smallest determinant to use in content calc.

   integer :: i

   allocate( ebt )
   ebt%sell => sell
   nullify(ebt%root)

   ebt%n_pt = n_pt_0
   if( present(n_pt) ) ebt%n_pt = n_pt
   allocate( ebt%be_pt( ebt%n_pt ) )

   do i = 1, ebt%n_pt
      nullify( ebt%be_pt(i)%be )
   end do

   if( present(idlist) ) then
      ebt%idlist => idlist
   else
      nullify( ebt%idlist)
   endif

   ebt%prop = 0.d0
   ebt%prop(4) = huge(1.d0)

   if( present(check) ) then
      call ebt_param_set( ebt, check )
   else
      ebt%check = 0
   endif

   if( present(force_search) ) then
      ebt%force_search = force_search 
   else
      ebt%force_search = 0
   endif

   if( present(pair_sep_qual) ) then
      ebt%pair_sep_qual = pair_sep_qual
   else
      ebt%pair_sep_qual = 0.9d0
   endif

   if( present(det_lim) ) then
      ebt%det_lim = det_lim
   else
      ebt%det_lim = tiny(1.d0)*1.d9
   endif

   return
end subroutine ebt_initialize  !------------------------------------------------------------

subroutine ebt_destroy( ebt )

!  Destroy and deallocate EBT ebt.
!  All BEs should first be destroyed by a call to ebt_purge

   type (ebt_type), pointer :: ebt    ! (in) data structure for EBT

   if( .not.associated(ebt) ) then
      write(0,*) 'ebt_destroy: ebt not associated; stopping'
	  stop
   endif

   deallocate( ebt%be_pt )
   deallocate( ebt )

   return
end subroutine ebt_destroy  !------------------------------------------------------------

subroutine ebt_add( ebt, id, ante_up, pair_cover, cut_update )  !  add ELL (ID=id) to EBT

! Based on the center of ELL, EBT is traversed to identify BE_SIBLING.
! A new BE node is added (in place of BE_SIBLING) whose 1st and 2nd 
! children are BE and BE_SIBLING, respectively.

   implicit none

   type (ebt_type), pointer :: ebt        ! EBT
   integer, intent(in)      :: id         ! ID of ELL
   logical, intent(in)      :: ante_up    ! .true. to update antecedent BEs
   integer, intent(in)      :: pair_cover ! ={1,2,3} algorithm to use
   integer, intent(in)      :: cut_update ! >=0 to update cutting planes

   type (ell_type), pointer :: ell
   integer                  :: ksp

   type (be_type), pointer  :: be, be_parent 
   type (be_type), pointer  :: be_sibling
   integer                  :: id_status, ids_assigned, id_max, nx, ng, ng_end
   logical                  :: intersect
   real(k_d)                :: det

   if( .not.associated(ebt) ) then  !  check input
      write(0,*) 'ebt_add: ebt not associated; stopping'
	  stop
   endif

   if( ebt%check > 1  .and.  associated( ebt%idlist ) ) then

	  call id_query( ebt%idlist, id, id_status, ids_assigned, id_max )
	  if( id_status /= 1 ) then
		 write(0,*) 'ebt_add: id has not been assigned, id = ', id
	     stop
	  endif
   endif

   do while ( id > ebt%n_pt )  !  increase size of be_pt as needed
      call ebt_be_pt_double( ebt ) 
   end do

   nx     = ebt%sell%nx
   ng     = ebt%sell%ng
   ng_end = nx + ng

!  initialize BE for new leaf
   call be_create( ebt, 2, be )
  
!  set geometry of BE from ELL

   ell => ebt%sell%ell_pt(id)%ell
   be%geom => ell%geom

!  set BE props 
   call ell_chol_det( nx, be%geom(nx+1:ng_end), det )

   be%prop(1) = 1.d0                            !  one leaf
   be%prop(2) = 1.d0 / max( det, ebt%det_lim )  !  content
   be%prop(3) = be%geom(nx+ng+2)                !  r_out

   ebt%prop(4) = min( det, ebt%prop(4) ) 
!  set other properties

   be%ebt  => ebt
   be%id   = id
   ebt%be_pt(id)%be => be

! traverse EBT to leaf = be_sibling
   call ebt_traverse( ebt, be%geom(1:ebt%sell%nx), be_sibling )

! special treatment for first BE in EBT
   if( .not. associated( be_sibling ) ) then
      ebt%root => be   !  BE is root of EBT
	  return
   endif

!  set  be  to 1st child; be_sibling to second child of be_parent

!  allocate and initialize BE = be_parent for new node
   call be_create( ebt, 1, be_parent )

   if( associated( be_sibling%parent ) ) then  !  be_sibling was not root
      be_parent%parent => be_sibling%parent
	  ksp              =  be_sibling%k
      be_parent%k      =  ksp
	  be_sibling%parent%child(ksp)%be => be_parent
   else                                        !  be_sibling was root
	  ebt%root => be_parent       !  be_parent is now root
   endif

   be_parent%ebt         => ebt
   be_parent%child(1)%be => be
   be_parent%child(2)%be => be_sibling

!  re-set be_sibling

   be_sibling%parent => be_parent
   be_sibling%k      =  2

!  set  be

   be%parent => be_parent
   be%k      =  1

!  set be_parent's cutting plane if it will not be set in ebt_update
      if( cut_update < 0 ) then

	     call ell_pair_separate( nx, be%geom(1:nx), be%geom(nx+1:ng_end), &
		                         be_sibling%geom(1:nx), be_sibling%geom(nx+1:ng_end), &
								 ebt%pair_sep_qual, 0, be_parent%p, be_parent%v, intersect )							 
      endif

!  update BEs of all be's antecedents

   if( ante_up ) call ebt_update( ebt, be, pair_cover, cut_update )

   nullify( be )
   nullify( be_sibling )
   nullify( be_parent )

   return
end subroutine ebt_add  !----------------------------------------------------------

subroutine ebt_remove( be, pair_cover, cut_update )  !  remove leaf BE from EBT
   implicit none

   type (be_type), pointer   :: be      ! given E
   integer, intent(in)       :: pair_cover ! algorithm to use
   integer, intent(in)       :: cut_update ! >=0 to update cutting planes

   integer                   :: k, ks, id
   type (ebt_type), pointer  :: ebt
   type (be_type), pointer   :: be_sibling, be_parent, be_grandparent 

   if( .not.associated(be) ) then  !  check input
      write(0,*) 'ebt_remove: be not associated; stopping'
	  stop
   endif

   if( be%id < 1 ) then
      write(0,*)'ebt_remove: be is not leaf'
	  stop
   endif

   ebt => be%ebt
   id  = be%id

   if( .not.associated(be%parent) ) then  !  BE is root of EBT
      nullify( ebt%root )                 !  EBT is empty
	  call be_destroy( be )
	  nullify( ebt%be_pt( id )%be )
	  nullify( ebt )
	  return
   endif

   be_parent  => be%parent     ! identify BE's parent and sibling
   k          =  be%k 
   ks         =  3 - k
   be_sibling => be_parent%child(ks)%be

   if( .not.associated(be_parent%parent) ) then  !  be_parent is root of EBT
      be%ebt%root => be_sibling                  !  be's sibling is sole surviving BE
	  be_sibling%k = 0
	  nullify( be_sibling%parent )

   else  !  grandparent exists

      be_grandparent    => be_parent%parent   !  new parent of sibling
      be_sibling%parent => be_grandparent
	  k                 = be_parent%k         !  index of parent's child (old and new)
      be_sibling%k      = k
      be_grandparent%child(k)%be => be_sibling
	  call ebt_update( ebt, be_sibling, pair_cover, cut_update ) 
	  nullify( be_grandparent )
   endif
  
   call be_destroy( be )
   call be_destroy( be_parent )
   nullify( be_sibling )
   nullify( ebt%be_pt( id )%be )
   nullify( ebt )

   return
end subroutine ebt_remove   !-------------------------------------------------------

subroutine ebt_purge( ebt )   !  remove all BEs from EBT

   implicit none
   type (ebt_type), pointer :: ebt   ! given EBT

   type (be_type), pointer  :: be, parent
   integer                  :: i

   if( .not.associated(ebt) ) then  !  check input
      write(0,*) 'ebt_purge: ebt not associated; stopping'
	  stop
   endif

!  start at root (if it exists)

   if( .not.(associated(ebt%root) ) ) return  !  already purged

   be => ebt%root  !  root
   nullify( parent )            !  ...has no parent

   do   !  on each iteration of loop,  be  and  parent  must be defined (parent may be null)

      if( be%id < 1 ) then  !  if BE is a node, first remove any children
	   
	   if( associated(be%child(1)%be) ) then      !  go to child(1) if it exists
	     parent => be
		 be     => be%child(1)%be
		 nullify( parent%child(1)%be )
		 cycle

       elseif( associated(be%child(2)%be) ) then  !  go to child(2) if it exists
	     parent => be
		 be     => be%child(2)%be
		 nullify( parent%child(2)%be )
		 cycle
       endif

      endif

!  be  is childless (leaf or node): destroy and go to parent (if it exists)
	     
	  call be_destroy( be )                   !  destroy  be

	  if( .not.associated( parent ) ) exit    !  all done

	  be => parent                            !  go to parent

	  if( associated(be%parent) ) then        !  identify parent (if any)
	     parent => be%parent
	  else
	     nullify( parent )
	  endif

   end do

! re-set related components of EBT

   nullify( ebt%root )

   do i = 1, ebt%n_pt
      nullify( ebt%be_pt(i)%be )
   end do

   ebt%prop = 0.d0  ! re-initialize ebt%prop
   ebt%prop(4) = huge(1.d0)

   return

end subroutine ebt_purge  !------------------------------------------------------------------

subroutine ebt_mark_set( ebt, be0, n_mark )  !  update nodes of EBT based on BE0

!  Mark all antecedent nodes of leaf BE0  in EBT.

   implicit none

   type (ebt_type), pointer :: ebt         ! given EBT
   type (be_type), pointer  :: be0         ! starting BE leaf
   integer, intent(out)     :: n_mark      ! number of nodes marked

   type (be_type), pointer :: be

   if( .not.associated(ebt) ) then  !  check input
      write(0,*) 'ebt_mark_set: ebt not associated; stopping'
	  stop
   endif

   if( .not.associated(be0) ) then 
      write(0,*) 'ebt_mark_set: be0 not associated; stopping'
	  stop
   endif

   if( be0%id < 1 ) then
      write(0,*) 'ebt_mark_set: be0 not a leaf; stopping'
	  stop
   endif

   n_mark = 0
   be     => be0
   do  !  loop over antecedent nodes

   	  if( .not. associated( be%parent ) ) return  !  be is the root: all done

	  be => be%parent  !  go to parent
	  be%marked = .true.
	  n_mark    = n_mark + 1 

   end do

   return
end subroutine ebt_mark_set  !------------------------------------------------------------------

subroutine ebt_mark_update( ebt, pair_cover, cut_update, n_update )  !  update marked nodes of EBT

!  Update the geometry of all marked nodes in EBT.
!  The structure of EBT is unaffected.

   implicit none

   type (ebt_type), pointer :: ebt         ! given EBT
   integer, intent(in)      :: pair_cover  ! ellipsoid pair covering algorithm to be used
   integer, intent(in)      :: cut_update  ! <0  - cutting planes are not updated
!                                          ! >=0 - cutting planes are updated.
   integer, intent(out)     :: n_update    ! number of nodes updated

!  pair_cover = 1, 2, 3 : spheroid, covariance, iterative algorithms, resp.

!  For cut_update >= 0,  the cutting planes are updated, using at most  cut_update
!  iterations in ell_pair_separate.

   integer   :: k, nx, ng, ng_end, j
   logical   :: intersect
   real(k_d) :: c(ebt%sell%nx), c1(ebt%sell%nx), c2(ebt%sell%nx), &
                gg(ebt%sell%ng), gg1(ebt%sell%ng), gg2(ebt%sell%ng), r_in, r_out

   type (be_type), pointer :: be

   if( .not.associated(ebt) ) then  !  check input
      write(0,*) 'ebt_mark_update: ebt not associated; stopping'
	  stop
   endif

   n_update = 0
   if( .not.associated(ebt%root) ) return  !  nothing to do if no root
   if( ebt%root%id > 0 )      return  !  nothing to do if EBT consists of a single leaf
   if( .not.ebt%root%marked ) return  !  nothing to do if root (and hence all nodes) not marked

!  set indices to local variables
   nx     = ebt%sell%nx   
   ng     = ebt%sell%ng   
   ng_end = nx + ng

!--------  traverse EBT

   be => ebt%root
   traverse: do

      do k = 1, 2  !  go to either child node if marked
         if( be%child(k)%be%marked ) then
		    be => be%child(k)%be
			cycle traverse
		 endif
	  end do

! no marked children: update this node

      n_update  = n_update + 1
	  be%marked = .false.

!  get BEs of both children
      c1  = be%child(1)%be%geom(1:nx)
      c2  = be%child(2)%be%geom(1:nx)
      gg1 = be%child(1)%be%geom(nx+1:ng_end)
      gg2 = be%child(2)%be%geom(nx+1:ng_end)

!  form bounding ellipsoid

	  call ell_pair_cover( pair_cover, nx, c1, gg1, c2, gg2, c, gg )

	  be%geom(1:nx)        = c
	  be%geom(nx+1:ng_end) = gg

      call ell_radii( nx, gg(1:ng), r_in, r_out )
	  be%geom(nx+ng+1) = r_in
	  be%geom(nx+ng+2) = r_out

!  form cutting plane
      if( cut_update >= 0 ) &
	     call ell_pair_separate( nx, c1, gg1, c2, gg2, ebt%pair_sep_qual, cut_update, &
		                         be%p, be%v, intersect )

!  project children onto normal line
	  call ell_line_proj( nx, c1, gg1, be%p, be%v, be%smin(1), be%smax(1) )
	  call ell_line_proj( nx, c2, gg2, be%p, be%v, be%smin(2), be%smax(2) )

!  properties of BE's sub-tree
      do j=1,3
         be%prop(j) = be%child(1)%be%prop(j) + be%child(2)%be%prop(j) 
	  end do
! end of updating:

   	  if( .not. associated( be%parent ) ) return  !  be is the root: all done

	  be => be%parent  !  go to parent

   end do traverse

   return
end subroutine ebt_mark_update  !------------------------------------------------------------------

subroutine ebt_update( ebt, be0, pair_cover, cut_update )  !  update nodes of EBT based on BE0

!  Update the geometry of all antecedent nodes of BE0  in EBT.
!  BE0 may be a leaf or a node.   Prior to this call, BE0 has been added 
!  to the EBT, or its geometry  has been modified.
!  The structure of EBT is unaffected.

   implicit none

   type (ebt_type), pointer :: ebt         ! given EBT
   type (be_type), pointer  :: be0         ! starting BE (leaf or node)
   integer, intent(in)      :: pair_cover  ! ={1,2,3} algorithm to use
   integer, intent(in)      :: cut_update  ! <  0 - cutting planes are not updated
!                                          ! >= 0 - cutting planes are updated.

!  For cut_update >= 0,  the cutting planes are updated, using at most  cut_update
!  iterations in ell_pair_separate.

   integer   :: nx, ng, ng_end, j
   logical   :: intersect
   real(k_d) :: c(ebt%sell%nx), c1(ebt%sell%nx), c2(ebt%sell%nx), &
                gg(ebt%sell%ng), gg1(ebt%sell%ng), gg2(ebt%sell%ng), r_in, r_out

   type (be_type), pointer :: be

   if( .not.associated(ebt) ) then  !  check input
      write(0,*) 'ebt_update: ebt not associated; stopping'
	  stop
   endif

   if( .not.associated(be0) ) then 
      write(0,*) 'ebt_update: be0 not associated; stopping'
	  stop
   endif

!  set indices to local variables
   nx     = ebt%sell%nx   
   ng     = ebt%sell%ng   
   ng_end = nx + ng

   be => be0
   do  !  loop over antecedent nodes

   	  if( .not. associated( be%parent ) ) then  !  be is the root: all done
         nullify( be )
		 return
      endif

	  be => be%parent  !  go to parent

!  get BEs of both children
      c1  = be%child(1)%be%geom(1:nx)
      c2  = be%child(2)%be%geom(1:nx)
      gg1 = be%child(1)%be%geom(nx+1:ng_end)
      gg2 = be%child(2)%be%geom(nx+1:ng_end)

!  form bounding ellipsoid

	  call ell_pair_cover( pair_cover, nx, c1, gg1, c2, gg2, c, gg )
	  be%geom(1:nx)        = c
	  be%geom(nx+1:ng_end) = gg

      call ell_radii( nx, gg(1:ng), r_in, r_out )
	  be%geom(nx+ng+1) = r_in
	  be%geom(nx+ng+2) = r_out

!  form cutting plane
      if( cut_update >= 0 ) &
	     call ell_pair_separate( nx, c1, gg1, c2, gg2, ebt%pair_sep_qual, cut_update, &
		                         be%p, be%v, intersect )

!  project children onto normal line
	  call ell_line_proj( nx, c1, gg1, be%p, be%v, be%smin(1), be%smax(1) )
	  call ell_line_proj( nx, c2, gg2, be%p, be%v, be%smin(2), be%smax(2) )

!  properties of BE's sub-tree
      do j=1,3
         be%prop(j) = be%child(1)%be%prop(j) + be%child(2)%be%prop(j) 
	  end do

   end do

   nullify( be )
   return
end subroutine ebt_update  !------------------------------------------------------------------

subroutine ebt_rebuild( ebt, mode, pair_cover )   !  rebuild EBT

   implicit none
   type (ebt_type), pointer :: ebt        !  EBT to be rebuilt
   integer, intent(in)      :: mode       !  method of rebuilding
   integer, intent(in)      :: pair_cover ! algorithm to use

! mode = 1 - retain the existing topology
! mode = 2 - purge EBT, and then re-add ELLs in a random order

   integer              :: i, id, ids_assigned, id_max
   integer, allocatable :: id_array(:)

   if( .not.associated(ebt) ) then  !  check input
      write(0,*) 'ebt_rebuild: ebt not associated; stopping'
	  stop
   endif

   if( mode == 1 ) then  !  update antecedent BEs 

      do id = 1, ebt%n_pt  !  mark all leaves
	     if( associated( ebt%be_pt(id)%be ) )  &
		    call ebt_mark_set( ebt, ebt%be_pt(id)%be, i )
	  end do

	  call ebt_mark_update( ebt, pair_cover, 0, id_max )  ! update

   elseif( mode ==2 ) then  !  purge and re-add ELLs in random order

	  call id_stats( ebt%idlist, ids_assigned, id_max )
	  allocate( id_array(ids_assigned) )
	  call id_rand_array( ebt%idlist, id_array )
      call ebt_purge( ebt ) 

      id_max = 0
	  do i = 1, ids_assigned
	     id = id_array(i)
	     if( associated(ebt%sell%ell_pt(id)%ell) )  then
		    call ebt_add( ebt, id, .false., pair_cover, -1 )
			id_max = id_max + 1
	     endif
	  end do

      id_max = 0
	  do id = 1, ebt%n_pt  !  mark all leaves
	     if( associated( ebt%be_pt(id)%be ) )  &
		    call ebt_mark_set( ebt, ebt%be_pt(id)%be, i )
			id_max = id_max + i
	  end do

      i = id_max
	  call ebt_mark_update( ebt, pair_cover, 0, id_max )  ! update

	  deallocate( id_array )

   else
      write(0,*)'ebt_rebuild: invalid mode = ', mode
	  stop
   endif

   return
end subroutine ebt_rebuild  !------------------------------------------------------------------

subroutine ebt_traverse( ebt, x, be )  !  traverse EBT based on x

!  Traverse EBT based on x to indentify leaf BE.
!  If the root of the EBT does not exist, nullify(be).

   implicit none

   type (ebt_type), pointer  :: ebt            ! given EBT
   real(k_d), intent(in)     :: x(ebt%sell%nx) ! query point
   type (be_type), pointer   :: be             ! BE

   real(k_d)                 :: s

   if( .not.associated( ebt ) ) then  ! check that EBT exists
      write(0,*)'ebt_traverse; ebt not associated'
	  stop
   endif

   if( .not.associated( ebt%root) ) then 
      nullify( be ) !  indicate that root does not exist
	  return
   endif

   be => ebt%root  !  start traverse at root

   do  !  loop over BEs  ------------------------------------------

      if( be%id > 0 ) return  !  leaf has been reached - all done

	  s = sum( be%v * (x-be%p) )  !  compute s
      if( s < 0.d0 ) then   !  go to child 1 or 2 depending on s
	     be => be%child(1)%be
	  else
	     be => be%child(2)%be
	  endif

   end do

end subroutine ebt_traverse  !---------------------------------------------------------------

subroutine ebt_query( ebt, x, mode, max_tests, be, tests )  !  find ELL containing x

!  Search EBT for the next BE leaf (ELL) containing x (if any).
!  For a given search, on the first call to ebt_query,  be  must be
!  nullified, and then the search starts of the root of the EBT.
!  On subsequent calls (for the same search), on entry  be  points to the leaf
!  identified on the previous call, and the search resumes from that leaf.
!  Use mode=mode in search (see ebt_query_which for explanation).
!  Perform no more than  max_tests  ellipsoid-covering tests.
!  Return BE in be, or nullify(be) if x is not in any (further) BE.

   implicit none

   type (ebt_type), pointer :: ebt             ! EBT
   real(k_d), intent(in)    :: x(ebt%sell%nx)  ! query point
   integer, intent(in)      :: mode            ! search mode ( 1 <= |mode| <= 8 )
   integer, intent(in)      :: max_tests       ! maximum number of leaf tests allowed
   type (be_type), pointer  :: be              ! BE
   integer, intent(out)     :: tests           ! number of leaf tests performed

   real(k_d) :: s
   integer   :: nx, ngeom, ka, kb
   logical   :: from_above, in

   if( .not.associated(be) ) then        !  fresh search
      if( associated( ebt%root) ) then   !  start at root (if it exists)
	     be => ebt%root
		 from_above = .true.
	  else
		 write(0,*)'ebt_query: no root'  !  root does not exist
		 return
	  endif

   else                                  !  resume previous search
      if( be%id > 0 ) then
		 from_above = .false.
	  else
	  	 write(0,*)'ebt_query: starting leaf not found'
		 stop
	  endif

   endif

   nx    = ebt%sell%nx
   ngeom = ebt%sell%ngeom

   tests = 0
   be_loop: do  !  loop over BEs  

!  At the start of the loop,  be  and from_above must be set.
!  If  from_above==.true., BE has been entered from above, and therefore for the first time.
!  If  from_above==.false., then move from  be  to  be%parent, which has been entered before,
!     at which time, be%k_first and be%done(1:2) had been set.

   if( from_above ) then  !  have entered BE from above (and therefore for the first time) 

      call sell_ell_pt_in( nx, ngeom, be%geom, x, in )

      if( be%id > 0 ) then  ! BE is a leaf
	  	 tests = tests + 1
	     if( in ) return    ! x is in leaf BE -- all done
		 if( tests >= max_tests ) then
		    nullify(be)     ! search failed
			return
		 endif
      endif 

	  if( .not.in ) then      !  x is not in BE
		 from_above = .false. !  move up to parent
		 cycle be_loop
	  endif

!  x is in BE which is not a leaf

      s          = sum( be%v * (x-be%p) )  !  compute s

      if( ebt%force_search == 0 ) then
!  set be%done(k) = .false. if child k exists, and x lies between its supporting hyperplanes
         be%done    = .true.         !  assume x not in children until found otherwise
	     be%k_first = 1              !  assume, reset as necessary

         if( s >= be%smin(1)  .and.  s <= be%smax(1) ) then
	        be%done(1) = .false.

            if( s >= be%smin(2)  .and.  s <= be%smax(2) ) then
		       be%done(2) = .false.  !  both directions are possible
		       if( mode == 2 ) then
			      if( s >= 0.d0 ) be%k_first = 2 
			   else
			      call ebt_query_which( be, s, mode, be%k_first )
			   endif
		    endif

	     elseif( s >= be%smin(2)  .and.  s <= be%smax(2) ) then
	        be%done(2) = .false.
            be%k_first = 2 
	     endif

      else  ! force_search > 0
	    be%done = .false.
	    call ebt_query_which( be, s, mode, be%k_first )
      endif


   else   !  .not.from_above:  move from be to be%parent

      if( .not.associated(be%parent) ) then   !  BE is root
		 nullify(be)                          !  at root:  no ELL containing x found 
	     return
	  endif

	  be%parent%done( be%k ) = .true. !  BE has been processed
	  be => be%parent                 !  moved to parent

   endif

!  x is in BE which is a node:  determine whether to process a child or go to parent 

   ka = be%k_first
   kb = 3 - ka
   if( .not.be%done(ka) ) then      !  go down to child ka
      from_above = .true.
	  be => be%child(ka)%be

   elseif( .not.be%done(kb) ) then  !  go down to child kb
      from_above = .true.
	  be => be%child(kb)%be

   else                             !  both children processed: go to parent
	  from_above = .false.
   endif

   end do be_loop 

end subroutine ebt_query  !-------------------------------------------------------------------

subroutine ebt_query_pair( sell_f, x_f, ebt_a, x_a, mode, max_f_tests, max_be_tests, &
                           be, f_tests, a_tests, be_tests )  

!  Search for an ELL in SELL_F which covers x_f.

!  SELL_F is a set of ellipsoids in the "full" space, and x_f is a given point in this space.
!  SELL_A is the corresponding set of ellipsoids projected onto a lower dimensional "affine" space,
!  and x_a is the projection of x_f.  EBT_A is an ellipsoidal binary tree constructed on SELL_A.
!  EBT_A is searched to find an ellipsoid in SELL_A which covees x_a; and then the corresponding 
!  ellipsoid in SELL_F is tested to determine if it covers x_f.
!  Perform no more than max_be_tests tests on BEs in EBT_A, and no more than max_f_tests in SELL_F.
!  There may be 0, 1 or more ELLs satisfying the query.  
!  On entry, if BE is not associated, start the search at the root of EBT_A.  
!  Otherwise, start at ebt_a%be (for ebt_a%be being a leaf identified on the previous call).  
!  Return BE in be, or nullify(be) if x_f is not in any (further) ELL.

   implicit none

   type (sell_type), pointer :: sell_f             ! SELL_F
   real(k_d), intent(in)     :: x_f(sell_f%nx)     ! query point in full space
   type (ebt_type),  pointer :: ebt_a              ! EBT defined on SELL_A
   real(k_d), intent(in)     :: x_a(ebt_a%sell%nx) ! query point in affine space
   integer, intent(in)       :: mode               ! search mode ( 1 <= |mode| <= 8 )
   integer, intent(in)       :: max_f_tests        ! maximum number of ELL tests allowed in SELL_F
   integer, intent(in)       :: max_be_tests       ! maximum number of BE tests allowed in EBT_A
   type (be_type), pointer   :: be                 ! BE in EBT_A at which to start search (or nullified)
   integer, intent(out)      :: f_tests            ! number of ELL tests performed in SELL_F
   integer, intent(out)      :: a_tests            ! number of ELL tests performed in SELL_A
   integer, intent(out)      :: be_tests           ! number of BE tests performed in EBT_A

   type (ell_type), pointer  :: ell
   real(k_d) :: s
   integer   :: nx, ngeom, ka, kb
   logical   :: from_above, in

   a_tests  = 0
   f_tests  = 0
   be_tests = 0
   if( max_be_tests <=0  .or.  max_f_tests <= 0 ) return !  fast return if no testing allowed

   if( .not.associated(sell_f) ) then  !  check input
      write(0,*)'ebt_query_pair: sell_f not associated'
	  stop
   endif

   if( .not.associated(ebt_a) ) then
      write(0,*)'ebt_query_pair: ebt_a not associated'
	  stop
   endif

!  determine BE at which to start

   if( .not.associated(be) ) then            !  fresh search
      if( associated( ebt_a%root) ) then     !  start at root (if it exists)
	     be => ebt_a%root
		 from_above = .true.
	  else
		 write(0,*)'ebt_query_pair: no root' !  root does not exist
		 return
	  endif

   else                                      !  resume previous search
      if( be%id > 0 ) then
		 from_above = .false.
	  else
	  	 write(0,*)'ebt_query_pair: BE is not a leaf'
		 stop
	  endif

   endif

   nx    = ebt_a%sell%nx
   ngeom = ebt_a%sell%ngeom

   be_loop: do  !  loop over BEs  

!  At the start of the loop,  be  and from_above must be set.
!  If  from_above==.true., BE has been entered from above, and therefore for the first time.
!  If  from_above==.false., then move from  be  to  be%parent, which has been entered before,
!     at which time, be%k_first and be%done(1:2) had been set.

   if( from_above ) then  !  have entered BE from above (and therefore for the first time) 

      call sell_ell_pt_in( nx, ngeom, be%geom, x_a, in )
	  be_tests = be_tests + 1

      if( be%id > 0 ) then  ! BE is a leaf
	  	 a_tests = a_tests + 1

	     if( in ) then      ! in ELL_A: test for x_f covered by ELL_F
	        ell => sell_f%ell_pt(be%id)%ell
			call sell_ell_pt_in( sell_f%nx, sell_f%ngeom, ell%geom, x_f, in )
	        f_tests = f_tests + 1
			ebt_a%prop(5) = ebt_a%prop(5) + 1.d0  ! leaf tests
			if( .not.in ) ebt_a%prop(6) = ebt_a%prop(6) + 1.d0  !  false positives
	     endif

	     if( in ) return    !  covering ellipsoid found

	     if( be_tests >= max_be_tests  .or.  f_tests >= max_f_tests ) then
		    nullify(be)     ! search failed
			return
		 endif
      endif 

	  if( .not.in ) then      !  x is not in BE
		 from_above = .false. !  move up to parent
		 cycle be_loop
	  endif

!  x is in BE which is not a leaf
!  set be%done(k) = .false. if child k exists, and x lies between its supporting hyperplanes

      be%done    = .true.    !  assume x_a not in children until found otherwise
	  be%k_first = 1         !  assume, reset as necessary
      s          = sum( be%v * (x_a-be%p) )  !  compute s

      if( s >= be%smin(1)  .and.  s <= be%smax(1) ) then
	     be%done(1) = .false.

         if( s >= be%smin(2)  .and.  s <= be%smax(2) ) then
		    be%done(2) = .false.  !  both directions are possible
		    if( mode == 2 ) then
			   if( s >= 0.d0 ) be%k_first = 2 
			else
			   call ebt_query_which( be, s, mode, be%k_first )
			endif
		 endif

	  elseif( s >= be%smin(2)  .and.  s <= be%smax(2) ) then
	     be%done(2) = .false.
         be%k_first = 2 
	  endif

	  if( ebt_a%force_search >= 1 ) be%done = .false.

   else   !  .not.from_above:  move from be to be%parent

      if( .not.associated(be%parent) ) then   !  BE is root
		 nullify(be)                          !  at root:  no ELL containing x_f found 
	     return
	  endif

	  be%parent%done( be%k ) = .true. !  BE has been processed
	  be => be%parent                 !  moved to parent

   endif

!  x_a is in BE which is a node:  determine whether to process a child or go to parent 

   ka = be%k_first
   kb = 3 - ka
   if( .not.be%done(ka) ) then      !  go down to child ka
      from_above = .true.
	  be => be%child(ka)%be

   elseif( .not.be%done(kb) ) then  !  go down to child kb
      from_above = .true.
	  be => be%child(kb)%be

   else                             !  both children processed: go to parent
	  from_above = .false.
   endif

   end do be_loop 

end subroutine ebt_query_pair  !-------------------------------------------------------------------

subroutine ebt_query_list( ebt, x, max_tests, be, n_id_max, n_id, id, tests )  !  find ELL containing x

!  Given a query point x, search the SELL, using EBT, to attempt to find 
!  up to n_id_max ELLs which covers x. Perform no more than max_tests tests in this attempt.
!  There may be 0, 1 or more ELLs found to satisfy the query.  
!  The number of covering ELLs found is returned in n_id and their IDs in id(:).
!  On entry, if be is not associated, start the search at the root.  
!  Otherwise, start at ebt%be (for ebt%be being a leaf identified in the previous call).  
!  On exit, be points to the last leaf visited (or is nullified if no leaf encountered).

   implicit none

   type (ebt_type), pointer :: ebt             ! EBT
   real(k_d), intent(in)    :: x(ebt%sell%nx)  ! query point
   integer,   intent(in)    :: max_tests       ! maximum number of leaf tests allowed
   type (be_type), pointer  :: be              ! leaf in EBT at which to start search (or nullified)
   integer,   intent(in)    :: n_id_max        ! maximum number of IDs to be returned
   integer,   intent(out)   :: n_id            ! number of covering ELLs found
   integer,   intent(out)   :: id(n_id_max)    ! array of IDs of covering Ells
   integer,   intent(out)   :: tests           ! number of leaf tests performed

   type (be_type), pointer  :: be_last_leaf
   real(k_d) :: s
   integer   :: nx, ngeom, k, ka, kb
   logical   :: from_above, in

   n_id  = 0
   id    = 0
   tests = 0

   if( max_tests <= 0 ) return  !  quick return if no tests allowed

   if( .not.associated( ebt ) ) then
      write(0,*)'ebt_query_list: ebt not associated'
	  stop
   endif

   if( .not.associated(be) ) then        !  fresh search
      if( associated( ebt%root) ) then   !  start at root (if it exists)
	     be => ebt%root
		 from_above = .true.
		 nullify(be_last_leaf)
	  else
		 write(0,*)'ebt_query_list: no root'  !  root does not exist
		 return
	  endif

   else                                  !  resume previous search
      if( be%id > 0 ) then
		 from_above = .false.
		 be_last_leaf => be
	  else
	  	 write(0,*)'ebt_query_list: starting leaf not found'
		 stop
	  endif

   endif

   nx    = ebt%sell%nx
   ngeom = ebt%sell%ngeom

   be_loop: do  !  loop over BEs  

!  At the start of the loop,  be  and from_above must be set.
!  If  from_above==.true., BE has been entered from above, and therefore for the first time.
!  If  from_above==.false., then move from  be  to  be%parent, which has been entered before,
!     at which time, be%k_first and be%done(1:2) had been set.

   if( from_above ) then  !  have entered BE from above (and therefore for the first time) 

      call sell_ell_pt_in( nx, ngeom, be%geom, x, in )
	  tests = tests + 1

      if( be%id > 0 ) then  ! BE is a leaf
	     be_last_leaf => be
	     if( in ) then
		    n_id     = n_id + 1
			id(n_id) = be%id
			if( n_id >= n_id_max  .or. tests >= max_tests ) exit be_loop

		    from_above = .false. !  move back up to parent
		    cycle be_loop
		 endif

		 if( tests >= max_tests ) exit be_loop
      endif 

	  if( .not.in ) then      !  x is not in BE
		 from_above = .false. !  move up to parent
		 cycle be_loop
	  endif

!  x is in BE which is not a leaf
!  set be%done(k) = .false. if child k exists, and x lies between its supporting hyperplanes

      be%done  = .true.           !  assume x not in children until found otherwise

      s = sum( be%v * (x-be%p) )  !  compute s
	  do k = 1, 2                 !  x may be in child if between supporting hyperplanes
		if( s >= be%smin(k)  .and.  s <= be%smax(k) ) be%done(k) = .false.
	  end do

	  if( ebt%force_search >= 1 ) be%done = .false.

	  if( s < 0.d0 ) then  !  select child to be tested first, based on s
	     be%k_first = 1
	  else
	     be%k_first = 2
	  endif

   else   !  .not.from_above:  move from be to be%parent

      if( .not.associated(be%parent) ) exit be_loop   !  BE is root

	  be%parent%done( be%k ) = .true. !  BE has been processed
	  be => be%parent                 !  moved to parent

   endif

!  x is in BE which is a node:  determine whether to process a child or go to parent 

   ka = be%k_first
   kb = 3 - ka
   if( .not.be%done(ka) ) then      !  go down to child ka
      from_above = .true.
	  be => be%child(ka)%be

   elseif( .not.be%done(kb) ) then  !  go down to child kb
      from_above = .true.
	  be => be%child(kb)%be

   else                             !  both children processed: go to parent
	  from_above = .false.
   endif

   end do be_loop 

   if( associated(be_last_leaf) ) then  !  set BE to last leaf encountered (or nullify)
      be => be_last_leaf
   else
      nullify( be )
   endif

   return

end subroutine ebt_query_list  !-------------------------------------------------------------------

subroutine ebt_query_pair_list( sell_f, x_f, ebt_a, x_a, max_f_tests, max_be_tests, be,  &
                                n_id_max, n_id, id, f_tests, a_tests, be_tests ) 

!  Search for ELLs in SELL_F which covers x_f.

!  SELL_F is a set of ellipsoids in the "full" space, and x_f is a given point in this space.
!  SELL_A is the corresponding set of ellipsoids projected onto a lower dimensional "affine" space,
!  and x_a is the projection of x_f.  EBT_A is an ellipsoidal binary tree constructed on SELL_A.
!  EBT_A is searched to find ellipsoids in SELL_A which cover x_a; and then the corresponding 
!  ellipsoid in SELL_F is tested to determine if it covers x_f.
!  Perform no more than max_be_tests tests in EBT_A, and no more than max_f_tests in SELL_F.
!  There may be 0, 1 or more ELLs satisfying the query.  
!  The number of covering ELLs found is returned in n_id and their IDs in id(:).
!  On entry, if BE is not associated, start the search at the root of EBT_A.  
!  Otherwise, start at ebt_a%be (for ebt_a%be being a leaf identified on the previous call).  

   implicit none

   type (sell_type), pointer :: sell_f             ! SELL_F
   real(k_d), intent(in)     :: x_f(sell_f%nx)     ! query point in full space
   type (ebt_type),  pointer :: ebt_a              ! EBT defined on SELL_A
   real(k_d), intent(in)     :: x_a(ebt_a%sell%nx) ! query point in affine space
   integer, intent(in)       :: max_f_tests        ! maximum number of ELL tests allowed in SELL_F
   integer, intent(in)       :: max_be_tests       ! maximum number of BE tests allowed in EBT_A
   type (be_type), pointer   :: be                 ! BE in EBT_A at which to start search (or nullified)
   integer,   intent(in)     :: n_id_max           ! maximum number of IDs to be returned
   integer,   intent(out)    :: n_id               ! number of covering ELLs found
   integer,   intent(out)    :: id(n_id_max)       ! array of IDs of covering Ells
   integer, intent(out)      :: f_tests            ! number of ELL tests performed in SELL_F
   integer, intent(out)      :: a_tests            ! number of ELL tests performed in SELL_A
   integer, intent(out)      :: be_tests           ! number of BE tests performed in EBT_A

   type (ell_type), pointer  :: ell
   type (be_type), pointer   :: be_last_leaf
   real(k_d) :: s
   integer   :: nx, ngeom, k, ka, kb
   logical   :: from_above, in

   n_id     = 0
   id       = 0
   f_tests  = 0
   a_tests  = 0
   be_tests = 0

   if( max_f_tests <= 0  .or.  max_be_tests <= 0 ) return  !  quick return if no tests allowed

   if( .not.associated(sell_f) ) then  !  check input
      write(0,*)'ebt_query_pair_list: sell_f not associated'
	  stop
   endif

   if( .not.associated(ebt_a) ) then
      write(0,*)'ebt_query_pair_list: ebt_a not associated'
	  stop
   endif

!  determine BE at which to start

   if( .not.associated(be) ) then            !  fresh search
      if( associated( ebt_a%root) ) then     !  start at root (if it exists)
	     be => ebt_a%root
		 from_above = .true.
	  else
		 write(0,*)'ebt_query_pair_list: no root' !  root does not exist
		 return
	  endif

   else                                      !  resume previous search
      if( be%id > 0 ) then
		 from_above = .false.
	  else
	  	 write(0,*)'ebt_query_pair_list: BE is not a leaf'
		 stop
	  endif

   endif

   nx    = ebt_a%sell%nx
   ngeom = ebt_a%sell%ngeom

   be_loop: do  !  loop over BEs  

!  At the start of the loop,  be  and from_above must be set.
!  If  from_above==.true., BE has been entered from above, and therefore for the first time.
!  If  from_above==.false., then move from  be  to  be%parent, which has been entered before,
!     at which time, be%k_first and be%done(1:2) had been set.

   if( from_above ) then  !  have entered BE from above (and therefore for the first time) 

      call sell_ell_pt_in( nx, ngeom, be%geom, x_a, in )
	  be_tests = be_tests + 1

      if( be%id > 0 ) then  ! BE is a leaf
         a_tests = a_tests + 1
	     be_last_leaf => be

	     if( in ) then      ! in ELL_A: test for x_f covered by ELL_F
	        ell => sell_f%ell_pt(be%id)%ell
			call sell_ell_pt_in( sell_f%nx, sell_f%ngeom, ell%geom, x_f, in )

	        f_tests = f_tests + 1
			ebt_a%prop(5) = ebt_a%prop(5) + 1.d0  ! leaf tests
			if( .not.in ) ebt_a%prop(6) = ebt_a%prop(6) + 1.d0  !  false positives
	     endif

	     if( in ) then
		    n_id     = n_id + 1
			id(n_id) = be%id
			if( n_id >= n_id_max  .or. f_tests >= max_f_tests  .or.  be_tests >= max_be_tests  ) exit be_loop

		    from_above = .false. !  move back up to parent
		    cycle be_loop
		 endif

		 if( f_tests >= max_f_tests  .or.  be_tests >= max_be_tests ) exit be_loop
      endif 

	  if( .not.in ) then      !  x is not in BE
		 from_above = .false. !  move up to parent
		 cycle be_loop
	  endif

!  x is in BE which is not a leaf
!  set be%done(k) = .false. if child k exists, and x lies between its supporting hyperplanes

      be%done  = .true.             !  assume x not in children until found otherwise

      s = sum( be%v * (x_a-be%p) )  !  compute s
	  do k = 1, 2                   !  x_a may be in child if between supporting hyperplanes
		if( s >= be%smin(k)  .and.  s <= be%smax(k) ) be%done(k) = .false.
	  end do

	  if( ebt_a%force_search >= 1 ) be%done = .false.

	  if( s < 0.d0 ) then  !  select child to be tested first, based on s
	     be%k_first = 1
	  else
	     be%k_first = 2
	  endif

   else   !  .not.from_above:  move from be to be%parent

      if( .not.associated(be%parent) ) exit be_loop   !  BE is root

	  be%parent%done( be%k ) = .true. !  BE has been processed
	  be => be%parent                 !  moved to parent

   endif

!  x_a is in BE which is a node:  determine whether to process a child or go to parent 

   ka = be%k_first
   kb = 3 - ka
   if( .not.be%done(ka) ) then      !  go down to child ka
      from_above = .true.
	  be => be%child(ka)%be

   elseif( .not.be%done(kb) ) then  !  go down to child kb
      from_above = .true.
	  be => be%child(kb)%be

   else                             !  both children processed: go to parent
	  from_above = .false.
   endif

   end do be_loop 

   if( associated(be_last_leaf) ) then  !  set BE to last leaf encountered (or nullify)
      be => be_last_leaf
   else
      nullify( be )
   endif

   return

end subroutine ebt_query_pair_list  !-------------------------------------------------------------------

subroutine ebt_query_which( be, s, mode, k )  !  decide of which child to visit first in search

   implicit none
   type (be_type), pointer :: be   ! node BE
   real(k_d), intent(in)   :: s    ! distance from cutting plane
   integer, intent(in)     :: mode ! criterion to be used (see below)
   integer, intent(out)    :: k    ! child to visit first (1 or 2)

!  The child (k=1 or k=2) to visit first is decided by the following criteria, 
!  depending on the value of mode.  For 1 <= mode <= 8:
 
!  1 - k=1
!  2 - k=1 for s<0, k=2 for s>=0
!  3 - k which maximizes min( |s-smin(k)|, |s-smax(k)| )
!  4 - k which maximizes the number of leaves in the sub-tree be%child(k)
!  5 - k which maximizes the content of leaves in the sub-tree be%child(k)
!  6 - k which maximizes the "length" of leaves in the sub-tree be%child(k)
!  7 - k such that be%child(k) was the path of the last successful search throuh BE
!  8 - k such that be%child(k) is the most successful direction in previous searches

!  For mode < 0, 1 <= |mode| <= 8, the opposite of mode=|mode| is used.

   integer :: m

   m = abs( mode )  !  Note: for efficiency, m is not checked
   k = 2  !  assume k=2, reset to k=1 as necessary

   if( m <= 4 ) then
      if( m <= 2 ) then
	     if( m ==1 ) then  ! m=1
            k=1
		 else              ! m=2
            if( s < 0.d0 ) k=1
		 endif
	  else  ! m > 2,  m <= 4
	     if( m ==3 ) then  ! m=3
		    if( min( s-be%smin(1) ,  be%smax(1)-s ) >  &
                min( s-be%smin(2) ,  be%smax(2)-s )  )  k=1
		 else              ! m=4
            if( be%child(1)%be%prop(1) > be%child(2)%be%prop(1) ) k=1
		 endif
	  endif
   else ! m > 4
      if( m <= 6 ) then
	     if( m ==5 ) then  ! m=5
            if( be%child(1)%be%prop(2) > be%child(2)%be%prop(2) ) k=1
		 else              ! m=6
            if( be%child(1)%be%prop(3) > be%child(2)%be%prop(3) ) k=1
		 endif
	  else  ! m > 6
	     if( m ==7 ) then  ! m=7
            if( be%prop(4) < 1.5d0 ) k=1
		 else              ! m=8
            if( be%prop(5) > be%prop(6) ) k=1
		 endif
	  endif
   endif

   if( mode < 0 ) k = 3 - k

   return
end subroutine ebt_query_which  !-------------------------------------------------------------------

subroutine ebt_query_success( ebt, be0 )  !  update node props of EBT based on BE0

!  Update the properties be%prop(4:6) of all antecedent nodes of BE0  in EBT.

   implicit none

   type (ebt_type), pointer :: ebt         ! given EBT
   type (be_type), pointer  :: be0         ! starting BE (leaf or node)

   integer   :: k
   type (be_type), pointer :: be

   if( .not.associated(ebt) ) then  !  check input
      write(0,*) 'ebt_query_success: ebt not associated; stopping'
	  stop
   endif

   if( .not.associated(be0) ) then 
      write(0,*) 'ebt_query_success: be0 not associated; stopping'
	  stop
   endif

   be => be0
   do  !  loop over antecedent nodes

   	  if( .not. associated( be%parent ) ) then  !  be is the root: all done
         nullify( be )
		 return
      endif

      k  =  be%k
	  be => be%parent  !  go to parent

	  be%prop(4)   = k
	  be%prop(4+k) = be%prop(4+k) + 1.d0

   end do

   nullify( be )
   return
end subroutine ebt_query_success  !------------------------------------------------------------------

subroutine ebt_status( ebt, status, leaves, max_depth )  !  check EBT for problems
   implicit none

   type (ebt_type), pointer :: ebt       ! EBT
   integer, intent(out)     :: status    ! =0 for a consistent rooted tree
   integer, intent(out)     :: leaves    ! number of leaves
   integer, intent(out)     :: max_depth ! maximum depth

!  Note: for status/=0, leaves and max_depth are not accurate

   type (be_type), pointer :: be
   integer :: i, depth

   status    = 0  !  assume consistent until found otherwise
   leaves    = 0
   max_depth = 0

   if( .not.associated(ebt) ) then
      status = -1  !  no EBT
	  return
   endif

   if( .not.associated(ebt%sell) ) then
      status = -2  !  no SELL
	  return
   endif

   if( .not.associated(ebt%root) ) then
      status = -3  !  no root
	  return
   endif

   do i = 1, ebt%n_pt  !  loop over leaves

      if( .not.associated( ebt%be_pt(i)%be ) ) cycle
	  leaves = leaves + 1
!  check leaf properties
      be => ebt%be_pt(i)%be
      if( .not.associated( be%ebt, ebt ) )then
	     status = -4
		 return
	  elseif( associated( be%child ) ) then
         status = -5
		 return
      elseif( be%id /= i ) then
         status = -6
		 return
      endif

!  loop over antecedent nodes
      depth = 1
	  do
	     if( .not.associated (be%parent) ) then  ! root
		    if( .not.associated(be,ebt%root) ) then
			   status = -7
			   return
			endif

            exit  !  root reached: proceed to next leaf

		 else  ! parent exists
		    if( .not.associated( be, be%parent%child( be%k )%be ) ) then
			   status = -8
			   return
		    endif
		 endif

		 be => be%parent
		 depth = depth + 1

		 if( .not.associated( be%ebt, ebt ) )then
	        status = -9
		    return
	     elseif( .not.associated( be%child(1)%be )  .or.  .not.associated( be%child(2)%be ) ) then
            status = -10
		    return
         elseif( be%id > 0 ) then
            status = -11
		    return
         endif
	  end do
	  max_depth = max( depth, max_depth )

   end do

end subroutine ebt_status  !-------------------------------------------------------------------

subroutine ebt_param_set( ebt, check, force_search, pair_sep_qual, det_lim ) 

! set parameters

   type (ebt_type), pointer        :: ebt    ! data structure for EBT. 
   integer, intent(in), optional   :: check  ! level of checking [0,1,2]
   integer, intent(in), optional   :: force_search  
   real(k_d), intent(in), optional :: pair_sep_qual  !  quality parameter (less than 1.0) 
   real(k_d), intent(in), optional :: det_lim  !  min. det.


   if( .not.associated(ebt) ) then
      write(0,*) 'ebt_param_set: ebt not associated; stopping'
	  stop
   endif

   if( present(check) ) then
      if( check < 0  .or.  check >2 ) then
	     write(0,*) 'ebt_param_set: invalid  check = ', check
	     stop
	  endif

      ebt%check = check
   endif

   if( present(force_search) )  ebt%force_search  = force_search
   if( present(pair_sep_qual) ) ebt%pair_sep_qual = pair_sep_qual
   if( present(det_lim) )       ebt%det_lim       = det_lim

   return
end subroutine ebt_param_set  !--------------------------------------------------------------

subroutine ebt_param_get( ebt, check, force_search, pair_sep_qual, det_lim ) 

! get parameters

   type (ebt_type), pointer         :: ebt   ! data structure for EBT. 
   integer, intent(out), optional   :: check  ! level of checking 
   integer, intent(out), optional   :: force_search  
   real(k_d), intent(out), optional :: pair_sep_qual  !  quality parameter (less than 1.0) 
   real(k_d), intent(out), optional :: det_lim  


   if( .not.associated(ebt) ) then
      write(0,*) 'ebt_param_get: ebt not associated; stopping'
	  stop
   endif

   if( present(check) ) check = ebt%check
   if( present(force_search) )  force_search  = ebt%force_search
   if( present(pair_sep_qual) ) pair_sep_qual = ebt%pair_sep_qual
   if( present(det_lim) )       det_lim       = ebt%det_lim

   return
end subroutine ebt_param_get  !--------------------------------------------------------------

subroutine ebt_prop( ebt, prop ) 

! get properties of EBT

   type (ebt_type), pointer :: ebt       ! data structure for EBT.   
   real(k_d), intent(out)   :: prop(10)  ! properties

! 1 - number of leaves in tree
! 2 - sum of content (volume) of leaves in tree
! 3 - sum of radii of circumscribed balls of leaves in tree
! 4 - det_min - minimum determinant of G encountered
! 5 - number of leaf tests in ebt_query_pair routines
! 6 - number of false positives in ebt_query_pair routines

   prop = 0.d0

   if( .not.associated(ebt) ) then
      write(0,*) 'ebt_prop: ebt not associated; stopping'
	  stop
   endif

   if( associated(ebt%root) ) prop(1:3) = ebt%root%prop(1:3)
   prop(4:6) = ebt%prop(4:6)

   return
end subroutine ebt_prop  !--------------------------------------------------------------

subroutine ebt_write( ebt, lu, n_leaves )  !  write EBT to  logical unit LU

   type (ebt_type), pointer  :: ebt       !  EBT to be written
   integer, intent(in)       :: lu        !  LU: must be open for writing
   integer, intent(out)      :: n_leaves  ! number of leaves written

   type (be_type), pointer   :: be
   logical                   :: opened, from_above
   integer                   :: i, id, id_parent

   if( .not.associated(ebt) ) then 
      write(0,*)'ebt_write: ebt not associated'
	  stop
   endif

   inquire( lu, opened=opened )
   if( .not.opened ) then
      write(0,*)'ebt_write: lu not opened ', lu
	  stop
   endif

   write(lu, err=100 ) ebt%n_pt, ebt%check, ebt%prop, ebt%force_search, &
                       ebt%pair_sep_qual, ebt%det_lim

   n_leaves = 0
   id       = 0
   do i = 1, ebt%n_pt  !  write leaf BEs
      if( associated(ebt%be_pt(i)%be) ) then
	     n_leaves = n_leaves + 1
		 id       = ebt%be_pt(i)%be%id
         write(lu, err=110 ) id
		 call be_write( ebt%be_pt(i)%be, lu )
	  else
	     write(lu) 0
	  endif
   end do

   if( n_leaves < 2 ) return  !  all done if no nodes

!  first traverse of EBT to assign unique negative id's to nodes, and to write node BEs
!  note that the root node has id = -1

   id         = 0
   from_above = .true.
   be         => ebt%root

   do
      if( from_above ) then  ! entering node or leaf for the first time

         if( be%id < 0 ) then   !  first time at node
	        id = id - 1
	        be%id = id          ! assign ID and write node BE
	        call be_write( be, lu )

			be%done(2) = .false.
		    be => be%child(1)%be  !  go to child 1, cycle
         else  !  leaf:  do nothing; return to parent
			from_above = .false.
			if( .not.associated(be%parent) ) exit
			be => be%parent  !  go to parent, cycle
		 endif

	  else ! from_below
	     if( .not.be%done(2) ) then
		    be%done(2) = .true.
			from_above = .true.
		    be => be%child(2)%be  !  go to child 2, cycle
         else  ! go to parent
			from_above = .false.
			if( .not.associated(be%parent) ) exit
			be => be%parent  !  go to parent, cycle
		 endif
	  endif
   end do

   if( -id /= n_leaves-1 ) then
      write(0,*)'ebt_write: node miscount, ', -id, n_leaves-1
	  stop
   endif

! second traverse of EBT nodes to write ID of BE, parent and children
   from_above = .true.
   be         => ebt%root

   do
      if( from_above ) then  !first time at node
	     if( associated(be%parent) ) then
		    id_parent = be%parent%id
		 else
		    id_parent = 0
		 endif

		 write(lu, err=120 ) be%id, id_parent, be%child(1)%be%id, be%child(2)%be%id

		 be%done(2) = .false.
		 if( be%child(2)%be%id > 0 ) be%done(2) = .true.

		 if( be%child(1)%be%id < 0 ) then
		    be         => be%child(1)%be
		 elseif( be%child(2)%be%id < 0 ) then
		    be%done(2) = .true.
		    be         => be%child(2)%be
		 else
		    be%done = .true.
		    from_above = .false.
			if( .not.associated(be%parent) ) exit
			be => be%parent
		 endif

      else !  from below

         if( .not.be%done(2) ) then
		    from_above = .true.
			be%done(2) = .true.
		    be => be%child(2)%be
		 else
			if( .not.associated(be%parent) ) exit
			be => be%parent
		 endif

      endif
   end do

   return

100   write(0,*)'ebt_write: error writing n_pt'
      stop

110   write(0,*)'ebt_write: error writing ebt%be_pt(i)%be%id'
      stop

120   write(0,*)'ebt_write: error writing be%id'
      stop


end subroutine ebt_write    !-----------------------------------------------

subroutine ebt_read( sell, ebt, idlist, lu, n_leaves )  !  read EBT from logical unit LU

!  EBT should not be allocated prior to call

   type (sell_type), pointer    :: sell
   type (ebt_type), pointer     :: ebt
   type (id_list_type), pointer :: idlist 
   integer, intent(in)          :: lu
   integer, intent(out)         :: n_leaves  ! number of leaves read

   type (be_type), pointer      :: be
   type (be_pointer), pointer   :: node_pt(:)
   integer   :: n_pt, check, force_search
   real(k_d) :: prop(10), pair_sep_qual, det_lim

   integer            :: i, id, n_nodes, id_node, id_parent, id_child1, id_child2
   logical            :: opened, from_above

   inquire( lu, opened=opened )
   if( .not.opened ) then
      write(0,*)'ebt_write: lu not opened ', lu
	  stop
   endif

   read(lu, err=100) n_pt, check, prop, force_search, pair_sep_qual, det_lim

   call ebt_initialize( sell, ebt, n_pt=n_pt, check=check, idlist=idlist, force_search=force_search, &
                        pair_sep_qual=pair_sep_qual, det_lim=det_lim )

   ebt%prop = prop

   n_leaves = 0
   do i = 1, ebt%n_pt  !  read leaf BEs
      read(lu, err=110 ) id
      if( id > 0 ) then
	     n_leaves = n_leaves + 1
		 call be_read( ebt, be, lu )
		 if( be%id /= id ) then
		    write(0,*)'ebt_read: id mismatch ', be%id, id
			stop
		 endif
         be%geom => sell%ell_pt(id)%ell%geom
		 ebt%be_pt(i)%be => be
	  else
	     nullify( ebt%be_pt(i)%be )
	  endif
   end do

   if( n_leaves < 2 ) then
      if( n_leaves == 1 ) ebt%root => be
	  return  !  all done if no nodes
   endif

   n_nodes = n_leaves - 1
   allocate( node_pt(n_nodes) )

   do id = 1, n_nodes  !  read node BEs
      call be_read( ebt, be, lu )
	  if( be%id /= -id ) then
	     write(0,*)'ebt_read: node ID mismatch ', be%id, -id
		 stop
	  endif
	  node_pt(id)%be => be
   end do

   ebt%root => node_pt(1)%be

!  traverse EBT to assign child(1:2) and parent

   from_above = .true.
   be         => ebt%root

   do
      if( from_above ) then  ! entering node or leaf for the first time
	     read(lu, err=120 ) id_node, id_parent, id_child1, id_child2

		 be => node_pt(-id_node)%be

		 if( id_parent < 0 ) then  !  set parent
		    be%parent => node_pt(-id_parent)%be
		 else
		    nullify( be%parent )
		 endif

         be%done = .false.
		 if( id_child1 < 0 ) then  !  set child 1
		    be%child(1)%be => node_pt(-id_child1)%be
		 else
            be%child(1)%be => ebt%be_pt(id_child1)%be 
			ebt%be_pt(id_child1)%be%parent => be
			be%done(1)     =  .true.
		 endif

		 if( id_child2 < 0 ) then  !  set child 2
		    be%child(2)%be => node_pt(-id_child2)%be
		 else
            be%child(2)%be => ebt%be_pt(id_child2)%be 
			ebt%be_pt(id_child2)%be%parent => be
			be%done(2)     =  .true.
		 endif

		 if( .not.be%done(1) ) then
		    be%done(1) = .true.
		    be         => be%child(1)%be
		 elseif( .not.be%done(2) ) then
		    be%done(2) = .true.
		    be         => be%child(2)%be
		 else
		    from_above = .false.
			if( .not.associated(be%parent) ) exit
			be => be%parent
		 endif

      else !  from below

         if( .not.be%done(2) ) then
		 	be%done(2) = .true.
		    from_above = .true.
		    be => be%child(2)%be
		 else
			if( .not.associated(be%parent) ) exit
			be => be%parent
		 endif

      endif
   end do

   deallocate( node_pt )
	    
   return

100   write(0,*)'ebt_read:  error reading n_pt'
      stop

110   write(0,*)'ebt_read:  error reading id'
      stop

120   write(0,*)'ebt_read:  error reading id_node'
      stop

end subroutine ebt_read    !-----------------------------------------------

subroutine ebt_be_pt_double( ebt ) 

!  double size of array ebt%be_pt

      type (ebt_type),   pointer :: ebt

      type (be_pointer), pointer :: pt_temp(:)
	  integer                    :: id

!  allocate temporary array pt_temp of twice size
   	  allocate( pt_temp(2*ebt%n_pt) )  

	  do id = 1, 2*ebt%n_pt
	     nullify( pt_temp(id)%be )
      end do

!  copy ell%pt to pt_temp
	  do id = 1, ebt%n_pt
	     if( associated( ebt%be_pt(id)%be ) ) pt_temp(id)%be => ebt%be_pt(id)%be
		 nullify( ebt%be_pt(id)%be )
	  end do

!  double size of ebt%be_pt
	  deallocate( ebt%be_pt )

      ebt%n_pt = 2 * ebt%n_pt
	  allocate( ebt%be_pt(ebt%n_pt) )  

	  do id = 1, ebt%n_pt
	     
	     if( associated( pt_temp(id)%be ) ) then
		    ebt%be_pt(id)%be => pt_temp(id)%be
	     else
		    nullify( ebt%be_pt(id)%be )
		 endif

		 nullify(pt_temp(id)%be) 
	  end do

	  deallocate( pt_temp )

	  return

end subroutine ebt_be_pt_double   !-----------------------------------------------

end module sell_ebt
