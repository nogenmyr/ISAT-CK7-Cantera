!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module sell_bt

!  Module for searching a set of ellipsoids (SELLs) using a binary tree (BT).

!  This code is based on sell_ebt, where nodes are referred to as BEs.
!  Here, to avoid name conflicts, BE is changed to ND (for node).

!  A BT has nodes and leaves, with an ND associated with each.  
!  The ND associated with a leaf is an ellipsoid (ELL) in the SELL.
!  Each node has two children (1 and 2).   

!  S.B. Pope  6/3/2006

!  Data structures
!     bt_type    - for BTs
!     nd_type    - for bounding ellipsoids (ND)
!     nd_pointer - ponter to ND

!================================================

! Public soubroutines

! bt_initialize  ! allocate and initialize BT
! bt_destroy     ! destroy and deallocate BT
! bt_add         ! add leaf ELL to BT
! bt_remove      ! remove leaf ELL from BT
! bt_purge       ! remove all nodes (NDs) from BT
! bt_rebuild     ! rebuild BT
! bt_traverse    ! traverse BT based on query point x
! bt_ball_pts    ! find points in ball 

! bt_param_set   ! set control parameters
! bt_param_get   ! get control parameters
! bt_prop        ! get BT properties
! bt_write       ! write BT to  logical unit LU (for checkpointing)
! bt_read        ! read BT from logical unit LU

! Private subroutines

! nd_create      ! allocate and initialize ND to null settings
! nd_destroy     ! deallocate ND
! nd_write       ! write ND to logical unit LU
! nd_read        ! create and read ND from logical unit LU
! bt_nd_pt_double! double pointer length of array

! bt_initialize( sell, bt, n_pt, check, idlist )
! bt_destroy( bt )
! bt_add( bt, id ) 
! bt_remove( bt, id )  
! bt_purge( bt )   
! bt_rebuild( bt, idlist, mode )  
! bt_traverse( bt, x, id ) 
! bt_ball_pts( bt, c, r, npts_max, npts, ids, cp_tests, in_tests ) 
! bt_param_set( bt, check ) 
! bt_param_get( bt, check ) 
! bt_prop( bt, prop ) 
! bt_write( bt, lu, n_leaves )  
! bt_read( sell, bt, idlist, lu, n_leaves ) 
! nd_create( bt, node_leaf, nd )  
! nd_destroy( nd ) 
! nd_write( nd, lu ) 
! nd_read( bt, nd, lu ) 
! bt_nd_pt_double( bt ) 

!================================================

use sell_m
implicit none

   private nd_create, nd_destroy, nd_write, nd_read, bt_nd_pt_double

   type :: bt_type  !  Binary Tree  (BT) !--------
      type (sell_type), pointer    :: sell     ! SELL to which BT belongs
      type(nd_type),    pointer    :: root     ! root ND of BT
	  type(nd_pointer), pointer    :: nd_pt(:) ! pointer to leaves: nd_pt(id) 
	  integer                      :: n_pt     ! dimension of nd_pt

      integer                      :: check    ! level of checking [0,1,2]
      type (id_list_type), pointer :: idlist   ! data structure used in id_list to generate IDs

      real(k_d)                    :: prop(10) ! properties of BT - see bt_prop
	  integer                      :: depth    ! depth of last traverse
   end type bt_type

   type :: nd_type   !  Bounding Ellipsoid (ND) ----------------
      type (bt_type),    pointer :: bt        ! BT to which ND belongs
      type (nd_type),    pointer :: parent    ! nd%parent  is the parent of ND
      type (nd_pointer), pointer :: child(:)  ! nd%child(k)%nd  is the k-th child of ND
	  integer                    :: id        ! nd%id is the ID of ND's leaf (if ND is a leaf of the bt)
	                                          ! nd%id < 0 indicates that ND is a node
	  integer                    :: k         ! ND is parent's k-th child (k = 1 or 2)

	  real(k_d), pointer         :: p(:)      ! point on cutting hyperplane
	  real(k_d), pointer         :: v(:)      ! cutting plane vector:  s(x) = v^T * (x-p)
	  integer                    :: k_first   ! first child to be processed
	  logical                    :: done(2)   ! flag for children having been treated
   end type nd_type  !------------------------------------------

   type :: nd_pointer  
      type (nd_type), pointer :: nd
   end type nd_pointer

contains  !===================================================================

subroutine bt_initialize( sell, bt, n_pt, check, idlist )

!  Allocate and initialize (empty) BT.  Optionally set check.

   type (sell_type), pointer              :: sell    ! (in) SELL on which bt is to be formed
   type (bt_type),  pointer               :: bt     ! (out) data structure for BT

!  Optional
   integer, intent(in), optional          :: n_pt    ! initial dimension of pointer array
   type (id_list_type), pointer, optional :: idlist  ! (in) data structure used in id_list to generate IDs

   integer, intent(in), optional          :: check   ! level of checking  =0 minimal; =1 some; =2 maximal

   integer :: i

   allocate( bt )
   bt%sell => sell
   nullify(bt%root)

   bt%n_pt = n_pt_0
   if( present(n_pt) ) bt%n_pt = n_pt
   allocate( bt%nd_pt( bt%n_pt ) )

   do i = 1, bt%n_pt
      nullify( bt%nd_pt(i)%nd )
   end do

   if( present(idlist) ) then
      bt%idlist => idlist
   else
      nullify( bt%idlist)
   endif

   if( present(check) ) then
      call bt_param_set( bt, check )
   else
      bt%check = 0
   endif

   bt%prop  = 0.d0
   bt%depth = 0

   return
end subroutine bt_initialize  !------------------------------------------------------------

subroutine bt_destroy( bt )

!  Destroy and deallocate BT bt.
!  All NDs should first be destroyed by a call to bt_purge

   type (bt_type), pointer :: bt    ! (in) data structure for BT

   if( .not.associated(bt) ) then
      write(0,*) 'bt_destroy: bt not associated; stopping'
	  stop
   endif

   deallocate( bt%nd_pt )
   deallocate( bt )

   return
end subroutine bt_destroy  !------------------------------------------------------------

subroutine bt_add( bt, id )  !  add ELL (ID=id) to BT

! BT is traversed to identify ND_SIBLING.
! A new ND node is added (in place of ND_SIBLING) whose 1st and 2nd 
! children are ND and ND_SIBLING, respectively.

   implicit none

   type (bt_type), pointer  :: bt         ! BT
   integer, intent(in)      :: id         ! ID of ELL

   type (nd_type), pointer  :: nd, nd_parent, nd_sibling
   integer                  :: id_status, ids_assigned, id_max, nx, ngeom, &
                               id_sibling, ksp, depth

   if( .not.associated(bt) ) then  !  check input
      write(0,*) 'bt_add: bt not associated; stopping'
	  stop
   endif

   if( bt%check > 1  .and.  associated( bt%idlist ) ) then
	  call id_query( bt%idlist, id, id_status, ids_assigned, id_max )
	  if( id_status /= 1 ) then
		 write(0,*) 'bt_add: id has not been assigned, id = ', id
	     stop
	  endif
   endif

   do while ( id > bt%n_pt )
      call bt_nd_pt_double( bt ) 
   end do

   nx    = bt%sell%nx
   ngeom = bt%sell%ngeom

!  initialize ND for new leaf
   call nd_create( bt, 2, nd )
   bt%prop(1) = bt%prop(1) + 1.d0

   nd%id  =  id
   bt%nd_pt(id)%nd => nd

! traverse bt to leaf = nd_sibling
   call bt_traverse( bt, bt%sell%ell_pt(id)%ell%geom(1:nx), id_sibling, depth )

! special treatment for first ND in BT
   if( id_sibling == 0 ) then
      bt%root => nd   !  ND is root of bt
	  return
   endif

   nd_sibling => bt%nd_pt(id_sibling)%nd

   bt%prop(2) = max( bt%prop(2), bt%depth+1.d0 )  !  update maximum depth

!  set  nd  to 1st child; nd_sibling to second child of nd_parent

!  allocate and initialize ND = nd_parent for new node
   call nd_create( bt, 1, nd_parent )

   if( associated( nd_sibling%parent ) ) then  !  nd_sibling was not root
      nd_parent%parent => nd_sibling%parent
	  ksp              =  nd_sibling%k
      nd_parent%k      =  ksp
	  nd_sibling%parent%child(ksp)%nd => nd_parent
   else                          !  nd_sibling was root
	  bt%root => nd_parent       !  nd_parent is now root
   endif

   nd_parent%child(1)%nd => nd
   nd_parent%child(2)%nd => nd_sibling

!  re-set nd_sibling

   nd_sibling%parent => nd_parent
   nd_sibling%k      =  2

!  set  nd

   nd%parent => nd_parent
   nd%k      =  1

!  set nd_parent's cutting plane (based on nd_sibling's ELL and nd's center)

   call ell_pt_hyper( nx, bt%sell%ell_pt(id_sibling)%ell%geom(1:nx), &
                          bt%sell%ell_pt(id_sibling)%ell%geom(nx+1:ngeom-2), &
                          bt%sell%ell_pt(id)%ell%geom(1:nx), nd_parent%p, nd_parent%v )

   nd_parent%v = - nd_parent%v ! s < 0 indicates parent%child(1) => nd

   return
end subroutine bt_add  !----------------------------------------------------------

subroutine bt_remove( bt, id )  !  remove leaf with ID=id from BT
   implicit none

   type (bt_type), pointer :: bt ! given BT
   integer, intent(in)     :: id ! ID of leaf to be removed

   integer                   :: k, ks, id_status, ids_assigned, id_max
   type (nd_type), pointer   :: nd, nd_sibling, nd_parent, nd_grandparent 

   if( .not.associated(bt) ) then  !  check input
      write(0,*) 'bt_remove: bt not associated; stopping'
	  stop
   endif

   if( id < 1 ) then  !  check that id corresponds to an extant leaf
      write(0,*)'bt_remove: id is non-positive'
	  stop
   elseif( bt%check > 1 ) then
      call id_query( bt%idlist, id, id_status, ids_assigned, id_max )
	  if( id_status /= 1 ) then
	     write(0,*)'bt_remove: id is not assigned ', id, id_status, ids_assigned, id_max
	     stop
	  endif
   endif

   nd         => bt%nd_pt(id)%nd   !  leaf to be removed
   bt%prop(1) = bt%prop(1) - 1.d0

   if( .not.associated(nd%parent) ) then  !  ND is root of bt
      nullify( bt%root )                  !  BT is empty
	  call nd_destroy( nd )
	  nullify( bt%nd_pt( id )%nd )
	  return
   endif

   nd_parent  => nd%parent     ! identify ND's parent and sibling
   k          =  nd%k 
   ks         =  3 - k
   nd_sibling => nd_parent%child(ks)%nd

   if( .not.associated(nd_parent%parent) ) then  !  nd_parent is root of BT
      nd%bt%root => nd_sibling                   !  nd's sibling is sole surviving ND
	  nd_sibling%k = 0
	  nullify( nd_sibling%parent )

   else  !  grandparent exists

      nd_grandparent    => nd_parent%parent   !  new parent of sibling
      nd_sibling%parent => nd_grandparent
	  k                 = nd_parent%k         !  index of parent's child (old and new)
      nd_sibling%k      = k
      nd_grandparent%child(k)%nd => nd_sibling
   endif
  
   call nd_destroy( nd )
   call nd_destroy( nd_parent )
   nullify( bt%nd_pt( id )%nd )

   return
end subroutine bt_remove   !-------------------------------------------------------

subroutine bt_purge( bt )   !  remove all NDs from BT

   implicit none
   type (bt_type), pointer  :: bt   ! given BT

   type (nd_type), pointer  :: nd, parent
   integer                  :: i

   if( .not.associated(bt) ) then  !  check input
      write(0,*) 'bt_purge: bt not associated; stopping'
	  stop
   endif

!  start at root (if it exists)

   if( .not.(associated(bt%root) ) ) return  !  already purged

   nd => bt%root  !  root
   nullify( parent )            !  ...has no parent

   do   !  on each iteration of loop,  ND  and  parent  must be defined (parent may be null)

      if( nd%id < 1  ) then  !  ND is a node
      
         if( associated(nd%child(1)%nd) ) then      !  go to child(1) if it exists
	        parent => nd
		    nd     => nd%child(1)%nd
		    nullify( parent%child(1)%nd )
		    cycle

         elseif( associated(nd%child(2)%nd) ) then  !  go to child(2) if it exists
	        parent => nd
		    nd     => nd%child(2)%nd
		    nullify( parent%child(2)%nd )
		    cycle
		 endif
		 
	  endif
	  
!  ND  is childless: destroy and go to parent (if it exists)
	     
	  call nd_destroy( nd )                   !  destroy  nd
	  if( .not.associated( parent ) ) exit    !  all done
	  nd => parent                            !  go to parent

	  if( associated(nd%parent) ) then        !  identify parent (if any)
		 parent => nd%parent
	  else
		 nullify( parent )
	  endif

   end do

   nullify( bt%root )

   do i = 1, bt%n_pt
      nullify( bt%nd_pt(i)%nd )
   end do

   bt%prop = 0.d0  ! re-initialize bt%prop

   return

end subroutine bt_purge  !------------------------------------------------------------------

subroutine bt_rebuild( bt, idlist, mode )   !  rebuild BT

   implicit none
   type (bt_type), pointer      :: bt     !  BT to be rebuilt
   type (id_list_type), pointer :: idlist !  IDLIST
   integer, intent(in)          :: mode   !  method of rebuilding

! mode = 1 - purge BT, and then re-add ELLs in a random order
! (only mode currently supported)

   integer              :: i, ids_assigned, id_max
   integer, allocatable :: id_array(:)

   if( .not.associated(bt) ) then  !  check input
      write(0,*) 'bt_rebuild: bt not associated; stopping'
	  stop
   endif

   call id_stats( idlist, ids_assigned, id_max )
   if( ids_assigned < 3 ) return  !  unique configuration

   if( mode == 1 ) then  !  purge and re-add ELLs in random order

	  allocate( id_array(ids_assigned) )
	  call id_rand_array( idlist, id_array )
      call bt_purge( bt ) 

	  do i = 1, ids_assigned
	     call bt_add( bt, id_array(i) )
	  end do

	  deallocate( id_array )

   else
      write(0,*)'bt_rebuild: invalid mode = ', mode
	  stop
   endif

   return
end subroutine bt_rebuild  !------------------------------------------------------------------

subroutine bt_traverse( bt, x, id, depth )  !  traverse BT based on x

!  Traverse BT based on x to indentify leaf with ID=id.
!  If the root of the BT does not exist, id=0 is returned.
!  bt%depth is set to the depth of the traverse.

   implicit none

   type (bt_type), pointer   :: bt            ! given BT
   real(k_d), intent(in)     :: x(bt%sell%nx) ! query point
   integer, intent(out)      :: id            ! ID of leaf
   integer, intent(out)      :: depth         ! depth of traverse path

   type (nd_type), pointer   :: nd
   real(k_d)                 :: s

   if( .not.associated( bt ) ) then  ! check that bt exists
      write(0,*)'bt_traverse; bt not associated'
	  stop
   endif

   if( .not.associated( bt%root) ) then 
      id       = 0 !  indicate that root does not exist
	  depth    = 0
	  bt%depth = 0
	  return
   endif

   nd         => bt%root            !  start traverse at root
   bt%prop(3) =  bt%prop(3) + 1.d0  !  increment number of traverses
   depth      =  1

   do  !  loop over NDs  ------------------------------------------

      if( nd%id > 0 ) then  !  leaf has been reached - all done
	     id         = nd%id
	     bt%prop(4) = bt%prop(4) + depth  !  increment sum of depths
		 bt%depth   = depth
	     return  
      endif

	  s = sum( nd%v * (x-nd%p) )  !  compute s
      if( s < 0.d0 ) then   !  go to child 1 or 2 depending on s
	     nd => nd%child(1)%nd
	  else
	     nd => nd%child(2)%nd
	  endif

	  depth = depth + 1

   end do

end subroutine bt_traverse  !---------------------------------------------------------------

subroutine bt_ball_pts( bt, c, r, npts_max, npts, ids, dsq, cp_tests, in_tests ) 

!  Given a ball, centered at  c  of radius r, traverse the binary tree BT
!  to identify up to npts_max points covered by the ball. 


   implicit none

   type (bt_type), pointer   :: bt            ! given BT
   real(k_d), intent(in)     :: c(bt%sell%nx) ! center of ball
   real(k_d), intent(in)     :: r             ! radius of ball
   integer, intent(in)       :: npts_max      ! max. no. of points to be returned
   integer, intent(out)      :: npts          ! no. of points found
   integer, intent(out)      :: ids(npts_max) ! IDs of points
   real(k_d), intent(out)    :: dsq(npts_max) ! square of distance of point from center
   integer, intent(out)      :: cp_tests      ! number of cutting plane tests
   integer, intent(out)      :: in_tests      ! number of in-ball tests

   type (nd_type), pointer   :: nd
   real(k_d)                 :: s, x(bt%sell%nx), rsq, xsq
   logical                   :: from_above
   integer                   :: nx, id

   if( .not.associated( bt ) ) then  ! check that bt exists
      write(0,*)'bt_traverse; bt not associated'
	  stop
   endif

   npts     = 0  !  null initialization
   cp_tests = 0
   in_tests = 0

   if( npts_max < 1 ) return

   if( .not.associated( bt%root) ) return 

   rsq        = r*r
   nx         = bt%sell%nx
   nd         => bt%root     !  start traverse at root
   from_above = .true.

   do  !  loop over NDs  ------------------------------------------
      if( from_above ) then
      
	     if( nd%id > 0 ) then  !  ND is a leaf - test for in ball
	        id  = nd%id
		    x   = bt%sell%ell_pt(id)%ell%geom(1:nx) - c
		    xsq = sum( x*x )
            in_tests = in_tests + 1

            if( xsq <= rsq ) then
               npts = npts + 1  !  in ball
			   ids(npts) = id
			   dsq(npts) = xsq
			   if( npts == npts_max ) return
		    endif

            if( .not.associated( nd%parent ) ) return  ! BT consists of single leaf
			   
			nd => nd%parent  !  go to parent
			from_above = .false.
			cycle

		 else  ! node - for the first time: determine which children need to be examined
		    nd%done = .false.

	        s = sum( nd%v * (c-nd%p) ) !  compute signed distance to center of ball
			cp_tests = cp_tests + 1

			if( s > r ) then
			   nd%done(1) = .true.  ! no need to examine child 1
			     
			elseif( s < -r ) then
			   nd%done(2) = .true.  ! no need to examine child 2
			endif
		 endif

	  endif

! at node; node%done has been set: go to child or parent

      if( .not.nd%done(1) ) then
	     nd%done(1) = .true.
	     nd         => nd%child(1)%nd
		 from_above = .true.

	  elseif( .not.nd%done(2) ) then
	     nd%done(2) = .true.
	     nd         => nd%child(2)%nd
		 from_above = .true.

	  else
	     if( .not.associated( nd%parent ) ) return  ! at root: all done
	     nd => nd%parent  !  go to parent
		 from_above = .false.
	  endif

   end do

end subroutine bt_ball_pts  !---------------------------------------------------------------

subroutine bt_param_set( bt, check ) 

! set parameters

   type (bt_type), pointer         :: bt    ! data structure for BT. 
   integer, intent(in), optional   :: check  ! level of checking [0,1,2]

   if( .not.associated(bt) ) then
      write(0,*) 'bt_param_set: bt not associated; stopping'
	  stop
   endif

   if( present(check) ) then
      if( check < 0  .or.  check >2 ) then
	     write(0,*) 'bt_param_set: invalid  check = ', check
	     stop
	  endif

      bt%check = check
   endif

   return
end subroutine bt_param_set  !--------------------------------------------------------------

subroutine bt_param_get( bt, check ) 

! get parameters

   type (bt_type), pointer         :: bt   ! data structure for bt. 
   integer, intent(out), optional   :: check  ! level of checking 

   if( .not.associated(bt) ) then
      write(0,*) 'bt_param_get: bt not associated; stopping'
	  stop
   endif

   if( present(check) ) check = bt%check

   return
end subroutine bt_param_get  !--------------------------------------------------------------

subroutine bt_prop( bt, prop ) 

! get properties of BT

   type (bt_type), pointer :: bt       ! data structure for BT.   
   real(k_d), intent(out)  :: prop(10) ! properties

! 1 - number of leaves in tree
! 2 - maximum depth of tree
! 3 - number of traverses performed
! 4 - sum of depths of traverses
! 5 - average depth of traverses

   prop = 0.d0

   if( .not.associated(bt) ) then
      write(0,*) 'bt_prop: bt not associated; stopping'
	  stop
   endif

   prop(1:4) = bt%prop(1:4)
   prop(5)   = prop(4) / max( prop(3), 1.d0 )

   return
end subroutine bt_prop  !--------------------------------------------------------------

subroutine bt_write( bt, lu, n_leaves )  !  write bt to  logical unit LU

   type (bt_type), pointer   :: bt        !  BT to be written
   integer, intent(in)       :: lu        !  LU: must be open for writing
   integer, intent(out)      :: n_leaves  ! number of leaves written

   type (nd_type), pointer   :: nd
   logical                   :: opened, from_above
   integer                   :: i, id, id_parent

   if( .not.associated(bt) ) then 
      write(0,*)'bt_write: bt not associated'
	  stop
   endif

   inquire( lu, opened=opened )
   if( .not.opened ) then
      write(0,*)'bt_write: lu not opened ', lu
	  stop
   endif

   write(lu, err=100 ) bt%n_pt, bt%check, bt%prop

   n_leaves = 0
   id       = 0
   do i = 1, bt%n_pt  !  write leaf NDs
      if( associated(bt%nd_pt(i)%nd) ) then
	     n_leaves = n_leaves + 1
		 id       = bt%nd_pt(i)%nd%id
         write(lu, err=110 ) id
		 call nd_write( bt%nd_pt(i)%nd, lu )
	  else
	     write(lu) 0
	  endif
   end do

   if( n_leaves < 2 ) return  !  all done if no nodes

!  first traverse of bt to assign unique negative id's to nodes, and to write node NDs
!  note that the root node has id = -1

   id         = 0
   from_above = .true.
   nd         => bt%root

   do
      if( from_above ) then  ! entering node or leaf for the first time

         if( nd%id < 0 ) then   !  first time at node
	        id = id - 1
	        nd%id = id          ! assign ID and write node ND
	        call nd_write( nd, lu )

			nd%done(2) = .false.
		    nd => nd%child(1)%nd  !  go to child 1, cycle
         else  !  leaf:  do nothing; return to parent
			from_above = .false.
			if( .not.associated(nd%parent) ) exit
			nd => nd%parent  !  go to parent, cycle
		 endif

	  else ! from_below
	     if( .not.nd%done(2) ) then
		    nd%done(2) = .true.
			from_above = .true.
		    nd => nd%child(2)%nd  !  go to child 2, cycle
         else  ! go to parent
			from_above = .false.
			if( .not.associated(nd%parent) ) exit
			nd => nd%parent  !  go to parent, cycle
		 endif
	  endif
   end do

   if( -id /= n_leaves-1 ) then
      write(0,*)'bt_write: node miscount, ', -id, n_leaves-1
	  stop
   endif

! second traverse of bt nodes to write ID of ND, parent and children
   from_above = .true.
   nd         => bt%root

   do
      if( from_above ) then  !first time at node
	     if( associated(nd%parent) ) then
		    id_parent = nd%parent%id
		 else
		    id_parent = 0
		 endif

		 write(lu, err=120 ) nd%id, id_parent, nd%child(1)%nd%id, nd%child(2)%nd%id

		 nd%done(2) = .false.
		 if( nd%child(2)%nd%id > 0 ) nd%done(2) = .true.

		 if( nd%child(1)%nd%id < 0 ) then
		    nd         => nd%child(1)%nd
		 elseif( nd%child(2)%nd%id < 0 ) then
		    nd%done(2) = .true.
		    nd         => nd%child(2)%nd
		 else
		    nd%done = .true.
		    from_above = .false.
			if( .not.associated(nd%parent) ) exit
			nd => nd%parent
		 endif

      else !  from below

         if( .not.nd%done(2) ) then
		    from_above = .true.
			nd%done(2) = .true.
		    nd => nd%child(2)%nd
		 else
			if( .not.associated(nd%parent) ) exit
			nd => nd%parent
		 endif

      endif
   end do

   return

100   write(0,*)'bt_write: error writing n_pt'
      stop

110   write(0,*)'bt_write: error writing bt%nd_pt(i)%nd%id'
      stop

120   write(0,*)'bt_write: error writing nd%id'
      stop


end subroutine bt_write    !-----------------------------------------------

subroutine bt_read( sell, bt, idlist, lu, n_leaves )  !  read BT from logical unit LU

!  BT should not be allocated prior to call

   type (sell_type), pointer    :: sell
   type (bt_type), pointer      :: bt
   type (id_list_type), pointer :: idlist 
   integer, intent(in)          :: lu
   integer, intent(out)         :: n_leaves  ! number of leaves read

   type (nd_type), pointer      :: nd
   type (nd_pointer), pointer   :: node_pt(:)
   real(k_d) :: prop(10)
   integer   :: n_pt, check

   integer            :: i, id, n_nodes, id_node, id_parent, id_child1, id_child2
   logical            :: opened, from_above

   inquire( lu, opened=opened )
   if( .not.opened ) then
      write(0,*)'bt_write: lu not opened ', lu
	  stop
   endif

   read(lu, err=100) n_pt, check, prop

   call bt_initialize( sell, bt, n_pt=n_pt, check=check, idlist=idlist )
   bt%prop = prop

   n_leaves = 0
   do i = 1, bt%n_pt  !  read leaf NDs
      read(lu, err=110 ) id
      if( id > 0 ) then
	     n_leaves = n_leaves + 1
		 call nd_read( bt, nd, lu )
		 if( nd%id /= id ) then
		    write(0,*)'bt_read: id mismatch ', nd%id, id
			stop
		 endif
		 bt%nd_pt(i)%nd => nd
	  else
	     nullify( bt%nd_pt(i)%nd )
	  endif
   end do

   if( n_leaves < 2 ) then
      if( n_leaves == 1 ) bt%root => nd
	  return  !  all done if no nodes
   endif

   n_nodes = n_leaves - 1
   allocate( node_pt(n_nodes) )

   do id = 1, n_nodes  !  read node NDs
      call nd_read( bt, nd, lu )
	  if( nd%id /= -id ) then
	     write(0,*)'bt_read: node ID mismatch ', nd%id, -id
		 stop
	  endif
	  node_pt(id)%nd => nd
   end do

   bt%root => node_pt(1)%nd

!  traverse bt to assign child(1:2) and parent

   from_above = .true.
   nd         => bt%root

   do
      if( from_above ) then  ! entering node or leaf for the first time
	     read(lu, err=120 ) id_node, id_parent, id_child1, id_child2

		 nd => node_pt(-id_node)%nd

		 if( id_parent < 0 ) then  !  set parent
		    nd%parent => node_pt(-id_parent)%nd
		 else
		    nullify( nd%parent )
		 endif

         nd%done = .false.
		 if( id_child1 < 0 ) then  !  set child 1
		    nd%child(1)%nd => node_pt(-id_child1)%nd
		 else
            nd%child(1)%nd => bt%nd_pt(id_child1)%nd 
			bt%nd_pt(id_child1)%nd%parent => nd
			nd%done(1)     =  .true.
		 endif

		 if( id_child2 < 0 ) then  !  set child 2
		    nd%child(2)%nd => node_pt(-id_child2)%nd
		 else
            nd%child(2)%nd => bt%nd_pt(id_child2)%nd 
			bt%nd_pt(id_child2)%nd%parent => nd
			nd%done(2)     =  .true.
		 endif

		 if( .not.nd%done(1) ) then
		    nd%done(1) = .true.
		    nd         => nd%child(1)%nd
		 elseif( .not.nd%done(2) ) then
		    nd%done(2) = .true.
		    nd         => nd%child(2)%nd
		 else
		    from_above = .false.
			if( .not.associated(nd%parent) ) exit
			nd => nd%parent
		 endif

      else !  from below

         if( .not.nd%done(2) ) then
		 	nd%done(2) = .true.
		    from_above = .true.
		    nd => nd%child(2)%nd
		 else
			if( .not.associated(nd%parent) ) exit
			nd => nd%parent
		 endif

      endif
   end do

   deallocate( node_pt )
	    
   return

100   write(0,*)'bt_read:  error reading n_pt'
      stop

110   write(0,*)'bt_read:  error reading id'
      stop

120   write(0,*)'bt_read:  error reading id_node'
      stop

end subroutine bt_read    !-----------------------------------------------

subroutine nd_create( bt, node_leaf, nd )  !  allocate and initialize ND to null settings

   type (bt_type), pointer  :: bt        !  (in)  BT to which ND belongs
   integer, intent(in)      :: node_leaf ! =1 for node; =2 for leaf
   type (nd_type),  pointer :: nd        !  (out) ND

   if( .not.associated( bt ) ) then  ! check that BT exists
      write(0,*)'nd_create; BT not associated'
	  stop
   endif

   if( node_leaf < 1  .or.  node_leaf > 2 ) then
      write(0,*)'nd_create; bad value of node_leaf = ', node_leaf
	  stop
   endif

   allocate( nd )
   nd%bt => bt
   nullify(  nd%parent )

   nd%id      = -1
   nd%k       = -1
   nd%k_first = -1
   nd%done    = .false.

   if( node_leaf == 1 ) then  ! node
      allocate( nd%child(2) )
      nullify(  nd%child(1)%nd )
      nullify(  nd%child(2)%nd )
      allocate( nd%p(bt%sell%nx) )
      allocate( nd%v(bt%sell%nx) )
	  nd%p       = -1.d0
      nd%v       = -1.d0
   else                        ! leaf
      nullify( nd%child )               
      nullify( nd%p )               
      nullify( nd%v )               
   endif

   return
end subroutine nd_create  !------------------------------------------------

subroutine nd_destroy( nd )  !   deallocate and nullify ND

   type (nd_type), pointer :: nd

   if( .not.associated( nd ) ) then  ! check that ND exists
      write(0,*)'nd_destroy; nd not associated'
	  stop
   endif

   if( nd%id < 1 ) then  !  node
      deallocate( nd%child )
      deallocate( nd%p )
      deallocate( nd%v )
   endif

   deallocate( nd )

   return
end subroutine nd_destroy  !------------------------------------------------

subroutine nd_write( nd, lu )  !  write ND to logical unit LU

   type (nd_type),  pointer :: nd 
   integer, intent(in)      :: lu

   integer :: nx

   nx    = nd%bt%sell%nx

! Note: nd%k_first, nd%done(1:2) are not written (a) because they are temporary, and not needed
!       for a perfect restart (b) to ensure repeatability
 
   write(lu, err=100 ) nd%id
   write(lu, err=110 ) nd%k

   if( nd%id < 0 ) write(lu, err=120) nd%p(1:nx), nd%v(1:nx)

   return

100  write(0,*)'nd_write: error writing id'
     stop

110  write(0,*)'nd_write: error writing k'
     stop

120  write(0,*)'nd_write: error writing p, v'
     stop

end subroutine nd_write  !------------------------------------------------

subroutine nd_read( bt, nd, lu )  !  create and read ND from logical unit LU

   type (bt_type), pointer  :: bt        !  (in)  BT to which ND belongs
   type (nd_type),  pointer :: nd        !  (out) ND
   integer, intent(in)      :: lu

   integer :: id, node_leaf, nx

   nx    = bt%sell%nx

   read(lu, err=100 ) id

   node_leaf = 1
   if( id > 0 ) node_leaf = 2

   call nd_create( bt, node_leaf, nd )

   nd%id = id
   read(lu, err = 110 ) nd%k

   if( nd%id < 0 ) read(lu, err = 120 ) nd%p(1:nx), nd%v(1:nx)

   return

100  write(0,*)'nd_read: error reading id'
     stop

110  write(0,*)'nd_read: error reading k'
     stop

120  write(0,*)'nd_read: error reading p, v'
     stop

end subroutine nd_read  !------------------------------------------------

subroutine bt_nd_pt_double( bt ) 

!  double size of array bt%nd_pt

      type (bt_type),   pointer :: bt

      type (nd_pointer), pointer :: pt_temp(:)
	  integer                    :: id

!  allocate temporary array pt_temp of twice size
   	  allocate( pt_temp(2*bt%n_pt) )  

	  do id = 1, 2*bt%n_pt
	     nullify( pt_temp(id)%nd )
      end do

!  copy ell%pt to pt_temp
	  do id = 1, bt%n_pt
	     if( associated( bt%nd_pt(id)%nd ) ) pt_temp(id)%nd => bt%nd_pt(id)%nd
		 nullify( bt%nd_pt(id)%nd )
	  end do

!  double size of bt%nd_pt
	  deallocate( bt%nd_pt )

      bt%n_pt = 2 * bt%n_pt
	  allocate( bt%nd_pt(bt%n_pt) )  

	  do id = 1, bt%n_pt
	     
	     if( associated( pt_temp(id)%nd ) ) then
		    bt%nd_pt(id)%nd => pt_temp(id)%nd
	     else
		    nullify( bt%nd_pt(id)%nd )
		 endif

		 nullify(pt_temp(id)%nd) 
	  end do

	  deallocate( pt_temp )

	  return

end subroutine bt_nd_pt_double   !-----------------------------------------------

end module sell_bt
