!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module sell_m  

!  Data structures and procedures for a ellipsoids (ELLs) and sets of 
!  ellipsoids (SELLs).

!  S.B. Pope  1/15/2006, 5/15/2006

!  The dimensionality of the space is nx >=1.

!  The ellipsoid ELL is defined by  {x | |G^T * (x-c)| <= 1}, 
!  where G is lower triangular (with ng=(nx*(nx+1))/2 elements).  
!  ELL's inscribed and circumscribed balls have radii squared 
!  rsq_in and rsq_out, respectively.
!  ELL's geometry is stored as
!  ell%geom(1:ngeom) =  [c G r_in r_out] where ngeom=nx+ng+2, and
!  G is in packed format.

!  In the context of ISAT, an ellipsoid pertains to a leaf.
!  The ellipsoid ID, ell%ell_id, has the same value as the ID of the 
!  corresponding ISAT leaf.

!  An ELL may belong to only one SELL.

!  Both ELLs and SELLs are independent of the search methods used in ISAT.

!  Data structures:
!    sell_type      for SELL
!    ell_type       for ELL
!    ell_pointer

!  Subroutines:
!    sell_create( sell, nx, check, idlist, n_pt )  ! create SELL
!    sell_destroy( sell )               ! destroy SELL 
!    sell_cov( sell, mean, cov )        !  evaluate mean and covariance of ellipsoid centers 
!    sell_stretch( sell, factor )       !  uniformly stretch all ELLs in SELL by factor
!    sell_set_check( sell, check )      ! set checking level 
!    sell_write( sell, lu )             ! write SELL (for checkpointing)
!    sell_read( sell, lu, idlist )      ! read SELL

!    sell_ell_create( sell, id )        ! create ELL
!    sell_ell_destroy( sell, id )       ! destroy ELL
!    sell_ell_geom_set( ell, c, g )     ! set geometry of ELL
!    sell_ell_rad_set( ell )            ! set radii (squared) of ELL

!    sell_ell_pt_in( ell, x, in )       ! determine whether ELL covers the point x
!    sell_ell_write( ell, lu )          ! write ELL (for checkpointing)
!    sell_ell_read( ell, lu )           ! read ELL

!    sell_ell_pt_double( sell )         !  for internal use only


!  Note: for sell%check=2, checking of ID is performed, assuming that id is 
!  set through id_get prior to calling sell_ell_create, and returned by id_return 
!  after calling sell_ell_destroy.

   use id_list
   implicit none

   private sell_ell_pt_double

   integer, parameter  :: k_d = kind(1.d0)  !  kind for double precision
   integer, parameter  :: n_pt_0 = 100      !  initial size of pointer array
   integer, parameter, private :: check_def =1      !  default level of checking

   type  :: sell_type  !  SELL_TYPE ------------------------------------
      integer :: nx       !  dimensionality of space
      integer :: n_ell    !  number of ELLs in this SELL
      integer :: ng       !  storage for Cholesky triangle = ((nx*(nx+1))/2
      integer :: ngeom    !  storage for a ELLs geometry = nx + ng + 2
      integer :: n_pt     !  dimension of array of pointers to ELLs, ell_pt
      integer :: check    !  level of checking: =0 - minimal, =1 - some, =2, maximal

	  type (ell_pointer),  pointer :: ell_pt(:) ! array of pointers to ELLs: 
!                                                 sell%ell_pt(ell_id)%ell => ell
      type (id_list_type), pointer :: idlist  ! data structure used in id_list to generate IDs
!                            This is optional, but necessary for checking (level 2).

! (The following counters are treated as doubles to prevent integer overflow)
      real(k_d) :: n_ell_add       !  number of ELLs added
      real(k_d) :: n_ell_remove    !  number of ELLs removed
	  real(k_d) :: ell_id_max      !  maximum value of ell_id encountered
   end type sell_type

   type  :: ell_type  ! ELL_TYPE ------------------------------------
      type (sell_type), pointer :: sell     !  SELL to which ELL belongs
      integer                   :: ell_id   !  ID of this ELL
	  real(k_d), pointer        :: geom(:)  !  ELL data {c, G, r_in, r_out}
	  real(k_d)                 :: used     !  number of times ELL has been "used"
	  logical                   :: tested   !  used as indicator in MRU/MFU 
   end type ell_type

   type :: ell_pointer   ! ELL_POINTER -------------------------------------
      type (ell_type), pointer :: ell
   end type ell_pointer  

contains  !====================================================================

subroutine sell_create( sell, nx, check, idlist, n_pt )

!  Create a SELL (initially containing no ELLs)

   implicit none

! Required input
   type (sell_type),  pointer :: sell  ! (output)  SELL created
   integer, intent(in)        :: nx    ! dimensionality of space

! Optional input
   type (id_list_type), pointer, optional :: idlist  ! data structure used in id_list to generate IDs
   integer, intent(in), optional          :: n_pt    ! initial size of pointer array
   integer, intent(in), optional          :: check   ! checking level [0,1,2]


   integer :: i

   if( nx < 1 ) then   !  check input
      write(0,*)'sell_create: invalid nx = ', nx
	  stop
   endif

   allocate( sell )
   sell%nx    = nx
   sell%n_ell = 0
   sell%ng    = (nx*(nx+1))/2
   sell%ngeom = nx + sell%ng + 2

   if( present(check) ) then
      call sell_set_check( sell, check )
   else
      sell%check = check_def
   endif

   if( present(idlist) ) then
      sell%idlist => idlist
   else
      nullify(sell%idlist)
   endif

   if( present( n_pt ) ) then
      if( n_pt <=0 ) then
	     write(0,*)'sell_create: invalid n_pt = ', n_pt
		 stop
	  endif
	  sell%n_pt = n_pt    !  use specified value of n_pt
   else
      sell%n_pt = n_pt_0  !  use default value of n_pt
   endif

   allocate( sell%ell_pt( sell%n_pt ) )

   do i = 1, sell%n_pt
      nullify(  sell%ell_pt(i)%ell )
   end do

   sell%n_ell_add    = 0.d0       
   sell%n_ell_remove = 0.d0   
   sell%ell_id_max   = 0.d0  
   
   return   

end subroutine sell_create  !-----------------------------------------------

subroutine sell_destroy( sell )

!  Destroy and deallocate SELL and all of its ELLs 
!  All ELLs must be removed from all datastructures prior to call

   implicit none
   type (sell_type), pointer :: sell  ! (in)  SELL to be destroyed

   integer :: id

   if( .not.associated( sell ) ) then
      write(0,*)'sell_destroy: SELL not associated'
	  stop
   endif

   do id = 1, sell%n_pt
      if( associated( sell%ell_pt(id)%ell ) ) call sell_ell_destroy( sell, id )
   end do

   deallocate( sell%ell_pt )
   deallocate( sell )
   
   return   

end subroutine sell_destroy  !-----------------------------------------------

subroutine sell_cov( sell, mean, cov )  !  evaluate mean and covariance of ellipsoid centers

   implicit none
   type (sell_type), pointer :: sell            ! SELL  - set of ellipsoids  
   real(k_d), intent(out)    :: mean(sell%nx)   !  mean of  ELL centers 
   real(k_d), intent(out)    :: cov(sell%nx,sell%nx)  !  covariance 

   integer   :: nx, id, i, j, n
   real(k_d) :: c(sell%nx)

   if( .not.associated(sell)  ) then 
      write(0,*)'sell_cov: sell not associated'
	  stop
   endif

   mean = 0.d0  !  mean of centers (to be formed)
   cov  = 0.d0  !  covariance matrix (to be formed)

   if( sell%n_ell < 1 ) return  !  no ELLs

   nx    = sell%nx

   n     = 0     !  count number of ELLs
   do id = 1, sell%n_pt
	  if( .not.associated( sell%ell_pt(id)%ell ) ) cycle
	  n    = n + 1
	  mean = mean + sell%ell_pt(id)%ell%geom(1:nx)
   end do

   if( n /= sell%n_ell ) then
      write(0,*)'sell_cov: bad ELL count ', n, sell%n_ell 
	  stop
   endif

   mean = mean / n

!  form upper triangle of cov(i,j)

   do id = 1, sell%n_pt
	  if( .not.associated( sell%ell_pt(id)%ell ) ) cycle
	  c = sell%ell_pt(id)%ell%geom(1:nx) - mean
      do j = 1, nx
	     cov(1:j,j) = cov(1:j,j) + c(1:j) * c(j)
	  end do
   end do

   do i = 2, nx  !  set lower triangle
      cov(i,1:i-1) = cov(1:i-1,i)
   end do

   cov  = cov / n

   return
end subroutine sell_cov  !-------------------------------------------------------

subroutine sell_stretch( sell, factor )  !  uniformly stretch all ELLs in SELL by factor

   implicit none
   type (sell_type), pointer :: sell     ! SELL  - set of ellipsoids  
   real(k_d), intent(in)     :: factor   ! factor by which ELLs are to be stretched

   type (ell_type), pointer  :: ell
   integer   :: ngs, ngf, nri, nro, id, n
   real(k_d) :: fac_inv

   if( .not.associated(sell)  ) then 
      write(0,*)'sell_stretch: sell not associated'
	  stop
   endif

   if( sell%n_ell < 1 ) return  !  no ELLs

   if( factor <= 0.d0 ) then
      write(0,*)'sell_stretch: factor <= 0 ', factor
	  stop
   endif

   fac_inv = 1.d0 / factor

   ngs = sell%nx + 1
   ngf = sell%nx + sell%ng
   nri = sell%ngeom - 1
   nro = sell%ngeom

   n     = 0     !  count number of ELLs
   do id = 1, sell%n_pt
	  if( .not.associated( sell%ell_pt(id)%ell ) ) cycle
	  n    = n + 1
	  ell => sell%ell_pt(id)%ell
	  ell%geom(ngs:ngf) = ell%geom(ngs:ngf) * fac_inv
	  ell%geom(nri:nro) = ell%geom(nri:nro) * factor
   end do

   if( n /= sell%n_ell ) then
      write(0,*)'sell_stretch: bad ELL count ', n, sell%n_ell 
	  stop
   endif

   return
end subroutine sell_stretch  !-------------------------------------------------------

subroutine sell_set_check( sell, check )

!  Set checking level

   implicit none
   type (sell_type), pointer :: sell  ! (in)  SELL 
   integer, intent(in)       :: check ! checking level [0,1,2]

   if( .not.associated( sell ) ) then
      write(0,*)'sell_set_check: sell not associated'
	  stop
   endif

   if( check < 0  .or.  check > 2  ) then
      write(0,*)'sell_set_check: stopping,  check out of range [0,2], check = ', check
	  stop
   endif

   sell%check = check
   
   return   

end subroutine sell_set_check  !-----------------------------------------------

subroutine sell_write( sell, lu )

!  Write SELL to checkpoint file on unit LU

   implicit none

   type (sell_type),  pointer :: sell  ! SELL to be written
   integer, intent(in)        :: lu    ! logical unit number (of open unformatted file)

   integer              :: id
   integer, allocatable :: ell_id(:)
   logical              :: opened

   if( .not.associated(sell) ) then  !  check input
     write(0,*)'sell_write: SELL not associated'
	 stop
   endif

   if( lu < 1 ) then
      write(0,*)'sell_write:  lu < 1, ', lu
	  stop
   endif

!  check file status
   inquire( lu, err=100, opened=opened )

   if( .not.opened ) then
      write(0,*)'sell_write:  file not opened', lu
	  stop
   endif

   write( lu, err=110 ) sell%nx, sell%n_ell
   write( lu, err=120 ) sell%ng, sell%ngeom, sell%n_pt, sell%check 

   if( associated( sell%idlist ) ) then
      write(lu) 1
   else
      write(lu) 0
   endif

   write( lu, err=130 ) sell%n_ell_add, sell%n_ell_remove, sell%ell_id_max

   allocate( ell_id(sell%n_pt) )
   ell_id = 0
   do id = 1, sell%n_pt
      if( associated( sell%ell_pt(id)%ell ) ) ell_id(id) = sell%ell_pt(id)%ell%ell_id
   end do
   write( lu, err=140 ) ell_id

   deallocate( ell_id )

   do id = 1,  sell%n_pt
      if( associated(sell%ell_pt(id)%ell) ) call sell_ell_write( sell%ell_pt(id)%ell, lu )
   end do

   return   

100   continue
      write(0,*)'sell_write:  error making file inquiry'
	  stop

110   continue
      write(0,*)'sell_write:  writing nx...'
	  stop

120   continue
      write(0,*)'sell_write:  writing ng...'
	  stop

130   continue
      write(0,*)'sell_write:  writing n_ell_add...'
	  stop
	  
140   continue
      write(0,*)'sell_write:  writing ell_id'
	  stop


end subroutine sell_write  !-----------------------------------------------

subroutine sell_read( sell, lu, idlist )

!  Allocate SELL and read properties from checkpoint file on unit LU.

   implicit none

   type (sell_type),  pointer   :: sell  ! SELL to be created and read
   integer, intent(in)          :: lu    ! logical unit number (of open unformatted file)
   type (id_list_type), pointer :: idlist

   integer              :: i, id, nx, n_ell, ng, ngeom, n_pt, check
   integer, allocatable :: ell_id(:)
   logical              :: opened
   real(k_d)            :: n_ell_add   

   if( lu < 1 ) then
      write(0,*)'sell_read:  lu < 1, ', lu
	  stop
   endif

!  check file status
   inquire( lu, err=100, opened=opened )

   if( .not.opened ) then
      write(lu,*)'sell_read:  file not opened', lu
	  stop
   endif

   if( .not.associated(idlist) ) then
      write(lu,*)'sell_read:  idlist not associated'
	  stop
   endif

   read( lu, err=110 ) nx, n_ell
   read( lu, err=120 ) ng, ngeom, n_pt, check 

   call sell_create( sell, nx, check=check, idlist=idlist, n_pt=n_pt )

   if( sell%ng /= ng  .or.  sell%ngeom /= ngeom ) then
      write(lu,*)'sell_read:  incompatible ng and/or ngeom ', sell%ng, ng, sell%ngeom, ngeom
	  stop
   endif

   read(lu) i  ! indicator for idlist being associated
   if( i ==0 ) nullify( sell%idlist )

   read( lu, err=130 ) n_ell_add, sell%n_ell_remove, sell%ell_id_max

   allocate( ell_id(sell%n_pt) )
   read( lu, err=140 ) ell_id

   do id = 1, sell%n_pt
      if( ell_id(id) > 0 ) then
	     call sell_ell_create( sell, id )
		 call sell_ell_read( sell%ell_pt(id)%ell, lu )
      else
	     nullify( sell%ell_pt(id)%ell )
	  endif	  
   end do

   if( sell%n_ell /= n_ell ) then
      write(lu,*)'sell_read:  n_ell mismatch ', sell%n_ell, n_ell
	  stop
   endif

   sell%n_ell_add = n_ell_add

   deallocate( ell_id )

   return   

100   continue
      write(0,*)'sell_read:  error making file inquiry'
	  stop

110   continue
      write(0,*)'sell_read:  reading nx...'
	  stop

120   continue
      write(0,*)'sell_read:  reading ng...'
	  stop

130   continue
      write(0,*)'sell_read:  reading n_ell_add...'
	  stop
	  
140   continue
      write(0,*)'sell_read:  reading ell_id'
	  stop

end subroutine sell_read  !-----------------------------------------------

subroutine sell_ell_create( sell, id )

!  Create and initialize ELL with ID=id (with null geometry) and add to SELL.

   implicit none
   type (sell_type), pointer :: sell  !  SELL to which ELL will belong
   integer, intent(in   )    :: id    !  ID 

   integer :: id_status, ids_assigned, id_max 

   if( .not.associated( sell ) ) then
      write(0,*)'sell_ell_create: SELL not associated'
	  stop
   endif

   if( sell%check > 0 ) then !  check ID
      if( id < 1 ) then
         write(0,*)'sell_ell_create: id < 1; id = ', id
         stop
      endif

	  if( sell%check > 1  .and.associated(sell%idlist) ) then
	     call id_query( sell%idlist, id, id_status, ids_assigned, id_max )
		 if( id_status /= 1 ) then
		    write(0,*)'sell_ell_create: stopping, id not assigned in idlist', id
            stop
         endif
	  endif
   endif

   do while ( id > sell%n_pt )
      call sell_ell_pt_double( sell ) 
   end do

   allocate( sell%ell_pt( id )%ell )
   allocate( sell%ell_pt( id )%ell%geom(sell%ngeom) )

   sell%ell_pt( id )%ell%sell   => sell
   sell%ell_pt( id )%ell%ell_id =  id
   sell%ell_pt( id )%ell%geom   =  -1.d0  !  indicate geometry not set
   sell%ell_pt( id )%ell%used   =  0.d0

   sell%n_ell      = sell%n_ell     + 1
   sell%n_ell_add  = sell%n_ell_add + 1.d0
   sell%ell_id_max = max( sell%ell_id_max, 1.d0*id )

   return

end subroutine sell_ell_create  !-----------------------------------------------

subroutine sell_ell_destroy( sell, id )

!  Remove ELL from SELL; destroy and deallocate ELL 
!  ELL must be removed from all datastructures prior to call

   implicit none
   type (sell_type), pointer :: sell  !  (in) ELL to be destroyed
   integer, intent(in)       :: id    !  ID of ELL to destroy

   integer :: id_status, ids_assigned, id_max

   if( .not.associated( sell ) ) then
      write(0,*)'sell_ell_destroy: SELL not associated'
	  stop
   endif

   if( sell%check > 0 ) then  !  check ID
      if( id < 1 ) then
         write(0,*)'sell_ell_destroy: id < 1; id = ', id
         stop
      endif

	  if( sell%check > 1  .and.associated(sell%idlist) ) then
	     call id_query( sell%idlist, id, id_status, ids_assigned, id_max )
		 if( id_status /= 1 ) then
		    write(0,*)'sell_ell_destroy: stopping, ID not assigned in idlist', id
            stop
         endif
	  endif
   endif

   if( .not.associated(sell%ell_pt(id)%ell) ) then
      write(0,*)'sell_ell_destroy: ell not associated'
      stop
   endif

   deallocate( sell%ell_pt(id)%ell%geom ) 
   deallocate( sell%ell_pt(id)%ell )

   sell%n_ell        = sell%n_ell        -1
   sell%n_ell_remove = sell%n_ell_remove +1.d0

   return

end subroutine sell_ell_destroy  !-----------------------------------------------

subroutine sell_ell_geom_set( ell, c, g )

!  Set (or reset) the geometrical properties of an ELL.
!  ELL must have been created prior to call.

   implicit none
   type (ell_type), pointer :: ell    !  (in/out) ELL which must be created prior to call
   real (k_d), intent(in)   :: c(ell%sell%nx), g(ell%sell%ng)  ! {c, G} for ELL

   integer    :: nx, ng

   if( .not.associated(ell) ) then
      write(0,*)'sell_ell_geom_set: ELL not associated'
      stop
   endif

   nx = ell%sell%nx
   ng = ell%sell%ng

   ell%geom(1:nx)       = c(1:nx)
   ell%geom(nx+1:nx+ng) = g(1:ng)

   call sell_ell_rad_set( ell )

   return
end subroutine sell_ell_geom_set   !-----------------------------------------------

subroutine sell_ell_rad_set( ell )

!  Set (or reset) the lower and upper radii of an ELL.

   implicit none
   type (ell_type), pointer :: ell    !  (in/out) ELL which must be created prior to call

   integer    :: nx, ng
   real (k_d) :: r_in, r_out

   if( .not.associated(ell) ) then
      write(0,*)'sell_ell_rad_set: ELL not associated'
      stop
   endif

   nx = ell%sell%nx
   ng = ell%sell%ng

   call ell_radii( nx, ell%geom(nx+1:nx+ng), r_in, r_out )
   ell%geom(nx+ng+1) = r_in
   ell%geom(nx+ng+2) = r_out

   return
end subroutine sell_ell_rad_set   !-----------------------------------------------

subroutine sell_ell_pt_in( nx, ngeom, geom, x, in )

!  Determine whether ELL covers the point x.

   implicit none
   integer, intent(in)    :: nx           !  dimensionality of space
   integer, intent(in)    :: ngeom        !  dimension of geom = nx + (nx*(nx+1))/2 + 2
   real (k_d), intent(in) :: geom(ngeom)  !  [c G r_in r_out]
   real (k_d), intent(in) :: x(nx)        !  query point
   logical, intent(out)   :: in           !  .true. if ELL covers x
   
   integer    :: i, is, ie
   real (k_d) :: r(nx), rsq, rsq_ell
    
   in      = .false.  !  assume x not covered by ELL until found otherwise
   
!-----------  form partial sum of rsq and test
   rsq     = 0.d0
   rsq_ell = geom(ngeom)**2  ! r_out^2

   do i = 1, nx
      r(i) = x(i) - geom(i)         ! r = x - c
	  rsq  = rsq + r(i)**2
	  if( rsq > rsq_ell ) return  !  x is outside circumscribed spheroid
   end do

   if( rsq <= geom(ngeom-1)**2 ) then  !  x is inside inscribed ball spheroid
      in = .true.
	  return
   endif

!  r = G^T * [x-c];  in  if  |r| < 1    
!------ partial sum, backward
   rsq = 0.d0
   ie  = ngeom - 2
      
   do i = nx, 1, -1
      is  = ie - nx + i
      rsq = rsq + sum( geom(is:ie)*r(i:nx) )**2
      if( rsq > 1.d0 ) return
      ie = is-1
   end do
   
   in = .true.
             
   return
end subroutine sell_ell_pt_in   !-----------------------------------------------

subroutine sell_ell_pt_in_test( nx, ngeom, geom, x, in )

!  Determine whether ELL covers the point x.
!  This version to test efficiency of different codings.

   implicit none
   integer, intent(in)    :: nx           !  dimensionality of space
   integer, intent(in)    :: ngeom        !  dimension of geom = nx + (nx*(nx+1))/2 + 2
   real (k_d), intent(in) :: geom(ngeom)  !  [c G r_in r_out]
   real (k_d), intent(in) :: x(nx)        !  query point
   logical, intent(out)   :: in           !  .true. if ELL covers x
   
   integer, save :: mode_r=0, mode_g

   integer    :: i, is, ie
   real (k_d) :: r(nx), rsq, rsq_ell
   
   if( mode_r == 0 ) then
      read(1,*) mode_r, mode_g
      write(0,*) 'mode_r, mode_g = ', mode_r, mode_g
   endif
      
   in      = .false.  !--------  assume x not covered by ELL until found otherwise
   
   if( mode_r == 1 ) then  !  no test on radii
      r = x - geom(1:nx)
      
   elseif( mode_r == 2 ) then  !---------  form rsq then test
      r = x - geom(1:nx)
      rsq = sum(r*r)
      if( rsq > geom(ngeom)**2 ) return
      if( rsq <= geom(ngeom-1)**2 ) then
         in = .true.
         return
      endif
      
   else  !-----------  form partial sum of rsq and test
      rsq     = 0.d0
      rsq_ell = geom(ngeom)**2  ! r_out^2

      do i = 1, nx
         r(i) = x(i) - geom(i)         ! r = x - c
	     rsq  = rsq + r(i)**2
	     if( rsq > rsq_ell ) return  !  x is outside circumscribed spheroid
      end do

      if( rsq <= geom(ngeom-1)**2 ) then    !  x is inside inscribed ball spheroid
         in = .true.
	     return
      endif
   endif

 !  r = G^T * [x-c];  in  if  |r| < 1
 
   if( mode_g == 1 ) then  !-----  form r directly
      call dtpmv( 'L', 'T', 'N', nx, geom(nx+1:ngeom-2), r, 1 ) 
      rsq = sum(r*r) 
      if( rsq > 1.d0 ) return
      
   elseif( mode_g == 2 ) then  !------ partial sum, forward
      rsq = 0.d0
      is  = nx + 1
      
      do i = 1, nx
         ie  = is + nx - i
         rsq = rsq + sum( geom(is:ie)*r(i:nx) )**2
         if( rsq > 1.d0 ) return
         is = ie+1
      end do
      
    else   !------ partial sum, backward
      rsq = 0.d0
      ie  = ngeom - 2
      
      do i = nx, 1, -1
         is  = ie - nx + i
         rsq = rsq + sum( geom(is:ie)*r(i:nx) )**2
         if( rsq > 1.d0 ) return
         ie = is-1
      end do
   endif
   
   in = .true.
             
   return
end subroutine sell_ell_pt_in_test   !-----------------------------------------------

subroutine sell_ell_pt_in_orig( nx, ngeom, geom, x, in )

!  Determine whether ELL covers the point x.

   implicit none
   integer, intent(in)    :: nx           !  dimensionality of space
   integer, intent(in)    :: ngeom        !  dimension of geom = nx + (nx*(nx+1))/2 + 2
   real (k_d), intent(in) :: geom(ngeom)  !  [c G r_in r_out]
   real (k_d), intent(in) :: x(nx)        !  query point
   logical, intent(out)   :: in           !  .true. if ELL covers x

   integer    :: i
   real (k_d) :: r(nx), rsq, rsq_ell

   in      = .false.  !  assume x not covered by ELL until found otherwise
   rsq     = 0.d0
   rsq_ell = geom(ngeom)**2  ! r_out^2

   do i = 1, nx
      r(i) = x(i) - geom(i)         ! r = x - c
	  rsq  = rsq + r(i)**2
	  if( rsq > rsq_ell ) return  !  x is outside circumscribed spheroid
   end do

   if( rsq <= geom(ngeom-1)**2 ) then    !  x is inside inscribed ball spheroid
      in = .true.
	  return
   endif

 !  r = G^T * [x-c]
   call dtpmv( 'L', 'T', 'N', nx, geom(nx+1:ngeom-2), r, 1 ) 
   rsq = sum(r*r) 
   if( rsq <= 1.d0 ) in = .true.

   return
end subroutine sell_ell_pt_in_orig   !-----------------------------------------------

subroutine sell_ell_write( ell, lu )

!  Write ELL to checkpoint file on unit LU

   implicit none

   type (ell_type),  pointer :: ell  ! ELL to be written
   integer, intent(in)       :: lu   ! logical unit number (of open unformatted file)

   logical :: opened

   if( .not.associated(ell) ) then  !  check input
     write(0,*)'sell_ell_write: SELL not associated'
	 stop
   endif

   if( lu < 1 ) then
      write(lu,*)'sell_ell_write:  lu < 1, ', lu
	  stop
   endif

!  check file status
   inquire( lu, err=100, opened=opened )

   if( .not.opened ) then
      write(lu,*)'sell_ell_write:  file not opened', lu
	  stop
   endif

   write( lu, err=110 ) ell%ell_id
   write( lu, err=120 ) ell%geom(1:ell%sell%ngeom)
   write( lu, err=130 ) ell%used
   
   return   

100   continue
      write(lu,*)'sell_ell_write:  error making file inquiry'
	  stop

110   continue
      write(lu,*)'sell_ell_write:  writing ell_id'
	  stop

120   continue
      write(lu,*)'sell_ell_write:  writing geom'
	  stop

130   continue
      write(lu,*)'sell_ell_write:  writing used'
	  stop

end subroutine sell_ell_write  !-----------------------------------------------

subroutine sell_ell_read( ell, lu )

!  Read ELL from checkpoint file on unit LU

   implicit none

   type (ell_type),  pointer :: ell  ! ELL to be written
   integer, intent(in)       :: lu   ! logical unit number (of open unformatted file)

   logical :: opened

   if( .not.associated(ell) ) then  !  check input
     write(0,*)'sell_ell_read: ELL not associated'
	 stop
   endif

   if( lu < 1 ) then
      write(lu,*)'sell_ell_read:  lu < 1, ', lu
	  stop
   endif

!  check file status
   inquire( lu, err=100, opened=opened )

   if( .not.opened ) then
      write(lu,*)'sell_ell_read:  file not opened', lu
	  stop
   endif

   read( lu, err=110 ) ell%ell_id
   read( lu, err=120 ) ell%geom(1:ell%sell%ngeom)
   read( lu, err=130 ) ell%used
   
   return   

100   continue
      write(lu,*)'sell_ell_read:  error making file inquiry'
	  stop

110   continue
      write(lu,*)'sell_ell_read:  reading ell_id'
	  stop

120   continue
      write(lu,*)'sell_ell_read:  reading geom'
	  stop

130   continue
      write(lu,*)'sell_ell_read:  reading used'
	  stop

end subroutine sell_ell_read  !-----------------------------------------------

subroutine sell_ell_pt_double( sell ) 

!  double size of array sell%ell_pt

      type (sell_type),   pointer :: sell

      type (ell_pointer), pointer :: pt_temp(:)
	  integer                     :: id

!  allocate temporary array pt_temp of twice size
   	  allocate( pt_temp(2*sell%n_pt) )  

	  do id = 1, 2*sell%n_pt
	     nullify( pt_temp(id)%ell )
      end do

!  copy ell%pt to pt_temp
	  do id = 1, sell%n_pt
	     if( associated( sell%ell_pt(id)%ell ) ) pt_temp(id)%ell => sell%ell_pt(id)%ell
		 nullify( sell%ell_pt(id)%ell )
	  end do

!  double size of sell%ell_pt
	  deallocate( sell%ell_pt )

      sell%n_pt = 2 * sell%n_pt
	  allocate( sell%ell_pt(sell%n_pt) )  

	  do id = 1, sell%n_pt
	     
	     if( associated( pt_temp(id)%ell ) ) then
		    sell%ell_pt(id)%ell => pt_temp(id)%ell
	     else
		    nullify( sell%ell_pt(id)%ell )
		 endif

		 nullify(pt_temp(id)%ell) 
	  end do

	  deallocate( pt_temp )

	  return

end subroutine sell_ell_pt_double   !-----------------------------------------------

end module sell_m 