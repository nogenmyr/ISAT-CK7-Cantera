!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!


subroutine ci_ice_pic_bound( r0, drds, sb, kfa, indic_rpos, indic_zpos, r_normal, iflag )
     
!  Determine the intersection of the line  r(s) = r0 + s * drds  with the boundary
!  of the reduced realizable region, with r_normal^T * drds >0, where r_normal is
!  the outward normal.    Linear programming is used to determine
!  sb, the maximum of s, subject to z = B^T( r + s * drds ) >= 0. 

!  r0 and drds are assumed to satisfy the normalization conditions:
!     w^T * r0 = 1, w^T * drds = 0,  where w = amolwt_n.

!  Input:
!    r0   - reduced composition {z^r, z^ue}
!    drds - tangent along line

!  Output:
!    sb   - distance from r0 to the intersection
!    kfa  - information about the facet intersected
!           kfa=k if r(k)=0 at the intersection (for some k, 1<= k <= nrc)
!           kfa=nrc+1  if all components of r are strictly positive at the intersection
!    indic_rpos - indicator of reduced compositions that can be positive on the facet
!    indic_zpos - indicator of species that can be positive on the facet
!    r_normal   - unit outward normal to the facet in r-space
!    iflag > 0 for success:
!            = 1 a component of r is zero on the boundary
!            = 2 a determined, unrepresented species is zero on the boundary
!            = 3 an undetermined, unrepresented species is zero on the boundary
!          < 0 for failure: 
!            = {-1,-2} no intersection
!            = -3 LP failure 
!            = {-10:-4} normal failure

!  This file also contains subroutine ci_ice_pic_bound_test, which can be called to 
!  test ci_ice_pic_bound.

  use ci_dat
  use ci_dat8
  use linprog
  use ci_utils
  use ci_stats

  implicit none
  
  logical, parameter    :: quit = .false.    ! set true to abort on failure; 
                                             ! otherwise return with  iflag < 0       
  logical, parameter    :: if_svd = .false.  ! set true to use SVD to check normal
  real(k_dp), parameter :: drds_lim = 0.d0   ! limit on small components

  real(k_dp), intent(in)  :: r0(nrc), drds(nrc)
  real(k_dp), intent(out) :: sb, r_normal(nrc)
  integer,    intent(out) :: kfa, indic_rpos(nrc), indic_zpos(ns), iflag
       
! local arrays
         
  integer       :: i, j, nx, nle, nge, neq, lhv(ne)
  
  real(k_dp), allocatable, save  :: aeq(:,:), beq(:),  f(:), x(:)
  real(k_dp), save   :: calls = 0.d0
  real(k_dp)    :: ale(1,1), age(1,1), ble(1), bge(1),  &
                   stats_lp(30), ss, tol_use, smin, smax,  &
                   sp, sint, drdsm(nrc), drds_min, &
                   adet, bdet, rr
   
   tol_use  = 1.d-12  ! tolerance used for LP and other purposes 
    
   call routine_start(i_ci_ice_pic_bound)

!  determine the intersections (s=smin, s=smax) of the line with the simplex r >=0
!  (which covers the realizable region)
 
  if( calls == 0.d0 ) then
     call determined  ! identify determined unrepresented species and other information
     nx = nzud+1
     allocate( aeq(nrud,nx) )
     allocate( beq(nrud) )
     allocate( f(nx) )
     allocate( x(nx) )
  endif
  
  calls = calls + 1.d0
  
  drdsm    = drds ! modified version
  drds_min = drds_lim * maxval( abs(drds) )

  smax =  huge(0.d0)
  smin = -huge(0.d0)
  
  do i = 1, nrc
     if( abs( drdsm(i) ) <= drds_min ) then
        drdsm(i) = 0.d0
        if( r0(i) < 0.d0 ) then
           call abort( 1, 'no intersection', isv=i )
           return
        endif
        
     else
        sint = - r0(i) / drdsm(i) ! intersection
        
        if( drdsm(i) > 0.d0 ) then  
           smin = max( smin, sint )  !  sint is left intersection
        else
           if( sint < smax ) then
              smax = sint            !  sint is right intersection
              kfa  = i
           endif
        endif
     endif
  end do
  
!  Consider determined, unrepresented species.  
!  (Normally k_zdet(i) <= 1)

  do i = 1, ns
     if( k_zdet(i) <= 1 ) cycle
     ! Need only consider species determined by two or more components of r,
     ! since those determined by a single component have already been treated (above).
     ! Along the line, z(i) = adet + bdet * s
     adet = dot_product( r0,   r2zdet(i,:) )
     bdet = dot_product( drds, r2zdet(i,:) )
     
     if( bdet == 0.d0 ) then
        if( adet < 0.d0 ) then
           call abort( 2, 'no intersection', isv=i )
           return
        endif
        
     else
        sint = - adet / bdet ! intersection
        
        if( bdet > 0.d0 ) then  
           smin = max( smin, sint )  !  sint is left intersection
        else
           if( sint < smax ) then
              smax = sint            !  sint is right intersection
              kfa  = -i              !  index of intersection
           endif
        endif
     endif
     
  end do
        
  if( smin > smax ) then
     call abort( 1, 'no intersection' )
     return
  endif
  
  sp = smax + sqrt(tol_use) * (1.d0 + abs(smax) ) ! define sp such that r0 + sp * drdsm is exterior
  
!  Now treat undetermined, unrepresented species, zud.
!  Intersection is at s = sp - ds,  ds > 0.  
!  Solve LP for ds:  minimize ds, subject to rud = rud0 + (sp-ds)*d(rud)/ds = BBud^T*zud, zud >=0.      
        
  nx  = nzud + 1 ! x = [zud; ds],  s = sp - ds,  ds > 0
  nle = 0
  nge = 0 
  neq = nrud    
  ale = 0.d0
  ble = 0.d0
  age = 0.d0   
  bge = 0.d0
       
 !  specify equality constraints: unrepresented elements 
  do i = 1, nzud
     do j = 1, nrud
        aeq(j,i) = BB(i_zud(i),i_rud(j)) 
     end do
  end do
  
  do j = 1, nrud
     i = i_rud(j) ! index of j-th undetermined component of r
     aeq(j,nx) = drdsm( i ) 
     beq(j) = r0(i) + sp * drdsm(i)
  end do
  
  rr  = maxval( abs(aeq(:,nx)) ) ! scale
  aeq(:,nx) = aeq(:,nx) / rr
       
  ss  = maxval( abs( beq) )  ! scale
  beq = beq / ss
       
  f        = 0.d0
  f(nx)    = 1.d0 ! minimize ds 
  stats_lp = 0.d0

  call routine_start(i_lp)

  call lp( nx, nle, nge, neq, ale, age, aeq, ble, bge, beq, &
           f, x, iflag, toler=tol_use, stats=stats_lp, &
           lhvars=lhv, t_pivot=1.d-20, t_z_stop=tiny(0.d0)  ) 

  call routine_stop(i_lp)

  x     = x * ss      !  unscale solution
  x(nx) = x(nx) / rr
  sb    = sp - x(nx)  !  location of intersection based on zue
           
  if( (iflag <= -2  .and.  iflag >= -4) .or. &
      (iflag == -1  .and.  maxval(aeq(:,nx)) < 0.d0 )  .or. & 
       sb > smax-tol_use ) then
       
!  Intersection is sb = smax (and not determined by LP solution) if 
!    (a) there is no feasible solution with ds >=0, or
!    (b) the LP solution is unbounded (which can happen if max(drds_ud)<0)
!    (c) sb (based on zud) is greater than smax
       
     sb         = smax
     indic_rpos = 1
     indic_zpos = 1
     
     if( kfa > 0 ) then  !  intersection with r(kfa) = 0
        iflag           = 1
        indic_rpos(kfa) = 0  !  zero r
        r_normal        = 0.d0
        r_normal(kfa)   = -1.d0 ! outward normal to facet r(kfa) = 0
        
        do i = 1, ns  !   identify zero z
           if( BB(i,kfa) > 0.d0 ) indic_zpos(i) = 0
        end do

     else  !  intersection with determined species z(-kfa) = 0
        iflag    = 2
        indic_zpos(-kfa) = 0
        
        r_normal = r2zdet(-kfa,:)
        r_normal = r_normal / norm( r_normal ) ! make unit outward
        if( dot_product( r_normal, drds ) < 0.d0 ) r_normal = -r_normal
        kfa      = nrc + 1
     endif

     call routine_stop(i_ci_ice_pic_bound)
     return  !  all done

  elseif( iflag /= 0 ) then  !  other LP failure (unbounded solution...)
     call abort( 4, 'LP failure, iflag =', isv = iflag )
     return
  endif
  
  if( sb < smin ) then
     call abort( 1, 'no intersection' )
     return
  endif
  
 !  intersection is determined by zud                           
  kfa   = nrc + 1
  iflag = 3
  
  call normal ! determine: r_normal, indic_zpos, indic_rpos
 
  return
  
  contains !-----------------------------------------------
  
    subroutine determined
  
!  Determine if any of the unrepresented species is completely determined by
!  the constraint on elements of unrepresented species.  
!  Determined species include:
!     1/ represented species
!     2/ unrepresented species that contain all of the unrepresented atoms of 
!          some element (e.g., noble gases, or N2 if no other unrepresented 
!          species contain N)
!     3/ possibly other unrepresented species, determined by the SVD of CEu, 
!          the element matrix of unrepresented species.
!  This routine is called just once, during initialization.

!  The following quantities (defined in ci_dat8) are set:
!    nzud - the number of undetermined species
!    nrud - the number of "undetermined" reduced compositions, i.e., the
!              number of reduced compositions (or, equiv., elements) in the
!              undetermined species
!    i_zud(i) - index of the i-th undetermined species (i=1:nzud)
!    i_rud(i) - index of the i-th undetermined reduced composition (i=1:nrud)
!    k_zdet(i) = 0 if species i is undetermined
!              = k > 0 if species i is determined by k components of r
!    r2zdet(ns,nrc) - matrix such that zdet = r2zdet * r yields z(i)=zdet(i) for
!                     determined species, and zdet(i)=0 for undetermined species

!  Caution: the unrepresented species are not necessarily contiguous in arrays. 
!           They can be identified by S2RS(i)==0.

   use ci_utils
      
   integer, save :: lu = -1  !  logical unit number of output, or lu<0 to suppress
   
   integer    :: i, ie, is, ieu, iuds, info, rank, ius, indic_zu_det(nus), &
                 indic_e_undet(ne), nzu_det
   real(k_dp) :: Eu(nus,ne), S(ne), U(nus,nus), VT(ne,ne), zue_to_zu(nus,ne), &
                 work(nus*nus+10*(nus+ne)), sum_comp
        
   allocate( k_zdet(ns) )
   allocate( r2zdet(ns,nrc) )
        
!  First phase: CEu(nus,ne) is the element matrix of the unrepresented species. 
!  In this phase zu(1:nus) are the unrepresented species (packed in order).                 
!  CEu is such that:  zue = CEu^T * zu.
!  The SVD of CEu is:  CEu = U * S * VT = [Ut Uh] * [St;0] * VT.
!  Thus:  zu = Ut * St^{-1} * VT * zue + Uh * b, where b is arbitrary. 
!  Hence, zu(j) is determined iff Uh(j,:)=0.
   if( lu >= 0 ) then
      do i = 1, nus
         write(lu,'(i4, 4x, 10i3)') i, (nint(CEu(i,j)), j=1,ne)
      end do
   endif

!  perform SVD of CEu
   Eu = CEu
   call dgesvd( 'A', 'A', nus, ne, Eu, nus, S, U, nus, VT, ne, &
                 work, nus*nus+10*(nus+ne), info )
      
   if (info /=0) call isat_abort('Determined', 1, mess= &
          'SVD failed, info = ', isv = info )
      
   rank = 0  !  check that rank of CEu is ne
   do i = 1, ne
      if( S(i) > epsilon_rank * S(1) ) rank = rank + 1
   end do
   
   if (rank /= ne ) call isat_abort('Determined', 2, mess= &
          'CEu is rank deficient, rank = ', isv = rank )
  
!  identify any determined species
   indic_zu_det = 0  !  indicator of determined unrepresented species
   nzu_det      = 0  !     number of determined unrepresented species
   do i = 1,nus
      if( maxval( abs( U(i,ne+1:nus) ) ) < epsilon_rank ) then
         nzu_det         = nzu_det + 1  !  species is determined
         indic_zu_det(i) = 1
      endif
   end do
   
   nzud = nus - nzu_det  ! number of undetermined, unrepresented species
   
   allocate( urud_to_ur(nzud) )
   
   ! determine urud_to_ur: 
   ! urud_to_ur(i) is the unrepresented-species index of the i-th undetermined, unrepresented species
   i = 0
   do ius = 1, nus
      if( indic_zu_det(ius) == 0 ) then
         i = i + 1
         urud_to_ur(i) = ius
      endif
   end do
   
   if( lu >= 0 ) then
      write(lu,*)' '
      write(lu,'(a,100i3)')' urud_to_ur = ', urud_to_ur
      write(lu,*)' '
   endif
   
   !  form zue_to_zu(1:nus,1:ne)  such that, for determined species: zu = zue_to_zu * zue
   S = 1.d0 / S
   zue_to_zu = matmul( U(:,1:ne), diagmul( S, VT ) )
   
   if( lu >= 0 ) then
      do i = 1, nus
         write(lu,'(i4,1p,10e13.4)') indic_zu_det(i), (zue_to_zu(i,j), j=1,ne)
      end do
   endif
   
   !write(lu,*)'nzu_det = ', nzu_det
   
   !  identify elements in undetermined species
   indic_e_undet = 0
   do ie = 1, ne
   
     sum_comp = 0.d0 ! sum of |zue_to_zu(:,ie)| for undetermined species
     do ius = 1, nus
        if( indic_zu_det(ius) == 0 ) sum_comp = sum_comp + abs( zue_to_zu(ius,ie) )
     end do
   
     if( sum_comp > epsilon_rank ) indic_e_undet(ie) = 1 ! element ie is in undet. spec.
   end do
   
   nrud = sum(indic_e_undet)  ! number of elements in undetermined, unrepresented species
   
   allocate( i_zud( max( nzud, 1 ) ) )
   allocate( i_rud( max( nrud, 1 ) ) )
   
   if( lu >= 0 ) then
      write(lu,*)' '
      write(lu,*)' nzud, nrud = ', nzud, nrud 
      write(lu,*)' '
   endif
   
   ieu = 0
   do ie = 1, ne
      if( indic_e_undet(ie) == 1 ) then
         !  element in undetermined species
         ieu = ieu + 1         ! index of undetermined element
         i_rud( ieu ) = nrs+ie ! corresponding index in r
      endif
   end do
   
   if( lu >= 0 ) then
      write(lu,*)'indic_e_undet = ', indic_e_undet
   
      write(lu,*)' '
      write(lu,'(a,100i3)')' i_rud = ', i_rud
      write(lu,*)' '
   endif
   
!  Phase 2: form needed quantities based on regular species ordering

   r2zdet = 0.d0
   k_zdet = 0
   ius    = 0  ! index of unrepresented species
   iuds   = 0  ! index of unrepresented determined species
   
   if( lu >= 0 ) then
      write(lu,*)' '
      write(lu,*)' is, k_zdet(is), r2zdet*'
      write(lu,*)' '
   endif
   
   do is = 1, ns
   
      if( S2RS(is) > 0 ) then 
         k_zdet(is) = 1    ! represented species
         r2zdet( is, S2RS(is) ) = 1.d0
         
      else
         ius = ius + 1 ! index of unrepresented species
         if( indic_zu_det(ius) == 1 ) then  
            ! determined unrepresented species
            r2zdet( is, nrs+1:nrc) = zue_to_zu(ius,1:ne)

            do ie = 1, ne  !  count number of non-zero reduced compositions
               if( abs( zue_to_zu(ius,ie) ) > epsilon_rank ) k_zdet(is) = k_zdet(is) + 1
            end do
         else
            k_zdet(is)    = 0 ! undetermined species  
            iuds          = iuds + 1 
            i_zud( iuds ) = is 
         endif
      endif
      
      if( lu >= 0 ) write(lu,'(2i3,100i2)') is, k_zdet(is), nint( min( 9.d0,1d4*abs(r2zdet(is,:)) ))
      
   end do
   
   if( lu >= 0 ) then
      write(lu,*)' '
      write(lu,'(a,100i3)')' i_zud = ', i_zud
      write(lu,*)' '
   endif
     
  end subroutine determined
  
  subroutine normal  !--------------------------------------
  
!  Determine the normal vector on the boundary of the realizable region using QR
!  for the case in which r > 0 on the boundary.
!  Optionally (for if_svd=.true.) check using SVD.
!  The quantities determined are: r_normal, indic_zpos and indic_rpos.

   integer    :: lwork, info, rank, iur
   real(k_dp) :: BTP(ne,ne), tau(ne), work(ne*(20+ne)), err, &
                 Bn(ns), BTPc(ne,ne), U(ne,ne), S(ne), VT(ne,ne), svd_norm(ne)
   
  lwork = ne*(20+ne)  !  work space for SVD and QR

! represented species can be positive
  indic_zpos = 1
      
! take unrepresented species to be zero, until determined otherwise                 
   do i = 1, nus
      indic_zpos( US(i) ) = 0  
   end do
   
!  assemble BTP = columns of CEu^T corresponding to left-hand variables
!  the size of BTP is BTP(ne,ne-1)

   BTP  = 0.d0
   j    = 0
   
   do i = 1, nrud
      if( lhv(i) > nzud+1 ) then
         call abort( 4, 'bad left-hand variable', isv = lhv(i) )
         return
      endif
      if( lhv(i) > nzud ) cycle  !  lhv corresponds to ds
      
      iur = urud_to_ur( lhv(i) )   ! unrepresented-species index of LHV
      indic_zpos( US(iur) ) = 1    ! take species corresponding to LH variables to be positive

      j = j + 1
      BTP(:,j) = CEu( iur, : )
   end do
   
   if( j /= nrud-1 ) then
      call abort( 5, 'too few basic species', isv = j )
      return
   endif
   
!  add columns for determined elements
   iur = 1
   do i = 1, ne
      if( i_rud(iur) == nrs+i ) then
         !  undetermined element
         iur = min( iur + 1, nrud )
      else
         !  determined
         j = j + 1
         BTP(i,j) = 1.d0
      endif
   end do
   
!  optionally perform SVD to check rank and normal

   if( if_svd ) then
      BTPc = BTP  ! make copy and form BTPc = U * diag(S) * VT
      call dgesvd( 'A', 'N', ne, ne, BTPc, ne, S, U, ne, VT, ne, &
                   work, lwork, info )
                   
      if( info /= 0 ) then
         call abort( 6, 'dgesvd failed', isv = info ) 
         return
      endif
      
      rank = 0  !  check that rank is ne-1
      do i = 1, ne
         if( S(i) > epsilon_rank * S(1) ) rank = rank + 1
      end do
   
      if( rank /= ne-1 ) then
         call abort( 7, 'rank deficiency: rank=', isv = rank ) 
         return
      endif

      svd_norm = U(:,ne)  !  normal to facet in r-space
   endif
   
!  form QR decomposition Q * R = BTP,  Q is ne x ne

   lwork = ne*(2+ne)
   call dgeqrf( ne, ne-1, BTP, ne, tau, work, lwork, info )
   
   if( info /= 0 ) then
      call abort( 8, 'QR dgeqrf failed', isv = info ) 
      return
   endif
                          
   call dorgqr( ne, ne, ne-1, BTP, ne, tau, work, lwork, info )
   
   if( info /= 0 ) then
      call abort( 9, 'QR dorgqr failed', isv = info )
      return
   endif 
           
   if( if_svd ) then  !  check normal from QR against that from SVD
      err = abs( dot_product( svd_norm, BTP(:,ne) ) ) - 1.d0
      if( err > tol_use ) then
         call abort( 10, 'error in normal' )
         return
      endif 
   endif
           
   r_normal = 0.d0                       
   r_normal(nrs+1:nrc) = BTP(:,ne)  !  normal is the last column of Q
   if( dot_product( drdsm, r_normal ) < 0.d0 ) r_normal = -r_normal ! ensure outward pointing
   
!  determine species that can be positive
   Bn = matmul( BB, r_normal )
   do i = 1, ns
      if( indic_zpos(i) == 1 ) cycle  ! represented and LH unrepresented already set
      if( abs(Bn(i)) < 1.d-6 ) indic_zpos(i) = 1 ! dz(i) does not change n^T * dr 
   end do
   
 ! deduce indic_rpos from indic_zpos
   indic_rpos = nint( matmul( BBT, 1.d0*indic_zpos ) )
   do i = 1, nrc
      indic_rpos(i) = min( indic_rpos(i), 1 )
   end do

  return
  end subroutine normal
  
 
  subroutine abort( loc, mess, isv )  !---------------------------------------------------
  
   integer, intent(in)                :: loc     ! location number
   character(*), intent(in)           :: mess    ! message
   integer, intent(in), optional      :: isv
         
   if( quit ) then
      if( present(isv) ) then
         call isat_abort( 'ci_ice_pic_bound', loc, mess=mess, isv=isv )
      else
         call isat_abort( 'ci_ice_pic_bound', loc, mess=mess )
      endif
      
   else
      iflag = -loc
      return
   endif
   
   end subroutine abort
 
  end subroutine ci_ice_pic_bound  !----------------------------------------------------------
  
subroutine ci_ice_pic_bound_test

!  Test ci_ice_pic_bound.
!  Generate an internal point and a second point, and use ci_ice_pic_bound
!  to determine the intersection between the line connecting them and the boundary.
!  For ext2 == .true., the second point is external, and the intersection has r>0.

use ci_dat
use ci_dat8
use ci_utils
use linprog

logical    :: ext2 = .false.  ! set true to make the second point external
integer    :: npt = huge(0)  ! number of test points

logical    :: ext, int
integer    :: it, i, kfa, indic_rp(nrc), indic_zp(ns), iflag, info, lhv(2)
             
real(k_dp) :: r(nrc), ri(nrc), re(nrc), drds(nrc), z(ns), zmm, tol=1.d-10, &
              zi(ns), sb, zmmb, zmm_max = 0.d0, rmin, rmin_max=0.d0, ale(1,1), &
              age(1,1), ble(1), bge(1), beq(2), f(ns), x(ns),  r_normal(nrc), &
              rnd, rnd_max=0.d0, zvec(ns), zv0, aeq2(2,ns), zpi_min = huge(0.d0), &
              z0i_max = -huge(0.d0)

do it = 1, npt    !  loop over test points
   ext = .false.  ! indicate whether points have been generated
   int = .false.
   
   do
      call random_number( r )
      r = r / dot_product( r, amolwt_n ) !  random point, r>=0, w^T*r=1, but not necessarily realizable 
      
      call ceq_maxmin( ns, nrc, BB, r, z, zmm, info)  ! get max-min composition
      if( info < 0 ) call isat_abort('bound_test',1,mess='ceq_maxmin failed, info = ', isv = info )
      
      if( zmm < -tol ) then       
         ext = .true.  !  point is not realizable
         re = r
         
      elseif( zmm > tol ) then
                       !  point is realizable
         if( int  .and.  .not.ext2 ) then  
            re = r     ! 2nd internal point
            exit
         endif  
                  
         int = .true.  ! 1st (internal) point
         ri = r
         zi = z
      endif
      
      if( int  .and.  ext ) exit  ! both points generated
   end do
   
! the two points, ri and re have been generated
   
   drds = re - ri
   drds = drds / norm(drds) ! unit tangent vector
         
   call ci_ice_pic_bound( ri, drds, sb, kfa, indic_rp, indic_zp, r_normal, iflag )
         
   if( iflag < 0 ) call isat_abort( 'ci_ice_pic_bound_test', 1, &
                          mess='ci_ice_pic_bound failed, iflag = ', isv = iflag )
                                          
!  successful call to ci_ice_pic_bound: verify that r = ri + sb * drds is a boundary point    
   r = ri + sb * drds
   call ceq_maxmin( ns, nrc, BB, r, z, zmmb, info)
   if( info < 0 ) call isat_abort( 'ci_ice_pic_bound_test', 2, &
                            mess='ceq_maxmin failed, info = ', isv = info )
      
   if( zmmb > zmm_max ) then
      write(0,'(a,i12,i4,1p,10e13.4)') 'bound_test: record max z on facet =', it, kfa, zmmb
      zmm_max = zmmb
   endif
      
   rmin = minval(r)
   if( kfa <= nrc .and. abs(rmin) > rmin_max ) then
      write(0,'(a,i12,i4,1p,10e13.4)') 'bound_test: record max r on facet =', it, kfa, rmin
      rmin_max = abs(rmin)
   endif
   
!  test n=r_normal: maximize n^T r = m^T * z, subject to z >= 0, w^T z=1, where m = B * n. 
   zvec = matmul(BB,r_normal) ! minimize -n^T * r = -m^T * z
   f    = -zvec
   beq  = 1.d0
   
   call lp( ns, 0, 0, 1, ale, age, amolwt, ble, bge, beq, &
            f, x, iflag )
   
   if( iflag <0 ) call isat_abort('bound_test',3,mess='LP failed, iflag = ', isv = iflag )
   
   ri  = matmul(BBT,x)
   rnd = dot_product( r_normal, r-ri ) / norm(r)
   if( abs(rnd) > rnd_max ) then
      write(0,'(a,i12,i4,1p,10e13.4)') 'bound_test: record error in norml =', it, kfa, rnd
      rnd_max = abs(rnd)
   endif
   
! loop over species:  maximize z(i), subject to w^T * z = 1, zvec^T * z = zv0

   zv0 = dot_product( zvec, x )
   aeq2(1,:) = amolwt
   aeq2(2,:) = zvec
   beq(1)    = 1.d0
   beq(2)    = zv0
   
   do i = 1, ns
      f    = 0.d0
      f(i) = -1.d0

      call lp( ns, 0, 0, 2, ale, age, aeq2, ble, bge, beq, &
            f, x, iflag, lhvars=lhv )
            
     if( iflag <0 ) call isat_abort('bound_test',4,mess='LP failed, iflag = ', isv = iflag )
     
     if( indic_zp(i) == 1 ) then
        if( x(i) < zpi_min ) then
           write(0,'(a,i12,i4,1p,10e13.4)') 'bound_test: record low pos. spec. =', it, kfa, x(i)
           zpi_min = x(i)
        endif
     else
        if( x(i) > z0i_max ) then
           write(0,'(a,i12,i4,1p,10e13.4)') 'bound_test: record hi zero. spec. =', it, kfa, x(i)
           z0i_max = x(i)
        endif
     endif

   end do
    
end do  !  end of loop over test points

end subroutine ci_ice_pic_bound_test
