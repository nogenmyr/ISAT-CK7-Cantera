!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!


      subroutine ci_ice_facet( Nbar_pre, z_pre, X_pre, Tlambda_pre, 
     1                         TNbar_pre, k_facet, alpha_min,
     2                         class, z0_index,
     3                         indic_pos, z_lin, ns_pos  )
        
! The ray r(alpha) = r_pre + alpha * dr/ds  intersects the boundary of the 
! realizable region in r-space on facet k_facet for alpha = alpha_min. 
! This routine returns k_facet and alpha_min.  
! Thus, r_j(alpha_min) >= 0 for all j, with the equality holding for j=k_facet.

! Input:
!    Nbar_pre:     total specific moles at pre-image point
!    z_pre:        specific moles of species at pre-image point
!    X_pre:        mole fractions at pre-image point
!    Tlambda_pre:  d(lambda)/ds on PIC

!  Output:
!    k_facet:      index of facet intersected by ray
!    alpha_min:    value of alpha at intersection

!  Output - temporary for consistency with original edge_face
!           class, z0_index, indic_pos, z_lin, alpha_lin, ns_pos

!  Note: it is not unusual for r_pre to be very close to a facet
!        (e.g., r_pre(i)<1d-18) and yet the intersection is with a
!        different facet. 

      use ci_dat
      use ci_dat8
      implicit none

	real(k_dp), intent(in)  :: Nbar_pre, z_pre(ns), X_pre(ns), 
     1                           Tlambda_pre(nrc), TNbar_pre

	integer, intent(out)    :: k_facet
      real(k_dp), intent(out) :: alpha_min

! temporary
	integer, intent(out):: class, z0_index, indic_pos(ns), ns_pos     
      real(k_dp), intent(out) :: z_lin(ns)   


! Local variables
      integer    :: i, is, ie
      real(k_dp) :: dzds(ns), drds(nrc), r_pre(nrc), alpha     
   
! Form the tangent vector dz/ds of the PIC in z-space.      
! Form the direction in the full space (vector v in Eq.48)
! dzds= Nbar_pre*diag_X_pre*BBF(1:ns,1:nrc)*Tlambda_pre+X_pre*TNbar_pre;
      
       dzds = X_pre * 
     1            ( Nbar_pre * matmul( BBF, Tlambda_pre ) + TNbar_pre )
     
!  Form the tangent vector dr/ds of the PIC in r space 
       drds = matmul( BBT, dzds )
       
       r_pre = matmul( BBT, z_pre )
       
       alpha_min = huge(1.d0)
       k_facet   = 0
       
       do i = 1, nrc
          if( drds(i) < 0.d0 )then
             alpha = -r_pre(i) / drds(i)
             if( alpha < alpha_min ) then
                k_facet   = i
                alpha_min = alpha
             endif
          endif
       end do
       
       if( k_facet == 0 ) call isat_abort('ci_ice_Facet', 1,
     1      mess='Intersection not found, drds=', rvar=drds )
     
     
!  evaluate temporary quantities ----------------------------

! z_lin may not be realizable, but r formed from it is
       z_lin = z_pre + alpha_min * dzds 
       
       if( k_facet <= nrs ) then  !  represented species is zero
          class               = 0
          z0_index            = CS(k_facet)
          indic_pos           = 1
          indic_pos(z0_index) = 0
          ns_pos              = ns - 1
          
       else  ! elements in unrepresented species is zero
          class               = 1
          ie                  = k_facet - nrs ! index of zero element
          
          indic_pos           = 1  !  set all species positive
          do is = 1, ns            !  find unrepresented species 
             if( indic_rs(is) == 0 ) then ! containing zero element
                if( CE(is,ie) > 0.d0 ) indic_pos(is) = 0
             endif
          end do
          
          z0_index            = 0
          ns_pos              = sum( indic_pos )
        endif
       
       return
       
       end subroutine ci_ice_facet
