!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ci_ice_pim_tan( A, S, T_PIM, svrat )

!  Determine an orthogonal basis (given by the columns of T_PIM) for the 
!  tangent space of the pre-image manifold (PIM).

!  Input:
!     A - reaction mapping sensitivity matrix
!     S - dz/dt at the pre-image point

!  Output:
!     T_PIM - span(T_PIM) is the tangent space of the pre-image manifold
!     svrat - ratio of singular values of TPI = [B^T * A  B^T * A * S]

use ci_dat
use ci_dat8
use ci_utils

real(k_dp), intent(in)  :: A(ns+1,ns+1), S(ns)
real(k_dp), intent(out) :: T_PIM(1:ns,1:ns+1-nrc), svrat

integer    :: lwork, info
real(k_dp) :: TPI(nrc,ns+1), Sig(nrc), U(1,1), VT(ns+1,ns+1), &
              work(10*(nrc+ns)+nrc*ns)

lwork = 10*(nrc+ns)+nrc*ns

!  For dx = [dz_g; dtau], infinitesimal increments in the PIM are given by:
!  TPI * dx = 0, where TPI = [B^T * A  B^T * A * S].  Thus, an orthogonal basis
!  for the PIM tangent space is given by T_PIM = V(1:ns,1:nrc+1:ns+1), 
!  where the SVD of TPI is:  TPI = U * Sig * V^T.

!  Assemble TPI
TPI(1:nrc,1:ns) = matmul( BBT, A(1:ns,1:ns) )
TPI(1:nrc,ns+1) = matmul( BBT, matmul(A(1:ns,1:ns), S) ) / norm(S)

!  Form SVD
call dgesvd( 'N', 'A', nrc, ns+1, TPI, nrc, Sig, U, 1, VT, ns+1, &
                 work, lwork, info )
                 
if( info /= 0 ) call isat_abort('ci_ice_pim_tan', 1, &
                                mess='SVD failed, info = ', isv = info )

svrat = Sig(nrc)/Sig(1)

T_PIM(1:ns,1:ns+1-nrc) = transpose( VT(nrc+1:ns+1,1:ns) )

return
end subroutine ci_ice_pim_tan
