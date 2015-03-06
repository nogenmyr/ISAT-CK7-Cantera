!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module isat_mpi_type

!  buffers and other data used in mpi message passing

use isat_prec
implicit none

!------------------------------------------------------------------------------
type :: mpi_type
   integer :: comm, nprocm, kxf, pend_max, &
              ngsb, gsb_j, gsb_tag, nasb, asb_j, asb_tag, asb_bytes

   real(k_xf), pointer :: gsb(:,:), gsb_to(:), gsb_from(:), gsb_to_recd(:), &
                          grb(:), asb_to(:), asb_from(:), asb_to_recd(:), &
			              asb(:,:), arb(:)
   integer, pointer    :: gsb_req(:), asb_req(:,:)
end type mpi_type
!------------------------------------------------------------------------------

end module isat_mpi_type
