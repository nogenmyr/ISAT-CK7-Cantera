!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the x2f_mpi      !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

module x2f_mpi_grpdata
!
!
! Store the important datas used in adaptive strategies:
!   grp_index  - indicate which group the processor belongs to
!       g      - number of processors in each group 
!      ng      - number of groups in each partition
!   pt_index   - indicate which partition the processor belongs to
!   mpicomm_pt - the communicator used within partition
!   intracomm  - the communicator used within group
!      ips     - the i-th pairing stage
!      mps     - total number of pairing stage
!   pairstg    - pairing stage or not
!
implicit none
integer, save        :: grp_index, g, ng, pt_index, mpicomm_pt
integer, save        :: intracomm, intercomm, ips, mps, rank, set_pt
real, save           :: t_C
logical, save        :: pairstg
logical, save        :: impinit = .false., pairinit = .false.
end module x2f_mpi_grpdata
