!
! Version: $Id$
!

subroutine iaga_set_mpi( comm_handle_ )
 
  use mpigroup
  implicit none
  include "mpif.h"

  integer, intent(in) :: comm_handle_

  integer inode, ierr

  comm_handle = comm_handle_

  call mpi_comm_rank(comm_handle,inode,ierr)
  call mpi_comm_rank(mpi_comm_world,irank,ierr)

  write(arank,'(I6)') inode
  arank = adjustl( arank )

end subroutine

subroutine iaga_call_escan( nbstates_ )
  use escan_comp_api, only: escancomp
  use eigenenergy
  use mpigroup
  implicit none
  include "mpif.h"

  integer, intent(in) :: nbstates_
  type ( escancomp ) ecp

  if( allocated( zebn ) )  deallocate( zebn )
  allocate( zebn( nbstates_ ) )

  ecp%comm_handle = comm_handle
  ecp%fileescaninput= "escan_input."//arank(1:len_trim(arank))
  ecp%escanfileonly=.TRUE.
  ecp%escandefaultprint=.false.
  ! The following do not need to be set for escanfileonly=.TRUE.
  ecp%filepot="pot.out"
  ecp%filewg_out = "wg.out"
  ecp%filewg_in = "wg.in"
  ecp%f_xatom = "atom.config"
  ecp%mx=1
  ecp%Eref=-2.43

  call escan_comp(ecp)

end subroutine

subroutine iaga_get_eigenvalues( states_, n_ )
  use eigenenergy
  use mpigroup
  implicit none
  include "mpif.h"

  integer, intent(in) :: n_
  real*8, intent(out), dimension(n_) :: states_
  integer i

  if( .not. allocated( zebn ) ) stop "Storage for eigenergies was never allocated."

  do i = 1, n_, 1
    states_(i) = zebn(i)*27.211396d0 ! goes to eV
  enddo

  deallocate( zebn )

end subroutine
