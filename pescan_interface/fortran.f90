!
! Version: $Id$
!
subroutine iaga_call_genpot( comm_handle_, rank_ )

  use mpigroup
  implicit none
  include "mpif.h"


  integer(KIND=MPI_OFFSET_KIND), intent(in) :: comm_handle_
  integer, intent(in) :: rank_

  irank = rank_
  write(arank,'(I6)') irank
  arank = adjustl( arank );
  comm_handle = comm_handle_

  call getVLarg()

end subroutine

subroutine iaga_call_escan()
  use escan_comp_api
  use mpigroup
  implicit none
  include "mpif.h"

  type ( escancomp ) ecp

  ecp%comm_handle = comm_handle
  ecp%fileescaninput= "escan_input."//arank(1:len_trim(arank))
  ecp%escanfileonly=.TRUE.
  ecp%escandefaultprint=.TRUE.
  ! The following do not need to be set for escanfileonly=.TRUE.
  ecp%filepot="pot.out"
  ecp%filewg_out = "wg.out"
  ecp%filewg_in = "wg.in"
  ecp%f_xatom = "atom.config"
  ecp%mx=1
  ecp%Eref=-2.43

  call escan_comp(ecp)

end subroutine
