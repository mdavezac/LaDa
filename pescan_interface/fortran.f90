!
! Version: $Id$
!
subroutine iaga_set_mpi( comm_handle_, rank_ )

  use mpigroup
  implicit none
  include "mpif.h"


        integer, intent(in) :: comm_handle_
        integer, intent(in) :: rank_

  irank = rank_
  comm_handle = comm_handle_

end subroutine
