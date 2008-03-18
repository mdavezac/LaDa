subroutine iaga_set_mpi( _comm_handle, _rank )
  implicit none

  include "mpif.h"

  use mpigroup

  integer( intent = in ) _comm_handle, _rank

  rank = _rank
  _comm_handle = _handle

end subroutine
