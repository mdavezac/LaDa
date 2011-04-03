subroutine iaga_set_mpi( in_comm_handle )
 
  use mpigroup, only : comm_handle, irank, arank
  use escan, only : params
  implicit none
  include "mpif.h"

  integer, intent(in) :: in_comm_handle

  integer inode, ierr, nnodes
  common /mpi_data/inode,nnodes

  comm_handle = in_comm_handle
  params%ecp%comm_handle = in_comm_handle

  call mpi_comm_rank(comm_handle,inode,ierr)
  call mpi_comm_rank(mpi_comm_world,irank,ierr)
  call mpi_comm_size(comm_handle,nnodes,ierr)

  write(arank,'(I6)') irank
  arank = adjustl( arank )

end subroutine

subroutine iaga_just_call_escan()
  use escan_comp_api, only: escancomp
  use eigenenergy
  use mpigroup
  implicit none
  include "mpif.h"

  type ( escancomp ) ecp

  ecp%comm_handle = comm_handle
  ecp%fileescaninput= "escan_input" !trim("escan_input.")//arank(1:len_trim(arank))
  ecp%escanfileonly=.TRUE.
  ecp%escandefaultprint=.true.
  ! The following do not need to be set for escanfileonly=.TRUE.
  ecp%filepot="pot.out"
  ecp%filewg_out = "wg.out"
  ecp%filewg_in = "wg.in"
  ecp%f_xatom = "atom.config"
  ecp%mx=1
  ecp%Eref=-2.43

  call escan_comp(ecp)

end subroutine

subroutine iaga_call_escan( in_nbstates, in_verbose, n, filename )
  use escan_comp_api, only: escancomp
  use eigenenergy
  use mpigroup
  implicit none
  include "mpif.h"

  integer, intent(in) :: in_nbstates
  integer, intent(in) :: in_verbose
  integer, intent(in) :: n ! size of filename string
  character(len=n), intent(in) :: filename
  type ( escancomp ) ecp

  if( allocated( zebn ) )  deallocate( zebn )
  allocate( zebn( in_nbstates ) )

  ecp%comm_handle = comm_handle
  ecp%fileescaninput= filename
  ecp%escanfileonly=.TRUE.
  ecp%escandefaultprint=.true.
  if( in_verbose == 0 ) ecp%escandefaultprint=.false.
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

  integer, intent(inout) :: n_
  real*8, intent(out), dimension(n_) :: states_
  integer i

  if( .not. allocated( zebn ) ) then
    n_ = 0
    return 
  endif 

  if( n_ .eq. 0 ) then
    n_ = size(zebn)
    return
  endif 

  i = min(n_, size(zebn))
  states_(1:i) = zebn(1:i)*27.211396d0 ! goes to eV from Hartree units.
  deallocate( zebn )

end subroutine


! prepares to read wavefunctions
subroutine escan_wfns_init(n, filename, latscale, smooth, &
                           kinscal, kpoint, pottype, comm)
  use wfns_module, only: init
  integer, intent(in) :: n                           ! length of filename
  character(len=n), intent(in) :: filename           ! filename of escan input
  real(kind=8), intent(in) :: latscale, smooth, kinscal, kpoint(3)
  integer, intent(in) :: comm, pottype            ! mpi communicator
  call init(filename, latscale, smooth, kinscal, kpoint, pottype, comm)
end subroutine escan_wfns_init
! cleanup module
subroutine escan_wfns_cleanup
  use wfns_module, only: wfns_cleanup
  call wfns_cleanup
end subroutine escan_wfns_cleanup

! gets array dimensions
subroutine escan_wfns_get_array_dimensions( n0, n2, g0 )
  use wfns_module, only: get_array_dimensions
  integer, intent(out) :: n0, n2, g0
  call get_array_dimensions(n0, n2, g0)
end subroutine escan_wfns_get_array_dimensions

! Reads wavefunction with given index
subroutine escan_wfns_read(n0, n1, n2, g0, indices, wfns, gvecs, projs, inverse)
  use wfns_module, only: read_wavefunctions
  implicit none

  integer, intent(in) :: n0, n1, n2, g0               ! all dimensions.
  integer, dimension(n1), intent(in) :: indices       ! indices to wavefunctions
  ! output wavefunctions
  complex(kind=8), dimension(n0, n1, n2), intent(out) :: wfns
  ! output g vectors.
  real(kind=8), dimension(g0,3), intent(out) :: gvecs
  ! output projects (g-space smooth cutoff)
  real(kind=8), dimension(g0), intent(out) :: projs
  ! output projects (g-space smooth cutoff)
  integer, dimension(g0), intent(out) :: inverse

  call read_wavefunctions(indices, wfns, gvecs, projs, inverse)

end subroutine escan_wfns_read

! Reads wavefunction with given index. No krammer degeneracy.
subroutine escan_wfns_read_nokram(n0, n1, n2, g0, indices, wfns, gvecs, projs)
  use wfns_module, only: read_wavefunctions
  implicit none

  integer, intent(in) :: n0, n1, n2, g0               ! all dimensions.
  integer, dimension(n1), intent(in) :: indices       ! indices to wavefunctions
  ! output wavefunctions
  complex(kind=8), dimension(n0, n1, n2), intent(out) :: wfns
  ! output g vectors.
  real(kind=8), dimension(g0,3), intent(out) :: gvecs
  ! output projects (g-space smooth cutoff)
  real(kind=8), dimension(g0), intent(out) :: projs
  ! output projects (g-space smooth cutoff)
  integer, dimension(1) :: inverse

  call read_wavefunctions(indices, wfns, gvecs, projs, inverse)

end subroutine escan_wfns_read_nokram

