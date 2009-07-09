!
! Version: $Id$
!

subroutine iaga_set_mpi( in_comm_handle )
 
  use mpigroup, only : comm_handle, irank, arank
  use MomentumDipole, only : params
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

subroutine iaga_call_escan( in_nbstates )
  use escan_comp_api, only: escancomp
  use eigenenergy
  use mpigroup
  implicit none
  include "mpif.h"

  integer, intent(in) :: in_nbstates
  type ( escancomp ) ecp

  if( allocated( zebn ) )  deallocate( zebn )
  allocate( zebn( in_nbstates ) )

  ecp%comm_handle = comm_handle
  ecp%fileescaninput= trim("escan_input.")//arank(1:len_trim(arank))
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
    states_(i) = zebn(i)*27.211396d0 ! goes to eV from Hartree units.
  enddo

  deallocate( zebn )

end subroutine

! Interface to subroutine in MomentumDipole module for easy C access.
subroutine momentum( in_inputfilename, in_fsize, &
                     in_dirvalence, in_dsizeA, &
                     in_dirconduction, in_dsizeB, &
                     in_indicesA, in_indicesB, & 
                     in_bandsA, in_bandsB, &
                     io_dipoles, in_dip2, in_dip3 )

  use MomentumDipole, only : module_momentum => momentum

  ! size of input filename
  integer, intent(in) :: in_fsize
  ! input filename. For folded spectrum calculations, only one input is
  ! needed, wether vbm or cbm. It is implicit that these calculations differ
  ! only by the reference energy. 
  character( len= in_fsize ), intent(in) :: in_inputfilename
  ! size of directory filename A.
  integer, intent(in) :: in_dsizeA
  ! wfn filename A. 
  character( len= in_dsizeA ), intent(in) :: in_dirvalence
  ! size of directory filename B.
  integer, intent(in) :: in_dsizeB
  ! wfn filename B.
  character( len= in_dsizeB ), intent(in) :: in_dirconduction

  ! number of bands A.
  integer, intent(in) :: in_bandsA
  ! band/spin indices A.
  integer, intent(in) :: in_indicesA( in_bandsA )
  ! number of bands B.
  integer, intent(in) :: in_bandsB
  ! band/spin indices B.
  integer, intent(in) :: in_indicesB( in_bandsB )
  ! number of valence dipole 
  integer, intent(in) :: in_dip2 
  ! number of conduction dipole 
  integer, intent(in) :: in_dip3 
  ! Dipole output.
  complex( kind=8 ),  intent(inout) :: io_dipoles( 3, in_dip2, in_dip3 )

  call module_momentum( in_inputfilename, in_fsize, &
                        in_dirvalence, in_dsizeA, &
                        in_dirconduction, in_dsizeB, &
                        in_indicesA, in_indicesB, & 
                        in_bandsA, in_bandsB, &
                        io_dipoles, in_dip2, in_dip3 )


end subroutine
