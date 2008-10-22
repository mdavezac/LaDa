!
! Version: $Id$
!

subroutine iaga_set_mpi( in_comm_handle )
 
  use mpigroup
  implicit none
  include "mpif.h"

  integer, intent(in) :: in_comm_handle

  integer inode, ierr

  comm_handle = in_comm_handle

  call mpi_comm_rank(comm_handle,inode,ierr)
  call mpi_comm_rank(mpi_comm_world,irank,ierr)

  write(arank,'(I6)') inode
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

! Sets stuff for supercell module.
subroutine iaga_set_cell( in_alat, in_cellin_ng1, in_ng2, in_ng3 )
 
  use supercell

  implicit none
  ! cell parameter.
  real*8, intent(in) :: in_alat
  ! cell vectors.
  real*8, intent(in) :: in_cell(1:3,1:3)
  ! grid sizes.
  integer, intent(in) :: in_ng1, in_ng2, in_ng3

  a1 = in_cell(1,1:3)
  a2 = in_cell(2,1:3)
  a3 = in_cell(3,1:3)
  ng1 = in_ng1
  ng2 = in_ng2
  ng3 = in_ng3

  call set_cell()

end subroutine iaga_set_cell 

! Computes the dipole elements between valence and conduction bands.
subroutine iaga_dipole_elements( in_dipole, in_bands, in_nval, &
                                 in_ncond, in_filename, in_length )

  use matelements
  use supercell
  
  implicit none
  ! Number of valence bands.
  integer, intent(in) :: in_nval
  ! Number of conduction bands.
  integer, intent(in) :: in_ncond
  ! in_dipole must have given format: 
  !   1:3 corresponds to the three cartesian coordinates.
  !   1:4 corresponds to the four possible spin arrangements
  !   1:in_nval corresponds to the valence bands
  !   1:in_ncond corresponds to the conduction bands 
  ! And it must be an array of 64bit fortan complexes.
  complex*16, intent(inout) :: in_dipole(1:3, 1:4, 1:in_nval, 1:in_ncond )
  ! Indices of the wavefunctions.
  integer, intent(in) :: in_bands(1:in_nval+in_ncond)
  ! Filename length
  integer, intent(in) :: in_length
  ! Wavefunction Filename.
  character, intent(in) :: in_filename(1:in_length)

  ! Locals. 
  complex*16 :: qdpkr(1:4)
  complex*16, allocatable :: psi_c(:,:)
  complex*16, allocatable :: psi_v(:,:)
  integer jc, iv

  allocate( psi_c(1:ngrid, 1:2) )
  allocate( psi_v(1:ngrid, 1:2) )
  
  ! Reads wavefunctions and performs calculations
  open (unit=4,file=filename,form='unformatted')
  do jc=1, in_ncond
     call read_wavefunction( in_bands( in_nval + jc ), psi_c, fileunit )
     do iv=1, in_nval
        call read_wavefunction( in_bands( iv ), psi_c, fileunit ) 
        call dqpole( psi_v, psi_c, in_dipole(:,:, iv,jc ), qdpkr )
     end do
  end do
  close(4)

  deallocate( psi_c )
  deallocate( psi_v )

  contains
    
    ! Reads wavefunction from file.
    subroutine read_wavefunction( in_n, in_psi, fileunit )

      use supercell

      implicit none
      ! band index.
      integer, intent(in) :: in_n
      ! Wavefunction to be read.
      complex*16, intent(inout) :: in_psi(1:ngrid, 1:2 )
      ! File unit from which to read.
      integer, intent(in) :: fileunit


      read (4, POS=2*ngrid*in_n+1) psi(:,1) 
      read (4) psi(:,2) 

    end subroutine read_wavefunction

end subroutine iaga_dipole_elements
