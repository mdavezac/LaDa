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

! Computes the dipole elements between valence and conduction bands.
subroutine iaga_dipole_elements( _dipole, _nval, _ncond )

  use matelements
  use supercell
  use wavefunction
  use strings
  
  implicit none
  ! _dipole must have given format: 
  !   1:3 corresponds to the three cartesian coordiantes.
  !   1:4 corresponds to the four possible spin arrangements
  !   1:_nval corresponds to the valence bands
  !   1:_ncond corresponds to the conduction bands 
  ! And it must be an array of 64bit fortan complexes.
  complex*16, intent(inout) :: _dipole(1:3, 1:4, 1:_nval, 1:_ncond )
  ! Number of valence bands.
  integer, intent(in) :: _nval
  ! Number of conduction bands.
  integer, intent(in) :: _ncond
  ! Locals. 
  complex*16 :: dipkr(1:3,1:4)
  complex*16 :: qdpkr(1:4)
  complex*16, allocatable :: psi_j(:,:)
  complex*16, allocatable :: psi_v(:,:)

  ! Filenames
  allocate( engv(1:_nval) )
  allocate( filev(1:_nval) )
  allocate( strgv(1:_nval) )
  allocate( engc(1:_ncond) )
  allocate( filec(1:_ncond) )
  allocate( strgc(1:_ncond) )

  open (unit=4,file='mxmat.d',status='old')

     read (4,*) alat
     read (4,*) a1(1),a1(2),a1(3) 
     read (4,*) a2(1),a2(2),a2(3) 
     read (4,*) a3(1),a3(2),a3(3) 
     read (4,*) ng1,ng2,ng3
     read (4,*) ngb1,ngb2,ngb3
     do iv=1, _nval
        read (4,*) engv(iv),filev(iv),strgv(iv)
     end do
     allocate (engc(1:ncon))
     allocate (filec(1:ncon))
     allocate (strgc(1:ncon))
     do ic=1, _ncon
        read (4,*) engc(ic),filec(ic),strgc(ic)
     end do

  close (unit=4)

  call set_cell()

  allocate( psi_i(1:ngrid, 1:2) )
  allocate( psi_j(1:ngrid, 1:2) )
  
  ! Reads wavefunctions and performs calculations
  do jc=1, _ncon
     call readwf ('c',jc,psi_j)
     do iv=1, _nval
        call readwf( 'v', iv, psi_i )
        call dqpole( psi_i(1), psi_j(1), _dipole(1,1, iv,jc ), qdpkr(1,1) )
     end do
  end do

  deallocate( psi_i )
  deallocate( psi_j )
  deallocate( engc )
  deallocate( filec )
  deallocate( strgc )
  deallocate( engv )
  deallocate( filev )
  deallocate( strgv )

end
