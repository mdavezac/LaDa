!
! Version: $Id$
!

subroutine iaga_set_mpi( in_comm_handle )
 
  use mpigroup, only : comm_handle, irank, arank
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
subroutine iaga_set_cell( in_alat, in_cell, in_ng1, in_ng2, in_ng3 )
 
  use supercell

  implicit none
  ! cell parameter.
  real*8, intent(in) :: in_alat
  ! cell vectors.
  real*8, intent(in) :: in_cell(1:3,1:3)
  ! grid sizes.
  integer, intent(in) :: in_ng1, in_ng2, in_ng3

  alat = in_alat
  a1 = in_cell(1,1:3)
  a2 = in_cell(2,1:3)
  a3 = in_cell(3,1:3)
  ng1 = in_ng1
  ng2 = in_ng2
  ng3 = in_ng3
  ngb1 = in_ng1
  ngb2 = in_ng2
  ngb3 = in_ng3

  call set_cell()

end subroutine iaga_set_cell 

! Computes the observable <r> between valence and conduction bands.
subroutine iaga_dipole_elements( in_dipole, in_val, in_cond, &
                                 in_nval, in_ncond, &
                                 in_valfilename, in_vallength, &
                                 in_condfilename, in_condlength )

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
  complex*16, intent(inout) :: in_dipole( 1:3, 1:4, 1:in_nval, 1:in_ncond )
  ! Indices of the wavefunctions.
  integer, intent(in) :: in_val(1:in_nval)
  ! Indices of the wavefunctions.
  integer, intent(in) :: in_cond(1:in_ncond)
  ! Valence Filename length
  integer, intent(in) :: in_vallength
  ! Conduction Filename length
  integer, intent(in) :: in_condlength
  ! Valence Wavefunction Filename.
  character(len=in_vallength), intent(in) :: in_valfilename
  ! Conduction Wavefunction Filename.
  character(len=in_condlength), intent(in) :: in_condfilename

  ! Locals. 
  complex*16 :: qdpkr(1:4)
  complex*16 :: psi_c(1:ngrid,1:2,1:in_ncond)
  complex*16 :: psi_v(1:ngrid,1:2,1:in_nval)
  complex*16 :: dipole(1:3, 1:4, 1:in_nval, 1:in_ncond )
  integer, parameter :: fileunit = 17
  integer :: jc, iv, i, j, k

  ! Reads wavefunctions and performs calculations
  open ( unit=fileunit, file=in_valfilename, &
         form='UNFORMATTED', position='REWIND', status='OLD' )
    rewind(fileunit)
    write(*,*) in_nval, in_val, ngrid
    call read_wavefunctions( in_val, in_nval, psi_v, fileunit )
    if( in_valfilename .ne. in_condfilename ) then
      close(fileunit)
      open ( unit=fileunit, file=in_condfilename, &
             form='UNFORMATTED', position='REWIND', status='OLD' )
      rewind(fileunit)
    endif
    call read_wavefunctions( in_cond, in_ncond, psi_c, fileunit )
  close(fileunit)

  ! Performs calculations
  k = 0
  do jc=1, in_ncond
     do iv=1, in_nval
       call dqpole( psi_v(:,:,iv), psi_c(:,:,jc), in_dipole(:,:, iv,jc ), qdpkr )
     end do
  end do

  contains
    
    ! read all wavefunctions.
    subroutine read_wavefunctions( in_bands, in_n, in_psi, in_fileunit )

      use supercell
      
      implicit none
    
      ! Number of wavefunctions to read
      integer, intent(in) :: in_n
      ! Index of band wavefunctions.
      integer, intent(in) :: in_bands(1:in_n)
      ! Wavefunctions.
      complex*16, intent(inout) :: in_psi(1:ngrid,1:2,1:in_n)
      ! Fileunit.
      integer, intent(in) :: in_fileunit

      ! locals
      integer :: j = 1, i = 0, k, u
      real*8 :: psi(1:ngrid,1:2) 

      do while( j <= in_n )
        read(in_fileunit) psi(:,1)
        read(in_fileunit) psi(:,2)
        in_psi(:,1,j) = cmplx(psi(:,1), psi(:,2))
        read(in_fileunit) psi(:,1)
        read(in_fileunit) psi(:,2)
        in_psi(:,2,j) = cmplx(psi(:,1), psi(:,2))
        i = i + 1
        if( i .eq. in_bands(j) ) j = j + 1
      enddo

    end subroutine

end subroutine iaga_dipole_elements
