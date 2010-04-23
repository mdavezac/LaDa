! A module to load/reload wavefunctions.
module Wfns_module
  use escan, only: t_Lattice, t_Escan
  implicit none
  include "mpif.h"
  
  ! Current lattice definition.
  type(t_Lattice), private :: lattice
  ! Current pescan definition.
  type(t_Escan), private :: params

  contains
    subroutine init(filename, mpicomm)
      use Escan, only: read_escaninput, read_potentialinput, set_common_blocks, &
                       prepare_fft_and_allocate_arrays, t_Lattice, t_Escan
  
      character(len=*), intent(in) :: filename            ! filename of escan input
      integer, intent(in) :: mpicomm                      ! mpi communicator
      ! Read parameters from file and sets common blocks
      call read_escaninput(filename, params, mpicomm)
      call read_potentialinput(params, lattice)
      call set_common_blocks(params, lattice)
      ! prepares wavefunctions.
      call prepare_fft_and_allocate_arrays(params, lattice)

    end subroutine init

    ! Reads wavefunctions from current directory.
    subroutine read_wavefunctions(indices, wfns, gvecs, projs)
      use Escan, only: read_escaninput, read_potentialinput, set_common_blocks, &
                       prepare_fft_and_allocate_arrays, pi
      use data, only : mg_nx, wg_n
      use load_data, only : ngtotnod, n1p_n, n2p_n, n3p_n
      ! indices of the wavefunctions.
      integer, dimension(:), intent(in) :: indices
      ! output wavefunctions
      complex(kind=8), dimension(:,:,:), intent(out) :: wfns
      ! output g vectors.
      real(kind=8), dimension(:,:), intent(out) :: gvecs
      ! output projectors (smooth cutoff)
      real(kind=8), dimension(:), intent(out) :: projs
     
      
      real(kind=8) ::  vol
      integer n1, n2, n3, ng, ng_n, mx, nr
      common /com123/n1,n2,n3,ng,ng_n,nr,mx,vol
      integer inode,nnodes
      common /mpi_data/inode,nnodes
      integer(kind=8), dimension(3) :: int_gpoints
      integer :: ig, nb_gpoints
 
      ! now reads actual data.
      if( params%with_spinorbit ) then
        call read_wg_comp( wfns(:,:,1), wfns(:,:,2), &
                           n1,n2,n3,size(indices),mg_nx, &
                           params%Ecutoff, lattice%rcell, &
                           params%ecp%filewg_out, size( indices ), &
                           indices(1), 1 )
      else
        call read_wg_comp( wfns(:,:,1), wfns(:,:,1), &
                           n1,n2,n3,size(indices),mg_nx, &
                           params%Ecutoff, lattice%rcell, &
                           params%ecp%filewg_out, size( indices ), &
                           indices, 0 )
      endif

      nb_gpoints = ngtotnod( inode )
      projs = wg_n(1:nb_gpoints) ! copy projector data
      ! Computes all gpoints.
      do ig = 1, nb_gpoints
        ! compute G point.
        int_gpoints(1) = n1p_n(ig) - 1
        int_gpoints(2) = n2p_n(ig) - 1
        int_gpoints(3) = n3p_n(ig) - 1
      
        if( int_gpoints(1) .gt. n1 / 2 ) int_gpoints(1) = int_gpoints(1) - n1
        if( int_gpoints(2) .gt. n2 / 2 ) int_gpoints(2) = int_gpoints(2) - n2
        if( int_gpoints(3) .gt. n3 / 2 ) int_gpoints(3) = int_gpoints(3) - n3
      
        gvecs(ig, 1) &
          = 2.d0 * pi * sum( lattice%kcell(1,:) * int_gpoints )
        gvecs(ig, 2) &
          = 2.d0 * pi * sum( lattice%kcell(2,:) * int_gpoints )
        gvecs(ig, 3) & 
          = 2.d0 * pi * sum( lattice%kcell(3,:) * int_gpoints )
      enddo ! ig
 
    end subroutine read_wavefunctions

    subroutine get_array_dimensions( n0, n2, g0 )
      use data, only: mg_nx
      use load_data, only: ngtotnod
      integer, intent(out) :: n0, n2, g0

      integer inode,nnodes
      common /mpi_data/inode,nnodes
      
      n0 = mg_nx
      n2 = 2
      if( .not. params%with_spinorbit ) n2 = 1
      g0 = ngtotnod(inode)
    end subroutine get_array_dimensions

    subroutine wfns_cleanup
      use escan, only: escan_cleanup
      call escan_cleanup
    end subroutine wfns_cleanup
end module Wfns_module

