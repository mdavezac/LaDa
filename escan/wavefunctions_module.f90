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
    subroutine init(filename, latscale, smooth, kinscal, kpoint, pottype, comm)
      use data, only : mg_nx, gkk_n, wg_n, inv_n
      use fft_data, only : fft_allocate, ncolx
      use load_data, only: nr1x, nr2x, nr3x
      use Escan, only: set_common_blocks, prepare_fft_and_allocate_arrays, pi
  
      character(len=*), intent(in) :: filename            ! filename of escan input
      integer, intent(in) :: comm, pottype             ! mpi communicator
      real(kind=8), intent(in) :: latscale, smooth, kinscal, kpoint(3)

      integer, parameter :: inputunit = 9
      real(kind=8) ::  AL(3,3), Ecut, vol, delta_k, volume
      integer ikpt, inode, nnodes, n1, n2, n3, ng, ng_n, mx, nr, ierr, totg
      integer if_so, mr, i
      common /com123/n1,n2,n3,ng,ng_n,nr,mx,vol
      common /mpi_data/inode,nnodes
      common /comAD/AL,Ecut
      common /comikpt/ikpt

      call iaga_set_mpi(comm)
      params%ecp%comm_handle = comm
      params%kpoint = kpoint * 2e0 * pi / latscale
      params%KineticScaling = kinscal
      params%PotentialType = pottype
      params%Gsmooth = smooth
      params%ecp%filewg_out = filename

      ! Read parameters from file and sets common blocks
      if(inode .eq. 0) then
        open(inputunit,file=filename,form='unformatted',status='old')
          rewind(inputunit)
          read(inputunit) (params%mesh(i), i=1,3), params%ecp%mx,if_so,ikpt
          read(inputunit) params%Ecutoff
          read(inputunit) lattice%rcell
        close(inputunit)
      endif
      call mpi_bcast(params%mesh, 3, MPI_INTEGER, 0, comm, ierr)
      call mpi_bcast(if_so, 1, MPI_INTEGER, 0, comm, ierr)
      call mpi_bcast(ikpt, 1, MPI_INTEGER, 0, comm, ierr)
      call mpi_bcast(params%ecp%mx, 1, MPI_INTEGER, 0, comm, ierr)
      call mpi_bcast(params%Ecutoff, 1, MPI_REAL8, 0, comm, ierr)
      call mpi_bcast(lattice%rcell, 9, MPI_REAL8, 0, comm, ierr)
      params%with_spinorbit = .true.
      if(if_so .eq. 0) params%with_spinorbit = .false.
      params%is_gamma = .false.
      if(ikpt .eq. 0) then 
        params%is_gamma = .true.
        params%kpoint = 0e0
      endif

      call set_common_blocks(params, lattice)
      call prepare_fft_and_allocate_arrays(params, lattice)

    end subroutine init

    ! Reads wavefunctions from current directory.
    subroutine read_wavefunctions(indices, wfns, gvecs, projs, inverse)
      use Escan, only: pi
      use data, only : mg_nx, wg_n, inv_n
      use load_data, only : ngtotnod, n1p_n, n2p_n, n3p_n
      ! indices of the wavefunctions.
      integer, dimension(:), intent(in) :: indices
      ! output wavefunctions
      complex(kind=8), dimension(:,:,:), intent(out) :: wfns
      ! output g vectors.
      real(kind=8), dimension(:,:), intent(out) :: gvecs
      ! output projectors (smooth cutoff)
      real(kind=8), dimension(:), intent(out) :: projs
      ! indices to -G components
      integer, dimension(:), intent(out) :: inverse
     
      
      real(kind=8) ::  vol
      integer n1, n2, n3, ng, ng_n, mx, nr
      common /com123/n1,n2,n3,ng,ng_n,nr,mx,vol
      integer inode,nnodes
      common /mpi_data/inode,nnodes
      integer(kind=8), dimension(3) :: int_gpoints
      integer :: ig, nb_gpoints
      real(kind=8) a
      integer ierr
 
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
      if(allocated(inv_n)) inverse = inv_n(1:nb_gpoints) - 1 ! copy -G indices (starting at 0)
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

