! A module to load/reload wavefunctions.
module Wfns_module

  ! Current lattice definition.
  type(t_Lattice), private :: lattice
  ! Current pescan definition.
  type(t_Escan), private :: params

  contains
    subroutine wfns_module_init(filename, mpicomm)
      use Escan, only: read_escaninput, read_potentialinput, set_common_blocks, &
                       prepare_fft_and_allocate_arrays, t_Lattice, t_Escan
  
      character(len=*), intent(in) :: filename           ! filename of escan input
      integer, intent(in) :: mpicomm                      ! mpi communicator
      ! Read parameters from file and sets common blocks
      call read_escaninput(filename, params, mpicomm)
      call read_potentialinput(params, lattice)
      call set_common_blocks(params, lattice)
      ! Sets number of spins
      nb_spins = 1
      if( params%with_spinorbit .eqv. .true. ) nb_spins = 2;
      ! prepares wavefunctions.
      call prepare_fft_and_allocate_arrays(params, lattice)

    end subroutine wfns_module_init

    ! Reads wavefunctions from current directory.
    subroutine read_wavefunctions(indices, wfns, gvecs, mpicomm)
      use Escan, only: read_escaninput, read_potentialinput, set_common_blocks, &
                       prepare_fft_and_allocate_arrays, t_Lattice, t_Escan
      use data, only : mg_nx
      use load_data, only : ngtotnod, n1p_n, n2p_n, n3p_n
      ! Filename of escan input (not wavefunction input!)
      character(len=*), intent(in) :: filename
      ! indices of the wavefunctions.
      integer, dimension(:), intent(in) :: indices
      ! output wavefunctions
      complex(kind=8), dimension(:,:,:), intent(out) :: wfns
      ! output g vectors.
      real(kind=8), dimension(:,:), intent(out) :: gvecs
      ! Group communicator.
      integer, intent(in) :: mpicomm
     
      
      real(kind=8) ::  vol
      integer n1, n2, n3, ng, ng_n, mx, nr
      common /com123/n1,n2,n3,ng,ng_n,nr,mx,vol
      integer inode,nnodes
      common /mpi_data/inode,nnodes
      integer(kind=8), dimension(3) :: int_gpoints
 
      ! now reads actual data.
      if( params%with_spinorbit ) then
        call read_wg_comp( wavefunctions%values(:,:,1), wavefunctions%values(:,:,2), &
                           n1,n2,n3,size(indices),mg_nx, &
                           params%Ecutoff, lattice%rcell, &
                           params%ecp%filewg_out, size( indices ), &
                           indices(1), 1 )
      else
        allocate( wavefunctions%values(mg_nx, size(indices), 1 ) )
        call read_wg_comp( wavefunctions%values(:,:,1), wavefunctions%values(:,:,2), &
                           n1,n2,n3,size(indices),mg_nx, &
                           params%Ecutoff, lattice%rcell, &
                           params%ecp%filewg_out, size( indices ), &
                           indices, 0 )
      endif

      nb_gpoints = ngtotnod( inode )
      allocate( wavefunctions%gpoints( nb_gpoints, 3 ) )
      ! Computes all gpoints.
      do ig = 1, nb_gpoints
        ! compute G point.
        int_gpoints(1) = n1p_n(ig) - 1
        int_gpoints(2) = n2p_n(ig) - 1
        int_gpoints(3) = n3p_n(ig) - 1
      
        if( int_gpoints(1) .gt. n1 / 2 ) int_gpoints(1) = int_gpoints(1) - n1
        if( int_gpoints(2) .gt. n2 / 2 ) int_gpoints(2) = int_gpoints(2) - n2
        if( int_gpoints(3) .gt. n3 / 2 ) int_gpoints(3) = int_gpoints(3) - n3
      
        gpoints(ig, 1) &
          = 2.d0 * pi * sum( lattice%kcell(1,:) * int_gpoints )
        gpoints(ig, 2) &
          = 2.d0 * pi * sum( lattice%kcell(2,:) * int_gpoints )
        gpoints(ig, 3) & 
          = 2.d0 * pi * sum( lattice%kcell(3,:) * int_gpoints )
      enddo ! ig
 
    end subroutine read_wavefunctions

end module Wfns_module
