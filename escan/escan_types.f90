! Module to hold a few escan general stuff
Module Escan

  use escan_comp_api, only: escancomp
  implicit none
  include 'mpif.h'

  public

  type t_Escan ! Could be an extension of escancomp... waiting for fortran 2003

    type( escancomp ) :: ecp 
    ! Energy Cutoff (in atomic units at end of read_input )
    real( kind=8 ) :: Ecutoff
    ! Method used in computation: complete diagonalization or folding.
    integer :: method
    ! G-vector smoothing parameter
    real( kind=8 ) :: Gsmooth
    ! Kinetic energy scaling factor
    real( kind=8 ) :: KineticScaling
    ! Kpoint. 
    real( kind=8 ) :: kpoint(3)
    ! Lattice scaling. 
    real( kind=8 ) :: lattice_scale
    ! Wether gamma calculations are intended.
    logical :: is_gamma 
    ! Type of the potential:
    !   1 = local, 2 = non-local, 3 = non-local + spin-orbit
    integer :: PotentialType
    ! Wether the potential has spin orbit. Equivalent to PotentialType = 3.
    logical :: with_spinorbit
    ! Mesh parameters.
    integer, dimension(3) :: mesh

  end type ! t_Escan

  type t_Lattice

    ! real-space cell vectors
    real(kind=8) :: rcell(3,3)
    ! reciprocal-space cell vectors.
    real(kind=8) :: kcell(3,3)

  endtype
  real( kind=8 ), parameter, public :: pi = 3.14159265358979323846264338327948d0
  ! Escan parameters.
  type( t_Escan ), public :: params 

  contains

    ! reads escan input necessary for momentum dipole computations
    subroutine read_escaninput( filename, escan_input, comm )
      ! input filename
      character(len=*), intent(in) :: filename
      ! Output parameters.
      type(t_Escan), intent(out) :: escan_input
      ! Handle to mpi communicator.
      integer, intent(in) :: comm


      ! dummy argument.
      integer linenumber
      ! Input unit number
      integer, parameter :: unitnumber = 9
      ! special kpoint  requested?
      integer ikpt, rank, ierr, width
      common /comikpt/ikpt

      escan_input%ecp%comm_handle    = comm
      escan_input%ecp%fileescaninput = filename;
      call mpi_comm_rank(comm, rank, ierr)
      if( rank == 0 ) then
        open(unitnumber,file=escan_input%ecp%fileescaninput, status='OLD')        
          rewind( unitnumber )
          read(unitnumber,*) linenumber, escan_input%ecp%filepot
          read(unitnumber,*) linenumber, escan_input%ecp%filewg_out
          read(unitnumber,*) linenumber, escan_input%method
          read(unitnumber,*) linenumber, escan_input%ecp%Eref, &
                                         escan_input%Ecutoff,&
                                         escan_input%Gsmooth,&
                                         escan_input%KineticScaling
          read(unitnumber,*) linenumber, escan_input%ecp%mx
          read(unitnumber,*) ! Minimizer parameters,.
          read(unitnumber,*) ! input wavefunctions.
          read(unitnumber,*) ! input wavefunctions.
          read(unitnumber,*) ! input wavefunctions.
          read(unitnumber,*) ! output wavefunction parameters.
          read(unitnumber,*) linenumber, ikpt, &
                             escan_input%kpoint(1:3), escan_input%lattice_scale
          read(unitnumber,*) linenumber, escan_input%PotentialType
        close(unitnumber)

        ! some extra logic. Done now so that nodes are exactly the same.
        escan_input%is_gamma = .false.
        if( ikpt .eq. 0 ) escan_input%is_gamma = .true.
        
        ! Consistency checks and unit mending
        if( escan_input%is_gamma .eqv. .true. ) then
          escan_input%kpoint = 0.d0
        else
          escan_input%kpoint = escan_input%kpoint * 2.d0 * pi / escan_input%lattice_scale
        endif
        ! to atomic units from Rydbergs(?)
        escan_input%Ecutoff = escan_input%Ecutoff * 0.5d0
        ! to atomic units from eV.
        escan_input%ecp%Eref = escan_input%ecp%Eref / 27.211396d0
      endif

      width = len(escan_input%ecp%filepot)
      call mpi_bcast(escan_input%ecp%filepot, width, MPI_CHARACTER, 0, comm, ierr)
      call mpi_bcast(escan_input%ecp%filewg_out, width, MPI_CHARACTER, 0, comm, ierr)
      call mpi_bcast(escan_input%method, 1, MPI_INTEGER, 0, comm, ierr)
      call mpi_bcast(escan_input%ecp%Eref, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      call mpi_bcast(escan_input%Ecutoff, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      call mpi_bcast(escan_input%Gsmooth, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      call mpi_bcast(escan_input%KineticScaling, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      call mpi_bcast(ikpt, 1, MPI_INTEGER, 0, comm, ierr)
      call mpi_bcast(escan_input%kpoint(1), 3, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      call mpi_bcast(escan_input%lattice_scale, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      call mpi_bcast(escan_input%PotentialType, 1, MPI_INTEGER, 0, comm, ierr)

      escan_input%is_gamma = .false.
      if( ikpt .eq. 0 ) escan_input%is_gamma = .true.
      if( escan_input%PotentialType .lt. 0 .or. escan_input%PotentialType .gt. 3 ) then
        print *, "Potential Type is not acceptable: ", escan_input%PotentialType 
        call exit( 1 )
      endif
      escan_input%with_spinorbit = .false.
      if( escan_input%PotentialType .eq. 3 ) escan_input%with_spinorbit = .true.

    end subroutine ! read_escan_input

    ! Reads in cell-vectors and proc numbers.
    subroutine read_potentialinput(escan_input, lattice)
      ! Escan parameters. Some are used and some are set here.
      type(t_Escan), intent(inout) :: escan_input
      ! Lattice parameters
      type(t_Lattice), intent(inout) :: lattice

      integer, parameter :: inputunit = 9
      integer nnodes, inode, ierr

      call mpi_comm_size(escan_input%ecp%comm_handle, nnodes, ierr)
      call mpi_comm_rank(escan_input%ecp%comm_handle, inode, ierr)
      inode = inode + 1

      ! Read in lattice and mesh information from potential input.
      if (inode .eq. 1) then
         open(unit=inputunit,file=escan_input%ecp%filepot,form='unformatted',status='old')
           rewind(inputunit)
           read(inputunit) escan_input%mesh, ierr ! number of nodes when written
!          if( ierr .ne. nnodes ) then
!            write (*,*) ierr, nnodes 
!            stop "Number of nodes used to write potential different from when reading it." 
!          endif
           read(inputunit) lattice%rcell
         close(inputunit)
      end if

      call mpi_bcast(escan_input%mesh(1),3,MPI_INTEGER,0, escan_input%ecp%comm_handle,ierr)
      call mpi_bcast(lattice%rcell, 9, MPI_REAL8, 0, escan_input%ecp%comm_handle,ierr)

    end subroutine ! read_potential_input

    ! Sets common block parameters using input arguments.
    subroutine set_common_blocks(escan_input, lattice)
      use data, only: mg_nx, mr_n
      use load_data, only: nr1x, nr2x, nr3x
      ! Input parameters to set to common block
      type(t_Escan), intent(out) :: escan_input
      ! Lattice parameters to set to common block
      type(t_Lattice), intent(out) :: lattice
      ! Common blocks from escan
      real(kind=8) ::  AL(3,3), Ecut, vol, delta_k
      integer ikpt, inode, nnodes, n1, n2, n3, ng, ng_n, mx, nr, ierr, totg
      common /com123/n1,n2,n3,ng,ng_n,nr,mx,vol
      common /mpi_data/inode,nnodes
      common /comAD/AL,Ecut
      common /comikpt/ikpt
 
      call mpi_comm_rank( escan_input%ecp%comm_handle,inode,ierr)
      call mpi_comm_size( escan_input%ecp%comm_handle,nnodes,ierr)
      inode = inode + 1

      Ecut = escan_input%Ecutoff
      AL(:,:) = lattice%rcell(:,:)
      n1 = escan_input%mesh(1)
      n2 = escan_input%mesh(2)
      n3 = escan_input%mesh(3)
      nr1x=n1
      nr2x=n2
      nr3x=n3+2
      vol =   al(3,1) * ( al(1,2)*al(2,3) - al(1,3)*al(2,2) ) &
            + al(3,2) * ( al(1,3)*al(2,1) - al(1,1)*al(2,3) ) &
            + al(3,3) * ( al(1,1)*al(2,2) - al(1,2)*al(2,1) )
      vol=abs(vol)
      delta_k=(2*pi)**3/vol
      totg=(0.5d0*4.0d0/3.0d0*pi*(sqrt(2.0d0*Ecut))**3)/delta_k
      mg_nx=2*(int(1.1*totg/nnodes)+100)    ! now we have the whole sphere
      mx = escan_input%ecp%mx

      nr=n1*n2*n3
      mr_n=2*n1*n2*(n3+2)/nnodes    
 
    end subroutine set_common_blocks

    subroutine prepare_fft_and_allocate_arrays(params, lattice)

      use fft_data, only : fft_allocate, ncolx
      use load_data, only : load_allocate, ngtotnod
      use data, only : mg_nx, gkk_n, wg_n, inv_n

      type(t_Escan),   intent(in) :: params
      type(t_Lattice), intent(in) :: lattice

      ! Escan common group stuff.
      real(kind=8) ::  vol
      integer inode,nnodes, n1, n2, n3, ng, ng_n, mx, nr
      common /mpi_data/inode,nnodes
      common /com123/n1,n2,n3,ng,ng_n,nr,mx,vol

      logical :: with_spinorbit = .false.

      if( params%PotentialType .eq. 3 ) with_spinorbit = .true.

      call fft_allocate(n1,n2,n3,nnodes)
      call load_allocate(ncolx,nnodes,mg_nx)
      allocate(gkk_n(mg_nx))
      allocate(wg_n(mg_nx))
      allocate(inv_n(mg_nx))

      if( params%is_gamma .eqv. .true. ) then 
        call gen_G_compk0( lattice%rcell, lattice%kcell, params%Ecutoff, &
                           params%Gsmooth, params%KineticScaling )
      else
        call gen_G_compk( lattice%rcell, lattice%kcell, params%Ecutoff, &
                          params%Gsmooth, params%KineticScaling, &
                          params%kpoint(1), params%kpoint(2), params%kpoint(3) )
      endif

      ! some mpi distribution checking.
      if(ngtotnod(inode)>mg_nx) then
         write(6,*)'Not allocated enough g points on node',inode
         write(6,*)'ngtot per node = ',ngtotnod(inode)
         write(6,*)'Allocated no. = ',mg_nx
         call abort()
      end if

      ! Set up index matrices for fft routinines 

      deallocate(gkk_n)
      if(      params%is_gamma .eqv. .false. &
          .or. params%with_spinorbit .eqv. .false. ) deallocate(inv_n)

    end subroutine ! prepare_fft_and_allocate_arrays

    subroutine escan_cleanup
      use fft_data, only: fft_deallocate
      use load_data, only: load_deallocate
      use data, only: gkk_n, wg_n, inv_n

      call fft_deallocate
      call load_deallocate
      if( allocated(gkk_n) ) deallocate(gkk_n)
      if( allocated(wg_n) )  deallocate(wg_n)
      if( allocated(inv_n) ) deallocate(inv_n)

    end subroutine escan_cleanup

end Module Escan

! gets dimension for real space wavefunctions.
subroutine escan_get_nr(n)
  integer, intent(out) :: n
  integer n1,n2,n3,ng,ng_n,nr,mx
  real*8 vol
  common /com123/n1,n2,n3,ng,ng_n,nr,mx,vol
  n = n1*n2*n3
end subroutine
! gets dimension for real space wavefunctions.
subroutine escan_get_n1_n2_n3(a1, a2, a3)
  integer, intent(out) :: a1, a2, a3
  real*8 vol
  integer n1,n2,n3,ng,ng_n,nr,mx
  common /com123/n1,n2,n3,ng,ng_n,nr,mx,vol
  a1 = n1
  a2 = n2
  a3 = n3
end subroutine
! gets dimension for real space wavefunctions.
subroutine escan_get_mr_n(n)
  use data, only: mr_n
  integer, intent(out) :: n
  n = mr_n
end subroutine

! gets cell vectors
subroutine escan_get_cell(a0, a1, a2)
  real(kind=8),dimension(3), intent(inout) :: a0
  real(kind=8),dimension(3), intent(inout) :: a1
  real(kind=8),dimension(3), intent(inout) :: a2
  real(kind=8) ::  AL(3,3)
  common /comAD/AL

  a0 = AL(:,1)
  a1 = AL(:,2)
  a2 = AL(:,3)
end subroutine

