!
! Version: $Id$
!

module MomentumDipole

  use escan_comp_api

  implicit none
  include 'mpif.h'

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

  end type ! t_Escan

  type t_Lattice

    ! real-space cell vectors
    real*8 :: rcell(3,3)
    ! reciprocal-space cell vectors.
    real*8 :: kcell(3,3)

  endtype

  ! Escan parameters.
  type( t_Escan ), public :: params 
  ! Lattice cell vectors.
  type( t_Lattice ), public :: lattice 

  private
  ! Pi parameter
  real( kind=8 ), parameter, private :: pi = 4.d0*datan(1.d0)

  public :: momentum

  contains

    ! reads escan input necessary for momentum dipole computations
    subroutine read_escaninput( in_filename, in_fnamelength )

      ! input filename length
      integer, intent(in) :: in_fnamelength
      ! input filename
      character(len=in_fnamelength), intent(in) :: in_filename


      ! dummy argument.
      integer linenumber
      ! Input unit number
      integer, parameter :: unitnumber = 9
      ! Common blocks from escan
      real(kind=8) ::  AL(3,3) 
      real(kind=8) ::  Ecut 
      integer ikpt
      common /comAD/AL,Ecut
      common /comikpt/ikpt

      params%ecp%fileescaninput = in_filename;
      
      ! all processors perform reading. filename may (should) be different.
      open(unitnumber,file=params%ecp%fileescaninput, status='OLD')        
        rewind( unitnumber )
        read(unitnumber,*) linenumber, params%ecp%filepot
        read(unitnumber,*) linenumber, params%ecp%filewg_out
        read(unitnumber,*) linenumber, params%method
        read(unitnumber,*) linenumber, params%ecp%Eref, &
                                       params%Ecutoff,&
                                       params%Gsmooth,&
                                       params%KineticScaling
        read(unitnumber,*) ! linenumber, params%ecp%mx
        read(unitnumber,*) ! Minimizer parameters,.
        read(unitnumber,*) ! input wavefunctions.
        read(unitnumber,*) ! input wavefunctions.
        read(unitnumber,*) ! input wavefunctions.
        read(unitnumber,*) ! output wavefunction parameters.
        read(unitnumber,*) linenumber, ikpt, &
                           params%kpoint(1:3), params%lattice_scale
        read(unitnumber,*) linenumber, params%PotentialType
      close(unitnumber)

      params%is_gamma = .false.
      if( ikpt .eq. 0 ) params%is_gamma = .true.

      ! Consistency checks and unit mending
      if( params%is_gamma .eqv. .true. ) then
        params%kpoint = 0.d0
      else
        params%kpoint = params%kpoint * 2.d0 * pi / params%lattice_scale
      endif
      if( params%PotentialType .lt. 0 .or. params%PotentialType .gt. 3 ) then
        print *, "Potential Type is not acceptable: ", params%PotentialType 
        call exit( 1 )
      endif
      params%with_spinorbit = .false.
      if( params%PotentialType .eq. 3 ) params%with_spinorbit = .true.
      ! to atomic units from Rydbergs(?)
      params%Ecutoff = params%Ecutoff * 0.5d0
      Ecut = params%Ecutoff
      ! to atomic units from eV.
      params%ecp%Eref = params%ecp%Eref / 27.211396d0


    end subroutine ! read_escan_input

    ! Reads in cell-vectors and proc numbers.
    subroutine read_potentialinput()
    
      use data, only : mg_nx

      ! input unit number
      integer, parameter :: inputunit = 9
      ! Escan common group stuff.
      real(kind=8) ::  AL(3,3) 
      real(kind=8) ::  vol,Ecut 
      integer inode,nnodes, n1, n2, n3, ng, ng_n, mx, nr, ierr
      common /mpi_data/inode,nnodes
      common /com123/n1,n2,n3,ng,ng_n,nr,mx,vol
      common /comAD/AL,Ecut
      ! mpi error integer
      real( kind = 8 ) delta_k, totg

      ! Read in lattice and mesh information from potential input.
      if (inode==1) then
         open(unit=inputunit,file=params%ecp%filepot,form='unformatted',status='old')
         rewind(inputunit)
         read(inputunit) n1,n2,n3 ! number of nodes when written
         read(inputunit) lattice%rcell
      end if

      call mpi_bcast(n1,1,MPI_INTEGER,0, params%ecp%comm_handle,ierr)
      call mpi_bcast(n2,1,MPI_INTEGER,0, params%ecp%comm_handle,ierr)
      call mpi_bcast(n3,1,MPI_INTEGER,0, params%ecp%comm_handle,ierr)
      call mpi_bcast(lattice%rcell, 9, MPI_REAL8, 0, params%ecp%comm_handle,ierr)

      ! escan common block stuff
      AL = lattice%rcell

      vol =   al(3,1) * ( al(1,2)*al(2,3) - al(1,3)*al(2,2) ) &
            + al(3,2) * ( al(1,3)*al(2,1) - al(1,1)*al(2,3) ) &
            + al(3,3) * ( al(1,1)*al(2,2) - al(1,2)*al(2,1) )
      vol=abs(vol)
      delta_k=(2*pi)**3/vol
      totg=(0.5d0*4.0d0/3.0d0*pi*(sqrt(2.0d0*Ecut))**3)/delta_k
      mg_nx=2*(int(1.1*totg/nnodes)+100)    ! now we have the whole sphere

    end subroutine ! read_potential_input

    subroutine prepare_fft_and_allocate_arrays

      use fft_data, only : fft_allocate, ncolx
      use load_data, only : load_allocate, ngtotnod
      use data, only : mg_nx, gkk_n, wg_n, inv_n

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
      call fftprep_comp(n1,n2,n3)

      deallocate(gkk_n)
      if(      params%is_gamma .eqv. .false. &
          .or. params%with_spinorbit .eqv. .false. ) deallocate(inv_n)

    end subroutine ! prepare_fft_and_allocate_arrays

    subroutine deallocate_data

      use fft_data, only : fft_deallocate
      use load_data, only : load_deallocate
      use data, only : wg_n, inv_n

      call load_deallocate()
      call fft_deallocate()
      deallocate(wg_n)
      if( params%is_gamma .and. params%with_spinorbit ) deallocate(inv_n)

    end subroutine ! deallicate_data

    subroutine read_wavefunctions( io_wfns, in_indices )

      use data, only : mg_nx

      ! Input wavefunctions. First index is gspace coord, then band coord, then
      ! spin index.
      complex(kind=8), dimension(:,:,:), intent(inout) :: io_wfns
      ! Indices of wavefunctions when written to disk.
      integer, dimension(:), intent(in) :: in_indices

      real(kind=8) ::  vol
      integer n1, n2, n3, ng, ng_n, mx, nr
      common /com123/n1,n2,n3,ng,ng_n,nr,mx,vol

      if( params%with_spinorbit ) then
        call read_wg_comp( io_wfns(:,:,1), io_wfns(:,:,2),n1,n2,n3,size(in_indices),mg_nx, &
                           params%Ecutoff, lattice%rcell, &
                           params%ecp%filewg_in, size( in_indices ), &
                           in_indices(1), params%with_spinorbit )
      else
        call read_wg_comp( io_wfns(:,:,1), io_wfns(:,:,1),n1,n2,n3,size(in_indices),mg_nx, &
                           params%Ecutoff, lattice%rcell, &
                           params%ecp%filewg_in, size( in_indices ), &
                           in_indices, params%with_spinorbit )
      endif

    end subroutine ! read_wavefunctions

    ! Performs actual momentum computation.
    subroutine compute_momentum( in_Awfns, in_Bwfns, io_dipoles, in_dip2, in_dip3 )

      use load_data, only : ngtotnod, n1p_n, n2p_n, n3p_n
      use data, only : wg_n

      ! Input wavefunctions ordered as in read_wavefunctions above.
      complex(kind=8), dimension(:,:,:), intent(in) :: in_Awfns
      ! Input wavefunctions ordered as in read_wavefunctions above.
      complex(kind=8), dimension(:,:,:), intent(in) :: in_Bwfns
      ! Dipoles with indices as follows:
      !   . 3 cartesian coordinates
      !   . first band index.
      !   . second band index
      ! Note that the system is either spin-degenerate ( eg
      ! params%with_spinorbit is false ) or the spin is not a good quantum
      ! number.
      ! When using spinors (params%with_spinorbit=true) and at gamma
      ! (params%is_gamma = true), computations are done using Kramer degeneracy.
      ! Hence the number of dipoles is twice the number of wfns in the arrays
      ! in_Awfns and in_Bwfnsm,  something which io_dipoles must reflect.
      ! See ortho_comp.fpp for relationship between calculated and implied
      ! wavefunctions.
      integer, intent(in):: in_dip2
      integer, intent(in):: in_dip3
      complex(kind=8), intent(inout) :: io_dipoles( 3, in_dip2, in_dip3 )

      real(kind=8) ::  vol
      integer n1, n2, n3, ng, ng_n, mx, nr
      common /com123/n1,n2,n3,ng,ng_n,nr,mx,vol

      integer inode,nnodes
      common /mpi_data/inode,nnodes

      ! number of band in A
      integer :: nb_Aband
      ! number of bands in B 
      integer :: nb_Bband

      real*8, dimension(:,:), allocatable :: gpoints
      integer*8, dimension(3) :: int_gpoint
      integer ::  nb_gpoints
      integer ierr, Aband, Bband, ig
      complex( kind=8 ) :: mpiresult
      complex(kind=8) realnb(5), a(5)

      nb_gpoints = ngtotnod( inode )
      nb_Aband = size( in_Awfns, 2 )
      nb_Bband = size( in_Bwfns, 2 )

      allocate( gpoints( nb_gpoints, 3 ) )
      ! Computes all gpoints.
      do ig = 1, nb_gpoints
        ! compute G point.
        int_gpoint(1) = n1p_n(ig) - 1
        int_gpoint(2) = n2p_n(ig) - 1
        int_gpoint(3) = n3p_n(ig) - 1

        if( int_gpoint(1) .gt. n1 / 2 ) int_gpoint(1) = int_gpoint(1) - n1
        if( int_gpoint(2) .gt. n2 / 2 ) int_gpoint(2) = int_gpoint(2) - n2
        if( int_gpoint(3) .gt. n3 / 2 ) int_gpoint(3) = int_gpoint(3) - n3

        gpoints(ig, 1) = 2.d0 * pi * sum( lattice%kcell(1,:) * int_gpoint ) &
                         * wg_n(ig) * wg_n(ig)
        gpoints(ig, 2) = 2.d0 * pi * sum( lattice%kcell(2,:) * int_gpoint ) &
                         * wg_n(ig) * wg_n(ig)
        gpoints(ig, 3) = 2.d0 * pi * sum( lattice%kcell(3,:) * int_gpoint ) &
                         * wg_n(ig) * wg_n(ig)
      enddo ! ig

      io_dipoles = 0d0
      ! now perform sums.
      do Aband = 1, nb_Aband
        do Bband = 1, nb_Bband

          ! With spin orbit and Kramer doubling.
          if(       params%is_gamma .eqv. .true. &
              .and. params%with_spinorbit .eqv. .true. ) then

            ! Kramer doubling means that half the states are computed
            ! implicitly using invariance through time reversal.
            ! We must still compute the dipole moments of these states.
            call sp_sum( in_Awfns(1:nb_gpoints, Aband, 1 ), & 
                         in_Bwfns(1:nb_gpoints, Bband, 1 ), & 
                         in_Awfns(1:nb_gpoints, Aband, 2 ), & 
                         in_Bwfns(1:nb_gpoints, Bband, 2 ), & 
                         io_dipoles(:, 2*(Aband-1)+1, 2*(Bband-1)+1 ) )
            call invA_sum( in_Awfns(1:nb_gpoints, Aband, 2 ), & 
                           in_Bwfns(1:nb_gpoints, Bband, 1 ), & 
                           in_Awfns(1:nb_gpoints, Aband, 1 ), & 
                           in_Bwfns(1:nb_gpoints, Bband, 2 ), & 
                           io_dipoles(:, 2*Aband, 2*(Bband-1)+1 ) )
            call invB_sum( in_Awfns(1:nb_gpoints, Aband, 1 ), & 
                           in_Bwfns(1:nb_gpoints, Bband, 2 ), & 
                           in_Awfns(1:nb_gpoints, Aband, 2 ), & 
                           in_Bwfns(1:nb_gpoints, Bband, 1 ), & 
                           io_dipoles(:, 2*(Aband-1)+1, 2*Bband ) )
            call invAinvB_sum( in_Awfns(1:nb_gpoints, Aband, 2 ), & 
                               in_Bwfns(1:nb_gpoints, Bband, 2 ), & 
                               in_Awfns(1:nb_gpoints, Aband, 1 ), & 
                               in_Bwfns(1:nb_gpoints, Bband, 1 ), & 
                               io_dipoles(:, 2*Aband, 2*Bband ) )

          ! With spin orbit and without Kramer doubling
          else if(       params%is_gamma .eqv. .false. &
                   .and. params%with_spinorbit .eqv. .true. ) then
           
            call sp_sum( in_Awfns(1:nb_gpoints, Aband, 1 ), & 
                         in_Bwfns(1:nb_gpoints, Bband, 1 ), & 
                         in_Awfns(1:nb_gpoints, Aband, 2 ), & 
                         in_Bwfns(1:nb_gpoints, Bband, 2 ), & 
                         io_dipoles(:, Aband, Bband ) )
          ! Without spin orbit and without Kramer doubling (eg non-spin
          ! polarized)
          else 

            call dipole_sum( in_Awfns(1:nb_gpoints, Aband, 1 ), & 
                             in_Bwfns(1:nb_gpoints, Bband, 1 ), & 
                             io_dipoles(:, Aband, Bband ) )
          endif
        enddo ! Bband
      enddo ! Aband
     
      deallocate( gpoints )

      ! normalize to volume.
      io_dipoles = io_dipoles * vol

      if( params%is_gamma .eqv. .true. .and. params%with_spinorbit .eqv. .true. ) then
        call mpi_allreduce( MPI_IN_PLACE, io_dipoles(1,1,1), &
                            nb_Aband * nb_Bband * 4 * 3 * 2 , &
                            MPI_DOUBLE_PRECISION, MPI_SUM, params%ecp%comm_handle, ierr )
      else
        call mpi_allreduce( MPI_IN_PLACE, io_dipoles(1,1,1), &
                            nb_Aband * nb_Bband * 3 * 2 , &
                            MPI_DOUBLE_PRECISION, MPI_SUM, params%ecp%comm_handle, ierr )
      endif

      ! success
      if( ierr .ne. MPI_SUCCESS ) then
        print *, "Error encountered while calling mpi_reduce.", ierr
        call mpierror( ierr )
        call abort()
      endif

      contains 
        subroutine sp_sum( in_wfnA_up, in_wfnB_up, in_wfnA_dw, in_wfnB_dw, io_dipole )
 
          complex( kind=8 ), intent(in) :: in_wfnA_up( nb_gpoints )
          complex( kind=8 ), intent(in) :: in_wfnB_up( nb_gpoints )
          complex( kind=8 ), intent(in) :: in_wfnA_dw( nb_gpoints )
          complex( kind=8 ), intent(in) :: in_wfnB_dw( nb_gpoints )
 
          complex( kind=8 ), intent(inout) :: io_dipole(3)
 
          integer k
          
          do k = 1, 3
            io_dipole(k) = sum &
                           (   &
                             gpoints(:,k) *  &
                             (    &
                                 conjg( in_wfnA_up )  * in_wfnB_up &
                               + conjg( in_wfnA_dw )  * in_wfnB_dw &
                             ) &
                           )  
          enddo
          print *, "sp_sum ", io_dipole * vol
        
        end subroutine 

        ! wfnA are the time reversed wavefunction of the computed wfn.
        ! sign of up = -down component comes from using spinors.
        subroutine invA_sum( in_wfnA_up, in_wfnB_up, in_wfnA_dw, in_wfnB_dw, io_dipole )
 
          use data, only : inv_n
 
          complex( kind=8 ), intent(in) :: in_wfnA_up( nb_gpoints )
          complex( kind=8 ), intent(in) :: in_wfnB_up( nb_gpoints )
          complex( kind=8 ), intent(in) :: in_wfnA_dw( nb_gpoints )
          complex( kind=8 ), intent(in) :: in_wfnB_dw( nb_gpoints )
 
          complex( kind=8 ), intent(inout) :: io_dipole(3)
 
          integer k
          
          do ig = 1, nb_gpoints
            do k = 1, 3
              io_dipole(k) = io_dipole(k) + gpoints(ig, k) &
                             * (   &
                                 - in_wfnA_up( inv_n(ig) ) * in_wfnB_up(ig) &
                                 + in_wfnA_dw( inv_n(ig) ) * in_wfnB_dw(ig) &
                               )
            enddo
          enddo
        
        end subroutine 

        ! wfnB are the time reversed wavefunction of the computed wfn.
        ! sign of up = -down component comes from using spinors.
        subroutine invB_sum( in_wfnA_up, in_wfnB_up, in_wfnA_dw, in_wfnB_dw, io_dipole )
 
          use data, only : inv_n
 
          complex( kind=8 ), intent(in) :: in_wfnA_up( nb_gpoints )
          complex( kind=8 ), intent(in) :: in_wfnB_up( nb_gpoints )
          complex( kind=8 ), intent(in) :: in_wfnA_dw( nb_gpoints )
          complex( kind=8 ), intent(in) :: in_wfnB_dw( nb_gpoints )
 
          complex( kind=8 ), intent(inout) :: io_dipole(3)
 
          integer k
          
          do ig = 1, nb_gpoints
            do k = 1, 3
              io_dipole(k) = io_dipole(k) + gpoints(ig, k) &
                             * (   &
                                 - conjg( in_wfnA_up(ig) ) * conjg( in_wfnB_up(inv_n(ig)) ) &
                                 + conjg( in_wfnA_dw(ig) ) * conjg( in_wfnB_dw(inv_n(ig)) ) &
                               )
            enddo
          enddo
        
        end subroutine 

        ! wfnA are the time reversed wavefunction of the computed wfn.
        ! wfnB are the time reversed wavefunction of the computed wfn.
        ! sign of up = -down component comes from using spinors.
        subroutine invAinvB_sum( in_wfnA_up, in_wfnB_up, in_wfnA_dw, in_wfnB_dw, io_dipole )
 
          use data, only : inv_n
 
          complex( kind=8 ), intent(in) :: in_wfnA_up( nb_gpoints )
          complex( kind=8 ), intent(in) :: in_wfnB_up( nb_gpoints )
          complex( kind=8 ), intent(in) :: in_wfnA_dw( nb_gpoints )
          complex( kind=8 ), intent(in) :: in_wfnB_dw( nb_gpoints )
 
          complex( kind=8 ), intent(inout) :: io_dipole(3)
 
          integer k
          
          do ig = 1, nb_gpoints
            do k = 1, 3
              io_dipole(k) = io_dipole(k) + gpoints(ig, k) &
                             * (   &
                                   in_wfnA_up(inv_n(ig)) * conjg( in_wfnB_up(inv_n(ig)) )&
                                 + in_wfnA_dw(inv_n(ig)) * conjg( in_wfnB_dw(inv_n(ig)) )&
                               )
            enddo
          enddo
        
        end subroutine 

        subroutine dipole_sum( in_wfnA, in_wfnB, io_dipole )
 
          complex( kind=8 ), intent(in) :: in_wfnA( nb_gpoints )
          complex( kind=8 ), intent(in) :: in_wfnB( nb_gpoints )
 
          complex( kind=8 ), intent(inout) :: io_dipole(3)
 
          integer k
          
          do k = 1, 3
            io_dipole(k) = sum( gpoints(:,k) * conjg( in_wfnA )  * in_wfnB )
          enddo
        
        end subroutine 

    end subroutine ! compute_dipole_moments.


    ! Computes momentum from complete diagonalization results (eg single wfn
    ! file ).
    subroutine momentum( in_inputfilename, in_fsize, &
                         in_dirvalence, in_dsizeA, &
                         in_dirconduction, in_dsizeB, &
                         in_indicesA, in_indicesB, & 
                         in_bandsA, in_bandsB, &
                         io_dipoles, in_dip2, in_dip3 )

      use data, only : mg_nx

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

      ! In case of spin orbit, an eigen wavefunction is a linear combination of
      ! Sz = up and Sz= down eigenstates. Sz (spin projected along z) itself is
      ! NOT a good quantum number. If there is no spin orbit, calculations are
      ! no spin polarized.
      ! Furthermore, if params%is_gamma is true, then Krammer degeneracy is
      ! used (see pescan's othornomalization procedure). With a particular set
      ! of input parameters, params%is_gamma may be false even-though
      ! params%kpoint is gamma. 
      ! The wave-functions below are arrays along G coords, then band indices,
      ! then Sz. The linear combination coefs are included within the
      ! wave-functions.
      complex( kind=8 ), dimension( :,:,: ), allocatable :: valence_wfns
      complex( kind=8 ), dimension( :,:,: ), allocatable :: conduction_wfns
      character, pointer :: cwd(:);
      integer cwd_size, dir_status, nb_spins

      integer inode, nnodes, ierr
      common /mpi_data/inode,nnodes

      real(kind=8) ::  vol
      integer n1, n2, n3, ng, ng_n, mx, nr
      common /com123/n1,n2,n3,ng,ng_n,nr,mx,vol

      call mpi_comm_rank( params%ecp%comm_handle,inode,ierr)
      call mpi_comm_size( params%ecp%comm_handle,nnodes,ierr)
      inode = inode+1

      ! These two subroutines are external "C++" routines. 
      cwd_size = 0
      call get_current_directory( cwd(1), cwd_size, dir_status ) 
      allocate( cwd( cwd_size +1 ) )
      cwd = ' '
      call get_current_directory( cwd(1), cwd_size, dir_status ) 

      
      call change_current_directory( in_dirvalence, in_dsizeA, dir_status )
      call read_escaninput( in_inputfilename, in_fsize )
      call read_potentialinput()

      nb_spins = 1
      if( params%with_spinorbit .eqv. .true. ) nb_spins = 2;

      allocate( valence_wfns(mg_nx, in_bandsA, nb_spins ) )
      allocate( conduction_wfns(mg_nx, in_bandsB, nb_spins ) )

      call prepare_fft_and_allocate_arrays()

      ! reads valence band wfns.
      params%ecp%filewg_in = params%ecp%filewg_out
      call read_wavefunctions( valence_wfns, in_indicesA )

      ! reads conduction input and wavefunctions
      call change_current_directory( cwd, cwd_size, dir_status )
      call change_current_directory( in_dirconduction, in_dsizeB, dir_status )
      call read_escaninput( in_inputfilename, in_fsize )
      params%ecp%filewg_in = params%ecp%filewg_out
      call read_wavefunctions( conduction_wfns, in_indicesB )


      call compute_momentum( valence_wfns, &
                             conduction_wfns, &
                             io_dipoles, in_dip2, in_dip3 )

      call deallocate_data()
      deallocate(valence_wfns)
      deallocate(conduction_wfns)
      call change_current_directory( cwd, cwd_size, dir_status )
      deallocate( cwd )

    end subroutine ! momentum

end module ! Escan
