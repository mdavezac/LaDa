program CLJ
!
! read ep.input & POSCAR_0
! units  Angstrom & eV

  use ep_param
  implicit none

  interface
     subroutine force_lj_ewald (natom, a, ityp, tau, f, stress, ener )
       use ep_param, only: dbl
       INTEGER, intent(in) :: natom
       ! Cell parameters.
       REAL(kind=dbl), intent(in) :: a(3,3)
       ! Occupation of each atomic site.
       integer, intent(in) :: ityp(natom)
       ! Coordinates of each atomic site.
       REAL(kind=dbl), intent(in) :: tau(3,natom)
       ! Forces on each atom.
       REAL(kind=dbl), intent(out) :: f(3,natom)
       ! Stress.
       REAL(kind=dbl), intent(out) :: stress(3,3)
       ! Energy
       real(kind=dbl), intent(out) :: ener
     end subroutine
  end interface



  integer :: i, j, k, natom, iind
  integer, allocatable :: mspecx(:), ityp(:)

  character(len=256) :: cdum

  real(kind=dbl) :: alat_ep, axis_ep(3,3), axis(3,3), stress(3,3)
  real(kind=dbl), allocatable :: tau(:,:), f(:,:) 
  real(kind=dbl) :: ener
  integer(kind=4) :: relaxer_handle;

  ! please search  calculation & relaxation

  ! input

  open(unit=35, file= "ep.input", status="old", form="formatted")

  read(35,*) nspec_tot                             !  number of atomic species

  allocate( id(nspec_tot), rad_ion(nspec_tot), charge_ion(nspec_tot), &
            mspecx(nspec_tot), ibond(2,nspec_tot*nspec_tot), &
            epslon(nspec_tot*nspec_tot), rsigma(nspec_tot*nspec_tot) )

  do i = 1, nspec_tot
     read(35,*) id(i), rad_ion(i), charge_ion(i)   !  id, rad_ion, charge_ion
  end do

  do iind = 1, nspec_tot*nspec_tot
     read(35,*) ibond(1:2,iind), epslon(iind), rsigma(iind)
  end do

  read(35,*) rcut_const0
  read(35,*) RCUT_ewaldA
 
  read(35,*) PEGS                                  !  PEGS   0.0-not  1.0-yes

  close(35) 


  open(unit=37, file = "POSCAR_0", status = "old", form = "formatted")

  read(37,*) cdum

  read(37,*) alat_ep
  read(37,*) axis_ep(:,1)
  read(37,*) axis_ep(:,2)
  read(37,*) axis_ep(:,3)

  read(37,*) mspecx(:)
  read(37,*) cdum

  natom = sum (mspecx(:))

  allocate( tau(3,natom), ityp(natom), f(3,natom) )

  do i = 1, natom
     read(37,*) tau(:,i)
  end do

  close(37)

  axis = alat_ep * axis_ep

  k = 0
  do i = 1, nspec_tot
     do j = 1, mspecx(i)
        k = k + 1
        ityp(k) = id(i)
     end do
  end do


  ! calculation & relaxation
  ! axis(3,3), ityp(natom), tau(3,natom), f(3,natom), stress(3,3)

  call create_relaxer( relaxer_handle, 8, "ep.input", force_lj_ewald )

  ener = 0e0
! call force_lj_ewald (natom, axis, ityp, tau, f, stress, ener )
  call call_relaxer ( relaxer_handle, natom, axis, ityp, tau, f, stress, ener )


  call release_relaxer( relaxer_handle )

  
  ! output

! open(unit=38, file = "CONTCAR", status = "unknown", form = "formatted")

  write(*,*) ener
  write(*,'("Crystal structure of individual")')

  write(*,'(f10.5)') 1.d0
  write(*,'(3(f15.5,2x))') axis(:,1)
  write(*,'(3(f15.5,2x))') axis(:,2)
  write(*,'(3(f15.5,2x))') axis(:,3)

  write(*,*) (mspecx(i), i=1, nspec_tot)
  write(*,'(a9)') 'Direct   '

  do i = 1, natom
     write(*,'(3(f15.5,2x))') tau(:,i)
  end do

! close(38)

! open(unit=39, file = "OSZICAR", status = "unknown", form = "formatted")

! write(39,'("   ep-energy")')
! write(39,'("   1 F= ", ES14.8, " d E= 0.000000E+00")') ener

! close(39)

end program CLJ
