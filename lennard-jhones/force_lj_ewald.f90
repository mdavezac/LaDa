subroutine force_lj_ewald (natom, a, ityp, tau, f, stress, ener )
  !  
  use ep_param    ! XZ

  implicit none
  
  !  Original LJ subroutine from Allen-Tildsley's book on
  !  Physics of Liquids
  
  INTEGER :: natom, iind, ij_typ    ! XZ
  
  REAL(kind=dbl) :: V, W 
  REAL(kind=dbl) :: sij(3), rij(3), a(3,3), axis_ewald(3,3)
  REAL(kind=dbl) :: tau(3,natom), f(3,natom)
  
  integer        :: ityp(natom)
  INTEGER        :: i,j, ii, jj, ispec, jspec
  INTEGER        :: icell1, icell2, icell3
  integer        :: ncell1, ncell2, ncell3
  
  REAL(kind=dbl) :: RCUTSQ, SIGSQ, EPS4, EPS24
  REAL(kind=dbl) :: FXI, FYI, FZI
  REAL(kind=dbl) :: radius_i, radius_j, sigma
  REAL(kind=dbl) :: r_ij, r_ij_inv, rsq_ij
  REAL(kind=dbl) :: rxij, ryij, rzij, fxij, fyij, fzij
  REAL(kind=dbl) :: SR2, SR6, SR12, V_LJ_IJ, vij, WIJ, FIJ, v_elec_ij
  real(kind=dbl) :: q_i, q_j, f_lj_ij, f_elec_ij, rsq_ij_inv
  
  real(kind=dbl) :: ener
  real(kind=dbl) :: alpha,alpha_inv, alpha_cube
  
  ! variables used to calculate the stress tensor
  
  REAL(kind=dbl) :: stress(3,3)
  REAL(kind=dbl) :: temp00,temp01,temp02
  REAL(kind=dbl) :: temp10,temp11,temp12
  REAL(kind=dbl) :: temp20,temp21,temp22

  integer :: IPR
  REAL(kind=dbl) :: EEWALD_ewald 
  REAL(kind=dbl) :: FEWA_ewald(3,natom) 
  REAL(kind=dbl) :: FEWAC_ewald(3,natom) 
  REAL(kind=dbl) :: sigma_ewald(6), ZZ_ewald(natom), RCUT_ewald
  INTEGER        :: MXDNAT_ewald
  REAL(kind=dbl) :: stress_ewald(3,3)

  REAL(kind=dbl) :: V_coul_IJ
  REAL(kind=dbl) :: rcut_const, rcut_sigma, vcut, AUE    ! XZ

  AUE = 27.211396d0            ! XZ
  rcut_const = rcut_const0     ! XZ
  alpha = 0.529177d0
  alpha_cube = alpha*alpha*alpha
  alpha_inv = 1.d0/0.529177d0

  stress(:,:) = 0.d0
  stress_ewald(:,:) = 0.d0
  
  ncell1 = 3
  ncell2 = 3
  ncell3 = 3
  
  V    = 0.d0
  V_coul_IJ = 0.d0
   
  ! Lennard-Jones part of energy, force, and stress tensor
  f(:,:) = 0.d0
  vcut = 0.d0

  do i = 1, natom
     do ii = 1, nspec_tot
        if (ityp(i) == id(ii)) then
           radius_i = rad_ion(ii)
           q_i = charge_ion(ii)
        end if
     end do
       
     do j = i+1, natom
        do jj = 1, nspec_tot
           if (ityp(j) == id(jj)) then
              radius_j = rad_ion(jj)
              q_j = charge_ion(jj)
           end if
        end do

        do iind = 1, nspec_tot*nspec_tot                                     ! XZ
           if (ityp(i) == ibond(1,iind) .and. ityp(j) == ibond(2,iind)) then
              ij_typ = iind
           end if
        end do
        
        sigma = (radius_i + radius_j) * rsigma(ij_typ)                       ! XZ
        sigsq = sigma * sigma

        do icell1 = -ncell1,ncell1
           do icell2 = -ncell2,ncell2
              do icell3 = -ncell3,ncell3
                 
                 sij(1) = tau(1,i) - tau(1,j)                   
                 sij(2) = tau(2,i) - tau(2,j)                   
                 sij(3) = tau(3,i) - tau(3,j) 
                 
                 sij(1) = sij(1) - anint (sij(1))
                 sij(2) = sij(2) - anint (sij(2))
                 sij(3) = sij(3) - anint (sij(3))
                 
                 sij(1) = sij(1) + dfloat(icell1)
                 sij(2) = sij(2) + dfloat(icell2)
                 sij(3) = sij(3) + dfloat(icell3)   
                 
                 rij = matmul(a,sij)
                 rsq_ij = sum(rij(:)*rij(:))
                 r_ij = dsqrt(rsq_ij)
                 
                 rxij = rij(1)
                 ryij = rij(2)
                 rzij = rij(3)
                 
                 rcut_sigma = rcut_const * sigma
                 if (r_ij < rcut_sigma ) then

                 ! -----------------------------
                 r_ij_inv = 1.d0/rcut_sigma
                 rsq_ij_inv = r_ij_inv * r_ij_inv
                  
                 sr2   = sigsq * rsq_ij_inv
                 sr6   = sr2 * sr2 * sr2
                 sr12  = sr6 * sr6
                 vij   = sr12 - sr6 * (1.d0-PEGS)             ! XZ

                 VIJ   = 4.d0*epslon(ij_typ) * VIJ            ! XZ

                 vcut = vcut - vij 

                 ! -----------------------------
                 r_ij_inv = 1.d0/r_ij
                 rsq_ij_inv = r_ij_inv * r_ij_inv
                 
                 sr2   = sigsq * rsq_ij_inv
                 sr6   = sr2 * sr2 * sr2
                 sr12  = sr6 * sr6
                 
                 vij   = sr12 - sr6 * (1.d0-PEGS)             ! XZ
                 WIJ   = vij + sr12
                 VIJ   = 4.d0*epslon(ij_typ) * VIJ            ! XZ
                 WIJ   = 4.d0*epslon(ij_typ) * 6.d0 * WIJ     ! XZ

                 V_LJ_IJ = V_LJ_IJ + VIJ

                 !V_coul_IJ = V_coul_IJ + v_elec_ij
                 !V  = V + VIJ + v_elec_ij

                 v  = v + vij
                 f_lj_ij =  -wij*r_ij_inv

                 ! f_elec_ij =   -v_elec_ij*(r_ij_inv + mu)
                 ! fij = f_lj_ij + f_elec_ij

                 ! contribution to forces from the pair potential

                 fij = f_lj_ij
                 
                 FXIJ  =  FIJ * RXIJ * r_ij_inv
                 FYIJ  =  FIJ * RYIJ * r_ij_inv
                 FZIJ  =  FIJ * RZIJ * r_ij_inv
                 
                 f(1,i) = f(1,i) - FXIJ
                 f(2,i) = f(2,i) - FYIJ
                 f(3,i) = f(3,i) - FZIJ

                 f(1,j) = f(1,j) + FXIJ
                 f(2,j) = f(2,j) + FYIJ
                 f(3,j) = f(3,j) + FZIJ
                 
                 !FEX(I) = FEX(I) - f_elec_ij * RXIJ * r_ij_inv
                 !FEY(I) = FEY(I) - f_elec_ij * RYIJ * r_ij_inv
                 !FEZ(I) = FEZ(I) - f_elec_ij * RZIJ * r_ij_inv

                 !FEX(J) = FEX(J) + f_elec_ij * RXIJ * r_ij_inv
                 !FEY(J) = FEY(J) + f_elec_ij * RYIJ * r_ij_inv
                 !FEZ(J) = FEZ(J) + f_elec_ij * RZIJ * r_ij_inv

                 ! stress tensor : components of the virial for pair potential

                 temp00 = -FXIJ*RXIJ 
                 temp11 = -FYIJ*RYIJ 
                 temp22 = -FZIJ*RZIJ 
                 
                 temp01 = -FXIJ*RYIJ 
                 temp02 = -FXIJ*RZIJ 
                 
                 temp10 = -FYIJ*RXIJ 
                 temp12 = -FYIJ*RZIJ 
                
                 temp20 = -FZIJ*RXIJ 
                 temp21 = -FZIJ*RYIJ
                 
                 stress(1,1) = stress(1,1) + temp00
                 stress(2,2) = stress(2,2) + temp11
                 stress(3,3) = stress(3,3) + temp22
                 
                 stress(1,2) = stress(1,2) + temp01
                 stress(1,3) = stress(1,3) + temp02
                 
                 stress(2,1) = stress(2,1) + temp10
                 stress(2,3) = stress(2,3) + temp12
                 
                 stress(3,1) = stress(3,1) + temp20
                 stress(3,2) = stress(3,2) + temp21

                 end if
                 
              end do  !  icell3
           end do  !  icell2
        end do  !  icell1
         
     end do  ! atom j
      
  end do  ! atom i
   
  do i = 1, natom
     do ii = 1, nspec_tot
        if (ityp(i) == id(ii)) then
           zz_ewald(i) = charge_ion(ii)
        end if
     end do
  end do

  RCUT_ewald = RCUT_ewaldA * alpha_inv    ! XZ
  
  axis_ewald = a * alpha_inv
  !axis_ewald = a 

  call EWALDF(0,EEWALD_ewald,FEWA_ewald,FEWAC_ewald,sigma_ewald, &
      natom,tau,ZZ_ewald,RCUT_ewald,axis_ewald,natom)

  FEWAC_ewald = FEWAC_ewald*0.5d0*alpha_inv

  ! sigma_ewald =  sigma_ewald * 0.5d0 * (alpha**3.d0)

  sigma_ewald =  sigma_ewald * 0.5d0

  ! write(6,'(2(f20.10,3x))') v , vcut                      ! XZ
  ! write(6,'(3(f10.8,3x))') f
  ! write(6,'(3(f10.8,3x))') stress

  stress_ewald(1,1) =  sigma_ewald(1) 
  stress_ewald(2,2) =  sigma_ewald(2)
  stress_ewald(3,3) =  sigma_ewald(3)

  stress_ewald(1,2) =  sigma_ewald(4)
  stress_ewald(2,1) =  sigma_ewald(4)

  stress_ewald(2,3) =  sigma_ewald(5)
  stress_ewald(3,2) =  sigma_ewald(5)

  stress_ewald(1,3) =  sigma_ewald(6)
  stress_ewald(3,1) =  sigma_ewald(6)

  stress = stress + stress_ewald * AUE                      ! XZ

  f = f + FEWAC_ewald * AUE                                 ! XZ
  
  ener = v + vcut + EEWALD_ewald* 0.5d0 * AUE               ! XZ
 
end subroutine force_lj_ewald
