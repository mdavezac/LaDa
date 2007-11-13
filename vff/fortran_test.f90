!
!  Version: $Id$
!
program main

  real*8 :: pos(3,8)
  real*8 :: cell(3,3)
  real*8 :: alphas(5);
  real*8 :: betas(5);
  real*8 :: energy;
  character (len=2) :: At(8)

  cell(:,:) = 0.d0
  cell(1,1) = 1.d0
  cell(2,2) = 1.d0
  cell(3,3) = 1.d0

  call vff_create()

  pos(1,1) = 0.00; pos(2,1) = 0.00; pos(3,1) = 0.00; At(1) =  "Ga"
  pos(1,2) = 0.25; pos(2,2) = 0.25; pos(3,2) = 0.25; At(2) =  "As"
  pos(1,3) = 0.50; pos(2,3) = 0.50; pos(3,3) = 0.00; At(3) =  "In"
  pos(1,4) = 0.75; pos(2,4) = 0.75; pos(3,4) = 0.25; At(4) =  "As"
  pos(1,5) = 0.50; pos(2,5) = 0.00; pos(3,5) = 0.50; At(5) =  "In"
  pos(1,6) = 0.75; pos(2,6) = 0.25; pos(3,6) = 0.75; At(6) =  "Sb"
  pos(1,7) = 0.00; pos(2,7) = 0.50; pos(3,7) = 0.50; At(7) =  "Ga"
  pos(1,8) = 0.25; pos(2,8) = 0.75; pos(3,8) = 0.75; At(8) =  "Sb"

  call vff_cell( cell(1,1) )
  call vff_scale( 5.65329836585d0 )
  call vff_atoms( 8, pos, At)

  alphas(:) = 0.d0; alphas(1) = 35.180; call vff_bond( "In-As", 2.622d0, alphas )
  alphas(:) = 0.d0; alphas(1) = 41.190; call vff_bond( "Ga-As", 2.448d0, alphas )
  alphas(:) = 0.d0; alphas(1) = 33.160; call vff_bond( "Ga-Sb", 2.640d0, alphas )
  alphas(:) = 0.d0; alphas(1) = 26.610; call vff_bond( "In-Sb", 2.805d0, alphas )

  betas(:) = 0.d0;
  betas(1) = 8.93823;  call vff_angle( "As-Ga-As", -0.3333333d0, 0.d0, betas )
  betas(1) = 8.08355;  call vff_angle( "As-Ga-Sb", -0.3333333d0, 0.d0, betas )
  betas(1) = 7.2289;   call vff_angle( "Sb-Ga-Sb", -0.3333333d0, 0.d0, betas )
  betas(1) = 5.48808;  call vff_angle( "In-As-In", -0.3333333d0, 0.d0, betas )
  betas(1) = 7.213155; call vff_angle( "Ga-As-In", -0.3333333d0, 0.d0, betas )
  betas(1) = 8.93823;  call vff_angle( "Ga-As-Ga", -0.3333333d0, 0.d0, betas )
  betas(1) = 5.48808;  call vff_angle( "As-In-As", -0.3333333d0, 0.d0, betas )
  betas(1) = 4.88615;  call vff_angle( "Sb-In-As", -0.3333333d0, 0.d0, betas )
  betas(1) = 4.2842;   call vff_angle( "Sb-In-Sb", -0.3333333d0, 0.d0, betas )
  betas(1) = 7.2289;   call vff_angle( "Ga-Sb-Ga", -0.3333333d0, 0.d0, betas )
  betas(1) = 4.2842;   call vff_angle( "In-Sb-In", -0.3333333d0, 0.d0, betas )
  betas(1) = 5.75655;  call vff_angle( "Ga-Sb-In", -0.3333333d0, 0.d0, betas )

  call vff_minimize( energy );

  print *, "result= ", energy
  call vff_print_structure()
! call vff_print_lattice()


  call vff_destroy()

end program
