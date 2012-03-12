from lada.crystal import read_poscar, make_surface
from numpy import array, all, abs
########################################################

#ZincBlende
structure = read_poscar(path='poscars/POSCAR_GaAs')
miller = array([ 1, 1, 2]) 
s=make_surface(structure=structure,miller=miller,nlayers=7,vacuum=15)

assert all(abs(s.cell-array([[ 7.19161544e-01, 1.29679552e-07, 1.80382592e-06],
                             [ 3.61568401e-24, 1.01704719e+00, 7.43418378e-06],
                             [ 5.23168124e-17, 1.45790583e-23, 7.68899277e+00]]))<1e-7)
#Wurtzite
structure = read_poscar(path='poscars/POSCAR_ZnO')
miller = array([ 2,-1, 0]) 
s=make_surface(structure=structure,miller=miller,nlayers=7,vacuum=15)

assert all(abs(s.cell-array([[ 1.63173500e+00,  8.30612023e-07, -1.20086204e-05],
                             [ 8.11412482e-23,  1.74969923e+00,  2.53831558e-05],
                             [ 5.29395592e-23, -1.30910525e-22,  1.16866921e+01]]))<1e-7)

#Rutile
structure = read_poscar(path='poscars/POSCAR_SnO2')
miller = array([ 1, 1, 0]) 
s=make_surface(structure=structure,miller=miller,nlayers=7,vacuum=15)

assert all(abs(s.cell-array([[ 3.24605173e+00,  1.64681605e-07,  2.07761230e-05],
                             [-4.23930064e-24,  6.82013915e+00,  1.71171464e-05],
                             [ 9.25408310e-24, -5.07751038e-16,  6.27410232e+01]]))<1e-7)

