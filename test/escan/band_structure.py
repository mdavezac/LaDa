from numpy import array
from boost.mpi import world
from lada.escan import read_input, exec_input, BPoints, ReducedBPoints
from lada.crystal import nb_valence_states  
from lada.crystal.binary import zinc_blende

X = array( [1,0,0], dtype="float64" )
G = array( [0,0,0], dtype="float64" )
L = array( [0.5,0.5,0.5], dtype="float64" )
W = array( [0, 0.5,1], dtype="float64" )

input = read_input('input.py')
kescan = exec_input(repr(input.escan).replace('Escan', 'KEscan')).functional

structure = zinc_blende().to_structure(subs={'A':'Si', 'B':'Si'}) 
structure.scale = 5.45

kescan.fft_mesh = 14, 14, 14
kescan.kpoints = ReducedBPoints(density=20) + (X, G) + (G, L)
result = kescan( structure, comm=world, outdir='results/kescan', 
                 nbstates = nb_valence_states(structure) + 4,
                 eref = None )
