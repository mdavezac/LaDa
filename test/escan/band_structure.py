from numpy import array
from pylada.escan import read_input, exec_input, ReducedBPoints
from pylada.crystal.binary import zinc_blende

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
result = kescan( structure, outdir='results/kescan', 
                 nbstates = len(structure.atoms) * 4 + 4,
                 eref = None )
