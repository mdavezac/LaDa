""" Input script for the calculation of Point Defects. """
from lada.crystal import A2BX4

lattice = A2BX4.b5()
""" Back-bone lattice. """
# changes species in lattice.
for site in lattice.sites:
  site.type = {"A": "Rh", "B": "Zn", "X": "O"}[site.type[0]]
lattice.scale = 8.506

supercell = array([[1, 0, 0],\
                   [0, 1, 0],\
                   [0, 0, 1]], dtype="float64" )
""" Supercell of defect structures. """

vasp = Vasp()
""" VASP functional """
vasp.kpoints    = "Automatic generation\n0\nMonkhorst\n2 2 2\n0 0 0"
vasp.precision  = "accurate"
vasp.ediff      = 1e-5
vasp.encut      = 1
vasp.lorbit     = 10
vasp.npar       = 2
vasp.lplane     = True
vasp.addgrid    = True
vasp.set_smearing   = "metal", 0.01
vasp.set_relaxation = "ionic"
vasp.set_symmetries = "off"

#                Symbol, directory of POTCAR, U parameters, max/min oxidation state, is magnetic
vasp.add_specie = "Rh", "pseudos/Rh", U("liechtenstein", "d", 3.3), 3, True
vasp.add_specie = "Zn", "pseudos/Zn", U("liechtenstein", "d", 6.0), 2, False
vasp.add_specie =  "O",  "pseudos/O", None, -2, True

######################## FAKE! #######################################
vasp.add_specie = "Li", "pseudos/Zn", None, 1
vasp.add_specie = "Na", "pseudos/Zn", None, 1
######################## FAKE! #######################################

# Parameters specific to the ionic relaxation.
relaxation_parameters = {}

# The following sets-up which point defects to do.
substitutions = {"Rh": ["Zn", None], "Zn": ["Rh", None], "O": [None]}

# The following sets-up which which interstitials to examine.
# below, we have:
# Li interstitial on 16c (0,0,0)
# Li interstitial on 32e (x,x,x) with x = 0.75.
# Na interstitial on 16c (0,0,0)
# Na interstitial on 32e (x,x,x) with x = 0.75.
# These are cartesian coordinates in units of lattice.scale (set above).
interstitials = { "Li": [(0,0,0, "16c"), (0.75,0.75,0.75, "32e_0.75")],\
                  "Na": [(0,0,0, "16c"), (0.75,0.75,0.75, "32e_0.75")] }
