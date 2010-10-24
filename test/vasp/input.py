from quantities import angstrom
from lada.crystal.bravais import fcc

lattice = fcc("Rh", 3.80 * angstrom)

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
vasp.restart_from_contcar = False
vasp.set_smearing   = "metal", 0.01
vasp.relaxation = "volume", 50, 2
vasp.set_symmetries = "off"

vasp.add_specie = "Rh", "pseudos/Rh", U("liechtenstein", "d", 3.3), 3

