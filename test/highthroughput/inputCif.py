from glob import iglob
from lada.crystal import icsd_cif
########################################

vasp = Vasp()
""" VASP functional """
vasp.precision  = "accurate"
vasp.ediff      = 1e-4
vasp.encut      = 340.0
vasp.npar       = 2
vasp.lplane     = True
vasp.addgrid    = True
vasp.set_smearing   = "gaussian", 0.01
vasp.relaxation = "ionic"
vasp.set_symmetries = "off"
vasp.kpoints        = "\n0\nAuto\n20"
vasp.lorbit         = 10
vasp.add_param      = "lmaxmix",4


# "Al" => specie symbol
# "pseudos/Al" => directory where relevant POTCAR is located

vasp.add_specie = "Ti", "/scratch/vnsteva/lada/pseudos/Ti", U("dudarev", "d", 3.0)
vasp.add_specie = "O", "/scratch/vnsteva/lada/pseudos/O"

vasp.species["Ti"].moment = [1.e0]

materials = {}

for name in iglob("icsd_structures/*.cif"):
    materials[name[name.index('/')+1:-4]]=icsd_cif(name)

#########################################################

first_trial = { "kpoints": "\n0\nAuto\n10", "encut": 0.9 }
""" parameter to override during first relaxation step. """
relaxation_dof = "volume ionic cellshape"
#relaxation_dof = "volume ionic"
#relaxation_dof = "ionic"
#relaxation_dof = "static"
""" Degrees of freedom to relax. """
relaxer = RelaxCellShape( vasp, relaxation_dof, first_trial, maxiter=5)
""" Cell shape relaxation algorithm. """

""" Materials to compute. """
nbantiferro = 3
""" Number of random anti-ferro trials. """

