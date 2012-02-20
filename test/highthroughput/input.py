from lada.crystal import ABX 
vasp = Vasp()
""" VASP functional """

vasp.precision      = "accurate"
vasp.ediff          = 1e-5 # precision per ATOM
vasp.encut          = 340.0    # if smaller than 3 than scale
vasp.npar           = 2
#vasp.lplane         = True
vasp.addgrid        = True
vasp.set_smearing   = "metal", 0.03
vasp.relaxation = "ionic"
vasp.set_symmetries = "off"
vasp.kpoints        = "\n0\nAuto\n20"
vasp.lorbit         = 10
# Magmom is a special parameter with specific implementation for point defects.
# We need to go through add_item to reset it as normal key/value parameter.
#vasp.add_param      = "magmom", None
vasp.add_param      = "lmaxmix",4

#################################################################################
vasp.add_specie = "H", "/global/homes/h/hwpeng/tools/vasp/pbe/H", None, 1

vasp.add_specie = "Mg", "/global/homes/h/hwpeng/tools/vasp/pbe/Mg", None, 2
vasp.add_specie = "Ca", "/global/homes/h/hwpeng/tools/vasp/pbe/Ca_pv", None, 2
vasp.add_specie = "Ba", "/global/homes/h/hwpeng/tools/vasp/pbe/Ba_sv", None, 2
vasp.add_specie = "Sr", "/global/homes/h/hwpeng/tools/vasp/pbe/Sr_sv", None, 2

vasp.add_specie = "Al", "/global/homes/h/hwpeng/tools/vasp/pbe/Al", None,  3
vasp.add_specie = "Ga", "/global/homes/h/hwpeng/tools/vasp/pbe/Ga_d", None, 3
vasp.add_specie = "In", "/global/homes/h/hwpeng/tools/vasp/pbe/In_d", None ,3

vasp.add_specie = "Si", "/global/homes/h/hwpeng/tools/vasp/pbe/Si", None, 4
vasp.add_specie = "Ge", "/global/homes/h/hwpeng/tools/vasp/pbe/Ge_d", None, 4
vasp.add_specie = "Sn", "/global/homes/h/hwpeng/tools/vasp/pbe/Sn_d", None, 4

vasp.add_specie = "O", "/global/homes/h/hwpeng/tools/vasp/pbe/O_s", None, -2

vasp.add_specie = "Sc", "/global/homes/h/hwpeng/tools/vasp/pbe/Sc_sv", U("dudarev", "d", 3.0),3
vasp.add_specie = "Ti", "/global/homes/h/hwpeng/tools/vasp/pbe/Ti_pv",U("dudarev", "d", 3.0), 4
vasp.add_specie = "V", "/global/homes/h/hwpeng/tools/vasp/pbe/V_pv", U("dudarev", "d", 3.0), 5
vasp.add_specie=  "Cr", "/global/homes/h/hwpeng/tools/vasp/pbe/Cr_pv", U("dudarev", "d", 3.0), 3
vasp.add_specie = "Mn", "/global/homes/h/hwpeng/tools/vasp/pbe/Mn",  U("dudarev", "d", 3.0 ), 2
vasp.add_specie = "Fe", "/global/homes/h/hwpeng/tools/vasp/pbe/Fe", U("dudarev", "d", 3.0), 2
vasp.add_specie = "Co", "/global/homes/h/hwpeng/tools/vasp/pbe/Co", U("dudarev", "d", 3.0), 3
vasp.add_specie = "Ni", "/global/homes/h/hwpeng/tools/vasp/pbe/Ni", U("dudarev", "d", 3.0), 2

vasp.add_specie = "Cu", "/global/homes/h/hwpeng/tools/vasp/pbe/Cu", U("dudarev", "d" , 5.0), 1
vasp.add_specie	= "Ag", "/global/homes/h/hwpeng/tools/vasp/pbe/Ag", U("dudarev", "d", 5.0), 1 

vasp.add_specie = "Zn", "/global/homes/h/hwpeng/tools/vasp/pbe/Zn", U("dudarev", "d", 6.0), 2
vasp.add_specie = "Cd", "/global/homes/h/hwpeng/tools/vasp/pbe/Cd", U("dudarev", "d" , 6.0), 2

#
vasp.add_specie = "Y", "/global/homes/h/hwpeng/tools/vasp/pbe/Y_sv", U("dudarev", "d", 3.0), 3
vasp.add_specie = "Zr", "/global/homes/h/hwpeng/tools/vasp/pbe/Zr_sv", U("dudarev", "d", 3.0), 4
vasp.add_specie = "Nb", "/global/homes/h/hwpeng/tools/vasp/pbe/Nb_pv", U("dudarev", "d", 3.0), 5
vasp.add_specie = "Mo", "/global/homes/h/hwpeng/tools/vasp/pbe/Mo_pv", U("dudarev", "d", 3.0), 3
# Tc
# Ru
vasp.add_specie = "Rh", "/global/homes/h/hwpeng/tools/vasp/pbe/Rh", U("dudarev", "d", 3.0), 3
# Pd

vasp.add_specie = "Hf", "/global/homes/h/hwpeng/tools/vasp/pbe/Hf_pv", U("dudarev", "d", 3.0), 2
vasp.add_specie = "Hf", "/global/homes/h/hwpeng/tools/vasp/pbe/Hf_pv", U("dudarev", "f", 6.0), 2
vasp.add_specie = "W", "/global/homes/h/hwpeng/tools/vasp/pbe/W", U("dudarev", "d", 3.0), 5
vasp.add_specie = "Ir", "/global/homes/h/hwpeng/tools/vasp/pbe/Ir", U("dudarev", "d", 3.0), 3
#
vasp.add_specie = "Ta", "/global/homes/h/hwpeng/tools/vasp/pbe/Ta_pv", U("dudarev", "d", 3.0),3

vasp.add_specie = "La", "/global/homes/h/hwpeng/tools/vasp/pbe/La", ("dudarev", "d", 3.0), 2  
vasp.add_specie = "La", "/global/homes/h/hwpeng/tools/vasp/pbe/La", ("dudarev", "f", 6.0), 2  

vasp.species["Cu"].moment = [ 8.0, 0.6 ] 
vasp.species["Cr"].moment = [ 8.0, 0.6 ]
vasp.species["Ni"].moment = [ 4.0, 0.6 ]
vasp.species["Co"].moment = [ 4.0, 0.6 ] 
vasp.species["Fe"].moment = [ 8.0, 0.6 ]
vasp.species["Mn"].moment = [ 8.0, 0.6 ]
vasp.species["Cr"].moment = [ 8.0, 0.6 ]
vasp.species["V"].moment =  [ 8.0, 0.6 ]
vasp.species["Ir"].moment = [ 4.0, 0.6 ]
vasp.species["Rh"].moment = [ 8.0, 0.6 ]
vasp.species["Ti"].moment = [ 5.0, 0.6 ]
vasp.species["Sc"].moment = [ 8.0, 0.6 ]
##################################################################################

#first_trial = { "kpoints": "\n0\nAuto\n40", "encut": 1.0 }
first_trial = { }
""" parameter to override during first relaxation step. """
relaxation_dof = "volume ionic cellshape"
""" Degrees of freedom to relax. """
relaxer = RelaxCellShape( vasp, relaxation_dof, first_trial, maxiter=5, keep_steps=False)
""" Cell shape relaxation algorithm. """

def scale(structure):
  """ Returns *guessed* scale (eg volume^(0.33)) for a given structure. """
  from numpy.linalg import det
  if "O" in [atom.type for atom in structure.atoms]:    spvol = 8.5**3/4e0
  elif "Se" in [atom.type for atom in structure.atoms]: spvol = 9.5**3/4e0
  elif "Te" in [atom.type for atom in structure.atoms]: spvol = 10.5**3/4e0

  nfu = float(len(structure.atoms)/7)*0.5 # 0.5 because 2 f.u. in spinel unit-cell.
  vol = det(structure.cell)
  return (nfu * spvol / vol)**(1e0/3e0)

#problem with Rh2TiO4 for some reason

""" Materials to compute. """
materials = [ "CrMnO"]

""" Number of random anti-ferro trials. """
nbantiferro = 8
nbrandom    = 3
do_ferro    = False
do_antiferro = False

lattices = [ ABX.s36() ]
#lattices = [ A2BX4.b5()]

 
#first_wave_dir = "/scratch/trpaude/lada/hte/"
