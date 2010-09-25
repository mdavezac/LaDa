from lada.crystal import A2BX4
vasp = Vasp()
""" VASP functional """
vasp.precision  = "accurate"
vasp.ediff      = 1e-5
vasp.encut      = 1
vasp.npar       = 2
vasp.lplane     = True
vasp.addgrid    = True
vasp.set_smearing   = "metal", 0.01
vasp.set_relaxation = "ionic"
vasp.set_symmetries = "off"
vasp.kpoints        = "\n0\nAuto\n30"

# "Al" => specie symbol
# "pseudos/Al" => directory where relevant POTCAR is located
vasp.add_specie = "Al", "pseudos/Al"
vasp.add_specie = "Mg", "pseudos/Mg"
vasp.add_specie =  "O",  "pseudos/O"
vasp.species["Mg"].moment = [5e0, 0.5e0]
vasp.species["Al"].moment = [5e0, 0.5e0]

first_trial = { "kpoints": "\n0\nAuto\n1", "encut": 0.9 }
""" parameter to override during first relaxation step. """
relaxation_dof = "volume ionic cellshape"
""" Degrees of freedom to relax. """
relaxer = RelaxCellShape( vasp, relaxation_dof, first_trial, maxiter=5)
""" Cell shape relaxation algorithm. """


def scale(structure):
  """ Returns *guessed* scale (eg volume^(0.33)) for a given structure. """
  from numpy.linalg import det
  if "O" in [atom.type for atom in structure.atoms]:    spvol = 8.5**3/4e0
  elif "Se" in [atom.type for atom in structure.atoms]: spvol = 9.5**3/4e0
  elif "Te" in [atom.type for atom in structure.atoms]: spvol = 10.5**3/4e0
  else: raise ValueError("Neither O, nor Se, nor Te atom found.")

  nfu = float(len(structure.atoms)/7)*0.5 # 0.5 because 2 f.u. in spinel unit-cell.
  vol = det(structure.cell)
  return (nfu * spvol / vol)**(1e0/3e0) 

materials = ["Al2MgO4"]
""" Materials to compute. """
nbantiferro = 3
""" Number of random anti-ferro trials. """

lattices = [ A2BX4.b1(),  A2BX4.b10()] #,  A2BX4.b10I(),  A2BX4.b11() ] #,
#            A2BX4.b12(),  A2BX4.b15(),  A2BX4.b16(),  A2BX4.b18(),
#            A2BX4.b19(), A2BX4.b1I(),  A2BX4.b2(),  A2BX4.b20(),  A2BX4.b21(),
#            A2BX4.b2I(), A2BX4.b33(),  A2BX4.b34(),  A2BX4.b35(),
#            A2BX4.b36(),  A2BX4.b37(), A2BX4.b38(),  A2BX4.b4(),  A2BX4.b4I(),
#            A2BX4.b5(),  A2BX4.b5I(),  A2BX4.b6(), A2BX4.b7(),  A2BX4.b8(),
#            A2BX4.b9(),  A2BX4.b9I(),  A2BX4.d1(),  A2BX4.d1I(), A2BX4.d3(),
#            A2BX4.d3I(),  A2BX4.d9(),  A2BX4.s1(),  A2BX4.s1I(),  A2BX4.s2(),
#            A2BX4.s2I(),  A2BX4.s3(),  A2BX4.s3I() ]
""" Lattice for which to create structures. """
