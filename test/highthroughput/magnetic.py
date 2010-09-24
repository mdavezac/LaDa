""" Some definitions of magnetic order used in high-throughput. """
def is_magnetic_system(structure, species):
  """ True if system is magnetic. """
  from lada.crystal import specie_list

  species = [u for name, u in species.items() if name in specie_list(structure)]
  return len([0 for u in species if hasattr(u, "moment")]) != 0

def ls_ferro(structure, species):
  """ Returns magmom VASP flag for low-spin ferromagnetic order. """
  from lada.crystal import specie_list

  magmom = ""
  for specie in specie_list(structure):
    assert specie in species,\
           KeyError("specie {0} not found in pseudo-potential dictionary.".format(specie))
    atoms = [atom for atom in structure.atoms if atom.type == specie]
    if hasattr(species[specie], "moment"):
      moment = species[specie].moment
      if hasattr(moment, "__iter__"): moment = min(moment)
      magmom += "{0}*{1:.2f} ".format(len(atoms), moment)
    else: magmom += "{0}*0 ".format(len(atoms), 0)
  return magmom

def hs_ferro(structure, species):
  """ Returns magmom VASP flag for high-spin ferromagnetic order. """
  from lada.crystal import specie_list

  magmom, has_both = "", False
  for specie in specie_list(structure):
    assert specie in species,\
           KeyError("specie {0} not found in pseudo-potential dictionary.".format(specie))
    atoms = [atom for atom in structure.atoms if atom.type == specie]
    if hasattr(species[specie], "moment"):
      moment = species[specie].moment
      if hasattr(moment, "__iter__"):
        moment = max(moment)
        has_both = True
      magmom += "{0}*{1:.2f} ".format(len(atoms), moment)
    else: magmom += "{0}*0 ".format(len(atoms), 0)
  return magmom if has_both else None

def sublatt_antiferro(structure, species):
  """ Anti ferro order with each cation type in a different direction. """
  from lada.crystal import specie_list

  magmom, sign, nb = "", 1e0, 0
  for specie in specie_list(structure):
    assert specie in species,\
           KeyError("specie {0} not found in pseudo-potential dictionary.".format(specie))
    atoms = [atom for atom in structure.atoms if atom.type == specie]
    if hasattr(species[specie], "moment"):
      moment = species[specie].moment
      if hasattr(moment, "__iter__"): moment = max(moment)
      magmom += "{0}*{1:.2f} ".format(len(atoms), moment * sign)
      nb += 1
      sign *= -1e0
    else: magmom += "{0}*0 ".format(len(atoms), 0)
  return None if nb < 2 else magmom

def random(structure, species):
  """ Random magnetic order. """
  from random import uniform
  from lada.crystal import specie_list

  magmom = ""
  for specie in specie_list(structure):
    assert specie in species,\
           KeyError("specie {0} not found in pseudo-potential dictionary.".format(specie))
    atoms = [atom for atom in structure.atoms if atom.type == specie]
    ps = species[specie]
    if hasattr(species[specie], "moment"): 
      moment = species[specie].moment
      if hasattr(moment, "__iter__"): moment = max(moment)
      for atom in atoms: magmom += "{0:.2f} ".format( uniform(-moment, moment) )
    else: magmom += "{0}*0   ".format(len(atoms))

  return magmom

