""" Some definitions of magnetic order used in high-throughput. """
def is_magnetic_system(structure, species):
  """ True if system is magnetic. """
  from lada.crystal import specie_list

  species = [u for name, u in species.items() if name in specie_list(structure)]
  return len([0 for u in species if hasattr(u, "moment")]) != 0



def ferro(structure, species):
  """ Returns magmom VASP flag for ferromagnetic order. """
  from lada.crystal import specie_list

  names = specie_list(structure)
  pseudos = [u for name, u in species.items() if name in names]
  numbers = [atom.type for atom in structure.atoms]
  numbers = [numbers.count(u) for u in names]

  magmom = ""
  for n, ps in zip(numbers, pseudos): 
    magmom +="%i*%.2f " % (n, ps.moment if hasattr(ps, "moment") else 0)
  return magmom

def sublatt_antiferro(structure, species):
  """ Anti ferro order with each cation type in a different direction. """
  from lada.crystal import specie_list

  names = specie_list(structure)
  pseudos = [u for name, u in species.items() if name in names]
  if len([0 for ps in pseudos if ps.magnetic]) < 2: return None
  numbers = [atom.type for atom in structure.atoms]
  numbers = [numbers.count(u) for u in names]

  magmom, sign = "", 1e0
  for n, ps  in zip(numbers, pseudos): 
    magmom +="%i*%.2f " % (n, sign * ps.moment if hasattr(ps, "moment") else 0)
    if ps.magnetic: sign *= -1e0
  return magmom

def random(structure, species):
  """ Random magnetic order. """
  from random import uniform
  from lada.crystal import specie_list

  names = specie_list(structure)
  pseudos = [u for name, u in species.items() if name in names]
  numbers = [atom.type for atom in structure.atoms]
  numbers = [numbers.count(u) for u in names]

  magmom = ""
  for n, ps  in zip(numbers, pseudos): 
    if hasattr(ps, "moment"):
      for i in range(n):
        magmom += "%.2f " % ( uniform(-ps.moment, ps.moment) )
    else: magmom += "%i*0   " % (n)
  return magmom

