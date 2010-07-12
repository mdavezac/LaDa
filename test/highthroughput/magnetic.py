""" Some definitions of magnetic order used in high-throughput. """
def is_magnetic_system(structure, species):
  """ True if system is magnetic. """
  from lada.crystal import specie_list 

  species = [u for name, u in species.items() if name in specie_list(structure)]
  return len([0 for u in species if u.magnetic]) != 0


def ferro(structure, species):
  """ Returns magmom VASP flag for ferromagnetic order. """
  from lada.crystal import specie_list

  names = specie_list(structure)
  pseudos = [u for name, u in species.items() if name in names]
  moments = [0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0] # high-spin, d elctrons.
  numbers = [atom.type for atom in structure.atoms]
  numbers = [numbers.count(u) for u in names]
  valence = [int(ps.valence) for ps in pseudos]

  magmom = ""
  for n, ps, val  in zip(numbers, pseudos, valence): 
    magmom +="%i*%f   " % (n, moments[val] if ps.magnetic else 0)
  return magmom

def sublatt_antiferro(structure, species):
  """ Anti ferro order with each cation type in a different direction. """
  from lada.crystal import specie_list

  names = specie_list(structure)
  pseudos = [u for name, u in species.items() if name in names]
  if len([0 for ps in pseudos if ps.magnetic]) < 2: return None
  moments = [0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0] # high-spin, d elctrons.
  numbers = [atom.type for atom in structure.atoms]
  numbers = [numbers.count(u) for u in names]
  valence = [int(ps.valence) for ps in pseudos]

  magom, sign = "", 1e0
  for n, ps, val  in zip(numbers, pseudos, valence): 
    magmom +="%i*%f   " % (n, sign * moments[ps.valence] if ps.magnetic else 0)
    if ps.magnetic: sign *= -1e0
  return magmom

def random(structure, species):
  """ Random magnetic order. """
  from random import uniform
  from lada.crystal import specie_list

  names = specie_list(structure)
  pseudos = [u for name, u in species.items() if name in names]
  moments = [0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0] # high-spin, d elctrons.
  numbers = [atom.type for atom in structure.atoms]
  numbers = [numbers.count(u) for u in names]
  valence = [int(ps.valence) for ps in pseudos]

  magmom = ""
  for n, ps, val  in zip(numbers, pseudos, valence): 
    if ps.magnetic:
      for i in range(n):
        magmom += "%f " % ( uniform(-moments[val], moments[val]) )
    else: magmom += "%i*0   " % (n)

