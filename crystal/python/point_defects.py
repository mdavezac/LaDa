# will be params:
#Aatoms = ["Ir", "Rh"]
#Batoms = ["Mg", "Zn"]
#defect_types = ["VA", "VB", "VX", "AonB", "BonA"]

Aatoms = ["Li"]
Batoms = ["Mg"]
defect_types = ["VA", "AonB"]

def inequivalent_sites(lattice, type):
  """ Yields sites occupied by type which are inequivalent.
  
      When creating a vacancy on, say, "O", or a substitution of "Al" by "Mg",
      there may be more than one site which qualifies. We want to iterate over
      those sites which are inequivalent only, so that only the minimum number
      of operations are performed. 
      @note: lattice sites can be defined as occupiable by more than one atomic type:
             C{lattice.site.type[i] = ["Al", "Mg"]}. These sites will be
             counted if C{type in lattice.site.type}, where type is the input
             parameter.
      @param lattice: Lattice for which to find equivalent sites.
      @type lattice: lada.crystal.Lattice
      @param type: Atomic specie for which to find inequivalent sites. 
      @return: indices of inequivalent sites.
  """
  from numpy.linalg import inv, norm
  from lada.crystal import fold_vector

  # all sites with occupation "type". 
  sites = [site for site in lattice.sites if type in site.type]
  site_indices = [i for i,site in enumerate(lattice.sites) if type in site.type]

  # inverse cell.
  invcell = inv(lattice.cell)
  # loop over all site with type occupation.
  i = 0
  while i < len(sites):
    # iterates over symmetry operations.
    for op in lattice.space_group:
      pos = op(site.pos)
      # finds index of transformed position, using translation quivalents.
      for t, other in enumerate(sites):
        if norm(fold_vector(pos, lattice.cell, invcell)) < 1e-12:
          print t
          break
      # removes equivalent site and index from lists if necessary
      if t != i and t < len(sites): 
        sites.pop(t)
        site_indices.pop(t)
    i += 1
  return site_indices

def vacancy(structure, lattice, type):
  """ Yields all inequivalent vacancies. 
  
      Loops over all equivalent vacancies.
      @param structure: structure on which to operate
      @param lattice: back-bone lattice of the structure.
      @param type: type of atoms for which to create vacancy.
  """
  from copy import deepcopy
  
  result = deepcopy(structure)
  for i in inequivalent_sites(lattice, type):

    # finds first qualifying atom
    for which, atom in enumerate(structure.atoms):
      if atom.site == i: break

    assert which < len(structure.atoms), RuntimeError("Site index not found.")

    # creates vacancy
    atom = result.atoms.pop(which)
    # returns structure with vacancy.
    yield deepcopy(result), deepcopy(atom), "vacancy_" + atom.type
    # removes vacancy
    result.atoms.insert(which, atom)

def substitution(structure, lattice, type, subs):
  """ Yields all inequivalent vacancies. 
  
      Loops over all equivalent vacancies.
      @param structure: structure on which to operate
      @param lattice: back-bone lattice of the structure.
      @param type: type of atoms for which to create vacancy.
      @param subs: substitution type
  """
  from copy import deepcopy
  
  result = deepcopy(structure)
  for i in inequivalent_sites(lattice, type):

    # finds first qualifying atom
    for which, atom in enumerate(structure.atoms):
      if atom.site == i: break

    assert which < len(structure.atoms), RuntimeError("Site index not found.")

    # creates substitution
    orig = result.atoms[which].type
    result.atoms[which].type = subs
    substituted = deepcopy(result.atoms[which])
    substituted.index = which
    # returns structure with vacancy.
    yield deepcopy(result), substituted, orig + "_on_" + subs
    # removes substitution
    result.atoms[which].type = orig

def oxidation(A=None, B=None):

  if A == None: A, B = B, A
  if A == None: # no oxidation
    yield None, "neutral"
    return
  if B == None:
    # max/min oxidation state
    max_oxidation = -A.oxidation if hasattr(A, "oxidation") else 0
    min_oxidation = 0
    if max_oxidation < min_oxidation: max_oxidation, min_oxidation = min_oxidation, max_oxidation
  else:
    # Finds maximum range of oxidation.
    maxA_oxidation = -A.oxidation if hasattr(A, "oxidation") else 0
    maxB_oxidation = -B.oxidation if hasattr(B, "oxidation") else 0
    max_oxidation = max(maxA_oxidation, maxA_oxidation - maxB_oxidation, 0)
    min_oxidation = min(maxA_oxidation, maxA_oxidation - maxB_oxidation, 0)

  for oxidation in range(min_oxidation, max_oxidation+1):
    # directory
    if   oxidation == 0:   oxdir = "neutral"
    elif oxidation > 0:    oxdir = "+" + str(oxidation) 
    elif oxidation < 0:    oxdir = str(oxidation) 
    yield oxidation, oxdir


if __name__ == '__main__':
  from sys import exit
  from os.path import join
  from lada.opt import read_input
  from lada.vasp import Vasp, Specie, specie
  from lada.vasp.incar import Standard, NElect
  from lada.crystal import fill_structure
  import jobs

  input_dict = {
                 "Specie": Specie, "Vasp":Vasp, "Standard":Standard,
                 "U": specie.U, "nlep": specie.nlep
               }
  input = read_input("input.py", input_dict)
  input.lattice.set_as_crystal_lattice()
  for site in input.lattice.sites:
    if   "A" in site.type: site.type[0] = input.species[0].symbol
    elif "B" in site.type: site.type[0] = input.species[1].symbol
    elif "X" in site.type: site.type[0] = input.species[2].symbol

# print "# input lattice. \n", input.lattice

  supercell = fill_structure(input.supercell)
# print "# supercell structure. \n", supercell

  # name of output directories
  outdir = "%s2%s%s4" % tuple([u.symbol for u in input.species])

  # vacancy
# for vac in input.species[:2]: # not doing oxygen right now.
#   for structure, atom in vacancy(supercell, input.lattice, vac.symbol):
#     vacdir = join(outdir, "vacancy_" + vac.symbol)
#     if not hasattr(vac, "oxidation"): 
#       jobs.current[vacdir].job["structure"] = structure
#       jobs.current[vacdir].job["vasp"] = input.vasp
#     else:
#       for i in range(0, vac.oxidation, 1 if vac.oxidation > 0 else -1): 
#         directory = join(vacdir, str(i) + "-" if vac.oxidation > 0 else "+")
#         jobs.current[directory].job["structure"] = structure
#         jobs.current[directory].job["nelect"] = NElect(-i)
#         jobs.current[directory].job["vasp"] = input.vasp
  
  jobs.current.update(jobs.load())
  for job in jobs.walk_through():
    print job.job["outdir"]
    if "nelect" in job.job: print job.job["nelect"]

# # substitutions.
# subs = [(input.species[0].symbol, input.species[1].symbol),\
#         (input.species[1].symbol, input.species[0].symbol)]
# for A, B in subs:
#  for structure in substitution(supercell, input.lattice, A, B):
#    directory = join(outdir, "%s_on_%s" % (A, B))
#    print directory, len([0 for atom in structure.atoms if atom.type == B])


