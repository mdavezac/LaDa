""" Point-defect helper functions. """
__docformat__ = "restructuredtext en"


def inequivalent_sites(lattice, type):
  """ Yields sites occupied by type which are inequivalent.
  
      When creating a vacancy on, say, "O", or a substitution of "Al" by "Mg",
      there may be more than one site which qualifies. We want to iterate over
      those sites which are inequivalent only, so that only the minimum number
      of operations are performed. 
      :note: lattice sites can be defined as occupiable by more than one atomic type:
             C{lattice.site.type[i] = ["Al", "Mg"]}. These sites will be
             counted if C{type in lattice.site.type}, where type is the input
             parameter.

      :Parameters:
          lattice : `lada.crystal.Lattice`
            Lattice for which to find equivalent sites.
          type : str 
            Atomic specie for which to find inequivalent sites. 

      :return: indices of inequivalent sites.
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
  
      Loops over all symmetrically inequivalent vacancies of given type.

      :Parameters:
        structure : `lada.crystal.Structure`
          structure on which to operate
        lattice : `lada.crystal.Lattice` 
          back-bone lattice of the structure.
        type : str
          type of atoms for which to create vacancy.

      :return: 
        a 3-tuple consisting of:
        - the structure with a vacancy.
        - the vacancy atom from the original structure. It is given an
          additional attribute, C{index}, referring to list of atoms of the
          original structure.
        - A suggested name for the vacancy: site_i, where i is the index of the
          vacancy in the original list of atoms.
  """
  from copy import deepcopy
  
  structure = deepcopy(structure)
  inequivs = inequivalent_sites(lattice, type)
  for i in inequivs:

    # finds first qualifying atom
    for which, atom in enumerate(structure.atoms):
      if atom.site == i: break

    assert which < len(structure.atoms), RuntimeError("Site index not found.")

    # name of this vacancy
    name = ("site_" + str(i)) if len(inequivs) > 1 else ""
    # creates vacancy and keeps atom for record
    atom = deepcopy(structure.atoms.pop(which))
    atom.index = i
    # structure 
    result = deepcopy(structure)
    result.name = structure.name + " %s vacancy on site %i" % (atom.type, i)
    # returns structure with vacancy.
    yield structure, atom, name
    # removes vacancy
    structure.atoms.insert(which, atom)

def substitution(structure, lattice, type, subs):
  """ Yields all inequivalent vacancies. 
  
      Loops over all equivalent vacancies.

      :Parameters:
        structure : `lada.crystal.Structure`
          structure on which to operate
        lattice : `lada.crystal.Lattice`
          back-bone lattice of the structure.
        type : str
          type of atoms for which to create vacancy.
        subs : str
          substitution type

      :return: a 3-tuple consisting of:
        - the structure with a substitution.
        - the substituted atom in the structure above. The atom is given an
          additional attribute, C{index}, referring to list of atoms in the
          structure.
        - A suggested name for the substitution: site_i, where i is the index
          of the substitution in the list of atoms.
  """
  from copy import deepcopy
  
  result = deepcopy(structure)
  inequivs = inequivalent_sites(lattice, type)
  for i in inequivs:

    # finds first qualifying atom
    for which, atom in enumerate(structure.atoms):
      if atom.site == i: break

    assert which < len(structure.atoms), RuntimeError("Site index not found.")

    # name of this substitution
    name = ("site_" + str(i)) if len(inequivs) > 1 else ""
    # creates substitution
    orig = result.atoms[which].type
    result.atoms[which].type = subs
    # substituted atom.
    substituted = deepcopy(result.atoms[which])
    substituted.index = which
    # returns structure with vacancy.
    yield result, substituted, name
    # removes substitution
    result.atoms[which].type = orig

def charged_states(A=None, B=None):
  """ Loops over charged systems. 

      :Types:
        - `A` : None or `lada.vasp.specie.Species`
        - `B` : None or `lada.vasp.specie.Species`

      - if only one of C{A} and C{B} is not None, then the accepted charge
        states are anything between 0 and C{-A.oxidation} included. This
        works both for negative and positive oxidation numbers. The "-" sign
        comes from the rational that an ion C{A} with oxidation
        C{A.oxidation} is removed from the system.
      - if both C{A} and C{B} are not None, than reviewed chared states are between
        C{max(-A.oxidation, B.oxidation-A.oxidation)} and
        C{min(-A.oxidation, B.oxidation-A.oxidation)}. 

      :return: Yields a 2-tuple:
        - Number of electrons to add to the system (not charge). 
        - a suggested name for the charge state calculation.
  """

  if A == None: A, B = B, A
  if A == None: # no oxidation
    yield None, "neutral"
    return
  if B == None:
    # max/min oxidation state
    max_charge = -A.oxidation if hasattr(A, "oxidation") else 0
    min_charge = 0
    if max_charge < min_charge: max_charge, min_charge = min_charge, max_charge
  else:
    # Finds maximum range of oxidation.
    maxA_oxidation = A.oxidation if hasattr(A, "oxidation") else 0
    maxB_oxidation = B.oxidation if hasattr(B, "oxidation") else 0
    max_charge = max(-maxA_oxidation, maxB_oxidation - maxA_oxidation, 0)
    min_charge = min(-maxA_oxidation, maxB_oxidation - maxA_oxidation, 0)

  for charge in range(min_charge, max_charge+1):
    # directory
    if   charge == 0:   oxdir = "neutral"
    elif charge > 0:    oxdir = "+" + str(charge) 
    elif charge < 0:    oxdir = str(charge) 
    yield -charge, oxdir


def band_filling_correction(defect, cbm):
  """ Returns band-filling corrrection. 

      :Parameters: 
        defect : return of `lada.vasp.Vasp.__call__`
          An output extraction object as returned by the vasp functional when
          computing the defect of interest.
        cbm : float 
          The cbm of the host with potential alignment.
  """
  from numpy import sum
  indices = defect.eigenvalues > cbm
  return sum( (defect.multiplicity * defect.eigenvalues * defect.occupations)[indices] - cbm )
  

def potential_alignment(defect, host, maxdiff=0.5):
  """ Returns potential alignment correction. 

      :Parameter:
        defect : return of `lada.vasp.Vasp.__call__`
          An output extraction object as returned by the vasp functional when
          computing the defect of interest.
        host : return of `lada.vasp.Vasp.__call__`
          An output extraction object as returned by the vasp functional when
          computing the host matrix.
        maxdiff : float
          Maximum difference between the electrostatic potential of an atom and
          the equivalent host electrostatic potential beyond which that atom is
          considered pertubed by the defect.

      :return: The potential alignment in eV (without charge factor).
  """
  from operator import itemgetter
  from numpy import average, array, abs
  from crystal import specie_list

  # first get average atomic potential per atomic specie in host.
  host_electropot = {}
  host_species = specie_list(host.structure)
  for s in host_species:
    indices = array([atom.type == s for atom in structure.atoms])
    host_electropot[s] = average(host.electropot[indices]) 
  
  # creates an initial list of unperturbed atoms. 
  types = array([atom.type for atom in defect.structure.atoms])
  unperturbed = array([(type in host_species) for type in types])

  # now finds unpertubed atoms according to maxdiff
  while any(unperturbed): # Will loop until no atom is left. That shouldn't happen!

    # computes average atomic potentials of unperturbed atoms.
    defect_electropot = {}
    for s in host_species:
      defect_electropot[s] = average( defect.electropot[unperturbed & (types == s)] )

    # finds atomic potential farthest from the average potentials computed above.
    discrepancies =   defect.electropot[unperturbed] \
                    - array([host_electropot[type] for type in types[unperturbed]])
    index, diff = max(enumerate(abs(discrepancies)), key=itemgetter(2))

    # if discrepancy too high, mark atom as perturbed.
    if diff > maxdiff: unperturbed[unperturbed][index] = False
    # otherwise, we are done for this loop!
    else: break

  # now computes potential alignment from unpertubed atoms.
  diff_from_host =    defect.electropot[unpertubed]  
                    - array([host_electropot[type] for type in types[unperturbed]])\
  return average(diff_from_host)
                    

def third_order_charge_correction(cell, n = 200):
  """ Returns energy of third order charge correction. 
  
      :Parameters: 
        cell : 3x3 numpy array
          The supercell of the defect in Angstroem(?).
        n 
          precision. Higher better.
      :return: energy of a single negative charge in supercell.
  """
  from math import pi
  from numpy import array, dot, det
  from ..physics import Rydberg

  def fold(vector):
    """ Returns smallest distance. """
    result = None
    for i in range(-1, 2):
      for j in range(-1, 2):
        for k in range(-1, 2):
          v = arrray([vector[0] + float(i), vector[1] + float(j), vector[2] + float(k)])
          v = dot(cell, v)
          m = dot(v,v)
          if result == None or result > m: result = m
    return result

  # chden = ones( (n, n, n), dtype="float64")
  result = 0e0
  for ix in range(n):
    for iy in range(n):
      for iz in range(n):
        vec = array([float(ix)/float(n)-0.5, float(iy)/float(n)-0.5, float(iz)/float(n)-0.5])
        result += fold(vec)
        
  return -result / float(n**3) * Rydberg("eV") * pi * 4e0 / 3e0 / det(cell)
