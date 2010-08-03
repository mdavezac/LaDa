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


def band_filling_correction(output, cbm):
  """ Returns band-filling corrrection. 

      :Parameters: 
        output : return of lada.vasp.Vasp.__call__
          an output extraction object as returned by the vasp functional when
          computing the defect of interest.
        cbm : float 
          The cbm of the host with potential alignment.
  """
  from numpy import sum
  indices = output.eigenvalues > cbm
  return sum( (output.multiplicity * output.eigenvalues * output.occupations)[indices] - cbm )
  
