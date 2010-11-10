""" Point-defect helper functions. """
__docformat__ = "restructuredtext en"
__all__ = [ 'inequivalent_sites', 'vacancy', 'substitution', 'charged_states', \
            'band_filling', 'potential_alignment', 'charge_corrections', \
            'magmom', 'low_spin_states', 'high_spin_states', 'magname', \
            'ExtractSingle', 'ExtractMaterial' ]

from extract import Single as ExtractSingle, Material as ExtractMaterial

def inequivalent_sites(lattice, type):
  """ Yields sites occupied by type which are inequivalent. 
  
      When creating a vacancy on, say, "O", or a substitution of "Al" by "Mg",
      there may be more than one site which qualifies. We want to iterate over
      those sites which are inequivalent only, so that only the minimum number
      of operations are performed. 

      :note:
        lattice sites can be defined as occupiable by more than one atomic type\:
        lattice.site.type[i] = ["Al", "Mg"]. These sites will be counted if
        type in lattice.site.type, where type is the input parameter.
 
      :Parameters:
          lattice : `lada.crystal.Lattice`
            Lattice for which to find equivalent sites.
          type : str 
            Atomic specie for which to find inequivalent sites. 

      :return: indices of inequivalent sites.
  """
  from numpy.linalg import inv, norm
  from .. import fold_vector

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


def non_interstitials(structure, indices, mods):
  """ Yields substitutions for given indices to structure. """
  from copy import deepcopy
  from regex import compile
  from re import compile
  if not hasattr(indices, '__len__'): indices = [indices]

  cation_regex  = compile(r'^\s*cations?\s*$')
  vacancy_regex = compile(r'^\s*vacanc(?:y|ies)\s*$')
  specie_regex  = compile(r'^\s*[A-Z][a-z]?(?:\d+)\s*$')

  # modify input to something which makes sense and check it.
  if mods == None: all_mods = [None]
  elif isinstance(mods, str): all_mods = [mods]
  else: all_mods = mods
  mods = []
  whatthehell = []
  for type in all_mods:
    if type == None: mods.append(None)
    elif cation_regex.match(type.lower()) != None:  mods.extend(_cationic_species(structure)) 
    elif vacancy_regex.match(type.lower()) != None: mods.append(None)
    elif specie_regex.match(type) != None: mods.append(type)
    else: whatthehell.append(type)
  assert len(whatthehell) == 0,\
         ValueError('Cannot understand following specie types: {0}.'.format(whatthehell))
  mods = list(set(mods))

  # loop over atoms to modify.
  for i, j in enumerate(indices):
    assert j < len(structure.atoms), RuntimeError("Site index not found.")
    atom = deepcopy(structure.atoms[j])
    atom.index = j
    # loop over atoms to modifications.
    for modif in mods:
      result = deepcopy(structure)
      if modif == atom.type: continue # can't substitute with self.
      if modif == None: 
        result.atoms.pop(j)
        name = "vacancy_{0}".format(atom.type)
        B = atom.type
        atom.type = 'None'
      else: 
        result.atoms[j].type = modif
        B = atom.type
        atom.type = modif
        name = "{0}_on_{1}".format(atom.type, B)
      if len(indices) > 1: name += "/site_{0}".format(i)
      result.name = name
      yield result, atom, B

def inequiv_non_interstitials(structure, lattice, type, mods):
  """ Loop over inequivalent non-interstitials. """
  inequivs = inequivalent_sites(lattice, type)
  indices = []
  for i in inequivs:
    # finds first qualifying atom
    for which, atom in enumerate(structure.atoms):
      if atom.site == i: break
    indices.append(which)
  for result in non_interstitials(structure, indices, mods): yield result


def interstitials(structure, lattice, interstitials):
  """ Yields interstitial. """
  from copy import deepcopy
  from numpy import dot
  for desc in interstitials:
    assert hasattr(desc, "__iter__"),\
           ValueError("For interstitials, desc should be a sequence: {0}".format(desc))
    assert len([u for u in desc]) == 3,\
           ValueError("For interstitials, desc should be a sequence of length 3: {0}".format(desc))
    type, position, name = tuple([u for u in desc])
    result = deepcopy(structure)
    result.add_atom = dot(lattice.cell, position), type
    result.name = "{0}_interstitial_{1}".format(type, name)
    defect = deepcopy(structure.atoms[-1])
    defect.index = -1
    yield result, defect, 'None'

def _cationic_species(structure):
  """ Returns list of cationic species. """
  return list(set([a.type for a in structure.atoms if a.type not in ['O', 'S', 'Se', 'Te']]))

def iterdefects(structure, lattice, defects):
  """ Iterates over all defects for any number of types and modifications. """
  from re import compile

  cation_regex       = compile(r'^\s*cations?\s*$')
  cation0_regex      = compile(r'^\s*cations?(\d+)\s*$')
  interstitial_regex = compile(r'^\s*interstitials?\s*$')
  type_regex         = compile(r'^\s*(?:vacanc(?:y|ies)|[A-Z][a-z]?(?:\d+)?)\s*$')
  already_cations    = False

  for key, value in defects.items():
    if key == None: keys = [None]
    elif interstitial_regex.match(key.lower()) != None: keys = [None]
    elif cation_regex.match(key.lower()) != None:
      assert already_cations == False, ValueError('Can only have one cations?\d* tag.')
      already_cations = True
      keys = _cationic_species(structure)
    elif cation0_regex.match(key.lower()) != None:
      assert already_cations == False, ValueError('Can only have one cations?\d* tag.')
      already_cations = True
      d = int(cation0_regex.match(key.lower()).group(1))
      keys = ['{0}{1}'.format(k, d) for k in _cationic_species(structure)]
    else: keys = [key]
    for type in keys:
      assert type == None or type_regex.match(type) != None,\
             ValueError("Cannot understand type {0}.".format(type))
      for result in any_defect(structure, lattice, type, value): yield result

def any_defect(structure, lattice, type, subs):
  """ Yields point-defects of a given type and different modifications. 
  
      Loops over all equivalent point-defects.

      :Parameters:
        structure : `lada.crystal.Structure`
          structure on which to operate
        lattice : `lada.crystal.Lattice`
          back-bone lattice of the structure.
        type : str or None or sequence
          type of atoms for which to create substitution.
          If None, will create a vacancy.
          If ``subs`` is None, then will create a vacancy. In that case, type
          should be a sequence describing the interstitials:

          >>> type = [ "Li", (0,0,0), "16c" ],\
          >>>        [ "Li", (0.75,0.75,0.75), "32e_0.75" ] 
           
          Each item in the sequence is itself a sequence where the first item is the
          specie, and the other items the positions and name of the
          interstitials for that specie. 
        subs : str or None
          substitution type. If None, will create an interstitial.

      :return: a 2-tuple consisting of:

        - the structure with a substitution.
        - the substituted atom in the structure above. The atom is given an
          additional attribute, C{index}, referring to list of atoms in the
          structure.
  """
  from re import compile
  from copy import deepcopy
  from numpy import dot

  single_regex    = compile(r'[A-Z][a-z]?')
  specie_regex    = compile(r'\b([A-Z][a-z]?)(\d+)\b')

  # Interstitials.
  if hasattr(type, 'rstrip'): type = type.rstrip()
  if hasattr(type, 'lstrip'): type = type.lstrip()
  if type == None or type.lower() in ['interstitial', 'interstitials', 'none']: 
    for result in interstitials(structure, lattice, subs): yield result
  # O\d: looking for specific atoms.
  elif specie_regex.match(type) != None:
    # looks for atom to modify
    found = specie_regex.match(type)
    type, index = found.group(1), int(found.group(2)) 
    if index < 1: index = 1
    for i, site in enumerate(lattice.sites):
      if type in site.type: continue
      index -= 1
      if index == 0: break
    assert index == 0, ValueError("Could not find {0}.".format(type))
    for index, atom in enumerate(structure.atoms):
      if atom.site == i: break
    assert atom.site == i, ValueError('Could not find atomic-site.')
    for result in non_interstitials(structure, index, subs): yield result
  # O, Mn ... but not O1: looking for inequivalent sites.
  elif single_regex.match(type) != None: 
    for result in inequiv_non_interstitials(structure, lattice, type, subs): yield result
  else: raise ValueError("Don't understand defect type {0}".format(type))

def charged_states(species, A, B):
  """ Loops over charged systems. 

      Charged states are given as A on B, where A and B are either None or an
      atomic-specie. If B is None, this indicates that A should be an
      interstitial. If A is None, than the defect is a vacancy.
      The charge states are as follows:
      
      - vacancy: between 0 and -B.oxidation.
      - interstitial: between 0 and A.oxidation.
      - substitution: between the maximum and minimum values of 
        A.oxidation - B.oxidation, -B.oxidation, 0.

      :Parameters:
        species : `lada.vasp.specie.Specie`
          A dictionary containing the description of the atomic species.
        A : None, or str
          If None, indicates that the charge states of an interstial are
          requested. If a string, it should be an atomic symbol with an entry
          in `species`. The defect is then either a vacancy or a substitution,
          depending on `B`. 
        B : None, or str
          If None, indicates that the charge states of a vacancy are
          requested. If a string, it should be an atomic symbol with an entry
          in `species`. The defect is then either an interstitial or a
          substitution, depending on `A`. 

      :return: Yields a 2-tuple:

        - Number of electrons to add to the system (not charge). 
        - a suggested name for the charge state calculation.
  """
  if A == 'None': A = None
  if B == 'None': B = None
  assert A != None or B != None, ValueError("Both A and B cannot be None")

  if A == None:   # vacancy! Charge states are in 0 to -B.oxidation.
    B = species[B]
    max_charge = -B.oxidation if hasattr(B, "oxidation") else 0
    min_charge = 0
  elif B == None: # interstial! Charge states are in 0, A.oxidation.
    A = species[A[0]]
    max_charge = A.oxidation if hasattr(A, "oxidation") else 0
    min_charge = 0
  else:           # substitution! Charge states are difference of A and B.
    A, B = species[A], species[B]
    Aox = A.oxidation if hasattr(A, "oxidation") else 0
    Box = B.oxidation if hasattr(B, "oxidation") else 0
    max_charge = max(Aox - Box, -Box, 0)
    min_charge = min(Aox - Box, -Box, 0)
    
  if max_charge < min_charge: max_charge, min_charge = min_charge, max_charge

  for charge in range(min_charge, max_charge+1):
    # directory
    if   charge == 0:   oxdir = "charge_neutral"
    elif charge > 0:    oxdir = "charge_" + str(charge) 
    elif charge < 0:    oxdir = "charge_" + str(charge) 
    yield -charge, oxdir


def band_filling(defect, host, **kwargs):
  """ Returns band-filling corrrection. 

      :Parameters: 

        defect 
          An output extraction object as returned by the vasp functional when
          computing the defect of interest.
        host 
          An output extraction object as returned by the vasp functional when
          computing the host matrix.
        kwargs 
          Parameters are passed on to potential alignment calculations.
         
      :return: Band-filling correction in eV.

      Accounts for Moss-Burnstein band-filling effects in the case of shallow
      donors and acceptors.
      The result of this call should be added to the energy of the point-defect.
  """
  from numpy import sum, multiply, newaxis
  from quantities import eV

  potal = potential_alignment(defect, host, **kwargs)

  cbm = host.cbm - potal
  if defect.eigenvalues.ndim == 3:
    dummy = multiply(defect.eigenvalues-cbm, defect.multiplicity[newaxis,:,newaxis])
  elif defect.eigenvalues.ndim == 2:
    dummy = multiply(defect.eigenvalues-cbm, defect.multiplicity[:, newaxis])
  dummy = multiply(dummy, defect.occupations)
  result = -sum(dummy[defect.eigenvalues > cbm])

  vbm = host.vbm - potal
  if defect.eigenvalues.ndim == 3:
    dummy = multiply(vbm-defect.eigenvalues, defect.multiplicity[newaxis,:,newaxis])
    dummy = multiply(dummy, 1e0-defect.occupations)
  elif defect.eigenvalues.ndim == 2:
    dummy = multiply(vbm-defect.eigenvalues, defect.multiplicity[:, newaxis])
    dummy = multiply(dummy, 2e0-defect.occupations)
  result -= sum(dummy[defect.eigenvalues < vbm])

  return -result.rescale(eV)
  

def potential_alignment(defect, host, maxdiff=0.5):
  """ Returns potential alignment correction. 

      :Parameters:

        defect 
          An output extraction object as returned by the vasp functional when
          computing the defect of interest.
        host 
          An output extraction object as returned by the vasp functional when
          computing the host matrix.
        maxdiff : float
          Maximum difference between the electrostatic potential of an atom and
          the equivalent host electrostatic potential beyond which that atom is
          considered pertubed by the defect.

          If this number is negative, than the algorithm is to exclude from the
          average the ``-maxdiff`` atoms furthest from the average (of the host). 

      :return: The potential alignment in eV (without charge factor).
  """
  from operator import itemgetter
  from numpy import mean, array, abs
  from quantities import eV
  from .. import specie_list

  if abs(defect.charge) < 1e-12: return 0 * eV

  # first get average atomic potential per atomic specie in host.
  host_electropot = {}
  host_species = specie_list(host.structure)
  for s in host_species:
    indices = array([atom.type == s for atom in host.structure.atoms])
    host_electropot[s] = mean(host.electropot[indices]) 
  
  # creates an initial list of unperturbed atoms. 
  types = array([atom.type for atom in defect.structure.atoms])
  unperturbed = array([(type in host_species) for type in types])

  # now finds unpertubed atoms according to maxdiff
  while any(unperturbed): # Will loop until no atom is left. That shouldn't happen!

    # computes average atomic potentials of unperturbed atoms.
    defect_electropot = {}
    for s in host_species:
      defect_electropot[s] = mean(defect.electropot[unperturbed & (types == s)])

    # finds atomic potential farthest from the average potentials computed above.
    by_type = array([host_electropot[type] for type in types[unperturbed]]) * eV
    discrepancies = defect.electropot[unperturbed] - by_type

    # if discrepancy too high, mark atom as perturbed.
    if maxdiff > 0e0: 
      index, diff = max(enumerate(abs(discrepancies)), key=itemgetter(1))
      if diff > maxdiff:
        unperturbed[ [i for i, u in enumerate(unperturbed) if u][index] ] = False
      # otherwise, we are done for this loop!
      else: break
    # mark n atoms with highest relative discrepancy.
    elif maxdiff < 0e0:
      maxdiff = int(0.01 - maxdiff)
      discrepancies = sorted(enumerate(discrepancies / by_type), key=itemgetter(1))[-maxdiff:]
      for i, v in discrepancies: unperturbed[i] = False
      break 
    # Do nothing, use all atoms to compute potential alignment.
    else: break

  # now computes potential alignment from unpertubed atoms.
  by_type = array([host_electropot[type] for type in types[unperturbed]]) * eV
  diff_from_host = defect.electropot[unperturbed]  - by_type
  return mean(diff_from_host).rescale(eV) * defect.charge
                    

def third_order_charge_correction(structure, charge = None, n = 200, epsilon = 1.0, **kwargs):
  """ Returns energy of third order charge correction. 
  
      :Parameters: 
        structure : `lada.crystal.Structure`
          Defect supercell, with cartesian positions in angstrom.
        n 
          precision. Higher better.
        charge 
          If no units are given, defaults to elementary charge. If None,
          defaults to 1 elementary charge.
        epsilon 
          Static dielectrict constant of the host. Most likely in atomic units
          (e.g. dimensionless), but I'm not sure.
      
      Taken as is from `Lany and Zunger, PRB 78, 235104 (2008)`__.
      Always outputs as eV. Not sure what the units of some of these quantities are. 

      .. __:  http://dx.doi.org/10.1103/PhysRevB.78.235104

      :return: third order correction  to the energy in eV. Should be *added* to total energy.
  """
  from .._crystal import third_order 
  from numpy import array
  from quantities import elementary_charge, eV, pi, angstrom, dimensionless
  from ...physics import a0, Ry

  if charge == None: charge = 1e0
  elif charge == 0: return 0e0 * eV
  if hasattr(charge, "units"):  charge  = float(charge.rescale(elementary_charge))
  if hasattr(epsilon, "units"): epsilon = float(epsilon.simplified)
  cell = structure.cell
  scale = structure.scale
  return - charge**2 / epsilon * third_order(cell, n/2) * (4e0*pi/3e0) \
         * (array(a0.rescale(angstrom))/scale)**3 * Ry.rescale(eV)


def first_order_charge_correction(structure, charge=None, epsilon=1e0, cutoff=15e1, **kwargs):
  """ First order charge correction of +1 charge in given supercell. 
  
      Units in this function are either handled by the module Quantities, or
      defaults to Angstroems and elementary charges.

      :Parameters:
        structure : `lada.crystal.Structure`
          Defect supercell, with cartesian positions in angstrom.
        charge 
          Charge of the point-defect. Defaults to 1e0 elementary charge. If no
          units are attached, expects units of elementary charges.
        epsilon 
          dimensionless relative permittivity.
        cutoff 
          Ewald cutoff parameter.

      :return: Electrostatic energy in eV.
  """
  from numpy.linalg import norm
  from quantities import elementary_charge, eV
  from .. import Structure
  from ...physics import Ry
  try: from ...pcm import Clj 
  except ImportError as e:
    print "Could not import Point-Charge Model package (pcm). \n"\
          "Cannot compute first order charge correction.\n"\
          "Please compile LaDa with pcm enabled.\n"
    raise

  if charge == None: charge = 1
  elif charge == 0: return 0e0 * eV
  if hasattr(charge, "units"): charge = float(charge.rescale(elementary_charge))

  clj = Clj()
  clj.charges["A"] = charge
  clj.ewald_cutoff = cutoff * Ry

  struc = Structure()
  struc.cell = structure.cell
  struc.scale = structure.scale
  struc.add_atom = ((0e0,0,0), "A")

  result = clj.ewald(struc).energy / epsilon
  return -result * eV

def charge_corrections(structure, **kwargs):
  """ Electrostatic charge correction (first and third order). """
  return   first_order_charge_correction(structure, **kwargs) \
         + third_order_charge_correction(structure, **kwargs) \

def magnetic_neighborhood(structure, defect, species):
   """ Finds magnetic neighberhood of a defect. 
   
       If the defect is a substitution with a magnetic atom, then the
       neighberhood is the defect alone. Otherwise, the neighberhood extends to
       magnetic first neighbors. An atomic specie is deemed magnetic if marked
       as such in `species`.

       :Parameters: 
         structure : `lada.crystal.Structure`
           The structure with the point-defect already incorporated.
         defect : `lada.crystal.Atom`
           The point defect, to which and *index* attribute is given denoting
           the index of the atom in the original supercell structure (without
           point-defect).
         species : dict of `lada.vasp.species.Specie`
           A dictionary defining the atomic species.

       :return: indices of the neighboring atoms in the point-defect `structure`.
   """
   from numpy import array
   from numpy.linalg import norm
   from . import Neighbors

   # checks if substitution with a magnetic defect.
   if hasattr(defect, "index") and defect.index < len(structure.atoms):
     atom = structure.atoms[defect.index]
     if species[atom.type].magnetic and norm(defect.pos - atom.pos) < 1e-12:
       return [defect.index]
   # now finds first neighbors. 12 is the highest coordination number, so
   # this should include the first shell.
   neighbors = [n for n in Neighbors(structure, 12, defect.pos)]
   # only take the first shell and keep indices (to atom in structure) only.
   neighbors = [n.index for n in neighbors if n.distance < neighbors[0].distance + 1e-1]
   # only keep the magnetic neighborhood.
   return [n for n in neighbors if species[structure.atoms[n].type].magnetic]

def equiv_bins(n, N):
  """ Generator over ways to fill N equivalent bins with n equivalent balls. """
  from itertools import chain
  from numpy import array
  assert N > 0
  if N == 1: yield [n]; return
  if n == 0: yield [0 for x in range(N)]
  for u in xrange(n, 0, -1):
    for f in  equiv_bins(n-u, N-1):
      result = array([x for x in chain([u], f)])
      if all(result[0:-1]-result[1:] >= 0): yield result

def inequiv_bins(n, N):
  """ Generator over ways to fill N inequivalent bins with n equivalent balls. """
  from itertools import permutations
  for u in equiv_bins(n, N):
    u = [v for v in u]
    history = []
    for perm in permutations(u, len(u)):
      seen = False
      for other in history:
        same = not any( p != o for p, o in zip(perm, other) )
        if same: seen = True; break
      if not seen: history.append(perm); yield [x for x in perm]

def electron_bins(n, atomic_types):
  """ Loops over electron bins. """
  from itertools import product 
  from numpy import zeros, array
  # creates a dictionary where every type is associated with a list of indices into atomic_types.
  Ns = {}
  for type in set(atomic_types):
    Ns[type] = [i for i,u in enumerate(atomic_types) if u == type]
  # Distributes electrons over inequivalent atomic types.
  for over_types in inequiv_bins(n, len(Ns.keys())):
    # now distributes electrons over each type independently.
    iterables = [ equiv_bins(v, len(Ns[type])) for v, type in zip(over_types, Ns.keys()) ] 
    for nelecs in product(*iterables):
      # creates a vector where indices run as in atomic_types argument.
      result = zeros((len(atomic_types),), dtype="float64")
      for v, (type, value) in zip(nelecs, Ns.items()): result[value] = array(v)
      yield result

def magmom(indices, moments, nbatoms):
  """ Yields a magmom string from knowledge of which moments are non-zero. """
  from operator import itemgetter
  s = [0 for i in range(nbatoms)]
  for i, m in zip(indices, moments): s[i] = m
  compact = [[1, s[0]]]
  for m in s[1:]:
    if abs(compact[-1][1] - m) < 1e-12: compact[-1][0] += 1
    else: compact.append( [1, m] )
    
  string = ""
  for n, m in compact:
    if n > 1: string +=" {0}*{1}".format(n, m)
    elif n == 1: string += " {0}".format(m)
    assert n != 0
  return string

def electron_counting(structure, defect, species, extrae):
  """ Enumerate though number of electron in point-defect magnetic neighberhood. 

      Generator over the number of electrons of each atom in the magnetic
      neighberhood of a point defect with `extrae` electrons. If there are no
      magnetic neighborhood, then `magmom` is set
      to None and the total magnetic moment to 0 (e.g. lets VASP figure it out).
      Performs a sanity check on integers to make sure things are correct.

      :Parameters:
        structure : `lada.crystal.Structure`
          Structure with point-defect already inserted.
        defect : `lada.crystal.Atom`
          Atom making up the point-defect.
          In addition, it should have an *index* attribute denoting the defect 
        species : dict of `lada.vasp.species.Specie`
          Dictionary containing details of the atomic species.
        extrae
          Number of extra electrons to add/remove.

      :return: yields (indices, electrons) where indices is a list of indices
        to the atom in the neighberhood, and electrons is a corresponding list of
        elctrons.
  """
  from numpy import array
  from ...physics import Z
  indices = magnetic_neighborhood(structure, defect, species)

  # no magnetic neighborhood.
  if len(indices) == 0: 
    yield None, None
    return

  # has magnetic neighberhood from here on.
  atoms = [structure.atoms[i] for i in indices]
  types = [a.type for a in atoms]
  nelecs = array([species[type].valence - species[type].oxidation for type in types])

  # loop over different electron distributions.
  for tote in electron_bins(abs(extrae), types):
    # total number of electrons on each atom.
    if extrae < 0:   tote = nelecs - tote
    elif extrae > 0: tote += nelecs

    # sanity check. There may be more electrons than orbitals at this point.
    sane = True
    for n, type in zip(tote, types):
      if n < 0: sane = False; break;
      z = Z(type)
      if (z >= 21 and z <= 30) or (z >= 39 and z <= 48) or (z >= 57 and z <= 80):  
        if n > 10: sane = False;  break
      elif n > 8: sane = False; break

    if not sane: continue

    yield indices, tote
 

def low_spin_states(structure, defect, species, extrae, do_integer=True, do_average=True):
  """ Enumerate though low-spin-states in point-defect. 

      Generator over low-spin magnetic states of a defect with
      `extrae` electrons. The extra electrons are distributed both as integers
      and as an average. All these states are ferromagnetic. In the special
      case of a substitution with a magnetic atom, the moment is expected to go
      on the substitution alone. If there are no magnetic neighborhood, then `magmom` is set
      to None and the total magnetic moment to 0 (e.g. lets VASP figure it out).
      
      :Parameters:
        structure : `lada.crystal.Structure`
          Structure with point-defect already inserted.
        defect : `lada.crystal.Atom`
          Atom making up the point-defect.
          In addition, it should have an *index* attribute denoting the defect 
        species : dict of `lada.vasp.species.Specie`
          Dictionary containing details of the atomic species.
        extrae
          Number of extra electrons to add/remove.

      :return: yields (indices, moments) where the former index the relevant
               atoms in `structure` and latter are their respective moments.
  """

  from numpy import array, abs, all, any

  history = []
  def check_history(*args):
    for i, t in history:
      if all(abs(i-args[0]) < 1e-12) and all(abs(t-args[1]) < 1e-12):
        return False
    history.append(args)
    return True


  if do_integer: 
    for indices, tote in electron_counting(structure, defect, species, extrae):
      if tote == None: continue # non-magnetic case
      indices, moments = array(indices), array(tote) % 2
      if all(abs(moments) < 1e-12): continue # non - magnetic case
      if check_history(indices, moments): yield indices, moments
  if do_average: 
    for indices, tote in electron_counting(structure, defect, species, 0):
      if tote == None: continue # non-magnetic case
      if len(indices) < 2: continue
      indices, moments = array(indices), array(tote) % 2 + extrae / float(len(tote))
      if all(abs(moments) < 1e-12): continue # non - magnetic case
      if check_history(indices, moments): yield indices, moments


def high_spin_states(structure, defect, species, extrae, do_integer=True, do_average=True):
  """ Enumerate though high-spin-states in point-defect. 

      Generator over high-spin magnetic states of a defect with
      `extrae` electrons. The extra electrons are distributed both as integers
      and as an average. All these states are ferromagnetic. In the special
      case of a substitution with a magnetic atom, the moment is expected to go
      n the substitution alone. If there are no magnetic neighborhood, then
      `magmom` is set to None and the total magnetic moment to 0 (e.g. lets
      VASP figure it out).

      :Parameters:
        structure : `lada.crystal.Structure`
          Structure with point-defect already inserted.
        defect : `lada.crystal.Atom`
          Atom making up the point-defect.
          In addition, it should have an *index* attribute denoting the defect 
        species : dict of `lada.vasp.species.Specie`
          Dictionary containing details of the atomic species.
        extrae 
          Number of extra electrons to add/remove.

      :return: yields (indices, moments) where the former index the relevant
               atoms in `structure` and latter are their respective moments.
  """
  from numpy import array, abs, all, any

  def is_d(t): 
    """ Determines whether an atomic specie is transition metal. """
    from ...physics import Z
    z = Z(t)
    return (z >= 21 and z <= 30) or (z >= 39 and z <= 48) or (z >= 57 and z <= 80) 

  def determine_moments(arg, ts): 
    """ Computes spin state from number of electrons. """
    f = lambda n, t: (n if n < 6 else 10-n) if is_d(t) else (n if n < 5 else 8-n)
    return [f(n,t) for n, t in zip(arg, ts)]

  history = []
  def check_history(*args):
    for i, t in history:
      if all(abs(i-args[0]) < 1e-12) and all(abs(t-args[1]) < 1e-12):
        return False
    history.append(args)
    return True
  
  if do_integer: 
    for indices, tote in electron_counting(structure, defect, species, extrae):
      if tote == None: continue # non-magnetic case
      
      types = [structure.atoms[i].type for i in indices]
      indices, moments = array(indices), array(determine_moments(tote, types))
      if all(moments == 0): continue # non - magnetic case
      if check_history(indices, moments):  yield indices, moments

  if do_average: 
    for indices, tote in electron_counting(structure, defect, species, 0):
      if tote == None: continue # non-magnetic case
      if len(indices) < 2: continue

      types = [structure.atoms[i].type for i in indices]
      indices = array(indices)
      moments = array(determine_moments(tote, types)) + float(extrae) / float(len(types))
      if all(abs(moments) < 1e-12): continue # non - magnetic case
      if check_history(indices, moments):  yield indices, moments

def magname(moments, prefix=None, suffix=None):
  """ Construct name for magnetic moments. """
  if len(moments) == 0: return "paramagnetic"
  string = str(moments[0])
  for m in moments[1:]: string += "_" + str(m)
  if prefix != None: string = prefix + "_" + string
  if suffix != None: string += "_" + suffix
  return string
