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
    name = "vacancy_{0}".format(type)
    if len(inequivs) > 1: name += "/site_{0}".format(i)
    # creates vacancy and keeps atom for record
    atom = deepcopy(structure.atoms.pop(which))
    atom.index = i
    # structure 
    result = deepcopy(structure)
    result.name = name
    # returns structure with vacancy.
    yield result, atom
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
          type of atoms for which to create substitution.
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

  # Case for vacancies.
  if subs == None: 
    for args in vacancy(structure, lattice, type):
      yield args
    return

  
  result = deepcopy(structure)
  inequivs = inequivalent_sites(lattice, type)
  for i in inequivs:

    # finds first qualifying atom
    for which, atom in enumerate(structure.atoms):
      if atom.site == i: break

    assert which < len(structure.atoms), RuntimeError("Site index not found.")

    # name of this substitution
    name = "{0}_on_{1}".format(subs, type) 
    if len(inequivs) > 1: name += "/site_{0}".format(i)
    # creates substitution
    orig = result.atoms[which].type
    result.atoms[which].type = subs
    # substituted atom.
    substituted = deepcopy(result.atoms[which])
    substituted.index = which
    result.name = name
    # returns structure with vacancy.
    yield result, substituted
    # removes substitution
    result.atoms[which].type = orig

def charged_states(species, A=None, B=None):
  """ Loops over charged systems. 

      :Types:
        - species: `lada.vasp.specie.Specie`
        - `A` : None or str index to `species`.
        - `B` : None or str index to `species`.

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
  else: A = species[A]
  if B == None:
    # max/min oxidation state
    max_charge = -A.oxidation if hasattr(A, "oxidation") else 0
    min_charge = 0
    if max_charge < min_charge: max_charge, min_charge = min_charge, max_charge
  else:
    B = species[B]
    # Finds maximum range of oxidation.
    maxA_oxidation = A.oxidation if hasattr(A, "oxidation") else 0
    maxB_oxidation = B.oxidation if hasattr(B, "oxidation") else 0
    max_charge = max(-maxA_oxidation, maxB_oxidation - maxA_oxidation, 0)
    min_charge = min(-maxA_oxidation, maxB_oxidation - maxA_oxidation, 0)

  for charge in range(min_charge, max_charge+1):
    # directory
    if   charge == 0:   oxdir = "charge_neutral"
    elif charge > 0:    oxdir = "charge_" + str(charge) 
    elif charge < 0:    oxdir = "charge_" + str(charge) 
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
  diff_from_host =    defect.electropot[unpertubed]  \
                    - array([host_electropot[type] for type in types[unperturbed]])
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


def first_order_charge_correction(cell, charge=1e0, cutoff=None):
  """ First order charge correction of +1 charge in given supercell.
  
      :Parameters:
        - `cell`: Supercell of the point-defect.
        - `charge`: Charge of the point-defect.
        - `cutoff`: Ewald cutoff parameter.
  """
  from numpy.linalg import norm
  from ..crystal import Structure
  try: from ..pcm import Clj 
  except ImportError as e:
    print "Could not import Point-Charge Model package (pcm). \n"\
          "Cannot compute first order charge correction.\n"\
          "Please compile LaDa with pcm enabled.\n"
    raise

  clj = Clj()
  clj.charge["A"] = charge
  if cutoff == None:
    clj.ewald = 10 * max( [norm(cell[:,i]) for i in range(3)] )
  else: clj.ewald_cutoff = cutoff

  structure = Structure(cell)
  structure.add_atom = ((0e0,0,0), "A")

  return clj.ewald(structure)

# only defined if vasp can be imported.
def first_neighbors(structure, origin):
   """ Finds fist-neighbors of origin. """
   from . import Neighbors

   # now finds first neighbors. 12 is the highest coordination number, so
   # this should include the first shell.
   neighbors = [n for n in Neighbors(structure, 12, substitution.pos)]
   # only take the first shell and keep indices (to atom in structure) only.
   neighbors = [n.index for n in neighbors if n.distance < neighbors[0].distance + 1e-1]
   return neighbors

def match_moment(structure, origin, species, moment, extrae):

  from ..physics import Z
  indices = first_neighbors(structure, origin)
  indices = filter(indices, lambda x: species[ structure.atoms[x] ].magnetic)
  atoms = [structure.atoms[i] for i in indices]
  types = [a.type for a in atoms]
  oxidations = [species[a].oxidation for a in types]

  maps = {} 
  for type in set(types):
    maps[type] = [(a, n, o) for a, n, t, o in zip(atoms, indices, types, oxidations) if t == type]
  nelecs = [species[type].valence - species[type].oxidation for type in types]

  for type in set(types):
    nelecs[type] = species[type].valence - species[type].oxidation
  
  def equiv_bins(n, N):
    """ Generator over ways to fill N equivalent bins with n equivalent balls. """
    from itertools import chain
    if N == 1: yield [n]; return
    if n == 0: yield [0 for x in range(N)]
    for u in xrange(n, 0, -1):
      for f in  equiv_bins(n-u, N-1):
        yield chain([u], f)

  def inequiv_bins(n, N):
  """ Generator over ways to fill N inequivalent bins with n equivalent balls. """
  from itertools import permutations
  for u in equiv_bins(n, N):
    u = [v for v in u]
    history = []
    for perm in permutations(u, len(u)):
      seen = False
      for other in history:
        same = False
        for p, o in zip(perm, other):
          if p != o: same = True; break
        if same == False: seen = True; break
      if not seen: history.append([x for x in perm]); yield perm

    nb_electrons = species[type].valence
    z = Z(type)
    if (z >= 21 and z <= 30) or (z >= 39 and z <= 48) or (z >= 57 and z <= 80):  
      nelecs[type] = species[type].valence - species[type].oxidation
    else:
      nelecs[type] = species[type].valence - species[type].oxidation
      nelecs = nelecs, nelecs + extrae
      nelecs = min(nelecs), max(nelecs)

      hs = [0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0]
      ls = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
      moments[type] = hs







try: from .. import vasp as _vasp
except ImportError as e: pass
else: 
  class Magmom(_vasp.incar._params.SpecialVaspParam):
    """ Creates a magmom configuration, whether low-spin or high-spin. """
    def __init__(self, value, config = None, indices = None):
      self._config, self._indices = config, indices
      super(Magmom, self).__init__(value)
  
    def _get_value(self):
      if self._config == None: return None
      if self._indices == None: return None
      elif len(self._indices) == 0: return None
      return self._config, self._indices
  
    def _set_value(self, value): 
      if value == None: self._indices, self._config = None, None
      elif isinstance(value, str):
        if value.lower() == "high":  self._config = value.lower()
        elif value.lower() == "low": self._config = value.lower()
        else: raise ValueError("Unkown value for magmom: " + str(value) + ".")
      elif hasattr(value, "__len__"): # sequence.
        if isinstance(value[0], str):
          self._config = value[0]
          value = value[1]
        if len(value) > 0: self._indices = sorted(u for u in value)
      else: raise ValueError("Unkown value for magmom: " + str(value) + ".")
    value = property( _get_value, _set_value, \
                      doc = """ ("low|high", [indices]") """ )
    
    @broadcast_result(key=True)
    def incar_string(self, vasp, *args, **kwargs):
      from ..crystal import specie_list
  
      magmom = ""
      all_types = [atom.type for atom in vasp._system.atoms]
      for specie in specie_list(vasp._system): # do this per specie.
        indices = [n for n in self._indices if vasp._system.atoms[n].type == specie]
        enum_indices = [i for i, n in enumerate(self._indices) if vasp._system.atoms[n].type == specie]
        if len(indices) == 0: # case where there are no magnetic species of this kind
          magmom += "%i*0 " % (all_types.count(specie))
          continue
        species = [vasp.species[ vasp._system.atoms[n].type ] for n in indices]
        extra_elecs = -vasp.nelect if vasp.nelect != None else 0
    
        # computes low or high spin configurations. 
        # Assumes that magnetic species are s^2 d^n p^0!
        per_d = extra_elecs / len(self._indices) # number of electrons added by oxidation per atom
        leftovers = extra_elecs % len(self._indices) # special occupations...
        d_elecs = [int(s.valence-2+0.5) for s in species] # d occupation.
        d_elecs = [  s + per_d + (0 if i < leftovers else 1)\
                     for i, s in zip(enum_indices, d_elecs) ] # adds extra e-/h+
        d_elecs = [0 if s > 10 else s for s in d_elecs] # checks for d^n with n > 10
        mag = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0] if self._config == "low" \
              else [0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0]
        d_elecs = [mag[d] for d in d_elecs] # d occupation.
    
        # now constructs magmom
        last_index = -1
        for index, valence in zip(indices, d_elecs):
          index = all_types[:index].count(specie)
          if index - last_index == 1:   magmom += "%f " % (valence)
          elif index - last_index == 2: magmom += "0 %f " % (valence)
          else:                         magmom += "%i*0 %f " % (index-last_index-1, valence)
          last_index = index
        # adds last memebers of specie.
        index = all_types.count(specie)
        if index - last_index == 1: pass
        elif index - last_index == 2: magmom += "0 "
        else:                         magmom += "%i*0 " % (index - last_index - 1)
    
      return "MAGMOM = %s" % (magmom)
  
    def __repr__(self):
      """ Returns a python script representing this object. """
      return "vasp = %s(%s, config=%s, indices=%s))\n" \
             % (self.__classname__, repr(self.value), repr(self._config), repr(self.indices))

    
  
          
