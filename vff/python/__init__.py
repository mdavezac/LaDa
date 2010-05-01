""" Valence Force Field functional for zinc-blende materials """

# imports C++ extension
import vff

def zinc_blend_lattice():
  """ Defines a default zinc-blende lattice (InGaNP). """
  from ..crystal import Lattice, Site
  from numpy import array

  lattice = Lattice()
  lattice.cell = array( [[0,0.5,0.5], [0.5,0,0.5], [0.5,0.5,0]], dtype="float64")
  lattice.scale = 1e0
  lattice.sites.append( Site(array([0,0,0]), ["None"]) )
  lattice.sites.append( Site(array([0.25,0.25,0.25]), ["None"]) )
  return lattice

# class Vff(Object): 
#   """ Valence Force Field functional for zinc-blende materials """
#   def __init__(self):
#     """ Initializes a valence force field functional. """
#     from ..minimizer import Minimizer

#     super(Minimizer, self).__init__()
#     self.lattice = zinc_blende_lattice()
#     """ Lattice for which to perform calculations.
#     
#         In practice, zinc-blende lattices can be defined for any number of parameterizations.
#         The one chosen here is such that atomic site 0 is at the origin of the
#         unit cell, and atomic site 1 at 0.25, 0.25, 0.25. The scale is
#         arbitrary at this point, and set to 4.5. The unit cell is [[0,0.5,0.5],
#         [0.5,0,0.5], [0.5,0.5,0]]. 
#     """

#     self.relax = True
#     """ Whether to relax strain or just get energy. """
#     self.direction = None
#     """ Epitaxial direction if any.
#     
#         If set, it should be a vector. In that case, relaxation will only
#         happen in that direction. Otherwise, everything will be relaxed. 
#     """
#     self.minimizer = Minimizer()
#     """ Minimizer to use with VFF """
#     self.bonds = {}
#     """ Bond parameters. """
#     self.angles = {}
#     """ Angle parameters. """
#     self.OUTCAR = "vffout" 
#     """ File where to redirect output. """
#     self.ERRCAR = "vfferr" 
#     """ File where to redirect errors. """

#   def _not_available(self): raise RuntimeError("Property cannot be gotten.")

#   def _add_bond(self, args):
#     """ Adds/Modifies the bond parameter dictionary.
#     
#         @param args: - the first element of arg is one endpoint of the bond (say "In"), 
#                      - the second element is another endpoint (say "N"),
#                      - the last element is a sequence with at most 5 elements
#                        defining the bond-stretching parameters (order 2 through 6).
#         @type args: "type1", "type2", numpy array
#     """
#     assert len(args) == 3, RuntimeError("Bonds should be set with 3 parameters: %s." % (args))
#     assert args[0] in [ u for u in site.type for site in self.lattice.sites ],\
#            RuntimeError( "Type %s is not in lattice." % (args[0]))
#     assert args[1] in [ u for u in site.type for site in self.lattice.sites ],\
#            RuntimeError( "Type %s is not in lattice." % (args[1]))

#     name = "%s-%s" % (args[0], args[1])
#     if args[0] > args[1]: name = "%s-%s" % (args[1], args[0])

#     params = [u for u in args]
#     assert len(params) <= 5 and len(params) > 0, RuntimeError("To few or too many parameters: %s." % (params))
#     if name in self.bond_stretching: # replaces none value with known value
#       for old, new in zip(self.bond_stretching[name], params):
#         if old == None: old = new
#     if len(params) < 5: params.extend( 0 for u in range(5-len(params)) )
#     
#     self.bonds[name] = params

#   def set_bond(self, A, B, d0 = None, a2 = None, a3 = None, a4 = None, a5 = None, a6 = None):
#     """ Adds/Modifies the bond parameter dictionary.
#     
#         @param A: Endpoint specie.
#         @type A: str
#         @param B: Middle specie.
#         @type B: str
#         @param C: Endpoint specie.
#         @type C: str
#         @param gamma: S{gamma} parameter. Can be a number or a string starting with "tet".
#         @param sigma: S{sigma} parameter.
#         @param a2: 2nd order bond-angle parameter.
#         @param a3: 3rd order bond-angle parameter.
#         @param a4: 4rd order bond-angle parameter.
#         @param a5: 5rd order bond-angle parameter.
#         @param a6: 6rd order bond-angle parameter.
#     """
#     self._add_bond( (A,B, d0, a2, a3, a4, a5, a6) )

#   add_bond = property(_not_available, _add_bond, _add_bond.__doc__)



#   def _add_angle(self, args):
#     """ Adds/Modifies the angle parameter dictionary.
#     
#         @param args: - the first element of arg is one endpoint of the angle (say "In"), 
#                      - the second element is middle of the angle (say "N"),
#                      - the third element is the other endpoint (say "Ga"),
#                      - the last element is a sequence with at most 7 elements
#                        defining the bond-angle parameters. The first is the
#                        gamma parameter, the second the sigma, and the others
#                        the bond bending proper.
#         @type args: "type1", "type2", numpy array
#     """
#     assert len(args) == 4, RuntimeError("Angle should be set with 4 parameters: %s." % (args))
#     assert args[0] in [ u for u in site.type for site in self.lattice.sites ],\
#            RuntimeError( "Type %s is not in lattice." % (args[0]))
#     assert args[1] in [ u for u in site.type for site in self.lattice.sites ],\
#            RuntimeError( "Type %s is not in lattice." % (args[1]))
#     assert args[2] in [ u for u in site.type for site in self.lattice.sites ],\
#            RuntimeError( "Type %s is not in lattice." % (args[2]))

#     name = "%s-%s-%s" % (args[0], args[1], args[2])
#     if args[0] > args[2]: name = "%s-%s" % (args[2], args[1], args[0])

#     params = tuple([u for u in args])
#     assert len(params) <= 7 and len(params) > 0, RuntimeError("To few or too many parameters: %s." % (params))
#     if name in self.bond_bending: # replaces none value with known value
#       for old, new in zip(self.bond_bending[name], params):
#         if old == None: old = new
#     if str(params[0]).lower()[:3] == "tet": params[0] = -1e0/3e0
#     if len(params) < 7: params.extend( 0e0 for u in range(7-len(params)) )

#     self.angles[name] = params
#     
#   def set_angle(self, A, B, gamma = None, sigma = None, a2 = None, a3 = None, a4 = None, a5 = None, a6 = None):
#     """ Adds/Modifies the angle parameter dictionary.
#     
#         @param A: Endpoint specie.
#         @type A: str
#         @param B: Middle specie.
#         @type B: str
#         @param C: Endpoint specie.
#         @type C: str
#         @param gamma: S{gamma} parameter. Can be a number or a string starting with "tet".
#         @param sigma: S{sigma} parameter.
#         @param a2: 2nd order bond-angle parameter.
#         @param a3: 3rd order bond-angle parameter.
#         @param a4: 4rd order bond-angle parameter.
#         @param a5: 5rd order bond-angle parameter.
#         @param a6: 6rd order bond-angle parameter.
#     """
#     self._add_angle( (A,B, C, gamma, sigma, a2, a3, a4, a5, a6) )

#   add_angle = property(_not_available, _add_angle, _add_angle.__doc__)

#   def __str__(self):
#     result  = "# Lattice definition: \n%s\n" % (self.lattice)
#     result += "# VFF definition:\n vff = Vff()\n vff.lattice = lattice"
#     result += "vff.lattice = lattice\n"
#     result += "vff.OUTCAR = %s\n" % (self.OUTCAR)
#     result += "vff.ERRCAR = %s\n" % (self.ERRCAR)
#     result += "vff.direction = %s\n" % (self.direction)
#     result += "vff.relax = %s\n" % (self.relax)
#     for name, params in self.bonds.items():
#       A, B = name.split("-")
#       result += "vff.add_bond = %s, %s, %s" % (A, B)
#       for u in params: result += ", %s" % (u)
#       result += "\n"
#     for name, params in self.angles.items():
#       A, B, C = name.split("-")
#       result += "vff.add_angle = %s, %s, %s" % (A, B, C)
#       for u in params: result += ", %s" % (u)
#       result += "\n"
#     return result


#   def __call__(self, structure, outdir, comm = None ):
#     """ Performs calculation """
#     from copy import deepcopy
#     from os.path import exists, isdir
#     
#     # make this functor stateless.
#     this      = deepcopy(self)
#     outdir    = deepcopy(outdir)

#     # if other keyword arguments are present, then they are assumed to be
#     # attributes of self, with value their expected value before launch. 
#     if len(kwargs) != 0: 
#       for key in kwargs.keys(): getattr(this, key).value = kwargs[key]

#     # First checks if directory outdir exists (and is a directory).
#     if exists(outdir):
#       if not isdir(outdir): raise IOError, "%s exists but is not a directory.\n" % (outdir)
#       # checks if it contains a successful run.
#       extract = Extract(comm = comm, directory = outdir)
#       if extract.success: return extract # in which case, returns extraction object.
#     
#     # Otherwise, performs calculation by calling base class functor.
#     
#     # checks if result was successful
#     extract = Extract(comm = comm, directory = outdir)
#     if not extract.success:
#       raise RuntimeError, "VASP calculation did not complete in %s.\n" % (outdir)

#     return extract
