""" Methods to write structures from file. """
def poscar(structure, file='POSCAR', vasp5=None, substitute=None):
  """ Writes a poscar to file. 

      :param structure:
          The structure to print out.
      :type structure:
          :py:class:`Structure`
      :param file:
          Object with a ''write'' method. If a string, then a file object is
          opened with that filename. The file is overwritten. If None, then
          writes to POSCAR in current working directory.
      :type file: str, stream, or None.
      :param bool vasp5:
          If true, include species in poscar, vasp-5 style.  Otherwise, looks
          for :py:data:`is_vasp_4 <pylada.is_vasp_4>` global config variable. 
          Defaults to False, in which case, does not print specie types.
      :param substitute:
          If present, will substitute the atom type in the structure. Can be
          incomplete. Only works with vasp5 = True (or :py:data:`is_vasp_4
          <pylada.is_vasp_4>` = True).
      :type substitute:
          dict or None
  
      >>> with open("POSCAR", "w") as file: write.poscar(structure, file, vasp5=True)

      Species in structures can be substituted for others (when using vasp5 format).
      Below, aluminum atoms are replaced by cadmium atoms. Other atoms are left unchanged.

      >>> with open("POSCAR", "w") as file:
      >>>   write.poscar(structure, file, vasp5=True, substitute={"Al":"Cd"})

      Selective dynamics are added to the POSCAR file if an atom in the
      structure has a freeze attribute (of non-zero length). It is expected
      that this attribute is a string and that it contains one of "x", "y",
      "z", corresponding to freezing the first, second, or third fractional
      coordinates. Combinations of these are also allowed.
  """
  from quantities import angstrom
  if file is None:
    with open('POSCAR', 'w') as fileobj: return poscar(structure, fileobj, vasp5, substitute)
  elif not hasattr(file, 'write'):
    with open(file, 'w') as fileobj: return poscar(structure, fileobj, vasp5, substitute)

  from numpy import matrix, dot
  from . import specieset

  if vasp5 is None:
    import pylada 
    vasp5 = not getattr(pylada, 'is_vasp_4', True)

  string = "{0}\n{1}\n".format(getattr(structure, 'name', ''),
                               float(structure.scale.rescale(angstrom)))
  for i in range(3):
    string += "  {0[0]} {0[1]} {0[2]}\n".format(structure.cell[:,i])
  species = specieset(structure)
  if vasp5: 
    if substitute is None: substitute = {}
    for s in species: string += " {0} ".format(substitute.get(s,s))
    string += "\n"
  for s in species: 
    string += "{0} ".format(len([0 for atom in structure if atom.type == s]))
  inv_cell = matrix(structure.cell).I
  selective_dynamics =\
      any([len(getattr(atom, 'freeze', '')) != 0 for atom in structure])
  if selective_dynamics: string += "\nselective dynamics\ndirect\n"
  else: string += '\ndirect\n'
  for s in species: 
    for atom in structure:
      if atom.type != s: continue
      string += "  {0[0]} {0[1]} {0[2]}"\
                .format(dot(inv_cell, atom.pos).tolist()[0])
      freeze = getattr(atom, 'freeze', '')
      if selective_dynamics:
        string += "  {1} {2} {3}\n"\
                    .format( 'T' if 'x' in freeze  != 0 else 'F', 
                             'T' if 'y' in freeze  != 0 else 'F', 
                             'T' if 'z' in freeze  != 0 else 'F' ) 
      else: string += '\n'
  if file == None: return string
  elif isinstance(file, str): 
    from ..misc import RelativePath
    with open(RelativePath(file).path, 'w') as file: file.write(string)
  else: file.write(string)

def castep(structure, file=None):
  """ Writes castep input. """
  from quantities import angstrom
  cell = structure.cell * float(structure.scale.rescale(angstrom))
  string = "%BLOCK LATTICE_CART\n" \
           "  {0[0]} {0[1]} {0[2]}\n" \
           "  {1[0]} {1[1]} {1[2]}\n" \
           "  {2[0]} {2[1]} {2[2]}\n" \
           "%ENDBLOCK LATTICE_CART\n\n"\
           "%BLOCK POSITIONS_ABS\n".format(*(cell.T))
  for atom in structure:
    pos = atom.pos * float(structure.scale.rescale(angstrom))
    string += "  {0} {1[0]} {1[1]} {1[2]}\n"\
              .format(atom.type, pos)
  string += "%ENDBLOCK POSITION_ABS\n"
  if file == None: return string
  elif isinstance(file, str): 
    from ..misc import RelativePath
    with open(RelativePath(file).path, 'w') as file: file.write(string)
  else: file.write(string)

def crystal( structure, file='fort.34',
             dimensionality=None, centering=None, type=None, spacegroup=None ):
  """ Writes structure as CRYSTAL's EXTPRT. 

      :param structure:
          The structure to print out.
      :type structure: :py:class:`Structure`
      :param file:
          Object with a ''write'' method. If a string, then a file object is
          opened with that filename. The file is overwritten. If None, then
          returns a string.
      :type file: str, stream, or None.
      :param int dimensionality:
          Dimensionality of the system as an integer between 0 and  3 included.
          If None, checks for a ''dimensionality'' attribute in ``structure``.
          If that does not fit the biil, defaults to 3.
      :param int centering:
          Centering in CRYSTAL_'s integer format. Is None, checks ``structure``
          for a ``centering`` integer attribute. If that does not exist or is
          not convertible to an integer, then defaults to 1.
      :param int type:
          Crystal type in CRYSTAL_'s integer format. Is None, checks
          ``structure`` for a ``type`` integer attribute. If that does not
          exist or is not convertible to an integer, then defaults to 1.
      :param spacegroup:
          The structure's space group as a sequence of 4x3 matrices. If this is
          None (default), then checks for ''spacegroup'' attributes. If that
          does not exist, uses :py:function:`~pylada.crystal.space_group`.
  """
  from StringIO import StringIO
  from numpy import zeros
  from quantities import angstrom
  from ..crystal.iterator import equivalence as equivalence_iterator
  from ..periodic_table import find as find_specie
  from . import space_group
  # makes sure file is a stream.
  # return string when necessary
  if file is None:
    file = StringIO()
    crystal(structure, file, dimensionality, centering, type, spacegroup)
    return file.getvalue()
  elif not hasattr(file, 'write'):
    with open(file, 'w') as fileobj:
      return crystal( structure, fileobj, dimensionality,
                      centering, type, spacegroup)
  
  # normalize input as keyword vs from structure vs default.
  try:
    if dimensionality is None:
      dimensionality = getattr(structure, 'dimensionality', 3)
    dimensionality = int(dimensionality)
    if dimensionality < 0 or dimensionality > 3: dimensionality = 3
  except: dimensionality = 3
  try:
    if centering is None: centering = getattr(structure, 'centering', 1)
    centering = int(centering)
  except: centering = 1
  try:
    if type is None: type = getattr(structure, 'type', 1)
    type = int(type)
  except: type = 1
  if spacegroup is None: spacegroup = getattr(structure, 'spacegroup', None)
  if spacegroup is None: spacegroup = space_group(structure)
  if len(spacegroup) == 0:
    spacegroup = zeros((1, 4, 3))
    spacegroup[0,0,0] = 1
    spacegroup[0,1,1] = 1
    spacegroup[0,2,2] = 1

  # write first line
  file.write('{0} {1} {2}\n'.format(dimensionality, centering, type))
  # write cell
  cell = structure.cell * float(structure.scale.rescale(angstrom))
  file.write( '{0[0]: > 18.12f} {0[1]: > 18.12f} {0[2]: > 18.12f}\n'           \
              '{1[0]: > 18.12f} {1[1]: > 18.12f} {1[2]: > 18.12f}\n'           \
              '{2[0]: > 18.12f} {2[1]: > 18.12f} {2[2]: > 18.12f}\n'           \
              .format( *(cell) ) )
  # write symmetry operators
  file.write('{0}\n'.format(len(spacegroup)))
  for op in spacegroup:
    file.write( '{0[0]: > 18.12f} {0[1]: > 18.12f} {0[2]: > 18.12f}\n'         \
                '{1[0]: > 18.12f} {1[1]: > 18.12f} {1[2]: > 18.12f}\n'         \
                '{2[0]: > 18.12f} {2[1]: > 18.12f} {2[2]: > 18.12f}\n'         \
                .format(*op[:3]) )
    file.write( '    {0[0]: > 18.12f} {0[1]: > 18.12f} {0[2]: > 18.12f}\n'     \
                .format(op[3]) )


  # figure out inequivalent atoms.
  groups = [u for u in equivalence_iterator(structure, spacegroup)]
  file.write('{0}\n'.format(len(groups)))
  for group in groups:
    atom = structure[group[0]]
    # Try figuring out atomic number.
    type = atom.type
    try: n = int(type)
    except: 
      try: n = find_specie(name=type)
      except:
        raise ValueError( 'Could not transform {0} to atomic number.'          \
                          .format(type) )
      else: type = n.atomic_number
    else: type = n
    pos = atom.pos * float(structure.scale.rescale(angstrom))
    file.write( '{0: >5} {1[0]: > 18.12f} {1[1]: > 18.12f} {1[2]: > 18.12f}\n' \
                .format(type, pos) )
