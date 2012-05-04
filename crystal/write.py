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
          returns a string.
      :type file: str, stream, or None.
      :param bool vasp5:
          If true, include species in poscar, vasp-5 style.  Otherwise, looks
          for :py:data:`is_vasp_4 <lada.is_vasp_4>` global config variable. 
          Defaults to False, in which case, does not print specie types.
      :param substitute:
          If present, will substitute the atom type in the structure. Can be
          incomplete. Only works with vasp5 = True (or :py:data:`is_vasp_4
          <lada.is_vasp_4>` = True).
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
  if file is None:
    with open('POSCAR', 'w') as fileobj: return poscar(structure, fileobj, vasp5, substitute)
  elif not hasattr(file, 'write'):
    with open(file, 'w') as fileobj: return poscar(structure, fileobj, vasp5, substitute)

  from numpy import matrix, dot
  from . import specieset

  if vasp5 is None:
    import lada 
    vasp5 = not getattr(lada, 'is_vasp_4', True)

  string = "{0}\n{1}\n".format(getattr(structure, 'name', ''), structure.scale)
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
  string = "%BLOCK LATTICE_CART\n" \
           "  {0[0]} {0[1]} {0[2]}\n" \
           "  {1[0]} {1[1]} {1[2]}\n" \
           "  {2[0]} {2[1]} {2[2]}\n" \
           "%ENDBLOCK LATTICE_CART\n\n"\
           "%BLOCK POSITIONS_ABS\n".format(*(structure.cell.T*structure.scale))
  for atom in structure:
    string += "  {0} {1[0]} {1[1]} {1[2]}\n"\
              .format(atom.type, atom.pos * structure.scale)
  string += "%ENDBLOCK POSITION_ABS\n"
  if file == None: return string
  elif isinstance(file, str): 
    from ..misc import RelativePath
    with open(RelativePath(file).path, 'w') as file: file.write(string)
  else: file.write(string)
  
