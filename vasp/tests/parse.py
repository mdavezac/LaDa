def test():
  from lada.vasp import parse_incar
  from lada.error import ValueError

  string = """ALGO = Fast\n"""\
           """ENCUT = 294.414\n"""\
           """EDIFF = 2e-12\n"""\
           """MAGMOM = -0*8 1.5*2\n"""\
           """ISPIN = 1\n"""\
           """ISMEAR = -1\n"""\
           """SIGMA = 0.001\n"""\
           """ISIF = 2\n"""\
           """NSW = 50\n"""\
           """IBRION = 2\n"""\
           """LMAXMIX = 4\n"""\
           """LCHARG = .TRUE.\n"""\
           """LVTOT = .FALSE.\n"""\
           """SYSTEM = Zinc-Blende\n"""


  def get_errors(found):
    errors = {}
    expected = { 'ALGO': 'Fast', 'ENCUT': '294.414', 'EDIFF': '2e-12',
                 'ISPIN': '1', 'MAGMOM': '0*8 1.5*2', 'ISMEAR': '-1',
                 'SIGMA': '0.001', 'ISIF': '2', 'NSW': '50', 'IBRION': '2',
                 'LMAXMIX': '4', 'LCHARG': '.TRUE.', 'LVTOT': '.FALSE',
                 'SYSTEM': 'Zinc-Blende' }
    for key in set(found.keys()) - set(expected.keys()): 
      errors[key] = found[key]
    for key in set(expected.keys()) - set(found.keys()): 
      errors[key] = None
    return errors
 
  result = parse_incar(string)
  assert len(get_errors(result)) == 0
  assert parse_incar(string.replace('\n', '\n#')) == [('ALGO', 'Fast')]
  assert len(get_errors(parse_incar(string.replace('\n', ';', 2)))) == 0
  assert len(get_errors(parse_incar(string.replace('=', '\\\n  =', 2)))) == 0
  assert len(get_errors(parse_incar(string.replace('=', '=  \\\n  ', 2)))) == 0

  try: parse_incar( string + "LVTOT = .TRUE.") 
  except ValueError: pass
  else: raise Exception()
  try: parse_incar( string + "   = .TRUE.") 
  except ValueError: pass
  else: raise Exception()

if __name__ == '__main__': 
  test()
   
    

