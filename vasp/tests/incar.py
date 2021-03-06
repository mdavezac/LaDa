def test(path):
  from shutil import rmtree
  from tempfile import mkdtemp
  from os.path import join
  from quantities import eV
  from pylada.vasp import Vasp, read_incar
  from pylada.crystal import Structure
  structure = Structure([[0, 0.5, 0.5],[0.5, 0, 0.5], [0.5, 0.5, 0]], scale=5.43, name='has a name')\
                       .add_atom(0,0,0, "Si")\
                       .add_atom(0.25,0.25,0.25, "Si")

  vasp = Vasp()
  vasp.kpoints    = "Automatic generation\n0\nMonkhorst\n2 2 2\n0 0 0"
  vasp.precision  = "accurate"
  vasp.ediff      = 1e-5
  vasp.encut      = 1
  vasp.ismear     = "metal"
  vasp.sigma      = 0.06
  vasp.relaxation = "volume"
  vasp.add_specie = "Si", "{0}/pseudos/Si".format(path)

  directory = mkdtemp()
  try: 
    vasp.write_incar(path=join(directory, 'INCAR'), structure=structure)
    other = read_incar(join(directory, 'INCAR'))
    assert abs(other.ediff - 1e-5)  < 1e-8
    assert abs(other.encut - 245.345) < 1e-8
    assert abs(other.sigma - 0.06 * eV) < 1e-8
    assert other.ibrion     == 2
    assert other.icharg     == 'atomic'
    assert other.isif       == 7
    assert other.ismear     == 'metal'
    assert other.istart     == 'scratch'
    assert other.lcharg     == False
    assert other.nsw        == 50
    assert other.relaxation == 'volume'
    assert other.system     == 'has a name'
    with open(join(directory, 'INCAR'), 'a') as file:
      file.write('\nSOMETHing = 0.5\n')
    other = read_incar(join(directory, 'INCAR'))
    assert abs(other.ediff - 1e-5)  < 1e-8
    assert abs(other.encut - 245.345) < 1e-8
    assert abs(other.sigma - 0.06 * eV) < 1e-8
    assert other.ibrion     == 2
    assert other.icharg     == 'atomic'
    assert other.isif       == 7
    assert other.ismear     == 'metal'
    assert other.istart     == 'scratch'
    assert other.lcharg     == False
    assert other.nsw        == 50
    assert other.relaxation == 'volume'
    assert other.system     == 'has a name'
    assert 'something' in other._input
    assert isinstance(other.something, float)
    assert abs(other.something - 0.5) < 1e-8
  finally: 
    rmtree(directory)
    pass

if __name__ == '__main__':
  from sys import argv
  test(argv[1])
