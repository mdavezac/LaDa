def test(path):
  from glob import glob
  from os.path import join
  from shutil import rmtree
  from tempfile import mkdtemp
  from numpy import all, abs
  from quantities import kbar, eV, angstrom
  from lada.crystal import Structure
  from lada.vasp import Vasp
  from lada.vasp.relax import relax, Relax

    
    
  structure = Structure([[0, 0.5, 0.5],[0.5, 0, 0.5], [0.5, 0.5, 0]], scale=5.43, name='has a name')\
                       .add_atom(0,0,0, "Si")\
                       .add_atom(0.25,0.25,0.25, "Si")

  vasp = Vasp()
  vasp.kpoints    = "Automatic generation\n0\nMonkhorst\n2 2 2\n0 0 0"
  vasp.precision  = "accurate"
  vasp.ediff      = 1e-5
  vasp.encut      = 1
  vasp.ismear     = "metal"
  vasp.sigma      = 0.01
  vasp.relaxation = "volume"
  vasp.add_specie = "Si", "{0}/pseudos/Si".format(path)
  directory = mkdtemp()
  try: 
    functional = Relax(copy=vasp)
    assert abs(functional.ediff - 1e-5) < 1e-8
    assert functional.precision == 'Accurate'
    result = functional(structure, outdir=directory, comm={'n': 2, 'ppn': 1},
                        relaxation="volume ionic cellshape")
    assert result.success
    def sortme(a): return int(a.split('/')[-1])
    dirs = sorted(glob(join(join(directory, '*'), '[0-9]')), key=sortme)
  # for previous, current in zip(dirs, dirs[1:]):
  #   assert len(check_output(['diff', join(previous, 'CONTCAR'), join(current, 'POSCAR')])) == 0
  # assert len(check_output(['diff', join(current, 'CONTCAR'), join(directory, 'POSCAR')])) == 0
    assert result.stress.units == kbar and all(abs(result.stress) < 1e0)
    assert result.forces.units == eV/angstrom and all(abs(result.forces) < 1e-1)
    assert result.total_energy.units == eV and all(abs(result.total_energy + 10.668652*eV) < 1e-2)

  finally: 
    if directory != '/tmp/test/relax': rmtree(directory)
    pass

if __name__ == "__main__":
  from sys import argv, path 
  if len(argv) > 2: path.extend(argv[2:])
  if len(argv) > 1: test(argv[1])
