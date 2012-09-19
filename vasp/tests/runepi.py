def test(path):
  from shutil import rmtree
  from os.path import exists
  from os import makedirs
  from tempfile import mkdtemp
  from numpy import abs
  from lada.crystal import Structure
  from lada.vasp import Vasp
  from lada.vasp.relax import epitaxial

    
    
  structure = Structure([[0, 0.5, 0.5],[0.5, 0, 0.5], [0.5, 0.5, 0]], scale=5.55, name='has a name')\
                       .add_atom(0,0,0, "Si")\
                       .add_atom(0.25,0.25,0.25, "Si")

  vasp = Vasp()
  vasp.kpoints    = "Automatic generation\n0\nMonkhorst\n2 2 2\n0 0 0"
  vasp.prec       = "accurate"
  vasp.ediff      = 1e-5
  vasp.encut      = 1.4
  vasp.ismear     = "fermi"
  vasp.sigma      = 0.01
  vasp.relaxation = "volume"
  vasp.add_specie = "Si", "{0}/pseudos/Si".format(path)
  directory = mkdtemp()
  if exists(directory) and directory == '/tmp/test': rmtree(directory)
  if not exists(directory): makedirs(directory)
  try: 
    result = epitaxial(vasp, structure, outdir=directory, epiconv=1e-5, comm={'n': 2, 'ppn': 1})
    assert result.success
    assert abs(result.stress[2,2]) < 1.0
  finally: 
    if directory != '/tmp/test': rmtree(directory)
    pass

if __name__ == "__main__":
  from sys import argv
  test(argv[1])
