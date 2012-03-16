def test(path):
  from subprocess import check_output, STDOUT, Popen, PIPE
  from glob import glob
  from os.path import join
  from shutil import rmtree
  from tempfile import mkdtemp
  from numpy import all, abs
  from quantities import kbar, eV, angstrom
  from lada.crystal import Structure
  from lada.vasp import Vasp
  from epirelax import epitaxial

  structure = Structure([[0, 0.5, 0.5],[0.5, 0, 0.5], [0.5, 0.5, 0]], scale=5.55, name='has a name')\
                       .add_atom(0,0,0, "Si")\
                       .add_atom(0.25,0.25,0.25, "Si")

  vasp = Vasp()
  vasp.kpoints    = "Automatic generation\n0\nMonkhorst\n2 2 2\n0 0 0"
  vasp.precision  = "accurate"
  vasp.ediff      = 1e-5
  vasp.encut      = 1.4
  vasp.set_smearing   = "metal", 0.01
  vasp.relaxation = "volume"
  vasp.add_specie = "Si", "{0}/pseudos/Si".format(path)
  directory = mkdtemp()
  try: 
    result = epitaxial(vasp, structure, outdir=directory, epiconv=1e-4, comm={'n': 2, 'ppn': 1})
    assert result.success
  finally: 
    rmtree(directory)
    pass

if __name__ == "__main__":
  from sys import argv, path 
  if len(argv) > 1:
    path.extend(argv[1:])
    test(argv[1])
