def test(path):
  from shutil import rmtree
  from tempfile import mkdtemp
  from pylada.crystal import Structure
  from pylada.vasp import Vasp
  from epirelax import epitaxial
  from pylada import default_comm

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
  try: 
    result = epitaxial(vasp, structure, outdir=directory, epiconv=1e-4, comm=default_comm)
    assert result.success
  finally: 
    rmtree(directory)
    pass

if __name__ == "__main__":
  from sys import argv
  test(argv[1])
