def test(path):
  from shutil import rmtree
  from os.path import exists
  from os import makedirs
  from tempfile import mkdtemp
  from lada.crystal import Structure
  from lada.vasp import Vasp
  from lada.vasp.emass import mass

    
    
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
  if exists(directory) and directory == '/tmp/test': rmtree(directory)
  if not exists(directory): makedirs(directory)
  try: 
    result = mass(vasp, structure, outdir=directory, comm={'n': 2, 'ppn': 1})
    result.inverse_mass_tensor
    assert result.success
  finally: 
    if directory != '/tmp/test': rmtree(directory)
    pass

if __name__ == "__main__":
  from sys import argv, path 
  if len(argv) > 2: path.extend(argv[2:])
  if len(argv) > 1: test(argv[1])
