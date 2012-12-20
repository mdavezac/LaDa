def test():
  """ Tests file setup. """
  from tempfile import mkdtemp
  from os.path import exists
  from shutil import rmtree
  from quantities import eV
  from lada.gulp import Functional
  from lada.crystal import Structure
  from lada.tools import create_directory

  functional = Functional()
  functional.qeq = True
  functional.morse.enabled     = True
  functional.morse['Ti', 'O']  = [1.0279493, 3.640737, 1.88265, 0, 25]
  functional.morse['Ti', 'Ti'] = [0.00567139, 1.5543, 4.18784, 0, 25]
  functional.morse['O', 'O']   = [0.00567139, 1.5543, 4.18784, 0, 25] 

  structure = Structure( 4.6391, 0, 0,
                         0, 4.6391, 0,
                         0, 0, 2.97938,
                         scale=1, symmgroup=136 )                              \
                .add_atom(0, 0, 0, 'Ti', asymmetric=True)                      \
                .add_atom(2.31955, 2.31955, 1.48969, 'Ti', asymmetric=False)   \
                .add_atom(1.42027, 1.42027, 0, 'O', asymmetric=True)           \
                .add_atom(-1.42027, -1.42027, 0, 'O', asymmetric=False)        \
                .add_atom(-0.899275, 0.899275, 1.48969, 'O', asymmetric=False) \
                .add_atom(0.899275, -0.899275, 1.48969, 'O', asymmetric=False)

  directory = mkdtemp()
  if directory == '/tmp/test/':
    if exists(directory): rmtree(directory)
    create_directory(directory)
  try:
    result = functional(structure, directory)
    assert result._finished
    assert result.success
    assert not result.optimize
    assert not result.conp
    assert abs(result.energy + 9.84502039*eV) < 1e-5
  finally:
    if directory != '/tmp/test/': rmtree(directory)

if __name__ == '__main__':
  test()
