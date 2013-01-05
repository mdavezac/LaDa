def test_doubleatomsymm():
  """ Bugs when calling atomsymm twice. """
  from tempfile import mkdtemp
  from shutil import rmtree
  from os.path import exists
  from pylada.dftcrystal import Crystal
  from pylada.misc import Changedir


  c = Crystal(136, 4.63909875, 2.97938395, 
              ifhr=0, 
              shift=0)                                  \
             .add_atom(0, 0, 0, 'Ti')                   \
             .add_atom(0.306153, 0.306153, 0, 'O')      \
             .append('ATOMSYMM')
  
  directory = '/tmp/test' #mkdtemp()
  if directory == '/tmp/test':
    if exists(directory): rmtree(directory)
    with Changedir(directory) as cwd: pass
  try: 
      c.append('ATOMSYMM')
      c.eval()
  finally:
    if directory != '/tmp/test': rmtree(directory)

if __name__ == '__main__':
  test_doubleatomsymm()
