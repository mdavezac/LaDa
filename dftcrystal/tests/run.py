def test():
  from tempfile import mkdtemp
  from shutil import rmtree
  from pylada.dftcrystal import Crystal, Functional, Shell
  from pylada import default_comm

  functional = Functional()
  functional.basis['Si'] = [
      Shell('s', a0=[16120.0, 0.001959],
                 a1=[2426.0, 0.01493], 
                 a2=[553.9, 0.07285],
                 a3=[156.3, 0.2461], 
                 a4=[50.07, 0.4859],
                 a5=[17.02, 0.325]),
      Shell('sp', a0=[292.7, -0.002781, 0.004438],
                  a1=[69.87, -0.03571, 0.03267],
                  a2=[22.34, -0.115, 0.1347],
                  a3=[8.15, 0.09356, 0.3287],
                  a4=[3.135, 0.603, 0.4496]), 
      Shell('sp', 4.0, a0=[1.22, 1.0, 1.0]),
      Shell('sp', 0.0, a0=[0.55, 1.0, 1.0]),
      Shell('sp', 0.0, a0=[0.27, 1.0, 1.0]) ]

  functional.dft.pbe0 = True
  functional.fmixing = 30
  functional.shrink = 8, 16
  functional.levshift = 5, True
  functional.maxcycle = 600

  crystal = Crystal(227, 5.43).add_atom(0.125, 0.125, 0.125, 'Si')
  directory = '/tmp/test/' # mkdtemp()
  try: 
     results = functional(crystal, outdir=directory, comm=default_comm)
     assert results.success
  finally: 
    if directory != '/tmp/test/': rmtree(directory)
    pass
if __name__ == '__main__':
  test()
