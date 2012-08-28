def test():
  from tempfile import mkdtemp
  from shutil import rmtree
  from lada.dftcrystal import Crystal, relax, Shell

  functional = relax.Relax()
  functional.basis['Si'] = [
      Shell( 's', 2.0, a0=[910.655, 0.0660823],
             a1=[137.336, 0.3862290], a2=[29.7601, 0.67238] ),
      Shell( 'sp', 8.0, a0=[36.6716, -0.1045110, 0.113355],
             a2=[8.31729, 0.10741, 0.4575780],
             a3=[2.21645, 0.951446, 0.607427] ),
      Shell( 'sp', 4.0, a0=[1.07913, -0.376108, 0.067103],
             a1=[0.302422, 1.25165, 0.956883] ),
      Shell( 'sp', 0.0, a0=[0.0933392, 1.0, 1.0] ) ]
  functional.dft.pbe0 = True
  functional.fmixing = 30
  functional.optgeom.fulloptg = True
  functional.shrink = 4, None
  functional.levshift = 5, True

  crystal = Crystal(227, 5.43).add_atom(0.125, 0.125, 0.125, 'Si')
  directory = '/tmp/test' #mkdtemp()
  try: 
     results = functional(crystal, outdir=directory, comm=None)
     assert results.success
  finally: 
#   if directory != '/tmp/test/': rmtree(directory)
    pass
if __name__ == '__main__':
  test()
