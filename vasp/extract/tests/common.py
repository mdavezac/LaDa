def test(directory=None):
  from os.path import join
  from numpy import array, all, abs
  from lada.vasp import Extract
  from quantities import eV

  a = Extract(directory=join(directory, 'COMMON'))
  assert a.algo   == "fast"
  assert a.is_dft == True
  assert a.is_gw  == False
  assert abs(a.encut - 245.3*eV) < 1e-8
  assert a.success == True
  assert repr(a.datetime) == 'datetime.datetime(2012, 3, 1, 9, 11, 26)'
  assert a.LDAUType is None
  assert len(a.HubbardU_NLEP) == 0
  assert a.pseudopotential == 'PAW_PBE'
  assert set(a.stoichiometry) == set([2])
  assert set(a.species) == set(['Si'])
  assert a.isif == 7
  assert a.nsw == 50
  assert a.ibrion == 2
  assert a.relaxation == "volume"
  assert a.ispin == 1
  assert a.name == 'has a name'
  assert a.system == "has a name"
  assert all(abs(array(a.ionic_charges) - [4.0]) < 1e-8)
  assert abs(a.nelect - 8.0) < 1e-8
  assert abs(a.nbands - 8) < 1e-8
  assert all(abs(1.00778*array([[0, 0.5, 0.5],[0.5, 0, 0.5],[0.5, 0.5, 0]])-a.structure.cell) < 1e-4) 
  assert all(abs(a.structure.scale-5.43) < 1e-4) 
  assert all(abs(a.structure.scale*a.structure.cell - a._grep_structure.cell) < 1e-4) 
  assert all(abs(a.structure[0].pos) < 1e-8)
  assert all(abs(a.structure[1].pos - 0.251945) < 1e-6)
  assert all([b.type == 'Si' for b in a.structure])
  assert abs(a.structure.energy + 10.665642 * eV) < 1e-6
  assert abs(a.sigma-0.2*eV) < 1e-6
  assert abs(a.ismear-1) < 1e-6
  assert abs(a.potim-0.5) < 1e-6
  assert abs(a.istart-0) < 1e-6
  assert abs(a.icharg-2) < 1e-6
  assert a.precision == "accurate"
  assert abs(a.ediff-2e-5) < 1e-8
  assert abs(a.ediffg-2e-5) < 1e-8

if __name__ == "__main__":
  from sys import argv, path 
  if len(argv) > 2: path.extend(argv[2:])
  
  if len(argv) > 1: test(argv[1])

