def test(path=None, dodate=True):
  from numpy import array, all, abs
  from lada.vasp import Extract
  from quantities import eV

  a = Extract(directory=path)
  assert a.success == True
  assert a.algo   == "Fast"
  assert a.is_dft == True
  assert a.is_gw  == False
  assert abs(a.encut - 245.3*eV) < 1e-8
  if dodate:
    assert repr(a.datetime) == 'datetime.datetime(2012, 3, 8, 21, 18, 29)'
  assert a.LDAUType is None
  assert len(a.HubbardU_NLEP) == 0
  assert a.pseudopotential == 'PAW_PBE'
  assert set(a.stoichiometry) == set([2])
  assert set(a.species) == set(['Si'])
  assert a.isif == 7
  assert abs(a.sigma - 0.2*eV) < 1e-8
  assert a.nsw == 50
  assert a.ibrion == 2
  assert a.relaxation == "volume"
  assert a.ispin == 1
  assert a.isym == 2
  assert a.name == 'has a name'
  assert a.system == "has a name"
  assert all(abs(array(a.ionic_charges) - [4.0]) < 1e-8)
  assert abs(a.nelect - 8.0) < 1e-8
  assert abs(a.extraelectron - 0e0) < 1e-8
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
  assert abs(a.ediffg-2e-4) < 1e-8
  assert a.smearing == 'metal'
  assert a.lorbit == 0
  assert abs(a.nupdown + 1.0) < 1e-8
  assert a.lmaxmix == 4
  assert abs(a.valence - 8.0) < 1e-8
  assert a.nonscf == False
  assert a.lwave == False
  assert a.lcharg
  assert a.lvtot == False
  assert a.nelm == 60
  assert a.nelmin == 2
  assert a.nelmdl == -5

if __name__ == "__main__":
  from sys import argv, path 
  from os.path import join
  if len(argv) > 2: path.extend(argv[2:])
  
  if len(argv) > 1: test(join(argv[1], 'COMMON'))

