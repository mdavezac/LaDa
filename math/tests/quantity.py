""" Checks atom methods and attributes. """
def test():
  """ Test atom initialization. """
  from _quantity import is_quantity, fromC, fromPy, get_angstrom
  from numpy import abs
  from quantities import meter, angstrom, eV# , Ry

  assert not is_quantity("a");
  assert not is_quantity(0);
  assert not is_quantity(0.4);
  assert is_quantity(eV)
  assert is_quantity(5*eV)

  assert abs(fromC(5, 'meter') - 5*meter) < 1e-8
  assert abs(fromC(float((5*meter).rescale(angstrom)), 'angstrom') - 5*meter) < 1e-8
  assert abs(fromC(float((4*meter).rescale(angstrom)), 'angstrom') - 5*meter) > 1e-8
  assert abs(fromPy(5, meter) - 5*meter) < 1e-8
  assert hasattr(fromPy(5, meter), 'rescale')
  assert abs(get_angstrom(5*meter) - float(5*meter.rescale(angstrom))) < 1e-8
  assert abs(get_angstrom(5) - 5) < 1e-8
  assert abs(get_angstrom(5.5) - 5.5) < 1e-8
  assert get_angstrom(4*eV) is None

if __name__ == "__main__":
  from sys import argv, path 
  if len(argv) > 0: path.extend(argv[1:])
  test()
