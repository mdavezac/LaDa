""" Package to compute fft meshes from structural data. """
class SmallCells(object):
  """ Computes FFT mesh for small structures. 
  
      The FFT mesh is adjusted 
  """
  def __init__(self): 
    """ Initializes FFT mesh values. """
    super(SmallCells, self).__init__()

  def __call__(self, functional, structure, comm):
    """ Computes parameters of FFT mesh. """
    from operator import itemgetter
    from numpy import pi, sqrt
    from numpy.linalg import norm
    from lada.physics import a0
    from quantities import angstrom

    para = structure.scale*2.0*sqrt(functional.cutoff)/pi/a0.rescale(angstrom).magnitude
    result = [int(norm(structure.cell[:,0]) * para + 5e-1), \
              int(norm(structure.cell[:,1]) * para + 5e-1), \
              int(norm(structure.cell[:,2]) * para + 5e-1)]

    assert result[0] * result[1] >= comm.size,\
           ValueError("Too many comms for a system this size.")

    if result[0] * result[1] % comm.size != 0:
      a = 0 if result[0] % comm.size  == 0 else comm.size - result[0] % comm.size
      b = 0 if result[1] % comm.size  == 0 else comm.size - result[1] % comm.size
      all_prods = [(i, j, i*j) for i in xrange(result[0], result[0]+a+1)\
                               for j in xrange(result[1], result[1]+b+1)\
                               if (i*j) % comm.size == 0]
      result[0], result[1], dummy = min(all_prods, key = itemgetter(2))

    return result, result, (0,0,0)

  def __repr__(self):
    """ Representation string of this object. """
    return "{0.__class__.__name__}()".format(self)


class Nanowire(object):
  """ Computes FFT mesh for Lijun's nanowires. """
  def ___init__(self):
    """ Initializes FFT mesh values. """
    super(SmallCells, self).__init__()

  def __call__(self, functional, structure, comm):
    """ Computes parameters of FFT mesh. """
    from operator import itemgetter
    from numpy import pi, sqrt
    from numpy.linalg import norm
    from lada.physics import a0
    from quantities import angstrom

    para = structure.scale*2.0*sqrt(functional.cutoff)/pi/a0.rescale(angstrom).magnitude
    result = [int(norm(structure.cell[:,0]) * para + 5e-1), \
              int(norm(structure.cell[:,1]) * para + 5e-1), \
              int(norm(structure.cell[:,2]) * para + 5e-1)]

    multiple = comm.size if result[0] <= 500 else 3 * comm.size 
    if result[0] * result[1] % multiple != 0:
      a = 0 if result[0] % multiple  == 0 else multiple - result[0] % multiple
      b = 0 if result[1] % multiple  == 0 else multiple - result[1] % multiple
      all_prods = [(i, j, i*j) for i in xrange(result[0], result[0]+a+1)\
                               for j in xrange(result[1], result[1]+b+1)\
                               if (i*j) % multiple == 0]
      result[0], result[1], dummy = min(all_prods, key = itemgetter(2))

    return result, (result[0]//3, result[1]//3, result[2]), (result[0]//6,result[1]//6,0)

  def __repr__(self):
    """ Representation string of this object. """
    return "{0.__class__.__name__}()".format(self)
