""" Standard A2BX4 lattices. """


def b1():
  """ Returns b1 lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (1.440000e-01, 2.880000e-01, 0.000000e+00),\
                     (1.440000e-01, 0.000000e+00, 2.880000e-01),\
                     (5.000000e-01, 0.000000e+00, 0.000000e+00)
  lattice.name = "b1"
  lattice.add_site = (2.880000e-01, 2.880000e-01, 3.500000e-01), "A", 
  lattice.add_site = (1.440000e-01, 1.440000e-01, 1.500000e-01), "A", 
  lattice.add_site = (0.000000e+00, 0.000000e+00, 0.000000e+00), "B", 
  lattice.add_site = (0.000000e+00, 1.440000e-01, 0.000000e+00), "X", 
  lattice.add_site = (1.440000e-01, 0.000000e+00, 0.000000e+00), "X", 
  lattice.add_site = (2.880000e-01, 2.880000e-01, 1.500000e-01), "X", 
  lattice.add_site = (1.440000e-01, 1.440000e-01, 3.500000e-01), "X", 
  lattice.find_space_group()
  return lattice


def b10():
  """ Returns b10 lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (1.000000e+00, 0.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 5.850985e-01, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 5.224911e-01)
  lattice.name = "b10"
  lattice.add_site = (0.000000e+00, 0.000000e+00, 0.000000e+00), "A", 
  lattice.add_site = (5.000000e-01, 0.000000e+00, 2.612456e-01), "A", 
  lattice.add_site = (0.000000e+00, 2.925493e-01, 0.000000e+00), "A", 
  lattice.add_site = (5.000000e-01, 2.925493e-01, 2.612456e-01), "A", 
  lattice.add_site = (2.177349e-01, 1.462746e-01, 2.666468e-01), "A", 
  lattice.add_site = (2.822651e-01, 4.388239e-01, 5.401262e-03), "A", 
  lattice.add_site = (7.822651e-01, 4.388239e-01, 2.558443e-01), "A", 
  lattice.add_site = (7.177349e-01, 1.462746e-01, 5.170898e-01), "A", 
  lattice.add_site = (3.956987e-01, 1.462746e-01, 1.828752e-02), "B", 
  lattice.add_site = (1.043013e-01, 4.388239e-01, 2.795331e-01), "B", 
  lattice.add_site = (6.043013e-01, 4.388239e-01, 5.042036e-01), "B", 
  lattice.add_site = (8.956987e-01, 1.462746e-01, 2.429580e-01), "B", 
  lattice.add_site = (4.064161e-01, 1.462746e-01, 3.566835e-01), "X", 
  lattice.add_site = (9.358386e-02, 4.388239e-01, 9.543795e-02), "X", 
  lattice.add_site = (5.935839e-01, 4.388239e-01, 1.658076e-01), "X", 
  lattice.add_site = (9.064161e-01, 1.462746e-01, 4.270532e-01), "X", 
  lattice.add_site = (5.598078e-01, 1.462746e-01, 1.173047e-01), "X", 
  lattice.add_site = (9.401922e-01, 4.388239e-01, 3.785502e-01), "X", 
  lattice.add_site = (4.401922e-01, 4.388239e-01, 4.051864e-01), "X", 
  lattice.add_site = (5.980776e-02, 1.462746e-01, 1.439409e-01), "X", 
  lattice.add_site = (3.320130e-01, 1.467676e-02, 1.401257e-01), "X", 
  lattice.add_site = (1.679870e-01, -1.467676e-02, 4.013712e-01), "X", 
  lattice.add_site = (6.679870e-01, 3.072260e-01, 3.823654e-01), "X", 
  lattice.add_site = (8.320130e-01, 2.778725e-01, 1.211199e-01), "X", 
  lattice.add_site = (6.679870e-01, -1.467676e-02, 3.823654e-01), "X", 
  lattice.add_site = (8.320130e-01, 1.467676e-02, 1.211199e-01), "X", 
  lattice.add_site = (3.320130e-01, 2.778725e-01, 1.401257e-01), "X", 
  lattice.add_site = (1.679870e-01, 3.072260e-01, 4.013712e-01), "X", 
  lattice.find_space_group()
  return lattice


def b11():
  """ Returns b11 lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (7.430000e-01, 0.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 5.730000e-01, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 1.000000e+00)
  lattice.name = "b11"
  lattice.add_site = (5.028624e-01, 1.432500e-01, 4.182000e-01), "A", 
  lattice.add_site = (6.116376e-01, 4.297500e-01, 9.182000e-01), "A", 
  lattice.add_site = (2.401376e-01, 4.297500e-01, 5.818000e-01), "A", 
  lattice.add_site = (1.313624e-01, 1.432500e-01, 8.180000e-02), "A", 
  lattice.add_site = (7.344555e-01, 1.432500e-01, 7.046000e-01), "A", 
  lattice.add_site = (3.800445e-01, 4.297500e-01, 2.046000e-01), "A", 
  lattice.add_site = (8.544500e-03, 4.297500e-01, 2.954000e-01), "A", 
  lattice.add_site = (3.629555e-01, 1.432500e-01, 7.954000e-01), "A", 
  lattice.add_site = (1.751994e-01, 1.432500e-01, 4.155000e-01), "B", 
  lattice.add_site = (1.963006e-01, 4.297500e-01, 9.155000e-01), "B", 
  lattice.add_site = (5.678006e-01, 4.297500e-01, 5.845000e-01), "B", 
  lattice.add_site = (5.466994e-01, 1.432500e-01, 8.450000e-02), "B", 
  lattice.add_site = (2.340450e-02, 1.432500e-01, 4.032000e-01), "X", 
  lattice.add_site = (3.480955e-01, 4.297500e-01, 9.032000e-01), "X", 
  lattice.add_site = (7.195955e-01, 4.297500e-01, 5.968000e-01), "X", 
  lattice.add_site = (3.949045e-01, 1.432500e-01, 9.680000e-02), "X", 
  lattice.add_site = (2.206710e-01, 1.432500e-01, 5.579000e-01), "X", 
  lattice.add_site = (1.508290e-01, 4.297500e-01, 5.790000e-02), "X", 
  lattice.add_site = (5.223290e-01, 4.297500e-01, 4.421000e-01), "X", 
  lattice.add_site = (5.921710e-01, 1.432500e-01, 9.421000e-01), "X", 
  lattice.add_site = (2.226771e-01, 2.349300e-02, 3.484000e-01), "X", 
  lattice.add_site = (1.488229e-01, 5.495070e-01, 8.484000e-01), "X", 
  lattice.add_site = (5.203229e-01, 3.099930e-01, 6.516000e-01), "X", 
  lattice.add_site = (5.941771e-01, 2.630070e-01, 1.516000e-01), "X", 
  lattice.add_site = (5.203229e-01, 5.495070e-01, 6.516000e-01), "X", 
  lattice.add_site = (5.941771e-01, 2.349300e-02, 1.516000e-01), "X", 
  lattice.add_site = (2.226771e-01, 2.630070e-01, 3.484000e-01), "X", 
  lattice.add_site = (1.488229e-01, 3.099930e-01, 8.484000e-01), "X", 
  lattice.find_space_group()
  return lattice


def b12():
  """ Returns b12 lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (5.840000e-01, 0.000000e+00, -5.500000e-02),\
                     (0.000000e+00, 7.340000e-01, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 1.000000e+00)
  lattice.name = "b12"
  lattice.add_site = (4.425800e-01, 6.172940e-01, 7.600000e-02), "A", 
  lattice.add_site = (-1.780800e-01, 2.502940e-01, 4.240000e-01), "A", 
  lattice.add_site = (8.642000e-02, 1.167060e-01, 9.240000e-01), "A", 
  lattice.add_site = (1.230800e-01, 4.837060e-01, 5.760000e-01), "A", 
  lattice.add_site = (1.226260e-01, 3.655320e-01, 2.020000e-01), "A", 
  lattice.add_site = (-4.421260e-01, 7.325320e-01, 2.980000e-01), "A", 
  lattice.add_site = (4.063740e-01, 3.684680e-01, 7.980000e-01), "A", 
  lattice.add_site = (3.871260e-01, 1.468000e-03, 7.020000e-01), "A", 
  lattice.add_site = (4.359910e-01, 2.025840e-01, 7.900000e-02), "B", 
  lattice.add_site = (-1.714910e-01, 5.695840e-01, 4.210000e-01), "B", 
  lattice.add_site = (9.300900e-02, 5.314160e-01, 9.210000e-01), "B", 
  lattice.add_site = (1.164910e-01, 1.644160e-01, 5.790000e-01), "B", 
  lattice.add_site = (3.810620e-01, 5.872000e-03, 4.300000e-01), "X", 
  lattice.add_site = (-1.165620e-01, 3.728720e-01, 7.000000e-02), "X", 
  lattice.add_site = (1.479380e-01, 7.281280e-01, 5.700000e-01), "X", 
  lattice.add_site = (6.156200e-02, 3.611280e-01, 9.300000e-01), "X", 
  lattice.add_site = (3.212280e-01, 1.688200e-01, 1.800000e-01), "X", 
  lattice.add_site = (-5.672800e-02, 5.358200e-01, 3.200000e-01), "X", 
  lattice.add_site = (2.077720e-01, 5.651800e-01, 8.200000e-01), "X", 
  lattice.add_site = (1.728000e-03, 1.981800e-01, 6.800000e-01), "X", 
  lattice.add_site = (5.714970e-01, 1.394600e-01, 1.530000e-01), "X", 
  lattice.add_site = (-3.069970e-01, 5.064600e-01, 3.470000e-01), "X", 
  lattice.add_site = (-4.249700e-02, 5.945400e-01, 8.470000e-01), "X", 
  lattice.add_site = (2.519970e-01, 2.275400e-01, 6.530000e-01), "X", 
  lattice.add_site = (3.262500e-01, 1.064300e-01, 9.700000e-01), "X", 
  lattice.add_site = (-1.167500e-01, 4.734300e-01, 5.300000e-01), "X", 
  lattice.add_site = (2.027500e-01, 6.275700e-01, 3.000000e-02), "X", 
  lattice.add_site = (6.175000e-02, 2.605700e-01, 4.700000e-01), "X", 
  lattice.find_space_group()
  return lattice


def b18():
  """ Returns b18 lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (0.000000e+00, 4.760000e-01, 4.760000e-01),\
                     (1.000000e+00, 0.000000e+00, 1.000000e+00),\
                     (7.980000e-01, 7.980000e-01, 0.000000e+00)
  lattice.name = "b18"
  lattice.add_site = (5.950000e-01, 1.250000e+00, 6.981316e-01), "A", 
  lattice.add_site = (3.570000e-01, 7.500000e-01, 8.978684e-01), "A", 
  lattice.add_site = (5.950000e-01, 1.250000e+00, 1.296868e+00), "A", 
  lattice.add_site = (3.570000e-01, 7.500000e-01, 2.991316e-01), "A", 
  lattice.add_site = (1.190000e-01, 2.500000e-01, 1.995000e-01), "B", 
  lattice.add_site = (8.330000e-01, 1.750000e+00, 1.396500e+00), "B", 
  lattice.add_site = (4.098286e-01, 1.066227e+00, 3.924881e-01), "X", 
  lattice.add_site = (5.421714e-01, 9.337735e-01, 1.203512e+00), "X", 
  lattice.add_site = (4.098286e-01, 4.337735e-01, 8.045119e-01), "X", 
  lattice.add_site = (5.421714e-01, 1.566227e+00, 7.914881e-01), "X", 
  lattice.add_site = (3.041714e-01, 1.066227e+00, 8.045119e-01), "X", 
  lattice.add_site = (6.478286e-01, 9.337735e-01, 7.914881e-01), "X", 
  lattice.add_site = (3.041714e-01, 4.337735e-01, 3.924881e-01), "X", 
  lattice.add_site = (6.478286e-01, 1.566227e+00, 1.203512e+00), "X", 
  lattice.find_space_group()
  return lattice


def b19():
  """ Returns b19 lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (3.165000e-01, 0.000000e+00, 6.330000e-01),\
                     (5.000000e-01, 0.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 7.710000e-01, 0.000000e+00)
  lattice.name = "b19"
  lattice.add_site = (6.330000e-01, 2.000000e-01, 1.927500e-01), "A", 
  lattice.add_site = (3.165000e-01, 3.000000e-01, 5.782500e-01), "A", 
  lattice.add_site = (3.165000e-01, 1.110223e-16, 0.000000e+00), "A", 
  lattice.add_site = (3.165000e-01, 1.110223e-16, 3.855000e-01), "A", 
  lattice.add_site = (3.165000e-01, 3.500000e-01, 1.927500e-01), "B", 
  lattice.add_site = (6.330000e-01, 1.500000e-01, 5.782500e-01), "B", 
  lattice.add_site = (8.102400e-01, 4.500000e-01, 1.927500e-01), "X", 
  lattice.add_site = (1.392600e-01, 5.000000e-02, 5.782500e-01), "X", 
  lattice.add_site = (4.557600e-01, 4.500000e-01, 1.927500e-01), "X", 
  lattice.add_site = (4.937400e-01, 5.000000e-02, 5.782500e-01), "X", 
  lattice.add_site = (6.330000e-01, 2.500000e-01, 4.394700e-01), "X", 
  lattice.add_site = (3.165000e-01, 2.500000e-01, 3.315300e-01), "X", 
  lattice.add_site = (3.165000e-01, 2.500000e-01, 5.397000e-02), "X", 
  lattice.add_site = (6.330000e-01, 2.500000e-01, 7.170300e-01), "X", 
  lattice.find_space_group()
  return lattice


def b2():
  """ Returns b2 lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (1.000000e+00, 0.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 1.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 5.900000e-01)
  lattice.name = "b2"
  lattice.add_site = (0.000000e+00, 5.000000e-01, 2.950000e-01), "A", 
  lattice.add_site = (5.000000e-01, 0.000000e+00, 2.950000e-01), "A", 
  lattice.add_site = (0.000000e+00, 0.000000e+00, 0.000000e+00), "B", 
  lattice.add_site = (2.324500e-01, 2.324500e-01, 0.000000e+00), "X", 
  lattice.add_site = (7.675500e-01, 7.675500e-01, 0.000000e+00), "X", 
  lattice.add_site = (7.675500e-01, 2.324500e-01, 0.000000e+00), "X", 
  lattice.add_site = (2.324500e-01, 7.675500e-01, 0.000000e+00), "X", 
  lattice.find_space_group()
  return lattice


def b20():
  """ Returns b20 lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (2.980000e-01, 2.980000e-01, 0.000000e+00),\
                     (-5.160000e-01, 5.160000e-01, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 1.000000e+00)
  lattice.name = "b20"
  lattice.add_site = (0.000000e+00, 0.000000e+00, 0.000000e+00), "A", 
  lattice.add_site = (0.000000e+00, 0.000000e+00, 5.000000e-01), "A", 
  lattice.add_site = (2.980000e-01, 1.720034e-01, 8.000000e-01), "A", 
  lattice.add_site = (2.980000e-01, -1.720034e-01, 2.000000e-01), "A", 
  lattice.add_site = (2.980000e-01, 1.720034e-01, 1.700000e-01), "B", 
  lattice.add_site = (2.980000e-01, -1.720034e-01, 6.700000e-01), "B", 
  lattice.add_site = (2.980000e-01, 1.720034e-01, 0.000000e+00), "X", 
  lattice.add_site = (2.980000e-01, -1.720034e-01, 5.000000e-01), "X", 
  lattice.add_site = (3.784600e-01, 3.147600e-01, 2.500000e-01), "X", 
  lattice.add_site = (1.341000e-01, 1.702800e-01, 2.500000e-01), "X", 
  lattice.add_site = (3.814400e-01, 3.096000e-02, 2.500000e-01), "X", 
  lattice.add_site = (2.175400e-01, -3.147600e-01, 7.500000e-01), "X", 
  lattice.add_site = (4.619000e-01, -1.702800e-01, 7.500000e-01), "X", 
  lattice.add_site = (2.145600e-01, -3.096000e-02, 7.500000e-01), "X", 
  lattice.find_space_group()
  return lattice


def b21():
  """ Returns b21 lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (3.840000e-01, 3.840000e-01, 0.000000e+00),\
                     (-6.640000e-01, 6.640000e-01, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 1.000000e+00)
  lattice.name = "b21"
  lattice.add_site = (3.840000e-01, 2.217760e-01, 5.000000e-02), "A", 
  lattice.add_site = (-3.840000e-01, -2.217760e-01, -5.000000e-02), "A", 
  lattice.add_site = (3.840000e-01, -2.217760e-01, 5.500000e-01), "A", 
  lattice.add_site = (-3.840000e-01, 2.217760e-01, -5.500000e-01), "A", 
  lattice.add_site = (0.000000e+00, 0.000000e+00, 2.500000e-01), "B", 
  lattice.add_site = (0.000000e+00, 0.000000e+00, 7.500000e-01), "B", 
  lattice.add_site = (3.840000e-01, 2.217760e-01, 2.500000e-01), "X", 
  lattice.add_site = (3.855360e-01, -2.191200e-01, 7.500000e-01), "X", 
  lattice.add_site = (2.561280e-01, -4.428880e-01, 0.000000e+00), "X", 
  lattice.add_site = (2.561280e-01, 4.428880e-01, 0.000000e+00), "X", 
  lattice.add_site = (-5.122560e-01, 0.000000e+00, 0.000000e+00), "X", 
  lattice.add_site = (-2.561280e-01, 4.428880e-01, 5.000000e-01), "X", 
  lattice.add_site = (-2.561280e-01, -4.428880e-01, 5.000000e-01), "X", 
  lattice.add_site = (5.122560e-01, 0.000000e+00, 5.000000e-01), "X", 
  lattice.find_space_group()
  return lattice


def b33():
  """ Returns b33 lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (5.770000e-01, 5.770000e-01, 5.770000e-01),\
                     (-1.000000e+00, 1.000000e+00, -3.330000e-01),\
                     (0.000000e+00, 0.000000e+00, 2.570000e-01)
  lattice.name = "b33"
  lattice.add_site = (6.791290e-01, 5.384090e-01, 5.833900e-02), "A", 
  lattice.add_site = (3.485080e-01, 3.194090e-01, 5.833900e-02), "A", 
  lattice.add_site = (7.033630e-01, 1.424090e-01, 5.833900e-02), "A", 
  lattice.add_site = (4.748710e-01, 1.285910e-01, 1.986610e-01), "A", 
  lattice.add_site = (8.054920e-01, 3.475910e-01, 1.986610e-01), "A", 
  lattice.add_site = (4.506370e-01, 5.245910e-01, 1.986610e-01), "A", 
  lattice.add_site = (6.693200e-01, -7.859110e-01, 6.861900e-02), "A", 
  lattice.add_site = (9.226230e-01, -2.691100e-02, 6.861900e-02), "A", 
  lattice.add_site = (1.390570e-01, -1.869110e-01, 6.861900e-02), "A", 
  lattice.add_site = (4.846800e-01, -5.470890e-01, 1.883810e-01), "A", 
  lattice.add_site = (8.083770e-01, -3.060890e-01, 1.883810e-01), "A", 
  lattice.add_site = (4.379430e-01, -1.460890e-01, 1.883810e-01), "A", 
  lattice.add_site = (4.760135e-01, 7.845307e-01, 1.932820e-01), "B", 
  lattice.add_site = (2.371989e-01, 2.060069e-02, 1.932820e-01), "B", 
  lattice.add_site = (1.017788e+00, 1.956207e-01, 1.932820e-01), "B", 
  lattice.add_site = (6.779865e-01, -1.175307e-01, 6.371801e-02), "B", 
  lattice.add_site = (3.398011e-01, -3.536007e-01, 6.371801e-02), "B", 
  lattice.add_site = (7.132124e-01, -5.286207e-01, 6.371801e-02), "B", 
  lattice.add_site = (9.676867e-01, -9.925000e-02, 1.927500e-01), "X", 
  lattice.add_site = (1.792739e-01, -1.114500e-01, 1.927500e-01), "X", 
  lattice.add_site = (5.840394e-01, -7.885500e-01, 1.927500e-01), "X", 
  lattice.add_site = (7.633133e-01, -2.337500e-01, 6.425000e-02), "X", 
  lattice.add_site = (3.977261e-01, -2.215500e-01, 6.425000e-02), "X", 
  lattice.add_site = (5.699606e-01, -5.444500e-01, 6.425000e-02), "X", 
  lattice.add_site = (3.886095e-01, 6.717502e-01, 1.929042e-01), "X", 
  lattice.add_site = (3.785120e-01, 1.250200e-03, 1.929042e-01), "X", 
  lattice.add_site = (9.638785e-01, 3.277502e-01, 1.929042e-01), "X", 
  lattice.add_site = (7.653905e-01, -4.750200e-03, 6.409580e-02), "X", 
  lattice.add_site = (1.984880e-01, -3.342502e-01, 6.409580e-02), "X", 
  lattice.add_site = (7.671215e-01, -6.607502e-01, 6.409580e-02), "X", 
  lattice.add_site = (5.580167e-01, 7.853993e-01, 7.656030e-02), "X", 
  lattice.add_site = (1.953145e-01, 9.099930e-02, 7.656030e-02), "X", 
  lattice.add_site = (9.776688e-01, 1.238993e-01, 7.656030e-02), "X", 
  lattice.add_site = (5.959833e-01, -1.183993e-01, 1.804397e-01), "X", 
  lattice.add_site = (3.816855e-01, -4.239993e-01, 1.804397e-01), "X", 
  lattice.add_site = (7.533312e-01, -4.568993e-01, 1.804397e-01), "X", 
  lattice.add_site = (5.922905e-01, 5.540972e-01, 2.034412e-01), "X", 
  lattice.add_site = (3.785120e-01, 2.365972e-01, 2.034412e-01), "X", 
  lattice.add_site = (7.601975e-01, 2.100972e-01, 2.034412e-01), "X", 
  lattice.add_site = (5.617095e-01, 1.129028e-01, 5.355880e-02), "X", 
  lattice.add_site = (7.754880e-01, 4.304028e-01, 5.355880e-02), "X", 
  lattice.add_site = (3.938025e-01, 4.569028e-01, 5.355880e-02), "X", 
  lattice.find_space_group()
  return lattice


def b34():
  """ Returns b34 lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (8.080000e-01, 0.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 2.810000e-01, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 1.000000e+00)
  lattice.name = "b34"
  lattice.add_site = (5.454000e-02, 7.025000e-02, 3.901000e-01), "A", 
  lattice.add_site = (3.494600e-01, 2.107500e-01, 8.901000e-01), "A", 
  lattice.add_site = (7.534600e-01, 2.107500e-01, 6.099000e-01), "A", 
  lattice.add_site = (4.585400e-01, 7.025000e-02, 1.099000e-01), "A", 
  lattice.add_site = (8.427440e-02, 7.025000e-02, 9.056000e-01), "A", 
  lattice.add_site = (3.197256e-01, 2.107500e-01, 4.056000e-01), "A", 
  lattice.add_site = (7.237256e-01, 2.107500e-01, 9.440000e-02), "A", 
  lattice.add_site = (4.882744e-01, 7.025000e-02, 5.944000e-01), "A", 
  lattice.add_site = (2.003032e-01, 7.025000e-02, 6.513000e-01), "B", 
  lattice.add_site = (2.036968e-01, 2.107500e-01, 1.513000e-01), "B", 
  lattice.add_site = (6.076968e-01, 2.107500e-01, 3.487000e-01), "B", 
  lattice.add_site = (6.043032e-01, 7.025000e-02, 8.487000e-01), "B", 
  lattice.add_site = (2.672056e-01, 7.025000e-02, 1.410000e-02), "X", 
  lattice.add_site = (1.367944e-01, 2.107500e-01, 5.141000e-01), "X", 
  lattice.add_site = (5.407944e-01, 2.107500e-01, 9.859000e-01), "X", 
  lattice.add_site = (6.712056e-01, 7.025000e-02, 4.859000e-01), "X", 
  lattice.add_site = (2.108072e-01, 7.025000e-02, 2.997000e-01), "X", 
  lattice.add_site = (1.931928e-01, 2.107500e-01, 7.997000e-01), "X", 
  lattice.add_site = (5.971928e-01, 2.107500e-01, 7.003000e-01), "X", 
  lattice.add_site = (6.148072e-01, 7.025000e-02, 2.003000e-01), "X", 
  lattice.add_site = (4.444000e-02, 7.025000e-02, 9.270000e-02), "X", 
  lattice.add_site = (3.595600e-01, 2.107500e-01, 5.927000e-01), "X", 
  lattice.add_site = (7.635600e-01, 2.107500e-01, 9.073000e-01), "X", 
  lattice.add_site = (4.484400e-01, 7.025000e-02, 4.073000e-01), "X", 
  lattice.add_site = (1.462480e-02, 7.025000e-02, 7.120000e-01), "X", 
  lattice.add_site = (3.893752e-01, 2.107500e-01, 2.120000e-01), "X", 
  lattice.add_site = (7.933752e-01, 2.107500e-01, 2.880000e-01), "X", 
  lattice.add_site = (4.186248e-01, 7.025000e-02, 7.880000e-01), "X", 
  lattice.find_space_group()
  return lattice


def b35():
  """ Returns b35 lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (8.100000e-01, 0.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 2.660000e-01, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 1.000000e+00)
  lattice.name = "b35"
  lattice.add_site = (1.174500e-01, 6.650000e-02, 6.200000e-02), "A", 
  lattice.add_site = (2.875500e-01, 1.995000e-01, 5.620000e-01), "A", 
  lattice.add_site = (6.925500e-01, 1.995000e-01, 9.380000e-01), "A", 
  lattice.add_site = (5.224500e-01, 6.650000e-02, 4.380000e-01), "A", 
  lattice.add_site = (7.800300e-01, 6.650000e-02, 6.140000e-01), "A", 
  lattice.add_site = (4.349700e-01, 1.995000e-01, 1.140000e-01), "A", 
  lattice.add_site = (2.997000e-02, 1.995000e-01, 3.860000e-01), "A", 
  lattice.add_site = (3.750300e-01, 6.650000e-02, 8.860000e-01), "A", 
  lattice.add_site = (2.559600e-01, 6.650000e-02, 3.340000e-01), "B", 
  lattice.add_site = (1.490400e-01, 1.995000e-01, 8.340000e-01), "B", 
  lattice.add_site = (5.540400e-01, 1.995000e-01, 6.660000e-01), "B", 
  lattice.add_site = (6.609600e-01, 6.650000e-02, 1.660000e-01), "B", 
  lattice.add_site = (6.520500e-01, 6.650000e-02, 7.280000e-01), "X", 
  lattice.add_site = (5.629500e-01, 1.995000e-01, 2.280000e-01), "X", 
  lattice.add_site = (1.579500e-01, 1.995000e-01, 2.720000e-01), "X", 
  lattice.add_site = (2.470500e-01, 6.650000e-02, 7.720000e-01), "X", 
  lattice.add_site = (3.434400e-01, 6.650000e-02, 1.840000e-01), "X", 
  lattice.add_site = (6.156000e-02, 1.995000e-01, 6.840000e-01), "X", 
  lattice.add_site = (4.665600e-01, 1.995000e-01, 8.160000e-01), "X", 
  lattice.add_site = (7.484400e-01, 6.650000e-02, 3.160000e-01), "X", 
  lattice.add_site = (1.830600e-01, 6.650000e-02, 4.920000e-01), "X", 
  lattice.add_site = (2.219400e-01, 1.995000e-01, 9.920000e-01), "X", 
  lattice.add_site = (6.269400e-01, 1.995000e-01, 5.080000e-01), "X", 
  lattice.add_site = (5.880600e-01, 6.650000e-02, 8.000000e-03), "X", 
  lattice.add_site = (4.446900e-01, 6.650000e-02, 5.950000e-01), "X", 
  lattice.add_site = (7.703100e-01, 1.995000e-01, 9.500000e-02), "X", 
  lattice.add_site = (3.653100e-01, 1.995000e-01, 4.050000e-01), "X", 
  lattice.add_site = (3.969000e-02, 6.650000e-02, 9.050000e-01), "X", 
  lattice.find_space_group()
  return lattice


def b36():
  """ Returns b36 lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (1.000000e+00, 0.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 4.020000e-01, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 2.140000e-01)
  lattice.name = "b36"
  lattice.add_site = (-1.277000e-01, -2.026080e-01, 5.350000e-02), "A", 
  lattice.add_site = (3.723000e-01, 4.036080e-01, 1.605000e-01), "A", 
  lattice.add_site = (1.277000e-01, 2.026080e-01, -5.350000e-02), "A", 
  lattice.add_site = (-3.723000e-01, -4.036080e-01, -1.605000e-01), "A", 
  lattice.add_site = (4.600000e-03, 3.233286e-01, 5.350000e-02), "A", 
  lattice.add_site = (5.046000e-01, -1.223286e-01, 1.605000e-01), "A", 
  lattice.add_site = (-4.600000e-03, -3.233286e-01, -5.350000e-02), "A", 
  lattice.add_site = (-5.046000e-01, 1.223286e-01, -1.605000e-01), "A", 
  lattice.add_site = (1.937000e-01, 3.505440e-02, 5.350000e-02), "B", 
  lattice.add_site = (6.937000e-01, 1.659456e-01, 1.605000e-01), "B", 
  lattice.add_site = (-1.937000e-01, -3.505440e-02, -5.350000e-02), "B", 
  lattice.add_site = (-6.937000e-01, -1.659456e-01, -1.605000e-01), "B", 
  lattice.add_site = (-2.212000e-01, -1.441974e-01, 5.350000e-02), "X", 
  lattice.add_site = (2.788000e-01, 3.451974e-01, 1.605000e-01), "X", 
  lattice.add_site = (2.212000e-01, 1.441974e-01, -5.350000e-02), "X", 
  lattice.add_site = (-2.788000e-01, -3.451974e-01, -1.605000e-01), "X", 
  lattice.add_site = (1.571000e-01, 2.712294e-01, 5.350000e-02), "X", 
  lattice.add_site = (6.571000e-01, -7.022940e-02, 1.605000e-01), "X", 
  lattice.add_site = (-1.571000e-01, -2.712294e-01, -5.350000e-02), "X", 
  lattice.add_site = (-6.571000e-01, 7.022940e-02, -1.605000e-01), "X", 
  lattice.add_site = (-7.790000e-02, -3.931158e-01, 5.350000e-02), "X", 
  lattice.add_site = (4.221000e-01, 5.941158e-01, 1.605000e-01), "X", 
  lattice.add_site = (7.790000e-02, 3.931158e-01, -5.350000e-02), "X", 
  lattice.add_site = (-4.221000e-01, -5.941158e-01, -1.605000e-01), "X", 
  lattice.add_site = (5.190000e-02, 1.361976e-01, 5.350000e-02), "X", 
  lattice.add_site = (5.519000e-01, 6.480240e-02, 1.605000e-01), "X", 
  lattice.add_site = (-5.190000e-02, -1.361976e-01, -5.350000e-02), "X", 
  lattice.add_site = (-5.519000e-01, -6.480240e-02, -1.605000e-01), "X", 
  lattice.find_space_group()
  return lattice


def b37():
  """ Returns b37 lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (5.168000e-02, 5.168000e-02, 0.000000e+00),\
                     (-8.951000e-02, 8.951000e-02, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 1.000000e+00)
  lattice.name = "b37"
  lattice.add_site = (0.000000e+00, 0.000000e+00, 4.320000e-01), "A", 
  lattice.add_site = (0.000000e+00, 0.000000e+00, 5.680000e-01), "A", 
  lattice.add_site = (5.168000e-02, -2.983673e-02, 7.653330e-01), "A", 
  lattice.add_site = (5.168000e-02, -2.983673e-02, 9.013330e-01), "A", 
  lattice.add_site = (5.168000e-02, 2.983673e-02, 9.866700e-02), "A", 
  lattice.add_site = (5.168000e-02, 2.983673e-02, 2.346670e-01), "A", 
  lattice.add_site = (0.000000e+00, 0.000000e+00, 0.000000e+00), "B", 
  lattice.add_site = (5.168000e-02, -2.983673e-02, 3.333330e-01), "B", 
  lattice.add_site = (5.168000e-02, 2.983673e-02, 6.666670e-01), "B", 
  lattice.add_site = (0.000000e+00, 0.000000e+00, 1.460000e-01), "X", 
  lattice.add_site = (0.000000e+00, 0.000000e+00, 8.540000e-01), "X", 
  lattice.add_site = (5.168000e-02, -2.983673e-02, 4.793330e-01), "X", 
  lattice.add_site = (5.168000e-02, -2.983673e-02, 1.873330e-01), "X", 
  lattice.add_site = (5.168000e-02, 2.983673e-02, 8.126670e-01), "X", 
  lattice.add_site = (5.168000e-02, 2.983673e-02, 5.206670e-01), "X", 
  lattice.add_site = (0.000000e+00, 0.000000e+00, 2.840000e-01), "X", 
  lattice.add_site = (0.000000e+00, 0.000000e+00, 7.160000e-01), "X", 
  lattice.add_site = (5.168000e-02, -2.983673e-02, 6.173330e-01), "X", 
  lattice.add_site = (5.168000e-02, -2.983673e-02, 4.933300e-02), "X", 
  lattice.add_site = (5.168000e-02, 2.983673e-02, 9.506670e-01), "X", 
  lattice.add_site = (5.168000e-02, 2.983673e-02, 3.826670e-01), "X", 
  lattice.find_space_group()
  return lattice


def b38():
  """ Returns b38 lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (5.194000e-02, 5.194000e-02, 0.000000e+00),\
                     (-8.997000e-02, 8.997000e-02, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 1.000000e+00)
  lattice.name = "b38"
  lattice.add_site = (0.000000e+00, 0.000000e+00, 1.660000e-01), "A", 
  lattice.add_site = (5.194000e-02, -2.999006e-02, 4.993330e-01), "A", 
  lattice.add_site = (5.194000e-02, 2.999006e-02, 8.326670e-01), "A", 
  lattice.add_site = (0.000000e+00, 0.000000e+00, 9.370000e-01), "A", 
  lattice.add_site = (5.194000e-02, -2.999006e-02, 2.703330e-01), "A", 
  lattice.add_site = (5.194000e-02, 2.999006e-02, 6.036670e-01), "A", 
  lattice.add_site = (0.000000e+00, 0.000000e+00, 3.960000e-01), "B", 
  lattice.add_site = (5.194000e-02, -2.999006e-02, 7.293330e-01), "B", 
  lattice.add_site = (5.194000e-02, 2.999006e-02, 6.266700e-02), "B", 
  lattice.add_site = (0.000000e+00, 0.000000e+00, 4.000000e-02), "X", 
  lattice.add_site = (5.194000e-02, -2.999006e-02, 3.733330e-01), "X", 
  lattice.add_site = (5.194000e-02, 2.999006e-02, 7.066670e-01), "X", 
  lattice.add_site = (0.000000e+00, 0.000000e+00, 2.940000e-01), "X", 
  lattice.add_site = (5.194000e-02, -2.999006e-02, 6.273330e-01), "X", 
  lattice.add_site = (5.194000e-02, 2.999006e-02, 9.606670e-01), "X", 
  lattice.add_site = (0.000000e+00, 0.000000e+00, 4.590000e-01), "X", 
  lattice.add_site = (5.194000e-02, -2.999006e-02, 7.923330e-01), "X", 
  lattice.add_site = (5.194000e-02, 2.999006e-02, 1.256670e-01), "X", 
  lattice.add_site = (0.000000e+00, 0.000000e+00, 8.720000e-01), "X", 
  lattice.add_site = (5.194000e-02, -2.999006e-02, 2.053330e-01), "X", 
  lattice.add_site = (5.194000e-02, 2.999006e-02, 5.386670e-01), "X", 
  lattice.find_space_group()
  return lattice


def b4():
  """ Returns b4 lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (1.000000e+00, 0.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 1.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 1.000000e+00)
  lattice.name = "b4"
  lattice.add_site = (0.000000e+00, 5.000000e-01, 5.000000e-01), "A", 
  lattice.add_site = (5.000000e-01, 0.000000e+00, 5.000000e-01), "A", 
  lattice.add_site = (0.000000e+00, 0.000000e+00, 0.000000e+00), "B", 
  lattice.add_site = (2.700000e-01, 2.700000e-01, 2.250000e-01), "X", 
  lattice.add_site = (-2.700000e-01, -2.700000e-01, 2.250000e-01), "X", 
  lattice.add_site = (-2.700000e-01, 2.700000e-01, -2.250000e-01), "X", 
  lattice.add_site = (2.700000e-01, -2.700000e-01, -2.250000e-01), "X", 
  lattice.find_space_group()
  return lattice

def b5D():
  """ Returns b5D lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (0.000000e+00, 5.000000e-01, 5.000000e-01),\
                     (5.000000e-01, 0.000000e+00, 5.000000e-01),\
                     (5.000000e-01, 5.000000e-01, 0.000000e+00)
  lattice.name = "b5D"
  lattice.add_site = (5.000000e-01, 5.000000e-01, 5.000000e-01), "A", 
  lattice.add_site = (5.000000e-01, 2.500000e-01, 2.500000e-01), "A", 
  lattice.add_site = (2.500000e-01, 5.000000e-01, 2.500000e-01), "A", 
  lattice.add_site = (8.750000e-01, 8.750000e-01, 8.750000e-01), "A", 
  lattice.add_site = (2.500000e-01, 2.500000e-01, 5.000000e-01), "B", 
  lattice.add_site = (1.250000e-01, 1.250000e-01, 1.250000e-01), "B", 
  lattice.add_site = (2.500000e-01, 2.500000e-01, 2.500000e-01), "X", 
  lattice.add_site = (2.500000e-01, 5.000000e-01, 5.000000e-01), "X", 
  lattice.add_site = (5.000000e-01, 2.500000e-01, 5.000000e-01), "X", 
  lattice.add_site = (5.000000e-01, 5.000000e-01, 2.500000e-01), "X", 
  lattice.add_site = (7.500000e-01, 7.500000e-01, 7.500000e-01), "X", 
  lattice.add_site = (7.500000e-01, 5.000000e-01, 5.000000e-01), "X", 
  lattice.add_site = (5.000000e-01, 7.500000e-01, 5.000000e-01), "X", 
  lattice.add_site = (5.000000e-01, 5.000000e-01, 7.500000e-01), "X", 
  lattice.find_space_group()
  return lattice

def b5(u=0.25):
  """ Returns b5 lattice."""
  from . import Lattice

  x = u
  y = 0.25 - x

  lattice = Lattice()
  lattice.scale = 1.0
  lattice.name = "b5"

  lattice.set_cell = (0.000000e+00, 5.000000e-01, 5.000000e-01),\
                     (5.000000e-01, 0.000000e+00, 5.000000e-01),\
                     (5.000000e-01, 5.000000e-01, 0.000000e+00)

  lattice.add_site = (5.000000e-01, 5.000000e-01, 5.000000e-01), "A",
  lattice.add_site = (5.000000e-01, 2.500000e-01, 2.500000e-01), "A",
  lattice.add_site = (2.500000e-01, 5.000000e-01, 2.500000e-01), "A",
  lattice.add_site = (2.500000e-01, 2.500000e-01, 5.000000e-01), "A",
  lattice.add_site = (8.750000e-01, 8.750000e-01, 8.750000e-01), "B",
  lattice.add_site = (1.250000e-01, 1.250000e-01, 1.250000e-01), "B",
  lattice.add_site = (     x,     x,     x), "X",
  lattice.add_site = (     x,     y,     y), "X",
  lattice.add_site = (     y,     x,     y), "X",
  lattice.add_site = (     y,     y,     x), "X",
  lattice.add_site = (    -x,    -x,    -x), "X",
  lattice.add_site = (    -x,    -y,    -y), "X",
  lattice.add_site = (    -y,    -x,    -y), "X",
  lattice.add_site = (    -y,    -y,    -x), "X",
  lattice.find_space_group()

  return lattice

def b6():
  """ Returns b6 lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (3.050000e-01, 6.100000e-01, 0.000000e+00),\
                     (3.050000e-01, 0.000000e+00, 6.100000e-01),\
                     (5.000000e-01, 0.000000e+00, 0.000000e+00)
  lattice.name = "b6"
  lattice.add_site = (3.050000e-01, 3.050000e-01, 0.000000e+00), "A", 
  lattice.add_site = (3.050000e-01, 0.000000e+00, 0.000000e+00), "A", 
  lattice.add_site = (4.575000e-01, 1.525000e-01, 2.500000e-01), "A", 
  lattice.add_site = (1.525000e-01, 1.525000e-01, 2.500000e-01), "A", 
  lattice.add_site = (6.100000e-01, 4.575000e-01, 1.250000e-01), "B", 
  lattice.add_site = (3.050000e-01, 4.575000e-01, 3.750000e-01), "B", 
  lattice.add_site = (6.100000e-01, 2.880420e-01, 2.589000e-01), "X", 
  lattice.add_site = (6.100000e-01, 6.269580e-01, 2.589000e-01), "X", 
  lattice.add_site = (1.694580e-01, 1.525000e-01, 8.900000e-03), "X", 
  lattice.add_site = (4.405420e-01, 1.525000e-01, 8.900000e-03), "X", 
  lattice.add_site = (3.050000e-01, 2.880420e-01, 2.411000e-01), "X", 
  lattice.add_site = (3.050000e-01, 6.269580e-01, 2.411000e-01), "X", 
  lattice.add_site = (7.455420e-01, 7.625000e-01, 4.911000e-01), "X", 
  lattice.add_site = (4.744580e-01, 7.625000e-01, 4.911000e-01), "X", 
  lattice.find_space_group()
  return lattice


def b7():
  """ Returns b7 lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (3.880000e-01, 7.760000e-01, 0.000000e+00),\
                     (3.880000e-01, 0.000000e+00, 7.760000e-01),\
                     (5.000000e-01, 0.000000e+00, 0.000000e+00)
  lattice.name = "b7"
  lattice.add_site = (4.074000e-01, 1.940000e-01, 1.250000e-01), "A", 
  lattice.add_site = (3.686000e-01, 5.820000e-01, 1.250000e-01), "A", 
  lattice.add_site = (5.820000e-01, 7.566000e-01, 3.750000e-01), "A", 
  lattice.add_site = (9.700000e-01, 7.954000e-01, 3.750000e-01), "A", 
  lattice.add_site = (0.000000e+00, 0.000000e+00, 0.000000e+00), "B", 
  lattice.add_site = (7.760000e-01, 3.880000e-01, 2.500000e-01), "B", 
  lattice.add_site = (2.172800e-01, 7.977280e-01, 1.200000e-01), "X", 
  lattice.add_site = (5.587200e-01, 7.542720e-01, 1.200000e-01), "X", 
  lattice.add_site = (4.097280e-01, 9.467200e-01, 3.800000e-01), "X", 
  lattice.add_site = (3.662720e-01, 6.052800e-01, 3.800000e-01), "X", 
  lattice.add_site = (5.587200e-01, 4.097280e-01, 1.300000e-01), "X", 
  lattice.add_site = (2.172800e-01, 3.662720e-01, 1.300000e-01), "X", 
  lattice.add_site = (7.542720e-01, 9.467200e-01, 3.700000e-01), "X", 
  lattice.add_site = (7.977280e-01, 6.052800e-01, 3.700000e-01), "X", 
  lattice.find_space_group()
  return lattice


def b8():
  """ Returns b8 lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (5.000000e-01, 1.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 1.026000e+00),\
                     (1.610000e-01, 0.000000e+00, 0.000000e+00)
  lattice.name = "b8"
  lattice.add_site = (1.260000e-01, 6.669000e-02, 0.000000e+00), "A", 
  lattice.add_site = (8.740000e-01, 9.593100e-01, 0.000000e+00), "A", 
  lattice.add_site = (1.260000e-01, 4.463100e-01, 0.000000e+00), "A", 
  lattice.add_site = (8.740000e-01, 5.796900e-01, 0.000000e+00), "A", 
  lattice.add_site = (3.840000e-01, 2.565000e-01, 0.000000e+00), "B", 
  lattice.add_site = (6.160000e-01, 7.695000e-01, 0.000000e+00), "B", 
  lattice.add_site = (5.000000e-01, 0.000000e+00, 0.000000e+00), "X", 
  lattice.add_site = (5.000000e-01, 5.130000e-01, 0.000000e+00), "X", 
  lattice.add_site = (6.000000e-02, 2.565000e-01, 0.000000e+00), "X", 
  lattice.add_site = (9.400000e-01, 7.695000e-01, 0.000000e+00), "X", 
  lattice.add_site = (2.360000e-01, 9.234000e-01, 0.000000e+00), "X", 
  lattice.add_site = (7.640000e-01, 1.026000e-01, 0.000000e+00), "X", 
  lattice.add_site = (2.360000e-01, 6.156000e-01, 0.000000e+00), "X", 
  lattice.add_site = (7.640000e-01, 4.104000e-01, 0.000000e+00), "X", 
  lattice.find_space_group()
  return lattice


def b9():
  """ Returns b9 lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (8.620000e-01, 0.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 1.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 2.820000e-01)
  lattice.name = "b9"
  lattice.add_site = (3.732460e-01, 6.100000e-01, 7.050000e-02), "A", 
  lattice.add_site = (-3.732460e-01, -6.100000e-01, -7.050000e-02), "A", 
  lattice.add_site = (8.042460e-01, -1.100000e-01, 7.050000e-02), "A", 
  lattice.add_site = (-8.042460e-01, 1.100000e-01, -7.050000e-02), "A", 
  lattice.add_site = (3.620400e-01, 1.080000e-01, 7.050000e-02), "A", 
  lattice.add_site = (-3.620400e-01, -1.080000e-01, -7.050000e-02), "A", 
  lattice.add_site = (7.930400e-01, 3.920000e-01, 7.050000e-02), "A", 
  lattice.add_site = (-7.930400e-01, -3.920000e-01, -7.050000e-02), "A", 
  lattice.add_site = (6.516720e-01, 6.540000e-01, 7.050000e-02), "B", 
  lattice.add_site = (-6.516720e-01, -6.540000e-01, -7.050000e-02), "B", 
  lattice.add_site = (1.082672e+00, -1.540000e-01, 7.050000e-02), "B", 
  lattice.add_site = (-1.082672e+00, 1.540000e-01, -7.050000e-02), "B", 
  lattice.add_site = (1.792960e-01, 1.620000e-01, 7.050000e-02), "X", 
  lattice.add_site = (-1.792960e-01, -1.620000e-01, -7.050000e-02), "X", 
  lattice.add_site = (6.102960e-01, 3.380000e-01, 7.050000e-02), "X", 
  lattice.add_site = (-6.102960e-01, -3.380000e-01, -7.050000e-02), "X", 
  lattice.add_site = (9.913000e-02, 4.770000e-01, 7.050000e-02), "X", 
  lattice.add_site = (-9.913000e-02, -4.770000e-01, -7.050000e-02), "X", 
  lattice.add_site = (5.301300e-01, 2.300000e-02, 7.050000e-02), "X", 
  lattice.add_site = (-5.301300e-01, -2.300000e-02, -7.050000e-02), "X", 
  lattice.add_site = (4.491020e-01, 7.840000e-01, 7.050000e-02), "X", 
  lattice.add_site = (-4.491020e-01, -7.840000e-01, -7.050000e-02), "X", 
  lattice.add_site = (8.801020e-01, -2.840000e-01, 7.050000e-02), "X", 
  lattice.add_site = (-8.801020e-01, 2.840000e-01, -7.050000e-02), "X", 
  lattice.add_site = (3.611780e-01, 4.240000e-01, 7.050000e-02), "X", 
  lattice.add_site = (-3.611780e-01, -4.240000e-01, -7.050000e-02), "X", 
  lattice.add_site = (7.921780e-01, 7.600000e-02, 7.050000e-02), "X", 
  lattice.add_site = (-7.921780e-01, -7.600000e-02, -7.050000e-02), "X", 
  lattice.find_space_group()
  return lattice


def d1():
  """ Returns d1 lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (1.000000e+00, 0.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 1.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 7.450000e-01)
  lattice.name = "d1"
  lattice.add_site = (1.400000e-01, 1.630000e-01, 0.000000e+00), "A", 
  lattice.add_site = (8.600000e-01, 8.370000e-01, 0.000000e+00), "A", 
  lattice.add_site = (8.370000e-01, 1.400000e-01, 3.725000e-01), "A", 
  lattice.add_site = (1.630000e-01, 8.600000e-01, 3.725000e-01), "A", 
  lattice.add_site = (3.600000e-01, 6.630000e-01, 0.000000e+00), "A", 
  lattice.add_site = (6.400000e-01, 3.370000e-01, 0.000000e+00), "A", 
  lattice.add_site = (6.630000e-01, 6.400000e-01, 3.725000e-01), "A", 
  lattice.add_site = (3.370000e-01, 3.600000e-01, 3.725000e-01), "A", 
  lattice.add_site = (0.000000e+00, 5.000000e-01, 1.862500e-01), "B", 
  lattice.add_site = (5.000000e-01, 0.000000e+00, 5.587500e-01), "B", 
  lattice.add_site = (0.000000e+00, 5.000000e-01, 5.587500e-01), "B", 
  lattice.add_site = (5.000000e-01, 0.000000e+00, 1.862500e-01), "B", 
  lattice.add_site = (6.710000e-01, 1.710000e-01, 1.862500e-01), "X", 
  lattice.add_site = (3.290000e-01, 8.290000e-01, 1.862500e-01), "X", 
  lattice.add_site = (8.290000e-01, 6.710000e-01, 5.587500e-01), "X", 
  lattice.add_site = (1.710000e-01, 3.290000e-01, 5.587500e-01), "X", 
  lattice.add_site = (3.290000e-01, 8.290000e-01, 5.587500e-01), "X", 
  lattice.add_site = (6.710000e-01, 1.710000e-01, 5.587500e-01), "X", 
  lattice.add_site = (1.710000e-01, 3.290000e-01, 1.862500e-01), "X", 
  lattice.add_site = (8.290000e-01, 6.710000e-01, 1.862500e-01), "X", 
  lattice.add_site = (9.600000e-02, 6.370000e-01, 0.000000e+00), "X", 
  lattice.add_site = (9.040000e-01, 3.630000e-01, 0.000000e+00), "X", 
  lattice.add_site = (3.630000e-01, 9.600000e-02, 3.725000e-01), "X", 
  lattice.add_site = (6.370000e-01, 9.040000e-01, 3.725000e-01), "X", 
  lattice.add_site = (4.040000e-01, 1.370000e-01, 0.000000e+00), "X", 
  lattice.add_site = (5.960000e-01, 8.630000e-01, 0.000000e+00), "X", 
  lattice.add_site = (1.370000e-01, 5.960000e-01, 3.725000e-01), "X", 
  lattice.add_site = (8.630000e-01, 4.040000e-01, 3.725000e-01), "X", 
  lattice.find_space_group()
  return lattice


def d3():
  """ Returns d3 lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (5.000000e-01, -7.400000e-01, 1.000000e+00),\
                     (1.495000e-01, 0.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 4.380000e-01, 0.000000e+00)
  lattice.name = "d3"
  lattice.add_site = (-1.850000e-01, 0.000000e+00, 1.095000e-01), "A", 
  lattice.add_site = (-5.550000e-01, 0.000000e+00, 3.285000e-01), "A", 
  lattice.add_site = (0.000000e+00, 0.000000e+00, 0.000000e+00), "B", 
  lattice.add_site = (-5.920000e-03, 0.000000e+00, 2.006040e-01), "X", 
  lattice.add_site = (2.659200e-01, 0.000000e+00, 2.373960e-01), "X", 
  lattice.add_site = (-3.759200e-01, 0.000000e+00, 4.196040e-01), "X", 
  lattice.add_site = (6.359200e-01, 0.000000e+00, 1.839600e-02), "X", 
  lattice.find_space_group()
  return lattice


def d9():
  """ Returns d9 lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (5.000000e-01, 1.000000e+00, 0.000000e+00),\
                     (5.000000e-01, 0.000000e+00, 1.000000e+00),\
                     (5.000000e-01, 0.000000e+00, 0.000000e+00)
  lattice.name = "d9"
  lattice.add_site = (2.500000e-01, 3.750000e-01, 0.000000e+00), "A", 
  lattice.add_site = (7.500000e-01, 1.250000e-01, 0.000000e+00), "A", 
  lattice.add_site = (1.000000e+00, 1.250000e+00, 3.750000e-01), "A", 
  lattice.add_site = (1.000000e+00, 7.500000e-01, 1.250000e-01), "A", 
  lattice.add_site = (3.750000e-01, 1.000000e+00, 2.500000e-01), "B", 
  lattice.add_site = (6.250000e-01, 5.000000e-01, 2.500000e-01), "B", 
  lattice.add_site = (8.300000e-02, 8.300000e-02, 8.300000e-02), "X", 
  lattice.add_site = (9.170000e-01, 4.170000e-01, 8.300000e-02), "X", 
  lattice.add_site = (9.170000e-01, 5.830000e-01, 4.170000e-01), "X", 
  lattice.add_site = (1.083000e+00, 9.170000e-01, 4.170000e-01), "X", 
  lattice.add_site = (3.330000e-01, 3.330000e-01, 3.330000e-01), "X", 
  lattice.add_site = (6.670000e-01, 1.167000e+00, 3.330000e-01), "X", 
  lattice.add_site = (3.330000e-01, 6.670000e-01, 1.670000e-01), "X", 
  lattice.add_site = (6.670000e-01, 8.330000e-01, 1.670000e-01), "X", 
  lattice.find_space_group()
  return lattice


def s1():
  """ Returns s1 lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (1.000000e+00, 0.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 6.010000e-01, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 9.940000e-01)
  lattice.name = "s1"
  lattice.add_site = (3.802000e-01, 5.356112e-01, 8.081220e-01), "A", 
  lattice.add_site = (6.198000e-01, 6.538880e-02, 3.111220e-01), "A", 
  lattice.add_site = (8.802000e-01, 6.538880e-02, 8.081220e-01), "A", 
  lattice.add_site = (1.198000e-01, 5.356112e-01, 3.111220e-01), "A", 
  lattice.add_site = (6.335000e-01, 5.191438e-01, 1.143100e-02), "A", 
  lattice.add_site = (3.665000e-01, 8.185620e-02, 5.084310e-01), "A", 
  lattice.add_site = (1.335000e-01, 8.185620e-02, 1.143100e-02), "A", 
  lattice.add_site = (8.665000e-01, 5.191438e-01, 5.084310e-01), "A", 
  lattice.add_site = (3.815000e-01, 2.279593e-01, 1.688806e-01), "A", 
  lattice.add_site = (6.185000e-01, 3.730407e-01, 6.658806e-01), "A", 
  lattice.add_site = (8.815000e-01, 3.730407e-01, 1.688806e-01), "A", 
  lattice.add_site = (1.185000e-01, 2.279593e-01, 6.658806e-01), "A", 
  lattice.add_site = (6.391000e-01, 2.193650e-01, 9.861474e-01), "A", 
  lattice.add_site = (3.609000e-01, 3.816350e-01, 4.891474e-01), "A", 
  lattice.add_site = (1.391000e-01, 3.816350e-01, 9.861474e-01), "A", 
  lattice.add_site = (8.609000e-01, 2.193650e-01, 4.891474e-01), "A", 
  lattice.add_site = (3.798000e-01, 5.336880e-01, 1.851822e-01), "B", 
  lattice.add_site = (6.202000e-01, 6.731200e-02, 6.821822e-01), "B", 
  lattice.add_site = (8.798000e-01, 6.731200e-02, 1.851822e-01), "B", 
  lattice.add_site = (1.202000e-01, 5.336880e-01, 6.821822e-01), "B", 
  lattice.add_site = (3.966000e-01, 2.325870e-01, 8.082214e-01), "B", 
  lattice.add_site = (6.034000e-01, 3.684130e-01, 3.112214e-01), "B", 
  lattice.add_site = (8.966000e-01, 3.684130e-01, 8.082214e-01), "B", 
  lattice.add_site = (1.034000e-01, 2.325870e-01, 3.112214e-01), "B", 
  lattice.add_site = (5.386000e-01, 3.593980e-01, 1.214668e-01), "X", 
  lattice.add_site = (4.614000e-01, 2.416020e-01, 6.184668e-01), "X", 
  lattice.add_site = (3.860000e-02, 2.416020e-01, 1.214668e-01), "X", 
  lattice.add_site = (9.614000e-01, 3.593980e-01, 6.184668e-01), "X", 
  lattice.add_site = (5.201000e-01, 3.618020e-01, 8.572256e-01), "X", 
  lattice.add_site = (4.799000e-01, 2.391980e-01, 3.602256e-01), "X", 
  lattice.add_site = (2.010000e-02, 2.391980e-01, 8.572256e-01), "X", 
  lattice.add_site = (9.799000e-01, 3.618020e-01, 3.602256e-01), "X", 
  lattice.add_site = (2.863000e-01, 5.288800e-01, 3.479000e-03), "X", 
  lattice.add_site = (7.137000e-01, 7.212000e-02, 5.004790e-01), "X", 
  lattice.add_site = (7.863000e-01, 7.212000e-02, 3.479000e-03), "X", 
  lattice.add_site = (2.137000e-01, 5.288800e-01, 5.004790e-01), "X", 
  lattice.add_site = (7.693000e-01, 5.336279e-01, 7.267134e-01), "X", 
  lattice.add_site = (2.307000e-01, 6.737210e-02, 2.297134e-01), "X", 
  lattice.add_site = (2.693000e-01, 6.737210e-02, 7.267134e-01), "X", 
  lattice.add_site = (7.307000e-01, 5.336279e-01, 2.297134e-01), "X", 
  lattice.add_site = (5.111000e-01, 6.851400e-02, 1.245482e-01), "X", 
  lattice.add_site = (4.889000e-01, 5.324860e-01, 6.215482e-01), "X", 
  lattice.add_site = (1.110000e-02, 5.324860e-01, 1.245482e-01), "X", 
  lattice.add_site = (9.889000e-01, 6.851400e-02, 6.215482e-01), "X", 
  lattice.add_site = (5.427000e-01, 8.233700e-02, 8.589154e-01), "X", 
  lattice.add_site = (4.573000e-01, 5.186630e-01, 3.619154e-01), "X", 
  lattice.add_site = (4.270000e-02, 5.186630e-01, 8.589154e-01), "X", 
  lattice.add_site = (9.573000e-01, 8.233700e-02, 3.619154e-01), "X", 
  lattice.add_site = (2.964000e-01, 2.229710e-01, 9.700446e-01), "X", 
  lattice.add_site = (7.036000e-01, 3.780290e-01, 4.730446e-01), "X", 
  lattice.add_site = (7.964000e-01, 3.780290e-01, 9.700446e-01), "X", 
  lattice.add_site = (2.036000e-01, 2.229710e-01, 4.730446e-01), "X", 
  lattice.add_site = (7.433000e-01, 2.187640e-01, 7.282044e-01), "X", 
  lattice.add_site = (2.567000e-01, 3.822360e-01, 2.312044e-01), "X", 
  lattice.add_site = (2.433000e-01, 3.822360e-01, 7.282044e-01), "X", 
  lattice.add_site = (7.567000e-01, 2.187640e-01, 2.312044e-01), "X", 
  lattice.find_space_group()
  return lattice


def s2():
  """ Returns s2 lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (9.870000e-01, 0.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 2.950000e-01, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 1.000000e+00)
  lattice.name = "s2"
  lattice.add_site = (8.808975e-01, 7.375000e-02, 5.834000e-01), "A", 
  lattice.add_site = (5.996025e-01, 2.212500e-01, 8.340000e-02), "A", 
  lattice.add_site = (1.061025e-01, 2.212500e-01, 4.166000e-01), "A", 
  lattice.add_site = (3.873975e-01, 7.375000e-02, 9.166000e-01), "A", 
  lattice.add_site = (3.406137e-01, 7.375000e-02, 2.890000e-01), "A", 
  lattice.add_site = (1.528863e-01, 2.212500e-01, 7.890000e-01), "A", 
  lattice.add_site = (6.463863e-01, 2.212500e-01, 7.110000e-01), "A", 
  lattice.add_site = (8.341137e-01, 7.375000e-02, 2.110000e-01), "A", 
  lattice.add_site = (1.302840e-01, 7.375000e-02, 8.420000e-02), "B", 
  lattice.add_site = (3.632160e-01, 2.212500e-01, 5.842000e-01), "B", 
  lattice.add_site = (8.567160e-01, 2.212500e-01, 9.158000e-01), "B", 
  lattice.add_site = (6.237840e-01, 7.375000e-02, 4.158000e-01), "B", 
  lattice.add_site = (7.718340e-01, 7.375000e-02, 7.620000e-01), "X", 
  lattice.add_site = (7.086660e-01, 2.212500e-01, 2.620000e-01), "X", 
  lattice.add_site = (2.151660e-01, 2.212500e-01, 2.380000e-01), "X", 
  lattice.add_site = (2.783340e-01, 7.375000e-02, 7.380000e-01), "X", 
  lattice.add_site = (2.408280e-01, 7.375000e-02, 4.740000e-01), "X", 
  lattice.add_site = (2.526720e-01, 2.212500e-01, 9.740000e-01), "X", 
  lattice.add_site = (7.461720e-01, 2.212500e-01, 5.260000e-01), "X", 
  lattice.add_site = (7.343280e-01, 7.375000e-02, 2.600000e-02), "X", 
  lattice.add_site = (5.221230e-01, 7.375000e-02, 6.190000e-01), "X", 
  lattice.add_site = (9.583770e-01, 2.212500e-01, 1.190000e-01), "X", 
  lattice.add_site = (4.648770e-01, 2.212500e-01, 3.810000e-01), "X", 
  lattice.add_site = (2.862300e-02, 7.375000e-02, 8.810000e-01), "X", 
  lattice.add_site = (4.599420e-01, 7.375000e-02, 1.190000e-01), "X", 
  lattice.add_site = (3.355800e-02, 2.212500e-01, 6.190000e-01), "X", 
  lattice.add_site = (5.270580e-01, 2.212500e-01, 8.810000e-01), "X", 
  lattice.add_site = (9.534420e-01, 7.375000e-02, 3.810000e-01), "X", 
  lattice.find_space_group()
  return lattice


def s3():
  """ Returns s3 lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (6.110000e-01, 0.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 1.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 3.470000e-01)
  lattice.name = "s3"
  lattice.add_site = (4.704700e-02, 3.190000e-01, 1.735000e-01), "A", 
  lattice.add_site = (5.639530e-01, 6.810000e-01, 1.735000e-01), "A", 
  lattice.add_site = (2.584530e-01, 8.190000e-01, 1.735000e-01), "A", 
  lattice.add_site = (3.525470e-01, 1.810000e-01, 1.735000e-01), "A", 
  lattice.add_site = (0.000000e+00, 0.000000e+00, 0.000000e+00), "B", 
  lattice.add_site = (3.055000e-01, 5.000000e-01, 0.000000e+00), "B", 
  lattice.add_site = (1.344200e-01, 5.000000e-02, 1.735000e-01), "X", 
  lattice.add_site = (4.765800e-01, 9.500000e-01, 1.735000e-01), "X", 
  lattice.add_site = (1.710800e-01, 5.500000e-01, 1.735000e-01), "X", 
  lattice.add_site = (4.399200e-01, 4.500000e-01, 1.735000e-01), "X", 
  lattice.add_site = (2.199600e-01, 3.100000e-01, 0.000000e+00), "X", 
  lattice.add_site = (3.910400e-01, 6.900000e-01, 0.000000e+00), "X", 
  lattice.add_site = (8.554000e-02, 8.100000e-01, 0.000000e+00), "X", 
  lattice.add_site = (5.254600e-01, 1.900000e-01, 0.000000e+00), "X", 
  lattice.find_space_group()
  return lattice


def b10I():
  """ Returns b10I lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (1.000000e+00, 0.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 5.850000e-01, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 5.220000e-01)
  lattice.name = "b10I"
  lattice.add_site = (9.294000e-02, 1.462500e-01, 2.262713e-01), "A", 
  lattice.add_site = (4.070600e-01, 4.387500e-01, 4.872713e-01), "A", 
  lattice.add_site = (9.070600e-01, 4.387500e-01, 2.957287e-01), "A", 
  lattice.add_site = (5.929400e-01, 1.462500e-01, 3.472866e-02), "A", 
  lattice.add_site = (0.000000e+00, 0.000000e+00, 0.000000e+00), "A", 
  lattice.add_site = (5.000000e-01, 0.000000e+00, 2.610000e-01), "A", 
  lattice.add_site = (0.000000e+00, 2.925000e-01, 0.000000e+00), "A", 
  lattice.add_site = (5.000000e-01, 2.925000e-01, 2.610000e-01), "A", 
  lattice.add_site = (2.731900e-01, 1.462500e-01, 5.188941e-01), "B", 
  lattice.add_site = (2.268100e-01, 4.387500e-01, 2.578941e-01), "B", 
  lattice.add_site = (7.268100e-01, 4.387500e-01, 3.105900e-03), "B", 
  lattice.add_site = (7.731900e-01, 1.462500e-01, 2.641059e-01), "B", 
  lattice.add_site = (9.051000e-02, 1.462500e-01, 4.124635e-01), "X", 
  lattice.add_site = (4.094900e-01, 4.387500e-01, 1.514635e-01), "X", 
  lattice.add_site = (9.094900e-01, 4.387500e-01, 1.095365e-01), "X", 
  lattice.add_site = (5.905100e-01, 1.462500e-01, 3.705365e-01), "X", 
  lattice.add_site = (4.334300e-01, 1.462500e-01, 1.257863e-01), "X", 
  lattice.add_site = (6.657000e-02, 4.387500e-01, 3.867863e-01), "X", 
  lattice.add_site = (5.665700e-01, 4.387500e-01, 3.962137e-01), "X", 
  lattice.add_site = (9.334300e-01, 1.462500e-01, 1.352137e-01), "X", 
  lattice.add_site = (1.631800e-01, 1.005030e-02, 1.349370e-01), "X", 
  lattice.add_site = (3.368200e-01, 5.749497e-01, 3.959370e-01), "X", 
  lattice.add_site = (8.368200e-01, 3.025503e-01, 3.870630e-01), "X", 
  lattice.add_site = (6.631800e-01, 2.824497e-01, 1.260630e-01), "X", 
  lattice.add_site = (8.368200e-01, 5.749497e-01, 3.870630e-01), "X", 
  lattice.add_site = (6.631800e-01, 1.005030e-02, 1.260630e-01), "X", 
  lattice.add_site = (1.631800e-01, 2.824497e-01, 1.349370e-01), "X", 
  lattice.add_site = (3.368200e-01, 3.025503e-01, 3.959370e-01), "X", 
  lattice.find_space_group()
  return lattice


def b1I():
  """ Returns b1I lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (2.880000e-01, 0.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 2.880000e-01, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 1.000000e+00)
  lattice.name = "b1I"
  lattice.add_site = (0.000000e+00, 0.000000e+00, 0.000000e+00), "A", 
  lattice.add_site = (1.440000e-01, 1.440000e-01, 5.000000e-01), "A", 
  lattice.add_site = (0.000000e+00, 0.000000e+00, 3.500000e-01), "A", 
  lattice.add_site = (0.000000e+00, 0.000000e+00, -3.500000e-01), "A", 
  lattice.add_site = (1.440000e-01, 1.440000e-01, 8.500000e-01), "B", 
  lattice.add_site = (-1.440000e-01, -1.440000e-01, -8.500000e-01), "B", 
  lattice.add_site = (0.000000e+00, 1.440000e-01, 0.000000e+00), "X", 
  lattice.add_site = (1.440000e-01, 0.000000e+00, 0.000000e+00), "X", 
  lattice.add_site = (1.440000e-01, 0.000000e+00, 5.000000e-01), "X", 
  lattice.add_site = (0.000000e+00, 1.440000e-01, 5.000000e-01), "X", 
  lattice.add_site = (0.000000e+00, 0.000000e+00, 1.500000e-01), "X", 
  lattice.add_site = (0.000000e+00, 0.000000e+00, -1.500000e-01), "X", 
  lattice.add_site = (1.440000e-01, 1.440000e-01, 6.500000e-01), "X", 
  lattice.add_site = (-1.440000e-01, -1.440000e-01, -6.500000e-01), "X", 
  lattice.find_space_group()
  return lattice


def b2I():
  """ Returns b2I lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (1.000000e+00, 0.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 1.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 5.900000e-01)
  lattice.name = "b2I"
  lattice.add_site = (0.000000e+00, 0.000000e+00, 0.000000e+00), "A", 
  lattice.add_site = (0.000000e+00, 5.000000e-01, 2.950000e-01), "A", 
  lattice.add_site = (5.000000e-01, 0.000000e+00, 2.950000e-01), "B", 
  lattice.add_site = (2.324500e-01, 2.324500e-01, 0.000000e+00), "X", 
  lattice.add_site = (7.675500e-01, 7.675500e-01, 0.000000e+00), "X", 
  lattice.add_site = (7.675500e-01, 2.324500e-01, 0.000000e+00), "X", 
  lattice.add_site = (2.324500e-01, 7.675500e-01, 0.000000e+00), "X", 
  lattice.find_space_group()
  return lattice


def b4I():
  """ Returns b4I lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (1.000000e+00, 0.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 1.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 1.000000e+00)
  lattice.name = "b4I"
  lattice.add_site = (0.000000e+00, 0.000000e+00, 0.000000e+00), "A", 
  lattice.add_site = (0.000000e+00, 5.000000e-01, 5.000000e-01), "A", 
  lattice.add_site = (5.000000e-01, 0.000000e+00, 5.000000e-01), "B", 
  lattice.add_site = (2.700000e-01, 2.700000e-01, 2.250000e-01), "X", 
  lattice.add_site = (-2.700000e-01, -2.700000e-01, 2.250000e-01), "X", 
  lattice.add_site = (-2.700000e-01, 2.700000e-01, -2.250000e-01), "X", 
  lattice.add_site = (2.700000e-01, -2.700000e-01, -2.250000e-01), "X", 
  lattice.find_space_group()
  return lattice


def b5I(u=0.25):
  """ Returns b5I lattice."""
  from . import Lattice
  
  lattice = b5(u)
  nb_B = 2
  for site in lattice.sites:
    if "A" in site.type and nb_B > 0: 
      site.type = ["B"]
      nb_B -= 1
    elif "B" in site.type: site.type = ["A"]
  lattice.name = "b5I"
  lattice.find_space_group()
  return lattice


def b9I():
  """ Returns b9I lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (8.620000e-01, 0.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 1.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 2.820000e-01)
  lattice.name = "b9I"
  lattice.add_site = (6.516720e-01, 6.540000e-01, 7.050000e-02), "A", 
  lattice.add_site = (-6.516720e-01, -6.540000e-01, -7.050000e-02), "A", 
  lattice.add_site = (1.082672e+00, -1.540000e-01, 7.050000e-02), "A", 
  lattice.add_site = (-1.082672e+00, 1.540000e-01, -7.050000e-02), "A", 
  lattice.add_site = (3.732460e-01, 6.100000e-01, 7.050000e-02), "A", 
  lattice.add_site = (-3.732460e-01, -6.100000e-01, -7.050000e-02), "A", 
  lattice.add_site = (8.042460e-01, -1.100000e-01, 7.050000e-02), "A", 
  lattice.add_site = (-8.042460e-01, 1.100000e-01, -7.050000e-02), "A", 
  lattice.add_site = (3.620400e-01, 1.080000e-01, 7.050000e-02), "B", 
  lattice.add_site = (-3.620400e-01, -1.080000e-01, -7.050000e-02), "B", 
  lattice.add_site = (7.930400e-01, 3.920000e-01, 7.050000e-02), "B", 
  lattice.add_site = (-7.930400e-01, -3.920000e-01, -7.050000e-02), "B", 
  lattice.add_site = (1.792960e-01, 1.620000e-01, 7.050000e-02), "X", 
  lattice.add_site = (-1.792960e-01, -1.620000e-01, -7.050000e-02), "X", 
  lattice.add_site = (6.102960e-01, 3.380000e-01, 7.050000e-02), "X", 
  lattice.add_site = (-6.102960e-01, -3.380000e-01, -7.050000e-02), "X", 
  lattice.add_site = (9.913000e-02, 4.770000e-01, 7.050000e-02), "X", 
  lattice.add_site = (-9.913000e-02, -4.770000e-01, -7.050000e-02), "X", 
  lattice.add_site = (5.301300e-01, 2.300000e-02, 7.050000e-02), "X", 
  lattice.add_site = (-5.301300e-01, -2.300000e-02, -7.050000e-02), "X", 
  lattice.add_site = (4.491020e-01, 7.840000e-01, 7.050000e-02), "X", 
  lattice.add_site = (-4.491020e-01, -7.840000e-01, -7.050000e-02), "X", 
  lattice.add_site = (8.801020e-01, -2.840000e-01, 7.050000e-02), "X", 
  lattice.add_site = (-8.801020e-01, 2.840000e-01, -7.050000e-02), "X", 
  lattice.add_site = (3.611780e-01, 4.240000e-01, 7.050000e-02), "X", 
  lattice.add_site = (-3.611780e-01, -4.240000e-01, -7.050000e-02), "X", 
  lattice.add_site = (7.921780e-01, 7.600000e-02, 7.050000e-02), "X", 
  lattice.add_site = (-7.921780e-01, -7.600000e-02, -7.050000e-02), "X", 
  lattice.find_space_group()
  return lattice


def d1I():
  """ Returns d1I lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (1.000000e+00, 0.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 1.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 7.450000e-01)
  lattice.name = "d1I"
  lattice.add_site = (0.000000e+00, 5.000000e-01, 1.862500e-01), "A", 
  lattice.add_site = (5.000000e-01, 0.000000e+00, 5.587500e-01), "A", 
  lattice.add_site = (0.000000e+00, 5.000000e-01, 5.587500e-01), "A", 
  lattice.add_site = (5.000000e-01, 0.000000e+00, 1.862500e-01), "A", 
  lattice.add_site = (1.400000e-01, 1.630000e-01, 0.000000e+00), "A", 
  lattice.add_site = (8.600000e-01, 8.370000e-01, 0.000000e+00), "A", 
  lattice.add_site = (8.370000e-01, 1.400000e-01, 3.725000e-01), "A", 
  lattice.add_site = (1.630000e-01, 8.600000e-01, 3.725000e-01), "A", 
  lattice.add_site = (3.600000e-01, 6.630000e-01, 0.000000e+00), "B", 
  lattice.add_site = (6.400000e-01, 3.370000e-01, 0.000000e+00), "B", 
  lattice.add_site = (6.630000e-01, 6.400000e-01, 3.725000e-01), "B", 
  lattice.add_site = (3.370000e-01, 3.600000e-01, 3.725000e-01), "B", 
  lattice.add_site = (6.710000e-01, 1.710000e-01, 1.862500e-01), "X", 
  lattice.add_site = (3.290000e-01, 8.290000e-01, 1.862500e-01), "X", 
  lattice.add_site = (8.290000e-01, 6.710000e-01, 5.587500e-01), "X", 
  lattice.add_site = (1.710000e-01, 3.290000e-01, 5.587500e-01), "X", 
  lattice.add_site = (3.290000e-01, 8.290000e-01, 5.587500e-01), "X", 
  lattice.add_site = (6.710000e-01, 1.710000e-01, 5.587500e-01), "X", 
  lattice.add_site = (1.710000e-01, 3.290000e-01, 1.862500e-01), "X", 
  lattice.add_site = (8.290000e-01, 6.710000e-01, 1.862500e-01), "X", 
  lattice.add_site = (9.600000e-02, 6.370000e-01, 0.000000e+00), "X", 
  lattice.add_site = (9.040000e-01, 3.630000e-01, 0.000000e+00), "X", 
  lattice.add_site = (3.630000e-01, 9.600000e-02, 3.725000e-01), "X", 
  lattice.add_site = (6.370000e-01, 9.040000e-01, 3.725000e-01), "X", 
  lattice.add_site = (4.040000e-01, 1.370000e-01, 0.000000e+00), "X", 
  lattice.add_site = (5.960000e-01, 8.630000e-01, 0.000000e+00), "X", 
  lattice.add_site = (1.370000e-01, 5.960000e-01, 3.725000e-01), "X", 
  lattice.add_site = (8.630000e-01, 4.040000e-01, 3.725000e-01), "X", 
  lattice.find_space_group()
  return lattice


def d3I():
  """ Returns d3I lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (1.000000e+00, 0.000000e+00, -7.400000e-01),\
                     (0.000000e+00, 2.990000e-01, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 4.380000e-01)
  lattice.name = "d3I"
  lattice.add_site = (0.000000e+00, 0.000000e+00, 0.000000e+00), "A", 
  lattice.add_site = (5.000000e-01, 1.495000e-01, 0.000000e+00), "A", 
  lattice.add_site = (-1.850000e-01, 0.000000e+00, 1.095000e-01), "A", 
  lattice.add_site = (1.850000e-01, 0.000000e+00, -1.095000e-01), "A", 
  lattice.add_site = (3.150000e-01, 1.495000e-01, 1.095000e-01), "B", 
  lattice.add_site = (-3.150000e-01, -1.495000e-01, -1.095000e-01), "B", 
  lattice.add_site = (-5.920000e-03, 0.000000e+00, 2.006040e-01), "X", 
  lattice.add_site = (5.920000e-03, 0.000000e+00, -2.006040e-01), "X", 
  lattice.add_site = (4.940800e-01, 1.495000e-01, 2.006040e-01), "X", 
  lattice.add_site = (-4.940800e-01, -1.495000e-01, -2.006040e-01), "X", 
  lattice.add_site = (-3.759200e-01, 0.000000e+00, 4.196040e-01), "X", 
  lattice.add_site = (3.759200e-01, 0.000000e+00, -4.196040e-01), "X", 
  lattice.add_site = (1.240800e-01, 1.495000e-01, 4.196040e-01), "X", 
  lattice.add_site = (-1.240800e-01, -1.495000e-01, -4.196040e-01), "X", 
  lattice.find_space_group()
  return lattice


def s1I():
  """ Returns s1I lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (1.000000e+00, 0.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 6.010000e-01, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 9.940000e-01)
  lattice.name = "s1I"
  lattice.add_site = (3.798000e-01, 5.336880e-01, 1.851822e-01), "A", 
  lattice.add_site = (6.202000e-01, 6.731200e-02, 6.821822e-01), "A", 
  lattice.add_site = (8.798000e-01, 6.731200e-02, 1.851822e-01), "A", 
  lattice.add_site = (1.202000e-01, 5.336880e-01, 6.821822e-01), "A", 
  lattice.add_site = (3.966000e-01, 2.325870e-01, 8.082214e-01), "A", 
  lattice.add_site = (6.034000e-01, 3.684130e-01, 3.112214e-01), "A", 
  lattice.add_site = (8.966000e-01, 3.684130e-01, 8.082214e-01), "A", 
  lattice.add_site = (1.034000e-01, 2.325870e-01, 3.112214e-01), "A", 
  lattice.add_site = (3.802000e-01, 5.356112e-01, 8.081220e-01), "A", 
  lattice.add_site = (6.198000e-01, 6.538880e-02, 3.111220e-01), "A", 
  lattice.add_site = (8.802000e-01, 6.538880e-02, 8.081220e-01), "A", 
  lattice.add_site = (1.198000e-01, 5.356112e-01, 3.111220e-01), "A", 
  lattice.add_site = (6.335000e-01, 5.191438e-01, 1.143100e-02), "A", 
  lattice.add_site = (3.665000e-01, 8.185620e-02, 5.084310e-01), "A", 
  lattice.add_site = (1.335000e-01, 8.185620e-02, 1.143100e-02), "A", 
  lattice.add_site = (8.665000e-01, 5.191438e-01, 5.084310e-01), "A", 
  lattice.add_site = (3.815000e-01, 2.279593e-01, 1.688806e-01), "B", 
  lattice.add_site = (6.185000e-01, 3.730407e-01, 6.658806e-01), "B", 
  lattice.add_site = (8.815000e-01, 3.730407e-01, 1.688806e-01), "B", 
  lattice.add_site = (1.185000e-01, 2.279593e-01, 6.658806e-01), "B", 
  lattice.add_site = (6.391000e-01, 2.193650e-01, 9.861474e-01), "B", 
  lattice.add_site = (3.609000e-01, 3.816350e-01, 4.891474e-01), "B", 
  lattice.add_site = (1.391000e-01, 3.816350e-01, 9.861474e-01), "B", 
  lattice.add_site = (8.609000e-01, 2.193650e-01, 4.891474e-01), "B", 
  lattice.add_site = (5.386000e-01, 3.593980e-01, 1.214668e-01), "X", 
  lattice.add_site = (4.614000e-01, 2.416020e-01, 6.184668e-01), "X", 
  lattice.add_site = (3.860000e-02, 2.416020e-01, 1.214668e-01), "X", 
  lattice.add_site = (9.614000e-01, 3.593980e-01, 6.184668e-01), "X", 
  lattice.add_site = (5.201000e-01, 3.618020e-01, 8.572256e-01), "X", 
  lattice.add_site = (4.799000e-01, 2.391980e-01, 3.602256e-01), "X", 
  lattice.add_site = (2.010000e-02, 2.391980e-01, 8.572256e-01), "X", 
  lattice.add_site = (9.799000e-01, 3.618020e-01, 3.602256e-01), "X", 
  lattice.add_site = (2.863000e-01, 5.288800e-01, 3.479000e-03), "X", 
  lattice.add_site = (7.137000e-01, 7.212000e-02, 5.004790e-01), "X", 
  lattice.add_site = (7.863000e-01, 7.212000e-02, 3.479000e-03), "X", 
  lattice.add_site = (2.137000e-01, 5.288800e-01, 5.004790e-01), "X", 
  lattice.add_site = (7.693000e-01, 5.336279e-01, 7.267134e-01), "X", 
  lattice.add_site = (2.307000e-01, 6.737210e-02, 2.297134e-01), "X", 
  lattice.add_site = (2.693000e-01, 6.737210e-02, 7.267134e-01), "X", 
  lattice.add_site = (7.307000e-01, 5.336279e-01, 2.297134e-01), "X", 
  lattice.add_site = (5.111000e-01, 6.851400e-02, 1.245482e-01), "X", 
  lattice.add_site = (4.889000e-01, 5.324860e-01, 6.215482e-01), "X", 
  lattice.add_site = (1.110000e-02, 5.324860e-01, 1.245482e-01), "X", 
  lattice.add_site = (9.889000e-01, 6.851400e-02, 6.215482e-01), "X", 
  lattice.add_site = (5.427000e-01, 8.233700e-02, 8.589154e-01), "X", 
  lattice.add_site = (4.573000e-01, 5.186630e-01, 3.619154e-01), "X", 
  lattice.add_site = (4.270000e-02, 5.186630e-01, 8.589154e-01), "X", 
  lattice.add_site = (9.573000e-01, 8.233700e-02, 3.619154e-01), "X", 
  lattice.add_site = (2.964000e-01, 2.229710e-01, 9.700446e-01), "X", 
  lattice.add_site = (7.036000e-01, 3.780290e-01, 4.730446e-01), "X", 
  lattice.add_site = (7.964000e-01, 3.780290e-01, 9.700446e-01), "X", 
  lattice.add_site = (2.036000e-01, 2.229710e-01, 4.730446e-01), "X", 
  lattice.add_site = (7.433000e-01, 2.187640e-01, 7.282044e-01), "X", 
  lattice.add_site = (2.567000e-01, 3.822360e-01, 2.312044e-01), "X", 
  lattice.add_site = (2.433000e-01, 3.822360e-01, 7.282044e-01), "X", 
  lattice.add_site = (7.567000e-01, 2.187640e-01, 2.312044e-01), "X", 
  lattice.find_space_group()
  return lattice


def s2I():
  """ Returns s2I lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (9.870000e-01, 0.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 2.950000e-01, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 1.000000e+00)
  lattice.name = "s2I"
  lattice.add_site = (1.302840e-01, 7.375000e-02, 8.420000e-02), "A", 
  lattice.add_site = (3.632160e-01, 2.212500e-01, 5.842000e-01), "A", 
  lattice.add_site = (8.567160e-01, 2.212500e-01, 9.158000e-01), "A", 
  lattice.add_site = (6.237840e-01, 7.375000e-02, 4.158000e-01), "A", 
  lattice.add_site = (8.808975e-01, 7.375000e-02, 5.834000e-01), "A", 
  lattice.add_site = (5.996025e-01, 2.212500e-01, 8.340000e-02), "A", 
  lattice.add_site = (1.061025e-01, 2.212500e-01, 4.166000e-01), "A", 
  lattice.add_site = (3.873975e-01, 7.375000e-02, 9.166000e-01), "A", 
  lattice.add_site = (3.406137e-01, 7.375000e-02, 2.890000e-01), "B", 
  lattice.add_site = (1.528863e-01, 2.212500e-01, 7.890000e-01), "B", 
  lattice.add_site = (6.463863e-01, 2.212500e-01, 7.110000e-01), "B", 
  lattice.add_site = (8.341137e-01, 7.375000e-02, 2.110000e-01), "B", 
  lattice.add_site = (7.718340e-01, 7.375000e-02, 7.620000e-01), "X", 
  lattice.add_site = (7.086660e-01, 2.212500e-01, 2.620000e-01), "X", 
  lattice.add_site = (2.151660e-01, 2.212500e-01, 2.380000e-01), "X", 
  lattice.add_site = (2.783340e-01, 7.375000e-02, 7.380000e-01), "X", 
  lattice.add_site = (2.408280e-01, 7.375000e-02, 4.740000e-01), "X", 
  lattice.add_site = (2.526720e-01, 2.212500e-01, 9.740000e-01), "X", 
  lattice.add_site = (7.461720e-01, 2.212500e-01, 5.260000e-01), "X", 
  lattice.add_site = (7.343280e-01, 7.375000e-02, 2.600000e-02), "X", 
  lattice.add_site = (5.221230e-01, 7.375000e-02, 6.190000e-01), "X", 
  lattice.add_site = (9.583770e-01, 2.212500e-01, 1.190000e-01), "X", 
  lattice.add_site = (4.648770e-01, 2.212500e-01, 3.810000e-01), "X", 
  lattice.add_site = (2.862300e-02, 7.375000e-02, 8.810000e-01), "X", 
  lattice.add_site = (4.599420e-01, 7.375000e-02, 1.190000e-01), "X", 
  lattice.add_site = (3.355800e-02, 2.212500e-01, 6.190000e-01), "X", 
  lattice.add_site = (5.270580e-01, 2.212500e-01, 8.810000e-01), "X", 
  lattice.add_site = (9.534420e-01, 7.375000e-02, 3.810000e-01), "X", 
  lattice.find_space_group()
  return lattice

def b16():
  """ Returns a b16 Lattice. """
  from . import Lattice
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (1.000000e+00, 0.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 6.000000e-01, 0.000000e+00),\
                     (-3.144352e-01, 0.000000e+00, 9.735098e-01)
  lattice.name = "b16"
  lattice.add_site = (2.050000e-01, 3.492000e-01, 3.006070e-01), "A", 
  lattice.add_site = (-2.050000e-01, -3.492000e-01, -3.006070e-01), "A", 
  lattice.add_site = (7.050000e-01, -4.920000e-02, 1.433895e-01), "A", 
  lattice.add_site = (-7.050000e-01, 4.920000e-02, -1.433895e-01), "A", 
  lattice.add_site = (4.550000e-01, 3.492000e-01, -2.137926e-02), "A", 
  lattice.add_site = (-4.550000e-01, -3.492000e-01, 2.137926e-02), "A", 
  lattice.add_site = (9.550000e-01, -4.920000e-02, -1.785968e-01), "A", 
  lattice.add_site = (-9.550000e-01, 4.920000e-02, 1.785968e-01), "A", 
  lattice.add_site = (3.190000e-01, 3.660000e-02, 1.430726e-01), "B", 
  lattice.add_site = (-3.190000e-01, -3.660000e-02, -1.430726e-01), "B", 
  lattice.add_site = (8.190000e-01, 2.634000e-01, -1.414493e-02), "B", 
  lattice.add_site = (-8.190000e-01, -2.634000e-01, 1.414493e-02), "B", 
  lattice.add_site = (4.920000e-01, -2.520000e-02, 8.867536e-02), "X", 
  lattice.add_site = (-4.920000e-01, 2.520000e-02, -8.867536e-02), "X", 
  lattice.add_site = (9.920000e-01, 3.252000e-01, -6.854222e-02), "X", 
  lattice.add_site = (-9.920000e-01, -3.252000e-01, 6.854222e-02), "X", 
  lattice.add_site = (1.860000e-01, -2.520000e-02, 3.789253e-02), "X", 
  lattice.add_site = (-1.860000e-01, 2.520000e-02, -3.789253e-02), "X", 
  lattice.add_site = (6.860000e-01, 3.252000e-01, -1.193250e-01), "X", 
  lattice.add_site = (-6.860000e-01, -3.252000e-01, 1.193250e-01), "X", 
  lattice.add_site = (2.800000e-01, -2.520000e-02, 3.023356e-01), "X", 
  lattice.add_site = (-2.800000e-01, 2.520000e-02, -3.023356e-01), "X", 
  lattice.add_site = (7.800000e-01, 3.252000e-01, 1.451181e-01), "X", 
  lattice.add_site = (-7.800000e-01, -3.252000e-01, -1.451181e-01), "X", 
  lattice.add_site = (3.190000e-01, 2.202000e-01, 1.430726e-01), "X", 
  lattice.add_site = (-3.190000e-01, -2.202000e-01, -1.430726e-01), "X", 
  lattice.add_site = (8.190000e-01, 7.980000e-02, -1.414493e-02), "X", 
  lattice.add_site = (-8.190000e-01, -7.980000e-02, 1.414493e-02), "X", 
  lattice.find_space_group()
  return lattice



def b15():
  """ Returns a b15 Lattice. """
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (1.000000e+00, 0.000000e+00, 0.000000e+00),\
  		   (0.000000e+00, 9.754710e-01, 0.000000e+00),\
  		   (0.000000e+00, 0.000000e+00, 5.509780e-01)
  lattice.name = "b15"
  lattice.add_site = (4.540000e-01, 4.672506e-01, 3.024870e-01), "A", 
  lattice.add_site = (4.600000e-02, 5.082204e-01, 2.699792e-02), "A", 
  lattice.add_site = (5.460000e-01, 9.549861e-01, 5.239801e-01), "A", 
  lattice.add_site = (9.540000e-01, 2.048489e-02, 2.484911e-01), "A", 
  lattice.add_site = (2.223000e-01, 1.902168e-01, 1.101956e-02), "A", 
  lattice.add_site = (2.777000e-01, 7.852542e-01, 2.865086e-01), "A", 
  lattice.add_site = (7.777000e-01, 6.779524e-01, 2.644695e-01), "A", 
  lattice.add_site = (7.223000e-01, 2.975187e-01, 5.399585e-01), "A", 
  lattice.add_site = (3.670000e-01, 2.038734e-01, 3.179143e-01), "B", 
  lattice.add_site = (1.330000e-01, 7.715976e-01, 4.242531e-02), "B", 
  lattice.add_site = (6.330000e-01, 6.916089e-01, 5.085527e-01), "B", 
  lattice.add_site = (8.670000e-01, 2.838621e-01, 2.330637e-01), "B", 
  lattice.add_site = (5.170000e-01, 1.277868e-01, 3.945003e-01), "X", 
  lattice.add_site = (9.830000e-01, 8.476843e-01, 1.190113e-01), "X", 
  lattice.add_site = (4.830000e-01, 6.155223e-01, 4.319668e-01), "X", 
  lattice.add_site = (1.700000e-02, 3.599488e-01, 1.564777e-01), "X", 
  lattice.add_site = (2.890000e-01, 1.082773e-01, 1.680484e-01), "X", 
  lattice.add_site = (2.110000e-01, 8.671937e-01, 4.435374e-01), "X", 
  lattice.add_site = (7.110000e-01, 5.960128e-01, 1.074407e-01), "X", 
  lattice.add_site = (7.890000e-01, 3.794582e-01, 3.829297e-01), "X", 
  lattice.add_site = (2.310000e-01, 2.594753e-01, 4.440883e-01), "X", 
  lattice.add_site = (2.690000e-01, 7.159957e-01, 1.685993e-01), "X", 
  lattice.add_site = (7.690000e-01, 7.472108e-01, 3.823787e-01), "X", 
  lattice.add_site = (7.310000e-01, 2.282602e-01, 1.068897e-01), "X", 
  lattice.add_site = (4.130000e-01, 3.804337e-01, 2.259010e-01), "X", 
  lattice.add_site = (8.700000e-02, 5.950373e-01, 5.013900e-01), "X", 
  lattice.add_site = (5.870000e-01, 8.681692e-01, 4.958802e-02), "X", 
  lattice.add_site = (9.130000e-01, 1.073018e-01, 3.250770e-01), "X"
  lattice.find_space_group()
  return lattice

def s3I():
  """ Returns s3I lattice."""
  from . import Lattice
  
  lattice = Lattice()
  lattice.scale = 1.000000e+00
  lattice.set_cell = (6.110000e-01, 0.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 1.000000e+00, 0.000000e+00),\
                     (0.000000e+00, 0.000000e+00, 3.470000e-01)
  lattice.name = "s3I"
  lattice.add_site = (0.000000e+00, 0.000000e+00, 0.000000e+00), "A", 
  lattice.add_site = (3.055000e-01, 5.000000e-01, 0.000000e+00), "A", 
  lattice.add_site = (4.704700e-02, 3.190000e-01, 1.735000e-01), "A", 
  lattice.add_site = (5.639530e-01, 6.810000e-01, 1.735000e-01), "A", 
  lattice.add_site = (2.584530e-01, 8.190000e-01, 1.735000e-01), "B", 
  lattice.add_site = (3.525470e-01, 1.810000e-01, 1.735000e-01), "B", 
  lattice.add_site = (1.344200e-01, 5.000000e-02, 1.735000e-01), "X", 
  lattice.add_site = (4.765800e-01, 9.500000e-01, 1.735000e-01), "X", 
  lattice.add_site = (1.710800e-01, 5.500000e-01, 1.735000e-01), "X", 
  lattice.add_site = (4.399200e-01, 4.500000e-01, 1.735000e-01), "X", 
  lattice.add_site = (2.199600e-01, 3.100000e-01, 0.000000e+00), "X", 
  lattice.add_site = (3.910400e-01, 6.900000e-01, 0.000000e+00), "X", 
  lattice.add_site = (8.554000e-02, 8.100000e-01, 0.000000e+00), "X", 
  lattice.add_site = (5.254600e-01, 1.900000e-01, 0.000000e+00), "X", 
  lattice.find_space_group()
  return lattice


