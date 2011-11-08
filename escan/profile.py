""" Computes real-space profiles of the wavefunctions. """
__docformat__ = "restructuredtext en"

def linear_profile(extract, direction=(0,0,1), nbpoints=20, sigma=None, indices=None):
  """ Computes profile for given direction of each wavefunction.
  
      :Parameters:
        extract 
          Extraction object from an escan calculation. 
        nbpoints 
          Number of points to consider along the direction.
        sigma 
          Smearing to use when creating profile. If None, defaults to 3/4 of
          distance between two points along the profile.
        indices 
          Indices of the wavefunctions for which to compute the profile.
          Can also be None, in which case the profile of each wavefunction is
          computed.

      To get nice plots, we compute the convolution of the real-space density
      and a gaussian in fourier space. Comes down to an interpolation in
      real-space.
  """
  from numpy import dot, max, arange, array, tensordot, exp, multiply, sum
  from numpy.linalg import inv, norm
  from quantities import angstrom
  from lada.physics import a0

  direction = array(direction, dtype="float64")
  assert norm(direction) > 1e-12, ValueError("Direction cannot be null.")
  assert hasattr(extract, "gwfns"), ValueError("extract does not seem to be an escan extraction object.")
  if indices is None: indices = range(extract.gwfns)
  assert nbpoints > 1, ValueError("The number of points should be strictly larger than 1.")

  # first computes intersection of direction and cell.  At the end, we should
  # end up in a0, knowing the structure cell should be in angstrom.
  udir = direction / norm(direction)
  dir = dot(inv(extract.structure.cell), direction)
  direction = direction/max(abs(dir)) * extract.structure.scale * angstrom.rescale(a0)

  # then creates array of points along direction.
  zpoints = arange(nbpoints+1, dtype="float64").reshape(nbpoints+1, 1) / float(nbpoints)
  zpoints = multiply(zpoints, direction.magnitude) * direction.units
  x =  sum(zpoints*udir, axis=1).rescale(angstrom) 

  # computes sigma in g-space.
  if sigma is None: sigma = 0.75 * norm(zpoints[1] - zpoints[0]) * zpoints.units
  if not hasattr(sigma, "units"): sigma *= angstrom.rescale(a0)
  gsigma = 1e0/sigma

  # finds gvectors with zero components in directions perpendicular to direction.
  # dot does not work well with units.
  gvectors = extract.gvectors
  nonzero = (gvectors - dot(gvectors.magnitude, udir).reshape(gvectors.shape[0], 1) * udir * gvectors.units)**2
  nonzero = sum(nonzero, axis=1) < 1e-12
  gvectors = extract.gvectors[nonzero, :]

  # computes all exponentials exp(-i r.g), with r in first dim, and g in second.
  # tensordot does not conserve units. Must do it by hand.
  units = (zpoints.units * gvectors.units).simplified
  translations = exp(1j * tensordot(zpoints, gvectors, ((1),(1))) * units)

  # creates gaussian x translations matrix.
  gpoints = (dot(gvectors.magnitude, udir) * gvectors.units / gsigma).simplified
  gaussians = multiply(exp(-gpoints**2), translations)

  # fourrier operator to get auto-correlation function in g-space.
  units = (extract.rvectors.units * gvectors.units).simplified
  fourier = exp(-1j * tensordot(extract.rvectors, gvectors, ((1),(1))) * units)

  # finally computes profile in fourier space.
  results = []
  for i, wfn in enumerate(extract.rwfns):
    if i not in indices: continue
    # compute fourier transform of real-space density.
    autocorr = dot(fourier.T, wfn.density)
    # now performs summation with gaussian in fourier space.
    results.append(dot(gaussians, autocorr))

  return x, array(results)

def linear_rprofile(extract, direction=(0,0,1), nbpoints=20, sigma=0.2, indices=None):
  """ Computes profile for given direction of each wavefunction.
  
      :Parameters:
        extract 
          Extraction object from an escan calculation. 
        nbpoints 
          Number of points to consider along the direction.
        sigma 
          Smearing to use when creating profile. 
        indices 
          Indices of the wavefunctions for which to compute the profile.
          Can also be None, in which case the profile of each wavefunction is
          computed.
      
      Mostly for debugging purposes. Makes much uglier plots than linear_profile.
  """
  from numpy import dot, max, arange, array, exp, multiply, sum
  from numpy.linalg import inv, norm
  from quantities import angstrom
  from lada.physics import a0

  direction = array(direction, dtype="float64")
  assert norm(direction) > 1e-12, ValueError("Direction cannot be null.")
  assert hasattr(extract, "gwfns"), ValueError("extract does not seem to be an escan extraction object.")
  if indices is None: indices = range(extract.gwfns)
  assert nbpoints > 1, ValueError("The number of points should be strictly larger than 1.")

  # first computes intersection of direction and cell.  At the end, we should
  # end up in a0, knowing the structure cell should be in angstrom.
  direction = array(direction)
  udir = direction / norm(direction)
  dir = dot(inv(extract.structure.cell), direction)
  direction = direction/max(abs(dir)) * extract.structure.scale * angstrom.rescale(a0)

  # then creates array of points along direction.
  zpoints = arange(nbpoints+1, dtype="float64").reshape(nbpoints+1, 1) / float(nbpoints)
  zpoints = multiply(zpoints, direction) 
  x =  sum(zpoints*udir, axis=1).rescale(angstrom) 
  if sigma is None: sigma = 0.75 * norm(zpoints[1] - zpoints[0]) * zpoints.units
  if not hasattr(sigma, "units"): sigma *= angstrom.rescale(a0)

  zpoints = sum(zpoints*udir, axis=1) 
  gaussians = array([exp(-(sum(extract.rvectors*udir, axis=1)-z)**2/sigma**2) for z in zpoints])
  results = []
  for i, wfn in enumerate(extract.rwfns):
    if i not in indices: continue
    results.append( dot(gaussians, wfn.density) )
    results[-1][-1] = results[-1][0]
  return x, array(results)
