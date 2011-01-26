import pickle
import matplotlib.pyplot as plt 
from lada.escan import ExtractBS

def compute_projections(extract, filename, alpha = 1e0, **kwargs):
  """ Computes projected densities around each atom for a calculation. """
  from os.path import exists
  from sys import exit, getrefcount
  if exists(filename): 
    from pickle import load
    with open(filename, "r") as file: return load(file)

  from pickle import dump
  from numpy import zeros, dot, array, exp, sum
  from numpy.linalg import norm
  from quantities import angstrom
  from lada import periodic_table as table
  from lada.crystal import gaussian_projector
  species = set([u.type for u in extract.structure.atoms])
  species = sorted(list(species))
  result = {}
  for key in species:
    if (key not in kwargs) and (key not in table.__dict__): continue
    radius = kwargs.get(key, getattr(table, key)).atomic_radius
    sigma = -alpha / radius / radius

    # create projection operator.
    proj = zeros(extract.rvectors.shape[:-1])
    cell = extract.structure.cell * extract.structure.scale 
    for atom in extract.structure.atoms:
      if atom.type != key: continue
      pos = atom.pos * extract.structure.scale 
      proj += gaussian_projector(extract.rvectors, pos * angstrom, cell * angstrom, sigma )
    result[key] = [(w.eigenvalue, w.expectation_value(proj)) for w in extract.rwfns]
  
  n = norm(array(result.values()))
  for key in result.keys(): result[key] /= n
  with open(filename, "w") as file: dump(result, file)
  return result

def compute_bs():
  import pickle
  from sys import exit
  from os.path import join
  from numpy import matrix, array
  from numpy.linalg import norm
  from boost.mpi import world
  from lada.opt import read_input
  from lada.escan import Escan, soH, band_structure
  from lada.vff import Vff
  from lada.crystal import fill_structure, sort_layers, FreezeCell, nb_valence_states  

  # reads input file.
  global_dict={"Vff": Vff, "Escan": Escan, "nb_valence_states": nb_valence_states, "soH": soH}
  input = read_input("input.py", global_dict=global_dict)

  # creating unrelaxed structure.
  structure = input.vff.lattice.to_structure()
  structure.atoms[0].type = "Si"
  structure.atoms[1].type = "Ge"
  structure.scale = 5.65

  # some kpoints + associated name
  X = array( [1,0,0], dtype="float64" )
  G = array( [0,0,0], dtype="float64" )
  L = array( [0.5,0.5,0.5], dtype="float64" )
  W = array( [0, 0.5,1], dtype="float64" )

  # Each job is performed for a given kpoint (first argument), at a given
  # reference energy (third argument). Results are stored in a specific directory
  # (second arguement). The expected eigenvalues are given in the fourth argument.
  kpoints = [ (X, G), (G, L) ]
  density = 20 / min( norm(X), norm(L), norm(W) )

  result = band_structure( input.escan, structure, kpoints, density, 
                           outdir = "results",
                           eref   = None, 
                           nbstates = nb_valence_states(structure) + 4,
                           pools = 4)
    
  if world.rank == 0:
    with open(join("results", "pickle"), "w") as file:
      pickle.dump(result, file) 


def create_data():
  from lada.escan import ExtractBS
  from os.path import join
  compute_bs()
  extract_bs = ExtractBS("results")
  for key, value in extract_bs.iteritems():
    filename = join(value.directory, "PROJECT_BS")
    a = compute_projections(value, filename)


def plot_bands(extractor, tolerance=1e-6, **kwargs):
  """ Tries and plots band-structure. """
  from os.path import join
  from numpy import dot, array, min, max, sqrt
  from numpy.linalg import norm

  bandcolor = kwargs.pop('bandcolor', 'blue')
  edgecolor = kwargs.pop('edgecolor', 'red')
  edgestyle = kwargs.pop('edgestyle', '--')
  awidth    = kwargs.pop('awidth', 5.)
  type      = kwargs.pop('type', 'Si')
  alpha     = kwargs.pop('alpha', 1.)

  # first finds breaking point.
  kpoints = extractor.kpoints
  delta = kpoints[1:] - kpoints[:-1]
  norms = [norm(delta[i,:]) for i in range(delta.shape[0])]
  bk = []
  for i, d in enumerate(norms[1:]):
    if abs(norms[i]-d) > 1e-6: bk.append(i+1)

  # then create array of mean x, y values.
  x = array([sum(norms[:i]) for i in range(len(norms)+1)])
  y = array(extractor.eigenvalues.values())

  # gets projected stuff
  si_projs = array([ compute_projections(value, join(value.directory, "PROJECT_BS"), alpha=alpha)["Si"]
                     for value in extractor.itervalues() ] )
  si_projs = si_projs[:,:,1]
  ge_projs = array([ compute_projections(value, join(value.directory, "PROJECT_BS"), alpha=alpha)["Ge"]
                     for value in extractor.itervalues() ] )
  ge_projs = ge_projs[:,:,1]

  # then loop over each band.
  for band, si, ge in zip(y.T, si_projs.T.real, ge_projs.T.real):
    width = (si if type == "Si" else ge)  / sqrt(si * si + ge * ge)
    plt.fill_between(x, band + awidth * width, band - awidth * width,
                     color=bandcolor, **kwargs)
  
  # plots break-lines.
  for i in bk: plt.axvline(x[i], color='black', **kwargs)

  # then plot vbm and cbm.
  kwargs.pop('linestyle', None) 
  plt.axhline(extractor.vbm, color=edgecolor, linestyle=edgestyle, **kwargs)
  plt.axhline(extractor.cbm, color=edgecolor, linestyle=edgestyle, **kwargs)



  plt.xlim((x[0], x[-1]))
  ylims = min(y) - (max(y) - min(y))*0.05, max(y) + (max(y) - min(y))*0.05
  plt.ylim(ylims)
  plt.show()

from lada.escan import ExtractBS
extract_bs = ExtractBS("results")
# p = compute_projections(extract_bs.values()[0], "shit", alpha=5.)
plot_bands(extract_bs, awidth=1)
