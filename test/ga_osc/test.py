""" Creates job-dictionary of GA jobs. """
from copy import deepcopy
import IPython.ipapi

from lada.ga.escan.elemental.functional import Darwin as Functional
from lada.ga.escan.elemental.evaluator  import Dipole as Evaluator
from lada.ga.escan.elemental            import LayeredConverter as Converter
from lada.jobs import JobDict
from lada.escan import Escan, soH
from lada.vff import Vff
from lada.opt import read_input

dictionary = {"Vff": Vff, "Escan": Escan, "soH": soH}
paths = ["workdir", "outdir"]
input = read_input("input.py", dictionary, paths=paths)


jobdict = JobDict()
for scale in input.scales:
  for p in input.periods: 
    for trial in input.trials:
      supercell = input.supercell(p)

      escan = deepcopy(input.escan)
      escan.fft_mesh = input.fftmesh(supercell)
      escan.vff.lattice.scale = scale

      converter = Converter(supercell, lattice=escan.vff.lattice)
      
      evaluator = Evaluator(converter = converter, escan = escan, **input.kwargs)
      evaluator._outdir.envvar = "$SCRATCH"

      functional = Functional(evaluator)
      functional.popsize     = input.population_size
      functional.max_gen     = input.max_generations
      functional.rate        = input.offspring_rate
      functional.mean_conc   = input.mean_conc
      functional.stddev_conc = input.stddev_conc
      functional.pools       = input.pools

      gajob = jobdict / "scale_{0:.2f}/period_{1}/trial_{2}".format(scale, p, trial)
      gajob.functional = functional 

ip = IPython.ipapi.get()
ip.user_ns["current_jobdict"] = jobdict
