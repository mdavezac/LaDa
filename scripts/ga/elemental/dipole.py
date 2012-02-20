""" Creates job-dictionary of GA jobs. """
def gajobs(path, inputpath = "input.py"):
  from copy import deepcopy
  import IPython.ipapi

  from lada.ga.escan.elemental.functional import Darwin as Functional
  from lada.ga.escan.elemental.evaluator  import Dipole as Evaluator
  from lada.ga.escan.elemental            import Converter
  from lada.jobs import JobDict
  from lada.escan import read_input

  # reads input file.
  input = read_input(inputpath)

  jobdict = JobDict()
  for scale in input.scales:
    for range in input.ranges: 
      nmin, nmax = min(range), max(range)
      for trial in input.trials:
        escan = deepcopy(input.escan)
        escan.vff.lattice.scale = scale

        converter = Converter(growth=input.growth_direction, lattice=escan.vff.lattice)
        
        kwargs = { "popsizse": input.population_size, "rate": input.offspring_rate,
                   "max_gen": input.max_generations, "pools": input.pools,
                   "crossover_rate": input.crossover_rate, "swap_rate": input.swap_rate, 
                   "growth_rate": input.growth_rate, "nmin": nmin, "nmax": nmax,
                   "dosym": input.dosym, "rootworkdir": input.rootworkdir }
        kwargs.update(getattr(input, 'kwargs', {}))
        evaluator = Evaluator( converter = converter, 
                               functional = escan, 
                               **getattr(input, 'kwargs', {}) )
        functional = Functional(evaluator, **kwargs)

        gajob = jobdict / "{4[0]}{4[1]}{4[2]}/scale_{0:.2f}/{1}_{2}/trial_{3}"\
                          .format(scale, nmin, nmax, trial, input.growth_direction)
        gajob.functional = functional 

  ip = IPython.ipapi.get()
  ip.user_ns["current_jobdict"] = jobdict
  ip.magic("savejobs " + path)
