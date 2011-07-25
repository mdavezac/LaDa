""" Creates job-dictionary of GA jobs. """
def gajobs(path, inputpath = "input.py"):
  from copy import deepcopy
  import IPython.ipapi

  from lada.ga.ce.functional import Darwin as Functional
  from lada.ga.ce.evaluator  import LocalSearchEvaluator as Evaluator
  from lada.jobs import JobDict
  from lada.escan import MassExtract
  from lada.ce import read_input, cluster_factory

  # reads input file.
  input = read_input(inputpath)

  extractor  = [extract for extract in MassExtract(input.directory) if extract.success]
  structures = [extract.structure for extract in extractor]
  dos        = array([extract.dos(extract, input.sigma) for extract]
  evaluator  = Evaluator(lattice, structure, dos[0], **input.clusters)




  jobdict = JobDict()
  for energy in input.energies:
    for trial in input.trials:
      escan = deepcopy(input.escan)
      escan.vff.lattice.scale = scale

      converter = Converter(growth=input.growth_direction, lattice=escan.vff.lattice)
      
      kwargs = { "popsizse": input.population_size, "rate": input.offspring_rate,
                 "max_gen": input.max_generations, "pools": input.pools,
                 "crossover_rate": input.crossover_rate, "swap_rate": input.swap_rate, 
                 "growth_rate": input.growth_rate, "nmin": nmin, "nmax": nmax,
                 "dosym": input.dosym, "rootworkdir": "$SCRATCH" }
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
