""" Creates job-dictionary of GA jobs. """
def gajobs(path, inputpath = "input.py", collector=None):
  from weakref import proxy

  import IPython.ipapi
  from numpy import array, arange, any, isnan
  from quantities import eV

  from lada.ga.functional import maximize, minimize
  from lada.ga.ce.functional import Darwin as Functional
  from lada.ga.ce.evaluator  import LocalSearchEvaluator as Evaluator
  from lada.jobs import JobDict
  from lada.escan import MassExtract
  from lada.ce import read_input

  # reads input file.
  input = read_input(inputpath, {"maximize": maximize, "minimize": minimize, "eV": eV, "arange": arange})

  if collector == None: collector = MassExtract(input.directory)
  extractor  = [extract for extract in collector.values() if extract.success]
  structures = [extract.input_structure for extract in extractor]
  for e in extractor: assert not any(isnan(e.eigenvalues))
  dos_values = array([extract.dos(input.energies, input.sigma) for extract in extractor]).T
  assert not any(isnan(dos_values))
  evaluator  = Evaluator(input.lattice, structures, energies=dos_values[0], **input.clusters)


  jobdict = JobDict()
  for energy, energies in zip(input.energies, dos_values):
    for trial in input.trials:
      kwargs = { "popsize": input.population_size, "rate": input.offspring_rate,
                 "max_gen": input.max_generations, "crossover_rate": input.crossover_rate, 
                 "history": input.history, "histlimit": input.histlimit,
                 "comparison": input.comparison, "rootworkdir": "$SCRATCH",
                 "alwayson": input.alwayson, "alwaysoff": input.alwaysoff }
      functional = Functional(evaluator, **kwargs)
      evaluator.darwin = proxy(functional)

      gajob = jobdict / "dos_{0}/trial_{1}".format(energy, trial)
      print "dos_{0}/trial_{1}".format(energy, trial)
      gajob.functional = functional 
      gajob.jobparams["energies"] = energies

  ip = IPython.ipapi.get()
  ip.user_ns["current_jobdict"] = jobdict
  ip.magic("savejobs " + path)
